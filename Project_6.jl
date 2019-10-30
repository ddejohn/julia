using Plots
using CurveFit
import GR
gr()

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Main programs

function analytic2d(n; gp=50, V=1., L=1.)
    h = L/(gp-1)
    
    srf = zeros(gp, gp)
    
    for i in 1:gp
        srf[i,:] = [(V*4/pi)*sum((sin(k*pi*(j-1)*h/L)*sinh(k*pi*(i-1)*h/L))/(k*sinh(k*pi))
                        for k in 1:2:n) for j in 1:gp]
    end
    
    return (srf, [x for x in 0:h:L], [y for y in 0:h:L])
end

function analytic3d(; n, z, gp=50, V=1., L=1.)
    h = L/(gp-1)    
    vol = zeros(gp,gp)
    
    # Fourier double sum
    for p in 1:2:n
        for q in 1:2:n
            for i in 1:gp
                for j in 1:gp
                    cnm = 16*V*sech(pi*sqrt(p^2 + q^2)/2)/(p*q*pi^2)
                    xx = sin(p*pi*(i-1)*h/L)
                    yy = sin(q*pi*(j-1)*h/L)
                    zz = cosh(pi*sqrt(p^2 + q^2)*(z - L/2)/L)
                    vol[i,j] += cnm*xx*yy*zz
                end
            end
        end
    end
    
    return (vol, [x for x in 0:h:L], [y for y in 0:h:L])
end
    
# Calculate potential in a cube.
# mthd = integrator, gp = # of grid points, int = interior values function.
function laplace(; mthd, gp=50, V=1., L=1., intr=zeroint)
    h = L/(gp-1)
    
    # Initialize our cube with points defined by an intr function. Defaults to all zeros.
    vol = intr(gp)
    
    # Boundary conditions (the other boundary conditions are satisfied by intr).
    vol[:, :, 1] .= V
    vol[:, :, end] .= V
    
    # For plotting.
    xx = [x for x in 0:h:L]
    yy = [y for y in 0:h:L]
    zz = [z for z in 0:h:L]
    
    # Calculate the potential using given mthd.
    v, imax = mthd(vol)
    return (v, xx, yy, zz, uppercase(string(mthd)), imax)
end

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Helper functions

# Zeros on the interior points.
zeroint(gp) = zeros(gp,gp,gp)

# Interior points randomized.
function randint(gp)
    vol = rand(gp,gp,gp)
    for z in 2:(gp-1)
        vol[1,:,z] .= 0.
        vol[:,1,z] .= 0.
        vol[end,:,z] .= 0.
        vol[:,end,z] .= 0.
    end
    return vol
end

# Interior points given by first 5 terms of analytic soln.
function analyticint(gp)
    vol = zeros(gp,gp,gp)
    for i in 2:(gp-1)
        zz, xx, yy = analytic3d(n=5, z=(i/gp), gp=gp)
        vol[:,:,i] = zz
        vol[1,:,i] .= 0.
        vol[:,1,i] .= 0.
        vol[end,:,i] .= 0.
        vol[:,end,i] .= 0.
    end
    return vol
end

# Volume in potential cube at z=lvl. Helper function for gradient.
function dim3(args...)
    gp, V, h, lvl = args
    
    # Initialize our cube with points defined by an intr function.
    vol = analyticint(gp)
    
    # Boundary conditions.
    vol[:, :, 1] .= V
    vol[:, :, end] .= V
    
    # z-slice of volume.
    zlvl = Int(round(lvl/h))
    
    # Potential.
    v, = sor(vol)
    
    return v[:,:,zlvl]
end

# Volume on square surface. Helper function for gradient. Constant potential along edge.
function dim2(args...)
    gp, V, = args
    s = zeros(gp,gp)
    s[2:(end-1), end] .= V
    
    return sor_srf(s)
end

# Volume on square surface. Helper function for gradient. Single charge on edge.
function poisson_dir(args...)
    gp, V, h, = args
    e0 = 8.8452e-12
    s = zeros(gp,gp)
    p = zeros(gp,gp)
    p[Int(round(gp/2)), Int(round(gp/2))] = 3e-7
    
    w = 2/(1+sin(pi/gp))
    div = (gp-2)^2
    
    for _ in 1:gp^2
        c_sum = 0
        
        for i in 2:(gp-1)
            for j in 2:(gp-1)
                temp = w/4*(s[i+1,j] + s[i-1,j] + s[i,j+1] + s[i,j-1] + (h^2/e0)*p[i,j]) + (1-w)*s[i,j]
                c_sum += abs(1-s[i,j]/temp)
                s[i,j] = temp
            end
        end

        if c_sum/div < 1e-4
            break
        end
    end
    
    return s
end

function poisson_per(args...)
    gp, V, h, = args
    e0 = 8.8452e-12
    s = zeros(gp,gp)
    p = zeros(gp,gp)
    p[Int(round(gp/2)), Int(round(gp/2))] = 3e-7
    
    w = 2/(1+sin(pi/gp))
    div = (gp-2)^2
    
    for _ in 1:gp^2
        c_sum = 0
        
        for i in 1:gp
            for j in 1:gp
                lft, rt = mod(j-2,gp)+1, mod(j,gp)+1
                up, dwn = mod(i,gp)+1, mod(i-2,gp)+1
                temp = w/4*(s[up,j] + s[dwn,j] + s[i,rt] + s[i,lft] + (h^2/e0)*p[i,j]) + (1-w)*s[i,j]
                c_sum += abs(1-s[i,j]/temp)
                s[i,j] = temp
            end
        end

        if c_sum/div < 1e-5
            break
        end
    end
    
    return s
end
    
# lvl = z height of desired vector field.
function gradient(; dim, lvl=0, scale=0.1)
    # Grid spacing.
    gp = 101
    V = 1.
    L = 1.
    h = L/(gp-1)
    
    v = dim(gp, V, h, lvl)
    
    # x- and y-coordinates of every gridpoint.
    stp = 10h:10h:(L-h)
    xy = [(x,y) for x in stp for y in stp]
    
    # Initialize arrays for x and y steepest descents (size is gp-2 because we iterate only over interior points)
    grad = []
  
    # Do every other grid point to cut down on redundant calcs.
    for i in 10:10:(gp-10)
        for j in 10:10:(gp-10)
            vij = v[i, j]   # V(i,j)
            
            vpi = v[i+1, j] # V(i+1, j)
            vmi = v[i-1, j] # V(i-1, j)
            
            vpj = v[i, j+1] # V(i, j+1)
            vmj = v[i, j-1] # V(i, j-1)
                        
            dxs = [(vpi - vij)/h, -(vmi - vij)/h]
            dys = [(vpj - vij)/h, -(vmj - vij)/h]
                        
            # Take the largest change in potential along x- and y-directions.
            dx = -scale*dxs[argmax(abs.(dxs))]
            dy = -scale*dys[argmax(abs.(dys))]

            push!(grad, (dx, dy))
        end
    end
    
    p = plot(size=(900,900), title="", label="",
                xlabel="x", ylabel="y", xlims=(0.,1.), ylims=(0.,1.),
                xticks=[i for i in 0.0:0.1:1.0], yticks=[i for i in 0.0:0.1:1.0])
                
    contour!(p, [x for x in 0:h:L], [y for y in 0:h:L], v', 
                title="Electric field vectors plotted over equipotential curves", 
                levels=[i for i in 0.:0.1:1.], st=:contourf,
                color=:Spectral_r, colorbar=false, contour_labels=true)
    
    quiver!(p, xy, quiver=(grad), lw=4, seriescolor=:white)
    p
end


#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Integrators

function jacobi(v)
    gp = size(v,1)
    div = (gp-2)^2
    imax = 0
    c = copy(v)
    
    for iter in 1:gp^2
        c_sum = 0
            
        for k in 2:(gp-1)
            for i in 2:(gp-1)
                for j in 2:(gp-1)
                    c[i,j,k] = 1/6*(v[i+1,j,k] + v[i-1,j,k] + 
                                    v[i,j+1,k] + v[i,j-1,k] + 
                                    v[i,j,k+1] + v[i,j,k-1])
                    
                    c_sum += abs(1-v[i,j,k]/c[i,j,k])
                end
            end
        end
        
        v = copy(c)
        imax = iter
       
        if c_sum/div < 1e-4
            break
        end
    end
    return (v,imax)
end

function gauss(v)
    gp = size(v,1)
    div = (gp-2)^2
    imax = 0
    
    for iter in 1:gp^2
        c_sum = 0
    
        for k in 2:(gp-1)
            for i in 2:(gp-1)
                for j in 2:(gp-1)
                    temp = 1/6*(v[i+1,j,k] + v[i-1,j,k] + 
                                v[i,j+1,k] + v[i,j-1,k] + 
                                v[i,j,k+1] + v[i,j,k-1])
                    
                    c_sum += abs(1-v[i,j,k]/temp)
                    v[i,j,k] = temp
                end
            end
        end
        
        imax = iter
        if c_sum/div < 1e-4
            break
        end
    end
    
    return (v,imax)
end

function sor(v)
    gp = size(v,1)
    w = 2/(1+sin(pi/gp))
    div = (gp-2)^2
    imax = 0
    
    for iter in 1:gp^2
        c_sum = 0
        
        for k in 2:(gp-1)
            for i in 2:(gp-1)
                for j in 2:(gp-1)
                    temp = w/6*(v[i+1,j,k] + v[i-1,j,k] + 
                                v[i,j+1,k] + v[i,j-1,k] + 
                                v[i,j,k+1] + v[i,j,k-1]) + (1-w)*v[i,j,k]
                    
                    c_sum += abs(1-v[i,j,k]/temp)
                    v[i,j,k] = temp
                end
            end
        end
        
        imax = iter
        if c_sum/div < 1e-4
            break
        end
    end
    
    return (v,imax)
end
                
function sor_srf(s)
    gp = size(s,1)
    w = 2/(1+sin(pi/gp))
    div = (gp-2)^2
    
    for _ in 1:gp^2
        c_sum = 0
        
        for i in 2:(gp-1)
            for j in 2:(gp-1)
                temp = w/4*(s[i+1,j] + s[i-1,j] + s[i,j+1] + s[i,j-1]) + (1-w)*s[i,j]
                c_sum += abs(1-s[i,j]/temp)
                s[i,j] = temp
            end
        end

        if c_sum/div < 1e-4
            break
        end
    end
    
    return s
end

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Plotting

function srf_plot(dat; ttl="", zslc="", zlm=(0.,1.))
    zz, xx, yy = dat
    l = @layout [a{0.1h} ; b{0.5w} c{0.5w}]
    st = string("V(x, y", zslc, ")")
        
    srf=surface(xx, yy, zz, title=st,
                xlabel="x", ylabel="y", zlabel="",
                xticks=0:0.5:1., yticks=0:0.5:1.,
                zlims=zlm,
                color=:Spectral_r, colorbar=false,
                st=:surface, camera=(45,45))

    con=contour(xx, yy, zz, title=string("Equipotential level curves"),
                xlabel="x", ylabel="y", 
                xticks=0:0.5:1., yticks=0:0.5:1.,
                levels=[i for i in 0.:0.1:1.],
                color=:Spectral_r, colorbar=false,
                st=:contourf, contour_labels=true)
    
    pt=plot(annotation = (0, 0, text(string(ttl,"\n"), 20, :left)), showaxis=:false, grid=:false)

    plot(pt, srf, con, layout=l, size=(950,550))
end

function vol_plot(dat; lvl)
    vol, xx, yy, zslices, mthd, imax = dat
    gp = size(vol,3)
    zz = vol[:,:,Int(round(gp*lvl)+1)]
    ttl = "Potential in a cube at z = $(lvl) using $(mthd), N = $(imax+1) iterations"
    srf_plot((zz,xx,yy), ttl=ttl, zslc=", z = $(lvl)")
end

# function vol_gif(dat)
#     vol, xx, yy, zslices, mthd, imax = dat
#     gp = size(vol,3)

#     anim = @animate for i in 2:(size(vol,1)-1)
#         zz = vol[:,:,i]
#         zi = "$(round(zslices[i], digits=2))"
#         if length(zi) < 4
#             zi = rpad(zi, 4, "0")
#         end
        
#         ttl = "Potential in a cube at z = $(zi); $(mthd) method, $(imax+1) iterations"
#         srf_plot((zz,xx,yy), ttl=ttl, zslc=", z = $zi")
#     end every 1
    
#     gif(anim, string("laplace_", mthd, ".gif"))
# end

# vol_gif(laplace(mthd=jacobi, gp=100))
# vol_gif(laplace(mthd=gauss, gp=100))
# vol_gif(laplace(mthd=sor, gp=100))

# Maximum emission projection suggested by Anshul Singhvi on Julia's Slack server.
function mep(; mthd, gp=50)
    v = zeros(gp, gp, gp)
    v[:,:,1] .= 1.
    v[:,:,end] .= 1.

    pv, numiter = mthd(v) # plot volume

    GR.title("Potential in a cube using $(uppercase(string(mthd))), with n=$(numiter+1)")
    GR.volume(pv, algorithm = 0, colormap = 47, size=(900,900))
end

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Exercise-specific methods

# 8.1
# TODO: calculate Fourier series only once;
# Pull relevant data at n=11, n=21, n=51, n=151 instead of running each of those sums individually.
function one_percent()
    nums=[11, 21, 51, 151]
    gp = 1000
    lft = Int(round(gp*0.4))
    rt = gp-lft+1
    
    hz = [i/gp for i in (lft-1):rt]
    v0 = [1. for i in (lft-1):rt]
    vp = [1.005 for i in (lft-1):rt]
    vm = [0.995 for i in (lft-1):rt]

    p = plot(title="Partial sum to n terms on $gp grid points",
             size=(900,900), yticks=0.9:0.01:1.1, legend=:bottomright)
    
    plot!(hz, v0, l=(3,:dash,:green), alpha=0.35, label="V_0")
    plot!(hz, vp, l=(3,:dash,:red), alpha=0.35, label="+0.5%")
    plot!(hz, vm, l=(3,:dash,:blue), alpha=0.35, label="-0.5%")

    for n in nums
        dat = analytic2d(n, gp=gp)
        srf, xx, yy = dat
        plot!(p, xx[lft:rt], srf[end, lft:rt], lw=3, label="n=$n")
    end

    p
end

# 8.3
function system_size(; mthd, V=1., L=1., intr)
    grid_sizes = [i for i in 10:5:50]
    num_iters = []
    for gp in grid_sizes
        vol = intr(gp)
        vol[:, :, 1] .= V
        vol[:, :, end] .= V
        v, imax = mthd(vol)
        push!(num_iters, imax)
    end
    
    return (grid_sizes, num_iters)
end

# 8.3
function sys_plot(; ttl, interior=zeroint)
    x1, y1 = system_size(mthd=jacobi, intr=interior)
    x2, y2 = system_size(mthd=gauss, intr=interior)
    x3, y3 = system_size(mthd=sor, intr=interior)
    
    a1, b1 = power_fit(x1, y1)
    a2, b2 = power_fit(x2, y2)
    a3, b3 = power_fit(x3, y3)
    
    xx = [x for x in 10.:0.1:50.]    
    p = plot(title=ttl, xlabel="Grid Size", ylabel="Iterations",
             xticks=[x for x in 10:5:50], yticks=[y for y in 100:200:2500], legend=:topleft, size=(900,900))
    
    plot!(p, xx, [a1*x^b1 for x in xx], lw=3, label="y = $(round(a1,digits=4)) * x^($(round(b1,digits=4)))")
    plot!(p, xx, [a2*x^b2 for x in xx], lw=3, label="y = $(round(a2,digits=4)) * x^($(round(b2,digits=4)))")
    plot!(p, xx, [a3*x^b3 for x in xx], lw=3, label="y = $(round(a3,digits=4)) * x^($(round(b3,digits=4)))")
    
    scatter!(p, x1, y1, markersize=6, label="Jacobi")
    scatter!(p, x2, y2, markersize=6, label="Gauss-Seidel")
    scatter!(p, x3, y3, markersize=6, label="SOR")
    
    p
end

# 8.3 SOR method only
function sor_plot(; ttl)
    x1, y1 = system_size(mthd=sor, intr=zeroint)
    x2, y2 = system_size(mthd=sor, intr=randint)
    x3, y3 = system_size(mthd=sor, intr=analyticint)
    
    a1, b1 = power_fit(x1, y1)
    a2, b2 = power_fit(x2, y2)
    a3, b3 = power_fit(x3, y3)
    
    xx = [x for x in 10.:0.1:50.]    
    p = plot(title=ttl, xlabel="Grid Size", ylabel="Iterations",
             xticks=[x for x in 10:5:50], yticks=[y for y in 0:10:300], legend=:topleft, size=(900,900))
    
    plot!(p, xx, [a1*x^b1 for x in xx], lw=3, label="y = $(round(a1,digits=4)) * x^($(round(b1,digits=4)))")
    plot!(p, xx, [a2*x^b2 for x in xx], lw=3, label="y = $(round(a2,digits=4)) * x^($(round(b2,digits=4)))")
    plot!(p, xx, [a3*x^b3 for x in xx], lw=3, label="y = $(round(a3,digits=4)) * x^($(round(b3,digits=4)))")
    
    scatter!(p, x1, y1, markersize=6, label="zeros")
    scatter!(p, x2, y2, markersize=6, label="random")
    scatter!(p, x3, y3, markersize=6, label="Fourier")
    
    p
end;