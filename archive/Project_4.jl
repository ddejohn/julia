using Plots
gr()

# Length of rod and diffusivity constants.
const L = 1.
const diff = 1.

# Dirichlet Forward Time Centered Space (ends fixed at 0°C).
function dftcs(x...)
    temps, n, steps, k = x

    for i in 2:steps
#         temps[i,1] = 0
#         temps[i,n] = 0
        for j in 2:(n-1)
            temps[i,j] = temps[i-1,j] + k*(temps[i-1,j-1] + temps[i-1,j+1] - 2*temps[i-1,j])
        end
    end
    return temps
end

# Neumann Forward Time Centered Space (zero flux at ends).
function nftcs(x...)
    temps, n, steps, k = x
    
    for i in 2:steps
        for j in 2:(n-1)
            temps[i,j] = temps[i-1,j] + k*(temps[i-1,j-1] + temps[i-1,j+1] - 2*temps[i-1,j])
        end
        temps[i,1] = temps[i,2]
        temps[i,n] = temps[i,n-1]
    end
    return temps
end

# Periodic Forward Time Centered Space (flux equal at both ends).
function pftcs(x...)
    temps, n, steps, k = x
    
    for i in 2:steps
        temps[i,1] = temps[i-1,1] + k*(temps[i-1,n] + temps[i-1,2] - 2*temps[i-1,1])
        
        for j in 2:(n-1)
            temps[i,j] = temps[i-1,j] + k*(temps[i-1,j-1] + temps[i-1,j+1] - 2*temps[i-1,j])
        end
        
        temps[i,n] = temps[i-1,n] + k*(temps[i-1,n-1] + temps[i-1,1] - 2*temps[i-1,n])
    end
    return temps
end

# DuFort-Frankel Method (not working).
function dffm(x...)
    temps, n, steps, k, t, h = x
    s = 2*diff*t/h^2
    
    # First iteration (DFFM starter).
    for j in 2:(n-1)
        temps[2,j] = temps[1,j] + (diff*1e-6/h^2)*(temps[1,j-1] + temps[1,j+1] - 2*temps[1,j])
    end
    
    # DFFM
    for i in 3:steps
        for j in 2:(n-1)
            temps[i,j] = ((1-s)/(1+s)) * temps[i-2,j] + (s/(1+s)) * (temps[i-1,j+1] + temps[i-1,j-1])
        end
    end
    
    return temps
end

# Compute temperature in 1D rod of length 'const L' and diffusivity 'const diff'.
# Where t = timestep, n = number of grid points, mthd = dftcs, nftcs, pftcs, dffm.
# (Dirichlet, Neumann, periodic, DuFort-Frankel)
function heat_eq(t, n; mthd, offset=1/2)
    # Grid spacing.
    h = L/(n-1)

    # FTCS coefficient.
    k = diff*t/h^2

    # Initial conditions; integrate to t=0.025 seconds.
    steps = Int(round(0.025/t))
    xi = Int(ceil(offset*n))
    temps = zeros(steps, n)
    temps[1, n-xi] = 1.0/h

    # x- and time-values for plotting
    xvals = [-L/2 + h*i for i in 0:(n-1)]
    tvals = [t*i for i in 0:(steps-1)]

    return (mthd(temps, n, steps, k, t, h), xvals, tvals, n, t)
end

function plot_temps(dat, inc=0.1)
    temp, x, t, nstep, tstep = dat
    ttl = string("Diffusion of a delta spike for dt = ", string(tstep), " s, N = ", string(nstep))
    p1=surface(x, t, temp, title=ttl,
            xlabel=:"x", ylabel=:"time", zlabel=:"T(x,t)", 
            color=:Spectral_r, colorbar=false,
            st=:surface, camera=(45,45))
    p2=contour(x, t, temp, title="Temperature levels ($(inc)°C increments)",
            xlabel=:"x", ylabel=:"time", 
            levels=[i^2 for i in 0.:inc:sqrt(maximum(temp))],
            color=:Spectral_r, colorbar=false,
            st=:contourf, contour_labels=true)
    plot(p1, p2, layout=2, size=(950,450))
end


#plot_temps(heat_eq(timestep, gridpoints, method), contour_plot_increment)
# EX:
plot_temps(heat_eq(1e-5, 61, mthd=dftcs), 0.15)

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#

# # GIF-MAKING STUFF
# function rand_eq(t, n; mthd)
#     # Grid spacing.
#     h = L/(n-1)

#     # FTCS coefficient.
#     k = diff*t/h^2

#     # Initial conditions; integrate to t=0.025 seconds.
#     steps = Int(round(0.025/t))
#     temps = zeros(steps, n)

#     # x- and time-values for plotting
#     xvals = [-L/2 + h*i for i in 0:(n-1)]
#     tvals = [t*i for i in 0:(steps-1)]
    
#     for (i,x) in enumerate(xvals)
#         temps[1, i] = 6*cos(7*x)*sin(-17*x) / (cos(13*x) + 3*sin(-37*x)^3 + 5) + 3
#     end

#     return (mthd(temps, n, steps, k, t, h), xvals, tvals)
# end

# function all_three()
#     ddat = heat_eq(1e-5, 41, mthd=dftcs, offset=1/8)
#     ndat = heat_eq(1e-5, 41, mthd=nftcs, offset=1/8)
#     pdat = heat_eq(1e-5, 41, mthd=pftcs, offset=1/8)

#     dtemp, dx, t = ddat
#     ntemp, nx, nt = ndat
#     ptemp, px, pt = pdat
#     max_T = maximum(dtemp)/2
#     min_x = minimum(dx)
#     max_x = maximum(dx)

#     anim = @animate for i=1:5:Int(round(size(dtemp,1)))
#         dt = "$(round(t[i], digits=4))"
#         if length(dt) < 6
#             dt = rpad(dt, 6, "0")
#         end
#         ttl = string("Temperature from a delta spike at L/8 at time t = ", dt, " s")
#         plot(size=(900,450), title=ttl, xlabel=:"x", ylabel=:"T(x)", xlims=(min_x, max_x), ylims=(0.,max_T))
#         plot!(dx, dtemp[i, :], lw=4, label=:"Dirichlet")
#         plot!(nx, ntemp[i, :], lw=4, label=:"Neumann")
#         plot!(px, ptemp[i, :], lw=4, label=:"Periodic")
#     end every 5
#     gif(anim, "all_three.gif")
# end

# function rand_temps()
#     temp1, x1, t1 = rand_eq(1e-6, 321, mthd=dftcs)
#     temp2, x2, t2 = rand_eq(1e-6, 321, mthd=nftcs)
#     temp3, x3, t3 = rand_eq(1e-6, 321, mthd=pftcs)
    
#     max_T = 1.25*maximum(temp1)
#     min_T = minimum(temp1)
#     min_x = minimum(x1)
#     max_x = maximum(x1)

#     anim = @animate for i=1:30:Int(round(size(temp1,1)))
#         dt = "$(round(t1[i], digits=4))"
#         if length(dt) < 6
#             dt = rpad(dt, 6, "0")
#         end
#         ttl = string("Sinusoidal temperature distribution at time t = ", dt, " s")
#         plot(size=(900,450), xlabel=:"x", ylabel=:"T(x)", xlims=(min_x, max_x), ylims=(min_T,max_T), title=ttl)
#         plot!(x1, temp1[i, :], lw=4, label=:"Dirichlet")
#         plot!(x2, temp2[i, :], lw=4, label=:"Neumann")
#         plot!(x3, temp3[i, :], lw=4, label=:"Periodic")
#         plot!(x1, temp1[1, :], lw=4, alpha=0.5, label=:"Initial")
#     end every 5
#     gif(anim, "rand_temps.gif")
# end;