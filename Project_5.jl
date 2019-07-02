# TODO: convert mod() to mod1()


using Plots
gr()

const L = 1. # System size
const c = 1. # wave speed

# t = timestep, n = number of grid points, mthd = integration method, tt = total time,
# pbc = periodic BCs (true/false), fq = frequency (for Dirichlet).
function advect(t, n; mthd, tt=1.0, tx=false, pbc=false, fq=10*pi)
    # Grid spacing.
    h = L/n

    k = c*t/(2*h) # Coefficient for all schemes.
    klw = 2*k^2 # Coefficient for Lax-Wendroff.

    # Initial conditions; integrate to total time.
    steps = Int(round(tt/t))
    
    sig = 0.1 # Width of Gaussian pulse.
    wn = pi/sig # Wave number.

    # x- and time-values for plotting.
    xvals = [h*i - L/2 for i in 0:(n-1)]
    tvals = [t*i for i in 0:(steps-1)]
    
    wv = zeros(steps, n)
    
    # Initial wave pulse if periodic BCs selected.
    if pbc 
        for (i,x) in enumerate(xvals)
            wv[1,i] = exp(-x^2/(2*sig^2)) * cos(wn*x)
        end
    end
    
    if tx
        xi = Int(ceil(n/8))
        wv[1, n-xi] = 1.0/(2*h)
    end

    return (mthd(wv, n, steps, k, pbc, fq, t, klw, h), xvals, tvals, n, t)
end

function ftcs(x...)
    wv, n, steps, k, pbc, fq, t = x

    for i in 2:steps
        if pbc
            # Periodic BCs
            for j in 1:n
                lft, rt = mod1(j-1,n), mod1(j+1,n)
                wv[i,j] = wv[i-1,j] - k*(wv[i-1,lft] - wv[i-1,rt])
            end
        else
            # Dirichlet BCs
            wv[i,1] = sin(t*i*fq)
            for j in 2:(n-1)
                wv[i,j] = wv[i-1,j] - k*(wv[i-1,j-1] - wv[i-1,j+1])
            end
            wv[1,n] = 0
        end
    end
    
    return wv
end

function lax(x...)
    wv, n, steps, k, pbc, fq, t = x

    for i in 2:steps
        if pbc
            # Periodic BCs
            for j in 1:n
                lft, rt = mod1(j-1,n), mod1(j+1,n)
                wv[i,j] = 0.5*(wv[i-1,rt] + wv[i-1,lft]) - k*(wv[i-1,rt] - wv[i-1,lft])
            end
        else
            # Dirichlet BCs
            wv[i,1] = sin(t*i*fq)
            for j in 2:(n-1)
                wv[i,j] = 0.5*(wv[i-1,j+1] + wv[i-1,j-1]) - k*(wv[i-1,j+1] - wv[i-1,j-1])
            end
            wv[1,n] = 0
        end
    end
    
    return wv
end

function wdf(x...)
    wv, n, steps, k, pbc, fq, t, klw = x

    for i in 2:steps
        if pbc
            # Periodic BCs
            for j in 1:n
                lft, rt = mod1(j-1,n), mod1(j+1,n)
                wv[2,j] = wv[1,j] - k*(wv[1,rt] - wv[1,lft]) + klw*(wv[1,rt] + wv[1,lft] - 2*wv[1,j])
            end
        else
            # Dirichlet BCs
            wv[i,1] = sin(t*i*fq)
            for j in 2:(n-1)
                wv[i,j] = wv[i-1,j] - k*(wv[i-1,j+1] - wv[i-1,j-1]) + klw*(wv[i-1,j+1] + wv[i-1,j-1] - 2*wv[i-1,j])
            end
            wv[i,n] = 0
        end
    end
    
    return wv
end

function upwind(x...)
    wv, n, steps, k, pbc, fq, t, klw = x
    
    for i in 2:steps
        if pbc
            # Periodic BCs
            for j in 1:n
                wv[i,j] = wv[i-1,j] - 2*k*(wv[i-1,j] - wv[i-1,mod(j-2,n)+1])
            end
        else
            # Dirichlet BCs
            wv[i,1] = sin(t*i*fq)
            for j in 2:(n-1)
                wv[i,j] = wv[i-1,j] - 2*k*(wv[i-1,j] - wv[i-1,j-1])
            end
            wv[i,n] = 0
        end
    end
    
    return wv
end

function leapfrog(x...)
    wv, n, steps, k, pbc, fq, t, klw = x
    
    # One iteration of Lax-Wendroff to get leapfrog started.
    for j in 1:n
        lft, rt = mod1(j-1,n), mod1(j+1,n)
        wv[2,j] = wv[1,j] - k*(wv[1,rt] - wv[1,lft]) + klw*(wv[1,rt] + wv[1,lft] - 2*wv[1,j])
    end
    
    for i in 3:steps
        if pbc
            # Periodic BCs
            for j in 1:n
                lft, rt = mod1(j-1,n), mod1(j+1,n)
                wv[i,j] = wv[i-2,j] - 2*k*(wv[i-1,rt] - wv[i-1,lft])
            end
        else
            # Dirichlet BCs
            wv[i,1] = sin(t*i*fq)
            for j in 2:(n-1)
                wv[i,j] = wv[i-2,j] - 2*k*(wv[i-1,j+1] - wv[i-1,j-1])
            end
            wv[i,n] = 0
        end
    end
    
    return wv
end

function transport(x...)
    wv, n, steps, k, pbc, fq, t, klw, h = x
    
    # Coefficients for discretized transport eqn.
    k1, k2 = -c*t/h, k*t/h^2
    p, q = k2-k1, k1-2*k2+1
    
    s1, s2 = (c*t/h)^2, (2*k*t)/h^2
    
    if s1 <= s2 && s2 <= 1
        println("Solution expected to be stable:")
        println("(c*t/h)^2 = $s1 , (2*k*t)/h^2 = $s2")
    else
        println("Solution expected to be unstable")
        println("(c*t/h)^2 = $s1 , (2*k*t)/h^2 = $s2")
    end    
    
    for i in 2:steps
        for j in 1:n
            lft, rt = mod1(j-1,n), mod1(j+1,n)
            a, b, c = wv[i-1,rt], wv[i-1,j], wv[i-1,lft]
            wv[i,j] = a*p + b*q + c*k2
        end
    end
    
    return wv
end     

function plt_dat(dat)
    wv, x, t, nstep, tstep = dat
    ttl = string("Advection for dt = ", string(tstep), " s, N = ", string(nstep))
    p1=surface(x, t, wv, title=ttl,
            xlabel=:"x", ylabel=:"time", zlabel=:"a(x,t)", 
            color=:Spectral_r, colorbar=false,
            st=:surface, camera=(25,50))
    p2=plot(x, wv[1,:], label="Initial")
    plot!(p2, x, wv[end,:], label="Final")
    plot(p1, p2, layout=2, size=(950,450))
end

function dir_gif(dat)
    wv, x, t, nstep, tstep = dat
    
    max_amp = 1.25
    min_amp = -1.25
    max_x = maximum(x)
    min_x = minimum(x)

    anim = @animate for i=1:2:Int(round(size(wv,1)))
        dt = "$(round(t[i], digits=4))"
        if length(dt) < 6
            dt = rpad(dt, 6, "0")
        end
        ttl = string("Sinusoidal dirichlet BC at time t = ", dt, " s")
        plot(size=(900,500), xlabel=:"x", ylabel=:"a(x)", xlims=(min_x, max_x), ylims=(min_amp, max_amp), title=ttl)
        plot!(x, wv[i, :], lw=4, label=:"")
    end every 5
    gif(anim, "advect.gif")
end

# Plots for exercise 7.2
plt_dat(advect(0.015, 50, mthd=ftcs))
plt_dat(advect(0.02, 50, mthd=ftcs))
plt_dat(advect(0.03, 50, mthd=ftcs))

plt_dat(advect(0.015, 50, mthd=lax))
plt_dat(advect(0.02, 50, mthd=lax))
plt_dat(advect(0.03, 50, mthd=lax))

plt_dat(advect(0.015, 50, mthd=wdf))
plt_dat(advect(0.02, 50, mthd=wdf))
plt_dat(advect(0.03, 50, mthd=wdf))

# Plots for exercise 7.4
dat = advect(0.002, 50, tt=1., tx=true, mthd=transport)
wv, x, t = dat
plt_dat(dat)
contour(x, t, wv, title="Contour levels for advection of delta spike",
            xlabel=:"x", ylabel=:"time", size=(900,900),
            levels=[i^2 for i in 0.:0.1:sqrt(maximum(wv))],
            color=:Spectral_r, colorbar=false,
            st=:contourf, contour_labels=true)

plt_dat(advect(0.002825, 50, tt=1., tx=true, mthd=transport))

# Plots for exercise 7.5
plt_dat(advect(0.015, 50, mthd=upwind))
plt_dat(advect(0.02, 50, mthd=upwind))
plt_dat(advect(0.03, 50, mthd=upwind))

# Plots for exercise 7.6
plt_dat(advect(0.015, 50, mthd=leapfrog))
plt_dat(advect(0.02, 50, mthd=leapfrog))
plt_dat(advect(0.03, 50, mthd=leapfrog))