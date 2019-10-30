using Plots
using Printf
gr()

#——————————————————————————————————————————————————————————————————————————————————————————————————#
# Constants

# const μ = 1.3271244004193938e20 # JPL ephermerides, August 2018
# const AU = 1.495978707e11 # JPL ephermerides, August 2018
# const m_h = 2.2e14 # Halley's comet mass (wikipedia)

const μ = 4*pi^2

#——————————————————————————————————————————————————————————————————————————————————————————————————#
# Orbital state vector

struct State
    r::AbstractVector{Float64}
    v::AbstractVector{Float64}
end

# Define scalar multiplication and addition on state vectors.
Base.:*(k::Number, S::State) = State(k*S.r, k*S.v)
Base.:+(A::State, B::State) = State(A.r + B.r, A.v + B.v)
Base.:-(A::State, B::State) = State(A.r - B.r, A.v - B.v)
Base.abs(s::State) = State([abs(i) for i in s.r], [abs(i) for i in s.v])

# Just a container for holding orbital data for plotting purposes.
struct Orbit
    ρ; θ; K; U; E; Ers; t
end

struct Taus_Vs 
    v; τ
end

#——————————————————————————————————————————————————————————————————————————————————————————————————#
# Helper functions

function nrm(r::AbstractVector{Float64})
    return sqrt(sum(r.^2))
end

function nrm(s::State)
    return sqrt(sum(vcat(s.r, s.v).^2))
end

function gf(X::State)
    dr = X.v
    dv = -(μ/nrm(X.r)^3) * X.r
    return State(dr, dv)
end

function df(X::State)
end

function pf(X::State)
    r = nrm(X.r)
    v = nrm(X.v)
    α = 0.0005
    
    dr = X.v
    dv = -(μ/r^3)*(1-α/r)* X.r
    return State(dr, dv)
end

function energy(X::State)
    k = 1/2*nrm(X.v)^2
    u = -μ/nrm(X.r)
    return (K=k, U=u, E=k+u)
end

function abs_error(exp, theo)
    return abs(100*(exp - theo)/theo)
end

# function diff(X_s::State, X_b::State)
#     return abs((nrm(X_s.r) + nrm(X_s.v)) - (nrm(X_b.r) + nrm(X_b.v)))
# end

# Returns an interpolated polynomial lambda.
function intrp(xx, yy; n=5)
    k = length(xx)
    if k != length(Set(xx))
        println("Cannot interpolate from duplicate x-values!")
    end

    function lbp(j, x)
        p = 1
        for m in (1:k)[1:end .!= j]
            p *= (x - xx[m])/(xx[j] - xx[m])
        end
        return p
    end

    return x -> round(sum(yy[j]*lbp(j, x) for j in 1:k), digits=n)
end

#——————————————————————————————————————————————————————————————————————————————————————————————————#
# Integrators

function euler(X::State, τ::Float64)
    v = X.v -(τ*μ/nrm(X.r)^3) * X.r
    r = X.r + τ*X.v
    return State(r, v)
end

function cromer(X::State, τ::Float64)
    v = X.v -(τ*μ/nrm(X.r)^3) * X.r
    r = X.r + τ*v
    return State(r, v)
end

# function rk(X::State, τ::Float64)
#     return X + τ * f(X + τ/2 * f(X))
# end

function rk4(X::State, τ::Float64, cf=gf)
    k1 = cf(X)
    k2 = cf(X + τ/2*k1)
    k3 = cf(X + τ/2*k2)
    k4 = cf(X + τ*k3)

    return X + τ/6 * (k1 + k4 + 2*(k2 + k3))
end

function rka(X::State, τ::Float64, Δi::Float64, cf=gf)
    s1 = 0.9
    s2 = 4.
    eps = 1e-15
    
    for i in 1:100
        # do two half-tau RK4 steps
        X_s = rk4(rk4(X, τ/2, cf), τ/2, cf)

        # do one full-tau RK4 step
        X_b = rk4(X, τ, cf)
        
        # get diff
        scale = Δi/2 * (abs(X_s) + abs(X_b))
        scale = vcat(scale.r, scale.v)
        diff = vcat(X_s.r - X_b.r, X_s.v - X_b.v)        
        error_ratio = maximum(@. abs(diff)/(scale + eps))
        
        τ_old = τ
        τ = s1*τ_old*error_ratio^(-0.2)
        τ = max(τ, τ_old/s2)
        τ = min(τ, s2*τ_old)
        
        if error_ratio < 1.
            return (X_s, τ)
        end
    end
    
    error("Error, adaptive Runge-Kutta failed")
end

#——————————————————————————————————————————————————————————————————————————————————————————————————#
# Plotting

function adaptive_orbit(s::State, τ=0.005, Δi=0.003; cf=gf, orb=2)
    # Initial total energy, against which we compare total energy through trajectory.
    E_i = energy(s).E
    
    r = nrm(s.r)
    v = nrm(s.v)
    
    # Find period of orbit.
    a = r*μ/(2*μ - r*v^2)
    T = 2*pi*sqrt(a^3/μ)    
    time = 0.
    
    # Preallocate arrays for orbit data.
    ρ = []
    θ = []
    K = []
    U = []
    E = []
    Ers = []
    t = []
    
    # Do orb number of orbits.
    while time < orb*T
        time += τ
        nrg = energy(s)
        err = abs_error(nrg.E, E_i)
        
        push!(ρ, nrm(s.r))
        push!(θ, atan(s.r[2], s.r[1]))
        push!(K, nrg.K)
        push!(U, nrg.U)
        push!(E, nrg.E)
        push!(Ers, err)
        push!(t, time)
        rk_res = rka(s, τ, Δi, cf)
        s = rk_res[1]
        τ = rk_res[2]
    end
    
    return Orbit(ρ, θ, K, U, E, Ers, t)
end

function get_orbit(s::State, method=rk4, τ=0.005; cf=gf, orb=2.)   
    # Initial total energy, against which we compare total energy through trajectory.
    E_i = energy(s).E
    
    r = nrm(s.r)
    v = nrm(s.v)
    
    # Find period of orbit.
    a = r*μ/(2*μ - r*v^2)
    T = 2*pi*sqrt(a^3/μ)
    time = 0.

    # Divide period of orbit by tau to get evenly spaced steps; do orb number of orbits.
    steps = convert(Int64, orb*div(T, τ))
    
    # Preallocate arrays for orbit data.
    ρ = Array{Float64}(undef, steps)
    θ = Array{Float64}(undef, steps)
    K = Array{Float64}(undef, steps)
    U = Array{Float64}(undef, steps)
    E = Array{Float64}(undef, steps)
    Ers = Array{Float64}(undef, steps)
    t = Array{Float64}(undef, steps)
    
    # Begin trajectory calcs.
    for i in 1:steps
        time += τ
        nrg = energy(s)
        err = abs_error(nrg.E, E_i)
        ρ[i], θ[i], K[i], U[i], E[i], Ers[i], t[i] = nrm(s.r), atan(s.r[2], s.r[1]), nrg.K, nrg.U, nrg.E, err, time  
        s = method(s, τ)
    end
    
    return Orbit(ρ, θ, K, U, E, Ers, t)
end

function get_best_tau(method, τ; cf=gf)
    vels = []
    taus = []
    
    function propagate(s::State, v, τ, cf)
        # Find period of orbit.
        r = 35.
        a = r*μ/(2*μ - r*v^2)
        T = 2*pi*sqrt(a^3/μ)
        
        # Divide period of orbit by tau to get evenly spaced steps; do one full orbit.
        steps = div(T, τ)
        
        # Save initial state if we get a bad trajectory (err > 1%).
        S_i = s

        # Initial total energy, against which we compare total energy through trajectory.
        E_i = energy(s).E
        Ers = []
        
        # Initialize the error term.
        err = 0.
        
        # Flag used to recall if bad trajectory is found.
        bad_traj = false
        
        # Begin trajectory calcs.
        for i in 1:steps
            err = abs_error(energy(s).E, E_i)

            # Short-circuit if we get a bad trajectory.
            if err > 1.
                bad_traj = true
                break
            else
                push!(Ers, err)
                s = method(s, τ)
            end
        end

        if bad_traj
            # Start over with 99.95% of original tau if bad trajectory found.
            propagate(S_i, v, τ*0.995)
        else
            return τ
        end
    end
        
    for v in 0.5:0.01:1.0
        s = State([35.,0.], [0., v])
        
        # Push a good tau into the array.
        push!(taus, propagate(s, v, τ, cf))
        push!(vels, v)
    end
    
    return Taus_Vs(vels, taus)
end

function plot_orbit(o::Orbit)
    max_err = maximum(o.Ers)
    @printf("Maximum error: %.4f %%\n", max_err)
    
    p1 = plot([0], [0], label=:"", proj=:polar)
    p2 = plot(o.t, o.K, label=:"kinetic", 
                yaxis=:"Energy (M AU^2/yr^2)", xaxis=:"Time (yr)",
                legend=:bottomright)

    plot!(p1, o.θ, o.ρ, label=:"")
    plot!(p2, o.t, o.U, label=:"potential")
    plot!(p2, o.t, o.E, label=:"total")
    plot(p1, p2, size=(950,475))
end

function scatter_orbit(o::Orbit)
    max_err = maximum(o.Ers)
    @printf("Maximum error: %.4f %%\n", max_err)
    
    p1 = plot([0], [0], label=:"", proj=:polar)
    p2 = plot(o.t, o.K, label=:"kinetic", 
                yaxis=:"Energy (M AU^2/yr^2)", xaxis=:"Time (yr)",
                legend=:bottomright)

    scatter!(p1, o.θ, o.ρ, m=:+, label=:"")
    plot!(p2, o.t, o.U, label=:"potential")
    plot!(p2, o.t, o.E, label=:"total")
    plot(p1, p2, size=(950,475))
end

function plot_taus(T1::Taus_Vs, T2::Taus_Vs)
    hspace = [v for v in T1.v]
    p = scatter(size=(900,900),
            xaxis=(:log, "Initial Velocity (AU/yr)"),
            yaxis=(:log, "Timestep (yr)"),
            legend=:bottomright)
    plot!(p, hspace, T1.τ, label=:"")
    plot!(p, hspace, T2.τ, label=:"")
    scatter!(p, hspace, T1.τ, label=:"cromer")
    scatter!(p, hspace, T2.τ, label=:"rk4")
    plot(p)
end;

plot_orbit(get_orbit(State([1, 0], [0, 2*pi]), euler, 0.02, orb=4.))
plot_orbit(get_orbit(State([1, 0], [0, 2*pi]), cromer, 0.02, orb=4.))
plot_orbit(get_orbit(State([1, 0], [0, pi]), cromer, 0.005))

# Start Cromer with a timestep of 3 years.
cr_taus = get_best_tau(cromer, 3.)

# Start RK with a timestep of 15 years.
rk_taus = get_best_tau(rk4, 15.)

plot_taus(cr_taus, rk_taus)

hc = State([35, 0], [0, 0.19])
plot_orbit(get_orbit(hc, rk4, 0.01585, orb=4.)) # Plot Halley's Comet over 4 orbits.

s = State([1,0], [0, pi/2])
scatter_orbit(adaptive_orbit(s, 0.1, orb=1.5))