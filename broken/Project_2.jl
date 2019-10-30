using Plots
using Printf
gr()

#--------------------------------------------------------------------------------------------------------------------------#
# Constants.

const g = 9.80665 # Sea-level gravitational acceleration.
const m_bb = 0.145 # Baseball mass given by book.
const A_xs = pi*(0.037)^2 # Cross-sectional area of baseball, diameter 7.4 cm (r = 0.037 m)
const rho = 1.2041 # Dry air density, 20°C at sea-level pressure (wikipedia).
const cd = 0.35 # Baseball drag coefficient given by book.
const drag_mult = -(A_xs*cd*rho)/(2*m_bb) # Air drag multiplier.

#--------------------------------------------------------------------------------------------------------------------------#
# Helper functions.

# Returns an interpolated polynomial lambda.
function intrp(xx, yy)
    k = length(xx)
    if k != length(Set(xx))
        println("Cannot interpolate from duplicate x-values!")
    end

    function lbp(j, x)
        p = 1
        for m in 1:(k-1)
            p *= (x - xx[mod((j+m-1), k) + 1])/(xx[j] - xx[mod((j-m-1), k) + 1])
        end
        return p
    end

    return x -> sum(yy[j]*lbp(j, x) for j in 1:k)
end

# Returns an acceleration vector, given velocity and a drag boolean (defaults to false).
function acc(r, v, db = false)
    if db
        return drag_mult*sqrt(sum(v.^2))*v + [0, 0, -g]
    else
        return [0, 0, -g]
    end
end

# Analytic equations of motion.
function analytic_traj(v, t, h=0, inc=45)
    x = v*sind(inc)*t
    z = h + v*cosd(inc)*t - 1/2*g*t^2
    return [x, 0, z]
end

# Analytic solution to airless projectile range.
function airless_range(v, inc=45)
    return 2*v^2*cosd(inc)*sind(inc)/g
end

# Analytic soltion to airless projectile time of flight.
function tof(v, h=0, inc=45)
    return (v*cosd(inc) + sqrt((v*cosd(inc))^2 + 2*g*h))/g
end

# ¯\_(ツ)_/¯
function yeswind(r, ws = 0)
    return [ws, 0, 0]
end

# Interpolates between the last three points on a trajectory and 
# returns the x-value where the height above the ground is zero.
function getrange(traj)
    xx = traj[1][(end-2):end]
    yy = traj[3][(end-2):end]
    i = intrp(yy, xx)
    return i(0)
end

function find_onepercent(errs)
    findlast(x -> abs(x) < 1.0, errs)
end

function theoretical_error(exp, theo)
    return 100*(exp - theo)/theo
end

#--------------------------------------------------------------------------------------------------------------------------#
# Projectile motion methods.

# Euler
function euler(rn, vn, tau, ws, wind, db = false)
    a = acc(rn, vn - wind(rn, ws), db)
    v = vn + a*tau
    r = rn + vn*tau
    return (r,v)
end

# Euler-Cromer
function cromer(rn, vn, tau, ws, wind, db = false)
    a = acc(rn, vn - wind(rn, ws), db)
    v = vn + a*tau
    r = rn + v*tau
    return (r,v)
end

# Midpoint
function midpoint(rn, vn, tau, ws, wind, db = false)
    a = acc(rn, vn - wind(rn, ws), db)
    v = vn + a*tau
    r = rn + vn*tau + 1/2*a*tau^2
    return (r,v)
end

# Returns an array of arrays of position vector components along the trajectory corresponding
# to the given initial conditions, computed using the desired method, where:
# r and v are vectors in R3
# db is the drag boolean (true for drag, defaults as false for airless)
# wind is a function of position, 
# f is the method to use to compute the trajectory (euler, cromer, or midpoint)   
function trajectory(r, v, f, tau; ws = 0, wind = yeswind, db = false)
    # Initial position
    pos = [[r[1]], [r[2]], [r[3]]]
    
    while true
        # Have we struck the ground?
        if r[3] < 0
            break
        else
            r,v = f(r, v, tau, ws, wind, db)
            push!(pos[1], r[1])
            push!(pos[2], r[2])
            push!(pos[3], r[3])
        end
    end

    return pos # Output is of type: [[x_vals], [y_vals], [z_vals]], for plotting purposes.
end

#--------------------------------------------------------------------------------------------------------------------------#
# Plotting functions.

function plot_euler()
    vi = 15. # Initial speed.
    inc = 45. # Inclination angle (angle above xy plane)
    tau = 0.1 # Timestep

    # Initial conditions (taking z to be vertical)
    r1 = [0., 0., 0.]
    v1 = [vi*sind(45), 0, vi*cosd(45)]

    tf = tof(vi)
    tt = [t for t in 0.0:0.001:tf]
    tx = [analytic_traj(vi, t)[1] for t in tt]
    tz = [analytic_traj(vi, t)[3] for t in tt]

    euler_traj = trajectory(r1, v1, euler, tau)

    th_range = airless_range(vi)
    eu_range = getrange(euler_traj)
    eu_error = theoretical_error(eu_range, th_range)

    @printf("theory x-range:   %-6.2f m\n", th_range)
    @printf("euler x-range:    %-6.2f m. Error: %-2.2f %%\n", eu_range, eu_error)

    plot(tx, tz, w=4, label=:"theory", size=(900,300), aspect_ratio=1)
    scatter!(euler_traj[1], euler_traj[3], label=:"euler")
end

function plot_all_three()
    vi = 15. # Initial speed.
    inc = 45. # Inclination angle (angle above xy plane)
    tau = 0.1 # Timestep

    # Initial conditions (taking z to be vertical)
    r1 = [0., 0., 0.]
    v1 = [vi*sind(45), 0, vi*cosd(45)]

    tf = tof(vi)
    tt = [t for t in 0.0:0.001:tf]
    tx = [analytic_traj(vi, t)[1] for t in tt]
    tz = [analytic_traj(vi, t)[3] for t in tt]

    euler_traj = trajectory(r1, v1, euler, tau)
    cromer_traj = trajectory(r1, v1, cromer, tau)
    mp_traj = trajectory(r1, v1, midpoint, tau)

    th_range = airless_range(vi)
    eu_range = getrange(euler_traj)
    cr_range = getrange(cromer_traj)
    mp_range = getrange(mp_traj)

    eu_error = theoretical_error(eu_range, th_range)
    cr_error = theoretical_error(cr_range, th_range)
    mp_error = theoretical_error(mp_range, th_range)

    @printf("theory x-range:   %-6.2f m\n", th_range)
    @printf("euler x-range:    %-6.2f m. Error: %-2.2f %%\n", eu_range, eu_error)
    @printf("cromer x-range:   %-6.2f m. Error: %-2.2f %%\n", cr_range, cr_error)
    @printf("midpoint x-range: %-6.2f m. Error: %-2.2f %%\n", mp_range, mp_error)

    plot(tx, tz, w=4, label=:"theory", size=(900,300), aspect_ratio=1)
    plot!(euler_traj[1], euler_traj[3], w=4, label=:"euler")
    plot!(cromer_traj[1], cromer_traj[3], w=4, label=:"cromer")
    plot!(mp_traj[1], mp_traj[3], w=4, label=:"midpoint")
end

function vi_vs_err()
    vels = [v for v in 10.0:10.0:500.0]
    eu_errs = []
    ec_errs = []
    mp_errs = []
    
    r = [0., 0., 0.]
    tau = 0.1
    
    for vi in vels
        v = [vi*sind(45), 0, vi*cosd(45)]
        th_range = airless_range(vi)
        
        eu = trajectory(r, v, euler, tau)
        ec = trajectory(r, v, cromer, tau)
        mp = trajectory(r, v, midpoint, tau)

        eu_range = getrange(eu) 
        ec_range = getrange(ec)
        mp_range = getrange(mp)

        eu_err = theoretical_error(eu_range, th_range)
        ec_err = theoretical_error(ec_range, th_range)
        mp_err = theoretical_error(mp_range, th_range)

        push!(eu_errs, eu_err)
        push!(ec_errs, ec_err)
        push!(mp_errs, mp_err)
    end

    plot(title=:"Range error vs initial velocity", xlabel=:"Initial velocity (m/s)", ylabel=:"Error (%)", size=(900,500))
    plot!(vels, eu_errs, w=4, label=:"euler")
    plot!(vels, ec_errs, w=4, label=:"cromer")
    plot!(vels, mp_errs, w=4, label=:"midpoint")
end

# Euler and Euler-Cromer only.
function tau_vs_err()
    timestep = [t for t in 0.01:0.001:0.08]
    eu_errs = []
    ec_errs = []

    r = [0., 0., 0.]
    v = [50*sind(45), 0, 50*cosd(45)]
    th_range = airless_range(50)
    
    for tau in timestep
        eu = trajectory(r, v, euler, tau)
        ec = trajectory(r, v, cromer, tau)

        eu_range = getrange(eu) 
        ec_range = getrange(ec)

        eu_err = theoretical_error(eu_range, th_range)
        ec_err = theoretical_error(ec_range, th_range)

        push!(eu_errs, eu_err)
        push!(ec_errs, ec_err)
    end
    
    last_eu = find_onepercent(eu_errs)
    last_ec = find_onepercent(ec_errs)
    
    println("Euler largest tau:  $(timestep[last_eu]) s")
    println("Cromer largest tau: $(timestep[last_ec]) s")
    
    plot(title=:"Range error vs time step", legend=:bottomleft, xlabel=:"Timestep (s)", ylabel=:"Error (%)", size=(900,600))
    plot!(timestep, eu_errs, w=4, label=:"euler")
    plot!(timestep, ec_errs, w=4, label=:"cromer")
    scatter!([timestep[last_eu]], [eu_errs[last_eu]], markercolor = :black, markersize=8, label=:"last tau within 1%")
    scatter!([timestep[last_ec]], [ec_errs[last_ec]], markercolor = :black, markersize=8, label=:"last tau within 1%")
end

# Midpoint only.
function tau_vs_mp_err()
    timestep = [t for t in 0.01:0.001:0.5]
    mp_errs = []

    r = [0., 0., 0.]
    v = [50*sind(45), 0, 50*cosd(45)]
    th_range = airless_range(50)
    
    for tau in timestep
        mp = trajectory(r, v, midpoint, tau)
        mp_range = getrange(mp)
        mp_err = theoretical_error(mp_range, th_range)
        push!(mp_errs, mp_err)
    end

    plot(title=:"Midpoint range error vs time step", xlabel=:"Timestep (s)", ylabel=:"Error (%)", size=(900,500))
    plot!(timestep, mp_errs, w=3, label=:"")
end

function range_vs_angle()
    angs = [a for a in 10.0:0.1:50.0]
    
    r = [0., 0., 1.]
    tau = 0.1

    rnges = []
    for ang in angs
        v = [50*sind(90.0-ang), 0, 50*cosd(90.0-ang)]
        traj = trajectory(r, v, midpoint, tau, db=true)
        xrange = getrange(traj)
        push!(rnges, xrange)
    end

    ind = argmax(rnges)
    maxr = rnges[ind]
    opta = angs[ind]
    println("Maximum range: $maxr")
    println("Optimum angle: $opta")

    plot(size=(900,900), title=:"Range vs launch angle", xlabel=:"Angle (°)", ylabel=:"Range (m)")
    scatter!(angs, rnges, markershape = :circle,
                           markersize = 3,
                           markeralpha = 0.6,
                           markercolor = :red,
                           markerstrokewidth = 0,
                           markerstrokealpha = 0.,
                           markerstrokecolor = :red, label=:"")
end

function range_vs_wind()
    wind = [w for w in -300.0:10.0:300.0]
    
    r = [0., 0., 1.]
    v = [100*sind(55), 0, 100*cosd(55)]
    tau = 0.1

    eu_ranges = []
    ec_ranges = []
    mp_ranges = []
    for w in wind
        eu = trajectory(r, v, euler, tau, ws=w, db=true)
        ec = trajectory(r, v, cromer, tau, ws=w, db=true)
        mp = trajectory(r, v, midpoint, tau, ws=w, db=true)

        eu_range = getrange(eu) 
        ec_range = getrange(ec)
        mp_range = getrange(mp)
        
        push!(eu_ranges, eu_range)
        push!(ec_ranges, ec_range)
        push!(mp_ranges, mp_range)
    end

    plot(size=(900,900), title=:"Range vs windspeed", xlabel=:"Windspeed (m/s)", ylabel=:"Range (m)")
    plot!(wind, eu_ranges, w=2, label=:"euler")
    plot!(wind, ec_ranges, w=2, label=:"cromer")
    plot!(wind, mp_ranges, w=2, label=:"midpoint")
end

function windy_day()
    wind = [w for w in -60.0:5.0:-20.0]
    
    r = [0., 0., 0.]
    v = [50*sind(45), 0, 50*cosd(45)]
    tau = 0.01 # Timestep
    
    p = plot(size=(900,300), title=:"Windy Day Trajectories", xlabel=:"Range (m/s)", ylabel=:"Height (m)", aspect_ratio=1)
    
    for w in wind
        mp = trajectory(r, v, midpoint, tau, ws=w, db=true)
        plot!(p, mp[1], mp[3], w=2, label=:"")
    end
    
    p
end

#--------------------------------------------------------------------------------------------------------------------------#
# Fancy vector fields.

function helix(r, ws = 0)
    x,y,z = r # unpack the position vector
    return 200*[sin(y)^3, -cos(x)^2, -1/2 * exp(-(x/10)^2 -(y/10)^2)]
end

function sinks(r, ws = 0)
    x,y,z = r # unpack the position vector
    return 50*[2*(sin(y/2+z)), 2*(sin(y) + sin(z)), sin(y+x)]
end

function hole(r, ws = 0)
    x,y,z = r # unpack the position vector
    return [10 + 3*y - x^3, 2*z^3 - y^3, x^3 - z^3]
end

#--------------------------------------------------------------------------------------------------------------------------#
# Fancy plot of fancy vector fields.

function neato_plot()
    azi = 52. # Azimuthal angle (angle in xy plane)
    inc = 75. # Inclination angle (angle above xy plane)
    r = [0., -10., 0]
    v = [15*sind(90-inc)*cosd(azi), 15*sind(90-inc)*sind(azi), 15*cosd(90-inc)]
    tau = 0.0001 # Timestep

    hlx = trajectory(r, v, midpoint, tau, wind=helix, db=true)
    snk = trajectory(r, v, midpoint, tau, wind=sinks, db=true)
    hol = trajectory(r, v, midpoint, tau, wind=hole, db=true)

    p = plot(hlx[1], hlx[2], hlx[3], w=5, label=:"", size=(900,900))
    plot!(p, snk[1], snk[2], snk[3], w=5, label=:"")
    plot!(p, hol[1], hol[2], hol[3], w=5, label=:"")

    azi = 30
    inc = 40

    v = [15*sind(90-inc)*cosd(azi), 15*sind(90-inc)*sind(azi), 15*cosd(90-inc)]

    hlx = trajectory(r, v, midpoint, tau, wind=helix, db=true)
    snk = trajectory(r, v, midpoint, tau, wind=sinks, db=true)
    hol = trajectory(r, v, midpoint, tau, wind=hole, db=true)

    plot!(p, hlx[1], hlx[2], hlx[3], w=5, label=:"")
    plot!(p, snk[1], snk[2], snk[3], w=5, label=:"")
    plot!(p, hol[1], hol[2], hol[3], w=5, label=:"")

    # savefig(p, "test.html")
    p
end