using Plots
using Random
gr()

# MCI iteration values, in decadal-magnitudinal steps.
nn = [k*10^n for n in 3:9 for k in 1:9];

function mc_unit_circle(n)
    k = sum(1 for _ in 1:n if rand()^2 + rand()^2 < 1)
    return 4*k/n
end

function min_thresh(vals, upper, lower)
    last_lower = findlast(x -> x < lower, vals)
    last_upper = findlast(x -> x > upper, vals)
    return max(last_lower, last_upper) + 1
end

function mc_unit_sphere(n)
    k = sum(1 for _ in 1:n if rand()^2 + rand()^2 + rand()^2 < 1)
    return 8*k/n
end

function trunc_ellipsoid(xx,yy,zz,n)
    k = 0
    dx = xx[2] - xx[1]
    dy = yy[2] - yy[1]
    dz = zz[2] - zz[1]

    for _ in 1:n
        x = xx[1] + dx*rand()
        y = yy[1] + dy*rand()
        z = zz[1] + dz*rand()

        # By definition, these points are already within the truncation boundaries,
        # so we simply need to check if they also lie within the ellipsoid.
        if 2x^2 + 3y^2 + z^2 <= 25
            k += 1
        end
    end

    return dx*dy*dz*k/n
end

function rho_ellipsoid(xx,yy,zz,n)
    k = 0
    dx = xx[2] - xx[1]
    dy = yy[2] - yy[1]
    dz = zz[2] - zz[1]

    for _ in 1:n
        x = xx[1] + dx*rand()
        y = yy[1] + dy*rand()
        z = zz[1] + dz*rand()

        if 2x^2 + 3y^2 + z^2 <= 25
            k += x^2
        end
    end

    return dx*dy*dz*k/n
end

function theoretical_error(exp, theo)
    return 100*(exp - theo)/theo
end

function unbiased_error(x1, x2)
    return 200*abs(x1 - x2)/(x1 + x2)
end

# Unseeded.
yu = [mc_unit_circle(n) for n in nn]

# Seeded.
Random.seed!(1)
ys = [mc_unit_circle(n) for n in nn]

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Plot 1

# Base plot x-range.
xx = [minimum(nn), maximum(nn)]

# Upper and lower bounds for values that will round to 3.14.
upper = 3.145
lower = 3.135

# Minimum iterations to stay within threshold.
last_2digits = min_thresh(yu, upper, lower)

# Base plot.
plot(xx, [pi,pi], l=(1,:dash,:green), label=:pi, size=(900,500), xscale=:log10)

# Pi threshold.
plot!(xx, [upper, upper], l=(1,:dash,:red), label=:"3.145")
plot!(xx, [lower, lower], l=(1,:dash,:blue), label=:"3.135")

# MCI.
plot!(nn, yu, m=(:d, :grey), l=(1,:grey), label=:"unseeded mci")

# Minimum iterations.
scatter!([nn[last_2digits]], [yu[last_2digits]], m=(8, :d, :orange), label=:"minimum iterations to 3.14")

# Center the view on pi.
max_delta = max(pi - 1.001*maximum(yu), pi - 0.999*minimum(yu))
yaxis!((pi - max_delta, pi + max_delta))

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Plot 2

# Base plot x-range.
xx = [minimum(nn[last_2digits:end]), maximum(nn)]

# Upper and lower bounds for values that will round to 3.1416.
upper = 3.14165
lower = 3.14155

# Minimum iterations to stay within threshold.
last_3digits = min_thresh(yu, upper, lower)

# Base plot.
plot(xx, [pi,pi], l=(1,:dash,:green), label=:pi, size=(900,500), xscale=:log10)

# Pi threshold.
plot!(xx, [upper, upper], l=(1,:dash,:red), label=:"3.14165")
plot!(xx, [lower, lower], l=(1,:dash,:blue), label=:"3.14155")

# MCI.
plot!(nn[last_2digits:end], yu[last_2digits:end], m=(:d, :grey), l=(1,:grey), label=:"unseeded mci")

# Minimum iterations.
scatter!([nn[last_2digits]], [yu[last_2digits]], m=(8, :d, :red), label=:"minimum iterations to 3.14")
scatter!([nn[last_3digits]], [yu[last_3digits]], m=(8, :d, :orange), label=:"minimum iterations to 3.1416")

# Center the view on pi.
max_delta = max(pi - 1.001*maximum(yu[last_2digits:end]), pi - 0.999*minimum(yu[last_2digits:end]))
yaxis!((pi - max_delta, pi + max_delta))

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Plot 3

# Base plot x-range.
xx = [minimum(nn), maximum(nn)]

# Upper and lower bounds for values that will round to 3.14.
upper = 3.145
lower = 3.135

# Minimum iterations to stay within threshold.
last_s = min_thresh(ys, upper, lower)

# Base plot.
plot(xx, [pi,pi], l=(1,:dash,:green), label=:pi, size=(900,500), xscale=:log10)

# Pi threshold.
plot!(xx, [upper, upper], l=(1,:dash,:red), label=:"3.145")
plot!(xx, [lower, lower], l=(1,:dash,:blue), label=:"3.135")

# MCI
plot!(nn, ys, m=(:d, :grey), l=(1,:grey), label=:"seeded mci")

# Minimum iterations.
scatter!([nn[last_s]], [ys[last_s]], m=(10, :d, :orange), label=:"minimum iterations")

# Center the view on pi.
max_delta = max(pi - 1.001*maximum(yu), pi - 0.999*minimum(yu))
yaxis!((pi - max_delta, pi + max_delta))

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Plot 4

# Base plot x-range.
xx = [minimum(nn[last_s:end]), maximum(nn)]

# Upper and lower bounds for values that will round to 3.1416.
upper = 3.14165
lower = 3.14155

# Minimum iterations to stay within threshold.
last_3digits = min_thresh(ys, upper, lower)

# Base plot.
plot(xx, [pi,pi], l=(1,:dash,:green), label=:pi, size=(900,500), xscale=:log10)

# Pi threshold.
plot!(xx, [upper, upper], l=(1,:dash,:red), label=:"3.14165")
plot!(xx, [lower, lower], l=(1,:dash,:blue), label=:"3.14155")

# MCI.
plot!(nn[last_s:end], ys[last_s:end], m=(:d, :grey), l=(1,:grey), label=:"seeded mci")

# Minimum iterations.
scatter!([nn[last_s]], [ys[last_s]], m=(8, :d, :red), label=:"minimum iterations to 3.14")
scatter!([nn[last_3digits]], [ys[last_3digits]], m=(8, :d, :orange), label=:"minimum iterations to 3.1416")

# Center the view on pi.
max_delta = max(pi - 1.001*maximum(ys[last_s:end]), pi - 0.999*minimum(ys[last_s:end]))
yaxis!((pi - max_delta, pi + max_delta))

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Plot 5

# Base plot x-range.
xx = [minimum(nn), maximum(nn)]

# Base plot.
plot(xx, [pi,pi], l=(1,:dash,:green), label=:pi, size=(900,500), xscale=:log10)

# MCIs.
plot!(nn, yu, m=(:d, :grey), l=(1,:grey), label=:"unseeded mci")
plot!(nn, ys, m=(:d, :orange), l=(1,:orange), label=:"seeded mci")

# Minimum iterations.
scatter!([nn[last_2digits]], [yu[last_2digits]], m=(8, :white), label=:"unseeded 3.14")
scatter!([nn[last_s]], [yu[last_s]], m=(8, :white), label=:"seeded 3.14")

# Center the view on pi.
max_delta = max(pi - 1.001*maximum(yu), pi - 0.999*minimum(yu))
yaxis!((pi - max_delta, pi + max_delta))

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Plot 6

# Base plot x-range.
xx = [10^8, 9*10^9]

# Base plot.
plot(xx, [pi,pi], l=(1,:dash,:green), label=:pi, size=(900,500), xscale=:log10)

# MCIs.
plot!(nn[(end-17):end], yu[(end-17):end], m=(:d, :grey), l=(1,:grey), label=:"unseeded mci")
plot!(nn[(end-17):end], ys[(end-17):end], m=(:d, :orange), l=(1,:orange), label=:"seeded mci")

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# Begin ellipsoid calcs

# Not enough RAM for this.
# t = @elapsed iters = [[rand(), rand(), rand()] for _ in 1:10^9]

# Ellipsoid axes.
a = 5/sqrt(2)
b = 5/sqrt(3)
c = 5

# Truncated ellipsoid bounding box.
xb = [-1, a]
yb = [-b, b]
zb = [-2, 2]

# Whole ellipsoid analytic mass and MCI result.
a_vol = 4/3*pi*a*b*c
t1 = @elapsed w_vol = a*b*c*mc_unit_sphere(10^11)
println("Analytic mass: ", a_vol)
println("MCI mass: ", w_vol)
println("Error %: ", theoretical_error(w_vol, a_vol))
println("time: ", t1, "\n")

# Truncated "analytic mass" and MCI result.
at_vol = 82.84888725086089
t2 = @elapsed t_vol = trunc_ellipsoid(xb,yb,zb,10^11)
println("Analytic truncated mass: ", at_vol)
println("Truncated ellipsoid mass: ", t_vol)
println("Error %: ", theoretical_error(at_vol, t_vol))
println("time: ", t2, "\n")

# Variable density "analytic mass" and MCI result.
ar_vol = 187.3430623682275
t3 = @elapsed r_vol = rho_ellipsoid(xb,yb,zb,10^11)
println("Analytic variable density mass: ", ar_vol)
println("Variable density ellipsoid mass: ", r_vol)
println("Error %: ", theoretical_error(ar_vol, r_vol))
println("time: ", t3)