using Random
#const rng = MersenneTwister(1) # use if multiple sequences are needed
Random.seed!(1)

function mc_unit_circle(n)
    k = sum(1 for _ in 1:n if rand()^2 + rand()^2 <= 1)
    return 4*k/n
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

# Iteration vals, in decadal-magnitudal steps.
nn = [k*10^n for n in 1:10 for k in 1:9]

# Ellipsoid axes.
a = 5/sqrt(2)
b = 5/sqrt(3)
c = 5

# Truncated ellipsoid bounding box.
xx = [-1, a]
yy = [-b, b]
zz = [-2, 2]

whole_analytic_vol = 4/3*pi*a*b*c
# truncated_analytic_vol = 

whole_vols = [a*b*c*mc_unit_sphere(n) for n in nn]
truncated_vols = [trunc_ellipsoid(xx,yy,zz,n) for n in nn]