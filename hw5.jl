# Cubic spline solver, defaults to free boundary conditions
function cspl(; data, bc="free")
    xx, yy = map(y -> map(x -> x[y], data), 1:2)
    
    # number of splines, number of interior points, etc
    nsp = length(xx) - 1
    inter = nsp - 1
    n = 2*nsp

    # Initialize a matrix
    m = zeros(2*n, 2*n)

    # Partition the columns into sets of four
    p = [4*(j-1) .+ [i for i in 1:4] for j in 1:nsp]
    
    # Construct the solution vector
    zz = fill(0.0, 2*n)
    zz[1] = yy[1]
    zz[n] = yy[end]
    for i in 2:inter+1
        j = 2*(i-2) .+ [1,2]

        zz[j[1]+1] = yy[i]
        zz[j[2]+1] = yy[i]
    end
    
    # s, s', s'' coefficients
    coefs = [[1,1,1,1], [0,1,2,3], [0,0,2,6]]
    exps = [[0,1,2,3], [0,0,1,2], [0,0,0,1]]

    # s, s', s'' general forms as lambda expressions
    s = [x -> coefs[1][j]*(x-xi)^exps[1][j] for xi in xx[1:end-1], j in 1:4]
    ds = [x -> coefs[2][j]*(x-xi)^exps[2][j] for xi in xx[1:end-1], j in 1:4]
    dds = [x -> coefs[3][j]*(x-xi)^exps[3][j] for xi in xx[1:end-1], j in 1:4]
    
    # store s, s', s'' in a single matrix
    s = vcat(s, vcat(ds, dds))

    # Evaluate spline lambdas at a given x value
    function seval(s, x)
        return [round(f(x), digits=5) for f in s]
    end

    # this loop fills in the first 2*(number of splines) rows
    # where each pair of rows corresponds to a single spline
    # evaluated at its two end points
    j=1
    for i in 1:2:n
        j = mod1(j, nsp)
        m[i, p[j]] = seval(s[j,:], xx[j])
        m[i+1, p[j]] = seval(s[j,:], xx[j+1])
        j += 1
    end

    # this loop fills in the next 2*(number of interior points) rows
    # where each row consists of the first derivatives of adjacent splines
    # evaluated at the boundary between the two splines
    # and each row after the all the first derivatives then evaluates
    # all of the second derivatives
    for i in 1:inter
        q = mod1(i, nsp)

        # s' is evaluated at each spline transition
        # the second 'seval' is negative which corresponds to
        # subtracting the (i+1)th spline on the (i)th spline's side
        m[i+n, p[q]] = seval(s[i+nsp,:], xx[i+1])
        m[i+n, p[q+1]] = -seval(s[i+nsp+1,:], xx[i+1])

        # s'' is evaluated at each spline transition
        m[i+n+inter, p[q]] = seval(s[i+n,:], xx[i+1])
        m[i+n+inter, p[q+1]] = -seval(s[i+n+1,:], xx[i+1])
    end

    if bc == "free"
        # Free boundary conditions
        m[end-1, p[1]] = seval(s[1+n,:], xx[1]) 
        m[end, p[end]] = seval(s[3*nsp,:], xx[end])
    elseif bc == "nk"
        # Not-a-knot boundary conditions
        m[end-1,:] = vcat([0, 0, 0, 6, 0, 0, 0, -6], fill(0.0, 2*n-8))
        m[end,:] = vcat(fill(0.0, 2*n-8), [0, 0, 0, 6, 0, 0, 0, -6])
    else
        println("please specify valid boundary conditions!")
        return 0,0
    end

    # Solve for the coefficients
    z = m\zz
    
    # Partition coefficients into sets of four
    z =  [[z[4*(j-1)+i] for i in 1:4] for j in 1:nsp]
    soln = []

    # Construct the splines with the given coefficients
    # then tore the splines in a vector to return
    for i in 1:length(z)
        push!(soln, x -> sum(z[i][j]*(x-xx[i])^(j-1) for j in 1:4))
    end

    return soln
end


data = [(-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101)]
# data = [(0.1, log(0.1)), (1, log(1)), (2, log(2)), (2.9, log(2.9))]
# data = [(1,1), (2,1), (3,2), (4,6), (5,24)]

splines = cspl(data=data)

# s1, s2 = splines
# x1, x2 = collect(-0.5:0.01:-0.25), collect(-0.25:0.01:0)

# p = plot([],[],label="",size=(900,900))
# plot!(p, x1, s1.(x1))
# plot!(p, x2, s2.(x2))