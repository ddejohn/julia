# using LinearAlgebra


function cspl(; data)
    xx, yy = map(y -> map(x -> x[y], data), 1:2)
    
    nsp = length(xx) - 1
    inter = nsp - 1
    n = 2*nsp
    p = [4*(j-1) .+ [i for i in 1:4] for j in 1:nsp]
    
    m = zeros(2*n, 2*n)
    
    zz = fill(0.0, 2*n)
    zz[1] = yy[1]
    zz[n] = yy[end]
    for i in 2:inter+1
        j = 2*(i-2) .+ [1,2]

        print(j)
        zz[j[1]+1] = yy[i]
        zz[j[2]+1] = yy[i]
    end

    print(zz)
    
    # s, s', s'' coefficients
    coefs = [[1,1,1,1], [0,1,2,3], [0,0,2,6]]
    exps = [[0,1,2,3], [0,0,1,2], [0,0,0,1]]

    s = [x -> coefs[1][j]*(x-xi)^exps[1][j] for xi in xx[1:end-1], j in 1:4]
    ds = [x -> coefs[2][j]*(x-xi)^exps[2][j] for xi in xx[1:end-1], j in 1:4]
    dds = [x -> coefs[3][j]*(x-xi)^exps[3][j] for xi in xx[1:end-1], j in 1:4]
    
    s = vcat(s, vcat(ds, dds))

    function seval(s, x)
        return [round(f(x), digits=5) for f in s]
    end

    # Loop 1
    j=1
    for i in 1:2:n
        j = mod1(j, nsp)
        m[i, p[j]] = seval(s[j,:], xx[j])
        m[i+1, p[j]] = seval(s[j,:], xx[j+1])
        j += 1
    end

    # Loop 2
    for i in 1:inter
        q = mod1(i, nsp)

        m[i+n, p[q]] = seval(s[i+nsp,:], xx[i+1])
        m[i+n, p[q+1]] = -seval(s[i+nsp+1,:], xx[i+1])

        m[i+n+inter, p[q]] = seval(s[i+n,:], xx[i+1])
        m[i+n+inter, p[q+1]] = -seval(s[i+n+1,:], xx[i+1])
    end

    # Free boundary conditions
    m[end-1, p[1]] = seval(s[1+n,:], xx[1]) 
    m[end, p[end]] = seval(s[3*nsp,:], xx[end]) 

    # println()
    # for r in 1:size(m,1)
    #     println(m[r,:])
    # end

    # Solve for the coefficients
    z = m\zz

    z =  [[z[4*(j-1)+i] for i in 1:4] for j in 1:nsp]
    for r in z
        println(r)
    end
end


# data = [(-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101)]
data = [(0.1, log(0.1)), (1, log(1)), (2, log(2)), (2.9, log(2.9))]
# data = [(1,1), (2,1), (3,2), (4,6), (5,24)]

cspl(data=data)