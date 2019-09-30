# using LinearAlgebra


function cspl(; data)
    xx, yy = map(y -> map(x -> x[y], data), 1:2)
    nsp = length(xx) - 1
    inter = nsp - 1
    n = 2*nsp
    k = 4*nsp
    p = [4*(j-1) .+ [i for i in 1:4] for j in 1:nsp]

    m = zeros(k, k)

    # s, s', s'' coefficients
    coefs = [[1,1,1,1], [0,1,2,3], [0,0,2,6]]
    exps = [[0,1,2,3], [0,0,1,2], [0,0,0,1]]

    s = [x -> coefs[1][j]*(x-xi)^exps[1][j] for xi in xx[1:end-1], j in 1:4]
    ds = [x -> coefs[2][j]*(x-xi)^exps[2][j] for xi in xx[1:end-1], j in 1:4]
    dds = [x -> coefs[3][j]*(x-xi)^exps[3][j] for xi in xx[1:end-1], j in 1:4]
    
    s = vcat(s, vcat(ds, dds))

    # # Loop 1
    # for i in 2:n-2
    #     q = mod1(i-1,nsp)
    #     si = s[i-1,:]
    #     dsi = s[i+nsp-1,:]
    #     ddsi = s[i+2*nsp-1,:]

    #     println("i: $(i)\ni+n: $(i+n)\ni+n+inter: $(i+n+inter)")
    #     # println("p[q]: $(p[q])")

    #     m[i, p[q]] = [round(f(xx[i]), digits=5) for f in si]
    #     m[i+1, p[q]] = [round(f(xx[i]), digits=5) for f in si]

    #     m[i+n, p[q]] = [round(f(xx[i]), digits=5) for f in dsi]
    #     m[i+n+inter, p[q]] = [round(f(xx[i]), digits=5) for f in ddsi]
    # end

    # m[1, p[1]] = [round(f(xx[1]), digits=5) for f in s[1,:]]
    # m[end-4, p[end]] = [round(f(xx[end]), digits=5) for f in s[end,:]]

    # Loop 1
    j=1
    for i in 1:2:n
        j = mod1(j, nsp)
        m[i, p[j]] = [round(f(xx[j]), digits=5) for f in s[j,:]]
        m[i+1, p[j]] = [round(f(xx[j+1]), digits=5) for f in s[j,:]]
        j += 1
    end

    # Loop 2
    for i in 1:inter
        q = mod1(i, nsp)

        m[i+n, p[q]] = [round(f(xx[i+1]), digits=5) for f in s[i+nsp,:]]
        m[i+n, p[q+1]] = [round(-f(xx[i+1]), digits=5) for f in s[i+nsp+1,:]]

        m[i+n+inter, p[q]] = [round(f(xx[i+1]), digits=5) for f in s[i+2*nsp,:]]
        m[i+n+inter, p[q+1]] = [round(-f(xx[i+1]), digits=5) for f in s[i+2*nsp+1,:]]
    end

    # Boundary conditions
    m[end-1, p[1]] = [round(f(xx[1]), digits=5) for f in s[1+2*nsp,:]]
    m[end, p[end]] = [round(f(xx[end]), digits=5) for f in s[3*nsp,:]] 

    println()
    for r in 1:size(m,1)
        println(m[r,:])
    end
end


data = [(-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101)]
# data = [(0.1, log(0.1)), (1, log(1)), (2, log(2)), (2.9, log(2.9))]
# data = [(1,1), (2,1), (3,2), (4,6), (5,24)]

cspl(data=data)