# using LinearAlgebra




# # cubic spline constructor
# function csplines(; data)
#     xx, yy = map(y -> map(x -> x[y], data), 1:2)
#     nsp = length(xx) - 1

#     m = zeros(4*nsp, 4*nsp)

#     # s, s', s'' coefficients
#     coefs = [[1,1,1,1], [0,1,2,3], [0,0,2,6]]
#     exps = [[0,1,2,3], [0,0,1,2], [0,0,0,1]]
    
    
#     s = [x -> coefs[1][j]*(x-xi)^exps[1][j] for xi in xx[1:end-1], j in 1:4]
#     ds = [x -> coefs[2][j]*(x-xi)^exps[2][j] for xi in xx[1:end-1], j in 1:4]
#     dds = [x -> coefs[3][j]*(x-xi)^exps[3][j] for xi in xx[1:end-1], j in 1:4]
    
#     s = vcat(s, vcat(ds, dds))
#     # s = ["$(coefs[1][j])*(x-$(xi))^$(exps[1][j])" for xi in xx[1:end-1], j in 1:4]
#     # ds = ["$(coefs[2][j])*(x-$(xi))^$(exps[2][j])" for xi in xx[1:end-1], j in 1:4]
#     # dds = ["$(coefs[3][j])*(x-$(xi))^$(exps[3][j])" for xi in xx[1:end-1], j in 1:4]
    
#     j = collect(1:2:size(s,1)+1)
#     k = collect(1:4)
#     p = collect(1+2*nsp:2:size(m,1))

#     dsi = collect(2*nsp:2*nsp+nsp-1)

#     z = 2*nsp
#     # for r in 1:size(s,1)
#     #     println(join(s[r,:], " + "))
#     # end

#     # s(x) conditions
#     for i in 1:nsp
#         pl = 4*(mod1(i, nsp-1)-1) .+ k
#         pr = 4 .+ pl

#         dsi = mod1(i + 2*nsp, )
#         ddsi = mod1(i + 3*nsp, )






#         # println("$(p), $(r3)")
        
#         m[j[i], r1] = [round(f(xx[i]), digits=5) for f in s[i,:]]
#         m[j[i]+1, r1] = [round(f(xx[i+1]), digits=5) for f in s[i,:]]
    
#         m[j[i]+z, r2] = [round(f(xx[i+1]), digits=5) for f in s[i+nsp,:]]
#         m[j[i]+z+1, r2] = [round(f(xx[i+1]), digits=5) for f in s[i+z,:]]
        
#         m[j[i]+z, r3] = [round(-f(xx[i+1]), digits=5) for f in s[i+nsp+1,:]]
#         m[j[i]+z+1, r3] = [round(-f(xx[i+1]), digits=5) for f in s[i+z+1,:]]
#     end







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

    for i in 2:n-1
        q = mod1(i,nsp)
        si = s[i-1,:]
        dsi = s[i+nsp-1,:]
        ddsi = s[i+2*nsp-1,:]

        m[i, p[q]] = [round(f(xx[i]), digits=5) for f in si]
        m[i+n-1, p[q]] = [round(f(xx[i]), digits=5) for f in dsi]
        m[i+n+inter-1, p[q]] = [round(f(xx[i]), digits=5) for f in ddsi]

        # m[i+1, p[q]] = [round(f(xx[i+1]), digits=5) for f in si]
        # m[i+n, p[q]] = [round(f(xx[i+1]), digits=5) for f in dsi]
        # m[i+n+inter, p[q]] = [round(f(xx[i+1]), digits=5) for f in ddsi]
    end

    m[1, p[1]] = [round(f(xx[1]), digits=5) for f in s[1]]
    m[end-4, p[end]] = [round(f(xx[end]), digits=5) for f in s[end]]

    println()
    for r in 1:size(m,1)
        println(m[r,:])
    end
end







    # for i in 1:nsp+1
    #     q = i + nsp
    #     p = i + 2*nsp
    #     z = i+1

    #     println("$(vcat(r1,r2))")

    #     s1 = [round(f(xx[z]), digits=5) for f in s[q,:]]
    #     s2 = [round(f(xx[z]), digits=5) for f in s[q+1,:]]
    #     # s3 = [round(f(xx[z]), digits=5) for f in s[q+2,:]]
    #     # s4 = [round(f(xx[z]), digits=5) for f in s[q+3,:]]

    #     println("s1($(xx[z])) = ", s1)
    #     println("s2($(xx[z])) = ", s2)
    #     println()
    #     # println("s_[i]''($(xx[z])) = ", s3)
    #     # println("s_[i+1]''($(xx[z])) = ", s4)
    #     # println()
        
    #     m[p, r1] = [round(f(xx[z]), digits=5) for f in s[q,:]]
    #     m[p, r2] = [round(-f(xx[z]), digits=5) for f in s[q+1,:]]
        
        # for z in 2:nsp
            
            # m[p[i]+1, r1] = [round(f(xx[z]), digits=5) for f in s[q+2,:]]
            # m[p[i]+1, r2] = [round(-f(xx[z]), digits=5) for f in s[q+3,:]]
        # end
    # end
    
#     println()
#     for r in 1:size(m,1)
#         println(m[r,:])
#     end
# end



# # polynomial derivative transformation matrix
# function ddx(; polys)
#     pdr = 1:length(polys)
#     m = hcat([0 for i in pdr], diagm([i for i in pdr]))
# end


# # evaluate s_i(x)
# function seval(; xi)
#     s = x -> [(x-xi)^i for i in 0:]


# data = [(-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101)]
# csplines(data=data)

data = [(0.1, log(0.1)), (1, log(1)), (2, log(2)), (2.9, log(2.9))]
cspl(data=data)

# data = [(1,1), (2,1), (3,2), (4,6), (5,24)]
# csplines(data=data)
