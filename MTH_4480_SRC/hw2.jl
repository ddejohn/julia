function bisection(; f, intvl, tol=1e-2)
    a, b = intvl
    # vals = []
    
    if f(a)*f(b) >= 0
        return "initial inputs are not opposite sign"
    end

    x = (a+b)*0.5
    fx = f(x)
    n = 1
    # push!(vals, fx)

    while abs(fx) > tol
        sign(f(a)) == sign(fx) ? a = x : b = x

        x = (a+b)*0.5
        fx = f(x)
        n += 1
        # push!(vals, x)
    end
    
    println(
        "x: $(x)\nf(x) = $(f(x))\nnumber of iterations: $(n)\n",
    )
end;


println("Bisection:\n")
# bisection(f=x-> x^4 - 2x^3 - 4x^2 + 4x + 4, intvl=(-2, -1))
# bisection(f=x-> x^4 - 2x^3 - 4x^2 + 4x + 4, intvl=(-1, 0))
# bisection(f=x-> x^4 - 2x^3 - 4x^2 + 4x + 4, intvl=(0, 2))
# bisection(f=x-> x^4 - 2x^3 - 4x^2 + 4x + 4, intvl=(2, 3)) 


function fixed_point(; f, g, p0, max_iter=100, tol=1e-2)
    p = g(p0)
    n = 1

    while n < max_iter
        p = g(p0)
        if abs(p - p0) < tol
            break
        end    
        p0 = p
        n += 1
    end    

    if n == max_iter
        println("failed after $(n) iterations\n")
    else
        println(
            "p: $(p)\nf(p) = $(f(p))\nnumber of iterations: $(n)\n",
        )
    end    
end;


# println("Fixed point:\n")
# fixed_point(f=x-> x^3, g=x-> (18x + 19/x^2)/19, p0=1)
# fixed_point(f=x-> x^3, g=x-> x - (x^3 - 19)/(3x^2), p0=1, tol=1e-5)
# fixed_point(f=x-> x^3, g=x-> x - (x^4 - 19x)/(x^2 - 19), p0=-1)
# fixed_point(f=x-> x^3, g=x-> sqrt(19/x), p0=1)
# fixed_point(f=x-> x^3, g=x-> sqrt((19+x^2)/(1+x)), p0=1, tol=1e-5)

# g1 = x-> x - (x^3 - 19)/(3x^2)
# g2 = x-> sqrt((19+x^2)/(1+x))

# for i in 1:15
#     tol = 10.0^(-i)
#     g1_i = fixed_point(f=x-> x^3, g=g1, p0=1, tol=tol)
#     g2_i = fixed_point(f=x-> x^3, g=g2, p0=1, tol=tol)
#     println("tolerance: $(tol)\ng1: $(g1_i)\ng2: $(g2_i)")
# end