function bisection(; f, intvl, tol=1e-2)
    a, b = intvl
    
    if f(a)*f(b) >= 0
        return "initial inputs are not opposite sign"
    end

    x = (a+b)*0.5
    fx = f(x)
    n = 1

    while abs(fx) > tol
        sign(f(a)) == sign(fx) ? a = x : b = x

        x = (a+b)*0.5
        fx = f(x)
        n += 1
    end
    
    println(
        "x: $(x)\nf(x) = $(f(x))\nnumber of iterations: $(n)\n",
    )
end


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
end


function newton(; f, df, xn, max_iter=100, tol=1e-4)
    n=1

    while n <= max_iter
        x = xn - f(xn)/df(xn)

        if abs(x-xn) < tol || abs(f(x)) < tol
            xn = x
            break
        end
        
        n += 1
        xn = x
    end

    if (n == max_iter) & (abs(f(xn)) > tol)
        out = "failed to find root after $(max_iter) iterations.\n\n"
    else
        out = "root found near x = $(xn) in $(n) iterations.\n\n"
    end

    println(out)
end


function secant(; f, x0, x1, max_iter=100, tol=1e-4)
    n=1

    while n <= max_iter
        fx0, fx1 = f(x0), f(x1)
        x = x1 - fx1*(x1-x0)/(fx1-fx0)

        if abs(x-x1) < tol || abs(f(x)) < tol
            x1 = x
            break
        end
        
        n += 1
        x0 = x1
        x1 = x
    end

    if (n == max_iter) & (abs(f(x1)) > tol)
        out = "failed to find root after $(max_iter) iterations.\n\n"
    else
        out = "root found near x = $(x1) in $(n) iterations.\n\n"
    end

    println(out)
end