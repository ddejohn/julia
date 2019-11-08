using LinearAlgebra


function lam_iter(x)
    xn = A*x
    lam = maximum(xn)
    x = xn ./ lam
    return x, lam
end


function power_method(A, x, tol=1.0e-8)
    x0, mu = pow_help(A, x)
    n = 0

    while norm(x0-x, Inf) > tol
        x0 = x
        x, mu = pow_help(A, x)
        n += 1
    end

    println("eigenvector: $(x)\neigenvalue: $(mu)\niterations: $(n)")
end


function pow_help(A, x)
    x = x ./ norm(x, 2)
    xn = A*x
    x = xn ./ norm(xn, 2)
    return x, x'*xn
end


A = [4. 2. 1.; 0. 3. 2.; 1. 1. 4.]
x = [1., 2., 1.]

# A = [1. -1. 0.; -2. 4. -2.; 0. -1. 2.]
# x = [-1., 2., 1.]

power_method(A, x)