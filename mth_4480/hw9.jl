using LinearAlgebra


# helper functions to use in iteration matrix
jacobi_M(A) = inv(diagm(diag(A)))
gauss_M(A) = inv(tril(A))


# Perform Gauss-Seidel or Jacobi iteration on matrix A
function lin_iter(; mthd, A, b, stop, xi=false, tol=1.0e-7)
    if !xi
        xi = [0. for i in 1:size(A,1)]
    end

    M = mthd(A)
    T = I - M*A
    c = M*b

    n = 0
    nb = norm(b,1)

    # iteration equation and remainder equation
    f(x) = T*x + c
    r(x) = b - A*x

    # stopping criteria boolean check helper function
    diff(xi) = norm(f(xi)-xi, 1) > tol
    ratio(xi) = norm(r(xi), 1)/nb > tol

    # which stopping criteria to use
    if stop == "r"
        bool = ratio
    else
        bool = diff
    end

    while bool(xi)
        xi = f(xi)
        n += 1
    end

    return round.(xi,digits=8), n
end


# compare methods
function lin_compare(A, b)
    @time soln_rj, nrj = lin_iter(mthd=jacobi_M, A=A, b=b, stop="r")
    @time soln_dj, ndj = lin_iter(mthd=jacobi_M, A=A, b=b, stop="d")

    @time soln_rg, nrg = lin_iter(mthd=gauss_M, A=A, b=b, stop="r")
    @time soln_dg, ndg = lin_iter(mthd=gauss_M, A=A, b=b, stop="d")

    @time soln = A\b

    jerr = norm(soln-soln_rj, 1)
    gerr = norm(soln-soln_rg, 1)
    
    println("\nstopping criteria: ratio")
    println("jacobi # iter: $(nrj)")
    println("gauss # iter: $(nrg)\n")

    println("stopping criteria: diff")
    println("jacobi # iter: $(ndj)")
    println("gauss # iter: $(ndg)\n")

    println("jacobi error: $(jerr)")
    println("gauss error: $(gerr)\n")
    
    println("jacobi solution: $(soln_rj)")
    println("gauss solution: $(soln_rg)")
    println("actual solution: $(round.(soln, digits=8))")
end