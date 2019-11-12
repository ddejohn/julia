using LinearAlgebra
using Plots
gr()


function power_method(A, x, tol=1.0e-8)    
    n = 0
    mu = undef
    err_plot = []
    err_n = []
    errs = [1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8]

    while true
        x0 = x
        xn = A*normalize(x)
        x = normalize(xn)
        mu = x'*xn

        n += 1
        err = norm(x0-x, Inf)
        push!(err_plot, err)

        if err < errs[1]
            push!(err_n, n)
            popfirst!(errs)
        end

        if err < tol
            break
        end
    end

    println("eigenvector: $(x)\neigenvalue: $(mu)\niterations: $(n)")

    err_n_plot = [err_plot[i] for i in err_n]
    plot(collect(1:n), err_plot, yaxis=:log, lw=4, lc=:teal, label="",
        yticks=[1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8],
        xticks=[i for i in 0:5:n],
        ylabel="absolute error", xlabel="iterations", size=(900,900),
        title="Maximum eigenvalue error vs number of iterations")
    scatter!(err_n, err_n_plot, mc=:gold, ms=6, label="")
end


function qr_decomp(A, tol=1.0e-8)
    n = 0
    err_plot = []
    err_n = []
    errs = [1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8]

    while true
        A = qr(A)
        A = A.R*A.Q
        
        n += 1
        err = norm(diag(A, -1), 1)
        push!(err_plot, err)

        if err < errs[1]
            push!(err_n, n)
            popfirst!(errs)
        end

        if err < tol
            break
        end
    end

    println("eigenvalues: $(sort(diag(A)))\niterations: $(n)")

    err_n_plot = [err_plot[i] for i in err_n]
    plot(collect(1:n), err_plot, yaxis=:log, lw=4, lc=:teal, label="",
        yticks=[1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8],
        xticks=[i for i in 0:5:n],
        ylabel="absolute error", xlabel="iterations", size=(900,900),
        title="Off-diagonal error vs number of iterations")
    scatter!(err_n, err_n_plot, mc=:gold, ms=6, label="")
end;