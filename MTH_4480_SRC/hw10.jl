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

    println("iterations: $(n)\neigenvalues:")
    for e in sort(diag(A))
        println("    $(e)")
    end

    err_n_plot = [err_plot[i] for i in err_n]
    plot(collect(1:n), err_plot, yaxis=:log, lw=4, lc=:teal, label="",
        yticks=[1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8],
        xticks=[i for i in 0:5:n],
        ylabel="absolute error", xlabel="iterations", size=(900,900),
        title="Off-diagonal error vs number of iterations")
    scatter!(err_n, err_n_plot, mc=:gold, ms=6, label="")
end


function shift_qr(A, tol=1.0e-8)
    n = 0
    eigs = []
    err_plot = []

    while true
        K = A[end,end]*I(size(A,1))
        A = qr(A - K)
        A = A.R*A.Q + K
        n += 1
        
        if size(A, 1) == 1
            push!(eigs, A[end,end])
            push!(err_plot, 1.0e-8)
            break
        end

        err = abs(A[end-1,end])
        push!(err_plot, err)

        if (abs(A[end-1,end]) < tol) & (abs(A[end,end-1]) < tol)
            push!(eigs, A[end,end])
            A = A[1:end-1, 1:end-1]
        end
    end

    println("iterations: $(n)\neigenvalues:")
    for e in eigs
        println("    $(e)")
    end

    plot(collect(1:n), err_plot, yaxis=:log, lw=4, lc=:teal, label="",
        ylabel="absolute error", xlabel="iterations", size=(900,900),
        title="QR+shift off-diagonal error vs number of iterations")
    scatter!(collect(1:n), err_plot, mc=:gold, ms=6, label="")
end;