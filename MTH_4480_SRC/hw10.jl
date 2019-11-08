using LinearAlgebra
using Plots
gr()


function power_method(A, x, tol=1.0e-8)    
    n = 0
    mu = undef
    err_plot = []

    while true
        x0 = x
        normalize!(x)
        xn = A*x
        x = normalize(xn)
        mu = x'*xn

        n += 1
        err = norm(x0-x, Inf)
        push!(err_plot, err)
        
        if err < tol
            break
        end
    end

    println("eigenvector: $(x)\neigenvalue: $(mu)\niterations: $(n)")
    
    plot([x for x in 1:n], err_plot,
        yaxis=:log, lw=4, lc=:teal, label="",
        ylabel="error", xlabel="iterations",
        size=(900,900), title="Maximum eigenvalue error vs number of iterations")
    scatter!([x for x in 1:n], err_plot, mc=:gold, label="")
end;