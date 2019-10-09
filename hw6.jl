using Plots
gr()


# Forward, backward, and centered difference helper functions, vectorized
@. f_diff(; f, x, h) = (f(x+h) - f(x))/h
@. b_diff(; f, x, h) = (f(x) - f(x-h))/h
@. c_diff(; f, x, h) = (f(x+h) - f(x-h))/(2*h)


# Takes a function f, a list of x-values, and a step size h (defaults to 1.0e-3)
# and returns a numerically approximated derivative of f at the same x-values
# diffs calculated via array slicing and vectorized difference formula helper functions
function num_diff(; f, xx, h=1.0e-3)
    # initialize a vector of the same length as xx
    dx = Array{Float64}(undef, length(xx))
    
    # forward difference on first point
    dx[1] = f_diff(f=f, x=xx[1], h=h)
    
    # centered difference on interior points
    dx[2:end-1] = c_diff(f=f, x=xx[2:end-1], h=h)
    
    # backward difference on last point
    dx[end] = b_diff(f=f, x=xx[end], h=h)
    
    return dx
end


function error_plot(; f, df, xx, file_name)
    hh = [abs(10.0^i) for i in -15.0:0.0]
    yy = [abs(10.0^i) for i in -15.0:2.0]
    f_err = []
    b_err = []
    c_err = []

    for h in hh
        dx = num_diff(f=f, xx=xx, h=h)
        dfx = df.(xx)
        err = abs.(dfx .- dx)
        fe, ce, be = err
        
        println("h: $(h)\nfe: $(round(fe, digits=4))\nbe: $(round(be, digits=4))\nce: $(round(ce, digits=4))")
        
        push!(f_err, fe)
        push!(b_err, be)
        push!(c_err, ce)
    end

    p = plot([], [], xaxis=:log, yaxis=:log, xticks=hh, yticks=yy, label="",
        ylabel="Absolute error", xlabel="step size h", legend=:bottomleft,
        title="Absolute error as a function of step size", size=(900,900)
    )
    
    plot!(p, hh, f_err, lw=4, lc=:pink, label="forward error")
    plot!(p, hh, b_err, lw=4, lc=:gold, label="backward error")
    plot!(p, hh, c_err, lw=4, lc=:teal, label="centered error")
    
    savefig(p, file_name)
    
    p
end;