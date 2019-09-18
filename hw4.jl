using LinearAlgebra
using Plots
gr()


# Newton's divided difference polynomial
function ndd(; data)
    n = length(data)
    # initialize nxm matrix
    m = zeros(Float64, n, n)
    
    # unzips list of (x,y) coords into separate lists
    # set first column of DD matrix to f(x) for x in xx
    xx, m[:,1] = map(y -> map(x -> x[y], data), 1:2)

    # generate lower triangular NDD matrix
    for i in 2:n
        for j in 2:i
            m[i,j] = (m[i,j-1]-m[i-1,j-1])/(xx[i]-xx[i-j+1])
        end
    end
    
    # return a lambda in x: 
    # sum the diagonal entries as coefficients
    # on the generated Newton polynomial
    return x-> sum(m[i,i]*poly(i, x, xx) for i in 1:n), hcat(xx, m)
end


# Helper function to generate Newton polynomials
function poly(i, x, xx)
    p = 1
    for j in 1:i-1
        p *= (x-xx[j]) # (x - xj) for j from 1 to i-1
    end
    return p
end


# helper function to add new point to existing divided difference matrix
function ndd_add(; m::Matrix, p::Tuple)
    # concatenate m with a row of zeros at the bottom
    m = vcat(m, fill(0.0, size(m, 2))')
    
    # concatenate m with a column of zeros on the right
    m = hcat(m, fill(0.0, size(m, 1)))
    
    # unpack x,y into first two columns of m
    m[end,1], m[end,2] = p

    # perform divided differences on the last row
    for j in 3:size(m,2)
        m[end, j] = (m[end,j-1]-m[end-1,j-1])/(m[end, 1]-m[end-j+2,1])
    end

    # same return type as ndd
    return x-> sum(m[i,i+1]*poly(i, x, m[:, 1]) for i in 1:size(m,1)), m
end


# print the interpolating polynomial
function ndd_print(m::Matrix)
    xx = m[:,1]
    diag = [round(m[i,i+1], digits=3) for i in 1:size(m, 1)]
    ndd_repr = []
    for i in 1:length(diag)
        push!(ndd_repr, "$(diag[i])")
        for j in 1:i-1
            if xx[j] < 0
                ndd_repr[i] *=  "*(x+$(-xx[j]))"
            else
                ndd_repr[i] *=  "*(x-$(xx[j]))"
            end
        end
    end
    return "p(x) = " * join(ndd_repr, " + ")
end


# Lagrange polynomial interpolation
function lag(; data)
    n = length(data)
    xx, yy = map(y -> map(x -> x[y], data), 1:2)

    # Lagrange basis polynomial helper
    function lbp(j, x)
        p = 1
        for i in (1:n)[1:end .!= j]
            p *= (x - xx[i])/(xx[j] - xx[i])
        end
        return p
    end

    return x -> sum(yy[j]*lbp(j, x) for j in 1:n)
end


# monomial polynomial interpolation (silly name)
function mono(; data)
    n = length(data)
    xx, yy = map(y -> map(x -> x[y], data), 1:2)

    m = [x^p for x in xx, p in 0:n-1]
    cc = m\yy
    
    return x -> sum(cc[i]*x^(i-1) for i in 1:length(cc))
end


# plot any number of interpolated functions against f(x) 
function make_plot(; f, funcs, names, x_range::StepRangeLen)
    # x-range to plot
    xx = [x for x in x_range]
    
    p = plot(
        xx, f.(xx), label="f(x) = 2/(1-x)",
        lw=3, xlabel="x", ylabel="y",
        size=(900,900),
        legend=:bottomright,
        title="f(x) vs " * join(names, ", ")
    )
    
    for i in 1:length(funcs)
        # unpack functions into degree 3, 4
        f3, f4 = funcs[i]
        name = names[i]
        lbl = "$(name): degree "
        plot!(p, xx, f3.(xx), lw=3, label=lbl*"3")
        plot!(p, xx, f4.(xx), lw=3, label=lbl*"4")
    end
    
    scatter!(p, [-0.1, 0, 0.2, 0.3, 0.35], [1.81818, 2., 2.5, 2.85714, 3.07692],
        ms=6, mc="lightgrey", label="given data points")
    p
end


# error plotter
function error_plot(; f, funcs, names, x_range::StepRangeLen)
    # x-range to plot
    xx = [x for x in x_range]
    
    p = plot(
        [], [], label="",
        lw=3, xlabel="x", ylabel="absolute error",
        size=(900,900),
        legend=:topleft,
        title="Absolute error f(x) vs " * join(names, ", ")
    )
    
    for i in 1:length(funcs)
        # unpack functions into degree 3, 4
        f3, f4 = funcs[i]
        name = names[i]
        lbl = "$(name): degree "
        plot!(p, xx, abs.(f.(xx) - f3.(xx)), lw=3, label=lbl*"3")
        plot!(p, xx, abs.(f.(xx) - f4.(xx)), lw=3, label=lbl*"4")
    end
    
    scatter!(p, [-0.1, 0, 0.2, 0.3, 0.35], [0., 0., 0., 0., 0.],
        ms=6, mc="lightgrey", label="given data points")
    p
end


# subplots for degree 3 and degree 4 polynomials
function multi_plot(; f, funcs, names, x_range::StepRangeLen)
    # x-range to plot
    xx = [x for x in x_range]
    
    p3=plot([], [], label="", title="degree 3 polynomial", xlabel="x", ylabel="y")
    p4=plot([], [], label="", title="degree 4 polynomial", xlabel="x", ylabel="y")
    
    for i in 1:length(funcs)
        # unpack functions into degree 3, 4
        f3, f4 = funcs[i]
        name = names[i]
        plot!(p3, xx, f3.(xx), lw=3, label=name)
        plot!(p4, xx, f4.(xx), lw=3, label=name)
    end

    scatter!(p3, [-0.1, 0, 0.2, 0.3, 0.35], [1.81818, 2., 2.5, 2.85714, 3.07692],
        ms=6, mc="lightgrey", label="given data points")

    scatter!(p4, [-0.1, 0, 0.2, 0.3, 0.35], [1.81818, 2., 2.5, 2.85714, 3.07692],
        ms=6, mc="lightgrey", label="given data points")

    p=plot(p3, p4, layout=(1,2), size=(950,550), legend=:topleft)
    p
end;