using SpecialMatrices


# polynomial least squares
function lsfit(; data, n)
    xx, yy = map(y -> map(x -> x[y], data), 1:2)

    A = x.^collect(0:n)'
    coefs = (A'A)\(A'yy)

    return x -> sum(c*x^(i-1) for (i,c) in enumerate(coefs))
end


function oldfit(; data, n)
    n += 1
    xx, yy = map(y -> map(x -> x[y], data), 1:2)
    xrng = collect(minimum(xx):0.01:maximum(xx))
    
    A = Array{Float64}(undef, (n,n))
    b = [sum(yy[i]*xx[i]^j for i in 1:length(xx)) for j in 0:n-1]

    for i in 1:n
        for j in 1:n
            A[i,j] = sum(xx[k]^(i+j-2) for k in 1:length(xx))
        end
    end

    coefs = A\b

    return x -> sum(coefs[i]*x^(i-1) for i in 1:n)
end