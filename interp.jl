# Lagrange polynomial interpolation.
function intrp(xx, yy; n=5)
    k = length(xx)
    if k != length(Set(xx))
        println("Cannot interpolate from duplicate x-values!")
    end

    function lbp(j, x)
        p = 1
        for m in (1:k)[1:end .!= j]
            p *= (x - xx[m])/(xx[j] - xx[m])
        end
        return p
    end

    return x -> round(sum(yy[j]*lbp(j, x) for j in 1:k), digits=n)
end

x = [1, -2, 3, 2, -3, -1]
y = [1, -8, 27, 8, -27, -1]

i = intrp(x, y)
println(i(5))

x = [1, 2, 3, -2, -1, -3]
y = [1, 4, 9, 4, 1, 9]

i = intrp(x, y)
println(i(5))

x = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
y = [1.0, 0.99, 0.9604, 0.9120, 0.8463, 0.7652]
i = intrp(x, y)

for xi in [0.3, 0.9, 1.1, 1.5, 2.0]
    println(i(xi))
end