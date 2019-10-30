using Memoize
using BenchmarkTools

# Slightly less slow. Slightly less ugly.
function matmul(a,b)
    m = size(a, 1)
    n = size(a, 2)
    p = size(b, 1)
    q = size(b, 2)

    if n==p
        return [sum(a[i,k]*b[k,j] for k in 1:n) for i in 1:m, j in 1:q]
    end
end

# Slow and ugly.
function matmul2(a,b)
    m = size(a, 1)
    n = size(a, 2)
    p = size(b, 1)
    q = size(b, 2)

    if n==p
        c = zeros(m, q)
        for i in 1:m
            for j in 1:q
                for k in 1:n 
                    c[i,j] += a[i,k]*b[k,j]
                end
            end
        end
        return c
    end
end

# Slow and ugly.
function matmul3(a,b)
    m = size(a, 1)
    n = size(a, 2)
    p = size(b, 1)
    q = size(b, 2)

    if n==p
        c = zeros(m, q)
        for j in 1:q
            for k in 1:n 
                for i in 1:m
                    c[i,j] += a[i,k]*b[k,j]
                end
            end
        end
        return c
    end
end

function fac(n)
    if n < 2
        return 1
    end
    return n*fac(n-1)
end

@memoize function fib(n)
    if n < 2
        return n
    end
    return fib(n-1) + fib(n-2)
end

function slowfib(n)
    if n < 2
        return n
    end
    return slowfib(n-1) + slowfib(n-2)
end

# @btime fib(35)
# @btime slowfib(35)

# for i = 10:10:500
#     a = rand(-9:9, i, i)
#     b = rand(-9:9, i, i)
#     println(matmul3(a,b) == matmul2(a,b) == matmul(a,b) == a*b)
#     println()
#     @time matmul3(a,b)
#     @time matmul2(a,b)
#     @time matmul(a,b)
#     @time a*b
#     println()
# end

n = 500
a = rand(-9:9, n, n)
b = rand(-9:9, n, n)

t = @elapsed matmul(a,b)

println((2*n^3 - n^2)/t)