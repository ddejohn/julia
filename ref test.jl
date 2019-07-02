using LinearAlgebra
using BenchmarkTools

# Row-echelon form
function ref(a::Array{Float64,2})
    m, n = size(a,)
    for i in 1:m
        for j in 1:(m-i)
            a[i+j,:] = (a[i,i]/a[i+j,i])*a[i+j,:] - a[i,:]
            a[i+j,:] = round.(a[i+j,:], digits=5)
        end
    end
    a
end

# Reduced row-echelon form
# WARNING: LITERALLY THE WORST CODE I'VE EVER WRITTEN
function rref(a::Array{Float64,2})
    m, n = size(a,)
    # forward reduction
    for i in 1:m
        if abs(a[i,i]) > 1e-15
            for j in 1:(m-i)
                a[i+j,:] = (a[i,i]/a[i+j,i])*a[i+j,:] - a[i,:]
            end
        else a[i,i] = 0.
        end
    end
    # backward reduction
    for i in m:-1:1
        if abs(a[i,i]) > 1e-15
            for j in i-1:-1:1
                a[i-j,:] = (a[i,i]/a[i-j,i])*a[i-j,:] - a[i,:]
            end
        else a[i,i] = 0.
        end
    end
    for i in 1:m
        for j in 1:m
            if abs(a[i,j]) < 1e-10 a[i,j] = 0.
            end
        end
    end
    a
end

function test_det(a::Array{Float64,2})
    m, n = size(a,)
    if m == n
        res = 1.
        for i in 1:n
            res *= a[i,i]
        end
    end
    res
end

a = ref(100*rand(5, 5))
# a = rref([i*j for i in 1.:5., j in 1:5])

for i in 1:size(a,1)
    println(a[i,:])
end

t1 = @btime det(a)
t2 = @btime test_det(a)

println("LA: $t1 s")
println("My: $t2 s")