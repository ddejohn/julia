using Pkg
Pkg.add("Plots")
Pkg.add("GR")
using Plots
gr()

function matmul(a,b)
    m = size(a, 1)
    n = size(a, 2)
    p = size(b, 1)
    q = size(b, 2)

    if n == p
        return [sum(a[i,k]*b[k,j] for k in 1:n) for i in 1:m, j in 1:q]
    end
end

function matexp(a,n)
    if size(a,1) == size(a,2)
        if n > 0
            if n == 1
                return a
            end
            return matmul(a, matexp(a, n-1))
        else
            return a^0
        end
    else
        println("Cannot exponentiate asymmetric matrix!")
    end
end

function vdot(u, w)
    return sum(i*j for (i,j) in zip(u,w))
end

function f1(x)
    return exp(-x/4)*sin(x)
end

function f2(x)
    return abs(exp(-x)*sin(x))
end

# int p is number of petals, int w is 'pinch' factor; larger w = skinnier petals
function flower(x, p, w)
    return abs(sin(p/2*x)^w)
end

# 1.1 Exercises

#1
println("Hello, world!")

#2 - Basic matrix manipulation.

A = [1 2 ; 3 4]
B = [1 2 ; 3 4]

println("a) ", A+A)
println("b) ", B+B)
println("c) ", A+B)
println("d) ", A-A)
println("e) ", B-B)
println("f) ", 2*A)
println("g) ", 2*B)
println("h) ", matexp(A, 2))
println("i) ", matexp(B, 2))
println("j) ", B*B)
println("k) ", B*B)
println("l) ", B/B)

#3 - Basic plotting.

x = [i for i in 1.0:10.0]
y = x.^2
p1 = scatter(x, y, label="x^2")
plot(p1)

#4 - Reproduce plots from book.

# Decaying periodic function + decay envelope.
x = [convert(Float64, i) for i in 0:0.001:20]
y = [f1.(x), exp.(-x/4)]
p2 = plot(x, y, label=["exp(-x/4)*sin(x)" "exp(-x/4)"])
plot(p2)

# log y-axis plot.
x = [convert(Float64, i) for i in 0.001:0.25:18.85]
p3 = scatter(x, f2.(x), label="exp(-x/4)*sin(x))", m=:x, yscale=:log10)
plot(p3)

# Rhodonea curve on polar grid.
x = range(0, stop=2pi, length = 1000)
r = flower.(x, 6, 2)
p4 = plot(x, r, proj=:polar, label="|sin(3x)^2|")
plot(p4)
