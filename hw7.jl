function mid_point(; f, intvl)
    a, b = intvl
    a, b = minmax(a, b)
    return (b-a)*f((a+b)/2)
end


function trap(; f, intvl)
    a, b = intvl
    a, b = minmax(a, b)
    return (b-a)*(f(a)+f(b))/2
end


function simpson(; f, intvl)
    a, b = intvl
    a, b = minmax(a, b)
    m = (a+b)/2
    return (b-a)*(f(a)+4*f(m)+f(b))/6
end


function integrate(; f, intvl, mthd, n)
    a, b = intvl
    a, b = minmax(a, b)
    subs = collect(a:(b-a)/n:b)

    area = 0

    for i in 1:length(subs)-1
        area += mthd(f=f, intvl=(subs[i+1], subs[i]))
    end
    
    return area
end


f(x) = 4/(1+x^2)
# f(x) = sqrt(x)

# for n in [2,4,8,16,32]
#     println(integrate(f=f, intvl=(0,1), mthd=mid_point, n=n))
# end


println("n = 32\n")
for mthd in [mid_point, trap, simpson]
    area = integrate(f=f, intvl=(0,1), mthd=mthd, n=32)
    println("$(string(mthd)): $(round(area, digits=6))")
    println("error: $(round(abs(pi - area), digits=10))\n")
end