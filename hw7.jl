using Printf


# integration helper functions
mid_point(; f, a, b) = (b-a)*f((a+b)/2)
trap(; f, a, b) = (b-a)*(f(a)+f(b))/2
simpson(; f, a, b) = (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6


# integrate a function 'f' over a given interval 'intvl'
# using 'mthd' mid_point, trap, or simpson, and number
# of sub-intervals 'n'
function integrate(; f, mthd, intvl::Tuple, n::Int)
    a, b = intvl # unpack the tuple
    a, b = minmax(a, b) # rearrange so a, b correspond to the correct endpoints
    subs = collect(a:(b-a)/n:b) # construct an array of sub intervals

    # composite area sum using given 'mthd'
    return sum(mthd(f=f, a=subs[i], b=subs[i+1]) for i in 1:length(subs)-1)
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
    @printf "error: %e\n\n" abs(pi - area)
end