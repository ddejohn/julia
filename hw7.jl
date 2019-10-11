using Plots
gr()


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


# plot errors for all three methods
function error_plot(; f, exact, intvl)
    mdp_err = []
    tpz_err = []
    smp_err = []
    
    nvals = [2, 4, 6, 8, 16, 32, 64]
    
    p = plot([], [], label="", title="absolute error vs # subintervals", size=(900,900),
        xlabel="# subintervals", ylabel="absolute error", xticks=nvals, yaxis=:log)
    
    for n in nvals
        mdp = abs(exact - integrate(f=f, intvl=intvl, mthd=mid_point, n=n))
        tpz = abs(exact - integrate(f=f, intvl=intvl, mthd=trap, n=n))
        smp = abs(exact - integrate(f=f, intvl=intvl, mthd=simpson, n=n))
        
        push!(mdp_err, mdp)
        push!(tpz_err, tpz)
        push!(smp_err, smp)
    end
    
    plot!(p, nvals, mdp_err, lc=:teal, lw=4, label="midpoint")
    plot!(p, nvals, tpz_err, lc=:gold, lw=4, label="trapezoid")
    plot!(p, nvals, smp_err, lc=:pink, lw=4, label="simpson's")
    
    p
end;