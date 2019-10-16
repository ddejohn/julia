using Plots
gr()


# integration helper functions
mid_point(; f, a, b) = (b-a)*f((a+b)/2)
trap(; f, a, b) = (b-a)*(f(a)+f(b))/2
simpson(; f, a, b) = (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6


# Abscissas for gaussian_quadrature
legendre_roots = [
    [sqrt(3)/3, -sqrt(3)/3],
    [0, sqrt(15)/5, -sqrt(15)/5],
    [sqrt(525-70*sqrt(30))/35, -sqrt(525-70*sqrt(30))/35, 
        sqrt(525+70*sqrt(30))/35, -sqrt(525+70*sqrt(30))/35
    ],
    [0, sqrt(245-14*sqrt(70))/21, -sqrt(245-14*sqrt(70))/21,
        sqrt(245+14*sqrt(70))/21, -sqrt(245+14*sqrt(70))/21
    ]
]


legendre_weights = [
    [1, 1],
    [8/9, 5/9, 5/9],
    [(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36],
    [128/225, (322+13*sqrt(70))/900, (322+13*sqrt(70))/900,
        (322-13*sqrt(70))/900, (322-13*sqrt(70))/900
    ]
]


# integrate a function 'f' over a given interval 'intvl'
# using 'mthd' mid_point, trap, or simpson, and number
# of sub-intervals 'n'
function integrate(; f, mthd, intvl::Tuple, n::Int, deg=3)
    a, b = intvl # unpack the tuple
    a, b = minmax(a, b) # rearrange so a, b correspond to the correct endpoints
    subs = collect(a:(b-a)/n:b) # construct an array of sub intervals

    if string(mthd) == "gauss"
        return sum(mthd(f=f, a=subs[i], b=subs[i+1], n=deg) for i in 1:length(subs)-1)
    end

    # composite area sum using given 'mthd'
    return sum(mthd(f=f, a=subs[i], b=subs[i+1]) for i in 1:length(subs)-1)
end


function gauss(; f, a, b, n=3)
    x(t) = 0.5*(t*(b-a)+a+b)
    xx = x.(legendre_roots[n-1])
    ww = legendre_weights[n-1]

    return 0.5*(b-a)*sum(wi*f(xi) for (wi, xi) in zip(ww, xx))
end


# plot errors for all three methods
function error_plot(; f, exact, intvl)
    mdp_err = []
    tpz_err = []
    smp_err = []
    gss_err = []
    
    nvals = [2, 4, 6, 8, 16, 32, 64]
    
    p = plot([], [], label="", title="absolute error vs # subintervals", size=(900,900),
        xlabel="# subintervals", ylabel="absolute error",
        xticks=nvals, xaxis=:log, xformatter=x -> round(x), yaxis=:log)
    
    for n in nvals
        mdp = abs(exact - integrate(f=f, intvl=intvl, mthd=mid_point, n=n))
        tpz = abs(exact - integrate(f=f, intvl=intvl, mthd=trap, n=n))
        smp = abs(exact - integrate(f=f, intvl=intvl, mthd=simpson, n=n))
        gss = abs(exact - integrate(f=f, intvl=intvl, mthd=gauss, n=n, deg=3))
        
        push!(mdp_err, mdp)
        push!(tpz_err, tpz)
        push!(smp_err, smp)
        push!(gss_err, gss)
    end
    
    plot!(p, nvals, mdp_err, lc=:teal, lw=4, label="midpoint")
    plot!(p, nvals, tpz_err, lc=:gold, lw=4, label="trapezoid")
    plot!(p, nvals, smp_err, lc=:pink, lw=4, label="simpson's")
    plot!(p, nvals, gss_err, lc=:cyan, lw=4, label="gauss")
    
    p
end;