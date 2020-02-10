using Plots

mutable struct ThreeVec <: AbstractVector{Float64}
    x::Float64
    y::Float64
    z::Float64
end

Base.size(V::ThreeVec) = (3,)
Base.IndexStyle(::Type{<:ThreeVec}) = IndexLinear()
Base.getindex(V::ThreeVec, i::Int) = i==1 ? V.x : (i==2 ? V.y : (i==3 ? V.z : error("index out of bounds")))
Base.setindex!(V::ThreeVec, val::Float64, i::Int) = (i==1 ? V.x = val : (i==2 ? V.y = val : (i==3 ? V.z = val : error("index out of bounds"))))

mutable struct Projectile
    r::ThreeVec
    v::ThreeVec
    dt::Float64
end

function ODEstep!(p::Projectile; f, wind=nowind, db=false)
    a = acc(p.v - wind(p.r), db)
    p.r, p.v = f(p.r, p.v, a, p.dt)
end

function euler(rn::ThreeVec, vn::ThreeVec, an::ThreeVec, dt::Float64)
    v = vn + an*dt
    r = rn + vn*dt
    return(r, v)
end

function cromer(rn::ThreeVec, vn::ThreeVec, an::ThreeVec, dt::Float64)
    v = vn + an*dt
    r = rn + v*dt
    return(r, v)
end

function midpoint(rn::ThreeVec, vn::ThreeVec, an::ThreeVec, dt::Float64)
    v = vn + an*dt
    r = rn + vn*dt + 1/2*a*dt^2
    return(r, v)
end

# Returns an acceleration vector, given velocity and a drag boolean (defaults to false).
function acc(v::ThreeVec; db=false)
    if db
        return drag_mult*sqrt(sum(v.^2))*v + ThreeVec(0., 0., -g)
    else
        return ThreeVec(0., 0., -g)
    end
end

function nowind(r::ThreeVec)
    return ThreeVec(0.,0.,0.)
end

vi = 15.
inc = 45.
p = Projectile(ThreeVec(0.,0.,0.), ThreeVec(vi*sind(inc), 0., vi*cosd(inc)), 0.1)

plt = plot3d(1, xlim=(-25,25), ylim=(-25,25), zlim=(0,50), title = "Projectile Test", marker = 2)

# @gif for i in 1:20 # while true
#     ODEstep!(p, f = euler)
#     push!(plt, p.r.x, p.r.y, p.r.z)
#     # if p.r.z < 0
#     #     break
#     # end
# end every 1
