using Plots
gr()

const μ = 4*pi^2

#——————————————————————————————————————————————————————————————————————————————————————————————————#
# Cartesian three-vector

mutable struct ThreeVec
    x::Number
    y::Number
    z::Number
end

Base.:*(k::Number, v::ThreeVec) = ThreeVec(k*v.x, k*v.y, k*v.z)
Base.:+(u::ThreeVec, v::ThreeVec) = ThreeVec(u.x + v.x, u.y + v.y, u.z + v.z)

#——————————————————————————————————————————————————————————————————————————————————————————————————#
# Orbital state vector

mutable struct State
    r::ThreeVec
    v::ThreeVec
end

Base.:*(k::Number, s::State) = State(k*s.r, k*s.v)
Base.:+(a::State, b::State) = State(a.r + b.r, a.v + b.v)

#——————————————————————————————————————————————————————————————————————————————————————————————————#

function step!(s::State, t::Float64)
    a = cent_acc(s.r)
    s.r += t*s.v
    s.v += t*a
end

# ThreeVec norm
nrm(V::ThreeVec) = sqrt(V.x^2 + V.y^2 + V.z^2)

# Already mass-normalized (little m divided out)
cent_acc(r::ThreeVec) = -(μ/nrm(r)^3) * r

# X' = [dr, dv]
dX(X::State) = State(X.v, cent_acc(X.r))

function rk4!(X::State, τ::Float64)
    k1 = dX(X)
    k2 = dX(X + τ/2*k1)
    k3 = dX(X + τ/2*k2)
    k4 = dX(X + τ*k3)

    X.r += τ/6 * (k1.r + 2*k2.r + 2*k3.r + k4.r)
    X.v += τ/6 * (k1.v + 2*k2.v + 2*k3.v + k4.v)
end

function test_plot()
    r1 = ThreeVec(1., 0., 0.)
    v1 = ThreeVec(0., pi, 0.)

    s = State(r1, v1)

    x = []
    y = []

    for i in 1:200
        rk4!(s, 0.005)
        push!(x, s.r.x)
        push!(y, s.r.y)
    end
    p = plot([0], [0], proj=:polar, size=(900,900))
    plot!(p, x, y)
    p
end

test_plot()

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# mutable struct ThreeVec <: AbstractVector{Float64}

# mutable struct State <: AbstractVector{Float64}

# Base.size(S::State) = (2,3)
# Base.IndexStyle(::Type{<:State}) = IndexCartesian()

# Base.size(V::ThreeVec) = (3,)
# Base.IndexStyle(::Type{<:ThreeVec}) = IndexLinear()

# Base.getindex(S::State, i::Tuple{Int, 2}) = begin
#     if i[1] == 1
#         S.r[i[2]]
#     elseif i[1] == 2
#         S.v[i[2]]
#     else
#         error("Get error: index out of bounds")
#     end
# end

# Base.setindex!(S::State, val::Float64, i::Tuple{Int, 2}) = begin
#     if i[1] == 1
#         S.r[i[2]] = val
#     elseif i[1] == 2
#         S.v[i[2]] = val
#     else
#         error("Set error: index out of bounds")
#     end
# end


# Base.size(S::State) = (2,)
# Base.IndexStyle(::Type{<:State}) = IndexLinear()

# Base.getindex(S::State, i::Int) = begin
#     i==1 ? S.r :
#         (i==2 ? S.v :
#             error("Get error: index out of bounds"))
# end

# Base.setindex!(S::State, val::ThreeVec, i::Int) = begin
#     i==1 ? S.r = val :
#         (i==2 ? S.v = val :
#             error("Set error: index out of bounds"))
# end


# Base.getindex(V::ThreeVec, i::Int) = begin
#     i==1 ? V.x :
#         (i==2 ? V.y :
#             (i==3 ? V.z :
#                 error("Get error: index out of bounds")))
# end

# Base.setindex!(V::ThreeVec, val::Float64, i::Int) = begin
#     i==1 ? V.x = val :
#         (i==2 ? V.y = val :
#             (i==3 ? V.z = val :
#                 error("Set error: index out of bounds")))
# end
# s3 = s1 + s2
# s4 = 5*s3

# println(s3)

# println(s3.r)
# println(s3.v)
# println(s4.r)
# println(s4.v)

# nrm = sqrt(sum(s3.r.^2))
# println(nrm)
