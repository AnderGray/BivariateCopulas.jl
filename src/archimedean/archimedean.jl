abstract type ArchimedeanCopula <: AbstractCopula end

function rosenblatt(M::AbstractMatrix, c::ArchimedeanCopula)
    u1 = M[:, 1]
    u2 = M[:, 2]

    return [u1 (D¹φ.(φ⁻¹.(u1, c) + φ⁻¹.(u2, c), c) ./ D¹φ.(φ⁻¹.(u1, c), c))]
end

function inverse_rosenblatt(U::AbstractMatrix, c::ArchimedeanCopula)
    u2 = map(eachrow(U)) do u
        return find_zero(x -> rosenblatt([u[1] x], c)[2] - u[2], (0, 1), Roots.A42())
    end

    return [U[:, 1] u2]
end

function cdf(c::ArchimedeanCopula, u1::Real, u2::Real)
    return φ(φ⁻¹(u1, c) + φ⁻¹(u2, c), c)
end

function density(c::ArchimedeanCopula, u1::Real, u2::Real)
    if iszero(u1) || iszero(u2)
        return zero(u1)
    end
    return D²φ(φ⁻¹(u1, c) + φ⁻¹(u2, c), c) * D¹φ⁻¹(u1, c) * D¹φ⁻¹(u2, c)
end

function sample(c::ArchimedeanCopula, n::Int64)
    @assert n > 0
    return inverse_rosenblatt(rand(n, 2), c)
end

function (c::ArchimedeanCopula)(
    X::ContinuousUnivariateDistribution, Y::ContinuousUnivariateDistribution
)
    # Construct a joint distribution using the copula
    return Joint(X, Y, c)
end

include("independence.jl")
include("clayton.jl")
include("frank.jl")
