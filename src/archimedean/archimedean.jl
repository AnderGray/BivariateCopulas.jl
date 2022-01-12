abstract type ArchimedeanCopula <: AbstractCopula end

function rosenblatt(M::AbstractMatrix, c::ArchimedeanCopula)
    u1 = M[:, 1]
    u2 = M[:, 2]

    return [u1 (φ².(φ⁻¹.(u1, c) + φ⁻¹.(u2, c), c) ./ φ².(φ⁻¹.(u1, c), c))]
end

function cdf(c::ArchimedeanCopula, u1::Real, u2::Real)
    return φ(φ⁻¹(u1, c) + φ⁻¹(u2, c), c)
end

function density(c::ArchimedeanCopula, u1::Real, u2::Real)
    return φ²(φ⁻¹(u1, c) + φ⁻¹(u2, c), c) * ∇φ⁻¹(u1, c) * ∇φ⁻¹(u2, c)
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
