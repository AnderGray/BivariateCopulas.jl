struct Clayton <: ArchimedeanCopula
    ϑ::Real

    function Clayton(ϑ::Real)
        # bivariate clayton is defined on ϑ ∈ (0, ∞)
        @assert 0 <= ϑ <= Inf

        if ϑ == 0.0
            @warn "Clayton returns an independence copula for ϑ == 0.0"
            return Independence()
        elseif ϑ == Inf
            @warn "Clayton returns an M copula for ϑ == Inf"
            return M()
        end

        return new(ϑ)
    end
end

"""
    φ(x::Real, c::Clayton)

Generator of the Clayton copula.
"""
function φ(x::Real, c::Clayton)
    return (1 + x)^(-1 / c.ϑ)
end

"""
    D¹φ(x::Real, c::Clayton)

First derivative of the Clayton copula generator.
"""
function D¹φ(x::Real, c::Clayton)
    α = 1 / c.ϑ
    return -α * (1 + x)^(-1 - α)
end

"""
    D²φ(x::Real, c::Clayton)

Second derivative of the Clayton copula generator.
"""
function D²φ(x::Real, c::Clayton)
    α = 1 / c.ϑ
    return (α + α^2) * (1 + x)^(-2 - α)
end

"""
    φ⁻¹(x::Real, c::Clayton)

Inverse generator of the Clayton copula.
"""
function φ⁻¹(x::Real, c::Clayton)
    return x^(-c.ϑ) - 1
end

"""
    D¹φ⁻¹(x::Real, c::Clayton)

First derivative of the Clayton copula inverse generator.
"""
function D¹φ⁻¹(x::Real, c::Clayton)
    return -c.ϑ * x^(-c.ϑ - 1)
end

function rosenblatt(M::AbstractMatrix, c::Clayton)
    u = 1 .+ cumsum(φ⁻¹.(M, c); dims=2)
    a = -1 / c.ϑ - 1

    return [M[:, 1] u[:, 2] .^ a ./ u[:, 1] .^ a]
end

function inverse_rosenblatt(U::AbstractMatrix, c::Clayton)
    u1 = U[:, 1]
    u2 = U[:, 2]

    return [u1 (u1 .^ (-c.ϑ) .* (u2 .^ (-c.ϑ / (c.ϑ + 1)) .- 1) .+ 1) .^ (-1 / c.ϑ)]
end

function τ(c::Clayton)
    return c.ϑ / (c.ϑ + 2)
end
