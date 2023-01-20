struct Frank <: ArchimedeanCopula
    ϑ::Real

    function Frank(ϑ::Real)
        # bivariate frank is defined on ϑ ∈ (0, ∞)
        @assert 0 <= ϑ < Inf

        if ϑ == 0.0
            @warn "Frank returns an independence copula for ϑ == 0.0"
            return Independence()
        end

        return new(ϑ)
    end
end

"""
    φ(x::Real, c::Frank)

Generator of the Frank copula.
"""
function φ(x::Real, c::Frank)
    return -1 / c.ϑ * log(exp(-x) * (exp(-c.ϑ) - 1) + 1)
end

"""
    D¹φ(x::Real, c::Frank)

First derivative of the Frank copula generator.
"""
function D¹φ(x::Real, c::Frank)
    z = (1 - exp(-c.ϑ)) * exp(-x)
    return -(1 / c.ϑ) * z / (1 - z)
end

"""
    D²φ(x::Real, c::Frank)

Second derivative of the Frank copula generator.
"""
function D²φ(x::Real, c::Frank)
    z = (1 - exp(-c.ϑ)) * exp(-x)
    return (1 / c.ϑ) * z / (1 - z)^2
end

"""
    φ⁻¹(x::Real, c::Frank)

Inverse generator of the Frank copula.
"""
function φ⁻¹(x::Real, c::Frank)
    return -log((exp(-c.ϑ * x) - 1) / (exp(-c.ϑ) - 1))
end

"""
    D¹φ⁻¹(x::Real, c::Frank)

First derivative of the Frank copula inverse generator.
"""
function D¹φ⁻¹(x::Real, c::Frank)
    return (c.ϑ * exp(-c.ϑ * x)) / (exp(-c.ϑ * x) - 1)
end

function sample(c::Frank, n::Int64)
    @assert n > 0
    E = rand(Exponential(1), (n, 2))

    # Sample from the logarithmic distribution with α = 1 - exp(-ϑ)
    # Kemp (1981) Efficient generation of logarithmically distributed pseudo-random variables
    α = 1 - exp(-c.ϑ)
    h = log(1 - α)
    u = rand(n, 2)
    M = trunc.(1 .+ log.(u[:, 2]) ./ log.(1 .- exp.(u[:, 1] .* h)))

    return [φ.(E[:, 1] ./ M, c) φ.(E[:, 2] ./ M, c)]
end

function τ(c::Frank)
    return 1 + 4 * (sf_debye_1(c.ϑ) - 1) / c.ϑ
end
