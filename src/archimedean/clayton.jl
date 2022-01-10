struct Clayton <: ArchimedeanCopula
    ϑ::Real

    function Clayton(θ::Real)
        # bivariate clayton can go as low as -1
        # TODO: Return W and M for -1 and Inf
        @assert -1 < θ < Inf

        if θ == 0.0
            return Independence()
        end

        return new(θ)
    end
end

function φ(x::Real, c::Clayton)
    return (1 + x)^(-1 / c.ϑ)
end

function φ⁻¹(x::Real, c::Clayton)
    return x^(-c.ϑ) - 1
end

function φ²(x::Real, c::Clayton)
    α = 1 / c.ϑ
    return (α + α^2) * (1 + x)^(-(1 + α))
end

function ∇φ⁻¹(x::Real, c::Clayton)
    return -c.ϑ * x^(-c.ϑ - 1)
end

function inverse_rosenblatt(U::AbstractMatrix, c::Clayton)
    u1 = U[1, :]
    u2 = U[2, :]

    return [u1'; ((u1 .^ (-c.ϑ) .* (u2 .^ (-c.ϑ / (c.ϑ + 1)) .- 1) .+ 1) .^ (-1 / c.ϑ))']
end

function τ(c::Clayton)
    return c.ϑ / (c.ϑ + 2)
end
