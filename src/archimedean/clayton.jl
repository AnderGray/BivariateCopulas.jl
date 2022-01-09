struct Clayton <: ArchimedeanCopula
    ϑ::Real

    function Clayton(θ::Real)
        # bivariate clayton can go as low as -1
        @assert -1 <= θ <= Inf

        return new(θ)
    end
end

function gen(x::Real, c::Clayton)
    return (1 + x)^(-1 / c.ϑ)
end

function inverse_gen(x::Real, c::Clayton)
    return x^(-c.ϑ) - 1
end

function gen_derivative(x::Real, c::Clayton)
    α = 1 / c.ϑ
    return gamma(1 + α) / gamma(α) * (1 + x)^(-(1 + α))
end

function inverse_rosenblatt(U::AbstractMatrix, c::Clayton)
    u1 = U[1, :]
    u2 = U[2, :]

    return [u1'; ((u1 .^ (-c.ϑ) .* (u2 .^ (-c.ϑ / (c.ϑ + 1)) .- 1) .+ 1) .^ (-1 / c.ϑ))']
end

function τ(c::Clayton)
    return c.ϑ / (c.ϑ + 2)
end
