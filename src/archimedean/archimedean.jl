abstract type ArchimedeanCopula <: AbstractCopula end

function rosenblatt(M::AbstractMatrix, c::ArchimedeanCopula)
    u1 = M[1, :]
    u2 = M[2, :]

    return [
        u1'
        (
            gen_derivative.(inverse_gen.(u1, c) + inverse_gen.(u2, c), c) ./
            gen_derivative.(inverse_gen.(u1, c), c)
        )'
    ]
end

function sample(c::ArchimedeanCopula, n::Int64)
    @assert n > 0
    return inverse_rosenblatt(rand(2, n), c)
end

include("clayton.jl")
