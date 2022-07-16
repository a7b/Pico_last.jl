module QuantumLogic

export GATES
export apply
export ket_to_iso
export iso_to_ket

using LinearAlgebra

const GATES = Dict(
    :X => [0 1;
           1 0],

    :Y => [0 -im;
           im 0],

    :Z => [1 0;
           0 -1],

    :H => [1 1;
           1 -1]/√2,

    :CX => [1 0 0 0;
            0 1 0 0;
            0 0 0 1;
            0 0 1 0],
)

function apply(gate::Symbol, ψ::Vector{T} where T<:Number)
    @assert norm(ψ) ≈ 1.0
    @assert gate in keys(GATES) "gate not found"
    U = GATES[gate]
    @assert size(U)[2] == size(ψ)[1] "gate size does not match ket dim"
    return ComplexF64.(normalize(U * ψ))
end

ket_to_iso(ψ) = [real(ψ); imag(ψ)]

iso_to_ket(ψ̃) = ψ̃[1:div(length(ψ̃), 2)] + im * ψ̃[(div(length(ψ̃), 2) + 1):end]

end