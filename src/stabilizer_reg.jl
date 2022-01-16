"""
    StabilizerReg{Tzv, Tm} <: AbstractRegister{1}

A data struct for stabilizer states in Yao. To initialize a `StabilizerReg`, 
one can use [`zero_stabilizer_state`](@ref) and [`rand_stabilizer_state`](@ref).
"""
mutable struct StabilizerReg{Tzv, Tm} <: AbstractRegister{1}
    st::Stabilizer{Tzv, Tm}
end

Base.copy(reg::StabilizerReg{Tzv, Tm}) where {Tzv, Tm} = StabilizerReg{Tzv, Tm}(copy(reg.st))
function Base.show(io::IO, reg::StabilizerReg)
    println(io, typeof(reg), " with $(Yao.nqubits(reg)) qubits:")
    show(io, reg.st)
end

"""
    zero_stabilizer_state(n::Integer) -> StabilizerReg

Initialize a stabilizer state |0...0> of n-qubits.
"""
zero_stabilizer_state(n::Integer) = StabilizerReg(one(Stabilizer, n))

"""
    rand_stabilizer_state(n::Integer) -> StabilizerReg

Initialize a random stabilizer state of n-qubits.
"""
rand_stabilizer_state(n::Integer) = StabilizerReg(random_stabilizer(n))

function Yao.addbits!(reg::StabilizerReg, n::Integer)
    reg.st = reg.st ⊗ one(Stabilizer, n)
    return reg
end
Yao.fidelity(r1::StabilizerReg, r2::StabilizerReg) = LinearAlgebra.dot(r1.st, r2.st)
Yao.instruct!(reg::StabilizerReg, op...) = throw(MethodError(Yao.instruct!, (reg, op...)))
Yao.nactive(reg::StabilizerReg) = Yao.nqubits(reg)
Yao.nbatch(reg::StabilizerReg{Tzv, Tm}) where {Tzv, Tm} = 1
Yao.nqubits(reg::StabilizerReg) = QuantumClifford.nqubits(reg.st)
Yao.nremain(reg::StabilizerReg) = 0
function Yao.partial_tr(reg::StabilizerReg, locs)
    n = Yao.nqubits(reg)
    qubit_map = sort!(setdiff(1:n, locs))
    new_st = zero(Stabilizer, length(qubit_map))
    rank = 0
    temp_st = traceout!(copy(reg.st), locs)
    for r = 1:length(temp_st)
        rank += 1
        for i = 1:length(qubit_map)
            if !iszero(temp_st[r])
                new_st[rank, i] = temp_st[r, qubit_map[i]]
            end
        end
    end
    return StabilizerReg(new_st)
end 
function Yao.reorder!(reg::StabilizerReg, orders)
    reg2 = copy(reg)
    for i = 1:length(orders), r = 1:length(reg.st)
        reg.st[r, orders[i]] = reg2.st[r, i]
    end
    return reg
end
Yao.invorder!(reg::StabilizerReg) = Yao.reorder!(reg, Yao.nqubits(reg):-1:1)
Yao.invorder(reg::StabilizerReg) = Yao.invorder!(copy(reg))

PauliString{N} = KronBlock{N, N, <:NTuple{N,Yao.YaoBlocks.PauliGate}}
"""
    PauliString(paulis...)

Generate a `KronBlock` to represent a Pauli operator.
"""
PauliString(args::YaoBlocks.PauliGate...) = kron(args...)

# focus!
# insert_qubits!
# collapseto!
# density_matrix
# probs
# purify
# relax!
# state
# statevec
# collapseto!
# density_matrix
# von_neumann_entropy
# ρ
# basis