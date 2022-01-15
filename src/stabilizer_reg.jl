mutable struct StabilizerReg{Tzv, Tm} <: AbstractRegister{1}
    st::Stabilizer{Tzv, Tm}
end
# StabilizerReg(st::Stabilizer{Tzv, Tm}) where {Tzv, Tm} = StabilizerReg{Tzv, Tm}(st)

Base.copy(reg::StabilizerReg{Tzv, Tm}) where {Tzv, Tm} = StabilizerReg{Tzv, Tm}(copy(reg.st))
function Base.show(io::IO, reg::StabilizerReg)
    println(io, typeof(reg), " with $(Yao.nqubits(reg)) qubits:")
    show(io, reg.st)
end

zero_stabilizer_state(n::Integer) = StabilizerReg(one(Stabilizer, n))

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
    new_reg = copy(reg)
    traceout!(new_reg, locs)
    # TODO: delete locs
    return new_reg
end 
function Yao.reorder!(reg::StabilizerReg, orders)
    reg2 = copy(reg)
    for i = 1:length(orders), r = 1:length(reg.st)
        @show i, r
        @show QuantumClifford.getzbit(reg2.st, r, orders[i]), QuantumClifford.getzbit(reg2.st, r, i)
        QuantumClifford.setxbit(reg.st, r, orders[i], QuantumClifford.getxbit(reg2.st, r, i), orders[i]-i)
        QuantumClifford.setzbit(reg.st, r, orders[i], QuantumClifford.getzbit(reg2.st, r, i), orders[i]-i)
    end
    return reg
end
Yao.invorder!(reg) = Yao.reorder!(reg, nqubits(reg):-1:1)

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