const single_qubit_ops = (
    (:X, :sX), (:Y, :sY), (:Z, :sZ), (:H, :sHadamard),
    (:S, :sPhase), (:Sdag, :sInvPhase),
    (:I2, :sId1)
)
for (sym, op) in single_qubit_ops
    @eval function Yao.instruct!(r::StabilizerReg, ::Val{$(QuoteNode(sym))}, locs)
        for l in locs
            QuantumClifford.apply!(r.st, QuantumClifford.$(op)(l))
        end
        return r
    end
end

function Yao.instruct!(r::StabilizerReg, ::Val{:SWAP}, locs)
    length(locs) == 2 || error("SWAP gate should be applied on 2 qubits")
    QuantumClifford.apply!(r.st, QuantumClifford.sSWAP(locs[1], locs[2]))
    return r
end

function Yao.instruct!(r::StabilizerReg, ::Val{:X}, locs, ctrl_locs, ctrl_configs)
    length(ctrl_locs) == 1 || error("Gates with multiple control qubits are not supported")
    for l in locs
        ctrl_configs !== (1,) && QuantumClifford.apply!(r.st, QuantumClifford.sX(ctrl_locs))
        QuantumClifford.apply!(r.st, QuantumClifford.sCNOT(ctrl_locs[1], l))
        ctrl_configs !== (1,) && QuantumClifford.apply!(r.st, QuantumClifford.sX(ctrl_locs))
    end
    return r
end
function Yao.instruct!(r::StabilizerReg, ::Val{:Z}, locs, ctrl_locs, ctrl_configs)
    Yao.instruct!(r, Val(:H), locs)
    Yao.instruct!(r, Val(:X), locs, ctrl_locs, ctrl_configs)
    Yao.instruct!(r, Val(:H), locs)
    return r
end