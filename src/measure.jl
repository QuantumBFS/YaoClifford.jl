Yao.measure(reg::StabilizerReg; nshot::Integer = 1) = [measure!(copy(reg)) for _ = 1:nshot]
function Yao.measure!(reg::StabilizerReg)
    n = Yao.nqubits(reg)
    if n < 64
        meas_res = zero(Int64)
    elseif n < 128
        meas_res = zero(Int128)
    else
        meas_res = zero(BigInt)
    end
    idn = ['I' for _ = 1:n]
    for l in 1:n
        op_str = copy(idn)
        op_str[l] = 'Z'
        op = QuantumClifford._P_str(string(op_str...))
        _, anti_ind, proj_res = project!(reg.st, op)
        if proj_res === nothing
            proj_res = rand((0x00, 0x02))
            reg.st.phases[anti_ind] = proj_res
        end
        proj_res == 0x02 && (meas_res = flip(meas_res, one(meas_res) << (l-1)))
    end
    return BitStr{n, typeof(meas_res)}(meas_res)
end
