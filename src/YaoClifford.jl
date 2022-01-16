module YaoClifford

using Yao
using QuantumClifford
using YaoExtensions
using BitBasis
using LinearAlgebra

export StabilizerReg,
    zero_stabilizer_state,
    rand_stabilizer_state

include("stabilizer_reg.jl")
include("instruct.jl")
include("measure.jl")

end
