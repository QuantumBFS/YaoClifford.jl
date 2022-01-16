using QuantumClifford
using Yao
using LinearAlgebra
using Plots

N = 12
M = 1000
sts_clifford = [random_stabilizer(N) for _ = 1:M]
sts_haar = [rand_state(N) for _ = 1:M]
kernel_mat_clifford = ones(M, M)
for i = 1:M, j = i+1:M
    kernel_mat_clifford[i,j] = sts_clifford[i]⋅sts_clifford[j]
    @show i,j
end
for i = 1:M, j = i+1:M
    kernel_mat_clifford[j,i] = kernel_mat_clifford[i,j]
end

kernel_mat_haar = ones(M, M)
for i = 1:M, j = i+1:M
    kernel_mat_haar[i,j] = fidelity(sts_haar[i], sts_haar[j])
    @show i,j
end
for i = 1:M, j = i+1:M
    kernel_mat_haar[j,i] = kernel_mat_haar[i,j]
end

begin
    histogram(eigvals(kernel_mat_haar); label = "haar")
    histogram!(eigvals(kernel_mat_clifford); label = "clifford")
end