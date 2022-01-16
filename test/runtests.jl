using YaoClifford
using Yao
using Test

@testset "GHZ states" begin
    n = 3
    reg1 = zero_stabilizer_state(n)
    @test measure!(reg1) == 0
    circ = chain(put(n, 1 => Yao.H), chain(cnot(i, i+1) for i = 1:n-1))
    Yao.apply!(reg1, circ)
    reg2 = zero_stabilizer_state(n)
    fidelity(reg1, reg2) ≈ sqrt(1/2)
    @test measure!(reg1) in (0b000, 0b111)
    @test fidelity(reg1, reg2) in (0, 1)
end

@testset "Reorder" begin
    reg = zero_stabilizer_state(2)
    reg |> put(2, 1 => Yao.X)
    @test measure!(invorder!(reg)) == 0b10
end

@testset "Partial trace" begin
    reg = zero_stabilizer_state(2)
    reg |> chain(put(2, 1 => H), cnot(1, 2))
    reg2 = partial_tr(reg, [1])
    @test expect(PauliString(X), reg2) == expect(PauliString(Y), reg2) == expect(PauliString(Z), reg2) == 0
end

@testset "Expectation of Hamiltonian" begin
    reg = zero_stabilizer_state(2)
    reg |> chain(put(2, 1 => H), cnot(1, 2))
    Hmts = [PauliString(X, X), PauliString(Z, Z), PauliString(Y, Y), -2im*PauliString(I2, I2)]
    @test [expect(Hmt, reg) for Hmt in Hmts] == [1, 1, -1, -2im]
    rand_reg = rand_stabilizer_state(2)
    @test sum(expect(Hmt, rand_reg) for Hmt in Hmts) == expect(sum(Hmts), rand_reg)
end