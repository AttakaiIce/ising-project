if !("@stdlib" in LOAD_PATH)
    push!(LOAD_PATH, "@stdlib")
end

using Test
using IsingProject

@testset "IsingProject Initialization Tests" begin
    L = 10
    config = initialize_spins(L)

    # サイズの確認
    @test size(config) == (L, L)

    # 要素が +1 か -1 のみであることの確認
    @test all(s -> s in [-1, 1], config)

    # 要素の型が Int であることの確認
    @test eltype(config) == Int
end

@testset "IsingProject Magnetization Tests" begin
    @test magnetization([1 -1; 1 1]) == 2
    @test magnetization(fill(-1, 3, 4)) == -12
end

@testset "IsingProject Energy Tests" begin
    spins_all = fill(1, 2, 2)
    @test energy(spins_all) == -8
    @test energy(spins_all; h = 0.5) == -10.0

    spins_alt = [1 -1; -1 1]
    @test energy(spins_alt) == 8
end
