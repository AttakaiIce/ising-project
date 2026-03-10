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
