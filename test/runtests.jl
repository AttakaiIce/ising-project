if !("@stdlib" in LOAD_PATH)
    push!(LOAD_PATH, "@stdlib")
end

using Test
using IsingProject
using Random

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

@testset "IsingProject Metropolis Tests" begin
    @testset "delta_energy_flip matches energy difference" begin
        spins = [1 -1; -1 1]
        for i in 1:2, j in 1:2
            ΔE = delta_energy_flip(spins, i, j; J = 1, h = 0.2)
            spins2 = copy(spins)
            spins2[i, j] *= -1
            ΔE_from_energy = energy(spins2; J = 1, h = 0.2) - energy(spins; J = 1, h = 0.2)
            @test isapprox(ΔE_from_energy, ΔE; atol = 1e-12, rtol = 0)
        end
    end

    @testset "metropolis_step! accepts ΔE ≤ 0 without extra RNG draw" begin
        # 乱数で選ばれる(i,j)に対して ΔE≤0 になる配置を作り、必ず受容されることを確認する
        L = 3
        seed = 2026
        rng_probe = MersenneTwister(seed)
        i = rand(rng_probe, 1:L)
        j = rand(rng_probe, 1:L)

        spins0 = fill(-1, L, L)
        spins0[i, j] = 1 # 周囲が-1なので ΔE = 2*1*(1*(-4)+0) = -8
        @test delta_energy_flip(spins0, i, j) == -8

        spins = copy(spins0)
        rng = MersenneTwister(seed)
        accepted, ΔE = metropolis_step!(spins, 1.0; rng = rng)
        @test accepted
        @test ΔE == -8
        @test spins[i, j] == -1
    end

    @testset "metropolis_step! acceptance for ΔE > 0 is stochastic but reproducible with fixed RNG" begin
        spins0 = fill(1, 2, 2)
        beta = 0.1
        seed = 7

        rng_probe = MersenneTwister(seed)
        i = rand(rng_probe, 1:2)
        j = rand(rng_probe, 1:2)
        ΔE_expected = delta_energy_flip(spins0, i, j; J = 1, h = 0)
        @test ΔE_expected > 0
        u = rand(rng_probe)
        accept_expected = u < exp(-beta * ΔE_expected)

        spins = copy(spins0)
        rng = MersenneTwister(seed)
        accepted, ΔE = metropolis_step!(spins, beta; rng = rng)
        @test ΔE == ΔE_expected
        @test accepted == accept_expected
        if accepted
            @test spins[i, j] == -1
        else
            @test spins == spins0
        end
    end

    @testset "metropolis_sweep! accepts all steps at beta=0" begin
        spins = [1 -1; -1 1]
        rng = MersenneTwister(123)
        nacc = metropolis_sweep!(spins, 0.0; rng = rng)
        @test nacc == length(spins)
    end
end
