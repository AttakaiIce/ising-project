module IsingProject

using Random

export initialize_spins, magnetization, energy, delta_energy_flip
export metropolis_step!, metropolis_sweep!

"""
    initialize_spins(L::Int)

一辺の長さが `L` の2次元正方格子におけるイジングモデルの初期状態を生成します。
各スピンはランダムに +1 または -1 の値を持ちます。

# 引数
- `L::Int`: 格子の線形サイズ（全スピン数は L x L となります）

# 戻り値
- `Matrix{Int}`: +1 と -1 で構成される L x L の行列
"""
function initialize_spins(L::Int)::Matrix{Int}
    return rand([-1, 1], L, L)
end

"""
    magnetization(spins::AbstractMatrix{<:Integer}) -> Int

2次元イジング格子 `spins` の**全磁化** `M = sum(s_i)` を返します。
各スピン `s_i` は通常 ±1 を想定します。
"""
function magnetization(spins::AbstractMatrix{<:Integer})::Int
    return sum(spins)
end

"""
    energy(spins::AbstractMatrix{<:Integer}; J::Real=1, h::Real=0)

2次元イジングモデルのエネルギー

    E = -J * Σ_{⟨i,j⟩} sᵢ sⱼ - h * Σ_i sᵢ

を計算します。ここで Σ_{⟨i,j⟩} は最近接（上下左右）の結合の和です。
実装では**二重カウントを避けるため**、各格子点から「右」と「下」方向のみを足し上げます。
境界は**周期境界条件**（端は反対側に巻き戻る）です。
"""
function energy(spins::AbstractMatrix{<:Integer}; J::Real = 1, h::Real = 0)
    nrows::Int, ncols::Int = size(spins)
    T = promote_type(eltype(spins), typeof(J), typeof(h))
    e::T = zero(T)

    @inbounds for i in 1:nrows, j in 1:ncols
        s::T = spins[i, j]
        ip::Int = (i == nrows) ? 1 : (i + 1)
        jp::Int = (j == ncols) ? 1 : (j + 1)
        e -= T(J) * s * (spins[i, jp] + spins[ip, j])
        e -= T(h) * s
    end

    return e
end

"""
    delta_energy_flip(spins::AbstractMatrix{<:Integer}, i::Int, j::Int; J::Real=1, h::Real=0)

格子点 `(i, j)` のスピンを `s → -s` に**反転**したときのエネルギー差 `ΔE` を返します。
最近接（上下左右）相互作用と外場を

    E = -J * Σ_{⟨i,j⟩} sᵢ sⱼ - h * Σ_i sᵢ

で定義し、境界は周期境界条件（端は反対側に巻き戻る）とします。

このとき1点反転による差分は局所的に計算でき、

    ΔE = 2*s*(J*sum_neighbors + h)

（`sum_neighbors` は上下左右4近傍のスピン和）です。
"""
function delta_energy_flip(spins::AbstractMatrix{<:Integer}, i::Int, j::Int; J::Real = 1, h::Real = 0)
    nrows::Int, ncols::Int = size(spins)
    @boundscheck checkbounds(spins, i, j)

    T = promote_type(eltype(spins), typeof(J), typeof(h))
    s::T = spins[i, j]

    im::Int = (i == 1) ? nrows : (i - 1)
    ip::Int = (i == nrows) ? 1 : (i + 1)
    jm::Int = (j == 1) ? ncols : (j - 1)
    jp::Int = (j == ncols) ? 1 : (j + 1)

    @inbounds begin
        sum_neighbors::T = spins[im, j] + spins[ip, j] + spins[i, jm] + spins[i, jp]
        return 2 * s * (T(J) * sum_neighbors + T(h))
    end
end

"""
    metropolis_step!(spins::AbstractMatrix{<:Integer}, beta::Real; J::Real=1, h::Real=0, rng::AbstractRNG=Random.default_rng())

2次元イジング格子 `spins` に対して、メトロポリス法による**1回の更新試行**を行います。

- 1つの格子点 `(i, j)` をランダムに選ぶ
- 反転 `s → -s` を提案する
- エネルギー差 `ΔE` を `delta_energy_flip` で計算する
- `ΔE ≤ 0` なら受容、`ΔE > 0` なら確率 `exp(-beta*ΔE)` で受容

戻り値は `(accepted::Bool, ΔE)` です（棄却された場合は格子は変化しません）。
"""
function metropolis_step!(
    spins::AbstractMatrix{<:Integer},
    beta::Real;
    J::Real = 1,
    h::Real = 0,
    rng::AbstractRNG = Random.default_rng(),
)
    nrows::Int, ncols::Int = size(spins)
    i::Int = rand(rng, 1:nrows)
    j::Int = rand(rng, 1:ncols)

    ΔE = delta_energy_flip(spins, i, j; J = J, h = h)
    accepted::Bool = (ΔE <= 0) || (rand(rng) < exp(-beta * ΔE))

    if accepted
        @inbounds spins[i, j] = -spins[i, j]
    end

    return accepted, ΔE
end

"""
    metropolis_sweep!(spins::AbstractMatrix{<:Integer}, beta::Real; J::Real=1, h::Real=0, rng::AbstractRNG=Random.default_rng())

メトロポリス法の「**1スイープ**」を行います。
ここでは `size(spins,1) * size(spins,2)` 回だけ `metropolis_step!` を繰り返します。

戻り値は受容回数 `n_accepted` です。
"""
function metropolis_sweep!(
    spins::AbstractMatrix{<:Integer},
    beta::Real;
    J::Real = 1,
    h::Real = 0,
    rng::AbstractRNG = Random.default_rng(),
)::Int
    nsteps::Int = length(spins)
    n_accepted::Int = 0

    for _ in 1:nsteps
        accepted, _ = metropolis_step!(spins, beta; J = J, h = h, rng = rng)
        n_accepted += accepted ? 1 : 0
    end

    return n_accepted
end

end
