module IsingProject

using Random

export initialize_spins, magnetization, energy

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

end
