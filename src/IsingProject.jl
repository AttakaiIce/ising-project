module IsingProject

using Random

export initialize_spins

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

end
