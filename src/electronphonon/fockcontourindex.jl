"""
	FockContourIndex(j, branch, band)

electron-phonon（Fock 空间）partial IF 的一条腿：时间步 `j`、contour 分支
`branch`（`:τ`、`:+` 或 `:-`）、以及 band（声子模式）编号 `band`。

与 `GrassmannContourIndex` / `BosonicContourIndex` 一样是
[`QuAPI.AbstractContourIndex`](@ref) 的子类型，因此 `branch(x)` 直接可用；
它比 `BosonicContourIndex` 多带一个 `band` 字段，使 partial IF 可以按
(node, band) 逐条构造，而不必对整条腿合并所有 band。
"""
struct FockContourIndex <: AbstractContourIndex
	j::Int
	branch::Symbol
	band::Int
end

"""
	FockContourIndex(j; branch=:τ, band=1)

关键字构造器；默认对应虚时 `:τ` 分支、band 1。
"""
FockContourIndex(j::Int; branch::Symbol=:τ, band::Int=1) = FockContourIndex(j, branch, band)

Base.:(==)(a::FockContourIndex, b::FockContourIndex) =
	a.j == b.j && a.branch == b.branch && a.band == b.band

"""
	x[i::FockContourIndex, j::FockContourIndex]

读取（玻色）hybridization function 的矩阵元：
`branch(x; b1=branch(i), b2=branch(j))[i.j, j.j]`。普通声子 bath 与 band
无关，故 `band` 字段不参与取值。
"""
function Base.getindex(x::AbstractCorrelationFunction, i::FockContourIndex, j::FockContourIndex)
	return branch(x; b1=branch(i), b2=branch(j))[i.j, j.j]
end
