# 链级线性代数（后端：FiniteMPSAlgorithms）
#
# dot/norm/lmul! 与 +/-/⊙(Hadamard 积) 全部委托给 FMA 在 CanonicalMPS payload
# 上的同名函数（含 per-site scaling 约定，总 scaling = scaling^L）。

LinearAlgebra.dot(psiA::FockMPS, psiB::FockMPS) = LinearAlgebra.dot(psiA.parent, psiB.parent)
LinearAlgebra.norm(psi::FockMPS) = LinearAlgebra.norm(psi.parent)
distance(a::FockMPS, b::FockMPS) = _distance(a, b)
distance2(a::FockMPS, b::FockMPS) = _distance2(a, b)

LinearAlgebra.lmul!(f::Number, psi::FockMPS) = (LinearAlgebra.lmul!(f, psi.parent); psi)

Base.:*(psi::FockMPS, f::Number) = FockMPS(psi.parent * f)
Base.:*(f::Number, psi::FockMPS) = psi * f
Base.:/(psi::FockMPS, f::Number) = psi * (1/f)
Base.:(-)(psi::FockMPS) = FockMPS(-psi.parent)

# 精确（未压缩）乘积：物理指标共享的逐点（Hadamard）乘积，委托给 FMA 的 ⊙
function Base.:*(x::FockMPS, y::FockMPS)
	(length(x) == length(y)) || throw(DimensionMismatch())
	return FockMPS(⊙(x.parent, y.parent))
end

# 块对角直和：FMA 的实现会把两边的 scaling 折入数据
function Base.:+(x::FockMPS, y::FockMPS)
	(length(x) == length(y)) || throw(DimensionMismatch())
	return FockMPS(x.parent + y.parent)
end
Base.:-(x::FockMPS, y::FockMPS) = x + (-y)

# 站点内容重排：委托给 FMA 的链级 permute!（内部经 canonicalize! 与相邻 swap!）
function TK.permute!(x::FockMPS, perm::Vector{Int}; trunc::TruncationScheme=DefaultKTruncation)
	FMA.permute!(x.parent, perm; trunc=_fmatrunc(trunc))
	return x
end
TK.permute(x::FockMPS, perm::Vector{Int}; kwargs...) = permute!(deepcopy(x), perm; kwargs...)
