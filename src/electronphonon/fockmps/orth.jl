# 正交化 / 规范化（后端：FiniteMPSAlgorithms）
#
# GTEMPO 的 QR/SVD 与 TruncationScheme（Z2Tensors 体系，Grassmann 路径共用）
# 和 FMA 的对应类型是相互独立的；`_fmaorth`/`_fmatrunc` 把 GTEMPO 侧的算法
# 配置翻译为 FMA 的类型后，实际的 QR/SVD sweep 由 FMA 在 FockMPS 内层的
# CanonicalMPS payload 上完成。

_fmaorth(::QR) = FMA.QR()
_fmaorth(::SVD) = FMA.SVD()
_fmatrunc(::NoTruncation) = FMA.NoTruncation()
_fmatrunc(t::TruncationDimension) = FMA.TruncateDim(t.dim)
_fmatrunc(t::TruncateRelError) = FMA.TruncateRelError(t.ϵ)
_fmatrunc(t::TruncateDimCutoff) = FMA.TruncateDimCutoff(t.D, t.ϵ, t.add_back)
_fmaorthalg(alg::Orthogonalize) = FMA.Orthogonalize(_fmaorth(alg.orth), _fmatrunc(alg.trunc), alg.normalize, alg.verbosity)

TK.leftorth!(psi::FockMPS; alg::Orthogonalize = Orthogonalize()) = (FMA._leftorth!(psi.parent, _fmaorth(alg.orth), _fmatrunc(alg.trunc), alg.normalize, alg.verbosity); psi)
TK.rightorth!(psi::FockMPS; alg::Orthogonalize = Orthogonalize()) = (FMA._rightorth!(psi.parent, _fmaorth(alg.orth), _fmatrunc(alg.trunc), alg.normalize, alg.verbosity); psi)

function canonicalize!(psi::FockMPS; alg::Orthogonalize = Orthogonalize(trunc=DefaultITruncation, normalize=false))
	alg.normalize && @warn "canonicalize with renormalization not recommanded for FockMPS"
	FMA._canonicalize!(psi.parent; alg=_fmaorthalg(alg))
	return psi
end
