# scaling.jl：缩放标量（后端直用）。
#
# 缩放标量 = **GrassmannMPS.ScaledScalar 直接使用**（2026-11-16：旧 ScaledComplex
# 包装退役——GM.ScaledScalar 为其严格超集：参数化实/复相位（实链结果保持实型）、
# 同构造归一化、value ≡ 旧 unscale（同 OverflowError 守卫）、完整 scaled 算术
# （×/÷/+/−/conj/abs，log 域精确）。公共 API：ScaledScalar/value 由主文件
# using+export 透传；unscale = value 的历史名别名；scaled_ratio = 精确比值物化。
#
# BoundedCache 已删（2026-11-16 A4 裁决：GTEMPO 全树零调用点）。
const unscale = GM.value
scaled_ratio(a::GM.ScaledScalar, b::GM.ScaledScalar) = GM.value(a / b)
