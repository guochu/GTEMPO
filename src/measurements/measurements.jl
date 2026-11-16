# measurements.jl：观测量层总纲（D5 两轨）。
#
# 两轨并行测量方式（D5 裁决，值差 O(胞内接触项)，各自独立成文件）：
#   - inserted_k.jl：inserted-K 边界语义——「带插入重建 K → GM 引擎
#     integrate → /Z」，边界时刻观测量（occupation/greater/lesser/
#     matsubara/nn/电流族）；golden 锚定路径（双栈对拍机器精度）。
#   - operators.jl：算符/胞坐标语义——GM ExpectationCache 快照 +
#     SparseGMPS 单项式（旧栈 cached_gf/TwosideExpectationCache 同构），
#     胞坐标 0⁺ 点的量，与边界插入差 O(胞内接触项)（cached_gf docstring
#     原文）。
#
# `GrassmannContourIndex` 在 lattices.jl（路线 2 消费）；`OperatorInsertion`
# 在 sysdynamics.jl（K 构造 API，测量层只消费）。

include("inserted_k.jl")
include("operators.jl")

# 公开皮肤
export occupation, greater, lesser, matsubara, default_branch, nn,
       raw_lead_correlator, electriccurrent, heatcurrent
export electriccurrent_fast, cached_electriccurrent_fast
