# 新旧栈 alignment sweep 报告（2026-11-16，Phase 6 验收项）

对照 = backup/GTEMPO（旧 Z2Tensors 栈，`--project=../backup/GTEMPO`）vs 本栈
（GrassmannMPS 后端）。场景（双栈同参数）：

| 场景 | 物理量 | 结果 |
|---|---|---|
| S1 虚时单带 δ-bath（N=4, δτ=0.2, μ=0.4） | Z | **Δ = 4.2e-15** |
| S1 occupation(τ, j=1..4) | 4 点 | **Δ = 2.8e-16** |
| S3 双带 U=1（N=4, δτ=0.2） | Z | **Δ = 9.0e-14** |
| S3 occupation(band1, j=1..4) | 4 点 | **Δ = 1.1e-15** |
| S2 实时单带热初态（N=4, δt=0.05） | Z | Δ = 1.9e-3（热态构造差：旧 systhermalstate! vs 新 initial_state 截断档） |
| S2 greater(k, 0) | 4+1 点 | **逐点 1e-5 吻合，差一个节点约定**（见下） |

**节点约定差（S2）**：旧栈格点 k = N+1（含 τ/t=0 端点），greater 第 1 点 =
τ=0⁺ 接触值 ≈ 0；新栈 k = N（节点 = δt..Nδt 边界），greater 第 1 点 =
旧栈第 2 点。移位后逐点吻合 ≤ 1e-5（两侧各自 IF 截断档）。

**已知预期分歧（在案，非回归）**：旧栈 ExactTTI 的 Prony ~5% 缺陷（新实现
8.7e-4 距 ED，fixtures check_golden_exactttiif）；M=2 带互换 bug（本 snapshot
的旧栈 S3 未触发，修复历史见 changes.md 2026-11-14 批次）。

复跑：
```
julia --project=. validation/alignment/run_new.jl records/align_new.json
julia --project=../backup/GTEMPO validation/alignment/run_old.jl records/align_old.json
```
（对照脚本内嵌于本报告数值表；数据 = records/*.json，Serialization 格式）
