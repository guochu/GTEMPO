# singlemodebath — 两带 Anderson 模型验证（SecondOrderWeyl 完整版）

两带 AndersonIM（U=0.8, μ=0.7）+ 每带一个 Dirac-δ 单模浴（β=1, ω=1, α=0.5）。
算法模块 = `../src/SecondOrderWeyl.jl`；ED 参考 = `../../validation/ed.jl`
（`SingleModeED.observables` 返回逐带占据；总占据 = 2×逐带）。

## 运行

```bash
julia --project=<GTEMPO仓库> benchmark/secondorderif/singlemodebath/convergence.jl <仓库根>
python3 benchmark/secondorderif/plot_scaling.py   # 出图
```

## 结果（2026-11-05，`results/`）

| 观测量 | order 1 斜率 | order 2 斜率 | 最细误差 |
|--------|--------------|--------------|----------|
| 虚时 ⟨n₁+n₂⟩（N∈{4,8,16}） | 1.607 | **2.004** | 7.5e-5 (N=16) |
| 实时 ⟨n₁+n₂(0.4)⟩（N∈{8,16,32}） | 0.584 | **1.927** | 3.0e-6 (N=32) |

Anderson 实时的 U≠0 离散误差常数较大（err/h² → 1.9e-2，见
`anderson_realtime.csv`），纯二阶渐近区 N≥8；逐带 ⟨n_b⟩、等时两点 0⁺
的 ED 验证见 `../../validation/twopoint_band_check.jl`（四节全绿）。

图：`results/scaling_tau.png`、`results/scaling_realtime.png`。
数据：`results/anderson_tau.csv`、`results/anderson_realtime.csv`、
`results/anderson_slopes.csv`。

历史（2026-11-04/05 的部件级对齐、七轮修复史）见 `../HISTORY.md`；
旧单模特例实现已删除（本目录现在只含对 `../src` 完整版模块的验证）。
