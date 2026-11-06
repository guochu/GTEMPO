# toulouse — 无相互作用 Toulouse 模型验证（SecondOrderWeyl 完整版）

单带（U=0）SIAM + Dirac-δ 单模浴（β=1, μ=0.7, ω=1, α=0.5）。
算法模块 = `../src/SecondOrderWeyl.jl`；ED 参考 = `../../secondorderif/validation/ed.jl`。

## 运行

```bash
julia --project=<GTEMPO仓库> benchmark/secondorderif/toulouse/convergence.jl <仓库根>
python3 benchmark/secondorderif/plot_scaling.py   # 出图
```

## 结果（2026-11-05，`results/`）

| 观测量 | order 1 斜率 | order 2 斜率 | 最细误差 |
|--------|--------------|--------------|----------|
| 虚时 ⟨n⟩（N∈{5,10,20,40}） | 1.380 | **2.002** | 5.9e-6 (N=40) |
| 实时 ⟨n(0.4)⟩（N∈{2,4,8}） | 0.822 | **1.996** | 1.9e-6 (N=8) |

图：`results/scaling_tau.png`、`results/scaling_realtime.png`
（log-log 误差-步长，order-2 线平行于 ∝h² 参考线）。

数据：`results/toulouse_tau.csv`、`results/toulouse_realtime.csv`、
`results/toulouse_slopes.csv`（case,order,N,h,value,reference,abs_err）。
