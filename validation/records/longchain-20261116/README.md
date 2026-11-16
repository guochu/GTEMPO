# longchain 性能基线（2026-11-16，GrassmannMPS 后端分支）

`gtempo_new.csv`：kind(real/imag) × order(1/2) × 胞数 n ∈ {8,16,32} × reps=3，
bands=2、symmetry=:u1、δ-bath (ω=1, α=0.5)、β=1、μ=0.7、real δt=0.15 / imag δτ=0.125。
stage = K（sysdynamics[_weyl2order]）/ IF（hybriddynamics·PartialIF 逐带）/
integrate（integrate_scaled）/ occupation / greater（均插入重建精确路线）。
中位数（秒）：

| n  | K(o2,real) | IF(real,o2) | integrate(real,o2) | occ(real,o2) | IF(imag,o2) | integrate(imag,o2) |
|----|------------|-------------|--------------------|--------------|-------------|--------------------|
| 8  | 0.076      | 1.081       | 0.097              | 0.212        | 0.245       | 0.050              |
| 16 | 0.106      | 4.438       | 0.209              | 0.320        | 1.034       | 0.099              |
| 32 | 0.212      | 17.522      | 0.411              | 0.676        | 4.358       | 0.202              |

标度：IF 构建每胞加倍 ~4.3×（行因子逐带链积，O(n²·bond) 主导）；K 构建、
integrate、观测量均近线性（n=32 全部亚秒）。与 U1GTEMPO2 基线
（validation/records/port_new_longchain.csv：imag n=16 IF ≈ 1.02s、
full_green ≈ 0.14–0.55s）同量级；stage 语义两栈不同构（U1GTEMPO2 的
structured_kernel/spin_relabel 等为本栈所无），故本表为**本栈纵向基线**
（未来回归对照），跨栈对照按 IF/full_green 量级参考。

复跑：`julia --project=. validation/longchain_benchmark.jl [outdir]`
