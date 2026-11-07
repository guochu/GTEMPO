# U1GTEMPO2 τ 两点对齐判决（第十二轮补遗，2026-11-07）

## 结论

**结构尾部不存在于 U1GTEMPO2**——金标准自身干净：

| 模型（双带 Anderson，Dirac-δ 浴 α=0.5，β=1，order 2） | U1GTEMPO2 误差 | SecondOrderWeyl 误差 |
|---|---|---|
| U=0, N=8, G(τᵢ,0) 全列 | 4.8e-5 ~ 1.47e-4 | 3.2e-4 ~ **2.9e-3** |
| U=0.8, N=8, G(τᵢ,0) 全列 | 2.0e-5 ~ 1.65e-4 | 3.0e-4 ~ **2.6e-3** |

U1 在 U=0 的双带列与单带离散化地板**逐位相同**（4.805e-5/2.290e-5/…/1.468e-4
= 我们 M=1 链的误差序列）——U1 管线对多带 τ 两点完全无结构尾部。

## 根因最终定位（修正 7782c0b 的「收缩引擎架构固有」结论）

两个框架**不共用收缩引擎**：U1GTEMPO2 是独立栈（直接 `import TensorKit`，
自有 GMPS/`berezin_vertex` 配对），本基准走 GTEMPO 的 Z2Tensors 栈。差异 =
**键电荷分辨力**：

- U1GTEMPO2：`SpinSector = U1Irrep ⊠ U1Irrep ⊕ FermionParity`——键携带
  每带独立粒子数电荷（`qsub(fqs[s],fqs[t])⊖prefix`），(1,0) 与 (0,1) 是
  不同扇区，分级融合逐带精确。
- 本基准（GTEMPO 栈）：单 Z2 宇称 `paritybits(s−1)⊕paritybits(t−1)⊕F`——
  (1,0)/(0,1) 合并进同一扇区，跨带配置的融合符号错标 ⟹ ~1e-3 bilinear
  尾部（K₁-slots × IF₂-selections，第十一轮定性不变）。

K 侧符号不可修的证明（7782c0b）仍然成立——正确的修复不是符号而是
**grading 升级**。GTEMPO 核心 Z2Irrep 硬编码 36 处（integration/observables/
grassmannmps/lattices/sysdynamics），泛化 = 核心重构；完整参考实现 =
U1GTEMPO2 的 SpinSector 方案（`tensors.jl` 的 `split_cell` + `from_core`
电荷标签拆分、`local_weight` 的 graded⊗ 符号、显式 oddlater）。

## 工具

`diag_u1_tau_twopoint.jl`——三框架同模型对比（ED / U1GTEMPO2 /
SecondOrderWeyl），U ∈ {0, 0.8} × N ∈ {4, 8} 的 G(τᵢ,0) 列逐点误差。
运行：

```
julia --project=~/Documents/Missile/U1GTEMPO2 \
    benchmark/secondorderif/validation/diag_u1_tau_twopoint.jl <GTEMPO repo>
```
