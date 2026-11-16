# 多轨道（bands=4）端到端验证分诊脚本（2026-11-16）

norb=2 KanamoriIM / bands=4 / 对角 1:1 混杂的验证与分诊脚本（下一轮修复的
工作资产；结论与完整线索见仓库 changes.md「D1 修订批次」及
GrassmannMPS/bugreports/2026-11-16-*.md）：

- `bare_k_atomic.jl`：裸 K 链迹 vs 原子极限 ED——:z2/:u1/:u1u1 机器精度 ✓；
  **:z2z2 在四对称性顺序循环的第 4 个阶段间歇 segfault**（`unknown function
  (ip: nil)`；单跑不复现，本脚本单跑 Nτ≤6 全部精确）。
- `norb2_z_validation.jl`：norb=2 :u1u1 全链验证——occ ✓ 1.7e-4 且 Nτ 收敛；
  绝对 Z 乘性抑制（缺 ~一份 Z_B = (1+e^{−βω})^bands，比值观测量免疫；
  Nτ 6→12 仅 4.73→4.57 近饱和 ⟹ 非逐格衰减，疑常数级归一化缺口）。
  NU/NJ/NT 环境变量可调。
- `merged_u1_discriminant.jl`：:u1 合并形 2 链同现象（5.6×）⟹ 排除链数。
- `z2z2_segfault_bisect.jl`：:z2z2 分阶段定位（lattice/vacuumstate/
  sysdynamics/integrate；NN 环境变量）。
- `alpha_scaling_diag.jl`：单带 × 耦合强度/浴模数缩放（C1/C2/D：Σ 核 ≡ 等效
  单模强耦合，逐位一致）。
