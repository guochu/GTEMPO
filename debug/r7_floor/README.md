# R7 排查档案：Weyl2 order=2 "~1e-5 δt 无关地板"（2026-11-19 定谳）

**结论：地板 = ED 参照初态 bug，不是 TN/order=2 装配问题。**
`validation/bcs_realtime_scaling.jl` 与 `test/bcs/runtests.jl` 的 ED 参照
`greater_lesser_ed(H, Ab, Bb, H0dec, ts, β)` 传了手搓的**未配对浴**热态
（`H0dec = kron(Himpbare,I4) + ω·kron(I4, n_bath)`），而 TN 的 kernel 编码
**配对浴**热态（bcs_Δt 的 u²/v²/uv 权重即源于配对浴热态）——ED 与 TN 的
物理模型不一致 ⟹ O(v²) 的 δt⁰ 系统差：
- Δ=0 时消失（v²=0）——臂 B（bcsbath Δ=0）与臂 C（normal bath）斜率 2.000 干净；
- 两阶共有（order=1 被更大 O(δt) 误差掩盖，且污染斜率拟合——旧 order=1
  slope 1.07–1.09，修正后 1.016/0.959）；
- 只入 greater/lesser 插入测量——Z 恒等式探针（U=1 下 Z=1 精确）全程无地板。

**修复**：`greater_lesser_ed` 传 bcs_ed 返回的 H0（杂质 ⊗ 配对浴）。
修正后 order=2 渐近斜率 = 2.000（4.3e-5 → 6.8e-7，Nt=4→32）。

## 探针脚本（按排查顺序）

1. `r7_floor_isolation.jl` —— 四臂隔离（BCS Δ=0.6 / BCS Δ=0 / normal 2带 /
   normal 1带）+ 配对浴初态修正后复跑（最终四臂全 2.000）。
2. `r7_z_identity.jl` —— Z 恒等式探针（U=1 Anderson 下 Keldysh 轨迹 Z=1
   精确成立）：Δ=0.6 order=2 |Z−1| 6.9e-8@Nt=32 斜率 2.000 无地板 ⟹ 地板
   只在插入测量。
3. `r7_lambda_bisect.jl` —— 反常块幅度 λ-二分：λ=0 地板仍在 ⟹ 排除反常块
   本身（地板绑定 Δ≠0 的 uu/dd Bogoliubov 双支结构）。
4. `r7_time_profile.jl` —— 逐时刻误差分解：order=2 固定时刻误差随 Nt
   收敛、晚时 δt⁰ 累积分量露出。
5. `r7_bdga_analytic.jl` / `r7_bdga2.jl` —— U=0 解析 BdG 第三参照
   （Nambu 8×8，G⁺(t)=U(t)·⟨ΞΞ†⟩ 结构）：裁决 ED vs TN 谁错——最终把
   ② 中解析-vs-ED 的随 t 增长偏差追到 ED 初态缺浴配对（探针自身的
   H0dec 也犯同样错，修正后解析 ≡ ED）。

## 排查教训

- **ED 参照初态必须与 IF 的 bath 定义一致**：BCS bath 的初态 = 配对浴热态
  （bcs_ed 返回的 H0），不能用"数守恒退耦"的未配对浴（旧约定注释
  「数守恒退耦 H0」是错的——TN 模型的 bath 从来不是数守恒的）。
- **Z 恒等式（Z=1）对任意线性耦合 H_S 精确成立**——是比 ED 更干净的
  收敛探针（无需 ED 参照、不受初态口径影响）。
- 单点斜率拟合会被 δt⁰ 失配污染（旧 order=1 slope 1.07–1.09 其实 = 1 被失配
  抬高）——收敛标定前先确认参照一致性。
