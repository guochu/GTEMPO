# 奇算符插入 GF 缺陷诊断工具(2026-11-14/15 会话)

未决问题与完整诊断链见 src/measurement.jl 文件头 + changes.md 批次条目。
**裁决已定**:GM 引擎无罪(decisive.jl),port split_cell/from_core/pairize
之一为凶。下会话首步 = 双栈 cell-1 切分变量张量逐块 diff。

## 脚本清单(运行方式见各文件头注释)
- `decisive.jl` —— **决定性实验**(--project=U1GTEMPO2):参照栈奇插入 K' 链
  喂 GM.Zvalue,机器精度正确(0.66081−0.06945i)→ 引擎无罪。
- `u1g_instr.jl`(--project=U1GTEMPO2)—— 参照栈 builder 逐胞仪器化
  (nnz 18/144/72/36、prefix ζ1@cut7)。
- `port_faithful_u1.jl`(--project=GTEMPO)—— port 忠实复刻(:u1 整数荷,
  逐胞结构与参照一致)+ 迹对拍(~e^{-49} 垃圾)。
- `port_faithful.jl` —— :z2 版三模式(含 cut-0 折叠/prefix 机制复刻)。
- `op_calibrate.jl` —— golden gt/lt vs 算符路线标定(k≤N−1 diff ~1.4e-3;
  k=N 缺 a(N+1) 变量)。
- `ed_m2_accept.jl` / `ed_m2_convergence.jl` —— M=2 ED 验收与收敛归因
  (绝对 ω 约定铁律的出处)。

## 下会话首步操作
双栈各 dump cell-1 切分的两个变量张量全块(扇区键+值)逐块 diff:
port 侧 :u1 荷标签 Int 直接对应参照 ChargeCharge 首槽;
:u1 荷方向 = d:−q / ddag:+q(U1GTEMPO2 correlation_insertions 同款)。
