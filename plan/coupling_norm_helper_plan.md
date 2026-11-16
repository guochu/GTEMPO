# `coupling_norm(bath)` 助手函数计划：谱强度求和与连续谱数值积分

- 日期：2026-11-16
- 基线：GTEMPO HEAD `c01f446c` + 工作区
- 背景：GTEMPO 的 `HybridizationKernel.coupling_norm`（全带总耦合 Σ_b α_b）目前**只能由用户手写数值**（influence.jl:59 "由调用方传入"）；本计划提供 `coupling_norm` 计算
- 关联：`measurement_cache_refactor_plan.md` / `bridge_interface_simplification_plan.md` 独立并行，无文本冲突

## 0. 已裁决

- 连续 bath **直接数值积分**（用户 2026-11-16）；不做解析特例、不做缓存（一次性标量，积分代价可忽略）。
- **μ 不进 `coupling_norm`**（用户 2026-11-16）：Δ 装配里的 `spectrumshift(f0, μ)` 是核函数构造步骤（费米因子随 μ 变），谱权重本身与 μ 无关——离散谱求和天然 μ 不变；bath 层方法直接 `coupling_norm(bath.spectrum)`。
- **落点 = GTEMPO/src/correlationfunction.jl**（用户 2026-11-16）：与 `correlationfunction(bath, lattice)`（bath→QuAPI 包装，同文件 8-19 行）同族——都是"浴侧标量/关联函数助手"；单包改动，ImpurityModelBase 零改动。

## 1. 探查定谳（计划赖以成立的关键事实）

**F1【积分机器已按谱类型齐备——本计划零新积分代码】** `quadgkwrapper` 体系（ImpurityModelBase/spectrumfuncs/）对三类谱各有方法：
- `DeltaMultF`（单 δ 峰 × 函数）：`α·g(ω₀)`（prodfunc.jl:40-43）；
- `DeltasMultF`（DiscreteSpectrum × 函数）：`Σₙ αₙ·g(ωₙ)`（prodfunc.jl:81-87）——**fs = αₙ 的 δ-峰语义由此定谳**（非"线性插值密度"——discretespectrum.jl 头注释的插值描述与实际积分语义不符，以 prodfunc 积分规则为准）；
- 泛型 `AbstractBoundedFunction`：真数值积分 `_quadgk(m.f, lowerbound, upperbound)`（specfuncwrapper.jl:29，容差 1e-6）。

⟹ `coupling_norm(f) = quadgkwrapper(f × (ε→1))` 对三种谱同时给出 α / Σαₙ / ∫f dε。

**F2【`spectrumcouplings` 语义先例】** `spectrumcouplings(x::DiscreteSpectrum) = sqrt.(spectrumvalues(x))`（discretespectrum.jl:46-51，已导出）——"耦合强度 = √谱值"（vₙ = √αₙ，与 ED 约定 test/util.jl:159 `v = sqrt(α)` 一致）已是包内官方语义，`coupling_norm` 是其平方和的自然补全。

**F3【δ 谱不能直接积分】** `DiracDelta`/`DiscreteSpectrum` 没有活跃的 `quadgkwrapper` 直接方法（discretespectrum.jl:58-68 已注释、DiracDelta 同），且 `DiscreteSpectrum` 禁止按函数调用（:56 显式 error）——δ 谱必须经 `* (ε→1)` 乘法包装走 `DeltaMultF`/`DeltasMultF` 路径。边界：`DiracDelta` bounds = ±Inf（diracdeltaspectrum.jl:45-46），无出界截断问题。

**F4【落点文件现状】** `src/correlationfunction.jl`（43 行）= QuAPI wrapper：`correlationfunction(bath, lattice)` 三方法（Imag/Real/Mixed，:8-19）+ BCS 三方法（:38-43）+ 实时约定注记；GTEMPO.jl:51 include、:98 export。新方法同文件追加、同点导出。

**F5【GTEMPO 多带约定】** 每带一个 bath/kernel（`hybriddynamics(lat, bath, alg; band=b)`，Is... 变参每带一链）——Weyl 标量需要的总量 = Σ_b `coupling_norm(bath_b)`。

**F6【范围边界】** BCS 浴不立项（GTEMPO BCS 仅 :z2 一阶，coupling_norm 惰性——bcs_symmetry_guard；BCS 谱含 Δ 反常通道，norm 语义未定义，**不定义方法** = 天然拒绝）；玻色/费米粒子种类对谱强度无区别，不分派。

## 2. 设计（全部落 GTEMPO/src/correlationfunction.jl）

```julia
"""
    coupling_norm(f) -> Real
    coupling_norm(bath) -> Real

浴谱总强度 α（`HybridizationKernel.coupling_norm` = Σ_b α_b 的单浴因子）：
- `DiracDelta`：α；
- `DiscreteSpectrum`：Σₙ αₙ（fs 为 δ-峰权重，与 DeltasMultF 积分语义一致；
  vₙ = √αₙ = spectrumcouplings）；
- 连续 `BoundedFunction`：∫ f(ε)dε over [lb, ub]（quadgk 数值积分，容差同 QuadGKTolerance）；
- δ×函数乘积（DeltaMultF/DeltasMultF，如热流 ω 加权谱）：按各自积分规则取 g ≡ 1；
- `AbstractNormalBath`：直接取其谱——**coupling_norm 与 μ、β 无关**（它们只进
  Δ 装配的核因子，不进谱权重）。
"""
coupling_norm(f::DiracDelta) = f.α
coupling_norm(f::DiscreteSpectrum) = sum(spectrumvalues(f))
coupling_norm(f::BoundedFunction) = quadgkwrapper(f)                        # 连续谱：数值积分（F1 泛型方法现成）
coupling_norm(f::Union{DeltaMultF,DeltasMultF}) = quadgkwrapper(f, ϵ->1.0)  # 或 f*(ϵ->1.0) 包装
coupling_norm(bath::AbstractNormalBath) = coupling_norm(bath.spectrum)
```

- export 加 GTEMPO.jl:98 export 行。
- 直接公式法（DiracDelta/DiscreteSpectrum）与乘积包装法数值恒等（g≡1），双写是为避免"δ 谱不能直接积分"（F3）的陷阱暴露给调用方。
- 用户自定义谱子类：实现 `quadgkwrapper` 乘积路径或直接补 `coupling_norm` 方法——docstring 注明扩展点。
- `hybridization(lattice, bath; coupling_norm=nothing)` docstring 补一行：单带可 `coupling_norm=coupling_norm(bath)`；多带 = `sum(coupling_norm, baths)`。**不自动接线**——`nothing` 语义与 order=2 必填守卫（sysdynamics.jl:157）保持原样（自动缺省列为后续可选项，需另行裁决）。

## 3. 验收

1. **单元测试**（test/api/correlation.jl 追加）：
   - `DiracDelta(ω, α)` → α；`DiscreteSpectrum` → Σfs（含多点）；
   - 连续谱解析对照：矩形谱 ∫ = α·(ub−lb)、高斯谱 vs 解析值，quadgk 容差内；
   - `DeltaMultF` → α·g(ω₀) 抽查；`coupling_norm(bath)` 对 μ/β 变化不变（bath.spectrum 直取的回归锚）。
2. **端到端**：fewmodesbath realtime 测试把手写 `α=0.5` 换成 `coupling_norm(bath)`，传入值逐位相同（纯替换，测试预期不变）；golden 全套零改动（不接线即零语义变化）。
3. GTEMPO 全测试绿；ImpurityModelBase 零改动。

## 4. 不做清单

- 不自动接线 `hybridization`/`sysdynamics_weyl2order` 的缺省值（裁决面留给用户）；
- μ/β 依赖（N2 裁决）、BCS 浴方法（F6）、连续谱解析积分特例、结果缓存。

## 5. 工作量与顺序

约半天：correlationfunction.jl 实现 + docstring → 单元测试 + 端到端替换 → 全绿。单提交。
