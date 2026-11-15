> ⚠️ **API 迁移声明（grassmannmps-backend 分支，2026-11-15）**：本分支已换轨
> GrassmannMPS 后端——Grassmann ordering 家族已删（仅存 Chrono 三序）、
> `sysdynamics/hybriddynamics/greater/lesser/matsubara` 等签名全部更新
> （`symmetry` 手选、`coupling_norm` 随 K、插入显式带荷）。本文档与
> docs/tutorials 下的脚本以旧栈 API 写成，仅作旧接口参考；新 API 逐项
> 对照见仓库根 [changes.md](../../changes.md)，可运行的新 API 用法见
> `test/api/`（physics_battery/port_smoke/observables 即用法示例）。

# GTEMPO.jl

**GTEMPO** (Grassmann TEMPO) solves quantum impurity problems by representing
the influence functional (IF) of the bath as a matrix product state over
*Grassmann-valued* tensors. Compared with conventional TEMPO-style algorithms,
the Grassmann formulation removes the need for a huge local Hilbert space and
enables exact, compact representations of the IF on arbitrary time-contour
orderings.

## How it works

A calculation is organized around a few building blocks:

1. **Bath** — define a fermionic or bosonic bath with a spectral density
   (`fermionicbath`, `bosonicbath`).
2. **Lattice** — discretize the time contour into a `GrassmannLattice`
   (`contour = :imag`, `:real` or `:mixed`), possibly with several bands.
3. **Correlation function** — `correlationfunction(bath, lattice)` computes
   the discrete hybridization matrices on the chosen contour.
4. **Influence functional** — `hybriddynamics` builds the Grassmann MPS `I`
   of the IF, with four interchangeable construction algorithms
   (`PartialIF`, `XTRGIF`, `ExactTTIIF`, `TDVPIF`).
5. **Impurity dynamics** — `sysdynamics` builds the impurity propagator `K`;
   `boundarycondition!` closes the contour, `systhermalstate!` prepares a
   thermal initial state.
6. **Observables** — evaluate Green's functions (`gf`, `greater`, `lesser`),
   occupations (`occupation`, `cached_occupation`) and currents
   (`electriccurrent`, `electriccurrent_fast`, ...) either directly or
   through precomputed environment caches (`environments`).

A minimal real-time example:

```julia
using GTEMPO

bath  = fermionicbath(semicircular(5.0), β = 1.0, μ = 0.0)
model = ToulouseIM(ϵ_d = 1.25π)
lattice = GrassmannLattice(N = 10, δt = 0.05, contour = :real)
trunc = truncdimcutoff(D = 100, ϵ = 1e-10)

corr = correlationfunction(bath, lattice)
I    = hybriddynamics(lattice, corr, trunc = trunc)
K    = sysdynamics(lattice, model, trunc = trunc)
K    = boundarycondition!(K, lattice)

cache = environments(lattice, K, I)
gt = [-im * cached_greater(lattice, k, K, I; cache = cache) for k in 1:lattice.k]
```

## Contents

```@contents
Pages = ["manual.md", "practice.md", "grassmann_lattice.md", "tensor_and_grassmann_conventions.md", "internals.md", "tutorials.md", "api.md"]
Depth = 2
```
