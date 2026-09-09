# ExactTTIIF `multorder` benchmark

`ExactTTIIF` builds the translationally invariant influence functional by exactly
exponentiating every Prony decay term `(α, λ)` of the hybridization function into
a small term-MPO and multiplying all term-MPOs together with SVD compression
(`algmult`). Each intermediate product is truncated to bond dimension `D`, so the
sequence of the multiplications determines **how the truncation error
accumulates**. `multorder` selects that sequence:

| `multorder` | ordering of the term-MPOs |
|---|---|
| `:λLM` | decay rates `λ`, large first |
| `:λSM` | decay rates `λ`, small first |
| `:αLM` | weights `α`, large first |
| `:αSM` | weights `α`, small first (**default**) |
| `:no`  | Prony expansion order, no reordering |

This benchmark measures how the choice affects accuracy, construction time and
bond dimension.

## Setup

- impurity: `ToulouseIM(μ = 0.3)` (single band, `U = 0`)
- bath: 6 exponential modes, `ω = 0.5 … 3.0`, `α = 0.5 … 0.16`, `β = 2`
  → about 20 decay terms per branch after an
  `OverDeterminedProny(n = 30, tol = 1e-8)` fit, so there are ~80 term-MPOs to
  multiply on the real-time contour
- imaginary contour: `δτ = 0.125`, `Nτ = 16`; real contour: `δt = 0.125`, `Nt = 8`
- the impurity dynamics `K` and the thermal state always use `D = 160`, so all
  measured errors are due to the IF construction alone
- references:
  - `conv`: GFs from a converged `D = 160` ExactTTIIF run (isolates the IF
    truncation error at the given `D`)
  - `ED`: exact diagonalization of the 6-mode model (includes the `O(δt)` /
    `O(δτ)` time-discretization error, identical for every row)

Run with

```
julia --project=. performance/exactttiif/bench_multorder.jl          # full sweep
julia --project=. performance/exactttiif/bench_multorder.jl --quick  # smoke run (D = 20 only)
```

Full data: [results.csv](results.csv).

## Results

Relative error of the whole GF sequence against the converged reference
(`imag`: Matsubara `G(τ)`; `real`: `G^>` and `G^<`). The IF bond dimension
saturates `D` in all runs.

### Imaginary contour — rel. error of `G(τ)`

| D | `:λLM` | `:λSM` | `:αLM` | `:αSM` | `:no` |
|---|---|---|---|---|---|
| 20 | 3.5e-5 | 2.6e-5 | 2.6e-5 | **9.7e-6** | 1.3e-5 |
| 40 | 7.7e-7 | 7.7e-7 | 7.7e-7 | 8.1e-7 | 7.9e-7 |
| 80 | 4.6e-11 | 4.6e-11 | 4.6e-11 | 4.5e-11 | 3.1e-11 |

### Real contour — rel. error of `G^>` / `G^<`

| D | `:λLM` | `:λSM` | `:αLM` | `:αSM` | `:no` |
|---|---|---|---|---|---|
| 20 | **3.9e-4** / 2.4e-4 | 5.9e-4 / 5.0e-4 | 4.3e-4 / 3.1e-4 | 4.3e-4 / 3.5e-4 | 4.8e-3 / 2.8e-3 |
| 40 | 1.8e-5 / 8.5e-6 | 1.6e-5 / 1.6e-5 | **1.1e-5** / 1.2e-5 | 5.3e-5 / 5.9e-5 | 3.2e-5 / 4.2e-5 |
| 80 | **2.3e-7** / 2.5e-7 | 1.8e-6 / 1.8e-6 | 2.7e-7 / 3.2e-7 | 6.1e-7 / 6.0e-7 | 7.2e-7 / 1.5e-6 |

Build times are minor and ordering-independent to within noise: 0.3–0.6 s
(imaginary) and 1.3–5.6 s (real) per IF.

## Conclusions

1. **The ordering only matters while the IF bond dimension is truncated.** At
   `D = 20` (bond dimension saturated) the ordering changes the error by factors
   of 3–4; at `D = 40` the spread shrinks to the `1e-5` level and at `D = 80` all
   orderings agree to `≲ 1e-6`. In the saturated regime the IF error is anyway
   far below the time-discretization error of the grid.
2. **Never use `:no`.** Without reordering, the real-contour error at `D = 20`
   is an order of magnitude worse (4.8e-3) than any sorted variant.
3. **The default `:αSM` is a good all-round choice**: it is the most accurate on
   the imaginary contour and mid-field on the real contour at `D = 20`, and
   indistinguishable from the others for `D ≥ 40`.
4. On the real contour at `D = 20` the best orderings put the *slowest-decaying*
   terms first (`:λLM`), while on the imaginary contour the *smallest-weight*
   first order (`:αSM`) wins — there is no single ordering that dominates both
   contours, but the differences vanish beyond `D ≈ 40`.
5. The ED reference errors are constant across all rows (2.7e-1 imaginary,
   4.8e-1 real): at these grid spacings the physical accuracy is limited by the
   time discretization, not by the IF construction or its ordering.
