# Tensor and Grassmann conventions

This document explains how GTEMPO handles (fermionic) tensors and the
mathematical principles behind the two contraction interfaces, `@tensor`
(bosonic) and `@grassmann` (fermionic). Every rule is stated as a formula,
then demonstrated by a small piece of code whose **actual output is shown**.
Paste any code block into a Julia session that has GTEMPO loaded
(`using GTEMPO, Z2Tensors, TensorOperations`) and you should see the same
numbers.

---

## 0. Notation used in this document

| symbol | meaning |
|---|---|
| $\xi_\alpha$ | one Grassmann (G-) number on leg $\alpha$ |
| $p_\alpha \in \{0,1\}$ | Z2 sector ("parity") of leg $\alpha$: $0$ even, $1$ odd |
| $V$ | $Z2Space(0{\Rightarrow}1, 1{\Rightarrow}1)$: one even + one odd sector |
| $\bar V = V'$ | the dual space (the conjugate/bra space) |
| `isdual(space(t, i))` | how the code tells a bra leg ($\bar\xi$, dual) from a ket leg ($\xi$) |

The Z2-graded tensor space `Z2Space(0=>1, 1=>1)` is the smallest space that
contains a fermionic sector. All examples below use it (or its mixed cousin
`Z2Space(0=>1, 1=>2)`), so that every block is at most a few numbers and
everything is easy to check by hand.

---

## 1. The mathematical object: Grassmann tensors

The fundamental object handled by this package is a **Grassmann tensor**
(GT): a multi-dimensional coefficient array multiplied by a basis of
Grassmann (G-) numbers. A three-leg GT reads

$$
X_{ijk} \, \xi_1^{\,i} \xi_2^{\,j} \xi_3^{\,k} ,
\qquad i,j,k \in \{0,1\},
$$

where the string $\xi_1^{\,i}\xi_2^{\,j}\xi_3^{\,k}$ is the *basis element* of
the G-number algebra and $X_{ijk}$ the *coefficient*. G-numbers anticommute
and are nilpotent:

$$
\xi_\alpha \xi_\beta = -\xi_\beta \xi_\alpha
\quad (\alpha \neq \beta),\qquad
(\xi_\alpha)^2 = 0 ,\qquad
\xi_\alpha^{0} = 1 . \tag{1}
$$

Only the coefficient array $X$ is stored; the G-number string is *implicit in
the leg ordering*. Any algebraic operation must nevertheless keep track of
two G-number effects:

* the **exchange sign** $-1$ picked up whenever two odd G-numbers are
  permuted past each other (Eq. (1)); and
* the **twist** $-1$ associated with reordering a contracted pair out of its
  canonical order, or with closing a fermionic loop.

## 2. Storage: bosonic Z2-graded tensors (Z2Tensors)

The storage layer (Z2Tensors) is a bosonic Z2-graded tensor library: a
`TensorMap` is a block-sparse array over fusion trees of `Z2Irrep` sectors,
with the category structure of `Vec_Z2` (symmetric braiding with **trivial**
R-symbol and trivial twist). No Grassmann signs are applied at this level.
All sign handling lives in the `GrassmannBackend` described below.

Space duality is bookkept by a `dual` flag on `Z2Space`. The function
`space(t, i)` performs a "getindex dualization" so that the arrow
(ket/bra) character of each leg is exposed; `isdual(space(t, i))` is the
practical test used throughout GTEMPO. Example 1 shows it.

**Example 1 — sectors and dual flags.**

```julia
using GTEMPO, Z2Tensors

V    = Z2Space(0=>1, 1=>1)            # one even + one odd (fermionic) sector
Vbar = dual(V)
A = zeros(ComplexF64, V⊗V, Vbar⊗Vbar) # a (2,2) tensor:  codomain V⊗V, domain Vbar⊗Vbar
B = zeros(ComplexF64, Vbar⊗Vbar, V)   # a (2,1) tensor
println("A: isdual codomain = ", [isdual(space(A,i)) for i in 1:2],
        "   isdual domain = ", [isdual(space(A,i)) for i in 3:4])
println("B: isdual codomain = ", [isdual(space(B,i)) for i in 1:2],
        "   isdual domain = ", [isdual(space(B,i)) for i in 3:3])
```

Output:

```
A: isdual codomain = Bool[0, 0]   isdual domain = Bool[0, 0]
B: isdual codomain = Bool[1, 1]   isdual domain = Bool[1]
```

`A` and `B` are the two tensors used in the junction example below: the
contracted legs of `A` (its domain, a *ket*-flagged leg) meet the
corresponding legs of `B` (its codomain, a *bra*-flagged leg).

## 3. The two understandings of a Grassmann MPS

A GrassmannMPS (GMPS) represents a **ket**. It admits two equivalent
readings:

1. **GT reading.** Every site tensor is a GT; the MPS is a product of
   G-number strings with fused coefficient tensors. Any *ket–ket operation*
   (MPS multiplication, environment contraction, time evolution of kets)
   must be performed in this picture, with all fermionic signs.
2. **Coefficient reading.** The GMPS is just the coefficient array stored as
   an ordinary Z2-symmetric MPS, with no explicit G-numbers.

In the coefficient reading the inner product of a bra and a ket is an
ordinary coefficient pairing. Writing the ket coefficients as $\psi$ and the
bra coefficients as the conjugates $\bar\psi$,

$$
\langle \varphi | \psi \rangle \;=\; \sum_{i_1\cdots i_n}
\bar\varphi_{i_1\cdots i_n}\, \psi_{i_1\cdots i_n},
\tag{2}
$$

no G-number sign appears, because each pair is a coefficient product
$(\bar\xi^{\,i}\leftrightarrow \xi^{\,i})$, not a reordering of G-numbers.

GTEMPO deliberately switches between the two readings:

* **ket–ket operations** use reading 1 and the `@grassmann` interface
  (fermionic signs applied automatically).
* **bra–ket inner products** (e.g. $\langle\varphi|\psi\rangle$, observable
  sandwiches, the SVD re-combination $u \cdot v^\dagger$ inside `svdmult`)
  use reading 2: the contraction is an ordinary coefficient pairing and is
  **bosonic** — it must be written with `@tensor`, *not* `@grassmann`.
* **influence-functional construction** (`partialif`, `ttiif`) builds the
  equivalent bosonic MPO first, with the required Jordan–Wigner signs
  inserted manually into the MPO, and only then multiplies it onto the
  vacuum state. These are coefficient-world operations as well.

The practical rule of thumb, verified by the full test suite (§8):

> **`conj(...)` inside a contraction marks a ket–bra (coefficient)
> contraction and must therefore be bosonic.** A contraction like
> `data * conj(left)` is a *coefficient* contraction even when the tensors
> themselves carry odd sectors; only `@tensor` is correct there.
>
> ```julia
> # ket-bra (coefficient) contraction: bosonic, no fermionic twist
> @tensor u[-1 -2; -3] = twositemps[-1,-2,1,2] * conj(v[-3,1,2])
> ```
>
> As of 2026-09-04 there are no `@grassmann` expressions containing
> `conj(...)` anywhere in `src/`; new code must keep it that way.

## 4. `@tensor`: the bosonic interface

`@tensor` is the plain TensorOperations macro acting on Z2Tensors. Its
semantics are those of the bosonic category `Vec_Z2`: index permutations and
contractions never generate signs, no matter the sector parities. It is the
correct interface for

* the coefficient reading of GMPS (bra–ket contractions),
* the influence-functional MPO algebra (with JW signs inserted by hand),
* any operation on genuinely bosonic tensors.

Because even-parity sectors never carry fermionic signs, on tensors living
entirely in the even sector the two interfaces are identical. This is
Example 6 at the end of §5.

## 5. `@grassmann`: the fermionic interface

`@grassmann` accepts exactly the same expressions and keyword arguments as
`@tensor`, but injects the `GrassmannBackend` into every generated
`tensoradd!` / `tensortrace!` / `tensorcontract!` call:

```julia
# syntax only: A, B, C below stand for any conforming tensors
@grassmann C[a,b;c] := A[a,b,x,y] * B[x,y,c]   # same syntax as @tensor
@grassmann s = A[1,2] * B[2,1]                   # scalar output
@grassmann C[a,b;c] += A[a,b,x,y] * B[x,y,c]    # in-place addition
```

It is the correct interface for ket–ket operations on GTs. The following
subsections state, one by one, the sign rules it implements. The rules are
*fixed mathematical conventions* (identical to TensorKit's `FermionParity`
tensors and to the GrassmannTensors package), not heuristics; each is
demonstrated by a tiny runnable example.

### 5.1 G-string convention and the meaning of `isdual`

For a tensor $T: V_1\cdots V_m \leftarrow W_1\cdots W_n$ the canonical
G-number string is

$$
\xi(V_1)\cdots \xi(V_m)\; \bar\xi(W_n)\cdots \bar\xi(W_1),
\tag{3}
$$

i.e. codomain legs read left-to-right and domain legs read right-to-left.
The G-type of a leg is read off its space flag: a leg is a ket
$\xi$ ("$a$") when `isdual(space(t, i)) == false` and a bra $\bar\xi$
("$\bar a$") when it is `true` (Example 1 prints exactly these flags).

### 5.2 Index permutation: fermionic `gpermute`

Reordering the legs of a tensor reorders its G-string. Only the exchange of
two *odd* G-numbers costs a sign (Eq. (1)); moving an even sector
($\xi^{0}=1$) costs nothing. Hence the fermionic permutation sign is

$$
\operatorname{coef}(\pi) = \prod_{\substack{\text{adjacent swaps}\\
(i,\ i+1)\ \text{of } \pi}} (-1)^{p_i\, p_{i+1}},
\qquad p_i,p_{i+1}\in\{0,1\},
\tag{4}
$$

i.e. $(-1)$ for each odd–odd swap and $+1$ otherwise. This is the internal
`gpermute` (renamed from `f_permute`; not exported — use the `@grassmann`
macro, which applies these signs to every index permutation); the
bosonic `permute` never carries these signs. Example 2 isolates the effect.

**Example 2 — swapping two legs of a rank-(2,0) GT.**

```julia
using GTEMPO, Z2Tensors

V = Z2Space(0=>1, 1=>1)
T = zeros(ComplexF64, V⊗V)            # a rank-(2,0) Grassmann tensor
for (_, b) in blocks(T); b .= 1 + 0.7im; end
Tb = permute(T, (2,1), ())            # bosonic  permute
@grassmann Tf[2 1;] := T[1,2]         # fermionic permute (via the macro)
for (f1, f2) in fusiontrees(T)
    println("block ", f1.uncoupled, ":  T = ", vec(T[f1,f2]),
            "   perm^bos = ", vec(Tb[f1,f2]),
            "   perm^fer = ", vec(Tf[f1,f2]))
end
```

Output:

```
block (Z2Irrep(0), Z2Irrep(0)):  T = ComplexF64[1.0 + 0.7im]   perm^bos = ComplexF64[1.0 + 0.7im]   perm^fer = ComplexF64[1.0 + 0.7im]
block (Z2Irrep(1), Z2Irrep(1)):  T = ComplexF64[1.0 + 0.7im]   perm^bos = ComplexF64[1.0 + 0.7im]   perm^fer = ComplexF64[-1.0 - 0.7im]
```

Swapping the two legs rewrites the string $\xi_1\xi_2 \to \xi_2\xi_1 =
-\xi_1\xi_2$. The (even, even) block is untouched ($1 = \xi_1^{0}\xi_2^{0}$
permutes trivially), while the (odd, odd) block $\xi_1\xi_2$ picks up
$-1$: exactly Eq. (4) with $p_1 p_2 = 1$.

### 5.3 Contraction junctions: the canonical pair $(a,\bar a)$ and the twist

The contraction convention fixes the canonical order of a contracted pair as
$(a, \bar a)$: the *A*-side leg contributes the ket $a$, the *B*-side leg the
bra $\bar a$. The plain "chain" junction (A's domain region against B's
codomain region) realizes this order directly and needs no sign.

Whenever a contracted pair appears in the crossed order $(\bar a, a)$ — the
*A*-side leg carries a non-dual (ket) space against a dual (bra) space on
the B side, as always happens at the closure of a ring-like network — the two
G-numbers must first be permuted into the canonical order. That permutation
costs one fermionic swap, i.e. a **twist of $-1$ per odd sector** on the
pair, applied on the A side:

$$
\theta_{\mathrm{junction}} = (-1)^{\, \sum_k p_k},
\qquad \text{over the contracted A-side legs that are kets } a .
\tag{5}
$$

The rule is *unconditional*: it depends only on the G-types of the
contracted pair, not on which region the B-side legs live in. In the code
this is exactly the last step before `mul!` in `_contract!`
(`src/grassmanntensor/tensoroperations.jl`):

```julia
No = length(oindA)
inds = Tuple(No + k for k in eachindex(cindA) if !isdual(space(A, cindA[k])))
g_twist!(A′, inds)
```

`g_twist!` multiplies every fusion-tree block of `A′` by $(-1)^{\#\text{odd
sectors among the twisted legs}}$, i.e. implements Eq. (5).

**Example 3 — a two-leg junction.**  Fill *every* block of the two tensors
of Example 1 with $1 + 0.5\mathrm{i}$ and contract.

```julia
using GTEMPO, Z2Tensors

V = Z2Space(0=>1, 1=>1)
Vbar = dual(V)
A = zeros(ComplexF64, V⊗V, Vbar⊗Vbar)
B = zeros(ComplexF64, Vbar⊗Vbar, V)
for (_, b) in blocks(A); b .= 1 + 0.5im; end
for (_, b) in blocks(B); b .= 1 + 0.5im; end
@grassmann Rf[x,y;c] := A[x,y,a,b] * B[a,b,c]
@tensor    Rb[x,y;c] := A[x,y,a,b] * B[a,b,c]
println("first block entry of fermionic Rf = ", first(blocks(Rf))[2][1])
println("Rf == Rb (this open graph)        = ", Rf == Rb)
```

Output:

```
first block entry of fermionic Rf = 1.5 + 2.0im
Rf == Rb (this open graph)        = true
```

This is an *open* network: every fermion line starts and ends at the
boundary, so no closed loop exists and — after the fermionic bookkeeping —
the two interfaces agree on the same number $1.5 + 2\mathrm{i}$. The
fermionic convention only differs where a loop closes; see Examples 4 and 5.
(The identical number $1.5+2\mathrm{i}$ is produced by the independent
GrassmannTensors package on the same tensors, §7.)

### 5.4 Trace: the U-turn twist

A trace closes a loop through the planar "U-turn" between the ket and the
bra string. Closing the loop, **every traced codomain leg with a non-dual
space and an odd sector contributes a factor $-1$**, the topological twist
$\theta_f = -1$ of the fermionic sector:

$$
\theta_{\mathrm{trace}} = (-1)^{\, \#\{\text{traced non-dual, odd legs}\}} .
\tag{6}
$$

This is TensorKit's `_trace_permute!` rule — flip the sign once for every
traced subtree leg whose stored flag is not dual and whose sector is odd,
i.e. `θ = twist(uncoupled[i]) = -1` for odd sectors — translated to the
Z2Tensors tree convention, where the dual flags are recovered from the
spaces instead of stored on the fusion tree
(`src/grassmanntensor/tensoroperations.jl`):

```julia
# fermionic U-turn twist: closing the trace loop, every traced codomain leg
# with a non-dual space and odd sector contributes a factor -1
@inbounds for i in 2:length(g₁.uncoupled)
    (isdual(space(tsrc, q₁[i - 1])) || iseven(g₁.uncoupled[i].n)) && continue
    coeff = -coeff
end
```

**Example 4 — full trace of one tensor.**  A single (2,2) tensor closed on
itself is the simplest fermion loop. Its trace computed with `@grassmann`
differs from the bosonic `@tensor` result exactly by the U-turn signs of
Eq. (6).

```julia
using GTEMPO, Z2Tensors, Random

s = Z2Space(0=>1, 1=>2)
Random.seed!(7)
O = randn(ComplexF64, s⊗s, s⊗s)
@tensor    tb = O[1,2,1,2]
@grassmann tf = O[1,2,1,2]
println("bosonic   trace = ", tb)
println("fermionic trace = ", tf)
```

Output:

```
bosonic   trace = 1.1945052196905066 - 0.2070731926835736im
fermionic trace = -0.6772458932398924 + 1.1988470456294624im
```

### 5.5 Ring networks: the fermionic value is the physical one

The junction and U-turn twists are not independent decorations: they are the
*same* topological sign seen from two sides. A ring network contracts a
chain first and then closes it, so both rules act. Example 5 below builds
the two-tensor ring

$$
X: (a,b,p,q) \leftarrow (c,d) ,\qquad
Y: (c,d) \leftarrow (a,b) ,
$$

contracts the pair $(c,d)$ and then closes the $(a,b)$ legs. The
`@grassmann` result is *not* the bosonic one; the difference is the sign of
the closed fermion loop.

**Example 5 — the two-tensor ring.**

```julia
using GTEMPO, Z2Tensors, Random

s = Z2Space(0=>1, 1=>2)
Random.seed!(12345)
X = randn(ComplexF64, s⊗s⊗s⊗s, s⊗s)  # codomain (a,b,p,q) ← domain (c,d)
Y = randn(ComplexF64, s⊗s, s⊗s)       # codomain (c,d) ← domain (a,b)
@grassmann Mf[1,2,3,4;5,6] := X[1,2,3,4,7,8] * Y[7,8,5,6]  # chain over (c,d)
@grassmann Rf[3,4] := Mf[1,2,3,4,1,2]                        # close over (a,b)
@tensor    Mb[1,2,3,4;5,6] := X[1,2,3,4,7,8] * Y[7,8,5,6]
@tensor    Rb[3,4] := Mb[1,2,3,4,1,2]
println("fermionic ring value = ", Rf.data[1])
println("bosonic   ring value = ", Rb.data[1])
```

Output:

```
fermionic ring value = -2.070316090562936 + 3.335625309730484im
bosonic   ring value = -1.7323179751406874 - 1.6491679689279306im
```

(`Rf`/`Rb` are rank-(2,0) tensors with several blocks; `.data[1]` selects
the totally-even block, the same coefficient compared against TensorKit in
§7.) Which of the two is correct? The fermionic number is, and the reason is
given in §7: GTEMPO's `@grassmann` reproduces TensorKit's fermionic and
GrassmannTensors' results bit-for-bit, whereas the bosonic contraction has no
loop sign at all.

### 5.6 Addition

`tensoradd!` only permutes and sums blocks; it carries the odd–odd swap
signs of `gpermute` (Eq. (4)) but no junction twist (there is no contracted
pair). It is exercised by the in-place forms `=` / `+=` in the macro test of
§8.

### 5.7 Even-parity tensors are bosonic

If every sector in a contraction is even ($p = 0$ everywhere), Eqs. (4)–(6)
return $+1$ and `@grassmann` must coincide with `@tensor` **exactly**. This
is why all macro machinery (assignment forms, scalar outputs, arbitrary
index orders) can be tested against the plain `@tensor` on even tensors.

**Example 6.**

```julia
using GTEMPO, Z2Tensors, Random, LinearAlgebra

se = Z2Space(0=>2)
Random.seed!(1)
e1 = randn(ComplexF64, se⊗se, se⊗se)
e2 = randn(ComplexF64, se⊗se, se⊗se)
@grassmann cg[1,2;5,6] := e1[1,2,3,4] * e2[3,4,5,6]
@tensor    cb[1,2;5,6] := e1[1,2,3,4] * e2[3,4,5,6]
println("‖@grassmann − @tensor‖ = ", norm(cg - cb))
```

Output:

```
‖@grassmann − @tensor‖ = 0.0
```

## 6. Where the rules live

| operation | file | formula |
|---|---|---|
| index permutation | `src/grassmanntensor/grassmanntensor.jl` (`gpermute`, `add_gpermute!`) | Eq. (4), odd–odd swap signs |
| contraction junction | `src/grassmanntensor/tensoroperations.jl` (`_contract!`) | Eq. (5), twist $(-1)^{\sum p}$ on crossed A-side ket legs |
| trace closure | `src/grassmanntensor/tensoroperations.jl` (`trace_permute!`) | Eq. (6), U-turn twist $(-1)^{\#}$ of traced non-dual odd legs |
| block twist helper | `src/grassmanntensor/grassmanntensor.jl` (`g_twist!`) | multiplies blocks by $(-1)^{\#\text{odd twisted legs}}$ |

All of these are keyed on the trailing `backend::GrassmannBackend` argument,
which distinguishes them from the sign-free Z2Tensors operations with
otherwise similar signatures. Plain `@tensor` never sees them.

## 7. Relation to TensorKit and GrassmannTensors (cross-checks)

The `GrassmannBackend` rules are numerically identical to the fermionic
(`FermionParity`) behavior of TensorKit's `blas_contract!` /
`_trace_permute!` and of the GrassmannTensors package. Two implementation
differences are worth noting:

* TensorKit/GrassmannTensors fusion trees carry an `isdual` array; Z2Tensors
  trees do not, and the flags are recovered from the spaces. One can show
  the recovery formula `isdual(tree leg) ↔ isdual(space(tsrc, q₁[i-1]))`
  makes the trace rule identical.
* TensorKit's `blas_contract!` chooses (by a copy-cost heuristic) whether to
  twist A's or B's legs; the two are equivalent for the result, and
  `GrassmannBackend` always twists the A side.

The cross-checks were done with three independent implementations of the
*identical* tensor network:

**Junction probe** (Example 3's tensors). In GrassmannTensors the same
expression gives the same flags and the same number:

```julia
# (GrassmannTensors, an independent fermionic tensor package)
using GrassmannTensors, TensorOperations
s = GradedSpace(FermionParity(0)=>1, FermionParity(1)=>1); sp = dual(s)
A = zeros(ComplexF64, s⊗s, sp⊗sp); B = zeros(ComplexF64, sp⊗sp, s)
for (_, b) in blocks(A); b .= 1 + 0.5im; end
for (_, b) in blocks(B); b .= 1 + 0.5im; end
@tensor Rf[x,y;cc] := A[x,y,aa,bb] * B[aa,bb,cc]
println("A flags: ", [isdual(space(A,i)) for i in 1:4])
println("B flags: ", [isdual(space(B,i)) for i in 1:3])
println("grassmann Rf = ", first(blocks(Rf))[2][1])
```

Output (identical to GTEMPO, Example 3):

```
A flags: Bool[0, 0, 0, 0]
B flags: Bool[1, 1, 1]
grassmann Rf = 1.5 + 2.0im
```

**Ring probe.** The two-tensor ring of Example 5 can equally be contracted
in the opposite order (cross $(a,b)$ first, trace $(c,d)$ afterwards). A
consistent fermionic convention must give the same number regardless of the
intermediate contraction order. TensorKit's fermionic tensors give:

```
[all-odd] R1 = -0.7335205857637537 + 0.06011419991705536im
[all-odd] R2 = -0.7335205857637537 + 0.06011419991705536im
[mixed  ] R1 = -2.070316090562936 + 3.335625309730484im
[mixed  ] R2 = -2.0703160905629354 + 3.335625309730484im
```

GTEMPO's `@grassmann` reproduces these numbers bit-for-bit on the same
seeded data — in particular the mixed-space ring value
$-2.070316090562936 + 3.335625309730484\mathrm{i}$ of Example 5. Both
frameworks are therefore order-consistent *and* mutually identical, while
the bosonic value (Example 5) differs.

## 8. Unit and integration tests

The sign rules are covered both by direct unit tests and by the full
physical test suite.

**`test/grassmanntensor.jl`** tests the permutation signs against a
hand-written sign loop (bosonic result with the odd–odd signs of Eq. (4)
applied manually must equal the `@grassmann` permute), the chain contraction
(`@grassmann == @tensor` on an open graph, cf. Example 3), and all macro
forms on even tensors (cf. Example 6). Run it with

```bash
julia --project=<env-with-GTEMPO> -e 'using GTEMPO, Test, Random; include("test/grassmanntensor.jl")'
```

Output:

```
------------------------------------
|        Grassmann Tensor          |
------------------------------------
Test Summary:           | Pass  Total  Time
GrassmannTensor permute |    8      8  5.1s
Test Summary:            | Pass  Total  Time
GrassmannTensor contract |    2      2  5.0s
Test Summary:    | Pass  Total  Time
@grassmann macro |   10     10  9.1s
```

**Full suite** (`test/runtests.jl`, `RUN_HEAVY_TESTS=true`, single process):

```bash
OMP_NUM_THREADS=1 julia --project=<env-with-GTEMPO> test/runtests.jl
```

Result on 2026-09-04: **115 testsets, 0 failures**. The suite covers the
two readings of §3 end-to-end (GrassmannMPS multiplication and time
evolution via `@grassmann`; inner products, observables and influence
functional construction via `@tensor`), so the `conj(...) → bosonic` rule of
§3 is enforced at the level of physical results, not only of unit tests.
(The current version history and the original interface baseline are kept
under the project memory, `baseline_b57b240/`.)

## 9. Summary of the design principles

1. **One storage, two readings.** Coefficients live in bosonic Z2Tensors;
   the G-number string is virtual.
2. **ket–ket → `@grassmann`; anything with `conj` → `@tensor`.** The ket–bra
   contraction is a coefficient pairing (Eq. (2)) and is bosonic by
   convention.
3. **The canonical contracted pair is $(a,\bar a)$.** Crossed $(\bar a,a)$
   junctions pay a twist (Eq. (5)), automated in `GrassmannBackend.contract!`.
4. **Closed fermionic loops pay the U-turn twist** (Eq. (6)), automated in
   `GrassmannBackend.trace_permute!`.
5. **Manual signs are confined to the coefficient world** (JW strings inside
   the influence-functional MPOs), where they can be checked against known
   limits; the GT world needs no hand-inserted signs.
