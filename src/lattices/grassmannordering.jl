abstract type GrassmannOrdering end
abstract type ImagGrassmannOrdering <: GrassmannOrdering end
abstract type RealGrassmannOrdering <: GrassmannOrdering end
abstract type MixedGrassmannOrdering <: GrassmannOrdering end

abstract type ConjugationStyle end
struct AdjacentConjugation <: ConjugationStyle end
struct GeneralConjugation <: ConjugationStyle end
ConjugationStyle(x::GrassmannOrdering) = ConjugationStyle(typeof(x))

abstract type LayoutStyle end
struct TimeLocalLayout <: LayoutStyle end
struct BandLocalLayout <: LayoutStyle end
"""
	struct BranchLocalLayout

TimeTimeLocalLayout intra branch, but the two branches are separated
"""
struct BranchLocalLayout <: LayoutStyle end
LayoutStyle(x::GrassmannOrdering) = LayoutStyle(typeof(x))

### the auxiliary Grassmann numbers for boundary dynamics are assumed to be on the left boundary!

"""
	struct A1A1B1B1

GrassmannOrdering aābb̄ a₂ā₂b₂b̄₂ a₁ā₁b₁b̄₁
"""
struct A1A1B1B1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A1A1B1B1}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1B1B1}) = TimeLocalLayout()
const AABB = A1A1B1B1

"""
	struct A1B1B1A1

GrassmannOrdering abb̄ā a₂b₂b̄₂ā₂ a₁b₁b̄₁ā₁
This ordering is convenient to build the impurity dynamics
"""
struct A1B1B1A1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A1B1B1A1}) = GeneralConjugation()
LayoutStyle(::Type{A1B1B1A1}) = TimeLocalLayout()
const ABBA = A1B1B1A1

"""
	struct A2A2A1A1B2B2B1B1

GrassmannOrdering aābb̄ a₂ā₂a₁ā₁ b₂b̄₂b₁b̄₁
"""
struct A2A2A1A1B2B2B1B1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A2A2A1A1B2B2B1B1}) = AdjacentConjugation()
LayoutStyle(::Type{A2A2A1A1B2B2B1B1}) = BandLocalLayout()


# caption means forward, small means backword

"""
	struct A1A1a1a1B1B1b1b1

GrassmannOrdering aābb̄ a₂+ā₂+b₂+b̄₂+a₂-ā₂-b₂-b̄₂- a₁+ā₁+b₁+b̄₁+a₁-ā₁-b₁-b̄₁-
"""
struct A1A1B1B1a1a1b1b1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1A1B1B1a1a1b1b1}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1B1B1a1a1b1b1}) = TimeLocalLayout()
const AABBaabb = A1A1B1B1a1a1b1b1

"""
	struct A1A1a1a1B1B1b1b1

GrassmannOrdering aābb̄ a₂+ā₂+a₂-ā₂-b₂+b̄₂+b₂-b̄₂- a₁+ā₁+a₁-ā₁-b₁+b̄₁+b₁-b̄₁-
"""
struct A1A1a1a1B1B1b1b1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1A1a1a1B1B1b1b1}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1a1a1B1B1b1b1}) = TimeLocalLayout()
const AAaaBBbb = A1A1a1a1B1B1b1b1

"""
	struct A1a1B1b1b1B1a1A1

GrassmannOrdering abb̄ā a₂+a₂-b₂+b₂-b̄₂-b̄₂+ā₂-ā₂+ a₁+a₁-b₁+b₁-b̄₁-b̄₁+ā₁-ā₁+
This is a historical ordering, which may not be very useful
"""
struct A1a1B1b1b1B1a1A1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1a1B1b1b1B1a1A1}) = GeneralConjugation()
LayoutStyle(::Type{A1a1B1b1b1B1a1A1}) = TimeLocalLayout()
const AaBbbBaA = A1a1B1b1b1B1a1A1

"""
	struct A1a1B1b1b1B1a1A1

GrassmannOrdering aābb̄ a₂+b₂+ā₂-b̄₂-ā₂+b̄₂+a₂-b₂-  a₁+b₁+ā₁-b̄₁-ā₁+b̄₁+a₁-b₁-
This ordering is convenient to build the impurity dynamics for time local ordering
"""
struct A1B1ā1b̄1A1B1a1b1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1B1ā1b̄1A1B1a1b1}) = GeneralConjugation()
LayoutStyle(::Type{A1B1ā1b̄1A1B1a1b1}) = TimeLocalLayout()


"""
	struct A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1

GrassmannOrdering aābb̄ a₂+ā₂+a₁+ā₁+ a₂-ā₂-a₁-ā₁- b₂+b̄₂+b₁+b̄₁+  b₂-b̄₂-b₁-b̄₁-
"""
struct A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}) = AdjacentConjugation()
LayoutStyle(::Type{A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}) = BandLocalLayout()

"""
	struct A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2

GrassmannOrdering abb̄ā a₂+b₂+b̄₂+ā₂+a₁+b₁+b̄₁+ā₁+ a₁-b₁-b̄₁-ā₁-a₂-b₂-b̄₂-ā₂-
This ordering is convenient to build the impurity dynamics for band local ordering
"""
struct A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = GeneralConjugation()
LayoutStyle(::Type{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = BranchLocalLayout()

"""
	struct A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2

GrassmannOrdering aābb̄ a₂+ā₂+b₂+b̄₂+a₁+ā₁+b₁+b̄₁+ a₁-ā₁-b₁-b̄₁-a₂-ā₂-b₂-b̄₂-
This ordering is convenient to build the impurity dynamics
"""
struct A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}) = AdjacentConjugation()
LayoutStyle(::Type{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}) = BranchLocalLayout()


# mixed time
"""
	struct A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2

A1A1B1B1 for ImagGrassmannOrdering 
A1A1a1a1B1B1b1b1 for RealGrassmannOrdering, the real-time index is from small to large!!!
"""
struct A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2}) = TimeLocalLayout()
const AABB_AAaaBBbb = A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2

"""
	struct A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2

A1A1B1B1 for ImagGrassmannOrdering 
a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2 for RealGrassmannOrdering, the real-time index is from small to large!!!
"""
struct A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2}) = TimeLocalLayout()
const AABB_aaAAbbBB = A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2

"""
	struct A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2

A1B1B1A1 for ImagGrassmannOrdering 
A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2 for RealGrassmannOrdering
"""
struct A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = GeneralConjugation()
LayoutStyle(::Type{A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = BranchLocalLayout()
