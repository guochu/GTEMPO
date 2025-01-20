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
	struct A1Ā1B1B̄1

GrassmannOrdering aābb̄ a₂ā₂b₂b̄₂ a₁ā₁b₁b̄₁
"""
struct A1Ā1B1B̄1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A1Ā1B1B̄1}) = AdjacentConjugation()
LayoutStyle(::Type{A1Ā1B1B̄1}) = TimeLocalLayout()
const AĀBB̄ = A1Ā1B1B̄1

"""
	struct A1B1B̄1Ā1

GrassmannOrdering abb̄ā a₂b₂b̄₂ā₂ a₁b₁b̄₁ā₁
This ordering is convenient to build the impurity dynamics
"""
struct A1B1B̄1Ā1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A1B1B̄1Ā1}) = GeneralConjugation()
LayoutStyle(::Type{A1B1B̄1Ā1}) = TimeLocalLayout()
const ABB̄Ā = A1B1B̄1Ā1

"""
	struct A2Ā2A1Ā1B2B̄2B1B̄1

GrassmannOrdering aābb̄ a₂ā₂a₁ā₁ b₂b̄₂b₁b̄₁
"""
struct A2Ā2A1Ā1B2B̄2B1B̄1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A2Ā2A1Ā1B2B̄2B1B̄1}) = AdjacentConjugation()
LayoutStyle(::Type{A2Ā2A1Ā1B2B̄2B1B̄1}) = BandLocalLayout()


# caption means forward, small means backword

"""
	struct A1Ā1a1ā1B1B̄1b1b̄1

GrassmannOrdering aābb̄ a₂+ā₂+b₂+b̄₂+a₂-ā₂-b₂-b̄₂- a₁+ā₁+b₁+b̄₁+a₁-ā₁-b₁-b̄₁-
"""
struct A1Ā1B1B̄1a1ā1b1b̄1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1Ā1B1B̄1a1ā1b1b̄1}) = AdjacentConjugation()
LayoutStyle(::Type{A1Ā1B1B̄1a1ā1b1b̄1}) = TimeLocalLayout()
const AĀBB̄aābb̄ = A1Ā1B1B̄1a1ā1b1b̄1

"""
	struct A1Ā1a1ā1B1B̄1b1b̄1

GrassmannOrdering aābb̄ a₂+ā₂+a₂-ā₂-b₂+b̄₂+b₂-b̄₂- a₁+ā₁+a₁-ā₁-b₁+b̄₁+b₁-b̄₁-
"""
struct A1Ā1a1ā1B1B̄1b1b̄1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1Ā1a1ā1B1B̄1b1b̄1}) = AdjacentConjugation()
LayoutStyle(::Type{A1Ā1a1ā1B1B̄1b1b̄1}) = TimeLocalLayout()
const AĀaāBB̄bb̄ = A1Ā1a1ā1B1B̄1b1b̄1

"""
	struct A1Ā1B1B̄1b̄1B̄1ā1Ā1

GrassmannOrdering abb̄ā a₂+a₂-b₂+b₂-b̄₂-b̄₂+ā₂-ā₂+ a₁+a₁-b₁+b₁-b̄₁-b̄₁+ā₁-ā₁+
This is a historical ordering, which may not be very useful
"""
struct A1Ā1B1B̄1b̄1B̄1ā1Ā1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1Ā1B1B̄1b̄1B̄1ā1Ā1}) = GeneralConjugation()
LayoutStyle(::Type{A1Ā1B1B̄1b̄1B̄1ā1Ā1}) = TimeLocalLayout()
const AaBbb̄B̄āĀ = A1Ā1B1B̄1b̄1B̄1ā1Ā1

"""
	struct A1Ā1B1B̄1b̄1B̄1ā1Ā1

GrassmannOrdering aābb̄ a₂+b₂+ā₂-b̄₂-ā₂+b̄₂+a₂-b₂-  a₁+b₁+ā₁-b̄₁-ā₁+b̄₁+a₁-b₁-
This ordering is convenient to build the impurity dynamics for time local ordering
"""
struct A1B1ā1b̄1Ā1B̄1a1b1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1B1ā1b̄1Ā1B̄1a1b1}) = GeneralConjugation()
LayoutStyle(::Type{A1B1ā1b̄1Ā1B̄1a1b1}) = TimeLocalLayout()


"""
	struct A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1

GrassmannOrdering aābb̄ a₂+ā₂+a₁+ā₁+ a₂-ā₂-a₁-ā₁- b₂+b̄₂+b₁+b̄₁+  b₂-b̄₂-b₁-b̄₁-
"""
struct A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1}) = AdjacentConjugation()
LayoutStyle(::Type{A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1}) = BandLocalLayout()

"""
	struct A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2

GrassmannOrdering abb̄ā a₂+b₂+b̄₂+ā₂+a₁+b₁+b̄₁+ā₁+ a₁-b₁-b̄₁-ā₁-a₂-b₂-b̄₂-ā₂-
This ordering is convenient to build the impurity dynamics for band local ordering
"""
struct A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}) = GeneralConjugation()
LayoutStyle(::Type{A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}) = BranchLocalLayout()

"""
	struct A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2

GrassmannOrdering aābb̄ a₂+ā₂+b₂+b̄₂+a₁+ā₁+b₁+b̄₁+ a₁-ā₁-b₁-b̄₁-a₂-ā₂-b₂-b̄₂-
This ordering is convenient to build the impurity dynamics
"""
struct A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2}) = AdjacentConjugation()
LayoutStyle(::Type{A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2}) = BranchLocalLayout()


# mixed time
"""
	struct A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2

A1Ā1B1B̄1 for ImagGrassmannOrdering 
A1Ā1a1ā1B1B̄1b1b̄1 for RealGrassmannOrdering, the real-time index is from small to large!!!
"""
struct A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2}) = AdjacentConjugation()
LayoutStyle(::Type{A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2}) = TimeLocalLayout()
const AĀBB̄_AĀaāBB̄bb̄ = A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2

"""
	struct A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2

A1Ā1B1B̄1 for ImagGrassmannOrdering 
A1Ā1a1ā1B1B̄1b1b̄1a2a2A2A2b2b2B2B2 for RealGrassmannOrdering, the real-time index is from small to large!!!
"""
struct A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}) = AdjacentConjugation()
LayoutStyle(::Type{A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}) = TimeLocalLayout()
const AĀBB̄_aāAĀbb̄BB̄ = A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2

"""
	struct A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2

A1B1B̄1Ā1 for ImagGrassmannOrdering 
A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2 for RealGrassmannOrdering
"""
struct A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}) = GeneralConjugation()
LayoutStyle(::Type{A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}) = BranchLocalLayout()
