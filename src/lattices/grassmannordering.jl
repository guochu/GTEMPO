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

GrassmannOrdering a\bar{a}b\bar{b} a_2\bar{a}_2b_2\bar{b}_2 a_1\bar{a}_1b_1\bar{b}_1
"""
struct A1A1B1B1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A1A1B1B1}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1B1B1}) = TimeLocalLayout()
const AABB = A1A1B1B1

"""
	struct A1B1B1A1

GrassmannOrdering ab\bar{b}\bar{a} a_2b_2\bar{b}_2\bar{a}_2 a_1b_1\bar{b}_1\bar{a}_1
This ordering is convenient to build the impurity dynamics
"""
struct A1B1B1A1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A1B1B1A1}) = GeneralConjugation()
LayoutStyle(::Type{A1B1B1A1}) = TimeLocalLayout()
const ABBA = A1B1B1A1

"""
	struct A2A2A1A1B2B2B1B1

GrassmannOrdering a\bar{a}b\bar{b} a_2\bar{a}_2a_1\bar{a}_1 b_2\bar{b}_2b_1\bar{b}_1
"""
struct A2A2A1A1B2B2B1B1 <: ImagGrassmannOrdering end
ConjugationStyle(::Type{A2A2A1A1B2B2B1B1}) = AdjacentConjugation()
LayoutStyle(::Type{A2A2A1A1B2B2B1B1}) = BandLocalLayout()


# caption means forward, small means backword

"""
	struct A1A1a1a1B1B1b1b1

GrassmannOrdering a\bar{a}b\bar{b} a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-
"""
struct A1A1a1a1B1B1b1b1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1A1a1a1B1B1b1b1}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1a1a1B1B1b1b1}) = TimeLocalLayout()
const AAaaBBbb = A1A1a1a1B1B1b1b1

"""
	struct A1a1B1b1b1B1a1A1

GrassmannOrdering ab\bar{b}\bar{a} a_2^+a_2^-b_2^+b_2^-\bar{b}_2^-\bar{b}_2^+\bar{a}_2^-\bar{a}_2^+ a_1^+a_1^-b_1^+b_1^-\bar{b}_1^-\bar{b}_1^+\bar{a}_1^-\bar{a}_1^+
This is a historical ordering, which may not be very useful
"""
struct A1a1B1b1b1B1a1A1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A1a1B1b1b1B1a1A1}) = GeneralConjugation()
LayoutStyle(::Type{A1a1B1b1b1B1a1A1}) = TimeLocalLayout()
const AaBbbBaA = A1a1B1b1b1B1a1A1

"""
	struct A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1

GrassmannOrdering a\bar{a}b\bar{b} a_2^+\bar{a}_2^+a_1^+\bar{a}_1^+ a_2^-\bar{a}_2^-a_1^-\bar{a}_1^- b_2^+\bar{b}_2^+b_1^+\bar{b}_1^+  b_2^-\bar{b}_2^-b_1^-\bar{b}_1^-
"""
struct A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}) = AdjacentConjugation()
LayoutStyle(::Type{A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}) = BandLocalLayout()

"""
	struct A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2

GrassmannOrdering ab\bar{b}\bar{a} a_2^+b_2^+\bar{b}_2^+\bar{a}_2^+a_1^+b_1^+\bar{b}_1^+\bar{a}_1^+ a_1^-b_1^-\bar{b}_1^-\bar{a}_1^-a_2^-b_2^-\bar{b}_2^-\bar{a}_2^-
This ordering is convenient to build the impurity dynamics
"""
struct A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = GeneralConjugation()
LayoutStyle(::Type{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = BranchLocalLayout()

"""
	struct A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2

GrassmannOrdering a\bar{a}b\bar{b} a_2^+\bar{a}_2^+b_2^+\bar{b}_2^+a_1^+\bar{a}_1^+b_1^+\bar{b}_1^+ a_1^-\bar{a}_1^-b_1^-\bar{b}_1^-a_2^-\bar{a}_2^-b_2^-\bar{b}_2^-
This ordering is convenient to build the impurity dynamics
"""
struct A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2 <: RealGrassmannOrdering end
ConjugationStyle(::Type{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}) = AdjacentConjugation()
LayoutStyle(::Type{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}) = BranchLocalLayout()


# mixed time
"""
	struct A1A1B1B1_A1A1a1a1B1B1b1b1

A1A1B1B1 for ImagGrassmannOrdering 
A1A1a1a1B1B1b1b1 for RealGrassmannOrdering
"""
struct A1A1B1B1_A1A1a1a1B1B1b1b1 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1A1B1B1_A1A1a1a1B1B1b1b1}) = AdjacentConjugation()
LayoutStyle(::Type{A1A1B1B1_A1A1a1a1B1B1b1b1}) = TimeLocalLayout()
const AABB_AAaaBBbb = A1A1B1B1_A1A1a1a1B1B1b1b1

struct A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2 <: MixedGrassmannOrdering end
ConjugationStyle(::Type{A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = GeneralConjugation()
LayoutStyle(::Type{A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}) = BranchLocalLayout()
