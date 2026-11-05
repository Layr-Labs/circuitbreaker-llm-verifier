# One fixed uniform language with O(log(N)^28 loglog(N)^7) loss

For the single total language defined here there are absolute constants
B>0 and N0 such that, for EVERY integer actual length N>=N0,

 C_B2(F_N) >=4N-B*(log2 N)^28*(log2(log2 N))^7.

The native leading projection is (4,1). More explicitly, for the SAME
family and every fixed epsilon>0, eventually

 C_B2(F_N)>=4N-(log2 N)^28*(log2(log2 N))^(7+epsilon),

and hence also 4N-(log2 N)^(28+epsilon) and 4N-N^epsilon.
The coefficient B is asserted to exist as a fixed constant, not given an
unproved numerical value. The source guarantee and projection are asserted
jointly by both canonical claims.

The language is fully specified without knowing B or choosing a good scale.
It reads an integer scale from actual input bits and evaluates every scale.
For hardness we prove that ONE fixed scale value S0 supplies a hard
restriction at all sufficiently large lengths. Header and field data are
all counted in actual input length. Existence of this hard restriction is
not used to establish existence of the language or an efficient verifier.

## 1. The exact total algorithm

On input x of length N, family.py first checks its binary alphabet. If
N.bit_length()<1024 it returns the fixed legacy evaluator. Otherwise put
L=ceil(log2 N), z=ceil(log2 L). The first L bits, in little-endian order,
encode S-1. Thus 1<=S<=2^L<=2N. The program calls engine.PolyPlan(N,S).
If ready is false, it returns ZERO. If ready is true, it reads the stated
field data AFTER the header and runs engine.PolyEngine. Its payload is the
suffix after the ENTIRE header-plus-field prefix, padded with zeros to N
coordinates. The algorithm is fixed for every input and length; there is
no promised good header or externally supplied parameter.

Let P2(t) be the least power of two at least positive integer t. The exact
scales computed from N and the actual header S are

 k=P2(S^15*L^28*z^7), h=e=P2(S*L^3),
 ell=P2(S^3*L^6*z^3), e'=P2(S^7*L^15*z^3),
 d=P2(S^11*L^21*z^8), Y=P2(S^13*L^25*z^7).

Y is the initial row width. The plan compares log2 k with floor(log2 N)
before expanding k. Large rejected plans give zero; only the finite small-
length cutoff uses the legacy routine. Every operation is integer arithmetic.

## 2. All construction data are input coordinates

The first L bits are the scale header. They are followed by the six A1
packages for initial extraction,Q-side,X-side,E1,E2,E3; then the hash
modulus and the design A1 package. Initial dimension is 10000 FIXED, the
five internal dimensions are z, and the design dimension is 100. An A1
package contains exactly (D^2+1)s bits: a monic modulus's lower coefficients
and every entry of its D-by-D traversing matrix. Repeated calls share data
while using independent seed substrings. Given N and the actual header,
every size, ordering and bit position is computed deterministically.

The engine starts its field-data cursor after the header, checks it has
consumed exactly data_prefix bits, and removes header+data_prefix from the
payload. Its prefix variable is the TOTAL of those two lengths. The hard
restriction fixes the header and every field/matrix datum as actual input
bits. None is a witness, advice string, external parameter or oracle.

Good field tuples exist by A1's elementary finite-field and cyclic-group
arguments. One tuple works simultaneously for all sources at the stated
length and scale. Bad tuples remain ordinary inputs: monic quotient-ring
multiplication and matrix operations are bounded, so they still return a
Boolean value. Their goodness is not tested or promised.

The ambient-degree condenser modulus is COMPUTED by B1's deterministic
Artin--Schreier tower and relative-minimal-polynomial algorithm. It is not
an omitted order-N input datum. The supplied algebra proves every such
field and modulus and gives polynomial-time construction. Its multiplication
uses an explicit stack, without unbounded interpreter recursion.

## 3. Supplied extraction and gate arguments

A1 supplies the low-degree strong seeded extractor including its list-
reconstruction and binary-output proof. B1 supplies the polynomial graph
expansion/root argument with power-of-two exponents for binary linearity.
On an affine K-flat the neighbor count is sum_y2^rank(C_y), directly
bounding rank-deficient seeds. No unexplained modern affine extractor is
used.

B2 composes that rank condenser with fresh seeded calls on the SAME affine
input. Fixing earlier seeds and outputs preserves affine structure and
loses at most their output rank. The seed bound retains D explicitly.
For the internal maps, D=z, H=O(z),R=O(z), so their half-output seed cost
is O(z^3*e^2). The short E1 output is separately checked to fit ONE call.
Initial dimension stays fixed, giving d0<=10^9 L before any enumeration.

B3--B4 supply the complete adaptive affine look-ahead/index proof. Every
seed is fixed before an output is conditioned on. A-parts are fixed for
ALL streams, including the designated stream, before using its independent
B-part. Same-prefix next seeds are not fixed prematurely. Each index round
refreshes with E1 and E2, INCLUDING the last. The joint row error is at
most100bh^2*2^-e for b=ceil(d0/log2 h). All initial indices except at most
2M0/N^2 are good, M0=2^d0. On any h good indices their selected rows are
almost independent. Taking their first bits gives the source used in B6.

B5 supplies the general-min-entropy shifted-agreement generator from a
polynomial hash prefix and fresh A1 calls, with exceptional-conditioning
probabilities accounted for. B6 supplies the monotone partition function,
its exact integer threshold search, a proved Harris/Janson inequality,
and the two-level Bonferroni comparison. It works with arbitrary coupling
of the bad bits, including extra conditional randomness. Its proof yields
finite constants A_t,A_gamma for the requirements

 t>=A_t*r^3, gamma<=2^(-A_gamma*r^3).

They follow from its fixed O(r^3) degree and coefficient-count constants;
they are not an external resilient-function premise. The threshold search
uses only exact comparisons

 2^r*(2^(uB)-(2^u-J)^B)^M <=2^(uBM),

with polynomial bit cost in M0. No balancing advice or real oracle is used.

The finite full-B2 theorem in elimination.md Section12 gives4n-8d for
an affine disperser of dimension d. Its promised-domain representation,
all case reductions and four-distinct-object refinement are supplied.
It allows all16 binary truth tables, arbitrary sharing, constants, aliases,
unused inputs and all terminal cases. The older supplemental atomic gate
claims are retained with their separate invariants.

## 4. One fixed good SCALE for the hardness restriction

B7 proves the following estimates with constants independent of a FIXED
scale S, once N is large enough for that S. Internal dimensions D=z have
H=O(z), because log S/z tends to zero and the L,z exponents are fixed.
There are finite constants C_seed,C_advice such that half-output seed cost
is at most C_seed*z^3*e^2 plus lower terms, and E1's single-call advice is
at most C_advice*ell*z^2*e' plus lower terms. Initial extraction uses one
call eventually, since its advice has L exponent26.0028<28. Its cap
d0<=10^9 L is independent of fixed S; the threshold may depend on S.

Let B0=2*10^9, so eventually b<=B0 L/z. Rounding each power to P2 costs
at most two. The key bounds are

 z^3e^2 <=4S^2L^6z^3,       ell>=S^3L^6z^3,
 ell*z^2*e' <=4S^10L^21z^8, d>=S^11L^21z^8,
 h^3ell <=16S^6L^15z^3,     e'>=S^7L^15z^3,
 bhd <=4B0S^12L^25z^7,     Y>=S^13L^25z^7,
 hY <=4S^14L^28z^7,        k>=S^15L^28z^7.

Each required coefficient ratio gains one power of S. Thus choosing one
sufficiently large S makes every seed, row-rank and weak-seed inequality
strict. In particular S>=3200 leaves weak-seed room, S>=8000B0 leaves
Y>=1000bhd room, and S>=16 gives hY<=k/4. Other source losses and 10d
have smaller L powers than k. The code retains the actual entropy checks.

The layout has the exact bound r<(d0+1)*10^12<=c_r L with c_r=2*10^21.
Choose ONE integer S0 exceeding the finitely many bounds just proved,
including the seed/advice constants and

 S0>=A_t*c_r^3,
 S0>=2A_gamma*c_r^3,
 S0>=100c_r^3.

For this fixed S0, eventually h>=A_t r^3 and e>=2A_gamma r^3. The
polynomial prefactor in100bh^2*2^-e is negligible in the exponent, so
B6's gamma condition holds. Its coordinate condition follows from
uS_groups>M0/5,d0<=10^9 L and delta=10^-12:

 q/S_groups<=10u/N^2, S_groups^delta<=2^(.001L).

The product tends to zero. Thus every B6 hypothesis holds. The weaker
executable checks h>=r^3,e>=100r^3 are sufficient resource checks, not a
claim that all admitted headers certify dispersion.

This argument is not a nonconstructive definition of the language. The
fixed evaluator already computes ALL header values, including small scales
that may be easy. S0 is used only in the circuit restriction. It is a
fixed integer independent of N, and fits in L header bits at all sufficiently
large lengths. The same algorithms and language work for every length.
The verifier never receives S0 from anywhere except its ordinary input.

## 5. Eligibility for EVERY header and invalid field string

The scale parsed from the input satisfies S<=2N. Raw scalar products have
polynomially many bits. If k>N the plan rejects before expanding it. On a
remaining plan all source, output and precision parameters are polynomially
bounded in N; the explicit component checks reject incompatible dimensions.
The prefix+k<=N check bounds admitted actual data allocation. The d0 cap
bounds every enumerated seed uniformly over headers.

A1 streams only the FIRST n' coefficient monomials with polynomial delay
in D,H; it does not enumerate an entire D-dimensional grid. D=z grows
only logarithmically in log N and all operation counts are polynomial in
the displayed parameters, not N^D. Internal fields are evaluated by bounded
algebraic operations, never enumerated. Deterministic field dimensions are
polynomial in N. The robust loops have r<=10^12(d0+1)=O(log N), hence
fixed polynomial length. Its integer threshold arithmetic has polynomial
bit cost in M0, also polynomial in N.

Every large rejected plan returns zero, so infinitely many bad headers do
not invoke an unproved fallback. Invalid fields on an admitted plan still
give bounded ring/matrix computation. Only finitely many small cutoff
lengths invoke the fixed legacy algorithm, which halts by existence of
primes and its bounded subsequent loops. Those finitely many cases contribute
a finite constant to a polynomial-time bound. This proves a total uniform
polynomial-time language independently of any good S0 or good fields.

For NP define V(x,w) to reject every nonempty witness and otherwise accept
exactly when family.evaluate(x)=1. Witness length is zero; runtime is
polynomial, and both directions of membership hold. Header and construction
data belong to x, not w. No promise is imposed on x.

## 6. Actual-length full-B2 accounting

Fix the header to S0-1 and fix good field/matrix data at actual large length
N. Their total prefix length is

 P(N,S0)=L+O(S0^7L^15z^5)=o(k(N,S0)).

This includes the ENTIRE header. Once the header is fixed the following
field positions are static for N,S0. Good tuples work simultaneously for
all affine sources. Padding every payload flat injectively into N coordinates
and applying B1--B7 proves an affine disperser at dimension k on the
N-P genuine payload bits. The finite gate theorem yields

 C_B2(F_N)>=4(N-P)-8k=4N-R(N), R(N)=4P(N,S0)+8k(N,S0).

Restriction of a circuit cannot increase its gate count. Every input bit
fixed in this argument has therefore been charged, including the scale
header and every repeated-use matrix. Since P=o(k), eventually R<=9k and
k<2S0^15L^28z^7. As L/log2 N->1 and z/log2(log2 N)->1, we may take the
fixed constant B=32S0^15 in the opening bound, for all sufficiently large
integer actual N. B is not claimed numerically small or explicitly computed.

For each fixed epsilon>0, that constant is eventually absorbed by
(log2(log2 N))^epsilon. The SAME family therefore has the stated numerical
weaker bound with loglog exponent7+epsilon, and also the stated log and
N^epsilon bounds. The native projection remains exactly(4,1). The first
two canonical claims assert this entire declaration with the same evidence;
coefficient and exponent are not selected from incompatible families.

## 7. Evidence scope and provenance

The finite audit establishes only its named bounded predicates. It does
not execute the enormous complete advanced branch. Small full-family
fixtures use the legacy cutoff. Prefix-offset and scalar-plan mocks check
ordering and accounting, not universal entropy claims. Those universal
claims rest on the supplied arguments. Independent critique is advisory;
configured-provider and human judgments are separate.

Gate provenance is Mishra--Santhanam revision1 for the coefficient-four
architecture, followed by the supplied Section12 finite refinement.
Extraction provenance is specified in disperser.md. Retained baseline and
predecessor gate/code contributions are credited in the public note. Their
earlier acceptance is not a mathematical premise.
