# Weighted marked-subset selection and one uniform coefficient-four family

The new result is a smaller sufficient affine dimension for the supplied
algebraic payload, proved by weighted averaging over marked subsets.
Its full-B2 certificate uses the retained degree-sensitive elimination
bound. Exact global selection never increases the original supplied
certificate loss, and strictly decreases it at the explicitly checked
lengths below. The leading coefficient/exponent remains(4,1).
No leading coefficient greater than4 or actual circuit-size ratio is claimed.

## 1. Exact global family at every length

For every integer N>=0 enumerate all primes2<=k<N. Put m=N-k and
r=ceil(m/k), requiring r>=3 and k>r+1. Enumerate every integer2<=t<r.
For this pair define

    H=r(r+3)/2,
    T=(t-1)k-t(t+1)/2+1,
    U=rk-r(r+1)/2,
    d_old=H+max(T,floor(U/(t+1)))+1,
    B=t(t+3)/2.

If k>=2r+1 AND2k>=(t+1)^2, set

    h=ceil(2r/(t+1)),
    F=r*(2k-2r-1+h)/(2(t+1)),
    d=H+max(T,floor(F))+1.

Otherwise set d=d_old. Call(k,t) admissible when B<=d<m as well.
Its complete loss is ell_N(k,t)=4k+8d-4B-3. Select the lexicographically
least(ell_N(k,t),k,t). If no pair is admissible, output0 at the entire
length and define L_N=4N. Otherwise L_N is the selected loss. All bounds,
ceilings, floors and ties are exact integers/rationals, without real or
floating-point choices. In particular the empty search defines N=0.

On an admissible length read the first k actual bits as lower coefficients
of Q(Z)=Z^k+sum_(j=0)^(k-1)b_j Z^j overF_2. Require
Z^(2^k)=Z modQ and gcd(Q,Z^2+Z)=1, otherwise output0. For a passing
prefix use K=F_2[Z]/Q. Pad the m payload bits by rk-m zeros into r blocks
a_i in the polynomial basis, first bit the constant coordinate. Put

    y_i=a_i^(2^(r+2-i)-1),  1<=i<=r,

and output the constant coordinate of the elementary symmetric polynomial
of order t in these y_i. Initialize E_0=1 and other E_j=0; for each y_i,
update E_j<-E_j+y_i*E_(j-1) in descending j=min(i,t),...,1. Splitting
subsets by whether they contain i proves this division-free recurrence.
Every input, rejection and exceptional length is defined by family.py.

## 2. Prefix correctness, uniform P and NP

In an algebraic closure of F_2, Z^(2^k)-Z has 2^k distinct roots because
its derivative is 1. The roots are closed under addition, multiplication
and nonzero inversion, hence form a field of size 2^k. Frobenius z->z^2
permutes them and its k-th iterate is the identity. A Frobenius orbit of
length d gives an irreducible degree-d polynomial over F_2: its product
of linear factors is Frobenius-fixed, while a proper factor over F_2 would
require a proper invariant subset of the orbit. Conversely, an irreducible
factor's roots form one orbit. Thus factor degrees divide k.

Since k is prime, the orbits have size 1 or k. Only 0 and 1 are fixed,
and 2^k>2, so an irreducible degree-k polynomial exists for every such k.
No algorithm selecting it is needed in the uniform evaluator.

The first prefix test says Q divides the squarefree polynomial Z^(2^k)-Z;
the second excludes the factors Z and Z+1. Remaining irreducible factors
all have degree k. Since Q has degree k, it must itself be irreducible.
Conversely, every irreducible degree-k Q passes both tests. The test is
therefore exact. Compute its first part by k squarings modulo Q and its
second by polynomial Euclid, both in polynomial bit time.

All degree and order loops have at most N^2 iterations, with trial primality
in polynomial time. Parameters have O(log(N+1)) bits, and every operation
in selecting them is polynomial in N. The selected k,r,t are all <=N.
Schoolbook field multiplication costs O(k^2) bit operations; all powers
use O(r^2) multiplications and the symmetric recurrence O(rt). The exact
prefix test is polynomial as just proved. A coarse O((N+1)^6) bound covers
selection, field arithmetic and rejected inputs. No practical runtime for
all large N, or effective short prime-gap rate, is presumed.

Let A be this single evaluator and let L={x:A(x)=1}. Define the NP verifier
to accept exactly empty witnesses w and inputs with A(x)=1. If x is in L
the empty witness is accepted; conversely acceptance implies x in L. This
proves both directions with one uniform polynomial-time verifier. Witnesses
are not target-circuit inputs.


## 3. The new affine-dimension theorem

Sections1--7 of disperser.md reproduce the field argument, including
prime-degree image-polynomial zero patterns, trace-coefficient isolation,
triangular affine decomposition and weighted-capacity noncancellation.
Section9 converts the same isolated monomial into Boolean algebraic
degree at least B=t(t+3)/2. The old positive-deficit threshold remains a
valid fallback. The new Section10 proves the smaller threshold under the
explicit guards above; it changes only the subset-selection step.

For clarity, the new combinatorial lemma concerns dimensions d_i of
triangular field blocks and deficits delta_i=d_i-w_i, w_i=r+2-i.
Writing S={i:delta_i>0}, a suitable positive t-subset I is one with
sum_(i in I)delta_i+max_(i in I)d_i>k+1. When |S|>=t+1, mark one of the
highest-weight h_|S| indices, then mix t-subsets containing either one
or two of those high-weight indices. All deficits have equal expected
coefficient(t+1)/|S|, while the marked w has larger average than the
uniform mark. The complete probability distribution, boundary cases,
ceiling calculation and monotonicity proof are given in Section10.
Cases |S|<t and |S|=t are handled separately, using respectively T and
the guard2k>=(t+1)^2. No zero-deficit block is inserted into the source
coefficient-isolation lemma.

The resulting integer d strictly exceeds the applicable real threshold.
Both selection routes therefore produce the SAME nonzero coefficient
and degree guarantee on EVERY affine payload subspace of dimension d.
The guard k>r+1 ensures all exponent weights are below k. This applies
for any valid degree-k prefix Q and the nonzero constant-coordinate map.

## 4. Full B2 and the exact finite lower bound

For each selected pair choose any irreducible degree-k Q, whose existence
and exact prefix test were proved in Section2. Fixing this actual k-bit
prefix adds no gates. Padding the m payload bits is injective and affine,
so every dimension-d affine payload space maps to a dimension-d affine
space in K^r. The payload consequently has degree at least B on every
such space. The whole prefix language is not claimed to be an affine
disperser: invalid prefixes have zero regions.

Section6 of elimination.md proves, for1<=B<=d,

    C_B2(g)>=4m-8d+min(5B-1,4B+3),

from degree>=B on every d-dimensional affine subspace. It retains every
binary truth table, arbitrary sharing/fanout, constants, unary and repeated
pins, charged sparse affine-system rows, all terminal cases and boundary
values. Its effective-input parity argument and connectivity proof are
supplied, not treated as an external premise or a free affine basis.
For t>=2, B>=5 and the offset is4B+3. Also B<=H<d because t<r.
Thus a selected admissible pair yields

    C_B2(F_N)>=4(N-k)-8d+4B+3=4N-L_N.

For a length with no admissible pair, F_N=0 and L_N=4N, so the same
finite assertion reduces to C_B2(F_N)>=0 and remains true. Witnesses
are not target-circuit inputs and all paid prefix bits are counted in N.

## 5. Every sufficiently large actual integer length

Let b be the least positive integer with b^5>32N^3 and k the least prime
>=b, solely as an upper reference for the global search. Bertrand gives
k=Theta(N^(3/5)) at EVERY sufficiently large integer N. For m=N-k,
r=ceil(m/k) and t=ceil(sqrt(r)), we have r=Theta(N^(2/5)) and
t=Theta(N^(1/5)). Eventually k>r+1,r>=3,2<=t<r. The old threshold is
O(r^2+kt+rk/t)=O(N^(4/5)), while B=O(N^(2/5)) and m=Theta(N).
So this reference pair is old-admissible at every sufficiently large N.

Section10's identity gives F<=U/(t+1), hence d<=d_old under the guards;
fallback equality holds otherwise. B<=H implies B<=d, so every
old-admissible pair is new-admissible. The reference therefore proves
nonemptiness and L_N<=O(N^(4/5)) at EVERY sufficiently large N.
Admissible losses are positive since B<=H and d>=H+1, and the exceptional
loss4N is nonnegative. Thus L_N=O(N^(4/5))=o(N) eventually.
For each epsilon>0, choose N_epsilon such that L_N/N<=epsilon for all
integer N>=N_epsilon. Then C_B2(F_N)>=(4-epsilon)N at all those lengths,
proving the exact joint coefficient4/exponent1 declaration.
No PNT, unproved short-interval prime statement, sparse subsequence,
parameter limit optimization, or a family varying with epsilon is used.

## 6. Compare complete certificates, not actual complexities

Define the original reference search by the SAME prime k, r and order t
loops and guards, but using d_old for admissibility and loss
ell_old=4k+8d_old-4B-3. Whenever this old search is nonempty, let L_old(N)
be its minimum loss with the same lexicographic tie rule. The original
selection is fully specified; no prior accepted theorem is a premise.

Every old-admissible pair remains available in the new search and its
loss does not increase. Hence L_N<=L_old(N) whenever the old search is
nonempty, in particular at every sufficiently large actual N. These
are proved losses of two fully specified parameter-selected languages,
not an ordering of their actual minimum circuit sizes.

The exact finite-audit experiment enumerates EVERY prime k and EVERY
integer t for the following finite comparisons with independent literal
capacity sums and rational weighted thresholds:

    N       L_old       L_N       saving
    200       669        661          8
    1000     2533       2493         40
    10000   16809      16601        208

These finite oracle statements do not assert a growing strict saving
at every length or an asymptotic improvement of the coefficient. The
new general threshold in Section10, rather than finite evidence alone,
is the mathematical improvement. Its exact saving in the real second
threshold is r*(r-h)/(2(t+1)); max/floor and global selection may reduce
or eliminate that saving at a particular pair or length.

## 7. Provenance and evidence boundaries

The retained algebraic disperser and full-B2 degree-sensitive argument,
field arithmetic and global minimum selection come from source submission
8eef1799-c27f-46f3-8e61-3c4e0d28d145, frozen commit
b4a12ebc5267f63705e1c277565e2779c088c286. The compact presentation,
uniformity explanation and supplementary audit also build on our
2309de11-70f3-4d9d-99be-be46178f7983, frozen commit
e02200ac4e09415c05ce384d695752641568e129. Both are material sources.
Their judgments are not mathematical premises. Underlying literature is
BK's algebraic construction and Mishra--Santhanam ECCC TR26-260 revision1
(November3,2026); all unfamiliar arguments needed here are supplied.

The current new work is the weighted marked-subset distribution,
three-case positive-set proof, guarded exact dimension, integration with
the same degree witness, and the resulting smaller complete certificates.
We claim no new leading coefficient, P/poly separation or equality with
actual circuit complexity. Familiar standard field algebra, Gaussian
elimination, root counting and Bertrand are used as specified.
The candidate finite-audit checks exact finite predicates, including the
marked distribution and selection theorem on bounded dimension vectors.
It does not establish universal coverage or infinite asymptotics; those
rest on the written proofs. No simulated result or copied judgment is used.
