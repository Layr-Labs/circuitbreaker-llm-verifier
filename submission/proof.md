# Exact integer block selection and complete coefficient-four certificates

The new result exactly maximizes the integer dimension vectors for which
the retained positive-deficit isolation predicate fails. This gives a
smaller sufficient dimension, with no extra weighted-selection guards.
The same field witness and degree-sensitive full-B2 theorem apply. Global
selection never increases the original or guarded-weighted certified loss.
The leading coefficient/exponent remains(4,1); no actual circuit-size ratio
or new leading coefficient is claimed.

## 1. One exact family, at every length

At actual integer N>=0 enumerate every prime2<=k<N. Put m=N-k and
r=ceil(m/k), requiring r>=3 and k>r+1. Enumerate every2<=t<r. Define

    H=r(r+3)/2, T=(t-1)k-t(t+1)/2+1, B=t(t+3)/2.

For every positive-set size s=t,...,r and integer A=t-1,...,k-s-2 put

    X=min(floor(A/(t-1)),floor((k-t-A)/2),
          floor(((t-1)k-t(t-1)/2-tA)/(t-1))).

Skip choices with X<1. For the others let C=k-A,J=min(s,C-2X), and

    e_(s,A)=A+(J-t+1)X+floor((C-J-1)^2/4)-floor((C-s-1)^2/4).

Let E be the maximum of T and all these e values (empty ranges add
nothing), and d=H+E+1. A pair is admissible when B<=d<m. Select the
lexicographically least(4k+8d-4B-3,k,t). If none exists, output0 at the
whole length and set L_N=4N; otherwise L_N is the selected loss. In
particular N=0 is defined. All floors, loops and tie rules are exact.

The code's safe shortcut uses the previous guarded threshold: if it equals
H+T+1, Section11 proves E=T, so the full scan is unnecessary. This is an
equality shortcut, not a heuristic or an alternate family definition.
It also skips H+T+1>=m, which cannot be admissible since E>=T.

On an admissible length, read the first k actual bits as lower coefficients
of Q(Z)=Z^k+sum_(j=0)^(k-1)b_j Z^j overF_2. Require Z^(2^k)=Z modQ
and gcd(Q,Z^2+Z)=1; otherwise output0. For a passing prefix use the
field K=F_2[Z]/Q. Pad the m payload bits by rk-m zeros into r field
blocks a_i in the polynomial basis, first bit the constant coordinate.
Put y_i=a_i^(2^(r+2-i)-1), and output the constant coordinate of the
order-t elementary symmetric polynomial in the y_i. Its recurrence starts
E_0=1 and other E_j=0 and updates E_j<-E_j+y_i*E_(j-1) in descending
j=min(i,t),...,1. Splitting subsets by membership of i proves this
recurrence without division. family.py defines every input and rejection.

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

There are at most N^2 degree/order pairs, with trial primality in polynomial
time. For each pair the exact selection threshold scans at most r*k<=N^2
integer choices; its proven shortcut only reduces this work. Parameters have O(log(N+1)) bits, and every operation
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



## 3. Why the exact integer threshold is sufficient

The field proof in disperser.md Sections1--7 reduces noncancellation to
a positive t-subset I with sum_I delta_i+max_I d_i>k+1, where
w_i=r+2-i and delta_i=d_i-w_i. Section9 turns the same isolated
monomial into Boolean degree at least B=t(t+3)/2. The strict inequality
is the unchanged source predicate; no zero-deficit padding or additional
field premise is inserted.

The new Section11 proves that H+E is EXACTLY the largest dimension of
an integer vector0<=d_i<=k failing that predicate. The main reduction
moves positive deficits to the smallest weights and sorts them in the
opposite order. For a fixed deficit multiset, g_j=2delta_j plus the sum
of the largest t-1 other deficits is nondecreasing in delta_j. Pairing
larger g with smaller weights cannot increase the largest failed
inequality. Both dimension caps and the positive sum are preserved.

For s>=t positive indices, write A for the top t-1 sum and x for the
next deficit. All failed subset inequalities reduce to the explicit top
upper bounds and descending tail in Section11. Feasible integer top
vectors exist for exactly the stated range1<=x<=X, using a constructive
left-to-right fill. For fixed A the maximal tail increases with x, so X
is optimal. The floor-sum identity gives e_(s,A) exactly. Every listed
value has an attaining failed vector, as does T for s<t. The proof
therefore covers all integer vectors, not just sampled dimensions or
real LP points. The asserted optimality concerns ONLY this unchanged
selection predicate, not the polynomial, circuits, or all proof methods.

Every affine payload dimension at least d=H+E+1 supplies a good subset.
The field coefficient-isolation argument, triangular mixing and
noncancellation are unchanged; degree>=B follows on every such space.
The field guard k>r+1 ensures all exponent weights are below k.

## 4. Full B2, paid prefix and the finite inequality

Choose an irreducible degree-k Q, whose existence was proved in Section2,
and fix that actual k-bit prefix in an arbitrary circuit. Constant
normalization cannot add gates. Padding the remaining m actual bits is
an injective affine map, preserving the dimension of every payload space.
Thus the payload has degree>=B on every dimension-d affine space.
The whole prefix language is not asserted to be an affine disperser;
rejected prefixes have zero regions.

The supplied elimination.md Section6 proves that this degree hypothesis,
1<=B<=d, gives C_B2(g)>=4m-8d+min(5B-1,4B+3). Its argument covers
all sixteen binary tables, arbitrary sharing/fanout, constants, unary and
repeated pins, charged fair affine systems, terminal outputs, and boundary
cases. It does not supply an uncharged affine computational basis.
Here t>=2 gives B>=5 and offset4B+3, while t<r gives B<=H<d. Hence

    C_B2(F_N)>=4(N-k)-8d+4B+3=4N-L_N.

On an empty search F_N=0,L_N=4N and the same assertion is the true zero
bound. All prefix bits are counted in actual N; witnesses are not inputs
to target circuits.

## 5. Every sufficiently large actual integer length

Let d_old=H+max(T,floor(U/(t+1)))+1, U=rk-r(r+1)/2. Let d_weighted
be the guarded Section10 dimension, with its old fallback. Both are
proved sufficient for the same subset predicate. The attaining failed
vector of Section11 therefore gives H+E<d_old and H+E<d_weighted.
Thus d<=d_weighted<=d_old for EVERY eligible k,r,t. This uses the
supplied proofs, not a prior platform judgment.

For an upper reference at each N, take the least prime k at least the
least integer b with b^5>32N^3, and set t=ceil(sqrt(r)). Bertrand gives
k=Theta(N^(3/5)),r=Theta(N^(2/5)),t=Theta(N^(1/5)) at every sufficiently
large integer N. Eventually k>r+1,r>=3,2<=t<r, and
B=O(N^(2/5))<=d_old=O(N^(4/5))<N-k. Since B<=H<d<=d_old, that pair
is new-admissible. The search is nonempty at EVERY sufficiently large N,
and its minimum L_N is at most O(N^(4/5)). Losses are nonnegative because
B<=H and d>=H+1; hence L_N=O(N^(4/5))=o(N) eventually.
For each epsilon>0 some N_epsilon makes L_N/N<=epsilon at EVERY integer
N>=N_epsilon, proving C_B2(F_N)>=(4-epsilon)N there. One fixed total
family and verifier work throughout. No PNT, sparse subsequence or
parameter-limit optimization is needed for the joint pair(4,1).

## 6. Complete certificate comparisons

The original reference search uses the SAME prime/order loops and base
guards, but d_old in admissibility and loss4k+8d_old-4B-3. The weighted
reference similarly uses d_weighted from Section10. Both tie rules are
lexicographic in(loss,k,t), with all definitions given above and in
family.py. Write their minima L_old,L_weighted whenever the respective
search is nonempty. Each reference pair remains new-admissible with no
larger loss, so L_N<=L_weighted and L_N<=L_old whenever the respective
reference exists, including every sufficiently large integer length.
These compare certified losses of parameter-selected languages, not
actual circuit complexities.

Exact independent finite oracles enumerate every prime and order and give

    N      L_old   L_weighted   L_N
    500     1461       1437     1429
    1000    2533       2493     2485
    2000    4513       4449     4433
    10000  16809      16601    16577.

The declared finite-audit validates these finite comparisons and the exact
integer threshold against exhaustive small dimension vectors. The written
rearrangement, feasibility and attainment proof establishes the universal
threshold. No claim asserts an additional strict gain at every length,
a new leading coefficient, or an actual circuit-size difference.

## 7. Provenance and evidence scope

The algebraic construction, field-prefix evaluator, global selector and
full-B2 degree-sensitive theorem are retained from8eef1799-c27f-46f3-8e61-3c4e0d28d145
(commit b4a12ebc5267f63705e1c277565e2779c088c286). The compact presentation
and audit build on2309de11-70f3-4d9d-99be-be46178f7983 (commit
e02200ac4e09415c05ce384d695752641568e129). The guarded marked-subset
argument, its exact reference selector and finite audit are retained from
8d0d93b5-9da6-405e-98f3-aef47d9b9c09, whose supplied source is the
previous seven-file package. These are material providers, not premises
based on their judgments. Underlying literature is BK and Mishra--Santhanam
ECCC TR26-260 revision1, November3,2026; required unfamiliar arguments are
supplied. Familiar field algebra, Gaussian elimination, root counting and
Bertrand are used as stated.

The new contribution is Section11's exact integer reduction, feasible top
and tail search, attaining vectors, safe fast path, and integration into
one total family. It gives optimality ONLY of the unchanged positive
isolation predicate. It neither bounds the power of other field arguments
nor asserts optimality of actual circuits. Finite experiments establish
only their declared predicates. Universal eligibility and asymptotics
rest on these written proofs, without copied judgments or simulation.
