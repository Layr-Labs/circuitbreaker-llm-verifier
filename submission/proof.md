# Exact integer block selection and complete coefficient-four certificates

The new result permits zero-deficit blocks to complete the isolated
coefficient, and exactly maximizes the integer dimensions failing this
stronger predicate. This improves finite certified losses while retaining
the full-B2 theorem and leading coefficient/exponent(4,1).

## 1. One exact family, at every length

At actual integer N>=0 enumerate every prime2<=k<N. Put m=N-k and
r=ceil(m/k), requiring r>=3 and k>r+1. Enumerate every2<=t<r. Define

    H=r(r+3)/2, T=(t-1)k-t(t+1)/2+1, B=t(t+3)/2.

For every positive-set size s=t,...,r and integer A=t-1,...,k-s-2 put

    X=min(floor(A/(t-1)),floor((k-t-A)/2),
          floor(((t-1)k-t(t-1)/2-tA)/(t-1))).

Skip choices with X<1. For the others let C=k-A,J=min(s,C-2X), and

    e_(s,A)=A+(J-t+1)X+floor((C-J-1)^2/4)-floor((C-s-1)^2/4).

For each integer1<=s<t put

    a_s=min(k-s-1,floor((2sk-s(s+1))/(2(s+1)))),
    M_0=max(0,T-r+t-1).

Let E0 be the maximum of M_0, all a_s with a_s>=s, and all e_(s,A)
above. Empty ranges add nothing. Set d=H+E0+1. A pair is admissible when
B<=d<m. Select the lexicographically least(4k+8d-4B-3,k,t). If none
exists output0 at the whole length and set L_N=4N; otherwise L_N is the
selected loss. In particular N=0 is defined. All loops and floors are exact.
Section12 proves the skip H+M_0+1>=m is safe because E0>=M_0.
The code also skips a pair if its lower loss key
(4k+8(H+M_0+1)-4B-3,k,t) is at least the best already selected key.
The true key is no smaller, so it cannot improve that incumbent. This
changes neither the exact minimum nor its lexicographic tie rule.

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
integer choices; its proven skip only reduces this work. Parameters have O(log(N+1)) bits, and every operation
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

Section12 of disperser.md proves the new field premise. Select t blocks
with delta_i=d_i-w_i>=0, with at least one positive deficit. The required
inequality is sum_positive delta_i+max_positive d_i>k+1. Section3's
explicit isolation construction on positive sources kills every nonzero
trace shift. Each zero-deficit source has a full-row-rank matrix by the
linearized-polynomial root bound, so can contribute a nonzero minor at
the surviving shift. Triangular mixing, symmetric noncancellation and
Boolean degree transfer apply unchanged. In particular the maximum in
the predicate is over POSITIVE sources, not arbitrary selected sources.

Section12 also proves that H+E0 is exactly the largest total dimension
failing this enlarged predicate. Fewer than t eligible blocks have maximal
excess T-r+t-1. With at least t eligible but s<t positive blocks, all
positives can be padded by zeros; upper bounds on individual deficits and
their sum give exactly a_s, with a constructive integer fill proving
attainment. For s>=t, replacing selected zeros by unused positives shows
equivalence to the old positive-only predicate; Section11 supplies the
exact feasible top/tail values e_(s,A). All branches have explicit
attaining vectors, and together exhaust all integer vectors0<=d_i<=k.
Optimality is asserted ONLY for this enlarged selection predicate.

Consequently every affine payload dimension at least d=H+E0+1 supplies
a good selection, and has Boolean degree>=B=t(t+3)/2. The guard k>r+1
ensures every selected weight is below k. All unfamiliar field and
elimination arguments used by the theorem are supplied in the package.

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
proved sufficient for the same subset predicate. Section11 gives the old exact threshold H+E+1 no larger than either
reference. Section12 gives E0<=E because its good predicate includes the
old one. Thus d<=H+E+1<=d_weighted<=d_old for EVERY eligible k,r,t. This uses the
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

    N      L_old   L_weighted   L_positive-exact   L_N
    500     1461       1437          1429          1413
    1000    2533       2493          2485          2485
    2000    4513       4449          4433          4417
    10000  16809      16601         16577         16521.

The positive-exact reference uses H+E+1 from Section11, otherwise the same
complete loops and lexicographic tie rule. Section12 gives L_N no larger
whenever that reference exists. These are certified-loss comparisons,
not actual circuit complexities. The finite audit checks the comparisons
and bounded dimension vectors; the written proof establishes universal
thresholds. No strict improvement at every length is asserted.

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

Section11's exact positive-only optimization and its finite audit are
retained from673f938b-d06e-4423-9514-ed84afe6b13f. The new contribution
is Section12's zero-deficit coefficient extension, exhaustive integer
threshold for that predicate, and integration into this one total family.
The exactness assertion concerns this predicate alone, not actual circuit
optimality or limits of other field arguments. Finite checks prove only
their declared bounded predicates. Universal eligibility and asymptotics
rest on the supplied written proofs, without simulation or copied judgments.
