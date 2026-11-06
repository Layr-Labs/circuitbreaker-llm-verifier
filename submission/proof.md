# Five-bit field-prefix compression and complete coefficient-four certificates

The new result encodes the field polynomial using k-5 actual bits.
An elementary trace-root argument and reciprocal normalization prove
that every eligible degree k has an irreducible polynomial in the encoded
class. The zero-deficit field theorem and full-B2 elimination theorem are
retained. The leading coefficient/exponent remains(4,1); finite certified
losses improve, with no assertion about actual optimal circuit sizes.

## 1. One exact family, at every length

At actual integer N>=0 consider each integer saved=0,...,5. Enumerate
every prime2<=k<N+saved, put m=N+saved-k and r=ceil(m/k), and require
r>=3,k>r+1. Consequently every eligible k is an odd prime at least5.
Enumerate every2<=t<r. Put ignored=5-saved and prefix length
p=k-5+ignored=k-saved. Define

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
B<=d<m. Among all these choices select the least key
(4p+8d-4B-3,k,t,ignored), together with a constant-zero option whose
key is(4N,0,0,0). On the zero option output0 and set L_N=4N; otherwise
L_N is the selected loss. In particular N=0 is defined. The code realizes
this minimum by finding the previous zero-deficit minimum separately at
each reference length N+saved. Within a fixed saved value the reference
loss and new loss differ by the constant4saved, preserving its tie rule.
Its proven dimension and incumbent lower-bound pruning does not change
those reference minima. All loops and floors are exact.

On a nonzero-option length read k-5 actual bits b_j, j=0,...,k-6.
Set R(Z)=Z^(k-4)+sum_j b_j Z^j and Q(Z)=1+(Z^4+Z^2)R(Z).
Ignore the next ignored actual bits; they are still counted in p and N.
Require Z^(2^k)=Z modQ and gcd(Q,Z^2+Z)=1; otherwise output0.
For a passing prefix use the field K=F_2[Z]/Q. Pad the m payload bits by rk-m zeros into r field
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

There are O((N+1)^2) degree/order pairs across the six reference lengths, with trial primality in polynomial
time. For each pair the exact selection threshold scans at most r*k<=O((N+1)^2)
integer choices; its proven skip only reduces this work. Parameters have O(log(N+1)) bits, and every operation
in selecting them is polynomial in N. The selected k,r,t are all <=N+5. Decoding R and Q is linear in k.
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



### 2.1 Existence in the compressed polynomial class

Let k>=5 be odd prime, K a field of2^k elements as above, and
Tr(z)=sum_(i=0)^(k-1)z^(2^i). Tr takes values in F_2 and Tr(1)=1,
so its kernel has M=2^(k-1) elements. There exists nonzero alpha with
Tr(alpha)=Tr(alpha^-1)=0. Suppose otherwise. The monic polynomial

    T(X)=Tr(X)/X=sum_(i=0)^(k-1)X^(2^i-1)

has precisely the M-1 distinct nonzero trace-zero elements as roots.
All would then have inverse trace1. Thus T divides the polynomial

    U(X)=X^M(Tr(1/X)+1)=X^M+sum_(i=0)^(k-1)X^(M-2^i).

U is monic degree M, while T has degree M-1, so the quotient is X+c.
The coefficient of X^(M-1) forces c=1. But U has coefficient1 at
X^(M-2), while (X+1)T has coefficient0 there. Its exponents are powers
of2 or powers of2 minus1; M-2 is neither when M>=8. This contradiction
proves the claim using polynomial division and the elementary root bound.

Neither alpha nor alpha^-1 is1 since Tr(1)=1. Replace alpha by alpha^-1
if necessary to arrange also Tr(1/(alpha+1))=0. This is always possible:

    Tr(1/(alpha^-1+1))=Tr(alpha/(alpha+1))
                       =1+Tr(1/(alpha+1)).

The first two zero traces are preserved under inversion. Since k is
prime and alpha is not in F_2, its Frobenius orbit has size k. Its monic
irreducible polynomial Q has constant1, coefficient of Z^(k-1) equal
Tr(alpha)=0, and coefficient of Z equal Tr(alpha^-1)=0. These follow by
expanding the product over the conjugates; the product of conjugates is
nonzero in F_2, hence1. Also Q(1)=1 since Q has no root1. The logarithmic
derivative identity, obtained by differentiating that product, gives

    Q'(1)/Q(1)=sum_(i=0)^(k-1)1/(1+alpha^(2^i))
              =Tr(1/(1+alpha))=0.

Consequently Q+1 is divisible by Z^2 (its constant and linear coefficients
vanish) and by (Z+1)^2 (its value and derivative at1 vanish). These
coprime factors give

    Q=1+Z^2(Z+1)^2 R=1+(Z^4+Z^2)R.

R is monic degree k-4. Its coefficient at degree k-5 is0 because Q's
coefficient at degree k-1 is0. All remaining coefficients are exactly
k-5 free bits, as in the decoder. Hence every eligible k has at least
one encoded irreducible polynomial. Multiplication by Z^4+Z^2 is
injective in F_2[Z], so different bit prefixes encode different Q. The
prefix test from Section2 remains exact for every decoded monic Q.
This proof constructs no field polynomial by advice or an external
algorithm: it proves existence for restricting circuits; evaluation just
decodes its actual input and checks it in polynomial time.

## 3. Why the exact integer threshold is sufficient

For the selected prime k, the retained Section12 of disperser.md proves the field premise. Select t blocks
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

Choose an encoded irreducible degree-k Q using Section2.1,
and fix its actual k-5 encoded prefix bits and all ignored bits in an arbitrary circuit. Constant
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

    C_B2(F_N)>=4(N-p)-8d+4B+3=4N-L_N.

On an empty search F_N=0,L_N=4N and the same assertion is the true zero
bound. All prefix bits are counted in actual N; witnesses are not inputs
to target circuits.

## 5. Every sufficiently large actual integer length

Use the saved=0 branch to recover every zero-deficit certificate at
actual length N: its prefix has k-5 encoded bits and5 ignored bits, so
p=k and payload m=N-k. The zero-deficit argument supplies d no larger
than the old exact, guarded weighted and original sufficient thresholds.
For an all-length upper reference take the least prime k at least the
least integer b with b^5>32N^3, and t=ceil(sqrt(r)). Bertrand gives
k=Theta(N^(3/5)),r=Theta(N^(2/5)),t=Theta(N^(1/5)) for every sufficiently
large integer N. Eventually all guards hold and B<=H<d<=d_old=
O(N^(4/5))<N-k. Thus the selected loss is at most O(N^(4/5)). All
nonzero-option losses are nonnegative since B<=H<d and p>=0; the zero
option is nonnegative as well. Therefore L_N=o(N) over ALL integer
lengths. For every epsilon>0, C_B2(F_N)>=(4-epsilon)N for all large N.
The total evaluator and empty-witness NP verifier are fixed throughout.

## 6. Complete certificate comparisons

Let L0(n) denote the complete zero-deficit reference loss from the preceding family and disperser.md Section12: it reads an unconstrained k-bit
polynomial prefix, uses m=n-k, minimizes(loss,k,t), and uses4n on an
empty search. The original, guarded-weighted and positive-exact references
replace its d by their respective previously defined dimensions. All
reference searches enumerate every eligible prime and order.

The new loss has the exact identity

    L_N=min(4N, min_(s=0,...,5)(L0(N+s)-4s)).

Indeed a nonempty reference search at N+s has precisely the new payload
and field/order guards with p=k-s, and its loss decreases by4s. An empty
reference contributes the constant-zero bound4N. The explicit zero
option covers this value even if every reference search is nonempty.
Ties between nonzero choices are in(loss,k,t,ignored), and the zero option
wins any tie at4N. In particular L_N<=L0(N), and also no larger than each
older reference whenever that reference exists.

Exact complete searches give the following finite certified losses:

    N      original   weighted   positive-exact   L0(N)   L_N
    500      1461       1437          1429         1413    1393
    1000     2533       2493          2485         2485    2465
    2000     4513       4449          4433         4417    4397
    10000   16809      16601         16577        16521   16501.

These are certificate comparisons between parameter-selected languages,
not actual optimal circuit sizes. Finite checks reproduce the reference
searches and compression predicates; the written proofs establish
universal correctness. No strict improvement at every length is asserted.

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

Section11's exact positive-only optimization and finite audit are retained
from673f938b-d06e-4423-9514-ed84afe6b13f. Section12's zero-deficit field
extension, exact three-branch threshold and witness audit are retained
frombf5b6832-24d7-4ceb-9d64-abd59ca97414. Neither a prior verdict nor promotion is
used as a premise. Exactness concerns only the stated finite predicate,
not actual circuit optimality or the limits of other field arguments.
The five-bit compressed class, elementary trace-root existence proof,
affine polynomial decoder and six-length selector are the new contribution.
Finite experiments establish only their bounded predicates; universal
eligibility and asymptotics rest on the supplied written proofs.
