# One uniform family selected by the best supplied certificate

## Claim and scope

For one fixed total P language at actual input length N, we prove

    C_B2(F_N)>=4N-L_N,
    L_N/N^(4/5)->c_star=10*2^(1/5).

Thus its native joint coefficient/exponent is (4,1). All input bits,
including its field-polynomial prefix, are counted. The current contribution
is global exact parameter selection, a two-gate finite refinement, an
additional infinite-length comparison to the preceding two-prime selector,
and optimality of the LEADING CERTIFICATE among all admissible parameters.
This is not optimality of actual circuit complexities or other constructions.
No restricted basis, direct-sum premise, conjecture or exponent gain is used.

The supplied algebraic proof follows Ben-Sasson--Kopparty's construction.
The full-B2 argument follows Mishra--Santhanam, ECCC TR26-260 revision1,
November3,2026, with the supplied later refinements. All necessary unfamiliar
arguments appear in disperser.md and elimination.md. Previous platform
judgments are not mathematical premises. Section 9 identifies material
providers; the proof is self-contained apart from familiar standard facts
explicitly stated below.

## 1. Exact admissible parameters and evaluator

For each nonnegative integer N enumerate every prime integer 2<=k<N.
For that degree put m=N-k and r=ceil(m/k). Require r>=3 and k>r+1.
For every integer order 2<=t<r define

    H=r(r+3)/2,
    T=(t-1)k-t(t+1)/2+1,
    U=rk-r(r+1)/2,
    d=H+max(T,floor(U/(t+1)))+1,
    B_deg=t(t+3)/2.

Call (k,t) admissible if B_deg<=d<m as well. Define its COMPLETE loss

    ell_N(k,t)=4k+8d-4B_deg-3.

If there are no admissible pairs, output0 on every length-N input.
Otherwise choose the lexicographically least (ell_N(k,t),k,t), and write
its parameters as k,m,r,t,d_N,B_deg and its loss as L_N. This is a finite
exact integer minimization; no floating point or unspecified tie is used.

Read the first k bits b_0,...,b_(k-1) as the lower coefficients of

    Q(Z)=Z^k+sum_(j=0)^(k-1) b_j Z^j in F_2[Z].

Check Z^(2^k)=Z modulo Q and gcd(Q,Z^2+Z)=1. If either fails, output0.
For valid Q use K=F_2[Z]/Q. Zero-pad the remaining m payload bits to rk
bits, grouping them in r field blocks a_i, with the first bit of each
block its constant coordinate. For i=1,...,r put w_i=r+2-i and

    y_i=a_i^(2^(w_i)-1).

Output the constant field coordinate of sum_(|I|=t) product_(i in I)y_i.
Initialize E_0=1,E_1=...=E_t=0 and update E_j<-E_j+y_i*E_(j-1) for
j=min(i,t),...,1 in decreasing order. Splitting subsets by membership of i
proves this division-free recurrence. This defines every bit string and
length, including N=0, where the admissible set is empty. Q is input, not
advice or the output of an unspecified irreducible-polynomial algorithm.

The code computes the capacity sums with their clipped closed forms,
which equal the displayed formulas under k>r+1 and t<r. Explicitly, with
u=min(t,k-1),v=min(r,k-1), T=(u-1)k-u(u+1)/2+1 and U=vk-v(v+1)/2.
All candidate tests and scores therefore use exact polynomial-size integers.

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

## 3. A finite full-B2 lower bound for every admissible pair

The precise general lemma in disperser.md proves that for any prime k,
1<=t<r and nonzero F_2-linear pi:K->F_2,

    pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

has algebraic degree at least B_deg=t(t+3)/2 on every d-dimensional
affine subspace, where d is the Section 1 threshold. It supplies the
image-polynomial zero patterns, trace-coefficient isolation, full-hypothesis
subset averaging, triangular affine decomposition and weighted-capacity
noncancellation, then converts binary monomial weight to Boolean degree.
The admissibility guard k>r+1 also ensures all w_i<k.

For any admissible (k,t), fix an irreducible degree-k Q, whose existence
was proved in Section 2, in the actual k-bit prefix. Zero padding is an
injective linear map of the m actual payload bits. Thus its payload G_Q
has degree at least B_deg on every d-dimensional affine payload space.
The whole prefix language is not asserted to be an affine disperser;
invalid prefixes have zero regions. Fixing a valid prefix adds no gates
and every prefix bit remains charged in N.

Section 6 of elimination.md strengthens the supplied finite theorem to

    C_B2(f)>=4n-8d+min(5B-1,4B+3)

under degree(f|A)>=B for every d-dimensional affine A,1<=B<=d. It retains
all sixteen binary truth tables, arbitrary sharing and fanout, charged
fair affine-system rows, constants, unary/degenerate gates, repeated pins
and all output cases. The new effective-input parity refinement and its
boundary cases are supplied there; no basis or formula transfer is assumed.
Here t>=2 gives B_deg>=5, so the offset is 4B_deg+3. Consequently every
admissible branch has

    C_B2(F_(N;k,t))>=4(N-k)-8d+4B_deg+3
                     =4N-ell_N(k,t).

The same proof applies to the selected branch, giving C_B2(F_N)>=4N-L_N
whenever the admissible set is nonempty. This is a finite statement, not
only an asymptotic resource declaration.

## 4. All sufficiently large lengths have admissible parameters

For comparison define a reference upper prime k_plus as follows: b is
the least positive integer with b^5>4N^3, and k_plus is the least prime
>=b. Let r_plus=ceil((N-k_plus)/k_plus) and t_up=ceil(sqrt(r_plus)).
Bertrand's prime theorem gives k_plus=Theta(N^(3/5)), hence
r_plus=Theta(N^(2/5)),t_up=Theta(N^(1/5)). For every large integer N,
k_plus>r_plus+1, r_plus>=3,2<=t_up<r_plus. The Section 1 dimension is
O(N^(4/5)), B_deg=O(N^(2/5)) and m=Theta(N), so B_deg<=d<m eventually.
This reference pair is therefore admissible at EVERY sufficiently large
length. Its loss is O(N^(4/5)), proving the selected minimum exists and
L_N=O(N^(4/5)) at all such lengths. The finite theorem already implies
liminf C_B2(F_N)/N>=4. No sparse sequence is used for the principal bound.

For the exact leading constant we use the classical unconditional prime
number theorem in its standard form pi(x)~x/ln(x). This familiar theorem
is not re-proved; its exact statement is NIST DLMF27.2.3,
https://dlmf.nist.gov/27.2.E3 (Hadamard/de la Vallee Poussin). For fixed
eta>0 the difference pi((1+eta)x)-pi(x)~eta*x/ln(x) is positive eventually.
Thus the least prime >=x is x*(1+o(1)). Similarly the largest prime <x is
x*(1+o(1)), from the interval ((1-eta)x,x); endpoints change counts by
at most one. These are all-length statements and do not supply a uniform
short-interval prime theorem.

It follows that k_plus/N^(3/5)->a_star=2^(2/5). For a reference degree
k/N^(3/5)->a>0 with t=ceil(sqrt(r)), direct rounding gives

    d=r^2/2+kt+o(N^(4/5)),
    ell_N(k,t)/N^(4/5)->8F(a), F(a)=sqrt(a)+1/(2a^2).

Indeed T=kt+O(k+t^2), U/(t+1)=kt+O(k+r^(3/2)), and r=N/k+O(1).
Also the prefix4k and degree term4B_deg are lower order. Calculus gives
F'(a)=(a^(5/2)-2)/(2a^3), with unique global minimum at a_star. Its value
is(5/4)*2^(1/5). Therefore the reference above proves

    limsup L_N/N^(4/5)<=c_star=10*2^(1/5).

## 5. Global optimality of the leading certificate

We prove that NO admissible parameter sequence has normalized loss with
liminf below c_star, without assuming convergence of its parameter ratios.
This is only optimality of the explicit Section 1 CERTIFICATE, not circuit
complexity, all lower-bound methods, or all algebraic constructions.

Consider any admissible pairs along lengths tending to infinity whose
loss ell is bounded above by C*N^(4/5), for a fixed C. Other sequences
cannot yield a smaller finite liminf unless such a subsequence exists.
Since t<r, B_deg<=H: at t=r-1 their difference is r+1, and smaller t
only reduces B_deg. Write M=max(T,floor(U/(t+1))). Then

    ell=4k+8(H+M+1)-4B_deg-3>=4k+4H+8M+5.          (1)

In particular ell>=4k and ell>=8U/(t+1)-3. Using the exact T in (1),

    ell>=(8t-4)k+2r^2+6r-4t^2-4t+13>=2kt+2r^2. (2)

The last inequality uses t<r<k, so4t^2<=4kt, t>=2 so4t-4>=2t,
and6r-4t+13>0.

These inequalities force the parameter scale, rather than assuming it.
First k=O(N^(4/5))=o(N),r=O(N^(2/5)). Since N-k<=rk<N, rk~N,
hence k=Omega(N^(3/5)). Equation(2) then gives t=O(N^(1/5)). Moreover
U=rk-r(r+1)/2=N+O(k+r^2)~N; the bound ell>=8U/(t+1)-3 forces
t=Omega(N^(1/5)). Returning to(2) gives k=O(N^(3/5)), then
r=Omega(N^(2/5)). Thus

    a=k/N^(3/5), b=t/N^(1/5)

lie in a compact subset of (0,infinity)^2, with constants depending on C.
Uniformly on this set,

    T=kt+O(N^(3/5)),
    floor(U/(t+1))=N/t+O(N^(3/5)),
    d=r^2/2+max(kt,N/t)+O(N^(3/5)).

For the second estimate use U=N+O(k+r^2), divide by t+1, and note
N/t-N/(t+1)=O(N/t^2)=O(N^(3/5)). The remaining errors are no larger.
Prefix and degree terms are lower order, so

    ell/N^(4/5)=G(a,b)+o(1),
    G(a,b)=4/a^2+8max(ab,1/b)>=4/a^2+8sqrt(a)>=c_star. (3)

The first inequality uses max(ab,1/b)>=sqrt(a); equality requires
b=a^(-1/2). The unique minimum in a is the derivative computation in
Section 4. This proves the asserted liminf bound for every admissible
parameter sequence, including ones without convergent ratios.

Apply this to the selected global minimum, whose loss is bounded by the
reference in Section 4. Equations(3) and that reference yield

    L_N/N^(4/5)->c_star.

The minimum of G is unique. Compactness therefore also implies the
selected ratios converge: otherwise a subsequence a positive distance
from that minimizer would have a further compact limit with G>c_star.
Consequently

    k/N^(3/5)->a_star, t/N^(1/5)->a_star^(-1/2),
    r=Theta(N^(2/5)), t/sqrt(r)->1.

No additive estimate t=sqrt(r)+O(1) for the globally selected order is
claimed or needed. This distinction matters for the retained comparisons.

## 6. Fully specified earlier references and uniform comparison

Define F_dual, the preceding two-prime reference, with threshold b from
Section 4, k_plus as above and k_minus the largest prime <b if it exists.
A branch is eligible when k<N and ceil(sqrt(r))<r. If the k_plus branch
is ineligible output0. Otherwise select between all eligible branches and
t in {floor(sqrt(r)),ceil(sqrt(r))} by the lexicographically least
(ell_old(k,t),k,t), where

    ell_old(k,t)=4k+8d(k,t)-4B_deg(t)-1.

Use the same prime-degree prefix test and symmetric recurrence. This is
one fixed total P reference language. At all large lengths both neighboring
degrees are Theta(N^(3/5)), their adjacent orders are Theta(N^(1/5)),
and all the global admissibility filters hold. Thus the chosen reference
pair is included in our search with its loss smaller by exactly2:

    L_N<=L_dual(N)-2 for every sufficiently large integer N. (4)

This is a comparison of complete proved losses for DIFFERENT fixed
languages, not their actual minimum circuit complexities. We next prove
an additional unbounded strict saving, not just this finite two-gate gain.

For completeness define F_adj similarly with only k_plus and its adjacent
orders, and F_up with only k_plus and t_up. The older generic finite bound
in elimination.md Sections 4-5 independently proves their losses with the
same ell_old formula. Inclusion gives L_dual<=L_adj<=L_up at every large
length, so (4) also retains the prior weak comparisons. No monotonicity of
the currently selected global parameters is assumed.

## 7. An additional infinite-length saving over the two-prime reference

For a sufficiently large odd prime q let h<q<p be its preceding and next
primes, and put

    N_q=ceil((q^5/4)^(1/3)), c=4^(-1/3).

Equality4N_q^3=q^5 is impossible by powers of two. The O(1) rounding of
N changes4N^3 by O(q^(10/3)), smaller than(q+1)^5-q^5~5q^4. Therefore
b(N_q)=q+1 for all large q, and F_dual's two degrees are exactly q,p.
Also N_q/q=c*q^(2/3)+O(1/q).

We show infinitely many such actual lengths have the SAME block count
for all three degrees h,q,p. For q in[X,2X], the PNT gives
(1+o(1))*X/log(X) middle primes. All consecutive gaps relevant to those
triples telescope over an interval contained in[X/2,4X], by Bertrand
for the first predecessor and final successor; their sum is at most4X.
At most4X/(32log(X)) gaps exceed32log(X). Each bad gap excludes at most
two middle primes, so at most X/(4log(X)) middles are excluded. A positive
constant times X/log(X) triples therefore have both gaps <=32log(X).

Let epsilon_X=1000log(X)/X^(1/3), eventually less than1/4. Uniformly for
these short-gap triples, with h>=X/2 and q,p>=X,

    |N_q/q-f(q)|<=1/q<epsilon_X, f(x)=c*x^(2/3),
    0<=N_q/h-N_q/q<=epsilon_X,
    0<=N_q/q-N_q/p<=epsilon_X.

For example N_q<=c*(2X)^(5/3)+1, so the larger displacement is at most
2*(c*(2X)^(5/3)+1)*32log(X)/X^2; the constant1000 is ample eventually.
Only o(X/log(X)) INTEGER q in[X,2X] have dist(f(q),Z)<=2epsilon_X.
Indeed f' is bounded below by a positive constant times X^(-1/3), so
each inverse integer window has length O(epsilon_X*X^(1/3))=O(log(X)).
There are O(X^(2/3)) such windows meeting f([X,2X]) and its endpoints.
Counting at most length plus one integer points per interval gives
O(X^(2/3)log(X))=o(X/log(X)). This removes bad primes by counting all
bad integers, without a prime equidistribution or short-interval theorem.

Thus short-gap middle primes exist outside those windows in every large
dyadic range. For them all three ratios N_q/h,N_q/q,N_q/p lie between
the same consecutive integers: f(q) is more than2epsilon_X from endpoints,
and each total displacement is strictly less than2epsilon_X. Hence

    ceil(N_q/h)=ceil(N_q/q)=ceil(N_q/p),
    r_h=r_q=r_p=R.

These yield infinitely many distinct actual integer lengths N_q.

For any common adjacent order t and large such q, the capacity formulas
are unclipped and

    T_t(p)-T_t(q)=(t-1)(p-q), U(p)-U(q)=R(p-q).

Since R>=(t-1)^2, R/(t+1)>=t-3. The floor difference is at least the
floor of the difference, so both entries of the maximum in d_t fall by
at least(t-3)(p-q) when p is replaced by q. Including prefix costs gives

    ell_old(p,t)-ell_old(q,t)>=(8t-20)(p-q)>0

eventually. In particular the two-prime reference selects degree q.
Let t_dual be its selected adjacent order there. It is also an available
global order at degree h: h/q->1 on these triples, so the global
admissibility filters hold for h,t_dual at all sufficiently large N_q.
The same slope argument now compares q and h, while the current finite
theorem contributes its extra2 gates. Therefore

    L_dual(N_q)-L_Nq >=(8t_dual-20)(q-h)+2.                 (5)

Since q-h>=2 for large odd primes, and t_dual=Theta(N_q^(1/5)), the
right side is Omega(N_q^(1/5)) and positive on this infinite set.
No equality with actual circuit sizes, strict gain at every length, or
uniform quantitative prime-gap bound is claimed. The principal leading
lower bound remains valid at EVERY sufficiently large integer length.

## 8. Retained certificate comparisons

The positive-deficit threshold comparison can still be made at the
CURRENT globally selected parameters. Put A=kr/t,J=kt,R=r(r+1), and

    d_coarse=floor(A)+J+R+1,
    d_base=H+max(J,floor(A))+1.

Both are sufficient dimensions for this same current family: T<=kt,
U/(t+1)<kr/t and R>=H imply d_N<=d_base<=d_coarse. Section 5 supplies
t/sqrt(r)->1, so A/J->1, r=o(J), and

    2d_base-d_coarse=|A-J|+2r+O(1)=o(J).

Also d_base-d_N=O(k+t^2+r^2/t)=O(N^(3/5)) from the two maximum entries,
so it is o(J). Prefix4k and the degree correction are lower order. Thus

    L_N/(4k+8d_coarse)->1/2.

This is a same-family comparison of sufficient remainders, not actual
circuit complexities. It does not improve other constructions with
polylogarithmic remainders.

For the supplied earlier two-stage reference, define F_pre on N=0 to
output0, with d_pre=0 and L_pre=0 there. For N>=1 let k0 be the least
prime at least the least positive b0 with b0^5>32N^3, and k_pre the least
prime >=ceil(12k0/25). If k_pre>=N, output0 and set d_pre=0,L_pre=4N.
Otherwise put m_pre=N-k_pre,r_pre=ceil(m_pre/k_pre),t_pre=ceil(sqrt(r_pre)).
If t_pre>=r_pre, again output0 and set d_pre=0,L_pre=4N. On all remaining
lengths use the same certificate-prefix evaluator from Section 1 with
k_pre,r_pre,t_pre, returning0 on an invalid prefix, and define

    H_pre=r_pre(r_pre+3)/2,
    T_pre=sum_(j=2)^t_pre max(k_pre-j,0),
    U_pre=sum_(j=1)^r_pre max(k_pre-j,0),
    d_pre=H_pre+max(T_pre,floor(U_pre/(t_pre+1)))+1,
    L_pre=4k_pre+8d_pre-4.

This specifies its degree, dimension, guards, exceptional lengths, prefix
interpretation, irreducibility test and output recurrence completely. It is
one fixed total P reference with empty NP witnesses by Section 2's argument.
All sufficiently large lengths pass these guards and obey k_pre>r_pre+1,
d_pre<m_pre. The supplied general algebraic proof and the earlier finite
ledger in elimination.md Section3 independently give this complete certified
loss L_pre there; no prior acceptance is used.
The standard PNT
consequence gives k0/N^(3/5)->2 and k_pre/N^(3/5)->24/25, so

    L_pre/N^(4/5)->8F(24/25),
    L_N/L_pre->((5/4)*2^(1/5))/(sqrt(24/25)+625/1152)<19/20.

For an exact rational certificate, 2^(1/5)<23/20 because23^5>2*20^5,
and sqrt(24/25)>979/1000 because25*979^2<24*1000^2. Hence
(5/4)*2^(1/5)<23/16<(19/20)*(979/1000+625/1152)<(19/20)*F(24/25).
This compares complete certified losses of different fixed languages,
not their actual minimum circuit sizes or a prior accepted theorem.

## 9. Provenance, computation and limits

The certificate-prefix evaluator, exact prime-degree test and original
algebraic proof are retained from db5f9c6e-607f-424a-89f8-5a49d29d07df.
The compact full-B2 promised-pair method and four-object closure are
retained from28feddfa-dfc6-412a-a507-a525b15567fe. Positive-deficit
accounting and the half-loss comparison build on440a0962-3f37-4edb-8097-86c3b209b5c1;
the two-stage reference ondce46c96-6868-4a18-875f-290e1e06b980; full-hypothesis
averaging and exact dimension onc42988f5-ee13-407e-a8b0-4932794457fd;
the inert ledger and cleanup-C clarification on3b8a3eb4-444e-4250-a94f-b8e5578f7772;
the direct-degree/PNT reference and leading optimization on811df8a3-8fde-441f-aee7-1b18dfcadce0;
the degree-sensitive gate theorem on6fd4cdbe-2d59-4478-9d8c-81812124629f;
adjacent-order selection on5629bae8-dba0-4c3b-8569-722a079fe560;
and the two-prime reference and bad-integer-window strategy on
6a80eec7-0e26-45f3-ab57-06d231c2ce17. All relevant arguments are supplied
here. Publication status or favorable judgments are never proof premises.

The current new work is elimination.md Section6 and global selection,
its complete comparison and global leading-certificate optimality above.
This work does not claim a leading coefficient greater than4. Familiar
standard field algebra, Gaussian elimination, elementary calculus,
Bertrand's theorem and the explicitly stated PNT are used as specified;
all unfamiliar construction and model-transfer arguments are reproduced.

The declared finite-audit({}) is supplementary exact bounded computation.
It checks the predicates described in audit.py, including parameter
selection, arithmetic inequalities, legal promised-pair lines and executable
small families. It does not prove uniformity, asymptotics, all circuits,
prime distribution or the universal lemmas. Those rest on the written
arguments. Organizer checker fixtures and simulated runs are not candidate
proof checks. Independent advisory critique is separate from configured
provider review and human promotion. No prior review or score is assumed.
