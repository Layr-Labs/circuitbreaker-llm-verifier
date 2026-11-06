# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). Coefficient 4 and the evaluator use existing
mathematics; the contribution here is the sharper supplied block-selection
and dimension estimate, and its resulting certified remainder. The gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, revision 1, November 3, 2026,
Sections 2-4. The algebraic construction is based on Ben-Sasson
and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*, author
version March 12, 2010, Sections 3, 6 and 7. The necessary arguments are
supplied in our notation in `elimination.md` and `disperser.md`.

Our input-prefix wrapper removes the deterministic irreducible-polynomial
construction prerequisite of the repository's BK/Shoup seed. It explicitly
pays for all prefix bits. The **whole** function is not an affine disperser:
only a fixed-prefix restriction is asserted to be one. No conjecture,
restricted circuit basis, direct-sum assumption or exponent gain is used.
We assert no exponent above 1, endpoint, or optimality over all families.

## 1. One evaluator on every bit string

The first bit of each field block is its constant coefficient. On x of length N:

1. If N=0, output 0.
2. Let B be the smallest positive integer with B^5 > 32*N^3. Let k0 be the
   smallest prime at least B. Let k be the smallest prime at least
   ceil(12*k0/25). Use deterministic trial division for both scans.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k), t=ceil(sqrt(r)).
   If t>=r, output 0.
4. Read the first k bits b_0,...,b_(k-1) as

       Q(Z)=Z^k + sum_(j=0)^(k-1) b_j Z^j    in F_2[Z].

   Check Z^(2^k)=Z modulo Q and gcd(Q,Z^2+Z)=1. If either fails, output 0.
5. Use K=F_2[Z]/(Q). Append rk-m zero bits to the remaining m payload bits
   and interpret the r blocks of k bits as a_1,...,a_r in K. This is an
   injective F_2-linear map of the payload.
6. For i=1,...,r set e_i=r+1-i and y_i=a_i^(2^(e_i+1)-1). Output the
   constant coefficient of the elementary symmetric polynomial of order t
   in y_1,...,y_r, in its representative of degree less than k.

For explicit evaluation, initialize E_0=1 and E_1=...=E_t=0. For each i,
update E_j <- E_j+y_i*E_(j-1), in decreasing order j=min(i,t),...,1.
Splitting subsets by whether they contain i proves this division-free
recurrence. Binary powering takes O(e_i+1) multiplications per y_i.

This specifies every choice and every input length. Q is part of x, not
advice or a value supplied by an unspecified construction algorithm.

## 2. Exact prefix test and existence of hard prefixes

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

## 3. Uniform P and NP with empty witnesses

Bertrand's prime theorem gives k0<2B and k<2*ceil(12*k0/25), so

    2*N^(3/5) < k0 < 4*N^(3/5)+2,
    12*k0/25 <= k < 24*k0/25+2.

In particular k=Theta(N^(3/5)) and k is prime.

The prime scan and trial division are polynomial in N. B is computed by
integer binary search and fifth-power comparisons. All subsequent integer,
polynomial and field operations have polynomial bit cost. Field multiplication
can use schoolbook O(k^2) operations; the powers take O(r^2) multiplications
in total and the symmetric recurrence O(rt). A coarse O((N+1)^6) bound
suffices for the whole algorithm, including rejected inputs.

Let A be this fixed evaluator, L={x:A(x)=1}, W(N)=0, and let V(x,w)
accept exactly when w is empty and A(x)=1. If x is in L the empty witness
is accepted; if a permitted witness is accepted then A(x)=1 and x is in L.
This proves both directions with one uniform polynomial-time verifier.
Witnesses are not inputs to the target circuits.

## 4. All sufficiently large integer lengths pass the guards

We have k=Theta(N^(3/5)), m=N-k, r=Theta(N^(2/5)),
t=Theta(N^(1/5)). Consequently k<N and t<r for every sufficiently
large integer N. The former guard kt>=r(r+1) is removed: the sharpened
deficit-selection lemma does not need it. These estimates cover all large
integer lengths, not a sparse subsequence.

## 5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    H=r(r+3)/2, T=sum_(j=2)^t max(k-j,0), U=sum_(j=1)^r max(k-j,0),
    d_N = H+max(T,floor(U/(t+1)))+1.

`disperser.md` proves that, for prime field degree k and r>t>=1, the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
H+max(T,U/(t+1)), for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

Since kr=O(N), r=Theta(N^(2/5)) and t=Theta(N^(1/5)),

    d_N=O(N^(4/5)),          k=O(N^(3/5)).

We use this conservative sublinear disperser dimension, not later
logarithmic-dimension constructions. Either suffices for the leading pair.

Fixing prefix inputs of any circuit for F_N yields a circuit for G_Q with
no additional gates. `elimination.md` proves the full-B2 bound

    |C| >= 4m-8d

for a dimension-d affine disperser. It includes arbitrary sharing,
repeated pins, constants, degenerate gates and terminal-output cases.
Apply it to the restricted circuit (or its no-larger normalized form):

    C_B2(F_N) >= 4(N-k)-8d_N
               = 4N-(4k+8d_N).

Set R(N)=4k+8d_N eventually and extend it on the finitely many guarded
lengths. Then R(N)=O(N^(4/5))=o(N). The choice of Q in this lower-bound
restriction is not advice to the evaluator: every circuit must work for
all prefixes, while the uniform algorithm simply reads and tests its input.

For each epsilon>0, take N_epsilon large enough that the guards pass and
R(N)/N<=epsilon for every N>=N_epsilon. This gives exactly
C_B2(F_N)>=(4-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete gate argument. `audit.py` supplies
supplementary finite checks of Boolean identities, rational estimates and
the prefix test; these do not replace either universal proof. No claimed
review, checker success or prior acceptance is a mathematical premise.

Sources, for attribution and version identification:

* Mishra--Santhanam, ECCC TR26-260 revision 1, 2026-11-03, Sections 2-4,
  https://eccc.weizmann.ac.il/report/2026/260/download/
  The revision is identified by its date and title; the older version hash does not identify this revision.
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4,
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the protected full-B2 conventions. Their LP and lower-bound case
  analysis are not invoked by this submission.

## 7. Reused submission material

The evaluator, certificate-prefix family, algebraic disperser proof and bounded
audit are adapted from submission db5f9c6e-607f-424a-89f8-5a49d29d07df
(candidate commit 02e21edbaa931b85f664acd2da6873155dc57716).
The full gate argument in elimination.md, including its 4n-8d finite refinement,
is reproduced from submission 28feddfa-dfc6-412a-a507-a525b15567fe
(candidate commit a28162b3c8200880a6cabcdb07572205030bc491).
All premises are supplied here; neither source's review outcome is a premise.
The current evaluator retains the certificate-prefix design but uses the
smaller degree rule and removes the obsolete guard. The proved coefficient
is 4 rather than the original 7/2. This package does not claim
a larger coefficient than the other coefficient-4 candidates or a smaller
remainder than their polylogarithmic remainders. The immediate own predecessor 440a0962-3f37-4edb-8097-86c3b209b5c1
supplies the refined positive-deficit lemma and remainder comparison. Its
publication or review outcome is not a mathematical premise. Its complete
threshold argument is reproduced here. The smaller paid field degree and its uniform loss comparison are carried
from own predecessor dce46c96-6868-4a18-875f-290e1e06b980. The current
new contribution is use of the full independent-source hypothesis, giving
the sharper exact threshold and gate improvement in Section 10.

## 8. Sharper certified remainder for the same family

The deficit-selection proof in disperser.md Section 5 implies the former
refined dimension d_base=kt+H+1 is sufficient; Section 10 proves our new
sharper d_N is no larger and differs from d_base by O(k+r). Here H=r(r+3)/2.
For comparison with the original coarse bound, put d_new=d_base.
For this section, compare thresholds at the current k,r,t; Section 9
separately compares degrees across two evaluators. Define

    A=kr/t, B=kt, R=r(r+1), H=r(r+3)/2,
    d_old=floor(A)+B+R+1, d_new=B+H+1.

The former sufficient threshold is also proved by the supplied argument:
it is no smaller than the new threshold eventually, as shown below. No
prior acceptance is assumed for either assertion. Both are O(N^(4/5)).
As r tends to infinity, A/B=r/t^2 tends to 1. Moreover H=R/2+r and
r/B=r/(kt)<=sqrt(r)/k tends to zero, using r=O(N^(2/5)) and
k=Theta(N^(3/5)). Floors and constants are o(B). Therefore

    2*d_new-d_old = B-A+2r+O(1) = o(B),
    d_old>=B, and hence d_new/d_old -> 1/2.

This proof does not require convergence of the prime ratio k/N^(3/5).
In particular d_new<d_old for every sufficiently large integer length.
The base certified circuit losses are L_new=4k+8d_new and
L_old=4k+8d_old. Since k/B=1/t tends to zero,

    L_new/L_old -> 1/2 < 1.

This compares two proved lower-bound remainders for the SAME fixed family,
not actual minimum circuit complexities. It improves the earlier compact
algebraic package's certified remainder while retaining coefficient 4 and
exponent 1. It does not claim to improve the stronger polylogarithmic
remainders supplied for other families.

## 9. Smaller field degree with the same sharpened lemma

The reference family uses k0 as its field degree and the refined threshold
d0=k0*t0+r0(r0+3)/2+1, where m0=N-k0, r0=ceil(m0/k0),
t0=ceil(sqrt(r0)). Its complete certified loss is L0=4k0+8d0.
The prior submitted proof and unchanged algebra supplied here establish
this reference independently of any platform outcome. The current evaluator
uses the new prime k specified in Section 1. Put lambda=k/k0. Bertrand
gives 12/25<=lambda<=24/25+o(1). Removing the obsolete guard is
justified by the actual lemma, which needs only prime k and 1<=t<r.

Uniformly for these choices, elementary rounding and k=Theta(N^(3/5)) give

    d_base(k)=sqrt(k*N)+N^2/(2*k^2)+o(N^(4/5)),
    d_N(k)=d_base(k)+O(k+r),
    L(k)=4k+8*d_N(k)=8*d_base(k)+o(N^(4/5)).

Section 10 shows d_N differs from d_base by O(k+r), which is subleading.
Indeed r=N/k+O(1), t=sqrt(r)+O(1), so kt=sqrt(kN)+O(k),
and r(r+3)/2=N^2/(2k^2)+O(N/k). Both error terms are o(N^(4/5)).
Write a=k0/N^(3/5), so 2<a<4+o(1), and s=1/(2*a^(5/2)).
Then 0<s<=1/(8*sqrt(2))<177/2000 eventually, and the leading ratio is

    (sqrt(lambda)+s/lambda^2)/(1+s).

For fixed a the unnormalized expression as a function of lambda has a
single minimum, because its derivative has the sign of
lambda^(5/2)-4s. Hence its maximum on [12/25,24/25] is at an endpoint.
At either endpoint the ratio increases with s: lambda^(-2)>sqrt(lambda).
The lower endpoint has sqrt(12/25)<693/1000; the upper endpoint has
sqrt(24/25)<49/50. Exact rational bounds at s=177/2000 give

    ((693/1000)+(177/2000)*(25/12)^2)/(2177/2000) < 99/100,
    ((49/50)+(177/2000)*(25/24)^2)/(2177/2000) < 99/100.

The interval's o(1) endpoint perturbation and all rounding errors vanish
uniformly since a stays bounded above and away from zero. Thus

    limsup L(k)/L0 <= 99/100 < 1.

This is a comparison of complete certified losses, including field-prefix
bits. It need not compare actual circuit complexities, fields, or functions
between the current family and the reference family. Both are fixed total
uniform languages and both retain the native coefficient-four theorem.

## 10. Use the actual independent-source threshold

The new selection proof uses sum delta_i+max d_i>k+1, precisely the
hypothesis of the coefficient-isolation lemma. The preceding coarser
selection demanded sum delta_i>k. This avoids paying an unnecessary
maximum-block loss. At the SAME current field degree and the same
fixed family, let d_base=kt+H+1 and let d_N be the sharper dimension
of Section 5. The exact clipped-capacity sum T is polynomial-time
computable: with u=min(t,k-1),

    T=(u-1)k-u(u+1)/2+1.

For prime k>=2 and t>=1, u>=1, so this formula includes empty sums.
Since T<=k(t-1), and r<=t^2,

    U/(t+1) <= r(k-1)/(t+1) <= k*t^2/(t+1) = kt-kt/(t+1),
    T <= kt-k <= kt-kt/(t+1).

It follows, keeping the floor exactly, that

    d_base-d_N >= ceil(kt/(t+1)).

Consequently the actual circuit lower bound improves on the preceding
same-family certified bound by at least 8*ceil(kt/(t+1)) gates:

    C_B2(F_N) >= 4N-(4k+8d_N)
              >= 4N-(4k+8d_base)+8*ceil(kt/(t+1)).

This is an exact sufficient-remainder comparison for every sufficiently
large integer N, not a statement about actual minimum circuit-size ratios.
Also T is at least the unclipped sum kt-k-t(t+1)/2+1, so

    0 <= d_base-d_N <= k+t(t+1)/2 = O(k+r).

Thus d_N=O(N^(4/5)); the threshold-halving limit in Section 8 and the
99/100 cross-degree limit in Section 9 remain true for the sharper loss.
The evaluator is unchanged from the smaller-degree predecessor; only the
proved dimension and certified remainder improve.

This refinement builds on own submission dce46c96-6868-4a18-875f-290e1e06b980
for the smaller-degree evaluator and its resource comparison. The full
arguments are supplied here independently of any platform outcome.

### 10.1 Exact eventual saving

For this family k/r tends to infinity and k/r^2 tends to zero. Eventually
k>r+1 and k>t, so U=rk-r(r+1)/2 and
T=(t-1)k-t(t+1)/2+1. Direct subtraction gives

    (t+1)T-U
      = (t^2-1-r)k+r(r+1)/2-t(t+1)^2/2+(t+1).

Since r<=t^2, its first term is at least -k. Dividing this lower bound
by r^2 gives limit 1/2: k/r^2->0 and t=Theta(sqrt(r)). Hence it is
positive eventually at every integer length. Therefore the maximum
in d_N is eventually attained by T, and

    d_N=H+(t-1)k-t(t+1)/2+2,
    d_base-d_N=k+t(t+1)/2-1.

Thus the SAME family enjoys the exact eventual certified improvement
8k+4t(t+1)-8 gates over its supplied preceding base-dimension bound.
This remains a remainder improvement, with native projection (4,1).
