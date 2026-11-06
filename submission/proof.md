# A compact coefficient-four family

We prove one fixed total uniform P language F has C_B2(F_N)>=4N-O(N^(4/5))
at EVERY sufficiently large integer length N. Thus liminf C_B2(F_N)/N>=4.
This is a literature-based improvement of the checkout's 31/10 seed, not a
new leading coefficient beyond the existing coefficient-four submissions.

## 1. Exact definition at every length

For N>=0 let b be the least positive integer with b^5>32N^3, and let k
be the least prime >=b. If k>=N, output zero. Otherwise put m=N-k,
r=ceil(m/k), t=ceil(sqrt(r)). If 2<=t<r and k>r+1 do not both hold,
output zero. Set

    H=r(r+3)/2,
    T=(t-1)k-t(t+1)/2+1,
    U=rk-r(r+1)/2,
    d=H+max(T,floor(U/(t+1)))+1.

If d>=m, output zero. These guards define every exceptional length,
including N=0. All integer comparisons and square roots are exact.

On remaining lengths, read the first k actual input bits b_j as the
lower coefficients of Q(Z)=Z^k+sum_(j=0)^(k-1)b_j Z^j over F_2.
Require Z^(2^k)=Z mod Q and gcd(Q,Z^2+Z)=1; otherwise output zero.
For an accepted prefix use the field K=F_2[Z]/Q. Pad the m payload bits
with rk-m zeros and interpret the r blocks a_i in the polynomial basis,
first bit the constant coordinate. Put y_i=a_i^(2^(r+2-i)-1).
Output the constant coordinate of the elementary symmetric polynomial
of order t in y_1,...,y_r. To compute it without division, initialize
E_0=1, E_1=...=E_t=0 and update E_j <- E_j+y_i E_(j-1) in descending
j=min(i,t),...,1 for each successive y_i. The subsets containing i give
exactly the new term, proving this recurrence. Code is in family.py.
The prefix is part of the actual input, not advice or an NP witness.

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

The threshold scan, prime scan and trial primality have polynomially many
steps in N: Bertrand gives k=O((N+1)^(3/5)). Parameters have O(log(N+1)) bits, and every operation
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


## 3. Full-B2 finite bound and bit accounting

The supplied disperser.md Sections1--7 prove that for prime k and
1<=t<r, the payload polynomial above is nonconstant on every affine
subspace of K^r whose binary dimension exceeds H+max(T,U/(t+1)).
The guard k>r+1 makes the displayed unclipped T,U agree with the sums
in that lemma. Our integer d strictly exceeds this real threshold:
when T is the maximum, the added1 suffices; otherwise floor(U/(t+1))+1
strictly exceeds U/(t+1).

An irreducible degree-k prefix exists by Section2. Fix any such prefix.
Padding is an injective linear map from the m actual payload bits into
K^r, hence maps every dimension-d affine payload space to an affine
space of that same dimension. Its payload G_Q is therefore a
dimension-d affine disperser. The whole prefix language is NOT asserted
to be an affine disperser, because rejected prefixes give zero regions.

The supplied elimination.md Sections1--2 prove the finite full-B2 theorem
C_B2(g)>=4m-8d for any dimension-d affine disperser g on m bits.
It explicitly handles all16 binary truth tables, constant/unary and
repeated-pin simplification, final literals, arbitrary sharing/fanout,
and affine-system rows charged one for one with gates. Its legal
parameter restrictions preserve an injective affine map and correctness
on the promised pair domain; no direct-sum premise is used.

Fixing the k-bit prefix in any circuit for F_N and normalizing constants
cannot increase its gate count, including under the pinned full-B2
conventions. Thus for every length passing the guards,

    C_B2(F_N)>=C_B2(G_Q)>=4(N-k)-8d=4N-(4k+8d).

The paid prefix loss4k is included. Witnesses are not circuit inputs.
No claim equates potential, circuit size, or the finite certificate.

## 4. Every sufficiently large integer length

The threshold satisfies b=2N^(3/5)+O(1). Bertrand's theorem says every
integer b>=2 has a prime between b and2b. Consequently the least prime
k>=b obeys k=Theta(N^(3/5)). This holds at EVERY integer length,
without a prime number theorem or sparse-length assumption. Hence
m=N-k=Theta(N), r=ceil(m/k)=Theta(N^(2/5)) and
t=ceil(sqrt(r))=Theta(N^(1/5)). Eventually 2<=t<r, k>r+1 and k<N.
Further H=O(r^2)=O(N^(4/5)), T<=kt=O(N^(4/5)) and
U/(t+1)<=rk/(t+1)=O(N^(4/5)). Therefore d=O(N^(4/5))<m eventually,
so every guard passes at every sufficiently large integer N.
The complete loss4k+8d=O(N^(4/5))=o(N). In particular for each real
epsilon>0 some N_epsilon makes (4k+8d)/N<=epsilon for EVERY integer
N>=N_epsilon. The finite theorem gives C_B2(F_N)>=(4-epsilon)N there.
This proves the canonical coefficient4/exponent1 joint declaration.

## 5. Provenance and finite evidence

This submission reuses the mathematical arguments and evaluator from
submission8eef1799-c27f-46f3-8e61-3c4e0d28d145, frozen candidate commit
b4a12ebc5267f63705e1c277565e2779c088c286. Its original provenance includes
Ben-Sasson--Kopparty's affine disperser and Mishra--Santhanam ECCC TR26-260
revision1 (November3,2026), with the supplied finite4n-8d refinement.
The relevant arguments are reproduced in full; no prior judgment,
incidental source claim or URL is used as a proof premise. This package
removes global certificate optimization, degree-sensitive refinements,
PNT comparisons and historical auxiliary claims. It claims no improvement
on an already proved coefficient4. It is a compact route from the seed.

Familiar standard premises are elementary finite-field theory, univariate
root counting, Gaussian elimination and Bertrand's prime theorem. The
unfamiliar coefficient isolation and full circuit elimination arguments
are supplied in disperser.md and elimination.md. Their correctness,
not any favorable checker or review output, is needed for the theorem.
The finite-audit experiment checks bounded arithmetic, field-test and
small evaluator predicates. It does not prove universal circuit coverage,
asymptotics or eligibility. These rely on the supplied written arguments.
