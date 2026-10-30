# Published coefficient 31/10: construction and proof

This applies the cited LY-full, BK-preprint and Shoup-preprint results; it claims no new theorem. Source versions, locators and hashes are listed below.

## Fixed total family

Set f_0(empty)=0. For n>=1, let k be the least prime strictly above 2*n^(3/5),
r=ceil(n/k), t=ceil(sqrt(r)). Fix once a deterministic procedure from Section 3
realizing Shoup-preprint Theorem 3.2, with deterministic choices fixed by that procedure,
to produce a monic irreducible q_k over F_2 of degree k. Let K=F_2[Z]/(q_k), with
polynomial basis 1,Z,...,Z^(k-1). This fixes a single language, not per-length advice.

Append r*k-n zeros to x and split it into r blocks. Interpret each block in K with
its first bit the constant coefficient; this is an injective F_2-linear map phi.
Let pi take the constant coefficient, so pi is linear and pi(1)=1. For 1<=i<=r put
y_i=phi_i(x)^(2^(i+1)-1). Define f_n(x)=pi(Sym_r^t(y_1,...,y_r)), where Sym_r^t
sums the products indexed by t-element subsets. This uses the precise ascending
powers in BK-preprint Theorem 2.2, equation (4).

## Uniform P and NP membership

One deterministic algorithm computes this function at all lengths. Find k by scanning
and primality trial division, using q^5>32*n^3 for the strict threshold. Bertrand's
bound makes k=O(n^(3/5)), so this scan is polynomial in actual input length n.
Shoup's unconditional fixed-prime procedure is polynomial in k. Binary powering uses
only O(i+1) field multiplications for each exponent, whose bit length is i+1.
Field operations in the polynomial basis have polynomial bit cost.

Compute the symmetric polynomial without characteristic-two division: start E_0=1
and E_j=0 for 1<=j<=t; after reading y_i update E_j=E_j+y_i*E_(j-1) in descending
j from min(i,t) to 1. Splitting subsets by inclusion of i proves the recurrence;
it takes O(r*t) field operations. Output pi(E_t). Thus one polynomial-time evaluator
A computes the total family. Use V(x,w)=1 iff w is empty and A(x)=1, with witness
bound zero. Both directions are immediate, establishing L in P, hence NP. The
witness is not an input to the target circuits.

## Every integer length and the lower bound

BK-preprint Theorem 2.2 applies to arbitrary integer input length via the injective
map, not a divisibility assumption or a sparse subsequence. With p=2 and the above
k,r,t,phi,pi and powers, it is nonconstant on all affine subspaces of dimension
strictly greater than 6*n^(4/5). Set d(n)=floor(6*n^(4/5))+1. Eventually 1<=d(n)<=n
and d(n)/n<=6*n^(-1/5)+1/n tends to zero. Zero padding is inside the theorem's
injection; the circuit bound remains measured against n input bits.

Apply LY-full Theorem 1.1 to obtain C_B2(f_n)>=(31/10)*n-o(n). The P evaluator
establishes eligibility and places no uniformity restriction on circuits. The full
B2 model includes all 16 binary truth tables, unbounded sharing and depth, one
output and all input assignments; inputs/wires are not gates. The source allows
constant nodes during substitutions and handles repeated pins and unary effects by
normalization. Output cases must be covered; no free extra computational basis is
assumed. Sources: LY-full Sections 1.1, 2.1, 2.4 footnotes 3–4 and Section 3.3.

## Conservative remainder, with source discrepancy retained

LY-full Corollary 3.8 bounds potential by n/2+5*d/2, while its Theorem 3.9 displays
n/2+d in a finite expression. The finite expressions differ. The following direct derivation uses the larger
loss from Corollary 3.8, together with Definition 3.6 and Theorem 4.1:

    delta*(n-2*d-2) <= mu <= |C|+alpha_I*n+alpha_phi*(n/2+5*d/2).

At alpha_I=48/5, alpha_Q=9/5, alpha_phi=1/5 and delta=64/5 this gives

    |C| >= (31/10)*n-(261/10)*d-128/5.

The loss is o(n) for the specified d. This does not re-prove the paper's long case
analysis. The seed coefficient is exactly 31/10; no stronger remainder is asserted.

## Exact arithmetic and the unchanged LP

The six terms bounding delta-alpha_I in LY-full Theorem 4.1 are

    alpha_I/3, 2-2*phi+Q, 4-4*phi, 3+phi, 5-Q, (5-2*phi+Q)/2.

At the published point they are 16/5,17/5,16/5,16/5,16/5,16/5. The admissibility
conditions include 0<=Q<=alpha_I and 0<=phi<=1/2. The coefficient is
c=delta-alpha_I-phi/2=31/10. The declared optional checker inputs verify this arithmetic using exact fractions.

Writing u=delta-alpha_I, combine u<=3+phi, u<=5-Q and
u<=(5-2*phi+Q)/2 with weights 7/10,1/10,1/5. They sum to
u<=31/10+phi/2, hence c<=31/10 for this unchanged LP. The feasible point attains
it. Retuning this LP alone cannot advance the record; the certificate is not a
barrier to another proof. Passing the LP checker certifies only these arithmetic
claims, not the full circuit theorem, NP construction or exhaustive case analysis.

## Appendix: explicit remainder and limit check (editorial, no new claim)

Let d(n)=floor(6*n^(4/5))+1. The proved bound is |C| >= (31/10)n - R(n) with
R(n) = (261/10)*d(n) + 128/5. Since d(n) <= 6*n^(4/5)+1, we get
R(n) <= (1566/10)*n^(4/5) + 261/10 + 128/5 = 156.6*n^(4/5) + 51.5, so
R(n)/n <= 156.6*n^(-1/5) + 51.5/n -> 0. Hence for every epsilon>0 there is
N_epsilon with |C| >= (31/10 - epsilon)*n for ALL integers n >= N_epsilon,
which is exactly liminf C_B2(f_n)/n >= 31/10 over every sufficiently large length,
not a sparse subsequence. The bound is monotone-free: it uses only the pointwise
inequality at each n, so no averaging over lengths is needed. The exponent p=1 is
attained trivially and no logarithmic factor is claimed (q=r=0). This appendix
restates the existing derivation; the theorem and the family are unchanged.

## Sources

- LY-full: ECCC TR21-023, front-page date 2021-02-20; Theorem 1.1 p.2; Sections 1.1,2.1,2.4; Definition 3.6, Corollary 3.8 pp.14–15; Theorem 4.1 pp.18–19; https://eccc.weizmann.ac.il/report/2021/023/download/; SHA256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062
- BK-preprint: Author preprint 2010-03-12; Theorem 2.2, equation (4), printed pp.6–7 (PDF pp.7–8); https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf; SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf
- Shoup-preprint: Author preprint 1989-01-31; Section 3, Theorem 3.2, printed p.3 and fixed-prime consequence; https://www.shoup.net/papers/detirred.pdf; SHA256 00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644
