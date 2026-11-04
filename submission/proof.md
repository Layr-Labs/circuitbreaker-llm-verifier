# Coefficient 31/10 at exponent 1 for one fixed BK/Shoup family: construction and proof

This submission applies three published results (LY-full, BK-preprint, Shoup-preprint;
exact versions, page locators and SHA-256 digests in "Sources") and claims no new
circuit theorem. Compared with the seed it (i) shows that the disperser property
does not depend on which irreducible polynomial the fixed algorithm outputs,
(ii) derives the explicit size bound term by term from LY-full Definition 3.6 rather
than from the displayed finite expression of LY-full Theorem 3.9, and (iii) gives
exact thresholds for every asymptotic step. Notation: n is always the actual input
length |x|. BK-preprint writes m for this length and n for the prime field degree;
here that prime is called k.

## 1. The fixed total family

Set f_0(empty)=0. For n>=1 define:

- k = k(n): the least prime with k > 2*n^(3/5), equivalently k^5 > 32*n^3
  (both sides are positive, so the equivalence is exact integer arithmetic);
- r = ceil(n/k) and t = ceil(sqrt(r)), so 1 <= t <= r and r*k >= n;
- q_k: the monic irreducible polynomial of degree k over F_2 output by one fixed
  deterministic algorithm S realizing Shoup-preprint Theorem 3.2 (Section 3)
  with p=2, all internal choices fixed by S;
- K = F_2[Z]/(q_k), a field with 2^k elements, with basis 1,Z,...,Z^(k-1);
- phi: F_2^n -> K^r: append r*k-n zero bits to x, cut the result into r
  consecutive blocks of k bits, and read block i as the element of K whose
  coefficient of Z^j is bit j of the block (j=0 first);
- pi: K -> F_2, the coefficient of Z^0;
- y_i = phi_i(x)^(2^(i+1)-1) for 1<=i<=r;
- f_n(x) = pi(Sym_r^t(y_1,...,y_r)), where Sym_r^t is the sum over all
  t-element subsets I of {1..r} of the product of y_i over i in I.

L is the set of x with f_|x|(x)=1. It is one fixed language: it is exactly the
language accepted by the single deterministic algorithm A of Section 2, which
contains S. No per-length advice or choice is involved.

Properties used later. phi is F_2-linear (padding with fixed zeros and regrouping
bits are linear) and injective (the padded string, hence x, is recovered from the
block coefficients). pi is F_2-linear and nontrivial since pi(1)=1. For p=2,
1+p+...+p^i = 2^(i+1)-1, so these are exactly the powers of BK-preprint
Theorem 2.2, equation (4).

## 2. Uniform polynomial time, hence NP

Algorithm A on input x with n=|x|: if n=0 output 0. Otherwise:

1. Find k: test q = 2,3,4,... for primality by trial division and stop at the
   first prime with q^5 > 32*n^3. By Bertrand's postulate some prime lies in
   (2*n^(3/5), 4*n^(3/5)+1], so at most O(n^(3/5)) candidates are tried, each in
   polynomial time.
2. Run S to obtain q_k. Shoup-preprint Theorem 3.2 bounds the cost by
   O(p^(1/2)(log p)^3 k^(3+eps) + (log p)^2 k^(4+eps)) operations in F_p; for p=2
   this is poly(k) <= poly(n).
3. Compute phi(x) and each y_i by square-and-multiply: the exponent 2^(i+1)-1 has
   i+1 bits, so at most 2(i+1) multiplications in K, each O(k^2) bit operations
   followed by reduction mod q_k. With i <= r <= n this is polynomial.
4. Compute E_t = Sym_r^t(y_1..y_r) without division (characteristic 2 forbids
   the Newton-Girard divisions): put E_0=1, E_j=0 for 1<=j<=t; for i=1..r update
   E_j <- E_j + y_i*E_(j-1) for j from min(i,t) down to 1. Invariant: after step
   i, E_j = Sym_i^j(y_1..y_i). It holds at i=0, and the update is the identity
   Sym_i^j = Sym_(i-1)^j + y_i*Sym_(i-1)^(j-1) (split subsets by whether they
   contain i); descending j ensures E_(j-1) still holds its step-(i-1) value.
   This costs O(r*t) field operations.
5. Output pi(E_t).

So A runs in time polynomial in n and computes f_n at every length. NP membership:
take witness bound W(n)=0 and V(x,w)=1 iff w is empty and A(x)=1. If x is in L
then w=empty gives V=1; if V(x,w)=1 for some |w|<=0 then w is empty and A(x)=1, so
x is in L. Both directions hold, L is in P, hence in NP. Witnesses are not circuit
inputs; the circuits below read exactly the n bits of x.

## 3. The disperser property at every length

BK-preprint Definition 2.1 (printed p.6): f: F_p^m -> F_p is an affine disperser for
dimension d if |f(A)|>1 for every affine subspace A of dimension at least d.
BK-preprint Theorem 2.2 (printed pp.6-7): for every integer m, with the prime field
degree the least prime above 2*m^(3/5), r=ceil(m/degree), t=ceil(sqrt r), ANY
injective F_p-linear phi into (F_(p^degree))^r and ANY nontrivial F_p-linear pi,
function (4) satisfies |f(A)|>1 for every affine A with dim(A) > 6*m^(4/5).

Hypothesis check with p=2, m=n: the degree is our k; r and t match; phi and pi
satisfy the hypotheses by Section 1. The theorem is stated for the field F_(2^k),
not for a particular model of it. K is such a field for whatever irreducible q_k
the algorithm S returns, and phi, pi are linear injective / nontrivial with respect
to K. Hence the conclusion holds no matter which irreducible q_k is produced; the
fixed algorithm only makes L a single language. There is no divisibility assumption
(phi absorbs the padding) and no restriction to a subsequence of lengths.
Theorem 2.2 is printed for every integer m, but its proof (BK Sections 3, 6, 7) uses
asymptotic estimates. For example, Theorem 7.3 needs t<r, which fails when r=1.
So we use it only in the weaker form "for every m >= n0", with n0 an absolute
constant. Nothing below needs small lengths: every later statement is for
n >= N0 := max(n0, 2^20), which suffices for a liminf over all integers.

Let d(n) = floor(6*n^(4/5)) + 1 > 6*n^(4/5). Every affine subspace of F_2^n of
dimension d(n) has dimension > 6*n^(4/5), so f_n is non-constant on it. This is
exactly LY-full Definition 2.1 (p.7): f_n is an affine disperser for dimension d(n).
For n >= N0 >= 2^20 we have 6*n^(4/5) <= 6*n/16 < n-1, hence 1 <= d(n) <= n. Moreover
d(n)/n <= 6*n^(-1/5) + 1/n -> 0.

## 4. Circuit model

LY-full Section 1.1 works with single-output circuits whose gates have fan-in
two and compute arbitrary functions in B2 = all 16 maps F_2 x F_2 -> F_2, with
unrestricted fan-out and depth. Sections 2.1 and 2.4 with footnote 3 count only
gates, not constants or variables; footnote 4 describes rewiring after substitutions.
Section 3.3 gives the normalization rules that remove gates fed by constants,
degenerate or trivial gates, and handles repeated pins and unary effects. The proof
of Theorem 1.1 is stated for semicircuits, a superset of acyclic circuits. This is
the model fixed in rules.md: every acyclic B2 circuit with one output that is exact
on all 2^n inputs, counting every retained gate once, is covered. No extra basis
is assumed; the lower bound places no uniformity requirement on circuits.
Degenerate outputs are covered too. For n >= N0, f_n is non-constant on F_2^n itself
(dimension n > 6n^(4/5)) and on every hyperplane x_j = b. So f_n is neither a constant
nor a literal x_j or its negation, and every circuit computing it outputs from a gate.
LY Theorem 4.1 holds for every B2 circuit computing the disperser, with constants,
unary effects and repeated pins handled by its normalization rules (Section 3.3).

## 5. Explicit lower bound at each length

Fix n >= N0 and ANY B2 circuit C computing f_n, with g=|C| gates. Apply the
following with d = d(n).

(a) LY-full Definition 3.6 (p.14): for a packing P and rdq-source R,
    mu(C,P,R) = g + alpha_I*i + alpha_Q*q + alpha_phi*Phi(C,P),
where i is the number of influential inputs and q the number of quadratic
equations of R. Take R = empty (no quadratic variables, Definition 2.2, p.8), so
q = 0. Influential inputs are input variables, so 0 <= i <= n.

(b) LY-full Lemma 3.7 (p.14; FGHK16 Lemma 3) and Corollary 3.8 (p.15): for a
circuit computing an affine disperser for dimension d and any packing P,
Phi(C,P) <= n/2 + 5*d/2.

(c) LY-full Theorem 4.1 (p.18): for a circuit computing an affine disperser for
dimension d and any packing P, mu(C,P,empty) >= delta*(n-2*d-2) where
    delta = alpha_I + min{alpha_I/3, 2-2*alpha_phi+alpha_Q, 4-4*alpha_phi,
                         3+alpha_phi, 5-alpha_Q, (5-2*alpha_phi+alpha_Q)/2},
for the nonnegative weights with alpha_phi <= 1/2 used in the proof of
Theorem 1.1 (pp.18-19).

Use the published point alpha_I=48/5, alpha_Q=9/5, alpha_phi=1/5 (all >= 0,
alpha_phi <= 1/2). The six terms in the min are exactly
    16/5, 17/5, 16/5, 16/5, 16/5, 16/5
(48/5/3; 2-2/5+9/5; 4-4/5; 3+1/5; 5-9/5; (5-2/5+9/5)/2 = (32/5)/2), so the
minimum is 16/5 and delta = 48/5 + 16/5 = 64/5.

Combining, with any packing P (for example the empty packing):
    g = mu - alpha_I*i - alpha_Q*q - alpha_phi*Phi
      >= (64/5)(n-2d-2) - (48/5)n - 0 - (1/5)(n/2 + 5d/2)
      = (64/5 - 48/5 - 1/10) n - (128/5 + 1/2) d - 128/5
      = (31/10) n - (261/10) d - 128/5.
(Check: 128/10 - 96/10 - 1/10 = 31/10 and 256/10 + 5/10 = 261/10.)
Each step uses only an upper bound on a subtracted nonnegative term, so it holds
even if Theorem 4.1's right side is negative. Since C was arbitrary,

    C_B2(f_n) >= (31/10) n - (261/10) d(n) - 128/5   for every n >= N0.   (*)

Source discrepancy, retained explicitly. LY-full Theorem 3.9 (p.15) displays the
finite loss alpha_phi*(n/2 + d), smaller than the n/2 + 5d/2 of Corollary 3.8 that
its proof cites. (*) does not use Theorem 3.9's finite expression; it uses the
larger, Corollary-3.8 loss. Both versions are o(n) for sublinear d, so the
asymptotic statement of Theorem 3.9 and Theorem 1.1 is unaffected either way.
This submission relies on, and does not re-prove, the case analysis behind
Theorem 4.1 (Section 4.1) and Proposition 3.10.

## 6. The asymptotic statement for all large lengths

From (*) and d(n) <= 6*n^(4/5) + 1,
    C_B2(f_n) / ((31/10) n) >= 1 - (10/31)*[(261/10)(6 n^(-1/5) + 1/n) + (128/5)/n].
The bracket tends to 0 along all integers n. Hence for every eps>0 there is
N_eps (explicitly any N_eps >= N0 beyond which the bracket is at most
(31/10)eps) such that C_B2(f_n) >= (1-eps)(31/10) n for EVERY integer n >= N_eps.
That is, liminf over all n of C_B2(f_n)/g(n) >= 1 for g(n) = (31/10) n^1
(log exponent 0, loglog exponent 0), which is the declared polynomial_log guarantee.
Its finite projection (31/10, 1) follows: liminf C_B2(f_n)/n >= 31/10 > 0 with
exponent 1. The same fixed family, verifier and witness bound give both coordinates.

## 7. The unchanged LP cannot exceed 31/10 (unscored lemma)

The LP in the proof of LY Theorem 1.1 (p.19) requires
delta <= alpha_I + min{six terms}, together with alpha_phi, alpha_Q, alpha_I >= 0 and
alpha_phi <= 1/2. Write u = delta - alpha_I and c = u - alpha_phi/2 (the coefficient
Theorems 3.9/4.1 yield). Then u is at most each of the six terms. In particular u <= 3+alpha_phi, u <= 5-alpha_Q and
u <= (5-2*alpha_phi+alpha_Q)/2. Weighting by 7/10, 1/10, 1/5 (sum 1):
    u <= (7/10)(3+alpha_phi) + (1/10)(5-alpha_Q) + (1/5)(5-2*alpha_phi+alpha_Q)/2
       = 21/10 + 5/10 + 5/10 + alpha_phi*(7/10 - 2/10) + alpha_Q*(-1/10 + 1/10)
       = 31/10 + alpha_phi/2,
so c <= 31/10 at every feasible point, and Section 5's point attains it. The declared
`lp-certificate` check verifies this exact arithmetic only; it does not certify the
circuit theorem, the construction or the case analysis. It is no barrier to other
proofs.

## Sources

- LY-full: Jiatu Li and Tianqi Yang, "3.1n - o(n) circuit lower bounds for explicit
  functions", ECCC TR21-023, front-page date 2021-02-20. Used: Theorem 1.1 (p.2);
  Sections 1.1, 2.1, 2.4 with footnotes 3-4; Definitions 2.1 (p.7), 2.2 (p.8),
  3.6 (p.14); Lemma 3.7 (p.14); Corollary 3.8, Theorem 3.9 (p.15); Section 3.3;
  Theorem 4.1 and proof of Theorem 1.1 (pp.18-19).
  https://eccc.weizmann.ac.il/report/2021/023/download/ ;
  SHA-256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062.
- BK-preprint: Eli Ben-Sasson and Swastik Kopparty, "Affine dispersers from subspace
  polynomials", author preprint 2010-03-12. Used: Definition 2.1 and Theorem 2.2 with
  equation (4), printed pp.6-7 (PDF pp.7-8).
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf ;
  SHA-256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
- Shoup-preprint: Victor Shoup, "New algorithms for finding irreducible polynomials
  over finite fields", author preprint 1989-01-31 (Math. Comp. 54 (1990) 435-447).
  Used: Theorem 3.2 and its fixed-prime consequence, printed p.3.
  https://www.shoup.net/papers/detirred.pdf ;
  SHA-256 00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644.

The digests identify the exact versions consulted for this revision.
