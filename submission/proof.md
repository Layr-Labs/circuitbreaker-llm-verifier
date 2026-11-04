# Coefficient four: certificate prefix, affine disperser, gate elimination

This package proves one theorem about one fixed total language in P, hence in NP:
in the pinned full B2 model, for every integer length N >= 2^20,

    C_B2(f_N) >= 4N - R(N),   R(N) = 4k + 9(floor(kr/t) + kt + r(r+1) + 1),

with k, r, t defined below, and R(N) <= 58*N^(4/5). Hence
liminf_{N -> infinity} C_B2(f_N)/N >= 4. The declaration in `bound.json`
records the version-1 polynomial/log guarantee with coefficient 4, exponent 1,
log exponents 0, and the finite projection [4,1].

Two instruments are proved in the accompanying evidence files and are not
re-derived here: `disperser.md` proves the algebraic separation lemma (the
binary case of Ben-Sasson--Kopparty, reproduced with attribution), and
`elimination.md` proves the ternary-pair gate-elimination theorem (adapted
from Mishra--Santhanam, with attribution). This file fixes the family, proves
its eligibility, verifies the number-theoretic guards at every sufficiently
large integer length, and combines the two instruments.

## 1. The fixed total family

Let N = |x| be the actual input length. Define f_0(empty) = 0. For N >= 1:

1. Let B_0 be the least positive integer B with B^5 > 32*N^3, and let k be
   the least prime with k >= B_0. Since B_0 = floor(2*N^(3/5)) + 1, this is a
   deterministic function of N.
2. Read the first min(k,N) bits of x as q_0,...,q_{k-1}, with zeros appended
   when N < k, and set
       Q(Z) = Z^k + q_{k-1} Z^(k-1) + ... + q_1 Z + q_0  over F_2.
3. Run the fixed deterministic Frobenius/gcd irreducibility test of Section 2
   on Q. If Q is not irreducible of degree k, output 0.
4. Put m = N - k. If m <= 0, output 0. Put r = ceil(m/k) and
   t = ceil(sqrt(r)). If 1 <= t < r and k*t >= r(r+1) do not both hold,
   output 0.
5. Let K = F_2[Z]/(Q), with the polynomial basis 1,Z,...,Z^(k-1). Append
   r*k - m zeros to the payload bits x_k,...,x_{N-1} and split them into r
   consecutive blocks of k bits. Let a_i in K be the element whose polynomial
   basis coefficients are block i, for i = 1,...,r.
6. Let pi:K -> F_2 take the constant coefficient, and put w_i = r+2-i.
   Output
       f_N(x) = pi( sum_{|I|=t} prod_{i in I} a_i^(2^(w_i) - 1) ),
   the sum being over all t-element subsets I of {1,...,r}.

This defines one function f_N on {0,1}^N for every N, hence one total
language L = {x : f_{|x|}(x) = 1}. Nothing depends on advice, on the length
in a nonconstructive way, or on a witness.

## 2. Uniform polynomial time and NP membership

The prime scan is polynomial in N: B_0 <= 2*N^(3/5) + 1, and by Bertrand's
postulate a prime occurs before 2*B_0, so at most polynomially many candidates
are tested. Each candidate is tested by trial division up to its square root,
which is polynomial in N.

The Frobenius/gcd test is polynomial in k. Compute Z^(2^k) mod Q by k
successive squarings of polynomials of degree less than k over F_2, each using
O(k^2) bit operations. Because k is prime, the standard criterion says Q is
irreducible of degree k exactly when gcd(Q, Z^2 - Z) = 1 and
Z^(2^k) = Z mod Q: Q divides Z^(2^k) - Z, so its distinct irreducible factors
have degrees dividing the prime k; the gcd condition excludes degree-one
factors, leaving degree k. Both the squarings and the gcd are polynomial in k.

The symmetric polynomial is evaluated without characteristic-two division:
start E_0 = 1 and E_j = 0 for 1 <= j <= t, and after reading the value
b_i = a_i^(2^(w_i)-1) update E_j = E_j + b_i * E_(j-1) in descending j from
min(i,t) down to 1. Splitting subsets by inclusion of i proves the recurrence.
Each b_i is computed by binary powering with O(w_i) <= O(r) field
multiplications, and the recurrence uses O(r*t) field multiplications. Field
operations in the polynomial basis cost O(k^2) bit operations. Since
k = O(N^(3/5)), r = O(N^(2/5)) and t = O(N^(1/5)), the whole evaluator is
polynomial in N.

Thus one deterministic polynomial-time algorithm A computes every f_N. Define
V(x,w) = 1 iff w is empty and A(x) = 1. Both directions are immediate, so L
is in P and hence in NP with witness bound zero. Witness bits are never
circuit inputs.

## 3. Guards at every sufficiently large integer length

Let N >= 2^20. Then N^(1/5) >= 16, N^(2/5) >= 256 and N^(3/5) >= 4096.

(i) B_0 = floor(2*N^(3/5)) + 1 satisfies 2*N^(3/5) <= B_0 <= 2*N^(3/5)+1.
    Since k is prime and k >= B_0 > 2*N^(3/5), and Bertrand gives a prime
    below 2*B_0, we have
        2*N^(3/5) < k < 4*N^(3/5) + 2 < N/2,
    the last inequality because 8*N^(3/5) + 4 < N when N^(2/5) >= 256.
    Hence m = N - k > N/2 > 0.

(ii) r = ceil(m/k) satisfies m/k <= r < m/k + 1 <= N/k + 1
     < N/(2*N^(3/5)) + 1 = N^(2/5)/2 + 1 <= 0.504*N^(2/5).
     Also r >= m/k > (N/2)/(4*N^(3/5)+2) >= N^(2/5)/8.01 >= 31, so r >= 32.

(iii) t = ceil(sqrt(r)) satisfies 2 <= t <= sqrt(r) + 1 < r, because
      sqrt(r) + 1 < r for r >= 3. Hence 1 <= t < r.

(iv) kt >= k*sqrt(r) >= 2*N^(3/5) * sqrt(N^(2/5)/8.01)
     = (2/sqrt(8.01))*N^(4/5) > 0.70*N^(4/5), while
     r(r+1) <= (0.504*N^(2/5))(0.504*N^(2/5)+1) < 0.26*N^(4/5) + 1
     <= 0.27*N^(4/5) for N^(2/5) >= 256. Hence kt >= r(r+1).

Therefore every guard in Section 1 holds at every integer N >= 2^20. No
length is skipped, and the definitions at smaller lengths or on failed guards
are total but carry no asymptotic claim.

## 4. Each fixed-prefix restriction is an affine disperser

Fix N >= 2^20 and a valid prefix Q. The payload bits x' in F_2^m are mapped
by zero-padded block splitting to (a_1,...,a_r) in K^r. This map P is F_2-
linear and injective. The payload function g(x') = f_N(Q, x') is exactly the
function of the lemma in `disperser.md`:

  Lemma. Let K be the field of size 2^k with k prime. Suppose 1 <= t < r and
  k*t >= r(r+1). For any nonzero linear pi:K -> F_2 put w_i = r+2-i and
  f(a_1,...,a_r) = pi( sum_{|I|=t} prod_{i in I} a_i^(2^(w_i)-1) ). Then f is
  nonconstant on every affine subspace of K^r of dimension strictly greater
  than kr/t + kt + r(r+1).

The constant-coefficient functional pi is nonzero (pi(1) = 1). By Section 3
the hypotheses hold for every N >= 2^20. Set
d(N) = floor(kr/t) + kt + r(r+1). Let A be an affine subspace of F_2^m of
dimension strictly greater than d(N). Dimensions are integers and
floor(kr/t) + 1 > kr/t, so
dim(A) >= d(N) + 1 > kr/t + kt + r(r+1). The padding map P is F_2-linear and
injective, so P(A) is an affine subspace of K^r of the same dimension, and
the lemma puts f nonconstant on P(A). Hence g is nonconstant on A. Therefore:

  (*) g is nonconstant on every affine subspace of F_2^m of dimension
      strictly greater than d(N).

This is the only place where the algebraic separation lemma is used. It does
not depend on the circuit, on acceptance history, or on any unproved
construction: the prefix is fixed first, and the field, the exponents and the
functional are then fixed by the input.

## 5. Gate elimination on the restriction

`elimination.md` proves the following theorem in the pinned full B2 model.

  Theorem. Let d be a nonnegative integer and let h:F_2^n -> F_2 be nonconstant
  on every affine subspace of dimension strictly greater than d. Then every
  circuit for h in the full binary basis has at least 4n - 9(d+1) gates. This
  covers arbitrary sharing, repeated pins, constants, unary tables, degenerate
  gates and a single terminal output.

Apply the theorem with n = m and d = d(N) to the payload function g of
Section 4. By (*) its hypothesis holds, so

    C_B2(g) >= 4m - 9(d(N)+1).

Now let C be any circuit for f_N on N input bits. Fix its first k input bits to
the bits defining Q. The resulting circuit has exactly the same gates and
computes g on the remaining m bits: fixed inputs become constants, constant
nodes introduced by the substitution cost zero in the pinned model, and no
gate is added. Hence the number of gates of C is at least C_B2(g). Therefore

    gates(C) >= 4(N-k) - 9(floor(kr/t) + kt + r(r+1) + 1)
             = 4N - R(N),
    R(N) = 4k + 9(floor(kr/t) + kt + r(r+1) + 1).

The inequality holds for every circuit C for f_N, so it holds for C_B2(f_N).

## 6. Remainder and the all-length liminf

Using Section 3 and kr < m + k = N (because r = ceil(m/k) implies
kr < m + k):

- kr/t <= k*sqrt(r) = sqrt(k * kr) < sqrt(kN)
  < sqrt(4*N^(8/5) + 2N) < 2*N^(4/5) + N^(2/5).
- kt <= k*(sqrt(r)+1) < (4*N^(3/5)+2)(N^(1/5)/sqrt(2)+2)
  <= 2.83*N^(4/5) + 8*N^(3/5) + sqrt(2)*N^(1/5) + 4
  <= 3.34*N^(4/5), using N^(3/5) <= N^(4/5)/16 and N^(1/5) <= N^(4/5)/65536.
- r(r+1) < (N^(2/5)/2 + 1)(N^(2/5)/2 + 2)
  <= 0.26*N^(4/5) + 1 <= 0.27*N^(4/5).
- 4k < 16*N^(3/5) + 8 <= 1.01*N^(4/5).

Adding the three bracket terms and the constant one gives
d(N) + 1 <= 5.61*N^(4/5), so 9(d(N)+1) <= 50.5*N^(4/5), and therefore
R(N) <= 51.5*N^(4/5) <= 58*N^(4/5).

Consequently, for every integer N >= 2^20,

    C_B2(f_N)/N >= 4 - 58*N^(-1/5).

The right side tends to 4. Thus for every real epsilon > 0, every integer
N >= N_epsilon = max(2^20, ceil((58/epsilon)^5)) satisfies
C_B2(f_N) >= (4 - epsilon)*N. Therefore

    liminf_{N -> infinity} C_B2(f_N)/N >= 4,

which is exactly the version-1 polynomial/log guarantee (coefficient 4,
exponent 1, log exponents 0) in `bound.json`. The declared finite projection
[4,1] is implied by that guarantee.

## 7. Model accounting

The bound is in Li--Yang's full B2 model: deterministic acyclic circuits, one
output, two input pins per gate, all sixteen binary truth tables, unrestricted
depth and fan-out, exact computation on every input, and nonuniform circuits
allowed. Inputs are x alone; witness bits and the polynomial Q's validation
are not circuit inputs. Gates are counted once, including trivial, constant or
degenerate gates; inputs, wires and output selection cost zero. Repeated pins
and shared subgraphs are allowed. Edge negations are absorbed in the child
truth table, and the elimination theorem's normalizations delete unary and
constant gates without charging a new basis. The elimination inequality
bounds the original gate count, not a potential mistaken for one. The
restriction in Section 5 fixes input bits only; it never adds a gate and never
uses a witness. Every quantifier is over the actual integer input length, and
the threshold N_epsilon depends only on epsilon.

## 8. Sources and attribution

- `disperser.md` reproduces the binary Ben-Sasson--Kopparty separation
  argument, including the linearized-polynomial, permanent and block-
  decomposition lemmas, with the attribution recorded there. The cited
  construction and theorem locators are in `proof.md`'s companion evidence.
- `elimination.md` adapts the fair affine system, role exchanges and deletion
  machinery of Mishra--Santhanam (ECCC TR26-260, 2026-10-19, Sections 2-4) and
  proves the ternary-pair row table used in Section 5, with the attribution
  recorded there.
- This file claims no literature-seed coefficient and no result beyond the
  inequality of Sections 5 and 6. The two instruments are joined only by the
  fixed-prefix restriction of Section 4 and the polynomial-time evaluator of
  Sections 1 and 2.
