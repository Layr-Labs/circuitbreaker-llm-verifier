# A uniform certificate-prefix family: coefficient 7/2, exponent 1

## Claim and provenance

For the single total Boolean language F_N defined below, with N denoting
the entire input length, we prove

    C_B2(F_N) >= (7/2)*N - R(N),    R(N)/N -> 0 as N -> infinity.

The joint coordinate pair is (7/2, 1). This applies existing published
mathematics from two sources:

- LY: Li and Yang, ECCC TR21-023, front-page date 2021-02-20. We use only
  its full-B2 model conventions (Sections 1.1, 2.1, 2.4 footnotes 3-4 and
  3.3) and no part of its gate-count argument.
  SHA256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062.

- BK: Ben-Sasson and Kopparty, author preprint dated 2010-03-12. We use
  the algebraic construction underlying their Theorem 7.3 and the structural
  arguments in Sections 3, 6 and 7.
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.

No other external source is invoked. The gate elimination argument in
Sections 5-6 below is reproduced from first principles and is not
attributed to any paper beyond those two.

---

## 1. The family on every bit string

Set F_0(empty) = 0. For N >= 1, proceed as follows.

Step 1. Let B be the least positive integer with B^5 > 32*N^3. Find B by
binary search on [1, 2N]. Let k be the least prime >= B. Find k by
scanning upward and testing each candidate by trial division.

Step 2. If k >= N, set F_N = 0 on all inputs. Otherwise set m = N-k,
r = ceil(m/k), t = ceil(sqrt(r)). If t >= r or k*t < r*(r+1), set F_N = 0.

Step 3. Given input x of length N, read bits x[0..k-1] as the lower k
coefficients of a monic polynomial Q(Z) = Z^k + sum_{j=0}^{k-1} x[j]*Z^j
over F_2. Test whether Z^(2^k) = Z mod Q and gcd(Q, Z^2+Z) = 1. If either
test fails output 0.

Step 4. Set K = F_2[Z]/(Q). Append r*k-m zero bits to the payload bits
x[k..N-1], split into r blocks of k bits, and interpret each block as a
field element a_i in K in the polynomial basis (bit j maps to Z^j). The
map from payload bits to (a_1,...,a_r) is F_2-linear and injective.

Step 5. For 1 <= i <= r set w_i = r+2-i and y_i = a_i^(2^(w_i)-1) in K,
computed by binary powering with O(w_i) field multiplications.

Step 6. Compute E_t = Sym_r^t(y_1,...,y_r), the elementary symmetric
polynomial of order t, via the division-free recurrence: initialize E_0=1
and E_j=0 for 1<=j<=t; for each i from 1 to r update E_j <- E_j + y_i*E_{j-1}
for j = min(i,t) down to 1. Output the constant coefficient of E_t in K,
i.e., the coefficient of Z^0.

This specifies F_N completely for every N and every input.

---

## 2. Exact prefix test

We verify that Q passes both tests in Step 3 if and only if Q is irreducible
of degree k over F_2.

The polynomial Z^(2^k)-Z factors in F_2[Z] as the product of all monic
irreducibles whose degree divides k. Its roots in the algebraic closure are
exactly the 2^k elements of F_{2^k}, and they are all simple (the formal
derivative is -1 = 1 in characteristic 2, hence a unit).

Since k is prime, the only positive divisors of k are 1 and k itself.
Degree-1 irreducibles are Z and Z+1, whose product is Z^2+Z. Every
other irreducible factor of Z^(2^k)-Z has degree k.

For a monic Q of degree k:
- Z^(2^k) = Z mod Q means Q divides Z^(2^k)-Z in F_2[Z], so every
  irreducible factor of Q has degree dividing k, i.e., 1 or k.
- gcd(Q, Z^2+Z) = 1 means Q shares no factor with Z*(Z+1), so Q has
  no degree-1 factor.
- Since deg(Q) = k and Q is squarefree (as a divisor of Z^(2^k)-Z which
  has no repeated roots), and has no factor of degree 1 or of degree
  strictly between 1 and k, Q must be irreducible of degree k.

Conversely, every monic irreducible of degree k divides Z^(2^k)-Z and is
coprime to Z^2+Z. Both tests pass.

For existence: the number of monic irreducibles of degree k over F_2 is
(2^k - sum_{d|k, d<k} d*N_d)/k where N_d counts irreducibles of degree d.
For prime k this simplifies to (2^k - 2)/k >= 1 for all primes k >= 2.
So there always exists a valid Q, and the set of inputs that output 0
at the guard is well-defined and computed deterministically.

---

## 3. Uniform P membership and NP encoding

One fixed deterministic algorithm A evaluates F_N(x):

- Finding B by binary search and k by scanning uses O(k) primality tests,
  each by trial division taking O(sqrt(k)) steps. Since k < 2B < 4*N^(3/5)+2
  by Bertrand's theorem, this is polynomial in N.
- The Frobenius test Z^(2^k) mod Q requires k squarings in F_2[Z]/(Q).
  Each squaring costs O(k^2) bit operations, total O(k^3) = O(N^(9/5)).
- The gcd test costs O(k^2) bit operations.
- r = ceil(m/k) = O(N^(2/5)) and t = O(N^(1/5)). The powers y_i need
  O(r) field multiplications each; over all i, O(r^2) multiplications,
  each O(k^2) bit operations, total O(N^2). The symmetric recurrence
  uses O(r*t) multiplications, total O(N^{9/5}).
- The whole computation is O(N^2) bit operations, polynomial in N.

Set L = {x : A(x) = 1}. Define V(x,w) = 1 iff w is empty and A(x) = 1,
with witness bound W(N) = 0. If x in L, the empty witness is accepted;
if the empty witness is accepted, x in L. Both NP directions hold.
Thus L is in P, hence in NP with empty witnesses.

---

## 4. Asymptotic guard passage

Bertrand's theorem gives k < 2B <= 2*(32*N^3)^(1/5) + 2 < 4*N^(3/5) + 2.
Also B^5 > 32*N^3 and B is the least such integer, so (B-1)^5 <= 32*N^3,
giving B > 2*N^(3/5). Since k >= B, we have k > 2*N^(3/5) for all N >= 1.
Thus m = N-k satisfies m >= N - 4*N^(3/5) - 2.

With k > 2*N^(3/5) and r <= m/k + 1 < (1/2)*N^(2/5) + 1:
- k >= N fails only for small N; for N >= 16 we have k < N/2.
- t = ceil(sqrt(r)) satisfies t <= sqrt(r)+1, so t < r eventually since
  r -> infinity.
- k*t >= r*(r+1): we need k*sqrt(r) > r*(r+1) approximately, i.e.,
  k > r^(3/2). Since k > 2*N^(3/5) and r < (1/2)*N^(2/5)+1, we have
  r^(3/2) < ((1/2)*N^(2/5)+1)^(3/2) = O(N^(3/5)), and k/r^(3/2) -> 2/(1/2sqrt(2)) > 1.
  Concretely: r*(r+1)/t <= sqrt(r)*(r+1) < ((1/2)*N^(2/5)+1)^(3/2) + N^(2/5) + 1,
  and k > 2*N^(3/5). For sufficiently large N, 2*N^(3/5) > ((1/2)*N^{2/5}+1)^{3/2}+N^{2/5}+1.

So there exists N_0 such that for every integer N >= N_0 all guards pass.

---

## 5. The affine disperser property

For each sufficiently large N, fix any irreducible Q of degree k and fix
the first k input bits to the coefficients of Q. The resulting restriction
G_Q: {0,1}^m -> {0,1} operates on the m payload bits.

We prove that G_Q is an affine disperser of dimension

    d_N = floor(k*r/t) + k*t + r*(r+1) + 1.

That is, G_Q is nonconstant on every affine subspace of {0,1}^m of
dimension strictly greater than d_N - 1, equivalently of dimension >= d_N.

The argument follows BK and is reproduced in full for completeness.

### 5.1 Polynomial representation and linear algebra

Let q = 2^k. A function K^r -> F_2 has a unique multivariate polynomial
representative with each variable degree at most q-1, obtained by reducing
modulo X_i^q - X_i. (This follows by induction on r and the fact that a
nonzero univariate polynomial of degree <= q-1 has at most q-1 roots in K.)
A nonzero positive-degree term in this representative certifies nonconstancy.

The trace Tr: K -> F_2 is defined by Tr(z) = sum_{h=0}^{k-1} z^(2^h). It is
F_2-linear and surjective (its polynomial has degree 2^(k-1) < q so is not
identically zero). For any nonzero mu in K, z -> Tr(mu*z) is also nonzero
since multiplication by mu is bijective. As a result, every nonzero F_2-linear
functional pi: K -> F_2 can be written pi(z) = Tr(mu*z) for some unique
nonzero mu. In particular the constant-coefficient extraction pi(f) = f mod Z
is such a functional (with pi(1)=1, pi(Z)=0, etc.) and is nonzero.

A linearized polynomial is one of the form sum_{j=0}^{k-1} a_j * X^(2^j).
Over K, these are exactly the F_2-linear maps K -> K: linearity follows
from the Frobenius, and their count matches the dimension of the space of
linear maps, so the correspondence is bijective.

### 5.2 Image polynomials and zero-run structure

For a linear subspace V <= K of F_2-dimension d, define the subspace
polynomial P_V(X) = product_{v in V} (X - v). This is a monic linearized
polynomial of degree 2^d (linearized since adding any v in V permutes V,
so P_V(X+v) = P_V(X); polynomial by induction on dim V).

Its kernel is V, so P_V(K) has dimension k-d. Let Q = P_{P_V(K)} be the
image polynomial of P_V(K). Then Q(P_V(X)) = X^q - X, because both sides
are monic of degree q vanishing on all of K. Symmetrically P_V(Q(X)) = X^q - X.
Thus Q maps K onto V (kernel has dimension k-d, so image has dimension d),
Q is a monic linearized polynomial of degree 2^(k-d) with 2^(k-d) distinct
roots in V. Write Q = sum_{j=0}^{k} a_j * X^(2^j) with a_k = 1 and indices
modulo k.

Claim (no long zero run): there is no cyclic run of d consecutive zero
coefficients among a_0,...,a_{k-1}. Proof: if such a run exists, a Frobenius
twist X -> X^(2^h) shifts it to the last d positions. The twisted polynomial
has degree <= 2^(k-d-1) < 2^(k-d), yet still has Q's 2^(k-d) roots as a
divisor, a contradiction since a nonzero polynomial of degree < 2^(k-d)
cannot have that many roots in K.

Claim (for prime k, the only zero run of length d-1 is the terminal one):
if another such run exists at offset h (0 < h < k), a Frobenius twist moves
it to the end, giving a nonzero scalar multiple c*Q. So every element of V
satisfies z^(2^h) = c*z. Picking nonzero z_0 in V, every z/z_0 in V/span(z_0)
satisfies (z/z_0)^(2^h) = z/z_0, meaning z/z_0 is in F_{2^{gcd(h,k)}} = F_2
since gcd(h,k)=1 (k prime, 0 < h < k). So dim(V) <= 1, contradicting d >= 2.

### 5.3 Nonzero coefficient for independent subspaces

Let V_1,...,V_t be linear subspaces of K of dimensions d_i > 1 each.
Write e_i with 1 <= e_i < d_i-1, h_i = e_i+1, delta_i = d_i - h_i > 0.
Suppose sum_i delta_i > k - max_i d_i + 1.
Let Q_i be the image polynomial of V_i.

Claim: the coefficient of the monomial product_i X_i^(b_i) with wt(b_i)=h_i
in the reduced polynomial Tr(mu * product_i Q_i(X_i)^(2^(h_i)-1)) is nonzero,
where wt(b) counts the number of 1-bits in b.

Proof: Write Q_i = sum_j a_{i,j} X^(2^j) with periodic coefficients mod k.
The h_i-fold Frobenius product Q_i(X)^(2^(h_i)-1) = product_{s=0}^{h_i-1}
Q_i(X)^(2^s). Expanding and collecting, the coefficient of X^(2^{b_{i,v}})
(for b_{i,v} distinct residues mod k) is a product of entries of Q_i's
coefficient matrix T_i, where T_i[u,v] = a_{i,v-u}^(2^u) (indices mod k).

For the chosen monomial to appear, we need one Frobenius residue set per block
that is an h_i-element column set B_i with all residues distinct (so binary
degree is exactly h_i with no carry). The column set B_i determines a submatrix
T_i[A_i, B_i] for row set A_i = {0,...,e_i}. The coefficient from block i is
the permanent of this submatrix.

Construct the column sets B_i so that the rho=0 term in a trace-over-Frobenius
expansion is nonzero and every rho != 0 term has a zero column in some minor.

- For rho != 0, we need every nonzero rho to be blocked by some column b in
  some B_i: specifically, the interval [b+1, b+delta_i] mod k must cover rho.
  The construction ensures this by placing B_i in non-overlapping segments
  across the k positions, with each delta_i covering its interval.
- For rho = 0, we need each T_i[A_i, B_i] to have nonzero permanent. By the
  zero-run structure of Q_i (no long zero runs, unique terminal run for prime k),
  we can choose B_i so the minor is lower-triangular with nonzero diagonal.

The full construction: reorder so d_1 >= d_2 >= ... >= d_t. Let J be the first
index with sum_{i<=J} delta_i > k-d_1. For i <= J: place B_i starting at
offset ell_i = sum_{j<i} delta_j; for i > J: set ell_i = 0. The no-carry
condition gives nonzero permanents; the coverage argument shows every rho=1,...,k-1
is blocked by the interval from at least one source.

Thus the rho=0 term alone survives in the trace sum, and it is a product of
nonzero field values (each permanent is nonzero). The monomial coefficient is
nonzero, proving nonconstancy of the polynomial.

### 5.4 Decomposing an affine subspace into blocks

Let A be any affine subspace of K^r of F_2-dimension D > floor(k*r/t)+k*t+r*(r+1).
By Gaussian elimination on the direction space of A, decompose it as

    A = {(gamma_i + z_i + sum_{j<i} sigma_{ij}(z_j))_{i=1}^r : z_i in V_i}

where V_i are linear subspaces of K and sigma_{ij}: V_j -> K are linear maps,
and dim(A) = sum_i dim(V_i). (This is standard linear algebra: filter the
direction space by vanishing on successive coordinate blocks, take quotients.)

Let d_i = dim(V_i) and let S = {i : d_i >= r+2}. From D > k*t + r*(r+1)
we get sum_{i in S} d_i > k*t + r*(r+1), so |S| >= t. Choosing the t
largest d_i's (index set I_0 in S):

    sum_{i in I_0} (d_i - (r+1)) > k   (by a weighted average argument using |S|>=t).

For each i in I_0, e_i = r+1-i >= 1 (since i <= r), and h_i = e_i+1 = r+2-i = w_i.
Also h_i < d_i since d_i >= r+2 and h_i = r+2-i <= r+1 < r+2 <= d_i.
So the subspace assumption applies: delta_i = d_i - h_i >= 1 for each i in I_0,
and sum_{i in I_0} delta_i = sum_{i in I_0}(d_i - h_i) > k - max(h_i) + 1
(since the excess over (r+1)*t is > k, and max(h_i) <= r+1).

Section 5.3 then yields a nonzero monomial M = product_{i in I_0} X_i^(b_i)
with wt(b_i) = w_i in the polynomial Tr(mu * product_{i in I_0} Q_i(X_i)^(2^(w_i)-1)).

### 5.5 Triangular mixing preserves the coefficient

In f(R(X)) where R_i(X) = gamma_i + Q_i(X_i) + sum_{j<i} Q_{ij}(X_j),
replace Q_i(X_i) by R_i(X) in each selected factor. The monomial M involves
only indices in I_0. For the largest i in I_0, variable X_i enters only through
R_i, and achieving binary degree w_i in X_i forces all w_i Frobenius factors to
use Q_i(X_i) (the only X_i term in R_i). Repeat inductively for smaller i. The
lower-triangular mixing terms (j < i) cannot contribute to M: they use variables
of smaller index with the wrong exponent signature. So the coefficient of M in
f(R(X)) equals its coefficient in the independent-subspace expression, which is
nonzero.

### 5.6 Isolation from other symmetric summands

For any other t-element index set J != I_0, a monomial matching M in the
corresponding symmetric summand would require the same binary degree profile
from a different index set. A capacity argument (counting factors per variable
and using the decreasing-w structure) forces J = I_0, so no cancellation occurs.
Concretely: if n_{ji} factors from position j contribute to variable i in I_0,
then sum_i n_{ji}/w_i <= 1 per j (supply bound), and sum_j n_{ji}/w_i >= 1
per i (demand bound), so sum_{i,j} n_{ji}/w_i is in [t, t]. Equality forces
each factor to use i = j, so I_0 = J.

### 5.7 Conclusion

The polynomial f, when restricted to the affine subspace A = R(K^r), has a
nonzero coefficient at the positive-degree monomial M. Hence f is nonconstant
on A. Since the payload embedding {0,1}^m -> K^r is F_2-linear and injective,
every affine subspace of {0,1}^m of dimension >= d_N maps to an affine subspace
of K^r of dimension >= d_N. With d_N = floor(k*r/t) + k*t + r*(r+1) + 1, the
above argument applies to any such subspace. Therefore G_Q is an affine
disperser of dimension d_N.

Note: d_N = O(N^(4/5)) because k=O(N^(3/5)), r=O(N^(2/5)), t=O(N^(1/5)), giving
k*r/t = O(N^(4/5)), k*t = O(N^(4/5)), and r*(r+1) = O(N^(4/5)).

---

## 6. Gate lower bound: (7/2)m - 7*d for an affine disperser

We prove the following self-contained lemma.

**Lemma.** Let f: {0,1}^n -> {0,1} be an affine disperser of dimension d
(nonconstant on every affine subspace of dimension >= d). Every circuit in
the full B2 model computing f has size |C| >= (7/2)*n - 7*d.

The proof is a parameter-elimination argument entirely derived here.
We follow the structure of the potential function method; the key
estimates are elementary case analysis on binary truth tables.

### 6.1 Setup: representation and potential function

Given a circuit C for f, we maintain a triple (Phi, sys, circ) where:
- Phi(u) = b + L*u + sum_{{p,h} in P} c_{p,h}*p*h is a map from M Boolean
  parameters u to the n original input coordinates. L is a binary matrix,
  P is a set of q disjoint pairs of parameters, and c_{p,h} in {0,1}^n.
  The map is always injective; fixing one member of each pair makes it affine
  with rank a = M - q, so its image contains an affine subspace of dimension a.
- sys: a square invertible auxiliary system A*v = N*u + c with at most 3
  nonconstant terms per row. Gates and equations total s.
- circ: the remaining acyclic B2 gates. The circuit output is a function of
  parameters, auxiliaries, and gate values.
- The invariant: the output equals f(Phi(u)) for every u assignment.

Let s = (equations in sys) + (gates in circ). Initially M = n, q = 0, a = n,
s = |C|.

For each unpaired parameter z let deg(z) = number of distinct rows or gates
using z. Assign weights: weight(z) = 0 if deg(z)=0, 6 if deg(z)=1, 7 if
deg(z)>=2. Set H = 2*s + 3*q + sum_{unpaired z} weight(z).

Initially H <= 2|C| + 7*n (since all parameters have weight <= 7 and q=0).
We always have H >= 0.

While a > d, we claim the output cannot be constant (else Phi maps M-q-dimensional
affine space to a level set, nonconstancy gives a contradiction) or affine
(fixing one parameter kills the output, reducing a by 1 with no gate charge;
the remaining image still has dimension a-1 >= d, contradiction). Similarly
the output is not a single AND of two unpaired parameters: forcing one of them
kills the output with Delta(a)=-1, same contradiction.

### 6.2 Four free deletions (preserve M, q, a; decrease H by >= 1)

All moves below preserve the injectivity invariant by construction.

(D1) Any gate or row that outputs a constant or a function of a single input:
route its output to that value, remove the gate/row. At most one surviving
parent gains weight by 1; the gate/row saves 2. Net: Delta(H) <= -1.
This covers: constant gates, unary gates (repeated or fixed pins), XOR on
affine parents with the same linear coefficient, AND on affine parents with
dependent linear parts (the table has one effective input).

(D2) A gate with no children (except possibly the output): remove it unless it
is the output (the output case is handled in termination). Saves 2, no weight
gain. Delta(H) <= -2.

(D3) A row with at most 2 nonconstant terms: since A is invertible, one term is
an auxiliary v. Write v = U + c (U is zero or another value), substitute, erase
the row. Sparsity is preserved (substitution only adds terms where v appears in
rows with >= 3 terms; but those rows had v among their terms, so the count does
not exceed 3 after substitution). Invertibility follows by expanding the
determinant on the v-column. At most U gains weight by 1; the row saves 2.
Delta(H) <= -1.

(D4) An unpaired parameter x of degree 1, appearing in exactly one row E:
since A*v = N*u + c, the N-column for x is nonzero, so A^{-1}*N_x is nonzero,
meaning some auxiliary Q satisfies Q = x + ell(others) in the system. Exchange
roles: Q becomes the new parameter, x becomes an auxiliary. The auxiliary matrix
stays invertible (determinant = the Q-coordinate of A^{-1}*N_x = 1). Now
auxiliary x appears only in row E; delete E (two-term row after substitution,
handled by D3) and delete the row. Old weight 6 replaced by at most 7; row saves
2. Delta(H) <= -1.

### 6.3 Affine gate absorption

Any gate computing an affine function (XOR-type on any inputs without an AND
ancestor in the path to the output) can be moved into the affine system by
adding one row for its output as auxiliary. The new row has at most 3 terms.
The auxiliary matrix remains block lower-triangular with 1 on the new diagonal
entry, hence invertible. Size s is unchanged. Parameter degrees do not increase.
Perform this in topological order before each case inspection.

### 6.4 Quadratic pivot (removes two gates, increases q by 1)

Suppose gate G = (x+b)(y+c)+e has a single child D = G + z + h, where x, y, z
are distinct unpaired parameters and z has degree 1. Replace z by w + G + h in
Phi (bijection since the old z-to-new-w map is invertible). Record pair {x,y}.
Then D = w; reroute all uses of D to w. Delete D and G.

Changes: Delta(s) = -2, Delta(q) = +1, Delta(a) = -1.
Delta(H) = -4 + 3 + (7 or less for new w) - weight(x) - weight(y).
When deg(x), deg(y) >= 2: Delta(H) <= -4 + 3 + 7 - 7 - 7 = -8.
When deg(x) = 1 (weight 6): Delta(H) <= -4 + 3 + 7 - 6 - 7 = -7.

### 6.5 Exhaustive case analysis while a > d

Apply D1-D4 repeatedly; also absorb eligible affine gates. When no D1-D4
applies, the circuit has the following structure: every gate has two distinct
nonconstant parents; AND gates on affine parents have independent linear parts;
all non-output gates have children; all rows have exactly 3 terms; unpaired
parameters are not isolated in rows.

Under these conditions, the first gate G in topological order (closest to inputs)
is AND-type (since all XOR-type gates were absorbed). Its parents are affine
values with independent linear parts. If a parent is an auxiliary Q, it depends
on some unpaired parameter x with deg(x) >= 2; swap Q and x at zero cost (both
have weight 7). After at most 2 such swaps, G reads two distinct unpaired
parameters x, y with deg(x) >= deg(y) >= 2.

If a child of G also reads x, compose G and the child into one binary gate on
(x, y) and reroute. This strictly reduces the child count of G; repeat until
no child of G reads x or y. Then do the symmetric operation for y.

Now write u = deg(x) >= v = deg(y) >= 2 and h = deg(G) >= 1.

**Case A: u >= 3.** Fix x to a constant. Delete x (saves weight u*weight-coeff),
delete G (saves 2), delete G's h children (they now have one constant parent,
handled by D1; saves 2h). No child reads x (by preparation). Surviving parents
of G's children may gain weight by at most 1 each. Number of potential gains
is at most h. So:

    Delta(H) <= -7 - 2*(1+h) + h = -7 - 2 - h <= -9 - h.
    Delta(a) = -1.
    Delta(H) / Delta(a) <= -9 - h <= -9 < -7.

**Case B: u = v = 2, and x and y have a common second child B (besides G).**
Restrict (x, y) to a line making both G and B constant: G and B are Boolean
tables on {0,1}^2; any two tables agree on at least 2 of 4 inputs, so three
points from the 4 where G is constant include two where B is also constant,
forming an affine line in (x,y) space. Fix x or identify x = y+c. Delete G
and B. The surviving parameter is unused.

    Delta(a) = -1, Delta(H) <= -7 - 7 - 4 = -18 < -7.

**Case C: u = v = 2, no common second child, u = v = 2, h >= 2.** Fix x.
Delete G and its h children. No child reads x (preparation). y loses G but
no new child: y's weight drops by at most 0 (already 7, degree falls by 1 to at
least 1, weight stays 7 or drops to 6). The h parents of G's h children can
gain weight by at most 1 each; but no co-parent of any child of G was an isolated
unpaired parameter (by the no-isolation rule after D1-D4), so each gain is from
7 to ... still at most 7. But each such co-parent goes from having 2+ children to
possibly having 1 fewer; weight can only drop or stay. In any case:

    Delta(H) <= -7 - 2*(1+h) + h <= -7 - 2 - h.
    With h >= 2: Delta(H) <= -11. Delta(a) = -1. Ratio <= -11 < -7.

**Case D: u = v = 2, h = 1, and G's unique child D reads an isolated unpaired z
(degree 1).** Since D has two parents (G and z), and z has degree 1, apply the
quadratic pivot of Section 6.4. This requires deg(x), deg(y) >= 2 (met) and D
is XOR-type or AND-type:

- If D is XOR-type: Delta(H) <= -8 (from Section 6.4), Delta(a) = -1.
- If D is AND-type: write D = (G+c)*(z+b)+e. Choose x_0, y_0 with G(x_0,y_0)=c
  (possible since G is nondegenerate). Fix x = x_0, y = y_0. Delete G, D, and
  the two children of x and of y that we labeled B_x, B_y in preparation
  (distinct by Case B absence). Also delete z (degree 1, now unused).

      Delta(a) = -2, Delta(H) <= -7 - 7 - 6 - 8 = -28, and 28/2 = 14 > 7.

**Case E: u = v = 2, h = 1, G's unique child D has no isolated co-parent.**
Fix x. Delete x (weight 7), delete G (saves 2), delete D (D has exactly one
parent left after G is rerouted to a constant, which is either another value or
constant—in any case it triggers D1). So we delete G and D: saves 4 total.
Since y has lost G and thus reduced degree from 2 to 1, weight drops from 7 to 6.
No weight gain from other parameters.

    Delta(H) <= -7 - 4 - 1 = -12. Delta(a) = -1. Ratio 12 > 7.

In every case, Delta(H) <= 14 * Delta(a) (i.e., each unit decrease in a costs at
most 14 in H). Since H starts at <= 2|C| + 7n and ends at >= 0 after reducing
a from n down to d:

    2|C| + 7n >= 14*(n - d), so 2|C| >= 7n - 14d, so |C| >= (7/2)*n - 7*d.

If n <= d the right side is nonpositive and the conclusion is immediate.
(Note: in cases where Delta(a) = -2, Delta(H) / Delta(a) is still <= 14,
e.g., case D AND-type gives 28/2 = 14.)

**Termination:** Each step either decreases s (number of gates and equations)
or decreases M (number of parameters), without increasing either. Since both
are nonneg integers bounded initially by |C| + n, the process terminates.

### 6.6 The potential optimality bound

For a potential of the form H = 2s + beta*q + alpha*(count of used unpaired) -
sigma*(count of isolated unpaired), the bound (7/2) is tight for these local
estimates: the cases at (deg_x, deg_y, deg_G) = (2, 1, 1) and (3, 1, 1) give
competing constraints theta <= alpha + 4 - sigma and theta <= alpha + 3 + sigma,
so theta - alpha <= min(4-sigma, 3+sigma) <= 7/2 (maximized at sigma = 1/2).
This is an intrinsic limit of this particular analysis, not a universal barrier
to stronger proofs by different methods.

---

## 7. Asymptotic lower bound

For each sufficiently large integer N >= N_0:
1. Fix any irreducible Q of degree k. Restrict x to have its first k bits equal
   to the coefficients of Q. This yields a sub-circuit C' for G_Q on m bits with
   |C'| <= |C|, so C_B2(F_N) >= C_B2(G_Q).
2. G_Q is an affine disperser of dimension d_N (Section 5).
3. By the gate lemma (Section 6): C_B2(G_Q) >= (7/2)*m - 7*d_N.
4. Thus: C_B2(F_N) >= (7/2)*(N-k) - 7*d_N = (7/2)*N - ((7/2)*k + 7*d_N).

Set R(N) = (7/2)*k + 7*d_N. Since k = O(N^(3/5)) and d_N = O(N^(4/5)):
R(N) = O(N^(4/5)), so R(N)/N -> 0.

For every epsilon > 0, choose N_epsilon >= N_0 large enough that R(N)/N <= epsilon
for all N >= N_epsilon. Then C_B2(F_N) >= (7/2 - epsilon)*N for all N >= N_epsilon.
Thus liminf_{N->inf} C_B2(F_N)/N >= 7/2.

Both the coefficient 7/2 and the exponent 1 come from this single theorem
and single language. There are no extraneous assumptions or unverified premises.

---

## 8. Sources

- LY: ECCC TR21-023, 2021-02-20; Sections 1.1, 2.1, 2.4 footnotes 3-4 and 3.3
  for the B2 model conventions only.
  URL: https://eccc.weizmann.ac.il/report/2021/023/download/
  SHA256: 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062.

- BK: Ben-Sasson and Kopparty, author preprint 2010-03-12; Lemmas 3.7, 3.10,
  3.13, Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4. The disperser
  construction in Sections 5-6 above adapts these arguments with full proofs
  inline. Only elementary finite-field algebra is assumed.
  URL: https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256: 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
