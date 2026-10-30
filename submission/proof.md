# C_B2(f_n) >= (7/2) n - 56 n^{4/5} - (35/2) n^{3/5} - 7 for an explicit family in P: complete proof

Declared: c = 7/2, p = 1, log exponents 0. Proof files: this file (family, eligibility, model, assembly), `elimination.md` (Part E: the gate-elimination theorem 3.5n' - 7d for every affine disperser on n' bits, re-derived from Mishra–Santhanam in the pinned model) and `disperser.md` (Part D: the affine disperser property, re-derived from Ben-Sasson–Kopparty). Everything is proved in these three files; no citation carries a step. The field modulus is read from the input (input-prefix modulus construction, credited in Section 5), so no irreducible-polynomial construction algorithm is needed.

## 1. The family

**Length 0.** f_0() := 0.

**Length n >= 1.** Define, as functions of n only:
- k = k(n): the least prime k with k^5 > 32 n^3, i.e. the least prime exceeding 2 n^{3/5}. (k >= 3.)
- If n <= k: f_n(x) := 0 for all x. Otherwise the *payload length* is m := n - k >= 1, r = r(n) := ceil(m / k) >= 1, t = t(n) := ceil(sqrt(r)); 1 <= t <= r.
- For i in [1, r]: alpha_i := 2^{r+2-i} - 1.

Given x in F_2^n with n > k, write x = (b_0, ..., b_{k-1}, y) with the *prefix* b in F_2^k and the *payload* y in F_2^m. Let Q(Z) := Z^k + sum_{j=0}^{k-1} b_j Z^j in F_2[Z], a monic polynomial of degree k.
- **Prefix test:** Q is *valid* iff Z^{2^k} ≡ Z (mod Q) and gcd(Q, Z^2 + Z) = 1. If Q is not valid, f_n(x) := 0.
- If Q is valid, let K := F_2[Z]/(Q) (a field with 2^k elements by Lemma 2.1), with polynomial basis 1, Z, ..., Z^{k-1}, and Tr : K -> F_2, Tr(z) := sum_{s=0}^{k-1} z^{2^s} (absolute trace; values in F_2, see `disperser.md`, preamble). Append rk - m zeros to y, split into r consecutive blocks of k bits, and map block i = (c_0, ..., c_{k-1}) to y_i := sum_j c_j Z^j in K. Then
    f_n(x) := Tr( Sym_r^t(y_1^{alpha_1}, ..., y_r^{alpha_r}) ),   Sym_r^t(w_1..w_r) := sum_{I subset [1,r], |I| = t} prod_{i in I} w_i.
The map iota_Q : y -> (y_1, ..., y_r) is F_2-linear and injective, and for a valid Q the restriction g_{n,Q} := f_n(b, ·) : F_2^m -> F_2 equals F o iota_Q with F as in Theorem D of `disperser.md` (with the field K = F_2[Z]/(Q)).

This is one total function at every length; nothing depends on n except through the computable quantities above, and nothing is chosen by any construction algorithm: the modulus is part of the input. (The payload function is the family of Ben-Sasson–Kopparty Theorem 2.2 with p = 2, pi = Tr and the blocks in reverse order; the prefix wrapper is the input-prefix modulus construction credited in Section 5.)

## 2. Eligibility: the prefix test, existence of valid prefixes, one polynomial-time evaluator

**Lemma 2.1 (prime-degree irreducibility criterion).** Let k be a prime and Q in F_2[Z] monic of degree k. Then Q is irreducible iff Z^{2^k} ≡ Z (mod Q) and gcd(Q, Z^2 + Z) = 1. If Q is irreducible, F_2[Z]/(Q) is a field with 2^k elements.
Proof. The last sentence is standard (a quotient of a polynomial ring over a field by an irreducible polynomial is a field; it has 2^{deg Q} elements). (=>) Let Q be irreducible and F := F_2[Z]/(Q), a field with 2^k elements. For every a in F, a^{2^k} = a (a = 0, or a^{2^k - 1} = 1 by Lagrange's theorem in the multiplicative group of order 2^k - 1). Applied to the class of Z this gives Z^{2^k} ≡ Z (mod Q). Z^2 + Z = Z (Z + 1) is the product of the two monic polynomials of degree 1; Q is irreducible of degree k >= 2, so it is coprime to both. (<=) Let g be a monic irreducible factor of Q, of degree h >= 1, and F := F_2[Z]/(g), a field with 2^h elements, with z := class of Z. From Z^{2^k} ≡ Z (mod Q) and g | Q we get z^{2^k} = z in F. Every element of F is a polynomial in z with F_2-coefficients, and a -> a^{2^k} is a ring homomorphism of F fixing F_2 (characteristic 2: (a + b)^2 = a^2 + b^2) and fixing z, hence it fixes every element: a^{2^k} = a for all a in F. Also a^{2^h} = a for all a in F (Lagrange as above). Let M := { j >= 1 : a^{2^j} = a for all a in F }. M contains h and k, is closed under addition (a^{2^{j+j'}} = (a^{2^j})^{2^{j'}}) and, for j > j' in M, contains j - j' (a = a^{2^j} = (a^{2^{j'}})^{2^{j-j'}} = a^{2^{j-j'}}). By the Euclidean algorithm gcd(h, k) in M, so all 2^h elements of F are roots of the polynomial X^{2^{gcd(h,k)}} - X of degree 2^{gcd(h,k)}; hence 2^h <= 2^{gcd(h,k)}, i.e. h = gcd(h, k) divides k. As k is prime, h in {1, k}. If h = 1 then g in {Z, Z + 1} divides gcd(Q, Z^2 + Z), excluded. So every irreducible factor of Q has degree k = deg Q, and Q (monic) is irreducible.

**Lemma 2.2 (valid prefixes exist).** For every prime k there is a monic irreducible polynomial of degree k over F_2.
Proof (standard finite-field facts). Let E be a splitting field of X^{2^k} - X over F_2. This polynomial has derivative -1, so it has 2^k distinct roots in E; the roots are closed under addition, multiplication and inversion of nonzero elements (a -> a^{2^k} is a ring homomorphism), so they form a field F with 2^k elements. Pick a in F \ F_2 (possible as 2^k > 2). The subfield F_2(a) of F has 2^h elements for h := deg of the minimal polynomial of a over F_2, and F is a vector space over F_2(a), so 2^k is a power of 2^h and h | k. Since a is not in F_2, h > 1, so h = k, and the minimal polynomial of a is monic irreducible of degree k.

**Evaluator E.** On input x of length n:
1. If n = 0 output 0. Compute k by testing k = 3, 4, 5, ... for primality by trial division and for k^5 > 32 n^3 in integer arithmetic; by Bertrand's postulate k <= 4 n^{3/5} + 2, so this takes poly(n) time. If n <= k output 0.
2. Form Q from the first k bits. Compute Z^{2^k} mod Q by k successive squarings modulo Q and compare with Z; compute gcd(Q, Z^2 + Z) by the Euclidean algorithm in F_2[Z]. If the test fails output 0. (poly(k) bit operations.)
3. Arithmetic in K = F_2[Z]/(Q) (addition, multiplication mod Q, squaring) costs poly(k) bit operations. Form y_1..y_r by padding and splitting the payload. For each i compute y_i^{alpha_i} with alpha_i = 2^{r+2-i} - 1 by the recursion w^{2^{j}-1} = (w^{2^{j-1}-1})^2 · w (j - 1 squarings and multiplications, j <= r + 1).
4. Compute Sym_r^t by the recurrence E_0 := 1, E_j := 0 for 1 <= j <= t; for i = 1..r: for j = min(i, t) down to 1: E_j := E_j + y_i^{alpha_i} E_{j-1}. (Splitting the t-subsets of [1, i] by whether they contain i shows that after processing i, E_j is the j-th elementary symmetric polynomial of the first i values.) O(rt) field operations. Output Tr(E_t), computed with k - 1 squarings and additions.
All quantities are bounded by polynomials in n (k, r, t <= n + 3), so E is one deterministic polynomial-time algorithm computing f_n for every n, with no advice. The language L := {x : f_{|x|}(x) = 1} is in P, hence in NP: the verifier V(x, w) accepts iff w is the empty string and E(x) = 1, witness bound W(n) = 0; both directions of "x in L iff exists w, |w| <= 0, V(x, w) = 1" are immediate. The circuits in the theorem read x (all n bits, prefix included) only; no witness bits are inputs; the prefix bits are counted in n.

## 3. The pinned model and the theorem

The circuit model is the one pinned by the rules (Li–Yang full B2: fan-in two, all 16 binary truth tables, each retained gate counted once including constant and unary gates, repeated pins allowed, unrestricted fan-out and depth, constants free, single output which may be a gate or an input node, exact computation on every input). `elimination.md` E1 restates this model, embeds every such circuit as an initial state of the elimination process with the same gate count, and explains why the generalized states only make the class of objects larger, so the bound proved for states applies to every pinned circuit. Cycles: no cyclic circuit is ever formed; the linear rewiring lives in an affine system with an invertible auxiliary matrix, maintained explicitly (E2–E3). Constant, unary and terminal-output cases: deletion D1 and E2.2.

**Restriction to a fixed prefix.** Let n > k, let b in F_2^k encode a valid Q, and let C be a pinned circuit on inputs x_1..x_n computing f_n. Replace the input nodes x_1..x_k by the constant nodes b_0..b_{k-1}. The result C' is a pinned circuit on the n - k = m payload inputs with exactly the gates of C (gates now reading constants are still counted), and it computes g_{n,Q} = f_n(b, ·). So |C| = |C'| >= C_B2(g_{n,Q}). (Only the gate count is transferred; no simplification is claimed or needed.)

**Theorem (main).** For every integer n >= 10^7 and every circuit C in the pinned model computing f_n,
    |C| >= (7/2)(n - k) - 7 (floor(8 n^{4/5}) + 1) >= (7/2) n - 56 n^{4/5} - (35/2) n^{3/5} - 7.
Consequently liminf_{n -> infinity} C_B2(f_n) / n >= 7/2. Explicitly, for every epsilon in (0, 7/2) put N_epsilon := max(10^7, (168/epsilon)^5, 21/epsilon); for every integer n >= N_epsilon we have 56 n^{4/5} <= (epsilon/3) n (as n^{1/5} >= 168/epsilon), (35/2) n^{3/5} <= (epsilon/3) n (as n^{2/5} >= (168/epsilon)^2 = 28224/epsilon^2 >= 52.5/epsilon for epsilon < 7/2) and 7 <= (epsilon/3) n, hence C_B2(f_n) >= (7/2 - epsilon) n; for epsilon >= 7/2 the inequality is trivial.

**Proof.** Let n >= 10^7, k = k(n) <= 5 n^{3/5} (Bertrand, `disperser.md` D7.1), m = n - k and d := d(n) = floor(8 n^{4/5}) + 1 < m (D7.1(iii)). By Lemma 2.2 choose a monic irreducible Q of degree k; its coefficient vector b is a valid prefix (Lemma 2.1). By Corollary D7.2 of `disperser.md`, g_{n,Q} : F_2^m -> F_2 is nonconstant on every affine subspace of F_2^m of dimension d. By Theorem E of `elimination.md` (with n' = m and this d), every pinned circuit computing g_{n,Q} has at least 3.5 m - 7 d gates. By the restriction paragraph, |C| >= 3.5 m - 7 d = 3.5 (n - k) - 7 d >= 3.5 n - 17.5 n^{3/5} - 56 n^{4/5} - 7. QED.

**Remainder and coverage.** The bound holds at every integer length n >= 10^7 (not a subsequence); the remainder 56 n^{4/5} + 17.5 n^{3/5} + 7 is o(n); the threshold 10^7 comes only from the parameter arithmetic D7.1 and could be lowered, but no threshold is required. All k prefix bits are charged to n: the bound is stated against the full input length n. The declared guarantee is g(n) = (7/2) n with liminf C_B2(f_n)/g(n) >= 1, and the finite projection (7/2, 1); the declared coefficient is exactly what is proved. Note that f_n itself is not claimed to be an affine disperser (it is 0 on all inputs with an invalid prefix); only its restriction to one fixed valid prefix is, and that restriction is a restriction of the actual input, which is why the circuit bound transfers.

## 4. Where each obligation is discharged

| Obligation | Where |
|---|---|
| Family defined at every length (n = 0, n <= k, invalid prefix, valid prefix), one algorithm, no advice | Sections 1–2 |
| P evaluator; NP verifier with empty witness; both directions | Section 2 |
| Field at every length: modulus from the input, exact prime-degree test, existence of valid prefixes | Section 2 (Lemmas 2.1, 2.2) |
| Disperser property of the restriction, complete proof | `disperser.md` D0–D7 (Theorem D, Corollary D7.2) |
| Parameter arithmetic, every n >= 10^7, d(n) < payload length | `disperser.md` D7.1 |
| Gate elimination in the exact pinned model | `elimination.md` E1–E7 (Theorem E) |
| Degenerate gates, constants, repeated pins, terminal output | `elimination.md` E1, D1, E2.2 |
| Cycles / rewiring | `elimination.md` E2–E3 |
| Every case stated and checked | `elimination.md` E3–E6 (D1–D4, AR, QR, Cases 1–10) |
| Prefix bits charged to n; restriction transfer | Section 3 |
| o(n) remainder, explicit | Section 3 |

## 5. Sources and credit

- Input-prefix modulus construction from Yukon submission f4f2c464 (the modulus Q is read from the first k input bits, tested by Z^{2^k} = Z mod Q and gcd(Q, Z^2 + Z) = 1, invalid prefixes output 0, and the lower bound fixes one valid prefix by restriction, charging the k prefix bits to n). The proofs of Lemmas 2.1–2.2 above are our own.
- O. S. Mishra and R. Santhanam, "An Explicit 3.5n − o(n) Size Lower Bound against General Boolean Circuits", ECCC TR26-260, 2026, Sections 2–4, Theorem 4.1 (|C| >= 3.5n - 7d for affine dispersers for dimension d). Part E re-derives it; the potential W = 2s + 3q + sum wt and the ten cases are theirs.
- E. Ben-Sasson and S. Kopparty, "Affine dispersers from subspace polynomials", author preprint 2010-03-12 (SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf), SIAM J. Comput. 41(4), 2012: Theorem 2.2, Sections 3, 6, 7. Part D re-derives Theorem 7.3 (hence 2.2) for p = 2 and prime k.
- M. G. Find, A. Golovnev, E. A. Hirsch, A. S. Kulikov, ECCC TR15-166 rev. 2 (fair semi-circuits, affine substitutions), and J. Li, T. Yang, ECCC TR21-023 (the pinned model conventions, Sections 1.1, 2.1, 2.4, 3.3): background only; no step is taken from them.
- E. Demenkov and A. S. Kulikov, MFCS 2011 (ECCC TR11-026): origin of the affine-disperser method; not used.

## 6. Honesty statement

Re-derived in full and checked line by line: all of Part E (every deletion, reparametrisation and case, with the potential accounting), all of Part D (subspace polynomials, structure lemma, permanent formula, independent-source theorem with its two claims, echelon decomposition, the two coefficient lemmas and the combinatorial claim, and the parameter arithmetic), and Lemmas 2.1–2.2. Nothing is transcribed by citation. Known deviations from the sources: Part D specialises BK to p = 2, linear block subspaces and prime k, which removes BK's "no affine shift of a proper subfield" hypothesis (D2.4(2)); Part E restates MS without change of parameters; the prefix wrapper is taken from submission f4f2c464 with our own proofs. Doubts: the threshold 10^7 is a convenient, not optimal, value. The whole function f_n is not an affine disperser; the theorem is proved through the restriction to a fixed valid prefix, as stated.
