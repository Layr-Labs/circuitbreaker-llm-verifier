# Coefficient 5703/1250, exponent 1: a validated-field family and a cutwidth argument

## Claim

Let L be the single total language defined in §1, A its evaluator and L_N its slice at length N, with
**N denoting the entire input length**.

**Main Theorem** (the `theorem` of `bound.json`). A runs in polynomial time on every input, so L is in P; L is in
NP with witness bound W(N) = 0 and the fixed verifier V that accepts (x, w) iff w is empty and A(x) = 1. In the
full pinned B2 model, with N the entire input length: for every real ε > 0 there is N_ε such that
C_B2(L_N) >= (5703/1250 − ε)·N for every integer N >= N_ε; equivalently, liminf over all integers N of
C_B2(L_N)/N >= 5703/1250.

**Route.** Put

    x0 = (1 + 2√2)/4,   θ0 = arccos x0,   γ = (3/(2π))·θ0,   c* = 1 + 1/(2γ) = 1 + π/(3θ0).

We prove that for every integer N >= 256

    C_B2(L_N) >= c*·N − R(N),                                                         (★)

with the explicit remainder R(N) of §6.1, that R(N) <= 1000·N/log2 N for every N >= 2^16 (§6.3), and that
c* > 5703/1250 (§6.2); then N_ε = 2^max(16, ⌈1000/ε⌉) works (§6.4). The declared exact pair is
**c = 5703/1250 (= 4.5624), p = 1**, both coordinates of this one theorem for this one family.

**Attribution.**
* **Method.** The route circuit -> wire graph -> cutwidth ordering -> read-once nondeterministic program ->
  rectangle edge counting, and the Gaussian edge-score ordering of cubic graphs that gives the constant γ (so
  their coefficient is the same c*), are from S. McGuire, S. Schlesinger, R. Williams, "Explicit functions
  requiring larger than 4.5n size circuits", preprint, November 2026; the authors report a Lean 4 formalization
  (§8). Their hard function is a sumset extractor of Li. Their paper and code are credited, not used as
  premises: every step is re-proved in `circuit.md` and `cutwidth.md`.
* **Norm graphs.** The rectangle property of the hard function rests on a theorem of J. Kollár, L. Rónyai and
  T. Szabó (premise P1), proved in `krs.md` along the lines of their proof.
* **Validated prefix.** Making the field-defining polynomials part of the input, checked by an exact test, is
  taken from an earlier CircuitBreaker package (coefficient 4, certificate-prefix family), which credits it to
  mooselumph's package f4f2c464.
* **This package** is an independent write-up with a different explicit function: the indicator of z^e = 1,
  e = 1 + q + ... + q^{t−1} with q = 2^m, for the payload element z of a field with q^t elements whose two
  defining polynomials are validated input prefixes (the condition defining the norm graphs of Kollár, Rónyai
  and Szabó). The cutwidth bound is proved for connected cubic multigraphs with loops and parallel edges, with
  an explicit remainder, so the final remainder is explicit.

No verdict, review or acceptance of any other package is a premise. There are no declared experiments and no
code files; every numerical fact used is proved with rational witnesses (§6.2–§6.3, `cutwidth.md` §5).

## 1. The family: one evaluator on every bit string

### 1.1 Length parameters and input layout

For an integer N >= 16 define

    k  = the largest integer with 4^k <= N,      t = 2^k,
    m0 = ⌊(N − t)/(t + 1)⌋,   m = m0 if m0 is odd, m = m0 − 1 if m0 is even,
    ℓ  = N − t − m − t·m.

Thus t is a power of two, m is odd and gcd(m, t) = 1. An input x = x_1 … x_N is split, in this order, into the
**R-prefix** x_1..x_t, the **Q-prefix** x_{t+1}..x_{t+m}, the **payload** of t·m bits c_{u,v} (0 <= u < m,
0 <= v < t; c_{u,v} is bit number t + m + v·m + u + 1), and ℓ **leftover** bits. So N = t + m + t·m + ℓ exactly.

### 1.2 The evaluator A

On x ∈ {0,1}^N:

1. If N < 16, output 0.
2. Compute k, t, m, ℓ as above. If m < 1, output 0.
3. Let R(W) = W^t + Σ_{v<t} x_{v+1} W^v ∈ F_2[W] and Q(Z) = Z^m + Σ_{u<m} x_{t+u+1} Z^u ∈ F_2[Z].
   If IrrTest(R) or IrrTest(Q) rejects, output 0.
4. Let F = F_2[Z]/(Q) and E = F[W]/(R). Form z = Σ_{v<t} (Σ_{u<m} c_{u,v} Z^u) W^v ∈ E.
5. Let e = Σ_{j=0}^{t−1} 2^{m·j}. Compute z^e in E by square-and-multiply. Output 1 iff z^e = 1.

The leftover bits are never read. **IrrTest(P)**, for a monic P ∈ F_2[X] of degree d >= 1, accepts iff
X^{2^d} ≡ X (mod P) and, for every divisor d' of d with 1 <= d' < d, gcd(P, (X^{2^{d'}} − X) mod P) = 1, where
X^{2^j} mod P is computed by j successive squarings modulo P.

L = {x : A(x) = 1} and L_N = L ∩ {0,1}^N. This defines L_N at **every** length, including N < 16, the guard
m < 1 (which never fires for N >= 16, `circuit.md` Lemma 6.1(a)) and every invalid prefix. The polynomials R and
Q are input bits checked by the evaluator; they are not advice and not the output of a construction algorithm.
When both prefixes pass, E is a field with 2^{tm} elements (`circuit.md` Lemmas 2.2 and 2.4).

## 2. Uniform P and NP with empty witnesses

A runs in O(N^4) bit operations on every input (`circuit.md` §2.2), so L ∈ P. Let V(x, w) accept exactly when w
is the empty string and A(x) = 1, with witness bound W(N) = 0. If x ∈ L, the empty witness is accepted. If
V(x, w) = 1 for some permitted w, then A(x) = 1 and x ∈ L. This proves both directions with one fixed
polynomial-time verifier for all lengths. Witnesses are not circuit inputs.

## 3. Circuit model and coverage

The bound is for the pinned full B2 model (Li–Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3):
deterministic acyclic circuits with one output, two ordered input pins per gate, any of the 16 binary truth
tables per gate (constants, projections and negated projections included), unrestricted fan-out and depth,
repeated pins and shared subcircuits allowed, exact computation on every input of length N, nonuniform circuits.
Every gate is counted once; inputs, wires and output selection cost nothing. The formal model is `circuit.md` §1.

* **All gate types, repeated pins, sharing.** The only gate-dependent step is the local rule
  β(e_g) = φ_g(β(p_1), β(p_2)) of `circuit.md` Theorem 4.4, valid for every φ_g; degenerate gates are counted like
  all others. A gate with both pins fed by one node has two distinct pin edges (`circuit.md` Definition 4.1),
  counted twice in Lemma 4.2; their equal values are forced by the copy rules, not assumed.
* **No substitution into circuits.** The proof never fixes inputs of a circuit. It compiles the whole circuit
  (with its R-, Q-prefix and leftover inputs) into a read-once program and restricts the *program*
  (`circuit.md` Lemma 3.1, Corollary 4.5). Hence constant nodes created by substitutions never arise.
* **Constant nodes.** If the model admits free constant nodes at all, `circuit.md` Lemma 1.1(i) removes them by
  absorbing each constant into the truth table of the gate it feeds; a gate with two constant pins becomes a
  constant-type gate reading x_1 twice. No gate is added, no gate is removed except non-ancestors of the output
  (Lemma 1.1(ii)), and no uncharged basis element is introduced. If the model has no constant nodes this step is
  void.
* **Terminal cases.** L_N is neither constant nor a projection for N >= 256 (`circuit.md` Proposition 6.2(1)),
  so after Lemma 1.1 the output is a gate and every other gate has fan-out >= 1. Inputs of fan-out 0 get no
  vertex and are never read (Theorem 4.4); unused payload variables are counted by r (Proposition 6.2(2)). Wire
  graphs that are trees or unicyclic are handled in Proposition 6.2(4).
* **Every length, invalid prefixes.** A is defined on all strings and (★) holds at every integer N >= 256.
  Every circuit must compute L_N on all inputs, including invalid prefixes; the lower bound uses only the
  restriction of the compiled program to one valid prefix pair (`circuit.md` §2.4, Corollary 4.5), which exists
  at every N >= 16 (`circuit.md` Lemma 2.3).

## 4. Structure of the proof

* `circuit.md` **Theorem 1** (circuit side): if the graph hypothesis H(γ,h) of `circuit.md` §0 holds, then for
  every N >= 256, C_B2(L_N) >= (1 + 1/(2γ))N − R(N) with R as in §6.1. Chain: normal form (Lemma 1.1); the
  restricted payload function has property D(t, t! + 1) and e ones (Lemmas 2.1–2.7, with P1 in Lemma 2.6); edge
  lower bound for read-once programs (Lemmas 3.1–3.2, Theorem 3.3); wire graph and sweep program (Lemma 4.2,
  Theorem 4.4, Corollary 4.5); kernel (Lemmas 5.1–5.4, Proposition 5.5); parameters and the main inequality
  (Lemma 6.1, Proposition 6.2).
* `cutwidth.md` **Theorem 1** (graph side): H(γ,h) holds with γ = (3/(2π))·arccos((1 + 2√2)/4) and
  h(D) = 24D/log2(D+2). Chain: Gaussian facts (§1), reduction of cubic multigraphs to simple cubic graphs (§2),
  Theorem 2 for simple cubic graphs (§3), assembly (§4), numerical inequalities (§5).
* `krs.md`: premise P1 (Kollár–Rónyai–Szabó, Combinatorica 16(3) (1996) 399–406, Theorem 3.3; stated in
  `circuit.md` §0 and `krs.md` §0), proved for every field via Theorem A (`krs.md` §0, §3).
* §6 below combines them and does the arithmetic.

## 5. Facts used without proof

* Textbook commutative algebra F1–F6 (Atiyah–Macdonald; Milne) and elementary facts E1–E5: `krs.md` §2.
* Textbook algebra (S1)–(S6): `circuit.md` §0.
* Undergraduate analysis and probability: `cutwidth.md` §1.0; and, for the numerical witnesses of §6.2 and
  `cutwidth.md` §5, the Taylor series of cos, arctan and ln, the fact that partial sums of an alternating series
  with decreasing terms bracket its sum, the addition formulas for cos and tan, and monotonicity of cos on [0, π].

Nothing else is used without proof. The Li–Yang paper fixes the model conventions only. The sources credited in
§8 are not premises. No deterministic construction of irreducible polynomials, no prime number theorem and no
Bertrand postulate is used: the field-defining polynomials are input bits, and their existence (`circuit.md`
Lemma 2.3) is used only inside the lower-bound proof.

## 6. The main theorem and its derivation

### 6.1 Combining `circuit.md` and `cutwidth.md`

`cutwidth.md` Theorem 1 is the hypothesis H(γ,h) of `circuit.md` §0 with

    γ = (3/(2π))·arccos((1 + 2√2)/4),   h(D) = 24·D/log2(D + 2).

Both files use the same class (connected cubic multigraphs; loops and parallel edges allowed; a loop adds 2 to
the degree), the same cut convention (a loop never crosses a cut, parallel edges count separately) and the same
cutwidth (maximum over prefix cuts, minimum over orderings). `circuit.md` applies H(γ,h) only to the kernel
K of Lemma 5.2, which is a connected cubic multigraph in exactly this sense, on D = 2(s' − n_u) >= 2 vertices,
and evaluates h only at positive integers, where h >= 0 and h is nondecreasing. So `circuit.md` Theorem 1
applies and gives (★) for every integer N >= 256 with c* = 1 + 1/(2γ) and

    R(N) = 5√N + (1/2)√N·log2 N + 2
           + (c* − 1)·[ 7√N + √N·log2 N + 6 + 3·log2(2s_0) + h(2s_0) ],     s_0 = ⌈c*·N⌉.

(Here 1/(2γ) = c* − 1.) For N < 256 one may set R(N) = c*·N; then (★) holds trivially at every length.

### 6.2 The declared coefficient: c* > 5703/1250

The comparison c* > 5703/1250 is equivalent to 1/(2γ) > 4453/1250, i.e. to γ < 625/4453. Proof by hand:

* Let a = 0.293955. Since 0 < a < 1, the terms a^{2k}/(2k)! of the cosine series decrease, so the alternating
  partial sum ending with a positive term is an upper bound: cos a <= 1 − a²/2 + a⁴/24. Exactly
  a² = 0.086409542025 and a⁴ = 0.007466608952970241100625, so a⁴/24 < 0.000311110 and
  cos a < 1 − 0.0432047710125 + 0.000311110 = 0.9571063389875.
* 1.41421356² = 1.9999999932878736 < 2, so √2 > 1.41421356 and x0 = (1 + 2√2)/4 > 0.95710678.
* Hence cos a < x0 = cos θ0. Since cos is strictly decreasing on [0, π] and a, θ0 ∈ [0, π], θ0 < a.
* With π > 3.14159 (`cutwidth.md` (W1) gives 3.14159265358 < π): γ = 3θ0/(2π) < 3a/(2·3.14159), and
  3a·4453 = 3926.944845 < 3926.9875 = 625·2·3.14159. So γ < 625/4453 and c* > 1 + 4453/1250 = 5703/1250.
  The same numbers give c* > 1 + 3.14159/(3a) > 4.56243.

We also need **γ > 1/8**: cos(π/12) = (√6 + √2)/4, and √6 + √2 > 1 + 2√2 ⟺ √6 > 1 + √2 ⟺ 6 > 3 + 2√2 ⟺
3 > 2√2 ⟺ 9 > 8. So cos(π/12) > x0, hence θ0 > π/12 and γ > 3/(2·12) = 1/8. Consequently c* − 1 < 4 and c* < 5.

### 6.3 The remainder: R(N) <= 1000·N/log2 N for every N >= 2^16

Let N >= 2^16 and L = log2 N >= 16.

* L² <= √N = 2^{L/2}: equality at L = 16, and 2 ln L − (L/2) ln 2 decreases for L > 4/ln 2. Hence
  √N·L <= N/L and √N <= N/L² <= N/(16L).
* N/L is increasing, so N/L >= 2^16/16 = 4096 and any constant b satisfies b <= (b/4096)·N/L.
* c* < 5 gives s_0 = ⌈c*·N⌉ <= 5N, so 2s_0 <= 10N and 3·log2(2s_0) <= 3(L + 4). Moreover
  3(L + 4) <= 0.015·N/L, since (3L² + 12L)/2^L is 960/65536 < 0.015 at L = 16 and decreases for L >= 16 (the
  logarithmic derivative of 3L² + 12L is (6L + 12)/(3L² + 12L) <= 108/960 < ln 2).
* h is nondecreasing, so h(2s_0) <= h(10N) = 240N/log2(10N + 2) <= 240N/L.

Therefore the first line of R(N) is at most (5/16 + 1/2 + 2/4096)·N/L < 0.8131·N/L, the bracket is at most
(7/16 + 1 + 6/4096 + 0.015 + 240)·N/L < 241.454·N/L, and with c* − 1 < 4 (§6.2)

    R(N) < (0.8131 + 4·241.454)·N/L < 966.7·N/L < 1000·N/log2 N.

### 6.4 Conclusion

Let ε > 0 and N >= N_ε = 2^max(16, ⌈1000/ε⌉), any integer. Then log2 N >= 1000/ε, so R(N) <= 1000N/log2 N <= εN,
and by (★) and §6.2

    C_B2(L_N) >= (c* − ε)·N > (5703/1250 − ε)·N.

This holds at **every** integer length N >= N_ε, not on a subsequence, for the one fixed language L with its one
verifier and witness bound 0. This is the Main Theorem, with this N_ε. Since ε > 0 is arbitrary, moreover

    liminf over all integers N of C_B2(L_N)/N >= c* > 5703/1250,

which gives the declared polynomial_log guarantee with coefficient 5703/1250, exponent 1, log and loglog exponents
0, and its finite projection (5703/1250, 1).

## 7. Bit accounting

* N = t + m + t·m + ℓ exactly. The R- and Q-prefixes (t + m <= 3√N bits) and the leftover (ℓ <= 2√N + 1 bits)
  are real circuit inputs; they are charged twice in the exact remainder (`circuit.md` (6.3)): once in the
  linear term, since they are not payload inputs and do not contribute to n_u, and once in the bracket, since
  they do not contribute to log2 N1. The payload loses m bits of density (N1 = e, not 2^{tm}) and
  O(√N log N) bits for the rectangle parameters (r_max and log2((t − 1)·t!)). The other losses are
  (3·log2(2s_0) + 5 + h(2s_0))/(2γ), from the program size and the cutwidth remainder. All are o(N); no linear
  loss is hidden.
* No formula, monotone, bounded-depth, XOR-free, multi-output, randomized or average-case statement is used, so
  no transfer is needed. No direct-sum additivity is assumed and no potential is identified with a gate count:
  the bound counts the gates of an arbitrary circuit through the exact vertex and edge counts of its wire graph
  (`circuit.md` Lemma 4.2).

## 8. Sources

* S. McGuire, S. Schlesinger, R. Williams, "Explicit functions requiring larger than 4.5n size circuits",
  preprint, November 2026. The authors' Lean 4 formalization: github.com/SamuelSchlesinger/circuit-lb-4.56n-lean
  (commit b7c27d13c0545dedd78ec2614ca6792df160e8ed, 2026-10-13). Credit for the method; not a premise.
* J. Kollár, L. Rónyai, T. Szabó, "Norm-graphs and bipartite Turán numbers", Combinatorica 16(3) (1996)
  399–406, Theorem 3.3: premise P1, proved in `krs.md`.
* N. Alon, L. Rónyai, T. Szabó, "Norm-graphs: variations and applications", J. Combin. Theory Ser. B 76 (1999)
  280–290: P1 with the indices exchanged (`krs.md` Theorem A).
* M. F. Atiyah, I. G. Macdonald, Introduction to Commutative Algebra, Addison-Wesley, 1969; J. S. Milne, Fields
  and Galois Theory, version 5.10: the textbook facts F1–F6 of `krs.md` §2.
* Li–Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3: the pinned full-B2 conventions. Its
  lower-bound arguments are not used.
* An earlier CircuitBreaker package (coefficient 4, certificate-prefix family) and mooselumph's package
  f4f2c464, credited there: the validated-prefix device. Their family and arguments are not used.
* S. Gaspers, arXiv:2611.05146v1 (preprint, November 2026), Lemma 3.4: a bound of the same 1/R form as
  `cutwidth.md` Lemma 3.2; not used.
