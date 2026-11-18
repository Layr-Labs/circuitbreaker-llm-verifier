# A strict cubic cutwidth improvement and an eventual full-B2 circuit bound

For the single total language L in P given in Section 1, there exists an integer N_0>=256 such that for every integer length N>=N_0,

    C_B2(L_N) >= c_* N,    c_* = 482709141518790465583497877634923/100000000000000000000000000000000.

The circuit model counts every gate of the full binary basis and permits arbitrary sharing, as specified below. The declared remainder is identically zero. The proof first gives a uniform sublinear error at the strictly larger coefficient c_0=1+1/(2 gamma), where gamma=0.130647519423167665.

The new graph argument changes a Gaussian ordering when a seven-vertex sign witness permits a crossing matching to shrink. A quantitative covariance interpolation handles vertices where identifications in the base graph make this witness scarce. The asymptotic worst-case cutwidth ratio is at most Gamma-0.000000047, where Gamma is the limiting alternating-path coefficient of the earlier global cover method. The displayed rational gamma uses an upward rounding of Gamma, certified by the order-46 and order-47 path sums. The two complementary estimates are proved in certificate.md Sections 8 and 9. Section 7.6 of cutwidth.md combines them using nine-edge patches and two passes before the uniform multigraph reduction.

The immediate precursor is our quantitative retiming submission 60c7a1c3-9b1a-4c29-ad98-5de7c3ef351b, which supplies the two finite probability estimates. The publication and exact proof revision are https://github.com/Puthranps/cubic-cutwidth-retiming/tree/a901153f78d426225407c02dec82995e51192f6d. The first retiming submission is 4c133ffa-7452-42c6-9629-d7821ee890d4, whose hosted reviewer marked its 10^(-14) graph saving and eventual circuit coefficient 48270900384055950810649810365841/10000000000000000000000000000000 valid. Here the needed components are supplied in full and the new finite checks execute afresh. The preceding source is DarkWingMcQuack's submission 3f2be7ea-aa22-4f7e-98e9-23b40557dbcf at commit 737818f9911e00e9d4909c644c771e0724fb6c37, which supplies the quantified declaration format and the coefficient 630647566423167697/130647566423167697. Its immediate numeric predecessor 0514f3af-5d65-4556-82df-a89e2a3ae4e1 established the finite order-39 estimate. ac0dd6e1 and 38e961ae supplied fixed-point contraction and factorization respectively. The original global crossing-cover and finite path construction is due to galarnon42, submission fcf21648-1b1f-463c-b5bc-c50f53f63436. The mathematical route from circuits to layouts is due to McGuire, Schlesinger and Williams. Arthur Ostapenko supplied the retained norm-one family, Hellinger comparison and algebraic transfer, crediting mooselumph for the validated-prefix device and the earlier submission fd902316-f679-4965-aca9-419a9ee5ec22.

The single language submitted to CircuitBreaker is the validated-prefix norm-graph language described below. The original McGuire–Schlesinger–Williams paper uses Li's sumset extractor. The universal graph inequality can also be inserted in their circuit-to-kernel theorem. These are separate applications.

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
* **Every length, invalid prefixes.** A is defined on all strings and the lower-bound transfer holds at every integer N >= 256.
  Every circuit must compute L_N on all inputs, including invalid prefixes; the lower bound uses only the
  restriction of the compiled program to one valid prefix pair (`circuit.md` §2.4, Corollary 4.5), which exists
  at every N >= 16 (`circuit.md` Lemma 2.3).

## 4. From layouts to circuits

The finite computations in certificate.md certify the good-center and tilted edge-credit sums, compare Gamma_47 with 0.130647566423167665 and verify the rational sum for the pivotal event. Section 7 of cutwidth.md proves that for each fixed odd J the asymptotic cubic coefficient is at most Gamma_J-0.000000047. Theorem 1 of that file extends the bound to loops, parallel edges, bridges and disconnected components and supplies a uniform nonnegative nondecreasing h_0 with h_0(D)=o(D). Thus the graph hypothesis H(gamma_0,h_0) of circuit.md holds with gamma_0=26129503884633533/200000000000000000. The algebraic norm-graph rectangle bound is proved in krs.md. The sweep simulation, wire-graph kernel conversion and all-length transfer are proved in circuit.md.

## 5. Assembling the eventual theorem

Theorem 1 of circuit.md, applied to this h_0 and the fixed evaluator in Section 1,
gives for EVERY integer N>=256

    C_B2(L_N) >= c_0 N-R_0(N),
    c_0 = 1+1/(2 gamma_0),       s_0 = ceil(c_0 N),
    R_0(N) = 5sqrt(N)+(1/2)sqrt(N)log2(N)+2
       +[7sqrt(N)+sqrt(N)log2(N)+6+3log2(2s_0)+h_0(2s_0)]/(2 gamma_0).

Here gamma_0, c_0 and h_0 are fixed. Since s_0=O(N) and tends to infinity,

    h_0(2s_0)/N = [h_0(2s_0)/(2s_0)] * [2s_0/N] -> 0.

Every other displayed summand of R_0(N)/N also tends to zero. Thus R_0(N)/N->0
along ALL integer lengths, with the same R_0 for every circuit computing L_N.

Here c_0=126129503884633533/26129503884633533. Choose c_* = floor(10^32 c_0)/10^32 = 482709141518790465583497877634923/100000000000000000000000000000000. Their difference is the strictly positive rational 12576547382327041/2612950388463353300000000000000000000000000000000. Since R_0(N)/N tends to zero along all integer lengths, there exists N_0 such that R_0(N)/N <= c_0-c_* for every N>=N_0. Substitution gives C_B2(L_N)>=c_* N for every such length. This also proves the finite projection with coefficient c_* and exponent 1. The threshold is existential. All witnesses, prefixes and leftover bits are handled by the fixed definition above.

## 6. Sources and reproduction

- Immediate submitted source: https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/737818f9911e00e9d4909c644c771e0724fb6c37/submission

- Earlier dependency: submission 0514f3af-5d65-4556-82df-a89e2a3ae4e1,
  https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/d18ef5c1b86c58ca9f0bdb93085f9b6f2698c55d/submission

- Exact-integer dependency: submission ac0dd6e1-8a73-4682-b7aa-67502b2de08b,
  https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/77f07bfd2fd633cf98e3737575ba188d20e95c6a/submission
- Earlier rejected variant: submission 38e961ae-bd55-4566-85ad-17112333f84a,
  source commit f6d3da95b1a680e9bb986008bd6d0890081cdaa7; attribution only.

- Direct dependency: submission fcf21648-1b1f-463c-b5bc-c50f53f63436,
  https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/81470c9513639af985c88699dbcd2e2efbef7463/submission

- McGuire--Schlesinger--Williams: https://github.com/rrwilliams0/circuit-lb-4.6155n,
  especially the circuit-to-kernel reduction and median-edge ordering.
- CircuitBreaker source for the retained family/circuit/algebra exposition:
  https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/5fef109bf7f1b664295244a0dee6df066ed2232d/submission
- Gaspers, "Path Decompositions of Sparse Graphs via Gaussian Sweeps and Faster Exponential Algorithms", https://arxiv.org/abs/2611.05146.
- Kollar--Ronyai--Szabo, "Norm-graphs and bipartite Turan numbers", Combinatorica
  16(3) (1996), 399--406; the needed theorem is proved in `krs.md`.

The first retiming version is publicly archived at https://github.com/Puthranps/cubic-cutwidth-retiming/tree/7b0d2587b3da62877fea4fae2d4fa9a370f3919f.

From the submission directory, reproduce the numerical certificate with

    python3 -I -B certify.py < /dev/null

A successful run exits zero and prints one JSON object with `certified: true`.
