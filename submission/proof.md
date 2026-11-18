# An eventual circuit lower bound with coefficient 315323783211583832170874691/65323783211583832170874691

## Claim and scope

For the single fixed total language L in P defined in Section 1, there is an
integer N_0 >= 256 such that EVERY integer N >= N_0 satisfies

    C_B2(L_N) >= (315323783211583832170874691/65323783211583832170874691) N.

This is an eventual lower bound with identically zero subtracted remainder,
in the version-2 declaration schema. It also implies the declared finite
asymptotic projection with that coefficient and exponent 1. The full B2 model
permits all 16 binary gate types, repeated pins, constants and arbitrary
sharing. All prefix and leftover bits count toward the actual length N.
The evaluator and empty-witness NP verifier are fixed independently of N_0
and every graph approximation parameter. No explicit value of N_0 is claimed.

The proof first establishes the strictly stronger intermediate coefficient

    gamma_0 = 130647566423167664341749381/1000000000000000000000000000
            = 0.130647566423167664341749381,
    c_0 = 1+1/(2 gamma_0) = 630647566423167664341749381/130647566423167664341749381,

with a single uniform sublinear loss R_0(N). A fixed positive coefficient
margin then absorbs this loss for all sufficiently large integer lengths.
Section 5 proves that quantifier conversion explicitly.

## Provenance and contribution of this revision

The immediate source is CircuitBreaker submission
3f2be7ea-aa22-4f7e-98e9-23b40557dbcf at commit
737818f9911e00e9d4909c644c771e0724fb6c37. Its complete fixed-language,
full-B2, algebraic, graph and positive-slack eventual-bound proofs are retained
for review anew. This revision replaces its numerical estimates with the full
localized quadrature and triple-tail proofs in certificate.md, certifies 65
finite path terms, and specializes the same general graph theorem to the
smaller gamma_0 in this proof. No prior judgment is a mathematical premise.

The directly used earlier submissions are also declared for attribution:
0514f3af-5d65-4556-82df-a89e2a3ae4e1 (commit
d18ef5c1b86c58ca9f0bdb93085f9b6f2698c55d) supplied the Gaussian-moment
and positive-completion refinements; ac0dd6e1-8a73-4682-b7aa-67502b2de08b
(commit 77f07bfd2fd633cf98e3737575ba188d20e95c6a) supplied exact integer
contractions; 38e961ae-bd55-4566-85ad-17112333f84a (commit
f6d3da95b1a680e9bb986008bd6d0890081cdaa7) supplied the factored kernel,
with its original hosted Decimal timeout disclosed.

The original direct mathematical package is submission
fcf21648-1b1f-463c-b5bc-c50f53f63436, by galarnon42, at commit
81470c9513639af985c88699dbcd2e2efbef7463, also a declared dependency.
The circuit-to-cutwidth route is due to McGuire, Schlesinger and Williams.
Arthur Ostapenko's snapshot 5fef109bf7f1b664295244a0dee6df066ed2232d supplies
the retained norm-one family, circuit transfer, Kollar--Ronyai--Szabo exposition
and Hellinger comparison, and credits mooselumph's f4f2c464 submission for the
validated-prefix device. The original direct source also credits submission
fd902316-f679-4965-aca9-419a9ee5ec22 (source aeabecef2df58965408d952605327debee2476e4).
Matching-cut submodularity and canonical covers have earlier literature,
including Jeong--Saether--Telle, arXiv:1507.02384; related Gaussian sweeps appear
in Gaspers, arXiv:2611.05146. All needed nonstandard arguments are supplied.

The new quadrature proof retains unaffected conditional-density factors when
differentiating one coordinate. Integrating or quadraturing them costs at most
one half per step, instead of replacing them all by volume factors. Reversible
three-coordinate anchors strengthen the outlier tails; two-step positive
completions supply further decay. The first path probability is evaluated
exactly by an arctangent enclosure. The rescaled exact tensor state
4^(j-1) M_j prevents loss of resolution on the fixed 2^128 grid. Certificate.md
proves every one of these changes and the interval code implements them.

The production parameters are N=80, J=65, L=3, and the fixed graph constant is

    C_65 = 1+3*65*2^64 = 3597115094373362565121.

The language, circuit model, exponent and eventual quantifier are unchanged.
The exact increase over the immediate source is

    315323783211583832170874691/65323783211583832170874691
      - 630647566423167697/130647566423167697
      = 8164562654500000000000000000/8534393306148005579050468452335659466056627
      > 0.

The increase is approximately 9.56665853285e-16. This is a small numerical
improvement within the same method. No global optimality or infinite-series
evaluation is claimed as part of this theorem.

## Evidence and execution status

The only finite numerical premise is

    Gamma_65 <= 0.1306475664231676643417493807
             < 0.130647566423167664341749381
             = 130647566423167664341749381/1000000000000000000000000000 = gamma_0.

certificate.md proves the path density, finite sum, root and weight validation,
localized Gaussian-moment quadrature error, reversible triple-outlier tails,
sharp positive-completion factors, and the rescaled exact integer enclosure.
certify.py computes all 65 terms using only the Python standard library. Its sole declared input is {}, and its production defaults
are N=80, J=65, L=3, 100-digit Decimal initialization and 2^128-scale exact
integer contractions rescaled by 4^(j-1). It exits zero only after a strict
final comparison with gamma_0. This finite arithmetic premise is supplemented by the full graph,
algebraic, circuit and uniformity proofs; computation alone is insufficient.

The complete theorem is assembled below from cutwidth.md, circuit.md, krs.md
and certificate.md. Both scored claims contain the identical complete version-2
declaration and the same seven evidence files. The local verifier is pinned
to the protected workflow revision fa3c2044e179e919bfb09b890594f21458ef142a.
Docker is unavailable locally; completed host executions are distinct from
the required fresh hosted experiment. No prior judgment is reused.

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

## 4. Proof dependencies

All nonstandard graph, Gaussian, layout and numerical arguments used here are
proved in the supplied files. The retained algebraic proof states its textbook
commutative-algebra premises. The graph proof states explicitly the standard
finite bipartite matching/cover theorem used in its cover-lifting and Hall steps.

1. `krs.md` proves the finite-field product-system bound for the rectangle argument.
2. `circuit.md` proves P membership, the rectangle lower bound, wire-graph
   simulation and the transfer from graph hypothesis H(gamma,h).
3. `cutwidth.md` proves H(gamma,h) with gamma=130647566423167664341749381/1000000000000000000000000000 from the sole finite
   numerical inequality Gamma_65<130647566423167664341749381/1000000000000000000000000000.
4. `certificate.md` and `certify.py` establish that finite inequality.

The graph theorem applies to every graph, not just graphs of large girth.
Finite covers are used only for the expectation comparison, and the actual
finite-range random ordering lives on the original graph. The field radius,
noise and threshold mesh are fixed before graph size tends to infinity. The
resulting remainder is uniform over graphs.

## 5. Assembling the eventual circuit theorem

### 5.1 A strictly stronger intermediate coefficient

Put gamma_0=130647566423167664341749381/1000000000000000000000000000. The numerical inequality in
certificate.md is strict below this constant. Theorem 1 of cutwidth.md is
explicitly instantiated at gamma_0: it supplies one nonnegative nondecreasing
h_0 with h_0(D)/D->0 such that every cubic multigraph on D vertices has cutwidth
at most gamma_0 D+h_0(D). Therefore H(gamma_0,h_0) of circuit.md holds, including
its uniformity over all connected cubic multigraphs.

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
Exact arithmetic gives

    c_0 = 1+1000000000000000000000000000/261295132846335328683498762
        = 1+500000000000000000000000000/130647566423167664341749381
        = 630647566423167664341749381/130647566423167664341749381.

### 5.2 Absorbing the sublinear loss with fixed positive slack

Define the declared coefficient

    c_* = 315323783211583832170874691/65323783211583832170874691.

The exact difference from the independently proved intermediate coefficient is

    delta = c_0-c_*
          = 250000000000000000000000000/8534393306148003445689984947206164501196135877816271 > 0.

Both fractions are reduced: their numerator-denominator gcds are 1.
Since R_0(N)/N->0 and delta/2 is a fixed positive number, the definition of
convergence provides an integer N_0>=256 such that EVERY integer N>=N_0 satisfies

    R_0(N)/N <= delta/2.

For all those lengths the bound in Section 5.1 yields

    C_B2(L_N) >= c_0 N-R_0(N)
              >= (c_0-delta/2)N
               = (c_*+delta/2)N
              >= c_* N.

This proves precisely the version-2 eventual_lower_bound with leading term
c_* N and explicit remainder terms []. The threshold N_0 may be large and is
not evaluated; its existence follows from the single uniform convergence
statement just proved. The proof uses a strictly larger c_0 before absorbing
the loss, so it does not infer a zero remainder from a same-coefficient liminf.
No rate of convergence for h_0 or R_0 is assumed.

The family is unchanged throughout: every parameter is independent of N_0,
and the assertion holds at all integer N>=N_0, not a sparse subsequence.
Dividing the eventual inequality by N also proves

    liminf over all integer N tending to infinity C_B2(L_N)/N >= c_*.

Thus the finite projection has coefficient c_* and exponent 1. The evaluator
and witness bound in Section 2 prove NP eligibility for this same language.
Together these statements prove the entire joint declaration in bound.json,
including both score coordinates and the zero explicit remainder. QED.

## 6. Sources and reproduction

- Immediate dependency: submission 3f2be7ea-aa22-4f7e-98e9-23b40557dbcf,
  https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/737818f9911e00e9d4909c644c771e0724fb6c37/submission

- Earlier numerical dependency: submission 0514f3af-5d65-4556-82df-a89e2a3ae4e1,
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

From the submission directory, reproduce the numerical certificate with

    python3 -I -B certify.py < /dev/null

A successful run exits zero and prints one JSON object with `certified: true`.
