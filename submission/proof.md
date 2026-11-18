# A structural improvement to the eventual full-B2 circuit lower bound

## Claim and scope

For the single fixed total language L in P defined in Section 1, there is an
integer N_0 >= 256 such that EVERY integer N >= N_0 satisfies

    C_B2(L_N) >= (315323783211583831670874691/65323783211583831670874691) N.

The exact coefficient is approximately 4.82709003840530217799497600319.
The version-2 declaration is an eventual lower bound with identically zero
subtracted remainder, and its finite projection has this same coefficient
and exponent 1. The full B2 model permits all sixteen binary truth tables,
repeated pins, constants and arbitrary sharing. All prefix and leftover bits
count toward the actual input length N. The evaluator and empty-witness
verifier are fixed independently of graph approximation parameters and N_0.
No explicit value of N_0 is claimed.

The argument keeps the previously certified finite path sum, and proves a
new structural improvement to the graph ordering. Distinguish

    gamma_base = 130647566423167664341749381/1000000000000000000000000000
               = 0.130647566423167664341749381,
    delta_g = 1/1000000000000000000,
    gamma_0 = gamma_base-delta_g = 130647566423167663341749381/1000000000000000000000000000
            = 0.130647566423167663341749381,
    c_0 = 1+1/(2 gamma_0) = 630647566423167663341749381/130647566423167663341749381.

The finite path inequality is Gamma_65 < gamma_base. It is NOT asserted to
be below gamma_0. Section 11 of cutwidth.md supplies the additional delta_g
through buffered local majority improvements and strict Gaussian comparison.
It proves H(gamma_0,h_0) with one common nondecreasing h_0=o(D). The circuit
transfer gives the strictly larger intermediate coefficient c_0 with one
uniform sublinear remainder. Section 5 absorbs that remainder using fixed
positive coefficient slack.

## Provenance and contribution of this revision

The immediate source is CircuitBreaker submission
ae372f38-09d4-452a-af05-4954cdd2f168, commit
2a8b03fe0be18fdf4918dc5fc34da51534f1aca1. Its complete fixed-language,
full-B2, algebraic, baseline graph and finite-path proofs are reproduced
in the supplied evidence for review anew. Its N=80, J=65, L=3 interval
certificate remains the finite baseline. The present improvement is the
new structural theorem in cutwidth.md Section 11, with the complete
quantitative local probability proofs in certificate.md Sections 9 and 10.
No previous judgment is a mathematical premise.

All six directly used submissions are declared for attribution. The five
earlier sources are 3f2be7ea-aa22-4f7e-98e9-23b40557dbcf (commit
737818f9911e00e9d4909c644c771e0724fb6c37),
0514f3af-5d65-4556-82df-a89e2a3ae4e1 (commit
d18ef5c1b86c58ca9f0bdb93085f9b6f2698c55d),
ac0dd6e1-8a73-4682-b7aa-67502b2de08b (commit
77f07bfd2fd633cf98e3737575ba188d20e95c6a),
38e961ae-bd55-4566-85ad-17112333f84a (commit
f6d3da95b1a680e9bb986008bd6d0890081cdaa7), and
fcf21648-1b1f-463c-b5bc-c50f53f63436 by galarnon42 (commit
81470c9513639af985c88699dbcd2e2efbef7463). The fourth source's hosted
Decimal timeout is part of the diagnostic history, not an accepted premise.
The successive numerical changes include exact integer contractions,
Gaussian-moment errors, positive completions, localized quadrature,
triple-outlier tails and rescaling of the tensor state.

The circuit-to-cutwidth route is due to McGuire, Schlesinger and Williams.
Arthur Ostapenko's snapshot 5fef109bf7f1b664295244a0dee6df066ed2232d supplies
the retained norm-one family, circuit transfer, Kollar--Ronyai--Szabo
exposition and Hellinger comparison, and credits mooselumph's f4f2c464
submission for the validated-prefix device. The original direct source
also credits fd902316-f679-4965-aca9-419a9ee5ec22, source
aeabecef2df58965408d952605327debee2476e4. Matching-cut submodularity and
canonical covers have earlier literature, including Jeong--Saether--Telle,
arXiv:1507.02384; related Gaussian sweeps appear in Gaspers,
arXiv:2611.05146. All needed nonstandard arguments are supplied here.

The new proof combines two local effects. A buffered crossing-cover pattern
forces mixed difference -2 and gives quantitative loss when profile
covariances increase. Another buffered pattern forces a majority update to
save a cut edge. Hellinger stability along short paths connects these two
effects for arbitrary simple cubic graphs. The same deterministic independent
set of local centers works at every threshold, and its improvements persist
through the entire canonical-cover interpolation. Strict loss away from zero
and a fixed threshold mesh make the ordering bound uniform. The multigraph
and circuit transfers then apply at gamma_0.

The exact increase over the immediate predecessor is

    315323783211583831670874691/65323783211583831670874691 - 315323783211583832170874691/65323783211583832170874691
      = 125000000000000000000000000000000000/4267196653074001690183100900473057770952638524345481 > 0,

approximately 2.92932363241221944e-17. This is a small structural improvement.
No optimality claim, superlinear circuit theorem, or numerical evaluation of
an infinite path sum is part of the declaration.

## Evidence and execution contract

The finite path premise remains

    Gamma_65 <= 0.1306475664231676643417493807 < gamma_base.

The additional quantitative premises, proved analytically in certificate.md
and supported by exact rational arithmetic, are

    q=3/100000, theta=1/100, zeta_c=1/700000,
    p=1/100000000000000000, zeta_g=1/35000,
    eta=99/70000000, alpha=1/6370000,
    A > 47223/7000000000000000,
    B_* > 1429569/70000000000000000,
    min(A alpha,p/8,B_*) > delta_g.

The curvature probability is greater than q/2 after conditioning the two
root values to equal the threshold; the gain probability is greater than
p/2 without conditioning. Both estimates hold throughout |threshold|<=1/100
and the stated entrywise covariance tolerances, including singular full
covariance matrices. The root pair needed for conditioning is nonsingular.

certify.py runs one declared experiment on the sole input {}. It uses only
the Python standard library. Its path computation uses 80 Gauss nodes,
65 terms, truncation L=3, 100-digit Decimal initialization, and outward
exact integer contractions at scale 2^128 with rescaling 4^(j-1). Its new
local-event computation uses downward integer lower sums and exact rational
inequalities. Successful exit requires both finite certificates and the
exact structural and eventual-coefficient comparisons. A checker or these
computations alone do not prove the theorem: the full graph, algebraic,
family, circuit, uniformity and slack arguments are supplied separately.

The first two scored claims contain the identical complete declaration and
the same sorted seven evidence files. Local structural and candidate
execution records are kept outside submission/. The verifier is pinned to
the protected workflow revision fa3c2044e179e919bfb09b890594f21458ef142a.
Docker is unavailable locally; a limited host execution is distinct from
the required fresh hosted experiment and configured mathematical review.
Judgment reuse is disabled and no prior_results are declared.

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

## 4. Proof dependencies and uniformity

All nonstandard graph, Gaussian, layout and numerical arguments used by the
declaration are proved in the supplied evidence. The retained algebraic proof
states its textbook commutative-algebra premises. The graph proof gives the
finite bipartite matching/cover argument used in its cover lifting and Hall
steps. No external submission has to be retrieved to complete the proof.

1. krs.md proves the finite-field product-system bound for the rectangle
   argument.
2. circuit.md proves the fixed language's P membership, the rectangle lower
   bound, whole-circuit wire-graph simulation, and the transfer from any
   graph hypothesis H(gamma,h).
3. cutwidth.md Sections 0--10 prove the baseline graph theorem at gamma_base.
   Section 11 proves Theorem 2 at gamma_0=gamma_base-delta_g, including one
   common nonnegative nondecreasing sublinear remainder for all cubic
   multigraphs.
4. certificate.md Sections 1--8 prove the finite inequality Gamma_65<gamma_base;
   Sections 9--10 prove the quantitative curvature and gain lemmas. certify.py
   verifies their stated finite arithmetic with enclosed rounding errors.

The graph theorem applies to every cubic multigraph, including loops,
parallel edges, bridges and disconnected graphs. Finite covers are only
expectation devices; the actual finite-range field, independent set of
gain centers and final ordering live on the original simple graph before
the extension. The threshold convention for the gained cost is b_t^+=1[Z>t]
and prefixes are taken in decreasing key order. Only F, not the gain count,
uses complement symmetry.

The constants delta_g, q, p, theta, zeta_c, zeta_g and J=65 are fixed first.
For each requested error r>0, radius, positive coordinate noise and the
finite threshold mesh are fixed before graph size. The fixed deletion
constant is C_65=3597115094373362565121. Its size changes the unspecified
thresholds, but not the uniform limit. Letting r decrease through arbitrary
positive values preserves the full delta_g when defining the common h_0.
There is no exchange of graph-size and profile-radius limits.

## 5. Assembling the eventual circuit theorem

### 5.1 A strictly stronger intermediate coefficient

Put gamma_0=130647566423167663341749381/1000000000000000000000000000. Theorem 2 of cutwidth.md,
proved in its Section 11, supplies one nonnegative nondecreasing h_0 with
h_0(D)/D->0 such that every cubic multigraph on D vertices has cutwidth
at most gamma_0 D+h_0(D). This is exactly H(gamma_0,h_0) in circuit.md;
its hypotheses include all connected cubic multigraphs needed there.
The finite numerical bound is below gamma_base; the structural theorem
is the additional reason the smaller gamma_0 is available.

Theorem 1 of circuit.md, applied to this h_0 and the fixed evaluator in
Section 1, gives for EVERY integer N>=256

    C_B2(L_N) >= c_0 N-R_0(N),
    c_0 = 1+1/(2 gamma_0),       s_0 = ceil(c_0 N),
    R_0(N) = 5sqrt(N)+(1/2)sqrt(N)log2(N)+2
       +[7sqrt(N)+sqrt(N)log2(N)+6+3log2(2s_0)+h_0(2s_0)]/(2 gamma_0).

Here gamma_0, c_0 and h_0 are fixed. Since s_0=O(N) and tends to infinity,

    h_0(2s_0)/N = [h_0(2s_0)/(2s_0)] * [2s_0/N] -> 0.

Every other displayed summand of R_0(N)/N also tends to zero. Consequently
R_0(N)/N->0 along ALL integer lengths, with the same R_0 for every circuit
computing L_N. Exact arithmetic gives

    c_0 = 1+500000000000000000000000000/130647566423167663341749381
        = 630647566423167663341749381/130647566423167663341749381.

### 5.2 Absorbing the sublinear loss with fixed positive slack

Define the declared coefficient

    c_* = 1+1/[2(gamma_0+1/1000000000000000000000000000)]
        = 315323783211583831670874691/65323783211583831670874691.

The exact difference from the independently proved intermediate coefficient is

    delta_c = c_0-c_* = 250000000000000000000000000/8534393306148003315042418524038500659446754377816271 > 0.

All displayed rational coefficients are reduced, as checked by integer gcd.
Since R_0(N)/N->0 and delta_c/2 is a fixed positive number, the definition
of convergence provides one integer N_0>=256 such that EVERY integer
N>=N_0 satisfies R_0(N)/N<=delta_c/2. At all those lengths,

    C_B2(L_N) >= c_0 N-R_0(N)
              >= (c_0-delta_c/2)N
               = (c_*+delta_c/2)N
              >= c_* N.

This proves precisely the version-2 eventual_lower_bound with leading term
c_* N and explicit remainder terms []. The existence of N_0 follows from
one uniform convergence statement. No explicit rate or threshold is assumed.
The strictly larger coefficient c_0 is essential; the argument does not
infer a zero remainder from a same-coefficient liminf bound.

The evaluator and language are unchanged throughout, and the result holds
at every integer N>=N_0. Dividing by N proves

    liminf over all integer N tending to infinity C_B2(L_N)/N >= c_*.

Thus the declared finite projection has coefficient c_* and exponent 1.
The verifier and witness bound in Section 2 prove NP eligibility for this
same fixed language. These statements prove the entire joint declaration
in bound.json, including both score coordinates and zero explicit remainder.
QED.

## 6. Sources and reproduction

- Immediate dependency: submission ae372f38-09d4-452a-af05-4954cdd2f168,
  https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/2a8b03fe0be18fdf4918dc5fc34da51534f1aca1/submission

- Earlier numerical dependency: submission 3f2be7ea-aa22-4f7e-98e9-23b40557dbcf,
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
