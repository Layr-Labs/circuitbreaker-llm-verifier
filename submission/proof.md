# Certified cubic cutwidth refinement and coefficient 39415472909/8165472909

## Claim and scope

For every cubic multigraph G on D vertices,

    cw(G) <= (8165472909/62500000000)D + h(D),

where h is one nonnegative nondecreasing function with h(D)/D tending to zero.
Loops, parallel edges, bridges and disconnected graphs are included. For the
single fixed total language L in P defined in Section 1, this proves

    liminf over ALL integers N tending to infinity C_B2(L_N)/N
        >= 39415472909/8165472909.

The coefficient is 4.82709003486573198794259816948... . The full B2 model permits
all 16 binary gate types, repeated pins, constants, and arbitrary sharing.
The evaluator and empty-witness NP verifier are fixed independently of every
asymptotic parameter. No explicit lower-order rate or practical threshold is claimed.

## Provenance and contribution of this revision

The direct source is CircuitBreaker submission
`fcf21648-1b1f-463c-b5bc-c50f53f63436`, by galarnon42, at exact repository commit
`81470c9513639af985c88699dbcd2e2efbef7463`. Its supplied graph, circuit,
finite-field and interval arguments are retained here and are mathematical
premises to be reviewed anew. No previous acceptance judgment is a premise.
The dependency is declared separately to Yukon for attribution.

This revision increases the fixed odd path truncation from 13 to 23 and the
Gauss-Legendre node count from 72 to 88. It evaluates the same tensor rule
with quadratic rather than cubic storage, by factoring its positive transition
kernel before contraction. The new finite certificate proves a smaller graph
constant gamma=8165472909/62500000000=0.130647566544. The required finite graph
error constant is updated to C_23=1+3*23*2^22=289406977. The graph proof applies
to every fixed odd J, so this parameter change preserves its quantifiers.
The new contribution is this numerical refinement and its implementation;
the graph reduction, explicit language and algebraic theorem are inherited.

The source's prior credits are preserved. McGuire, Schlesinger and Williams
supply the circuit-to-cutwidth route. Arthur Ostapenko's snapshot
5fef109bf7f1b664295244a0dee6df066ed2232d supplies the retained validated norm-one
language, circuit transfer, Kollar--Ronyai--Szabo exposition and Hellinger
pushforward comparison, and credits mooselumph's submission f4f2c464 for the
validated-prefix device. The direct source also builds on submission
fd902316-f679-4965-aca9-419a9ee5ec22 (source aeabecef2df58965408d952605327debee2476e4).
These ancestor credits describe the source's provenance; the complete package
used for this revision was obtained from the direct source above.

In particular, global crossing covers, the Hall interpolation, exact finite-cover
transfer, and the Gaussian-wave graph argument are inherited from that source.
Matching-cut submodularity and canonical covers have earlier literature,
including Jeong--Saether--Telle, arXiv:1507.02384; related Gaussian sweeps appear
in Gaspers, arXiv:2611.05146. All needed nonstandard arguments are supplied in
this package, with the standard textbook premises explicitly stated.

Compared with the source's coefficient 12613/2613, the exact increase is

    39415472909/8165472909 - 12613/2613
        = 168990000/2370708967913 > 0.

This is about 0.0000712824738453, or 0.0014767391 percent. It is a small
improvement of the same asymptotic circuit theorem, not a new exponent.

## Evidence and execution status

The only numerical premise is the finite inequality

    Gamma_23 <= 0.130647566543255
             < 0.130647566544 = 8165472909/62500000000.

`certificate.md` proves the density, finite sum, tensor rule and all analytic
and arithmetic error bounds. `certify.py` computes rigorous intervals using
only the Python standard library. Its declared input is {} and its defaults
are 88 nodes, J=23, L=3 and 100-digit Decimal arithmetic. It exits zero and
prints `certified: true` only when the upper endpoint lies strictly below the
exact target. This establishes only that finite numerical premise.

The complete theorem is assembled below from the supplied supporting proofs
in `cutwidth.md`, `circuit.md`, `krs.md` and `certificate.md`. Both scored claims
reference the identical complete declaration and all seven evidence files.
Local exact-byte execution and advisory review are recorded outside the
submission. The locally installed verifier is pinned to
fa3c2044e179e919bfb09b890594f21458ef142a. Docker is unavailable in this workspace,
so local host execution is distinct from the digest-pinned hosted experiment.
No hosted acceptance or independent human approval is asserted here.

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
3. `cutwidth.md` proves H(gamma,h) with gamma=8165472909/62500000000 from the sole finite
   numerical inequality Gamma_23<8165472909/62500000000.
4. `certificate.md` and `certify.py` establish that finite inequality.

The graph theorem applies to every graph, not just graphs of large girth.
Finite covers are used only for the expectation comparison, and the actual
finite-range random ordering lives on the original graph. The field radius,
noise and threshold mesh are fixed before graph size tends to infinity. The
resulting remainder is uniform over graphs.

## 5. Assembling the circuit theorem

Put gamma=8165472909/62500000000. Theorem 1 of `cutwidth.md` supplies one nonnegative
nondecreasing h with h(D)/D->0, such that every connected cubic multigraph on D
vertices has cutwidth at most gamma D+h(D). Hence H(gamma,h) of `circuit.md` holds.

Theorem 1 of `circuit.md`, applied to the evaluator in Section 1, gives for all
integers N>=256,

    C_B2(L_N) >= c N-R(N),   c=1+1/(2gamma),   s_0=ceil(cN),
    R(N)=5sqrt(N)+(1/2)sqrt(N)log2(N)+2
         +[7sqrt(N)+sqrt(N)log2(N)+6+3log2(2s_0)+h(2s_0)]/(2gamma).

Since s_0=O(N), h(2s_0)=o(N), and every other term of R is o(N), R(N)/N->0
along ALL integer lengths. Therefore liminf C_B2(L_N)/N>=c. Exact arithmetic gives

    c=1+62500000000/16330945818=1+31250000000/8165472909
     =39415472909/8165472909.

The fraction is reduced since gcd(31250000000,8165472909)=1. In particular c>4.8270900348.
The uniform evaluator and witness bound in Section 2 give NP eligibility for
this same language. This proves the entire joint declaration in `bound.json`,
including both score coordinates. QED.

## 6. Sources and reproduction

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
