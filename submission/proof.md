# Cubic cutwidth 0.13065 and a 12613/2613 circuit lower bound

## Claim and scope

For every cubic multigraph G on D vertices,

    cw(G) <= (2613/20000)D + h(D),

where h is one nonnegative nondecreasing function with h(D)/D tending to zero.
Loops, parallel edges, bridges and disconnected graphs are included. For the
single total language L in P defined in Section 1, this gives

    liminf over ALL integers N tending to infinity C_B2(L_N)/N >= 12613/2613.

The coefficient is 12613/2613 = 4.827018752391886... . The full B2 model permits
all 16 binary gate types and arbitrary sharing. The evaluator and the empty
NP-witness verifier are fixed independently of every asymptotic parameter.
No explicit lower-order rate or practical size threshold is claimed.

## Attribution and the new contribution

The circuit-to-cutwidth reduction is due to McGuire, Schlesinger and Williams.
The validated norm-one language, its circuit transfer and the Kollar--Ronyai--Szabo
exposition are retained from Arthur Ostapenko's CircuitBreaker work in snapshot
5fef109bf7f1b664295244a0dee6df066ed2232d of Layr-Labs/circuitbreaker-llm-verifier.
The snapshot credits mooselumph's submission f4f2c464 for the validated-prefix
device. These credits are preserved in `circuit.md` and `krs.md`. The Hellinger
pushforward comparison also builds on Ostapenko's graph proof. Gaussian sweep
methods appear in Gaspers, arXiv:2611.05146. Maximum-matching cut
submodularity and monotone canonical covers are established ingredients; see
Jeong--Sæther--Telle, arXiv:1507.02384, Section 3. The new graph proof reproduces
the needed facts and does not claim these ingredients as new.

Our earlier 0.1341 construction jointly inserted complementary-cycle edges
around a Gaussian matching-edge order. The present proof uses a different
optimization: at each cut of a Gaussian VERTEX order, take a minimum vertex
cover of ALL crossing edges. It proves the deterministic conversion

    cw(G) <= linear maximum matching width(G) + 6

for loopless cubic graphs. The proof uses monotone least minimizers, a Hall
matching and a tree/cycle interpolation to turn a potentially large jump
between optimal covers into individual edge insertions at constant extra cost.

A second new ingredient transfers the expectation through finite graph covers:
the minimum crossing vertex-cover cost lifts exactly under a finite cover,
whereas Hellinger pushforward increases Gaussian covariances. This permits a
universal high-girth comparison without losing a boundary fraction. On the
limiting cubic-tree Gaussian wave, the eigenvector identity excludes crossing
stars of degree three. A finite alternating path count then bounds the minimum
cover. A finite-radius field, independent jitter and a variance argument supply
one ordering uniformly over all its cuts. These steps and the extension to
all cubic multigraphs are proved in `cutwidth.md`.

The numerical premise is a finite 13-term sum of Gaussian path probabilities.
Its exact AR(2) representation permits tensor Gauss--Legendre quadrature with
two live coordinates; directed Decimal arithmetic and analytic tail/quadrature
bounds certify the result. Neither Monte Carlo estimates nor an infinite
series identity enter the proof.

The new coefficient 0.13065 is 2.5727% below our preceding 0.1341 result; the
circuit coefficient rises from 6341/1341 to 12613/2613, an increase greater than
0.09845. It also improves the earlier frontier snapshot's cubic coefficient
2452701554003947/2^54 and circuit coefficient
11459900808744939/2452701554003947. These are exact rational comparisons.
No prior submission's acceptance or judgment is used as a mathematical premise.

## Evidence and execution status

This is a computer-assisted written proof. The only numerical premise is

    Gamma_13 < 2613/20000.

`certificate.md` proves the density, finite sum, quadrature rule, truncation and
all arithmetic error bounds; `certify.py` supplies the reproducible computation.
The experiment in `submission.json` runs that exact program with input {}.
It exits zero and emits `certified: true` only after its certified upper endpoint
is strictly below the exact rational target. The complete mathematical argument
also includes `cutwidth.md`, `circuit.md` and `krs.md`; a successful computation
alone does not prove those arguments.

Local execution and independent advisory review are recorded outside the
submission. The digest-pinned Linux/amd64 Docker replay is unavailable on this
Mac because Docker is absent. No hosted pass, configured-provider judgment,
human approval or platform acceptance of this candidate is claimed.

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
3. `cutwidth.md` proves H(gamma,h) with gamma=2613/20000 from the sole finite
   numerical inequality Gamma_13<2613/20000.
4. `certificate.md` and `certify.py` establish that finite inequality.

The graph theorem applies to every graph, not just graphs of large girth.
Finite covers are used only for the expectation comparison, and the actual
finite-range random ordering lives on the original graph. The field radius,
noise and threshold mesh are fixed before graph size tends to infinity. The
resulting remainder is uniform over graphs.

## 5. Assembling the circuit theorem

Put gamma=2613/20000. Theorem 1 of `cutwidth.md` supplies one nonnegative
nondecreasing h with h(D)/D->0, such that every connected cubic multigraph on D
vertices has cutwidth at most gamma D+h(D). Hence H(gamma,h) of `circuit.md` holds.

Theorem 1 of `circuit.md`, applied to the evaluator in Section 1, gives for all
integers N>=256,

    C_B2(L_N) >= c N-R(N),   c=1+1/(2gamma),   s_0=ceil(cN),
    R(N)=5sqrt(N)+(1/2)sqrt(N)log2(N)+2
         +[7sqrt(N)+sqrt(N)log2(N)+6+3log2(2s_0)+h(2s_0)]/(2gamma).

Since s_0=O(N), h(2s_0)=o(N), and every other term of R is o(N), R(N)/N->0
along ALL integer lengths. Therefore liminf C_B2(L_N)/N>=c. Exact arithmetic gives

    c=1+20000/5226=1+10000/2613=12613/2613.

The fraction is reduced since gcd(10000,2613)=1. In particular c>4.82701875.
The uniform evaluator and witness bound in Section 2 give NP eligibility for
this same language. This proves the entire joint declaration in `bound.json`,
including both score coordinates. QED.

## 6. Sources and reproduction

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
