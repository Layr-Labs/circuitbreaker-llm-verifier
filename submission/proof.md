# A 0.1341 cubic cutwidth bound and a 6341/1341 circuit lower bound

## Claim and scope

Let L be the one total language defined in Section 1, and let N be the entire input length.
The complete claim is

    L is in P, hence in NP with empty witnesses, and
    liminf over ALL integers N tending to infinity C_B2(L_N)/N >= 6341/1341.

Equivalently, for every epsilon>0 there is N_epsilon such that every full-B2 circuit
computing L_N has at least (6341/1341-epsilon)N gates for every integer N>=N_epsilon.
The coefficient is 6341/1341 = 4.72856077554064... . No explicit lower-order rate or
practical threshold is claimed. The family, evaluator and NP verifier are fixed
independently of epsilon.

The new graph theorem is: every cubic multigraph G on D vertices has

    cw(G) <= (1341/10000)D + h(D),

for one nonnegative nondecreasing h with h(D)/D->0. Loops, parallel edges, bridges
and disconnected graphs are included. `cutwidth.md` proves the graph theorem;
`certificate.md` and `certify.py` give its deterministic numerical certificate.
The exact circuit coefficient is 1+1/(2*(1341/10000))=6341/1341.

## Attribution and the new contribution

The circuit-to-cutwidth reduction is due to McGuire, Schlesinger and Williams.
The validated norm-one family, its circuit transfer and the Kollar--Ronyai--Szabo
write-up are reproduced from ArthurOstapenko's CircuitBreaker work, as contained
in the frontier snapshot 5fef109bf7f1b664295244a0dee6df066ed2232d of
Layr-Labs/circuitbreaker-llm-verifier. That snapshot is the source of the retained
`circuit.md`, `krs.md`, and the family/model exposition below. Its numerical
refinement by lamb356 is not used to establish this theorem. The validated-prefix
device is credited there to mooselumph's submission f4f2c464.

The Hellinger pushforward comparison is adapted from ArthurOstapenko's graph proof.
The new construction uses EDGE-centered squared tree profiles, jointly optimizes
the complementary 2-factor in bounded blocks, identifies its cost with transition
vertex cover, and handles short cycles by periodic unrolling. The resulting
universal bound reduces to five alternating-sign probabilities for six Gaussian
variables. All these ingredients are proved in `cutwidth.md`; the finite integral
certificate is new and independent of the previous frontier's certificate.
Gaussian sweep methods also appear in Gaspers, arXiv:2611.05146.

The previous frontier's cutwidth coefficient is 2452701554003947/2^54,
approximately 0.1361522868894542. The new rational coefficient 0.1341 is smaller.
Its circuit coefficient exceeds 11459900808744939/2452701554003947, approximately
4.672358440853541, by more than 0.0562. These comparisons are exact rational
arithmetic, not assumptions in the proof. No prior submission's judgment or
acceptance is used as a mathematical premise.

## Evidence and execution status

This is a computer-assisted written proof. The only numerical obligation is
gamma_5<1341/10000. The supplied standard-library Python source computes
outward-rounded interval enclosures for an explicitly derived upper majorant U.
The recorded run proves

    gamma_5 <= U <= 0.134085456508304113 < 1341/10000.

The displayed upper endpoint is rounded outward as a decimal rational.
`certificate.md` proves the majorant, the finite integral formulas and every
error bound, and records the source hash and local output. The initial experiment
in `submission.json` runs `certify.py` with the empty object as its input. A
successful run exits zero and emits `certified: true` only after comparing the
upper endpoint with an enclosure of the exact target 1341/10000.

The certificate has been independently reproduced locally in about eight seconds.
It uses only the Python standard library and fits the declared twenty-second
budget on the tested Mac. The pinned Linux/amd64 Docker replay has not been run
because Docker is unavailable locally. No hosted computational pass, configured
provider judgment, human approval or platform acceptance is claimed.

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

The argument is self-contained apart from the familiar textbook facts explicitly
listed in the component files and Petersen's perfect-matching theorem for
bridgeless cubic graphs. The latter applies only after the graph is reduced to a
simple 2-connected cubic kernel. All nonstandard graph and Gaussian comparison
steps used here are proved in `cutwidth.md`.

1. `krs.md` proves the finite-field product-system bound used in the rectangle
   argument, from its stated textbook commutative-algebra facts.
2. `circuit.md` proves P membership, the rectangle lower bound, the wire-graph
   simulation, and the transfer from its graph hypothesis H(gamma,h).
3. `cutwidth.md` proves H(gamma,h) for gamma=1341/10000, using only the numerical
   claim proved by `certificate.md` and the locally executed `certify.py`.

The graph argument does not assume that locally tree-like graphs are extremal.
It proves the needed comparison for every graph using the squared-profile
pushforward, including repeated terminal coordinates. Short cycles incur no
positive-density error. The block size, field radius, noise, threshold grid and
concentration threshold are chosen in that order before the graph; the resulting
remainder is uniform over graphs.

## 5. Assembling the circuit theorem

Put gamma=1341/10000. Theorem 1 of `cutwidth.md` supplies a nonnegative
nondecreasing h with h(D)/D->0 such that every connected cubic multigraph on D
vertices has cutwidth at most gamma D+h(D). Thus H(gamma,h) in `circuit.md` holds.

Theorem 1 of `circuit.md`, applied to the evaluator of Section 1, gives, for all
integers N>=256,

    C_B2(L_N) >= c N-R(N),   c=1+1/(2gamma),   s_0=ceil(cN),
    R(N)=5sqrt(N)+(1/2)sqrt(N)log2(N)+2
         +[7sqrt(N)+sqrt(N)log2(N)+6+3log2(2s_0)+h(2s_0)]/(2gamma).

Here s_0=O(N), so h(2s_0)=o(N). Every other term of R is also o(N). Therefore
R(N)/N->0 along ALL integer lengths and liminf C_B2(L_N)/N>=c. Finally

    c=1+10000/2682=1+5000/1341=6341/1341.

The fraction is reduced: gcd(6341,1341)=gcd(5000,1341)=1 (Euclid). The exact
arithmetic also gives c>4.72856. The uniform evaluator and witness bound from
Section 2 prove NP eligibility for this same language. This proves the entire
joint declaration in `bound.json`, both score coordinates together. QED.

## 6. Sources and reproduction

- McGuire--Schlesinger--Williams: https://github.com/rrwilliams0/circuit-lb-4.6155n,
  especially the circuit-to-kernel reduction and median-edge ordering.
- CircuitBreaker frontier source used for the retained family/circuit/algebra
  exposition: https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/5fef109bf7f1b664295244a0dee6df066ed2232d/submission
- Gaspers: https://arxiv.org/abs/2611.05146.
- Kollar--Ronyai--Szabo, "Norm-graphs and bipartite Turan numbers", Combinatorica
  16(3) (1996), 399--406; the needed theorem is proved in `krs.md`.

From the submission directory, reproduce the numerical certificate with

    python3 -I -B certify.py < /dev/null

A successful run exits 0 and prints one JSON object with `certified: true`.
This finite check does not replace the graph and circuit proofs above.
