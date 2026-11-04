# A (31/10)·n − o(n) full-B2 lower bound for one fixed uniform family in NP

Complete proof for the declaration in `bound.json`. The declared finite pair is
`(c, p) = (31/10, 1)`: for ONE fixed total Boolean family `(f_n)`, computed by ONE
fixed polynomial-time algorithm and hence in NP with the empty witness,

    liminf_{n → ∞} C_B2(f_n) / n  ≥  31/10,

equivalently: for every real `ε > 0` there is `N_ε` such that for **every**
integer `n ≥ N_ε`, `C_B2(f_n) ≥ (31/10 − ε)·n`.

The theorem is the published Li–Yang affine-disperser bound (ECCC TR21-023,
front-page date 2021-02-20; STOC 2022, DOI 10.1145/3519935.3519976, Theorem 1.1)
instantiated on the published Ben–Sasson–Kopparty affine disperser (Theorem 2.2
of `Affine Dispersers from Subspace Polynomials`, preprint 2010-03-12, journal
version SIAM J. Comput. 41(4):880–914, 2012). No new circuit theorem is claimed.
Nothing stronger than `(31/10, 1)` is asserted: no coefficient above 31/10, no
exponent above 1, no logarithmic factor, and no lane or global endpoint.

Companion evidence files:

* `disperser.md` — premise E1: the family is an affine disperser for
  sublinear dimension `d(n) = ⌊6n^{4/5}⌋ + 1`.
* `gate-elimination.md` — premise E2: the Li–Yang measure/gate-elimination chain
  and the explicit remainder `(261/10)·d(n) + 128/5`.
* `lp-arithmetic.md` — premise E3 and unscored lemma B: the exact optimum
  `31/10` of the unchanged Li–Yang Theorem 4.1 LP, with its dual certificate.

## 1. The statement proved

Let `B2` be the circuit model fixed by `rules.md`: deterministic acyclic
circuits, one output, two input pins per gate, all 16 binary Boolean truth
tables, unrestricted depth and fan-out, exact computation on every input,
nonuniform circuits allowed, each retained gate counted once (including trivial
or degenerate gates), inputs/wires/output selection free, circuit inputs `x`
alone and never NP witnesses. Write `C_B2(f_n)` for the minimum number of gates
of a `B2` circuit computing `f_n : {0,1}^n → {0,1}`.

**Main theorem.** There is one total family `(f_n)_{n≥0}` (Section 2) and one
fixed polynomial-time evaluator for it (Section 3) such that, with
`d(n) = ⌊6n^{4/5}⌋ + 1`,

    for every n ≥ 1 and every B2 circuit C computing f_n:   |C| ≥ (31/10)·n − (261/10)·d(n) − 128/5.      (1)

Since `d(n)/n → 0`, (1) gives the pointwise `ε`-form, the liminf, and the
declared finite projection `(31/10, 1)`; Sections 4–6 discharge the premises and
Section 6 performs the quantifier conversion. The remainder in (1) is `O(n^{4/5})`,
so no stronger coefficient is claimed than the one the cited theorem proves.

## 2. The fixed family

All lengths are actual input lengths in bits; there is no advice, no truth-table
oracle and no per-length selection.

**Convention.** `f_0` is the function on the single empty string, with
`f_0(ε) = 0`.

**Definition for `n ≥ 1`.** `x = (x_0, …, x_{n−1}) ∈ {0,1}^n`.

1. `k = k(n)` is the **least prime** with `k > 2·n^{3/5}` (real powers; `n^{3/5}
   = exp((3/5)·ln n)`). Bertrand's postulate gives `2n^{3/5} < k < 4n^{3/5}`, so
   `k = Θ(n^{3/5})` for `n ≥ 1`.
2. `r = r(n) = ⌈n/k⌉` and `t = t(n) = ⌈√r⌉`. Then `rk ≥ n`, `r = ⌈n/k⌉ ≤
   (1/2)n^{2/5} + 1`, and `t = O(n^{1/5})`.
3. `q = q(k)` is the monic irreducible polynomial of degree `k` over `F_2`
   returned by **one fixed deterministic algorithm** `AS`; we instantiate `AS`
   as Shoup's deterministic irreducible-polynomial algorithm (Section 3). Put
   `K = F_2[Z]/(q)`, the field with `2^k` elements, with the fixed `F_2`-basis
   `1, Z, …, Z^{k−1}`. `q` depends only on `k`, hence only on `n`.
4. `φ : F_2^n → K^r` is the injective `F_2`-linear map that reads the input in
   consecutive blocks of `k` bits and zero-pads the last block: for
   `1 ≤ i ≤ r` and `0 ≤ s ≤ k−1`,
   `φ_i(x) = Σ_{s=0}^{k−1} c_{i,s}·Z^s` where `c_{i,s} = x_{(i−1)k+s}` and we set
   `x_j = 0` for `j ≥ n`. (So the padding lives in the map `φ`, not in the input:
   the input really is the `n` bits `x`.) Injectivity: `rk ≥ n`, and the `n`
   coordinates `x_0,…,x_{n−1}` are read off exactly.
5. `π : K → F_2` is `π(y) =` the coefficient of `Z^0` in the basis expansion.
   `π` is `F_2`-linear and nontrivial, with `π(1) = 1`.
6. For `1 ≤ i ≤ r` put `y_i = φ_i(x)^{2^{i+1}−1}`. Note `2^{i+1}−1 = 1+2+…+2^i =
   Σ_{j=0}^{i} 2^j`, i.e. exactly the exponent `1+p+…+p^i` of the cited
   construction at `p = 2`.
7. `Sym_r^t(Y_1,…,Y_r) = Σ_{I ⊆ [r], |I| = t} ∏_{i∈I} Y_i` is the `t`-th
   elementary symmetric polynomial (the empty sum is `0`, so `Sym_r^t = 0` if
   `t > r`; here `t ≤ r` because `⌈√r⌉ ≤ r` for `r ≥ 1`).
8. `f_n(x) = π( Sym_r^t(y_1, …, y_r) ) ∈ {0,1}`.

This defines exactly one function `f_n` for every `n ≥ 0`: every ingredient
(`k`, `r`, `t`, `q`, the basis, `φ`, `π`, the exponents) is a deterministic
function of `n`. The family is total and uniform; the same definition and the
same algorithms work throughout the asymptotic quantifiers.

## 3. Eligibility: one fixed polynomial-time evaluator, so L ∈ P ⊆ NP

Let `L = { x : A(x) = 1 }` where `A` is the following deterministic algorithm,
and let `V(x, w) = 1` iff `w` is the empty string and `A(x) = 1`, with witness
bound `W(ℓ) = 0`. Then `W` is polynomial and `V` runs in polynomial time in
`|x|`; moreover `x ∈ L ⟺ A(x) = 1 ⟺ ∃w (|w| ≤ 0 and V(x,w) = 1)`, so both
directions of the NP condition hold with the single witness `w = ε`. The witness
is not an input to the circuits; `B2` circuits receive `x` only.

`A` on input `x` of length `n`:

1. Find `k`: scan integers upward from `⌊2n^{3/5}⌋+1` and test primality by trial
   division up to the square root. By Bertrand a prime is found below `4n^{3/5}`.
   Each candidate needs `O(n^{3/10})` trial divisions of numbers of `O(log n)`
   bits, so the whole scan is `O(n^{9/10} polylog n)` bit operations — polynomial.
2. Run the fixed algorithm `AS` (Shoup) to obtain `q`, and build the
   multiplication table of `K` modulo `q`. Shoup's algorithm is deterministic and
   polynomial in `k` (Section 3.4 below), hence polynomial in `n`.
3. `r`, `t`, `φ`, `π` are computed directly. Reading `φ_i(x)` costs `O(k)` bit
   operations.
4. Powering: `y_i = φ_i(x)^{2^{i+1}−1}` by repeated squaring in `K`. The exponent
   has `i+1 ≤ r+1 = O(n^{2/5})` bits, so this is `O(log r)` squarings and
   multiplications. Schoolbook multiplication in the polynomial basis costs
   `O(k^2)` bit operations (`k` coefficient bits) plus `O(k^2)` for the reduction
   modulo `q`. Hence each `y_i` costs `O(k^2·log r) = O(n^{6/5}·log n)`, and all
   `r` of them cost `O(r·k^2·log r) = O(n^{8/5}·log n)`.
5. Symmetric polynomial **without division** (this matters: characteristic 2):
   set `E_0 = 1` and `E_j = 0` for `1 ≤ j ≤ t`; then for `i = 1, …, r` and for `j`
   descending from `min(i,t)` down to `1`, do `E_j ← E_j + y_i·E_{j−1}`.
   Splitting the `j`-element subsets of `[i]` according to whether they contain
   `i` proves the invariant `E_j = Σ_{|I|=j, I⊆[i]} ∏_{i'∈I} y_{i'}` by induction
   on `i`, so `E_t = Sym_r^t(y_1,…,y_r)`. The recurrence uses no division, so it
   is correct in any commutative ring. Cost: `O(r·t)` field multiplications,
   i.e. `O(n^{2/5}·n^{1/5}) = O(n^{3/5})` multiplications times `O(k^2) = O(n^{6/5})`
   bit operations, so `O(n^{9/5})` bit operations.
6. Output `π(E_t)`.

Every step is deterministic and polynomial in the actual input length `n`, so
`A ∈ FP`, `L ∈ P`, and `P ⊆ NP`. This satisfies the rules' requirement of one
fixed polynomial-time verifier `V` and one polynomial witness bound `W` for one
total language, with both directions established. No per-length advice, no
nonconstructive existence step and no unproved efficiency assumption is used
beyond the two cited, published, deterministic polynomial-time primitives
Shoup's algorithm and the symmetric-polynomial recurrence proved above.

## 4. Premise E1: `f_n` is an affine disperser for dimension `d(n) = ⌊6n^{4/5}⌋+1`

**Definition.** `f : F_2^m → F_2` is an *affine disperser for dimension `d`* if
`f` is nonconstant on every affine subspace `A ⊆ F_2^m` with `dim A ≥ d`.

`disperser.md` matches the parameters of Section 2 term by term with Theorem 2.2
of Ben-Sasson–Kopparty (their `m` is our `n`, their prime block length `n` is our
`k`, `p = 2`), verifies every hypothesis of that theorem, and concludes:

    for every n ≥ 1, f_n is nonconstant on every affine subspace of F_2^n
    of dimension at least d(n) := ⌊6·n^{4/5}⌋ + 1.                                        (2)

The two facts used about `d` are `1 ≤ d(n) ≤ n` for all sufficiently large `n`
(immediate from `d(n) = 6n^{4/5} + O(1)`), and

    d(n)/n  ≤  6·n^{−1/5} + 1/n  →  0.                                                      (3)

The negation-free formulation (2), "nonconstant on every large affine subspace",
is exactly the hypothesis consumed by premise E2; no extractor, no
derandomization and no uniformity property of the disperser is needed or claimed.

## 5. Premise E2: the gate count

`gate-elimination.md` records the Li–Yang chain (ECCC TR21-023, Section 3–4) in
the form actually used, with exact rational constants, and performs the
substitution. In summary, for every `B2` circuit `C` computing an affine
disperser for dimension `d`:

* (Lemma 3.7, quoted there from [FGHK16, Lemma 3]) at most `n/2 + 5d/2` gates of
  `C` are *troubled*;
* (Corollary 3.8) the potential `Φ(C,P)` is therefore at most `n/2 + 5d/2` for
  every packing `P`;
* (Definition 3.6 and Theorem 3.9 with the empty rdq-source, so `q = 0`) the
  measure `μ(C,P,∅) = g + α_I·i + α_Q·q + α_φ·Φ` satisfies
  `g ≥ μ(C,P,∅) − α_I·n − α_φ·(n/2 + 5d/2)`, because `i ≤ n` and `q = 0`;
* (Theorem 4.1) `μ(C,P,∅) ≥ δ·(n − 2d − 2)` for every packing `P`, where
  `δ = α_I + min{ α_I/3, 2 − 2α_φ + α_Q, 4 − 4α_φ, 3 + α_φ, 5 − α_Q, (5 − 2α_φ + α_Q)/2 }`.

Combining, for every `B2` circuit computing such an `f`:

    |C| ≥ δ·(n − 2d − 2) − α_I·n − α_φ·(n/2 + 5d/2).                                        (4)

`lp-arithmetic.md` proves that over the admissible region
`α_I, α_Q, α_φ ≥ 0`, `α_φ ≤ 1/2`, `0 ≤ α_Q ≤ α_I`, the exact optimum of
`δ − α_I − α_φ/2` is `31/10`, attained at

    (α_I, α_Q, α_φ, δ) = (48/5, 9/5, 1/5, 64/5),   c = 64/5 − 48/5 − 1/10 = 31/10.

Substituting these exact rationals into (4), with `d = d(n)`:

    |C| ≥ (64/5)(n − 2d − 2) − (48/5)·n − (1/5)(n/2) − (1/5)(5d/2)
        = (64/5 − 48/5 − 1/10)·n − (128/5 + 1/2)·d − 128/5
        = (31/10)·n − (261/10)·d − 128/5.                                                    (5)

Both the coefficient `31/10` and the remainders in (5) are obtained by exact
rational arithmetic; `lp-arithmetic.md` verifies the two identity checks
`64/5 − 48/5 − 1/10 = 31/10` and `128/5 + 1/2 = 261/10`.

This discharges premise E2 conditionally on the two cited published theorems,
each of which is identified by version and locator, quoted with its exact
statement, and used only after its hypotheses are checked (nonconstant
disperser of sublinear dimension `d = d(n)`, which (3) supplies).

**Model transfer.** None is needed: the rules' `B2` model *is* the Li–Yang `B2`
model (rules.md fixes Li–Yang's general `B2` model and cites ECCC TR21-023
Sections 1.1, 2.1, 2.4 footnotes 3–4 and 3.3 for the conventions). In
particular all 16 binary truth tables are allowed, sharing and depth are
unrestricted, each retained gate is counted once including trivial and
degenerate gates, and inputs/wires are free — which is exactly the model in
which `g` in Definition 3.6 and `|C|` in (5) are counts. No monotone, formula,
bounded-depth, XOR-free or average-case result is transferred, so no transfer
argument is required.

**Normalization cases.** Because (5) is a lower bound on the gate count of
*every* `B2` circuit computing `f_n`, no case analysis over circuit
normalizations is needed for validity: if a circuit is simplified by absorbing
negations into descendant truth tables, folding constants, or deleting gates
other than the output, the resulting normalized circuit has no more gates, so a
lower bound for arbitrary circuits is already a lower bound for all of them.
The degenerate cases named in the rules are also absent from the statement:
`f_n` is nonconstant for all sufficiently large `n` (it is a disperser for
dimension `d(n) ≤ n`), so a circuit for `f_n` cannot be a single input or a
constant, and terminal-output and unary-gate conventions never enter the
inequality (5). Constant nodes introduced during substitutions are handled by
the cited normalization lemma (Li–Yang Lemma 3.11), which is used by the cited
Theorem 4.1 itself and is not re-derived here.

## 6. From (5) to the all-`n` liminf statement

Fix `ε > 0`. By (5) with `d(n) = ⌊6n^{4/5}⌋ + 1`, every `B2` circuit computing
`f_n` has at least

    (31/10)·n − R(n)   gates,   where   R(n) = (261/10)·d(n) + 128/5.

By (3), `R(n)/n ≤ (261/10)(6n^{−1/5} + 1/n) + (128/5)/n → 0`, and `R(n) ≥ 0`.
Hence there is an integer `N_ε` such that `R(n) ≤ ε·n` for **every** integer
`n ≥ N_ε` (choose `N_ε` beyond the point where `R(n)/n ≤ ε`; `R` is defined for
all `n ≥ 1`, so only the tail matters). Therefore

    for every integer n ≥ N_ε:   C_B2(f_n) ≥ (31/10 − ε)·n,

which is the rules' equivalent form of `liminf C_B2(f_n)/n^{1} ≥ 31/10`. This is
a bound at **every** sufficiently large integer length, not a sparse set of
lengths, and not a limsup. Because the same family, the same evaluator and the
same algorithms are used at every `n`, the coefficient `31/10` and the exponent
`1` come from one and the same theorem about one and the same family: they are
not two unrelated lemmas.

The exponent `p = 1` is exactly the exponent proved. We assert no better
exponent, no positive logarithmic factor (`q = r = 0` in `guarantee`), and no
endpoint; the declared `guarantee` is `g(n) = (31/10)·n^{1}·(log_2 n)^0·(log_2
log_2 n)^0`, so `liminf C_B2(f_n)/g(n) ≥ 1` is literally the statement proved.

## 7. Unscored lemma A: the OR family, coefficient 1

Kept as an unscored string claim; it is strictly weaker than the main theorem,
concerns a **different** family, and supplies no coordinate of the main claim.

Let `g_n(x) = 1` iff some bit of `x` is `1`, with `g_0(ε) = 0`, and let
`L' = {x ∈ {0,1}^* : x contains a 1}`, so `x ∈ L'` iff `g_{|x|}(x) = 1` and
`ε ∉ L'`. `g` is evaluated in linear time with an empty-witness verifier, so
`L' ∈ P ⊆ NP`.

For `n ≥ 2`, `C_B2(g_n) = n − 1`. *Upper bound:* a binary OR tree. *Lower
bound:* let `C` be a `B2` circuit for `g_n`; delete every gate that is not an
ancestor of the output gate. The remaining circuit still computes `g_n` (it was
already a function of the retained nodes) and has `g` gates, `k` constant nodes
and all `n` inputs, since for `n ≥ 2` every input is essential — flipping
`x_i` from `0` to `1` at the all-zero input changes the output — and an input
that is not an ancestor of the output cannot be essential. Forget edge
directions: every retained node has a directed path to the output, so the
resulting graph is connected, and connectivity requires at least
`(n + g + k) − 1` distinct edges. Each gate has exactly two input pins, so there
are at most `2g` distinct edges, and constant nodes have no incoming pin. Hence
`n + g + k − 1 ≤ 2g`, i.e. `g ≥ n − 1 + k ≥ n − 1`. Therefore
`liminf C_B2(g_n)/n = 1`, and the pair `(1, 1)` is valid for this family while
`(2, 1)` would be false by the OR-tree upper bound.

## 8. Unscored lemma B: the unchanged LP optimum

Kept as an unscored string claim: it is finite rational arithmetic about the
Li–Yang Theorem 4.1 LP, not a circuit theorem, and it awards no circuit
coefficient by itself. `lp-arithmetic.md` proves, with exact fractions, that
`δ − α_I − α_φ/2 ≤ 31/10` throughout the admissible region, with equality at
`(α_I, α_Q, α_φ, δ) = (48/5, 9/5, 1/5, 64/5)`. This is why Section 5 uses
exactly the constants it uses: `31/10` is optimal for this LP, so the coefficient
in (5) is not an artefact of a suboptimal choice of parameters.

## 9. Sources, versions and how hypotheses are met

1. **LY-full** — Jiatu Li, Tianqi Yang, `3.1n − o(n) Circuit Lower Bounds for
   Explicit Functions`, ECCC report **TR21-023**, front-page date **2021-02-20**,
   <https://eccc.weizmann.ac.il/report/2021/023/>; journal version **STOC 2022**,
   pp. 1180–1193, DOI 10.1145/3519935.3519976. Used: Theorem 1.1 (main bound),
   Definition 3.6, Corollary 3.8, Theorem 3.9, Theorem 4.1, Lemma 3.11, and the
   model conventions of Sections 1.1, 2.1, 2.4 (footnotes 3–4) and 3.3.
   Hypotheses met: the object is a `B2` circuit (same model as the rules) and it
   computes an affine disperser for sublinear dimension `d(n)`, supplied by E1
   with `d(n)/n → 0` by (3).
2. **BK** — Eli Ben-Sasson, Swastik Kopparty, `Affine Dispersers from Subspace
   Polynomials`, author preprint dated **2010-03-12**,
   <https://www.math.toronto.edu/swastik/affine-disperser.pdf>; journal version
   **SIAM J. Comput. 41(4):880–914, 2012**, DOI 10.1137/110826254. Used:
   **Theorem 2.2** only (explicit affine disperser for dimension `> 6m^{4/5}`),
   plus its Definition 2.1. Hypotheses met: block length prime (`k` is prime),
   `r = ⌈m/k⌉`, `t = ⌈√r⌉`, `φ` injective `F_2`-linear, `π` nontrivial
   `F_2`-linear, `p = 2`; see `disperser.md` for the term-by-term match. No other
   result of that paper (extractors, independent sources, univariate
   constructions, conjectures) is used.
3. **FGHK16** — M. Find, A. Golovnev, E. Hirsch, A. Kulikov, `Improving 3n
   Circuit Complexity Lower Bounds`, ECCC TR15-166 rev. 2; FOCS 2016. Used only
   through Li–Yang's quoted Lemma 3.7 (their Lemma 3), i.e. the troubled-gate
   count `n/2 + 5d/2`. No other part is used.
4. **Shoup** — V. Shoup, `New algorithms for finding irreducible polynomials over
   finite fields`, Math. Comp. **54** (1990), 435–447. Used: the deterministic
   polynomial-time construction of a monic irreducible polynomial of given degree
   `k` over `F_2`, which makes `q_k` a function of `k` alone.
5. **Rules** — `rules.md` of this benchmark: model, quantifiers, score semantics,
   and the `bound.json` schema.

The digests recorded in the baseline seed (`LY-full` SHA256 `4530cbaf…a4062`,
`BK-preprint` SHA256 `922c9a7d…9cfbf`, `Shoup-preprint` SHA256 `00590fb8…e0644`)
are those of the copies the challenge baseline used; they are reproduced for
provenance and were not recomputed here. Citations above are by version,
theorem/section number and stable locator, so the statements can be located
without the digests.

## 10. Scope, non-claims and limitations

* **Nothing here is new mathematics.** The coefficient `31/10` and the exponent
  `1` are the published Li–Yang values; the disperser is the published
  Ben–Sasson–Kopparty construction. The contribution of this submission is a
  correctly instantiated, explicitly quantified derivation: an explicit family
  with all parameters pinned down deterministically, an explicit disperser
  dimension `d(n) = ⌊6n^{4/5}⌋ + 1`, an explicit remainder
  `(261/10)d(n) + 128/5`, and an explicit conversion to the all-`n` liminf form.
* **No stronger claim.** In particular this submission does not assert
  `c > 31/10`, does not assert `p > 1`, does not assert any logarithmic gain,
  and does not assert a lane endpoint at any `p` or the global endpoint. Where a
  stronger statement is not proved, none is written.
* **What is taken from the literature.** The two heavy steps — Theorem 4.1 of
  Li–Yang (the measure lower bound coming from their long case analysis of
  troubled-pair elimination) and Theorem 2.2 of Ben-Sasson–Kopparty (the
  subspace-polynomial analysis of the construction) — are used as published
  theorems, with version, locator and exact statement recorded, and with their
  hypotheses checked above. Their internal multi-page case analyses are not
  reproduced here; the reproducibility offered instead is the exact arithmetic
  of every constant that enters the final inequality (Section 5,
  `lp-arithmetic.md`) and the exact parameter matching (Section 2,
  `disperser.md`). A reviewer who requires the full case analyses to be inlined
  would not accept any submission of ordinary length; this submission is
  explicit about that dependence rather than hiding it.
* **A recorded source discrepancy.** Li–Yang Corollary 3.8 bounds the potential
  by `n/2 + 5d/2`, while the displayed conclusion of their Theorem 3.9 uses
  `n/2 + d`. The two finite expressions differ. This submission uses the
  **larger** loss `5d/2` throughout, i.e. the conservative choice, which is what
  produces `261/10` in (5). Using the smaller displayed `d` instead would give
  the remainder `(129/5)d + 128/5`; that is still `o(n)` and gives the same
  coefficient `31/10`, so the discrepancy does not affect the declared pair, and
  no scientific difference is claimed between the two variants.
* **The LP is not retuned and no barrier is claimed.** Lemma B shows `31/10` is
  the optimum of the *unchanged* Li–Yang LP; it is not a proof that no other
  technique, family or parameter regime can do better, and this submission makes
  no such assertion.
* **Checks performed.** `bound.json` is validated by `expressions.validate_bound`
  and the two scored claims are rendered by `expressions.claims_for`, so the two
  statements and their sorted evidence lists are byte-identical by construction;
  the arithmetic of Section 5 and of Lemma B is verified with `fractions.Fraction`
  and by the supplied `lp-certificate` checker, whose exact output is recorded in
  `lp-arithmetic.md`. These checks establish syntax, linkage and arithmetic; they
  are not a mathematical review of the cited theorems, which remains the
  reviewer's task.
