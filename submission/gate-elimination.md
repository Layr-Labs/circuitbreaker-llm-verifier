# Premise E2: the Li–Yang gate-elimination chain and the explicit remainder

This file records the exact statements of the Li–Yang chain that are used, checks
their hypotheses against the family of `proof.md` §2, and performs the exact
rational substitution producing

    |C| ≥ (31/10)·n − (261/10)·d(n) − 128/5.                                      (5)

Source: Jiatu Li, Tianqi Yang, `3.1n − o(n) Circuit Lower Bounds for Explicit
Functions`, ECCC report **TR21-023**, front-page date **2021-02-20**
(<https://eccc.weizmann.ac.il/report/2021/023/>); journal version **STOC 2022**,
pp. 1180–1193, DOI 10.1145/3519935.3519976. Section and theorem numbers below are
those of the ECCC version.

Notation: `n` is the actual number of circuit inputs; `C` is a `B2` circuit with
`g = |C|` gates, one output; `f` is an affine disperser for dimension `d`; `B2`
is the model fixed by `rules.md` (all 16 binary gates, two input pins per gate,
unrestricted depth and fan-out, exact on all inputs, gates counted once,
inputs/wires free).

## 1. Theorem 1.1 (main theorem), exact statement

> **Theorem 1.1.** Assume that `C` is a `B2` circuit computing an affine disperser
> for sublinear dimension, then `C` has size at least `3.1n − o(n)`.

`3.1 = 31/10` exactly. The paper obtains it from the chain of §2 below; we
reproduce the chain's statements and the substitution, and we make the `o(n)`
explicit, because the rules of this challenge require a proved
`c·n^p − r(n)` bound with `r(n)/n^p → 0` and require coverage of every
sufficiently large integer length.

## 2. The chain, with the exact constants

**Definition (affine disperser; Li–Yang Definition 2.1, after BK Definition 2.1).**
`f : F_2^n → F_2` is an affine disperser for dimension `d` if `f` is not constant
on any affine subspace of `F_2^n` of dimension `d`.

**Objects used by the measure.** A *packing* `P` of `C` and an *rdq-source* `R`;
`R = ∅` means the source has no quadratic equations, so `q = 0`. A gate is
*troubled* if it is an `∧`-type 1-gate fed by two 2-variables; `Φ(C,P)` is the
*potential* `= (#truly troubled gates) + |P|`. An input is *influential* if it is
a 1⁺-variable or protected; `i` is the number of influential inputs, so `i ≤ n`.

**Definition 3.6 (circuit complexity measure).**

    μ(C, P, R) = g + α_I·i + α_Q·q + α_φ·Φ(C, P),

with `α_I, α_Q, α_φ ≥ 0` constants chosen later.

**Lemma 3.7 (= Lemma 3 of [FGHK16]).** If `C` computes an affine disperser for
dimension `d`, then the number of troubled gates of `C` is at most `n/2 + 5d/2`.

**Corollary 3.8.** If `C` computes an affine disperser for dimension `d`, then
for all packings `P` of `C`,

    Φ(C, P) ≤ n/2 + 5d/2.

**Theorem 3.9.** Assume that computing `f : F_2^n → F_2` requires
`μ(C,P,∅) ≥ μ₀` for any circuit `C` and packing `P`. Then

    |C| ≥ μ₀ − α_I·n − α_φ·(n/2 + d),

and if `f` is an affine disperser for sublinear dimension with `μ₀ = δn − o(n)`,
then `|C| ≥ (δ − α_I − α_φ/2)·n − o(n)`.

**Theorem 4.1.** Let `C` be a circuit computing an affine disperser `f` for
dimension `d`. For all packings `P` of `C`,

    μ(C, P, ∅) ≥ δ·(n − 2d − 2),      where

    δ = α_I + min{ α_I/3,  2 − 2α_φ + α_Q,  4 − 4α_φ,  3 + α_φ,  5 − α_Q,
                   (5 − 2α_φ + α_Q)/2 }.

**Admissibility.** `α_I, α_Q, α_φ ≥ 0`, `α_φ ≤ 1/2`, `0 ≤ α_Q ≤ α_I`,
`δ ≥ 0`.

## 3. Substituting the constants, exactly

Take `(α_I, α_Q, α_φ, δ) = (48/5, 9/5, 1/5, 64/5)`; `lp-arithmetic.md` verifies
that this point is admissible, that five of the six terms of Theorem 4.1 equal
`16/5` there and the second equals `17/5` (so the minimum is `16/5` and the
defining identity for `δ` holds: `48/5 + 16/5 = 64/5`), and
that `31/10` is the *exact optimum* of `δ − α_I − α_φ/2` on the admissible
region.

Apply Theorem 4.1 with `μ₀ = δ(n − 2d − 2)`, then Theorem 3.9 with `R = ∅`
(hence `q = 0`) and with the potential bound of Corollary 3.8, i.e. the
**conservative** `Φ ≤ n/2 + 5d/2` (see §5). For every `B2` circuit `C` computing
`f_n`:

    g = |C| ≥ μ(C,P,∅) − α_I·i − α_Q·q − α_φ·Φ(C,P)
            ≥ δ(n − 2d − 2) − α_I·n − α_φ·(n/2 + 5d/2),

where `i ≤ n` and `q = 0`. Substituting the rationals:

    g ≥ (64/5)·n − (128/5)·d − 128/5 − (48/5)·n − (1/5)·(n/2) − (1/5)·(5d/2)
      = (64/5 − 48/5 − 1/10)·n − (128/5 + 1/2)·d − 128/5
      = (31/10)·n − (261/10)·d − 128/5.

Exact rational identities used, verified in `lp-arithmetic.md`:

    64/5 − 48/5 = 16/5,   16/5 − 1/10 = 32/10 − 1/10 = 31/10,
    128/5 + 1/2 = 256/10 + 5/10 = 261/10.

With `d = d(n) = ⌊6n^{4/5}⌋ + 1` from premise E1 this is (5), and since
`d(n)/n → 0` the remainder is `o(n)`; the coefficient is `31/10` and the
exponent is `1`.

## 4. Hypothesis checks for Theorem 1.1 / Theorem 4.1

1. **`C` is a `B2` circuit.** This is the model of `rules.md`, which fixes
   Li–Yang's general `B2` model (deterministic acyclic, one output, two pins per
   gate, all 16 binary truth tables, free sharing and depth, exact on every
   input, one gate count per retained gate, inputs/wires/output selection free).
   The circuit quantified over in (5) is such a circuit. No transfer from a
   restricted model is involved, so the rules' demand for a proved transfer does
   not arise.
2. **`C` computes an affine disperser for sublinear dimension.** By premise E1
   (`disperser.md`), `f_n` is nonconstant on every affine subspace of dimension
   `≥ d(n)`, and `d(n) = 6n^{4/5} + O(1) = o(n)`. So `f_n` is an affine disperser
   for the sublinear dimension `d(n)`, which is the hypothesis of Lemma 3.7,
   Corollary 3.8, Theorem 3.9 and Theorem 4.1.
3. **`n` is the number of circuit inputs.** Yes: `B2` circuits for `f_n` receive
   exactly the `n` bits of `x`; witnesses are not inputs (`proof.md` §3).
4. **`d ≤ n`, so the potential/floors are meaningful.** For all sufficiently
   large `n`, `d(n) = ⌊6n^{4/5}⌋ + 1 ≤ n`. Theorem 4.1's right side
   `δ(n − 2d − 2)` is then nonnegative; when it is negative the inequality
   `μ ≥ δ(n−2d−2)` still holds trivially (`μ ≥ 0`), so no restriction on `n` is
   hidden in (5).
5. **Constants are admissible.** `α_Q = 9/5 ≤ 48/5 = α_I`, `0 ≤ α_φ = 1/5 ≤ 1/2`,
   all nonnegative, `δ = 64/5 ≥ 0`.
6. **The `o(n)` is explicit.** Theorem 3.9's `o(n)` and Theorem 4.1's `−2d − 2`
   are both absorbed into the displayed remainder `(261/10)d(n) + 128/5`, which
   is `O(n^{4/5})` by premise E1. No asymptotic is left unquantified.

## 5. The recorded source discrepancy, and why the conservative branch is used

Li–Yang Corollary 3.8 bounds the potential by `n/2 + 5d/2`, whereas the displayed
conclusion of Theorem 3.9 uses `n/2 + d`. The two finite expressions therefore
differ. This submission consistently uses the **larger** loss `5d/2`, which is
the weaker (conservative) hypothesis, and thus yields the larger remainder
`(261/10)d + 128/5`. If the smaller displayed `d` is used instead, the same
derivation gives

    g ≥ (31/10)·n − (128/5 + 1/5)·d − 128/5 = (31/10)·n − (129/5)·d − 128/5,

which is still `o(n)`. Either way the coefficient is exactly `31/10` and the
exponent `1`; the declared pair is unaffected, and no scientific difference is
claimed between the two variants. Nothing in the submission depends on resolving
which of the two displays is the intended one.

## 6. Normalization, degenerate gates and the model's bookkeeping

The rules require that a normalized proof cover constant, unary and
terminal-output cases without adding an uncharged computational basis. This is
automatic here, because (5) is a lower bound on the gate count of *every* `B2`
circuit computing `f_n`, with each retained gate counted once including trivial
and degenerate gates:

* Any circuit can be simplified (constant folding, absorbing negations into
  descendant truth tables, deleting non-ancestors of the output) to a normalized
  circuit with **no more** gates computing the same function. A lower bound for
  all circuits is therefore already a lower bound for the normalized ones and
  vice versa; no case split is needed for the inequality.
* The same holds for constant nodes introduced during substitutions, repeated
  pins and shared subgraphs: those conventions are internal to the cited measure
  argument (Li–Yang Lemma 3.11 and the substitutions of Section 4), and this
  submission does not re-run them, so it neither adds nor removes an uncharged
  basis.
* The terminal-output and unary cases cannot occur for the functions bounded
  here: `f_n` is nonconstant for all sufficiently large `n` (it is a disperser
  for dimension `d(n) ≤ n`), so no circuit for it is a single input or a
  constant, and no output-selection or unary-gate convention enters (5).
* Output selection, wires and inputs cost zero in `B2`, and `g` in Definition
  3.6 counts gates only, so the accounting matches `rules.md` exactly.

## 7. Non-claims

* Theorem 4.1 (and with it the long case analysis of troubled-pair elimination)
  and Theorem 1.1 are used as published theorems, cited by version, theorem and
  locator, with the hypotheses checked in §4. Their internal case analyses are
  not reproduced; the exact constants that enter the final inequality are
  reproduced and checked in `lp-arithmetic.md`.
* No claim is made that `31/10` is optimal among all possible techniques; the LP
  optimality in `lp-arithmetic.md` concerns only the *unchanged* Li–Yang LP.
* No endpoint, no `p > 1`, no logarithmic gain and no coefficient above `31/10`
  is asserted anywhere in this submission.
