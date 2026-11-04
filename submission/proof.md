# A `4n - o(n)` lower bound for one fixed total uniform NP family in the full B2 model

**Claim declared in `bound.json`.** There is one fixed total Boolean family
`(f_n)` in NP such that

    liminf_{n -> infinity} C_B2(f_n) / n  >=  4,

i.e. coefficient `4` at exponent `1`, where `C_B2` counts gates in Li--Yang's
general B2 model at the actual input length `n = |x|`, for every sufficiently
large integer length.

The submission has three parts:

1. `family.md` — one family `(f_n)`, defined at every length by a fixed
   deterministic rule, with one fixed polynomial-time evaluator (so the language
   is in P, hence in NP with an empty witness), and with the published affine
   disperser property of the subspace-polynomial construction, at dimension
   `d(n) = floor(6 n^{4/5}) + 1 = o(n)`.
2. `gate-elimination.md` — `Theorem 10.1`: every B2 circuit computing an affine
   disperser of dimension `d` on `n` inputs has at least `4n - 10 d` gates. This
   is proved there from its own definitions, including the complete five-case
   gate-elimination analysis, the potential, the subspace lemma and the
   termination argument.
3. This file — the assembly: the theorem is applied to the family at every
   length, the `o(n)` loss is evaluated, and the limit and quantifiers are
   written out.

Together: `C_B2(f_n) >= 4n - 10 d(n) = 4n - O(n^{4/5})`, hence
`liminf C_B2(f_n)/n >= 4`.

The previously published best bound for this model is `3.1n - o(n)` (the seed in
the repository, from Li--Yang, STOC 2022, applied to the same
subspace-polynomial family). The improvement from `3.1` to `4` comes from the
gate-elimination theorem of part 2, not from the family: the family here is the
same one the seed uses.

## 1. The family (summary; full detail in `family.md`)

For `n >= 1`: let `k` be the least prime strictly larger than `2 n^{3/5}`, let
`r = ceil(n/k)` and `t = ceil(sqrt(r))`. Fix one deterministic polynomial-time
procedure producing a monic irreducible `q_k` of degree `k` over `F_2` (Shoup
1989, Section 3, Theorem 3.2, with a fixed tie-breaking rule), and put
`K = F_2[Z]/(q_k)`. Let `iota : F_2^n -> K^r` be the injective `F_2`-linear map
that pads `x` with `r*k - n` zeros, cuts it into `r` blocks of `k` bits, and
reads each block in the basis `(1, Z, ..., Z^{k-1})`. Let `pi : K -> F_2` be the
coefficient of `Z^0`. Set `y_i = iota_i(x)^{2^{i+1}-1}` and

    f_n(x) = pi( sum_{I subseteq [r], |I| = t} prod_{i in I} y_i ),

the sum of the products over `t`-element subsets (the `t`-th elementary
symmetric polynomial in `y_1, ..., y_r`, over `K`). Put `f_0(empty) = 0`, and let
`L = {x in {0,1}^* : f_{|x|}(x) = 1}`.

This is equation (4) of Ben-Sasson--Kopparty with `p = 2` and their field degree
written `k`, read through an injective `F_2`-linear map and a nontrivial
`F_2`-linear output functional, exactly as in their Theorem 2.2 (their
Proposition 3.3 relates that parameterised form to the non-parameterised form of
their Theorem 7.3, up to renaming variables).

## 2. Eligibility: one fixed polynomial-time verifier, both directions

`family.md` Lemma 3.1 gives a single deterministic algorithm `A`, the same at all
lengths, that computes `f_n(x)` in time polynomial in `n`: `k` is found by a scan
with trial-division primality tests (Bertrand gives `k <= 4 n^{3/5}`), the
irreducible `q_k` comes from the fixed deterministic procedure, the powers use
`O(r^2)` field multiplications by repeated squaring, and the elementary symmetric
polynomial is computed by the addition-only recurrence

    E_0 = 1,  E_j = 0 (j >= 1),  then for i = 1..r:  E_j <- E_j + y_i E_{j-1},
    for j descending from min(i,t) to 1,

which needs no division by `j` (characteristic two has no `1/j`), costs `O(r t)`
field operations, and is correct because the `t`-subsets of `[i]` split into
those avoiding `i` (old `E_t`) and those containing `i` (`y_i` times a
`(t-1)`-subset of `[i-1]`, old `E_{t-1}`).

Take the verifier `V(x,w) = 1` iff `w` is the empty word and `A(x) = 1`, with
witness bound `W(N) = 0`. Then `V` is polynomial time and, for every `x`,

    x in L  <=>  f_{|x|}(x) = 1  <=>  A(x) = 1
            <=>  there exists w, |w| <= 0, with V(x,w) = 1.

Both directions are immediate, since the only admissible witness is the empty
word. Hence `L` is in P, and therefore in NP. This is a language-level
statement: one language, one verifier, one witness bound, for all lengths. The
witness is not an input to the circuits; the circuits compute `f_n` on `x` alone.

## 3. The model and the absence of a transfer

`gate-elimination.md` Section 1 fixes the model used there to be exactly the
challenge model: deterministic acyclic circuits, one output, fan-in two, all
sixteen binary operations, unbounded depth and fan-out, exact computation on
every input, gates counted once with inputs, wires and output selection free,
and negations absorbed into the operations of the reading gates.

`Theorem 10.1` of `gate-elimination.md` is stated and proved for B2 circuits
directly. Its proof starts from an arbitrary B2 circuit `C_0` and uses fair
semi-circuits only as an analysis device whose *initial state is that circuit*
(empty affine system, acyclic part `C_0`). A semi-circuit is a *superset* of the
B2 model in the sense that everything the elimination produces is simulated by
the invariant `C(u) = C_0(Phi(u))` on `Z(P)`; the argument never assumes the
input circuit has fewer gates, restricted fan-out, restricted depth, or
restricted gate types. Therefore no restricted-model transfer (monotone,
formula, bounded depth, XOR-free, multiple output, randomized, average-case) is
invoked anywhere: the bound is proved in the full model, for exact computation
on all inputs.

## 4. Applying the theorem at every length

By `family.md` Corollary 4.2, `f_n` is an affine disperser for dimension
`d(n) = floor(6 n^{4/5}) + 1`, for every `n >= 1`. By `family.md` Lemma 5.1(5),
`d(n) < n` for every `n >= 16807 = 7^5` (indeed `n^{4/5} <= n/7` there, so
`6 n^{4/5} + 1 < n`; the exact first index `7781` is computed in `checks.py`).
So for every integer `n >= 16807`, `Theorem 10.1` applies with
`d = d(n)` and gives

    C_B2(f_n) >= 4n - 10 d(n) = 4n - 10 floor(6 n^{4/5}) - 10 >= 4n - 60 n^{4/5} - 10.

Lengths `0 <= n < 16807` are irrelevant to a limit and are covered by the same
definition of the family; no statement is made there beyond the bound being
allowed to be vacuous or negative.

**Theorem (main).** With `(f_n)` and `L` as in Sections 1--2,

    for every epsilon > 0 there is N_epsilon such that, for every integer
    n >= N_epsilon,   C_B2(f_n) >= (4 - epsilon) n.

*Proof.* Fix `epsilon in (0, 1]`. Choose `N_epsilon >= 16807` with
`60 n^{-1/5} + 10/n <= epsilon` for all `n >= N_epsilon`, which exists because
`n^{-1/5} -> 0` and `1/n -> 0`. Then for every integer `n >= N_epsilon`,

    C_B2(f_n) >= 4n - 60 n^{4/5} - 10 = (4 - 60 n^{-1/5} - 10/n) n >= (4 - epsilon) n.

`QED`

**Corollary.** `liminf_n C_B2(f_n)/n >= 4`, and `liminf_n C_B2(f_n)/n^p >= 4`
at `p = 1`. The declared finite projection is therefore `(coefficient, exponent)
= (4, 1)`, with the full guarantee expressed as the `polynomial_log` expression
`4 n^1 (log_2 n)^0 (log_2 log_2 n)^0`: `liminf C_B2(f_n)/g(n) >= 1` for
`g(n) = 4n`, which is the same statement.

## 5. Bit accounting and the quantifiers, explicitly

- `n = |x|` is the actual input bit length; `f_n : F_2^n -> F_2`.
- The blocks are exactly `k` bits and `K` is a field over `F_2`, so a field
  element is `k` bits and `pi` returns one bit. The padding `r*k - n` is inside
  the definition of `iota`, not an extra circuit input; the circuits are on `n`
  inputs and are charged against the same `n`.
- `k <= 4 n^{3/5}`, `r = ceil(n/k) <= n/k + 1`, `t = ceil(sqrt r)`, and
  `d(n) = floor(6 n^{4/5}) + 1`, so the additive loss is
  `10 d(n) = 60 n^{4/5} + O(1) = o(n)`; `d(n)/n -> 0` is verified in
  `family.md` Lemma 5.1(5).
- The quantifier structure is exactly the one the rules require: one family, one
  algorithm, one verifier, one witness bound, fixed before `epsilon` and `n` are
  quantified. No length is skipped: the bound is proved for every integer
  `n >= 16807`, not along a sparse subsequence, and not in the limsup sense.
- No step depends on a future choice of `epsilon`: `N_epsilon` is derived from
  the same explicit `O(n^{4/5})` remainder at every length.

## 6. Why this is stronger than the seed, and what it does not claim

The seed (`3.1n - o(n)`) and this bound use the same family. The difference is
the gate-elimination theorem: the seed's engine is the Li--Yang potential LP
whose optimum is exactly `31/10` (the repo's `lp-certificate` checker certifies
the corresponding feasible point and the dual bound `31/10`, and the seed's own
text notes that "retuning this LP alone cannot advance the record"). The
engine of Section 4 above is a different analysis, in which substitutions act on
a space of *parameters* subject to affine reparameterisation, so that two inputs
of a minimal AND-gate can be *paired* and the gate's value fixed without any
substitution being performed. That removes the alternative analysis and yields
the potential bound `dW <= 10 da`, hence `4n - 10d`; `proof.md` makes no use of
the Li--Yang LP.

What is not claimed:

- Not claimed: any improvement of the *family*; not claimed: `p > 1`, or any
  exponent beyond `1`. The bound is linear, so the declared exponent is exactly
  `1`.
- Not claimed: `4n` with a smaller loss. `10 d(n) = O(n^{4/5})` is what this
  route gives with the Ben-Sasson--Kopparty family; Li's `d = O(log n)`
  construction would give `4n - O(log n)`, still coefficient `4`.
- Not claimed: the constructive form. The theorem is a refutation argument;
  finding the counterexample input explicitly costs `poly(n, 2^d)`, which is
  exponential at `d = O(n^{4/5})`. Only the lower bound is asserted.
- Not claimed: that any part of the argument is new mathematics. Both halves are
  reproductions: the disperser property is proved in full in `disperser.md`
  (Parts I--VIII: subspace-polynomial structure, the monomial-coefficient
  criterion, the independent-sources theorem, Bourgain's decomposition, monomial
  invariance, and the assembly), and the conversion from "affine disperser of
  dimension `d`" to "`4n - 10d` gates" is proved in full in
  `gate-elimination.md`. `disperser.md` Part IX states the two places where it
  re-derives bookkeeping rather than transcribing it.
- Not claimed: that `4` is the best possible coefficient for B2 circuits, or
  that this family cannot have a better bound.

## 7. Sources

- **MS26** — Om Swostik Mishra, Rahul Santhanam, *An Explicit `4n - O(log n)`
  Size Lower Bound against General Boolean Circuits*, ECCC TR26-260,
  revision 1, 3 November 2026. Used for the gate-elimination theorem
  (`Theorem 1.1`, `Theorem 4.1`, `Lemma 4.2`, `Remark 4.3`).
  URL `https://eccc.weizmann.ac.il/report/2026/260/revision/1/download`,
  SHA256 `a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921`.
  The preceding version of the same report (TR26-260, 19 October 2026) proves
  `3.5n - o(n)` and is superseded by revision 1.
- **BK** — Eli Ben-Sasson, Swastik Kopparty, *Affine Dispersers from Subspace
  Polynomials*, author preprint dated 2010-03-12, `Theorem 2.2`, equation (4),
  printed pp. 6--7; journal version SIAM J. Comput. 41(4):880--914, 2012. Used
  for the disperser property of the family.
  URL `https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf`,
  SHA256 `922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf`.
- **Shoup** — Victor Shoup, *New algorithms for finding irreducible polynomials
  over finite fields*, author preprint dated 1989-01-31, Section 3,
  Theorem 3.2. Used only to make the choice of `q_k` canonical and
  polynomial-time.
  URL `https://www.shoup.net/papers/detirred.pdf`,
  SHA256 `00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644`.
- **FGHK**, **LY** — Find--Golovnev--Hirsch--Kulikov (FOCS 2016) and
  Li--Yang (STOC 2022, ECCC TR21-023) are the antecedents of MS26 and of the
  seed. They are cited for context; nothing in this submission depends on their
  statements. The seed's LP is not used here.

## 8. Checks performed, and what they establish

`checks.py` is exact-arithmetic verification of the parts that admit it:

1. the classification of the sixteen binary operations into 2 constant, 4
   unary, 2 XOR-type, 8 AND-type, and the identification of the 8 affine
   operations;
2. every row of the accounting table of `gate-elimination.md` Section 9,
   recomputed from its stated changes in `s`, `q`, `m` and in weights, checking
   `da = dm - dq` and `dW <= 10 da`; and that the rows needing `dW <= -10` are
   exactly the rows with `da = -1`;
3. the numerical part of `Theorem 10.1`: that
   `W_0 <= |C_0| + 6n < 10(n-d)` is equivalent to `|C_0| < 4n - 10d`, and that
   the summation `0 <= W_t <= W_0 + 10(a_t - n) < 10(a_t - d)` forces `a_t > d`,
   checked on an exact-rational grid;
4. the parameter arithmetic of `family.md` Lemma 5.1 for
   `64 <= n <= 10^5` (including `d(n) < n` from `16807`, `t < r` for `r >= 3`,
   `k` prime and `2n^{3/5} < k < 4n^{3/5}`), and the Ben-Sasson--Kopparty
   parameter hypothesis `k >= r(r+1)/t` together with
   `k r/t + k t + r(r+1) <= 6 n^{4/5} + o(n^{4/5})`;
5. the exact rationals used in `bound.json`: reduced `4/1` and `1/1`, and that
   the declared projection is implied by the declared guarantee under
   `expressions.py`.

These checks certify bookkeeping, not mathematics. In particular they do not
establish: that the five cases cover every configuration in which Cases 1--3
fail; the geometric claims in Case 2 and in Lemma 5.1; the preservation of the
invariant `(*)`; the disperser property; or the correctness of the case analysis
in `gate-elimination.md`. Those are proved in prose in
`gate-elimination.md` Sections 2--8 and must be read as an argument, not
verified as arithmetic. The repo's `lp-certificate` checker is not claimed and
not used: it certifies the seed's Li--Yang LP, which plays no role here.

## 9. Summary of the claim as declared

- Family: the single subspace-polynomial family of Section 1, at every length.
- Eligibility: in P by one fixed polynomial-time algorithm, hence in NP with the
  empty witness and witness bound `0`.
- Model: full B2, exact computation on all inputs, gates counted once.
- Guarantee: `C_B2(f_n) >= 4n - 60 n^{4/5} - 10` for every integer `n >= 16807`,
  hence `liminf C_B2(f_n)/n >= 4`.
- Declared finite projection: coefficient `4/1`, exponent `1/1`.
- Weaker, unscored statements that are also proved here and declared as string
  claims: (a) `C_B2(f_n) >= 4n - 10 d(n)` with
  `d(n) = floor(6 n^{4/5}) + 1`, for every integer `n >= 16807`; (b) the general
  theorem `C_B2(g) >= 4n - 10d` for every affine disperser `g` of dimension `d`
  on `n` inputs, `d < n`.
