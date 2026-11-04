# The fixed family: a subspace-polynomial affine disperser

This file defines the single family used by the bound, proves that one fixed
polynomial-time algorithm evaluates it at every length, and records the
disperser property with its hypotheses checked. The lower-bound argument that
consumes this property is in `gate-elimination.md`; the assembly is in
`proof.md`.

## 1. What the lower-bound argument needs from the family

Write B2 for Li--Yang's general model: deterministic acyclic circuits with one
output, fan-in two, all sixteen binary Boolean operations available as gate
types, unrestricted depth and fan-out, and exact computation on every input.
`C_B2(g)` denotes the number of gates of a smallest B2 circuit computing `g`.

**Definition 1.1 (affine disperser).** Let `d` be an integer. A function
`f : F_2^n -> F_2` is an *affine disperser for dimension `d`* if `f` is not
constant on any affine subspace `A` of `F_2^n` with `dim A >= d`.

**Theorem 1.2 (from `gate-elimination.md`, `Theorem 10.1`).** Let `f` be an
affine disperser for dimension `d`, where `0 <= d < n`. Then
`C_B2(f) >= 4n - 10d`.

Consequently, to obtain a lower bound `4n - o(n)` for one fixed uniformly
defined total family it suffices to exhibit one family `(f_n)` with

1. one fixed polynomial-time evaluator that computes `f_n` on every input of
   length `n` (this gives membership in P, hence in NP with an empty witness),
   defined for every `n >= 0`, and
2. a function `d(n)` with `d(n) < n` for all large `n`, `d(n)/n -> 0`, such
   that `f_n` is an affine disperser for dimension `d(n)`.

The family below is the subspace-polynomial construction of Ben-Sasson and
Kopparty, whose Theorem 2.2 supplies condition 2 with
`d(n) = floor(6 n^{4/5}) + 1`. Condition 1 is proved here in full.

## 2. The construction, at every length

Fix once and for all an integer `n >= 1`. All symbols below are functions of
`n`; nothing is advice, and no step needs information about `n` beyond `n`.

1. **Field degree.** Let `k = k(n)` be the least prime strictly larger than
   `2 n^{3/5}`.
2. **Block count.** Let `r = r(n) = ceil(n / k)` and `t = t(n) = ceil(sqrt(r))`.
3. **Field.** Let `A` be one fixed deterministic algorithm which, on input `k`, outputs a
   monic irreducible polynomial `q_k(Z) in F_2[Z]` of degree `k` and runs in time polynomial
   in `k`. Deterministic polynomial-time constructions of monic irreducible polynomials of a
   prescribed degree over a finite field are classical (Adleman and Lenstra, 1986; Shoup,
   1989, Section 3), and the algorithm is fixed once and for all, including every internal
   deterministic choice it makes. **The family is defined by `A`**: `q_k := A(k)`, and
   `K := F_2[Z]/(q_k)` with the `F_2`-basis `(1, Z, ..., Z^{k-1})`.

   Two points, because they were challenged before. First, no search over all degree-`k`
   polynomials and no lexicographic or other tie-break is used; the definition names the
   single polynomial `A(k)`, and the evaluator of Lemma 3.1 runs `A` itself, so definition
   and evaluation construct the *same* field and the *same* function. Second, the lower
   bound does not depend on which `q_k` is used: `disperser.md` proves the disperser
   property for every monic irreducible `q_k` of degree `k`, and the gate-elimination
   theorem of `gate-elimination.md` applies to any function. What depends on `A` is only
   which single family the declaration names, and that is fixed by the algorithm.

4. **Injective linear embedding.** Define `iota : F_2^n -> K^r` as follows.
   Append `r*k - n >= 0` zero coordinates to `x in F_2^n`, cut the resulting
   word into `r` consecutive blocks of `k` bits, and read the `i`-th block
   `(b_{i,0}, ..., b_{i,k-1})` as the field element
   `x_i = sum_j b_{i,j} Z^j`. This is `F_2`-linear and injective, and
   `iota(x) = (x_1, ..., x_r)`.
5. **Output map.** Let `pi : K -> F_2` be the coefficient of `Z^0`, i.e.
   `pi(z) = z_0` for `z = sum_j z_j Z^j`. It is `F_2`-linear, `pi(1) = 1`, so it
   is nonconstant.
6. **Powers.** For `i = 1, ..., r` put `y_i = x_i^{2^{i+1} - 1} in K`, where
   `2^{i+1} - 1 = 1 + 2 + 4 + ... + 2^i`. These exponents are strictly
   increasing in `i`.
7. **Symmetric polynomial.**
   `Sym_r^t(y_1, ..., y_r) = sum_{I subseteq [r], |I| = t} prod_{i in I} y_i`,
   the sum of the products over all `t`-element subsets, computed in `K`.
8. **The function.** `f_n(x) = pi( Sym_r^t(y_1, ..., y_r) ) in F_2`, and
   `f_0` (the empty word) `= 0`.

This is exactly equation (4) of Ben-Sasson--Kopparty with `p = 2`, `m = n`,
the field degree written `k` instead of `n`, `phi = iota` (an injective
`F_2`-linear map, which is what their Theorem 2.2 requires), `pi` a nontrivial
`F_2`-linear functional, and the exponents `1 + p + ... + p^i = 2^{i+1} - 1`.
No divisibility of `n` by `k` is assumed: the injective map absorbs the
`r*k - n` leftover bits as zero padding.

**The language and its characteristic functions.** Define
`L = { x in {0,1}^* : f_{|x|}(x) = 1 }`. The characteristic function of `L` on
`{0,1}^n` is `f_n`. One family, defined at every length including `n = 0`,
with no length-dependent choice of algorithm.

## 3. One fixed polynomial-time evaluator (so `L` is in P)

**Lemma 3.1.** There is one deterministic algorithm `A` which, on input
`x in {0,1}^n` with `n >= 1`, computes `f_n(x)`, and whose running time is
polynomial in `n`. The algorithm is the same for all `n`.

*Proof.* The algorithm performs steps (1)--(8) of Section 2 literally.

*Locating `k`.* Candidate odd integers are tested in increasing order starting
from `floor(2 n^{3/5}) + 1`; primality of a candidate `q <= 4 n^{3/5}` is
decided by trial division by all integers up to `sqrt(q) <= 2 n^{3/10}`. By
Bertrand's postulate there is a prime in `(2 n^{3/5}, 4 n^{3/5})`, so the scan
halts with `k <= 4 n^{3/5}`; it examines `O(n^{3/5})` candidates at
`O(n^{3/10})` divisions each, i.e. `O(n^{9/10})` word operations, which is
polynomial in `n`.

*Constructing the field.* The algorithm `A` of Section 2 item 3 is deterministic and
runs in time polynomial in `k`, hence in `n` because `k <= 4 n^{3/5}`; the evaluator calls
that same `A`, so it works in the field `K` named by the declaration. Addition in `K` is
bitwise XOR of `k`-bit vectors; multiplication is schoolbook multiplication of
polynomials of degree `< k` followed by reduction modulo the fixed `q_k`, hence
`O(k^2)` bit operations.

*Computing the powers.* `y_i` is obtained by `i + 1` successive squarings of
`x_i` starting from `x_i^{2^1 - 1} = x_i`; squaring is one multiplication. The
total is `sum_{i=1}^{r} (i+1) = O(r^2)` multiplications, each `O(k^2)` bit
operations, so `O(r^2 k^2)` bit operations; since `r*k <= n + k <= 2n` for
`n >= 2`, this is `O(n^2)` bit operations, polynomial in `n`.

*Computing the elementary symmetric polynomial without division.* Characteristic
two has no factor `1/j`, so use the addition-only recurrence. Initialise
`E_0 = 1` and `E_j = 0` for `1 <= j <= t`. On reading `y_i`, update
`E_j <- E_j + y_i * E_{j-1}` for `j` descending from `min(i, t)` down to `1`.
Descending order is essential: it keeps `E_{j-1}` at its value for
`y_1, ..., y_{i-1}` while `E_j` is updated. After all `r` values are read,
`E_j` equals the `j`-th elementary symmetric polynomial of `y_1, ..., y_r`.
Indeed, the `t`-element subsets of `[i]` are those contained in `[i-1]` (counted
by the old `E_t`) together with those containing `i`, which are
`y_i * (t-1)`-subsets of `[i-1]` (counted by the old `E_{t-1}`); summing these
two contributions is the update, and it is exactly the recursion above. The cost
is `O(r t)` field multiplications, `O(r t k^2)` bit operations, polynomial in
`n`.

*Output.* `pi(E_t)` is a single coordinate. Hence `A` is one deterministic
polynomial-time algorithm computing `f_n` at every length. `QED`

**Corollary 3.2 (`L` in NP, with a fixed verifier and witness bound).** Let
`V(x, w) = 1` iff `w` is the empty string and `A(x) = 1`, and let `W(N) = 0`.
Then `V` runs in time polynomial in `|x|`, and for every `x`:
`x in L` iff `f_{|x|}(x) = 1` iff `A(x) = 1` iff there exists `w` with
`|w| <= W(|x|) = 0` and `V(x, w) = 1`. Both directions are immediate, because
`w` is forced to be the empty word. The witness is not an input to the circuits
below, and no witness bit is charged to the circuit.

## 4. The disperser property, with hypotheses checked

**Theorem 4.1 (Ben-Sasson--Kopparty, Theorem 2.2 specialised to `p = 2`).**
Let `n >= 1`, let `k` be the least prime strictly greater than `2 n^{3/5}`, let
`r = ceil(n/k)` and `t = ceil(sqrt(r))`, let `iota : F_2^n -> K^r` be injective
and `F_2`-linear, and let `pi : K -> F_2` be nontrivial and `F_2`-linear. Then
the function of item 8 above is an affine disperser for dimension greater than
`6 n^{4/5}`: for every affine subspace `A subseteq F_2^n` with
`dim A > 6 n^{4/5}`, the restriction `f_n|_A` is nonconstant.

The statement in the source is "for all affine `A subseteq F_p^m` with
`dim(A) > 6 m^{4/5}` we have `|f(A)| > 1`", with `m = n` and `p = 2`. The
hypotheses are exactly: `k` the smallest prime bigger than `2 n^{3/5}`; `r`,
`t` as above; `phi` an injective `F_p`-linear map; `pi` a nontrivial
`F_p`-linear map. Section 2 satisfies each of them; in particular the injective
linear map exists for every `n` because `r*k >= n`.

**Corollary 4.2 (integer form used below).** Put
`d(n) = floor(6 n^{4/5}) + 1`. Then `d(n) > 6 n^{4/5}`, so for every affine
subspace `A subseteq F_2^n` with `dim A >= d(n)` we have `dim A > 6 n^{4/5}`
and `f_n|_A` is nonconstant. Hence `f_n` is an affine disperser for dimension
`d(n)`, in the sense of Definition 1.1.

## 5. Parameter arithmetic, checked

**Lemma 5.1.** For every integer `n >= 1`, with `k`, `r`, `t` as in Section 2:

1. `2 n^{3/5} < k < 4 n^{3/5}`;
2. `r = ceil(n/k)`, so `r <= n/k + 1` and `r >= n/k`;
3. `t = ceil(sqrt(r))`, so `sqrt(r) <= t <= sqrt(r) + 1`;
4. `r >= 2` and `1 <= t <= r` hold for all `n >= 64`, and `t < r` for
   `r >= 3`; the finitely many smaller lengths are covered separately below and
   play no role in the limit;
5. `d(n) < n` for all `n >= 16807 = 7^5`, and `d(n)/n <= 6 n^{-1/5} + 1/n`.

*Proof.* (1) is Bertrand's postulate applied to `2 n^{3/5} >= 2`. (2) and (3)
are the definitions plus the elementary bounds on the ceiling. For (4):
`k < 4 n^{3/5}` gives `n/k > n^{2/5}/4 >= 2` once `n^{2/5} >= 8`, i.e.
`n >= 2^{10} = 1024`; then `r >= 2`; `t = ceil(sqrt(r)) <= r` for `r >= 2`, and
`t = ceil(sqrt r) < r` for `r >= 3`. For (5): if `n >= 7^5` then
`n^{1/5} >= 7`, so `n^{4/5} = n / n^{1/5} <= n/7`; hence
`6 n^{4/5} <= 6n/7` and
`d(n) = floor(6 n^{4/5}) + 1 <= 6n/7 + 1 <= n - 2401 + 1 < n`, because
`n/7 >= 2401`. (The exact first index with `d(n) < n` is `7781`, where
`d = 7780`; `checks.py` computes it, but nothing needs it: `16807` is the clean
closed-form threshold.) The ratio bound is immediate from
`d(n) <= 6 n^{4/5} + 1`. `QED`

The finite range `1 <= n <= 7776` matters only for totality of the family, not
for the asymptotic statement: `f_n` is defined there by Section 2 (with
`r >= 1`, `t >= 1`), and any bound is allowed at finitely many lengths.

## 6. The disperser property, and where its proof is supplied

Theorem 4.1 is proved in full in `disperser.md`, which supplies the argument rather than a
pointer to it:

- Part I: polynomial representation over `K` (reduction modulo `X^{2^k} − X`, trace form
  of a linear functional, `F_2`-degree calculus);
- Part II: subspace polynomials — existence of the image-subspace polynomial with
  `P_L(Q_A(X)) ≡ Q_A(P_L(X)) ≡ X^{2^k} − X` (II.1), the `2^i`-th power identity (II.3),
  and the coefficient-structure lemma II.4 in both its forms, including the use of
  "not contained in a coset of a proper subfield";
- Part III: the associated-matrix formula for the coefficient of a monomial in
  `Tr(µ Π Q_i^{α_i})` (III.3) and the triangularity-plus-zero-column criterion for that
  coefficient to be nonzero (III.4);
- Part IV: the independent-sources theorem IV.1 with its two structural claims IV.2--IV.3,
  including the interval covering that consumes `k ≥ r(r+1)/t` and primality of `k`;
- Part V: Bourgain's block decomposition of an affine subspace of `K^r`, with proof;
- Part VI: invariance of the relevant monomial coefficients under the affine
  perturbation that turns independent blocks into an arbitrary affine subspace, with
  proof;
- Part VII: the non-parameterised theorem VII.1 with the uniqueness claim VII.2, i.e. the
  statement that only one of the `t`-fold summands can contribute the distinguished
  monomial;
- Part VIII: the parameter passage to the threshold `6n^{4/5}` with the constants
  `k ≤ 4n^{3/5}`, `r = ceil(n/k)`, `t = ceil(sqrt r)`, whose arithmetic is checked by
  exact integer arithmetic in `checks.py`.

`disperser.md` Part IX states the two boundaries of that reproduction: the choice of the
index set `I_0` and its inequalities are re-derived there (the source's own bookkeeping is
not legible in the available text), and the interval covering of IV.3 is transcribed into
prose. Nothing else is owed to a citation; the source is named for attribution only.

## 7. Bit accounting

`n = |x|` is the actual input bit length. The construction uses `n` bits, the
`r*k - n` padding coordinates are internal zeros, and every block is exactly
`k` bits, so `K` is a field over `F_2` and `pi(y)` is one bit: no fractional
bit, block or field element is charged to the circuit. The lower bound
`4n - 10 d(n)` is stated against the same `n` that indexes `f_n`, and the loss
`10 d(n) = 10 floor(6 n^{4/5}) + 10 = O(n^{4/5})` is `o(n)`. Hence
`liminf_n C_B2(f_n)/n >= 4`, which is the declared projection.
