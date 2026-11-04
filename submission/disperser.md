# Premise E1: the family is an affine disperser for dimension ⌊6n^(4/5)⌋+1

This file matches the family of `proof.md` Section 2 term by term against
Theorem 2.2 of Ben-Sasson–Kopparty, `Affine Dispersers from Subspace
Polynomials` (author preprint dated 2010-03-12,
<https://www.math.toronto.edu/swastik/affine-disperser.pdf>; journal version
SIAM J. Comput. **41**(4):880–914, 2012, DOI 10.1137/110826254).

## 1. The cited theorem, restated

**BK Theorem 2.2 (Affine dispersers for sublinear dimension).** *Given the
integer `m`, fix parameters `n, r, t` as follows. Let `n` be the smallest prime
bigger than `2·m^{3/5}`. Let `r = ⌈m/n⌉` and `t = ⌈√r⌉`. (So `n ≈ m^{3/5}`,
`r ≈ m^{2/5}`, `t ≈ m^{1/5}`.) Let `φ : F_p^m → (F_{p^n})^r` be an injective
`F_p`-linear map, where `φ(y) = (φ_1(y), …, φ_r(y))` and `φ_i(y) ∈ F_{p^n}`. Let
`π : F_{p^n} → F_p` be a nontrivial `F_p`-linear map. Then the function
`f : F_p^m → F_p` defined by*

    f(x) = π( Sym_t^r( (φ_1(x))^{1+p}, (φ_2(x))^{1+p+p^2}, …, (φ_r(x))^{1+p+…+p^r} ) )      (BK-4)

*with `Sym_t^r(Y_1,…,Y_r) = Σ_{I⊆[r], |I|=t} ∏_{i∈I} Y_i`, is an affine disperser
for dimension greater than `6·m^{4/5}`, i.e. for all affine `A ⊆ F_p^m` with
`dim(A) > 6·m^{4/5}` we have `|f(A)| > 1`.*

Here `p` is a prime and the field is `F_p`; the paper notes that this `f` is
computable in polynomial time in `p` and `m` (e.g. via Newton–Girard), and its
Definition 2.1 fixes the notion used: `f : F_p^m → F_p` is an affine disperser
for dimension `d` if `|f(A)| > 1` for every affine subspace `A ⊆ F_p^m` of
dimension at least `d`.

## 2. Setting `p = 2` and reading (BK-4) literally

With `p = 2`, the exponent attached to the `i`-th block is

    1 + p + … + p^i = 1 + 2 + … + 2^i = 2^{i+1} − 1,

so (BK-4) becomes `f(x) = π(Sym_t^r(y_1,…,y_r))` with `y_i = φ_i(x)^{2^{i+1}−1}`.
This is *exactly* the definition in `proof.md` Section 2, items 6–8: there is no
reindexing, no change of exponent and no extra factor. (In particular the first
block carries exponent `2^2−1 = 3 = 1+p`, the second `2^3−1 = 7 = 1+p+p^2`, and
so on.)

## 3. Parameter match

The ambiguity to remove is that the cited paper names its ambient dimension `m`
and its prime block length `n`, whereas this submission names the ambient
(actual input) length `n` and the prime block length `k`. With that renaming the
two parameter lists coincide:

| BK Theorem 2.2 | this submission (`proof.md` §2) | value used |
| --- | --- | --- |
| ambient dimension `m` | actual input length `n` | `n = \|x\|` (bits) |
| prime block length `n` | `k = k(n)` | least prime `> 2n^{3/5}` |
| `r = ⌈m/n⌉` | `r = ⌈n/k⌉` | `= O(n^{2/5})` |
| `t = ⌈√r⌉` | `t = ⌈√r⌉` | `= O(n^{1/5})` |
| field `F_{p^n}`, `p = 2` | `K = F_2[Z]/(q_k)` | `\|K\| = 2^k` |
| injective `F_p`-linear `φ` | block map of §2 item 4 | injective since `rk ≥ n` |
| nontrivial `F_p`-linear `π` | constant-coefficient functional | `π(1) = 1 ≠ 0` |
| `Sym_t^r` | `Sym_r^t` | same object, index order only |

## 4. Hypothesis checks

1. **`k` is prime and is the smallest prime above `2n^{3/5}`.** This is
   `proof.md` §2 item 1 verbatim: `k` is *defined* as the least prime with
   `k > 2n^{3/5}`. Primality is what the paper needs in order that `F_{2^k}`
   have no nontrivial subfields, so that the "not contained in an affine shift
   of a proper subfield" obstructions of its intermediate results cannot occur.
   For `k = 2` or `k` prime, `F_{2^k}` has exactly two subfields, `F_2` and
   itself; and `k > 2n^{3/5} ≥ 2` for `n ≥ 1`, so `k ≥ 3`. Existence of the least
   prime above `2n^{3/5}` is Bertrand's postulate (`2n^{3/5} < k < 4n^{3/5}`).
2. **`r = ⌈m/k⌉` with `m = n`.** `rk ≥ n` holds by definition of the ceiling,
   and `r ≥ 1`.
3. **`t = ⌈√r⌉`.** Verbatim.
4. **`φ : F_2^n → K^r` is an injective `F_2`-linear map.** The paper's statement
   allows *any* injective `F_2`-linear `φ`; the block map of §2 item 4 is such a
   map: it is `F_2`-linear (each coordinate of `φ_i` is a single input bit), and
   injective because `rk ≥ n` and the `n` input bits occupy the first `n`
   coordinate slots of the `r·k` slots, so they are recovered exactly from
   `(φ_1(x),…,φ_r(x))`. The paper's own reading of the construction ("partition
   the `m` coordinates into `r` blocks of `n` coordinates each, assuming `n`
   divides `m` by discarding a few field-elements if necessary") is the
   special case where the discarded tail is the zero-padding used here; using an
   injective map with `rk ≥ m` instead of discarding input coordinates is
   explicitly within the theorem's quantification, and it is also the variant
   that keeps the *actual* input length equal to `n` with no free extra
   coordinates.
5. **`π : K → F_2` is a nontrivial `F_2`-linear map.** `π(y)` is the coefficient
   of `Z^0` in the basis `1,Z,…,Z^{k−1}`; it is `F_2`-linear and `π(1) = 1`.
6. **`f` is Boolean-valued.** `π` maps into `F_2 = {0,1}`, so `f_n(x) ∈ {0,1}`.
7. **Totality and uniformity.** The theorem is applied for each `n ≥ 1` with the
   deterministically defined `k, r, t, q_k, φ, π`; nothing in the application
   depends on advice, on a random choice, or on a per-length search for a good
   construction. `q_k` is produced by the fixed deterministic algorithm `AS`
   (Shoup), so it is a function of `k` alone, hence of `n` alone.

## 5. Conclusion extracted

By Theorem 2.2 applied with `m = n` and `p = 2`, for every `n ≥ 1` the function
`f_n : F_2^n → F_2` of `proof.md` §2 satisfies

    f_n is nonconstant on every affine subspace A ⊆ F_2^n with dim(A) > 6·n^{4/5}.         (2)

Taking `d(n) = ⌊6·n^{4/5}⌋ + 1`, every affine subspace with
`dim(A) ≥ d(n)` has `dim(A) > 6n^{4/5}`, so `f_n` is nonconstant on it. This is
(2) of `proof.md` §4. Two elementary consequences used later:

* `d(n) ≤ n` for all sufficiently large `n` (as `6n^{4/5} + 1 ≤ n` eventually),
  so the whole space is among the subspaces covered and `f_n` is nonconstant;
* `d(n)/n ≤ 6n^{−1/5} + 1/n → 0`, i.e. `d(n) = o(n)`, so the disperser
  dimension is sublinear as premise E2 requires.

## 6. Computational cost of the deterministic instantiation

`f_n` must be computable by one fixed polynomial-time algorithm (`proof.md` §3).
The parts that are not literal bit operations are:

* **`q_k`.** Shoup's deterministic algorithm (V. Shoup, *New algorithms for
  finding irreducible polynomials over finite fields*, Math. Comp. **54** (1990),
  435–447) finds a monic irreducible polynomial of degree `k` over `F_2` in
  deterministic polynomial time in `k`, using a fixed enumeration of candidates
  so that its output is a well-defined function of `k`. Since `k = O(n^{3/5})`,
  this is polynomial in `n`.
* **Field arithmetic.** `K` is represented in the polynomial basis
  `1,Z,…,Z^{k−1}` over `F_2`. Addition is `O(k)` bit operations; schoolbook
  multiplication followed by reduction modulo the degree-`k` polynomial `q_k` is
  `O(k^2) = O(n^{6/5})` bit operations. No inversion and no discrete logarithm is
  needed anywhere.
* **Powering.** `y_i = φ_i(x)^{2^{i+1}−1}` costs `O(log(2^{i+1})) = O(i)` field
  multiplications by repeated squaring, and `i ≤ r = O(n^{2/5})`.
* **Symmetric polynomial.** The division-free recurrence of `proof.md` §3 step 5
  needs `O(r·t)` field multiplications and is valid in characteristic 2 because
  it never divides by an integer.

The total cost is polynomial in `n` (`proof.md` §3 gives the exponents), so `f_n`
lies in `P` through a single fixed algorithm, and the NP verifier with the empty
witness applies.

## 7. What is not used and not claimed

* Only **Theorem 2.2** of Ben-Sasson–Kopparty (and its Definition 2.1) is used.
  The paper's affine *extractors*, its results for independent sources, its
  univariate cubic/quartic constructions and its conjectures play no role.
* No claim is made that `f_n` is an extractor, that it is computable in a
  restricted circuit class, or that the disperser dimension `⌊6n^{4/5}⌋+1` is
  tight or optimal. Any smaller valid dimension bound would only weaken the
  remainder in `proof.md` (5), never the coefficient.
* No claim is made that the construction is the only one, or that this
  instantiation is the one used by Li–Yang. It is a published affine disperser
  for sublinear dimension, which is all premise E2 requires.
