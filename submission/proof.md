# A (4,1) lower bound: C_B2(f_N) >= 4N - 29 N^{4/5} for all N >= 2^20

**Theorem.** Let `(f_N)` be the single total uniform family defined in
Section 1. Then `f_N` is computable in polynomial time (hence its language
is in P, and in NP with empty witnesses), and for every integer
`N >= 2^20`, in the full Li–Yang B_2 model,

```
C_B2(f_N) >= 4N - R(N) >= 4N - 29 * N^{4/5},
```

where `R(N)` is the explicit loss of Section 4. In particular
`liminf_N C_B2(f_N) / N >= 4`. The proof combines a pair-domain gate
elimination with deletion charge 9 (`elimination.md`, Theorem E3) and the
Ben-Sasson–Kopparty affine disperser in its exact-threshold form
(`disperser.md`, Corollary D7).

## 1. The family

One deterministic evaluator `A` computes `f_N(x)` for every input length
`N = |x| >= 0`:

1. Let `B` be the least positive integer with `B^5 > 32 N^3`, and let `k`
   be the least prime `>= B`. (For `N = 0`: `B = 1`, `k = 2`.)
2. If `k >= N`, output `0`. Else set `m = N - k >= 1`.
3. Read the first `k` bits of `x` as the lower coefficients of
   `Q_x(Z) = Z^k + sum_{j<k} x_j Z^j` over `F_2` (a monic degree-`k`
   polynomial; the leading coefficient is the implicit `1`).
4. Test whether `Q_x` is irreducible over `F_2` by the exact test of
   Section 2. If it fails, output `0`.
5. Let `r = ceil(m/k)` and `t = ceil(sqrt(r))`. Split the remaining `m`
   bits into `r` blocks of `k` bits, padding with zeros; block `i` gives
   `a_i in K := F_2[Z]/(Q_x)`, first bit the constant coefficient. This is
   an injective `F_2`-linear map `phi : F_2^m -> K^r`.
6. Put `y_i = a_i^{2^{i+1} - 1}` (binary powering) and compute
   `E = Sym_r^t(y_1, ..., y_r)` by the division-free recurrence
   `E_0 = 1`, `E_j = 0` (`1 <= j <= t`); for each `i`, update
   `E_j <- E_j + y_i E_{j-1}` for `j = min(i, t)` down to `1`. Output the
   constant coefficient of `E_t`.

Every step is polynomial in `N`: `B` and `k` by scanning with trial
division (`k = O(N^{3/5})`); field arithmetic in the polynomial basis;
`O(i)` field multiplications for the power in step 6 (the exponent
`2^{i+1} - 1` has `i + 1` bits); `O(r t)` field operations for the
recurrence, which is correct by splitting the `t`-subsets according to
whether they contain `i` (no division, so it is valid in characteristic
2). Hence one polynomial-time evaluator computes the total family; with
`V(x, w) = 1` iff `w` is empty and `A(x) = 1` and witness bound `0`, both
directions of the NP condition are immediate, so the language is in P
(hence NP), with no per-length advice: the same algorithm and the same
formulas define every length, including the degenerate ones (step 2).

## 2. The irreducibility test is exact

**Lemma.** For prime `k >= 2` and monic `Q` of degree `k` over `F_2`, `Q`
is irreducible iff `Z^{2^k} ≡ Z (mod Q)` and `gcd(Q, Z^2 + Z) = 1`. Both
conditions are checkable in time polynomial in `k` (`k` repeated squarings
modulo `Q`; one polynomial gcd).

*Proof.* If `Q` is irreducible, `F_2[Z]/(Q) = F_{2^k}`, every element
satisfies `a^{2^k} = a`, giving the congruence; and `Q` is neither `Z` nor
`Z + 1` (degree `k >= 2`), so the gcd is `1`. Conversely, the congruence
says `Q | Z^{2^k} - Z`. The polynomial `Z^{2^k} - Z` is squarefree (its
derivative is `1`), so `Q` is squarefree and every irreducible factor `P`
of `Q` divides `Z^{2^k} - Z`; the roots of `P` then lie in `F_{2^k}`, so
`deg P` divides `k` (the field they generate is `F_{2^{deg P}} ⊆
F_{2^k}`). Since `k` is prime, `deg P in {1, k}`; the gcd condition
excludes the only monic degree-1 irreducibles `Z` and `Z + 1`. Hence
every factor has degree `k`; as `Q` is monic of degree `k` and
squarefree, `Q` is a single irreducible factor. ∎

For each `N` with `k < N` there is at least one monic irreducible of
degree `k` over `F_2`: their number is
`(1/k) * sum_{d | k} mu(d) 2^{k/d} = (2^k - 2)/k >= 1` for prime `k`
(standard Möbius count of monic irreducibles). Fix `Q_N` to be the
lexicographically least one (as a `k`-bit coefficient string). This is a
definite mathematical choice used only in the lower-bound argument; the
evaluator never computes `Q_N`.

## 3. The hard restriction and the disperser step

Fix `N >= 2^20` and a B_2 circuit `C` computing `f_N`. Restrict the first
`k` input pins of `C` to the coefficient bits of `Q_N`. Substituting
constants and repeatedly deleting gates that have become constant or
unary (rerouting their uses) never increases the gate count — it is the
deletion process (D1)/(D2) of `elimination.md` — so the restricted
function `g(z) = f_N(cert_N || z)` on `m = N - k` bits has a circuit of
size at most `|C|`. By Section 2 the evaluator's test accepts `Q_N`, so
`g` is exactly the function of Corollary D7 of `disperser.md` for the
field `K = F_2[Z]/(Q_N)`: `g(z) = pi(Sym_r^t(y_1, ..., y_r))`,
`y_i = phi_i(z)^{2^{i+1} - 1}`, `pi` the constant-coefficient functional.

The hypotheses of Corollary D7 hold for all `N >= 2^20`: with
`v = N^{1/5} >= 16` one has (Section 4) `2v^3 < k <= 2.0140 v^3`,
`0.4576 v^2 <= r <= 0.50391 v^2` and `0.6764 v <= t <= 0.77237 v`, whence
`kt >= 2 v^3 * 0.6764 v = 1.3528 v^4 > 0.25590 v^4 >= r(r + 1)`,
`1 <= t < r` (as `0.77237 v < 0.4576 v^2` for `v >= 16`), `r >= 2`, and
`k < N`. Therefore `g` is nonconstant on every affine subspace of
`F_2^m` of dimension at least

```
d = floor(kr/t) + kt + r(r + 1) + 1,
```

i.e. `g` is an affine disperser for dimension `d`.

## 4. Elimination and the explicit remainder

Theorem E3 of `elimination.md` applied to `g` gives
`C_B2(g) >= 4m - 9(d + 1)`, so

```
C_B2(f_N) >= 4(N - k) - 9(d + 1) = 4N - R(N),
R(N) := 4k + 9 * ( floor(kr/t) + kt + r(r + 1) + 2 ).
```

It remains to bound `R(N) <= 29 N^{4/5}` for all `N >= 2^20`. Write
`v = N^{1/5} >= 16`, so `N^{2/5} = v^2 >= 256`, `N^{3/5} = v^3 >= 4096`,
`N^{4/5} = v^4`.

* `B`: from `(B - 1)^5 <= 32 N^3 < B^5` and `(2 N^{3/5})^5 = 32 N^3`,
  `2 v^3 < B <= 2 v^3 + 1`.
* `k`: by Dusart's thesis result (for `x >= 3275` there is a prime in
  `(x, (1 + 1/(2 ln^2 x)) x]`; P. Dusart, PhD thesis, Université de
  Limoges, 1998), applicable since `B >= 8193 > 3275`, there is a prime
  in `(B, (1 + 1/(2 ln^2 B)) B]` (if `B` is itself prime it serves
  itself). Since `B > 2^13` and `ln` is increasing, `ln B > 13 ln 2 >
  9.0109`, so `2 ln^2 B > 162.39` and the interval factor satisfies
  `1 + 1/(2 ln^2 B) < 1.0062`. The least prime `>= B` therefore
  satisfies `2 v^3 < k <= 1.0062 * (2 v^3 + 1) = 2.0124 v^3 + 1.0062
  <= 2.0140 v^3`, using `1.0062 <= 1.0062 v^3 / 4096`.
* `r = ceil(m/k)` with `m = N - k`: `r <= m/k + 1 < v^5 / (2 v^3) + 1 =
  v^2/2 + 1 <= 0.50391 v^2`, and
  `r >= m/k >= (v^5 - 2.0140 v^3) / (2.0140 v^3) = (v^2 - 2.0140) /
  2.0140 >= 0.49261 v^2 >= 0.4576 v^2`, using `v^2 >= 256`.
* `t = ceil(sqrt(r))`: `0.6764 v <= sqrt(0.4576) v <= t <= sqrt(r) + 1
  <= 0.77237 v`.
* `kr/t + kt`: for `t >= sqrt(r)` the identity
  `kt + kr/t = 2 k sqrt(r) + k (t - sqrt(r))^2 / t` holds (expand the
  square). Since `t = ceil(sqrt(r))`, `0 <= t - sqrt(r) <= 1`, so the last
  term is at most `k / t <= k / sqrt(r)`. Now `r <= m/k + 1` gives
  `k^2 r <= k m + k^2 = k (m + k) = k N` (note `m + k = N` exactly), hence
  `k sqrt(r) <= sqrt(k N) <= sqrt(2.0140) v^4 <= 1.41916 v^4`, and
  `k / sqrt(r) <= 2.0140 v^3 / (0.6764 v) <= 2.9777 v^2 <= 0.01164 v^4`.
  Therefore `kr/t + kt <= 2 * 1.41916 v^4 + 0.01164 v^4 <= 2.84996 v^4`.
* `r(r + 1) <= 0.50391 v^2 * 0.50782 v^2 <= 0.25590 v^4`.
* Hence `d + 1 <= 2.84996 v^4 + 0.25590 v^4 + 2 <= 3.10590 v^4`
  (using `2 <= 2 v^4 / 65536`), and

```
R(N) <= 4 * 2.0140 v^3 + 9 * 3.10590 v^4
     <= 0.50350 v^4 + 27.95310 v^4 = 28.45660 v^4 < 29 N^{4/5},
```

using `4 * 2.0140 v^3 = 8.0560 v^3 <= 8.0560 v^4 / 16 = 0.50350 v^4`
(since `v >= 16`).

This proves the theorem. (The bound is not tight: the analytic loss is
`28.46 N^{4/5}`, and exact evaluation of `R(N)` gives a maximum ratio
`R(N) / N^{4/5}` of about `28.28` on `N in [2^20, 2^20 + 2000000]`. We
declare the safe integer constant `29`.)

## 5. Model and quantifiers

The bound is for the single family of Section 1, at actual input length
`N`, for *every* integer `N >= 2^20` (no sparse subsequence, no limsup).
`R(N) / N = 29 N^{-1/5} -> 0`, so the eventual bound implies
`liminf C_B2(f_N)/N >= 4`, the declared `(4, 1)` projection, with the
same family supplying both coordinates jointly.

The circuit model is Li–Yang's full B_2 model (ECCC TR21-023, Sections
1.1, 2.1, 2.4, 3.3): all 16 binary truth tables, two pins per gate,
repeated pins and shared subgraphs allowed, one output, inputs/wires/
output selection free, every retained gate counted once. The elimination
proof of `elimination.md` works directly in this model: it starts from
an arbitrary such circuit, counts degenerate, constant and unary gates
in `s` until they are deleted by (D1)/(D2), absorbs negations into gate
operations, and covers the terminal output forms (constant; literal of a
parameter, auxiliary or gate; literal of an AND gate with affine inputs)
in the preliminary check of Lemma E2 via Lemma E1. Constant nodes created
by restrictions cost zero and are handled by rerouting, exactly the
normalisation the model prescribes. No uniformity is assumed of the
circuits; no direct-sum additivity or potential-to-gate identification
is used — the potential `W` is only an accounting device bounded by
`|C_0| + 5m` at the start and by `0` throughout, and the conclusion is a
statement about `|C_0|` alone. Witness bits are never circuit inputs.

## Sources

* Li–Yang, ECCC TR21-023 (2021-02-20), Sections 1.1, 2.1, 2.4, 3.3:
  circuit model and normalisation conventions. SHA-256
  4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062.
  (Their Theorem 1.1, the `3.1n - o(n)` bound, is the seed this result
  improves; the LP-cap argument of their Section 4 applies to their
  potential only and not to the present one.)
* Ben-Sasson–Kopparty, *Affine dispersers from subspace polynomials*,
  SIAM J. Comput. 41(4):880–914, 2012; preprint 2010-03-12, SHA-256
  922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf:
  Theorem 2.2 (equation (4)), Lemma 3.10, Lemma 3.13, Theorem 6.1,
  Lemmas 7.1–7.2, Theorem 7.3, Claim 7.4 — re-derived in `disperser.md`.
* Find–Golovnev–Hirsch–Kulikov, FOCS 2016 / Comput. Complexity 32 (2023),
  and Mishra–Santhanam, ECCC TR26-260 Revision 1 (2026-11-03), SHA-256
  a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921,
  Sections 2–4: the semi-circuit / parameter-domain framework that
  `elimination.md` follows; the `(5, 3)` weighting and charge 9 there are
  derived in full in `elimination.md`.
* P. Dusart, *Autour de la fonction qui compte le nombre de nombres
  premiers*, PhD thesis, Université de Limoges, 1998: for `x >= 3275`
  there is a prime in `(x, (1 + 1/(2 ln^2 x)) x]`.
