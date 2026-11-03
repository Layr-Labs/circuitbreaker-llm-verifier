# Coefficient 7/2 for an explicit affine disperser: complete proof

**Theorem.** Let `F = (F_N)` be the total family of Section 1. Then `F` is
computable by one fixed deterministic polynomial-time algorithm, hence its
language is in P and in NP with empty witnesses, and in the full B2 model of
Section 0,

    C_B2(F_N) >= (7/2)*N - 120*N^(4/5)   for every integer N >= 2^22.

In particular `liminf_{N->oo} C_B2(F_N)/N >= 7/2` over all integer lengths
`N = |x|`: for every `eps > 0`, taking
`N_eps = max(2^22, ceil((120/eps)^5))` gives
`C_B2(F_N) >= (7/2 - eps)*N` for all `N >= N_eps`. This is the joint
declaration in `bound.json`: `polynomial_log` guarantee
`c = 7/2, p = 1, q = r = 0` with the matching finite projection `(7/2, 1)`.

The proof has three mathematical ingredients, each reproduced in full below
in the notation of this file; no ingredient is used as an unproved citation.

* (D) The Ben-Sasson–Kopparty construction (preprint of 2010-03-12, Sections
  3, 6, 7) is an affine disperser: Sections 5–7 prove that the payload
  function of Section 1 is nonconstant on every affine subspace of dimension
  `> kr/t + kt + r(r+1)` of its block space.
* (E) The Mishra–Santhanam gate-elimination theorem (ECCC TR26-260,
  2026-10-19, Sections 2–4): Sections 8–9 prove that every affine disperser
  for dimension `d` needs at least `(7/2)n - 7d` B2 gates.
* (A) Assembly and accounting: Sections 3–4 and 10 fix an irreducible prefix
  (a restriction of the real input, all of whose `k = O(N^(3/5))` bits are
  charged to `N`), apply (D) and (E), and check every sufficiently large
  integer length.

## 0. Model conventions

Circuits are as in Li–Yang, ECCC TR21-023 (2021-02-20), Sections 1.1, 2.1,
2.4 (footnotes 3–4) and 3.3, which is also Section 2.1 of Mishra–Santhanam:
deterministic acyclic circuits, one output, fan-in two, every gate labelled
by one of the 16 binary Boolean truth tables (the full basis B2), unbounded
fan-out and depth, exact computation on every input. Input vertices are the
variables `x_1..x_n` or constants `0,1`; inputs, wires and output selection
cost zero; the size `|C|` is the number of gates, each retained gate counted
once, including constant and unary (degenerate) gates. Both pins of a gate
may come from the same vertex; subcircuits may be shared. We identify
`{0,1}` with the field `F_2`; addition is XOR, negation is `+1`, and every
operation has a unique multilinear form `c + b_1 t_1 + b_2 t_2 + b_12 t_1 t_2`.
The two constant and four unary operations are *affine*; of the ten
operations depending on both inputs, two are *XOR-type*
(`t_1 + t_2 + c`) and eight are *AND-type*
(`(t_1 + a)(t_2 + b) + g`). An operation with `b_12 = 0` is called *affine*.
Negations on wires are absorbed into the reading gate's table at no cost, as
both sources allow; a retained negated-projection gate still costs one, and
the elimination argument of Sections 8–9 never deletes a gate without
charging it in the size `s`. Constant vertices appearing during
substitutions cost zero and are handled by deletion rule D1 and the
normalisation discussion in Section 9.1; no uncharged computational basis
is ever introduced. The output may be a gate, a constant or a literal;
Section 8 (Lemma L3) rules out the degenerate output forms whenever the
disperser hypothesis is active. Inputs to circuits are the bits `x` alone,
never NP witnesses.

## 1. The family

All parameters are functions of the actual input length `N = |x|` only.

* `F_0` is the constant `0`.
* For `N >= 1`, let `B` be the least positive integer with `B^5 > 32*N^3`,
  and let `k` be the least prime `>= B`. If `k >= N`, set `F_N = 0`.
* Otherwise read the first `k` input bits `c_0, ..., c_{k-1}` as the lower
  coefficients of the monic polynomial
  `Q(Z) = Z^k + c_{k-1} Z^{k-1} + ... + c_0` over `F_2`. If `Q` fails the
  irreducibility test T of Section 2, set `F_N(x) = 0`.
* If `Q` passes T (so `Q` is irreducible, by Section 2), let
  `K = F_2[Z]/(Q)`, a field with `2^k` elements, identified with `F_2^k` by
  the polynomial basis `1, Z, ..., Z^{k-1}`. Put `m = N - k`,
  `r = ceil(m/k)`, `t = ceil(sqrt(r))`. Append `rk - m` zero bits to the
  remaining `m` input bits and split them into `r` blocks
  `u_1, ..., u_r in K`, block `i` read in the polynomial basis with its
  first bit the constant coefficient. For `1 <= i <= r` put
  `e_i = r + 1 - i` and `y_i = u_i^{a_i}` with
  `a_i = 2^{e_i+1} - 1 = sum_{j=0}^{e_i} 2^j` (descending powers: block 1
  gets the largest exponent). Let `E_t` be the `t`-th elementary symmetric
  polynomial and `Tr` the trace `K -> F_2`. Define

      F_N(x) = Tr( E_t(y_1, ..., y_r) ).

This defines one total Boolean function at every length. The prefix bits
are part of the input: they are read and tested by the evaluator (Section
3) and are counted in `N` in the lower bound (Section 10). Nothing is
chosen by advice: the polynomial `Q` is supplied by the input itself, and
for the lower bound we merely fix the prefix to one irreducible `Q`, whose
existence is proved in Section 2; fixing input bits to constants never
increases circuit size.

## 2. Prime-degree irreducibility: test, correctness, existence

**Test T.** For monic `Q` of prime degree `k >= 2`, T accepts iff
(i) `Z^{2^k} = Z (mod Q)`, i.e. `Q | Z^{2^k} - Z`, and
(ii) `gcd(Q, Z^2 + Z) = 1` (over `F_2`, `Z^2 + Z = Z(Z+1)`).

**Lemma F1.** T accepts iff `Q` is irreducible.

*Proof.* The polynomial `Z^{2^k} - Z` has derivative `1` in characteristic 2
(the derivative of `Z^{2^k}` is `2^k Z^{2^k - 1} = 0`), hence is separable;
its roots in an algebraic closure form the field `F_{2^k}`, and it equals
the product of all monic irreducible polynomials over `F_2` whose degree
divides `k`. Since `k` is prime, those degrees are `1` and `k` only.

(<=) If `Q` is irreducible of degree `k`, it divides `Z^{2^k} - Z`, giving
(i). It is not `Z` or `Z + 1` (degree `k >= 2`), and `Z^2 + Z` has only the
monic irreducible factors `Z, Z + 1`, so (ii) holds.

(=>) Assume (i) and (ii). Suppose `Q` reducible, and let `U` be a monic
irreducible factor of `Q` of degree `d <= k/2 < k` (a reducible monic
polynomial has an irreducible factor of degree at most half its degree, and
`U | Q | Z^{2^k} - Z`). Every root `a` of `U` satisfies `a^{2^k} = a`, so
`a in F_{2^k}`; also `U | Z^{2^d} - Z`, so `a in F_{2^d}`. Hence
`a in F_{2^k} ∩ F_{2^d} = F_{2^{gcd(k,d)}} = F_2`, since `k` prime and
`d < k` force `gcd(k,d) = 1`. So `U` has a root in `F_2`, and being monic
irreducible, `U = Z` or `U = Z + 1`; thus `U | Z^2 + Z` and `U | Q`,
contradicting (ii). Therefore `Q` is irreducible. ∎

**Lemma F2.** For every prime `k >= 2` there is a monic irreducible
polynomial of degree `k` over `F_2`; hence some input prefix passes T.

*Proof.* `Z^{2^k} - Z` has `2^k` distinct roots forming `F_{2^k}`. Any root
`a notin F_2` (there are `2^k - 2 > 0`) has minimal polynomial of degree
dividing `[F_{2^k}:F_2] = k` and exceeding `1`, hence equal to `k` (`k`
prime). That minimal polynomial is monic irreducible of degree `k`, and
every monic degree-`k` polynomial occurs as a prefix. ∎

**Prime location.** We use Bertrand's postulate (for every integer
`B >= 1` there is a prime in `[B, 2B]`; classical, Erdős 1932) only to
locate `k`: there is a prime in `[B, 2B]`, so the least prime `k >= B`
satisfies `k <= 2B`. All other prime facts used are proved above; the only
property of `k` used in Sections 5–7 is primality (so that `F_{2^k}` has
no subfields other than `F_2` and `F_{2^k}`).


## 3. Eligibility: one polynomial-time evaluator, NP membership

Fix one deterministic algorithm `A` computing `F_N` at all lengths `N`.

1. Compute `B` by binary search on `O(log N)`-bit integers (`B^5` compared
   against `32 N^3`), then find `k` by testing `B, B+1, ...` for primality
   by trial division up to the square root; Bertrand guarantees a prime
   `<= 2B`. Cost polynomial in `N` (in fact in `N^{9/10}` elementary
   operations on `O(log N)`-bit numbers, using the bounds of Section 4).
2. If `k >= N`, output `0`. Else form `Q` from the first `k` bits and run
   T: `Z^{2^k} mod Q` by `k` repeated squarings of polynomials of degree
   `< k` reduced mod `Q`, and `gcd(Q, Z^2 + Z)` by the Euclidean algorithm;
   each step costs polynomially many bit operations in `k <= N`.
3. If T fails, output `0`. Else compute in `K = F_2[Z]/(Q)`: form the `r`
   blocks, compute each `y_i = u_i^{a_i}` by repeated squaring (`a_i` has
   bit length `e_i + 1 <= r + 1`, so `O(r)` field multiplications per
   block), then `E_t(y_1, ..., y_r)` by the division-free recurrence
   `E_0 = 1`, `E_j <- E_j + y_i * E_{j-1}` for `j = min(i,t)` down to `1`;
   splitting `t`-subsets by membership of `i` proves the recurrence, and
   it uses `O(r t)` field operations. Finally output `Tr(E_t)`, the sum of
   the `k` Frobenius conjugates. Field operations in the polynomial basis
   cost polynomially in `k`. Total: polynomial in `N`.

So `A` is a single fixed polynomial-time evaluator for the total family.
Let `L = {x : F_{|x|}(x) = 1}` and define the verifier
`V(x, w) = 1` iff `w` is empty and `A(x) = 1`, with witness bound
`W(|x|) = 0`. If `x in L` then the empty `w` is accepted; if `x notin L`
no `w` is accepted (nonempty `w` is rejected outright, and the empty `w`
gives `A(x) = 0`). Both directions hold, so `L in P ⊆ NP`. The witness is
never a circuit input: lower bounds concern circuits on `x` alone.

## 4. Parameter lemma

Write `N^{a/b}` for the real power. Recall `B` is least with
`B^5 > 32 N^3`, `k` the least prime `>= B`, `m = N - k`,
`r = ceil(m/k)`, `t = ceil(sqrt(r))`, and set

    Lambda = k*r/t + k*t + r*(r+1),   D = floor(k*r/t) + k*t + r*(r+1) + 1,

so `D` is an integer with `D > Lambda`.

**Lemma P.** For every integer `N >= 2^22`:

(a) `2 N^{3/5} < B <= 2 N^{3/5} + 1` and `2 N^{3/5} < k <= 5 N^{3/5}`.
(b) `N/2 <= m < N`.
(c) `N^{2/5}/10 <= r <= (3/5) N^{2/5}`.
(d) `N^{1/5}/4 <= t <= 2 N^{1/5}` and `t < r`.
(e) `k >= r(r+1)/t` (the hypothesis of Theorem D).
(f) `D <= 17 N^{4/5}` and `(7/2) k + 7 D <= 120 N^{4/5}`.

*Proof.* `N >= 2^22` gives `N^{1/5} > 21`, `N^{2/5} > 445`,
`N^{3/5} > 9362`, `N^{4/5} = N / N^{1/5}`.

(a) `(B-1)^5 <= 32 N^3 < B^5`; taking fifth roots,
`B - 1 <= 2 N^{3/5} < B`, i.e. the claim for `B`. Bertrand (Section 2)
gives `B <= k <= 2B <= 4 N^{3/5} + 2 <= 5 N^{3/5}`, the last step since
`N^{3/5} >= 2`.
(b) `m = N - k >= N - 5 N^{3/5} >= N/2` since `N^{2/5} > 445 >= 10`.
(c) `r = ceil(m/k)` lies between `m/k` and `m/k + 1`.
`m/k < N / (2 N^{3/5}) = N^{2/5}/2`, so
`r <= N^{2/5}/2 + 1 <= (3/5) N^{2/5}` since `1 <= N^{2/5}/10`.
`m/k >= (N/2) / (5 N^{3/5}) = N^{2/5}/10`.
(d) `sqrt(r) <= t <= sqrt(r) + 1`. Upper: `sqrt(r) <= N^{1/5}` (as
`r <= (3/5) N^{2/5}`), so `t <= N^{1/5} + 1 <= 2 N^{1/5}`.
Lower: `t >= sqrt(r) >= sqrt(N^{2/5}/10) >= N^{1/5}/4` (as `sqrt(10) < 4`).
`t < r`: `t <= 2 N^{1/5} < N^{2/5}/10 <= r` since `N^{1/5} > 21 > 20`.
(e) Since `t >= sqrt(r)`, `r(r+1)/t <= sqrt(r) (r+1) <= 2 r^{3/2}`
(using `r >= 1`). And `r <= N^{2/5}/2 + 1 <= 0.503 N^{2/5}` (as
`1 <= N^{2/5}/445`), so `2 r^{3/2} <= 2 * 0.503^{3/2} N^{3/5}
< 0.72 N^{3/5} < 2 N^{3/5} < k`.
(f) `k r <= m + k <= N + 5 N^{3/5}` (from `r <= m/k + 1`), hence
`k r / t <= (N + 5 N^{3/5}) / (N^{1/5}/4) = 4 N^{4/5} + 20 N^{2/5}
<= 5 N^{4/5}` (as `20 N^{2/5} <= (20/445) N^{4/5} < N^{4/5}`). Also
`k t <= 5 N^{3/5} * 2 N^{1/5} = 10 N^{4/5}`, and
`r(r+1) <= (3/5) N^{2/5} * ((3/5) N^{2/5} + 1) <= N^{4/5}` for
`N^{2/5} > 445`. Summing, `D <= 5 N^{4/5} + 10 N^{4/5} + N^{4/5} + 1
<= 17 N^{4/5}`. Finally `(7/2) k <= (35/2) N^{3/5} <= N^{4/5}` (as
`N^{1/5} > 21 > 35/2`), so `(7/2) k + 7 D <= (1 + 119) N^{4/5}
= 120 N^{4/5}`. ∎

The guard `k >= N` and degenerate cases `r < 2` or `t >= r` can occur only
for small `N < 2^22` (Lemma P(b,d) excludes them above that threshold);
the evaluator defines `F_N = 0` there, which affects only finitely many
lengths and not the asymptotic statement.

## 5. Subspace polynomials (Ben-Sasson–Kopparty, Section 3)

From here `p = 2` throughout; we write the general prime `p` only where
the source argument is identical. Let `K = F_{2^k}` with `k` prime, and
identify `K` with its polynomial-basis copy as in Section 1. For a
nonnegative integer `s`, `wt(s)` is the sum of the binary digits of `s`.
For a monomial `prod X_i^{b_i}`, its *`F_2`-degree in `X_i`* is `wt(b_i)`;
the `F_2`-degree of a polynomial is the maximum over its monomials.

**Proposition A1 (reduction).** Reducing a polynomial
`h(X_1, ..., X_r)` modulo `X_i^{2^k} - X_i` for all `i` does not change
its values on `K^r`; the reduced form `h0` satisfies `|h(K^r)| > 1` iff
`h0` has a monomial of positive degree.
*Proof.* Values are unchanged since `a^{2^k} = a` on `K`. The reduced
form is the *unique* representative of the function `h : K^r -> K` with
every variable of degree `< 2^k`: the monomials
`prod X_i^{b_i}`, `0 <= b_i < 2^k`, interpolate arbitrary functions on
`K^r` (tensor product of the univariate Lagrange bases for `K`), so two
such representatives of one function coincide. Hence if `h` is constant
as a function, `h0` is the constant polynomial; and if `h0` has a
positive-degree monomial, `h0` is not the constant polynomial, so `h`
is not constant as a function. ∎

**Proposition A2 (trace form).** Every `F_2`-linear map
`pi : K -> F_2` is `pi(x) = Tr(mu x)` for a unique `mu in K`, and `pi` is
nontrivial iff `mu != 0`. (Standard: the trace form is nondegenerate.)
We use `pi = Tr`, i.e. `mu = 1`.

A polynomial `P(X) = sum_{j=0}^{k-1} a_j X^{2^j}` is *linearized*; it is
exactly an `F_2`-linear map `K -> K` in polynomial form, and affine
linearized polynomials (`+ const`) are exactly the polynomials of
`F_2`-degree `<= 1` in one variable. For an `F_2`-linear subspace
`L ⊆ K` of dimension `d`, its *kernel-subspace polynomial* is
`P_L(X) = prod_{a in L} (X - a)`.

**Lemma A3.** `P_L` is a monic linearized polynomial of degree `2^d`.
*Proof (Ore).* Induction on `d`. For `d = 1`, `L = {0, a}` and
`P_L = X(X + a) = X^2 + aX` (characteristic 2), which is linearized and
monic of degree `2`. For the step write `L = L0 ∪ (a + L0)` with
`dim L0 = d - 1` and `a notin L0`. Let `R(Y) = Y^2 + P_{L0}(a) Y`, the
kernel polynomial of `{0, P_{L0}(a)}`. Then `P_L := R ∘ P_{L0}` is a
composition of linearized polynomials, hence linearized; it is monic of
degree `2 * 2^{d-1} = 2^d`; and `P_L(x) = 0` iff
`P_{L0}(x) in {0, P_{L0}(a)}` iff `x in L0 ∪ (a + L0) = L`, so it equals
`prod_{a in L} (X - a)` (both monic of the same degree with the same
`2^d` roots). The affine case is identical with a constant added. ∎

**Lemma A4 (image-subspace polynomial).** For `L` as above there is a
monic linearized `Q_L` of degree `2^{k-d}` with `Q_L(K) = L`, and
`P_L ∘ Q_L = Q_L ∘ P_L = X^{2^k} - X`. In particular `Q_L` has exactly
`2^{k-d}` roots in `K`.
*Proof (Berlekamp).* Let `L0 = P_L(K)`, a subspace of dimension `k - d`,
and put `Q_L = P_{L0}` (Lemma A3 applied to `L0`). Then
`Q_L ∘ P_L` is monic of degree `2^k` and vanishes on `K`, so equals
`X^{2^k} - X`. Hence `P_L(Q_L(P_L(X))) = P_L(X^{2^k} - X)
= P_L(X)^{2^k} - P_L(X)`; writing `g(Y) = P_L(Q_L(Y)) - (Y^{2^k} - Y)` we
get `g(P_L(X)) = 0`, and `deg g(P_L(X)) = (deg g)(deg P_L)` forces
`g = 0`. So `P_L ∘ Q_L = X^{2^k} - X`, the image of `Q_L` lies in the
kernel `L` of `P_L`, and equality holds by dimension count
(`|im Q_L| = 2^k / |ker Q_L| = 2^d = |L|`). ∎

**Claim A5.** If `a^{2^s} - c a = 0` with `s >= 1` and `a != 0`, then
`c != 0` and `a in b * F_{2^s}` for any `(2^s - 1)`-th root `b` of `c`
in an algebraic closure of `K`.
*Proof.* `a^{2^s - 1} = c != 0` (as `a != 0`), so such a root `b`
exists in the closure and `(a/b)^{2^s - 1} = 1`, i.e.
`a/b in F_{2^s}`. ∎

**Claim A6.** For linearized `Q(X) = sum a_j X^{2^j} + â`,
`(Q(X))^{2^s} mod (X^{2^k} - X) = sum_j (a_{(j-s) mod k})^{2^s} X^{2^j}
+ â^{2^s}`. (Direct expansion; Frobenius is `F_2`-linear and
`X^{2^k} = X` after reduction.)

**Lemma A7 (main structural lemma).** Let `L ⊆ K` be a linear subspace of
dimension `d` and `Q_L(X) = sum_{j=0}^{k-1} a_j X^{2^j}` its
image-subspace polynomial.
(1) Any `d` consecutive indices in `Z_k` contain some `j` with `a_j != 0`.
In particular `a_0 != 0` and `a_{k-d} != 0` (the top `d - 1` coefficients
`a_{k-d+1}, ..., a_{k-1}` vanish since `deg Q_L = 2^{k-d}`).
(2) If moreover `d >= 2` and `L` is not contained in any constant multiple
of a proper subfield of `K`, then for any `r0 != k - d + 1`, any `d - 1`
consecutive indices starting at `r0` contain some `j` with `a_j != 0`.

*Proof.* Rotations: by Claim A6, for any `s` the reduction
`Q0 := Q_L^{2^s} mod (X^{2^k} - X)` has coefficient at position `j`
equal to `a_{j-s}^{2^s}` (indices mod `k`); in particular `Q0 != 0`
(some `a_j != 0`, e.g. the leading coefficient) and `Q0` vanishes
exactly where `Q_L` does (Proposition A1: vanishing is preserved under
Frobenius and under reduction mod `X^{2^k} - X`), i.e. on the `2^{k-d}`
roots of `Q_L` in `K` (Lemma A4).
(1) Suppose `a_j = 0` on a window `J` of `d` consecutive indices starting
at `r0`. Rotate with `s := (k - d) - r0 (mod k)`, which sends `J` to the
top `d` positions `{k-d, ..., k-1}`; then the top `d` coefficients of
`Q0` vanish, so `deg Q0 <= 2^{k-d-1}`, but `Q0 != 0` has `2^{k-d}`
roots — impossible. The "in particular": apply (1) to the window
starting at `k - d`, i.e. `{k-d, k-d+1, ..., k-1}`: the last `d - 1`
coefficients vanish by degree, so `a_{k-d} != 0`; and to the window
starting at `k - d + 1`, i.e. `{k-d+1, ..., k-1, 0}`: again the first
`d - 1` vanish by degree, so `a_0 != 0`.
(2) Suppose `a_j = 0` on a window of `d - 1` consecutive indices starting
at `r0 != k - d + 1`. Rotate with `s := (k - d + 1) - r0 (mod k)`, which
sends the window to the top `d - 1` positions `{k-d+1, ..., k-1}`. First,
`s != 0`: `s = 0` would mean `r0 = k - d + 1 (mod k)`, excluded; and
`s <= k - 1` by construction, so `0 < s < k`. (This is exactly where the
excluded start matters: for `r0 = k - d + 1` the window already sits at
the top and the rotation is trivial.) The top `d - 1` coefficients of
`Q0` are precisely the rotated window, hence vanish, so
`deg Q0 <= 2^{k-d}`. Choose `c in K`
killing the coefficient of `X^{2^{k-d}}` in `Q0` if `deg Q0 = 2^{k-d}`
(`c` = ratio of leading coefficients), and `c = 0` if
`deg Q0 < 2^{k-d}`. Then `Q0 - c Q_L` has degree `< 2^{k-d}` and
vanishes on all `2^{k-d}` roots of `Q_L` (both `Q0` and `Q_L` do), so
`Q0 = c Q_L`. And `c != 0`: else `Q0 = 0`, contradicting the rotation
paragraph. Thus `Q_L(a)^{2^s} = c Q_L(a)` for all `a in K`; since
`im Q_L = L`, Claim A5 gives `L ⊆ b * F_{2^s}` for some `b != 0`. Pick
`0 != beta in L`; then `beta^{-1} L ⊆ F_{2^s} ∩ K = F_{2^{gcd(s,k)}}`.
Since `k` is prime and `0 < s < k`, `gcd(s,k) = 1`, so
`L ⊆ beta * F_2`, a set of size 2 — contradicting `dim L = d >= 2`.
Hence some `a_j != 0` in the window. ∎

**Associated matrices.** For linearized `Q` as above, its associated
matrix `M_Q` has `(i, j)`-entry `(a_{j-i})^{2^i}` (indices mod `k`), and
`M_Q^0` is its zero-one indicator: `(M_Q^0)_{i,j} = 1` iff `a_{j-i} != 0`.

**Claim A8.** The `(i, j)`-entry of `M_Q` is the coefficient of
`X^{2^j}` in `(Q(X))^{2^i} mod (X^{2^k} - X)`. (From Claim A6.)

**Lemma A9 (coefficient formula).** Let `Q_1, ..., Q_t` be linearized,
`A_i, B_i ⊆ {0, ..., k-1}` nonempty with `|A_i| = |B_i|`,
`al_i = sum_{j in A_i} 2^j`, `be_i = sum_{j in B_i} 2^j`, `mu != 0`, and

    R(X_1, ..., X_t) = Tr( mu * prod_i Q_i(X_i)^{al_i} )
                       mod (X_i^{2^k} - X_i).

Then the coefficient `c_M` of `M = prod_i X_i^{be_i}` in `R` is

    c_M = sum_{r=0}^{k-1} mu^{2^r} * prod_i Perm( M_i[A_i + r, B_i] )
        = sum_{r=0}^{k-1} mu^{2^r} * prod_i Perm( M_i[A_i, B_i - r] )^{2^r}.

*Proof.* Expand the trace as `sum_r (mu * prod_i Q_i(X_i)^{al_i})^{2^r}`.
The `r`-th summand's contribution factors over `i`; for fixed `i`,
expanding `prod_{j in A_i} (Q_i(X_i))^{2^{j+r}}` and collecting
`X_i^{be_i} = prod_{s in B_i} X_i^{2^s}`, each bijection
`h : B_i -> A_i + r` contributes, by Claim A8,
`prod_{s in B_i} (M_i)_{h(s), s}`; summing over bijections is exactly the
permanent. ∎

**Remark A10 (how Lemma A9 is used).** If each indicator minor
`M_i^0[A_i, B_i]` is, up to row/column reordering, triangular with nonzero
diagonal, the `r = 0` summand is nonzero. If additionally for every
`r in {1, ..., k-1}` some `M_i^0[A_i, B_i - r]` has an all-zero column,
every other summand vanishes, so `c_M != 0`.


## 6. Independent sources (Ben-Sasson–Kopparty, Section 6)

**Theorem B1.** Let `Y_1, ..., Y_t ⊆ K` be *linear* subspaces of
dimensions `d_1, ..., d_t > 1`, none contained in a constant multiple of
a proper subfield of `K`. (For `k` prime the only proper subfield is
`F_2`, whose multiples are 1-dimensional, so the hypothesis is automatic
from `d_i > 1`.) Let `Q_i` be the image-subspace polynomial of `Y_i`
(Lemma A4). Let `mu != 0`, `1 <= e_i < d_i - 1`,
`al_i = sum_{j=0}^{e_i} 2^j`, and

    R = Tr( mu * prod_i Q_i(X_i)^{al_i} )  mod (X_i^{2^k} - X_i).

If `sum_i (d_i - (e_i + 1)) > k - max_i d_i + 1`, then `R` has a monomial
`prod_i X_i^{be_i}` with `wt(be_i) = e_i + 1` and nonzero coefficient; in
particular `|R(K^t)| > 1` (Proposition A1).

*Proof.* Write `Q_i(X) = sum_j a_{i,j} X^{2^j}` (linearized, Lemma A4)
and let `M_i^0` be
the indicator matrix of `Q_i`. We choose sets `B_i` with
`|B_i| = e_i + 1` and apply Lemma A9 with index sets
`{0, ..., e_i}` in the role of `A_i`. Order
so `d_1 = max d_i`. Put `l_1 = 0` and
`l_{i+1} = l_i + d_i - (e_i + 1) (mod k)`.
Case split on `a_{i, l_i}`:

* If `a_{i,l_i} != 0`: `B_i = {l_i} ∪ {k - d_i + 1, ..., k - d_i + e_i}`.
* If `a_{i,l_i} = 0`: let `j_i` be the least index `> l_i` with
  `a_{i,j_i} != 0` and `ĵ_i` the greatest index `< l_i` with
  `a_{i,ĵ_i} != 0` (indices mod `k`; both exist by Lemma A7(1) applied to
  the windows adjacent to the top vanishing block). Let
  `q_i = j_i - ĵ_i - 1` (the zero gap, cyclic) and `s_i = min(q_i, e_i)`.
  Set `B_i = {j_i} ∪ {ĵ_i + 1, ..., ĵ_i + s_i}
  ∪ {k - d_i + s_i + 1, ..., k - d_i + e_i}`.

**Claim B2.** Each `M_i^0[{0..e_i}, B_i]` is lower triangular with nonzero
diagonal, for a suitable ordering of the columns `B_i`.
*Proof.* Recall the `(u, κ)` entry of the minor (rows `u in {0..e_i}`,
columns `κ in B_i`) is the indicator of `a_{i, κ-u}` (indices mod `k`).
*Case 1* (`a_{i,l_i} != 0`): order the columns as
`(l_i, k - d_i + 1, ..., k - d_i + e_i)`. Diagonal entries: at `(0, l_i)`
the index is `l_i`, nonzero by assumption; at `(v, v)`, `v >= 1`, the
column is `k - d_i + v` and the index is `k - d_i`, nonzero by
Lemma A7(1). Above the diagonal (`u < v`, `v >= 1`) the index is
`k - d_i + (v - u)` with `1 <= v - u <= e_i <= d_i - 2`, hence lies in
`{k - d_i + 1, ..., k - 2}`, where all coefficients vanish by degree.
*Case 2* (`a_{i,l_i} = 0`): order the columns as
`(j_i, ĵ_i + 1, ..., ĵ_i + s_i, k - d_i + s_i + 1, ..., k - d_i + e_i)`.
For the left `(s_i + 1) x (s_i + 1)` block: the `j_i` column has index
`j_i - u`, which for `1 <= u <= s_i <= q_i` lies in the zero gap
`{ĵ_i + 1, ..., j_i - 1}`, so only its diagonal entry `(0,0)`, index
`j_i`, is nonzero. The column `κ = ĵ_i + v` (`1 <= v <= s_i`) has index
`κ - u = ĵ_i + (v - u)`: on the diagonal (`u = v`) this is `ĵ_i`,
nonzero by choice; above the diagonal (`u < v`) it lies in the zero gap
`{ĵ_i + 1, ..., ĵ_i + s_i} ⊆ {ĵ_i + 1, ..., j_i - 1}`; below the
diagonal it is unconstrained. The remaining columns
`κ = k - d_i + w` (`s_i + 1 <= w <= e_i`) sit at positions `v = w`, with
diagonal index `κ - w = k - d_i` (nonzero by Lemma A7(1)) and
above-diagonal indices `κ - u`, `u < w`, in
`{k - d_i + 1, ..., k - d_i + e_i} ⊆ {k - d_i + 1, ..., k - 1}`, which
vanish by degree — including the rows `u <= s_i` of the left block,
since there `κ - u >= k - d_i + 1` as `w >= s_i + 1`. ∎

**Claim B3.** For every `r in {1, ..., k-1}` some
`M_i^0[{0..e_i}, B_i - r]` has an all-zero column.
*Proof.* Column `s` of the shifted minor is zero if the cyclic interval
`[s - r - e_i, s - r]` lies inside the top zero block
`{k - d_i + 1, ..., k - 1}` of `Q_i` (Lemma A7: precisely the indices
`k - d_i + 1 .. k - 1` may vanish, and `a_{k-d_i} != 0` bounds the block).
So it suffices that `[k+1, k + d_i - (e_i + 1)]`-type intervals cover:
for each `i` and each `κ in B_i`, all shifts
`r in [κ + 1, κ + d_i - (e_i + 1)] (mod k)` give a zero column in factor
`i`. We show `[1, k-1] ⊆ ⋃_i ⋃_{κ in B_i} [κ+1, κ + d_i - (e_i+1)]`.
For `i = 1`: `l_1 = 0 in B_1`-position gives coverage of
`[1, d_1 - (e_1 + 1)] ∪ [k - d_1 + 2, k - 1] = [l_1 + 1, l_2] ∪
[k - d_1 + 2, k - 1]`. Let `t0` be minimal with
`sum_{i <= t0} (d_i - (e_i + 1)) >= k - d_1 + 1` (exists by hypothesis).
No wrap mod `k` occurs up to `t0`: by minimality
`sum_{i < t0} (d_i - (e_i + 1)) <= k - d_1` (the sums are integers and
the threshold is first reached at `t0`), and the last term satisfies
`d_{t0} - (e_{t0} + 1) <= d_1 - 2` (as `d_{t0} <= d_1` and
`e_{t0} >= 1`), so `l_{t0 + 1} <= (k - d_1) + (d_1 - 2) = k - 2 < k`;
hence every `l_i` for `i <= t0 + 1` is a genuine integer in `[0, k)` and
every interval below is a genuine integer interval.
For `1 < i <= t0` we claim the union over `κ in B_i` covers
`[l_i + 1, l_{i+1}]`: if `a_{i,l_i} != 0` then `l_i in B_i` and we are
done by the displayed interval with `κ = l_i`. If `a_{i,l_i} = 0` and
`q_i < e_i`, the set `B_i` rewrites as the interval
`[ĵ_i + 1, ĵ_i + q_i + 1] ∪ {k - d_i + q_i + 1, ..., k - d_i + e_i}` and
again `l_i in B_i` (as `ĵ_i < l_i <= j_i - 1 <= ĵ_i + q_i`). If
`q_i >= e_i`, the elements of `B_i` in `[ĵ_i + 1, ĵ_i + e_i]` contribute
intervals `[κ + 1, κ + d_i - (e_i + 1)]` which are consecutive (each has
length `d_i - (e_i + 1) >= 1` and starts one past the previous) and
therefore cover `[ĵ_i + 2, ĵ_i + d_i - 1]`; and `j_i in B_i` contributes
the further interval `[j_i + 1, j_i + d_i - (e_i + 1)]`. Lemma A7(2) —
applicable because `Y_i` is not contained in a constant multiple of a
proper subfield and `ĵ_i + 1 != k - d_i + 1` (from
`ĵ_i < l_i <= k - d_1 <= k - d_i` and maximality of `d_1`) — gives
`j_i - ĵ_i <= d_i - 1`, i.e. `j_i + 1 <= ĵ_i + d_i`, so the `j_i`
interval starts no later than one past the end of the covered range and
the union covers `[ĵ_i + 2, j_i + d_i - (e_i + 1)]`, which contains
`[l_i + 1, l_{i+1}]`: at the bottom because `ĵ_i < l_i`, and at the top
because `j_i >= l_i` gives
`j_i + d_i - (e_i + 1) >= l_i + d_i - (e_i + 1) = l_{i+1}`.
Chaining these coverings from `i = 1` to `t0` covers `[1, k-1]`: the
intervals `[l_i + 1, l_{i+1}]` tile `[1, l_{t0+1}]` with
`l_{t0+1} >= k - d_1 + 1`, and the `i = 1` top block covers
`[k - d_1 + 2, k - 1]`. ∎

Claims B2, B3 and Remark A10 now give a monomial with nonzero
coefficient. ∎ (Theorem B1)

## 7. The disperser theorem (Ben-Sasson–Kopparty, Section 7)

**Lemma C1 (Bourgain decomposition).** Let `A ⊆ K^r` be `F_2`-affine,
`γ in A`. There are linear subspaces `Y_1, ..., Y_r ⊆ K` and `F_2`-linear
maps `σ_{ij} : Y_j -> K` (`j < i`) with
`A = { (x_i) : x_i = γ_i + y_i + sum_{j<i} σ_{ij}(y_j), y_i in Y_i }`
and `dim A = sum_i dim Y_i`.
*Proof.* Induction on the block index. Project `A` to the first `i`
blocks. The fibre of this linear projection over any fixed prefix is a
coset of one fixed linear space; call its `i`-th block component space
`Y_i`. Hence, once the earlier free components `y_1, ..., y_{i-1}` are
fixed, the `i`-th component has the form
`γ_i + y_i + (an F_2-affine function of y_1, ..., y_{i-1})`, and the
linear part of that function in `y_j` is `σ_{ij}`. Each step adds
exactly `dim Y_i` new degrees of freedom, so dimensions add. ∎

**Lemma C2 (top-weight transfer).** With `P_i, P_{ij}` linearized,
`γ ∈ K^r`, `I_0 ⊆ [r]`, `al_i` as in Theorem B1, let
`g = Tr( mu * prod_{i in I_0} (P_i(X_i) + sum_{j<i} P_{ij}(X_j) + γ_i)^{al_i} )`
and `g0 = Tr( mu * prod_{i in I_0} P_i(X_i)^{al_i} )`, both reduced mod
`(X_i^{2^k} - X_i)`. For every tuple `(be_i)_{i in I_0}` with
`wt(be_i) = e_i + 1`, the coefficients of `prod X_i^{be_i}` in `g` and
`g0` are equal.
*Proof.* Expand each factor as in the source: writing
`P̂_i = sum_{j<i} P_{ij}(X_j) + γ_i`, the `F_2`-degree of
`(P_i(X_i) + P̂_i)^{al_i}` in `X_i` splits as the pure term
`prod_{l=0}^{e_i} (P_i(X_i))^{2^l}` of `X_i`-degree exactly `e_i + 1`
plus a remainder `G_i` of `X_i`-degree `<= e_i` (each mixed term replaces
at least one linearized factor in `X_i` by one not mentioning `X_i`, and
linearized polynomials have `F_2`-degree 1). Order `I_0` decreasingly;
inductively, to reach `X_{i}^{be_i}`-degree `e_i + 1` in the largest
remaining index one must take the pure term there, since the cofactor
does not mention that variable. Hence only pure terms contribute, which
is exactly `g0`. ∎

**Theorem D (disperser, unparameterised).** Let `t < r`, let `k` be prime
with `k >= r(r+1)/t`. For `1 <= i <= r` let `e_i = r + 1 - i`,
`al_i = sum_{j=0}^{e_i} 2^j`, `mu != 0`, and

    f(x_1, ..., x_r) = Tr( mu * sum_{I ⊆ [r], |I| = t} prod_{i in I} x_i^{al_i} ).

Then `f` is nonconstant on every `F_2`-affine `A ⊆ K^r` with
`dim A > kr/t + kt + r(r+1)`.

*Proof.* Fix such `A`, `dim A =: dA`. Decompose `A` by Lemma C1 with
spaces `Y_i`, `d_i = dim Y_i`, `sum d_i = dA`. Let `Q_i` be the
image-subspace polynomial of `Y_i`, and choose linearized `Q_{ij}` with `Q_{ij}(x) = σ_{ij}(Q_j(x))`
(exists: extend `σ_{ij} : Y_j -> K` arbitrarily to an `F_2`-linear map
`K -> K`; the composite `σ_{ij} ∘ Q_j : K -> K` is then `F_2`-linear,
hence a linearized polynomial by the correspondence of Section 5, and
its values do not depend on the extension because `im Q_j = Y_j`).
Then `H(X) = (R_1, ..., R_r)`,
`R_i = Q_i(X_i) + sum_{j<i} Q_{ij}(X_j) + γ_i`, has image exactly `A`:
as `X_j` ranges over `K`, `Q_j(X_j)` ranges over `Y_j`, so `R_i` ranges
over `γ_i + Y_i + sum_{j<i} σ_{ij}(Y_j)`, which is Lemma C1's form.
By Proposition A1 it suffices to show `h := f ∘ H` (reduced) has a
positive-degree monomial with nonzero coefficient.

Let `S = {i : d_i > r + 1}`. Each `d_i <= k`, and since the terms
`d_i - (r+1)` for `i notin S` are `<= 0`,
`sum_{i in S} (d_i - (r+1)) >= dA - r(r+1) > kr/t + kt`. Also
`|S| >= t`: if `|S| <= t - 1`, then
`sum_{i in S} (d_i - (r+1)) <= (t-1)(k - r - 1) <= (t-1) k < kt
< kr/t + kt`, contradicting the displayed lower bound on the sum. Choose `I_0 ⊆ S` with
`|I_0| = t` maximising `sum_{i in I_0} (d_i - (r+1))`. Averaging over all
`t`-subsets of `S` (the average of this sum is
`(t/|S|) * sum_{i in S} (d_i - (r+1))`, and `|S| <= r`, so the average is
at least `(t/r)` times the total) gives

    sum_{i in I_0} (d_i - (r+1))
      > (t/r) * (kr/t + kt) = k + k t^2 / r >= k,

the inequality being strict because
`sum_{i in S} (d_i - (r+1)) > kr/t + kt` is strict. (No use of
`k >= r(r+1)/t` is needed at this step; that hypothesis is in any case
available from Lemma P(e).) Since `e_i + 1 <= r + 1` for all `i`,

    sum_{i in I_0} (d_i - (e_i + 1)) >= sum_{i in I_0} (d_i - (r+1)) > k.

Theorem B1's hypothesis is the *strict* inequality
`sum (d_i - (e_i + 1)) > k - max_i d_i + 1` with the maximum taken over
`I_0`; it follows because every `i in I_0 ⊆ S` has
`d_i > r + 1 >= 2`, so `k - max_i d_i + 1 <= k - 1 < k < sum`.
Moreover every `i in I_0` has `d_i > r + 1 >= e_i + 1`, so
`1 <= e_i < d_i - 1` holds, and `Y_i` (a
linear space of dimension `> 1` in the prime-degree field `K`, whose
only proper subfield is `F_2`) is not contained in a constant multiple
of a proper subfield. Theorem B1 applied to `(Y_i)_{i in I_0}` yields a
monomial `M = prod_{i in I_0} X_i^{be_i}`, `wt(be_i) = e_i + 1`, with
nonzero coefficient in
`g0 = Tr(mu prod_{i in I_0} Q_i(X_i)^{al_i})`. By Lemma C2 the same holds
in `g = Tr(mu prod_{i in I_0} R_i^{al_i})`, the `I_0`-summand of `h`.

No other summand contributes `M`: for `J != I_0`, `|J| = t`, expanding the
trace and the powers writes the `J`-summand as a product over
`j in J`, `l in [0, r - j + 1]` of affine linearized polynomials
`R_j^{2^{l+w}}` (each a sum of monomials of `F_2`-degree `<= 1`). If `M`
occurred, choosing for each `(j, l)` the variable index `λ_j(l) ∈ [1, j]`
of the monomial used, the `F_2`-degree of `M` in `X_i` gives
`sum_{j in J} |λ_j^{-1}(i)| >= r + 2 - i` for all `i in I_0` (since
`wt(be_i) = e_i + 1 = r + 2 - i`).

**Claim C3.** Under these inequalities, `J = I_0`.
*Proof.* For each `j`, `sum_{i in I_0} |λ_j^{-1}(i)| <= r + 2 - j`
(only `r + 2 - j` values of `l`). Summing the hypothesis over `i in I_0`
gives `sum_{I_0} (r + 2 - i) <= sum_J (r + 2 - j)`, so with
`|I_0| = |J| = t`, `sum_{I_0} i >= sum_J j`. More generally, for every
`s >= 0`, `sum_j sum_l (j^s - λ_j(l)^s) >= 0` (each `λ_j(l) <= j`),
combined with the hypothesis weighted by the multiplicities
`|λ_j^{-1}(i)| >= r + 2 - i`, gives
`sum_{I_0} i^{s+1} - sum_J j^{s+1} >= (r + 2)(sum_{I_0} i^s - sum_J j^s)`,
hence by induction on `s` the power sums are ordered:
`p_s(I_0) := sum_{I_0} i^s >= sum_J j^s =: p_s(J)` for all `s >= 0`.

We now compare the elements from the top down (this is the source's
argument; power-sum dominance alone does not give it). Let
`i_1 = max I_0`, `j_1 = max J`. If `i_1 < j_1`, then for large `s`,
`p_s(J) >= j_1^s` while `p_s(I_0) <= t * i_1^s`, contradicting
`p_s(I_0) >= p_s(J)`; so `i_1 >= j_1`. On the other hand the hypothesis
requires the value `i_1` to occur among the `λ_j(l)` at least
`r + 2 - i_1 >= 2` times (as `i_1 <= r`), but every `λ_j(l) <= j <= j_1`;
if `i_1 > j_1` it could not occur at all. Hence `i_1 = j_1`. Moreover
the occurrences of the value `i_1` can only come from `j = j_1` (for
`j < j_1`, `λ_j(l) <= j < i_1`), and `j_1` supplies exactly
`r + 2 - j_1 = r + 2 - i_1` values of `l`, all of which must therefore
equal `i_1`: `λ_{j_1}(l) = i_1` for all `l`.

Delete `i_1` from `I_0`, `j_1` from `J`, and all pairs `(j_1, l)` from
the `λ`-family. The remaining data satisfy the same hypotheses:
cardinalities drop to `t - 1` on both sides; the occurrence counts for
`i != i_1` are unchanged (the deleted pairs all had value `i_1`); the
bound `λ_j(l) <= j` persists; and the power sums stay ordered, since
`p_s(I_0 \ {i_1}) = p_s(I_0) - i_1^s >= p_s(J) - j_1^s
= p_s(J \ {j_1})` using `i_1 = j_1`. Repeating the argument `t` times
shows the sorted lists of `I_0` and `J` coincide element by element;
both are sets, so `J = I_0`. ∎

So `M` occurs in `h` with the nonzero coefficient from the `I_0`
summand alone; `h` is nonconstant, and `f` is nonconstant on `A`. ∎

**Corollary D1 (our parameters).** Fix an irreducible prefix `Q` in
`F_N` (`N >= 2^22`) and let `g` be `F_N` with the prefix bits fixed to
that prefix, a function of the `m` payload bits. Then `g` is nonconstant
on every affine subspace of `F_2^m` of dimension `>= D` (Section 4).
*Proof.* The zero-padding embedding `F_2^m -> K^r ≅ F_2^{rk}` is
`F_2`-linear and injective, so an affine `A ⊆ F_2^m` maps to an affine
subspace of the same dimension in `K^r`, and `g = f ∘ pad` there with
`f` as in Theorem D (`mu = 1`), since `Tr` and the exponents `al_i`
match Section 1 exactly (`al_i = 2^{r+2-i} - 1`). If `dim A >= D >
Lambda = kr/t + kt + r(r+1)`, Theorem D — whose hypotheses `t < r` and
`k >= r(r+1)/t` are Lemma P(d,e) — makes `f`, hence `g`, nonconstant on
`A`. ∎


## 8. Gate elimination setup (Mishra–Santhanam, Section 2–3)

We prove Theorem E (Section 9): *every affine disperser `f` on `n` bits
for dimension `d` (nonconstant on every affine subspace of dimension `d`,
hence on every larger one, which contains a `d`-dimensional affine
subspace) has `C_B2(f) >= (7/2) n - 7 d`.* Sections 8–9 reproduce the
source proof in full.

**Parameters and maps.** The algorithm maintains an injective map
`Phi : F_2^m -> F_2^n` of the form

    Phi(u) = b + L u + sum_{{p,r} in P} c_{p,r} * p * r,        (1)

where `u = (u_1, ..., u_m)` are *parameters*, `b in F_2^n`,
`L in F_2^{n x m}`, `P` is a set of disjoint parameter pairs, and
`c_{p,r} in F_2^n`. A parameter in a pair is *paired* (its mate is its
*partner*); others are *unpaired*. Write `q = |P|` and `a = m - q`.
Together with a circuit `C` on the parameters satisfying the invariant
`C(u) = f(Phi(u))` for all `u`. Initially `Phi = id`, `m = n`, `P = ∅`,
`C` the given circuit.

**Proposition E1.** `im(Phi)` contains an affine subspace of dimension
`a`; hence `f ∘ Phi` is nonconstant whenever `a >= d`.
*Proof.* Fix one member of each pair to a constant: every quadratic term
becomes affine, the restricted map stays injective (as `Phi` is), and its
image is affine of dimension `m - q = a`. If `a >= d` it contains (or is)
a `d`-dimensional affine subspace on which `f` is nonconstant. ∎

**Fair semi-circuits.** The circuit is kept as a *fair semi-circuit*: an
affine system `A v = N u + c` (`A` invertible `r x r`; every row of
`(A | N)` has at most three nonzero entries; `v` are *auxiliary
variables*) followed by an acyclic part (an ordinary B2 circuit whose
gates read constants, parameters, auxiliary variables or earlier gates),
with output a constant or a literal of a parameter, auxiliary variable or
gate. Fairness (`A` invertible) gives a unique extension of every `u` to
all *internal values* (parameters, auxiliary variables, gates). Size
`s` = (#equations) + (#gates); initially the system is empty and
`s = |C|`. An affine gate whose two inputs are parameters/auxiliary
variables/constants can be *moved into the system*: replace gate `J`
with the equation `J + b_1 t_1 + b_2 t_2 = c`; the new auxiliary matrix is
block-triangular with a `1` on the diagonal, hence invertible; `s` is
unchanged and the function is unchanged. We move every affine gate with
no AND-type ancestor into the system.

**Degrees, weights, potential.** The *children* of an internal value are
the gates/equations reading it (a double pin counts once); `deg` is the
number of children. A parameter is *unused* if `deg = 0`, *isolated* if
`deg = 1`. When a value becomes constant or a literal of another value we
*reroute* all its uses accordingly. Weights of unpaired parameters:

    wt(x) = 0 if deg(x) = 0;  6 if deg(x) = 1;  7 if deg(x) >= 2,

and the potential is `W = 2 s + 3 q + sum_{unpaired x} wt(x) >= 0`.
Paired parameters contribute only through `3q`. For an operation,
`DX = X' - X`.

**Parameter changes.** Four operations, each replacing `Phi` by
`Phi ∘ tau` (injectivity and form (1) preserved; the invariant is
maintained by applying `tau` to `C` as specified):
(AR) *Affine reparameterisation:* for unpaired `x` and affine `phi`,
`x <- w + phi(u_{-x})`, `w` fresh. `Dm = Dq = Da = 0`, `im Phi` unchanged.
(QR) *Quadratic reparameterisation:* for distinct unpaired `x, y, z`,
`z <- w + (x + a)(y + b) + g`, adding the pair `{x, y}` to `P`.
`Dm = 0`, `Dq = +1`, `Da = -1`, `im Phi` unchanged.
(UF) *Fixing an unpaired parameter:* `x = psi(u_{-x})` (`psi` affine;
we use only constants and `x = y + c`). `Dm = -1`, `Dq = 0`, `Da = -1`.
(PF) *Fixing a paired parameter:* for `{p, r} in P`, fix `p` to a constant
or impose `p = r + c`, removing `p` and the pair. `Dm = -1`, `Dq = -1`,
`Da = 0`.

**Implementation without new gates.** (UF)/(PF) substitute a constant or
a literal in `C` directly — no new gates. (AR) is implemented by *role
exchange*: if an auxiliary variable `Q` has solved form
`Q = x + phi(u_{-x})` with `x` unpaired, substitute
`x <- w + phi(u_{-x})` in `Phi` only, rename `Q` as the new parameter
`w`, and regard `x` as auxiliary. No equation or gate changes; only the
`x`-column of `N` and the `Q`-column of `A` are swapped. The new `A'` is
`A B` where `B` is the identity with column `j` (index of `Q`) replaced
by `t = A^{-1} N_x`; `det B = t_j = 1` because `Q` depends on `x`, so
`A'` is invertible, and every equation keeps at most three nonconstant
terms. Thus `Ds = Dm = Dq = Da = 0` and
`DW = wt'(w) - wt(x)` with `deg'(w) = deg(Q)`; `w` occurs in at least one
equation (its `A`-column is nonzero). We apply (AR) only when
`deg(x) >= 2` and `Q` has another child besides its equation occurrence
as needed, giving `wt(x) = wt'(w) = 7`, `DW = 0`; and inside D4 below,
where the possible `+1` is offset by a deletion. (QR) is applied only in
the configuration of Case 10 (Section 9), where it deletes two gates;
the affine system is untouched since `z` occurs in no equation there.
Rerouting preserves both invariants: substituting `v = U + c` (from an
equation with `<= 2` nonconstant terms) into other equations is a sequence
of elementary row operations on `(A | N | c)` followed by expansion along
the now-singleton column of `v`, preserving `det A != 0` and the
three-term bound; rerouting inside the acyclic part does not touch the
system (no wires go back).

**Lemma L3 (output forms).** If `a > d`, the output is neither an affine
function of the parameters nor a literal of an AND-type gate on two
unpaired parameters.
*Proof.* A constant output contradicts Proposition E1. If the output is
a nonconstant affine function, pick a parameter `x` with nonzero
coefficient; if `x` is paired, fix its partner to a constant (PF, `Da =
0`) making `x` unpaired; then a (UF) substitution for `x` makes the
output constant while `Da = -1`, so the restricted image still contains
an affine subspace of dimension `>= a - 1 >= d`, contradicting E1. If the
output is a literal of an AND gate on unpaired `x, y`, fix `y` to a
constant making the gate constant: again output constant with `Da = -1`,
contradiction. ∎

**Deletions** (no change to `Phi`; `Dm = Dq = Da = 0`, `Ds = -1` each).
(D1) A gate computing a constant or a literal of one input — every
constant/unary gate, every gate with a constant input or two equal
inputs, and every AND-type gate whose two inputs are
parameters/auxiliaries with dependent linear parts — is rerouted to that
constant/literal and deleted. An unpaired isolated parameter whose only
child was this gate can gain weight `6 -> 7` only if it acquires
`>= 2` children by the rerouting; hence `DW <= -2 - 6 + 7 = -1`, and
`DW <= -2` if rerouted to a constant.
(D2) A gate with no children that is not the output is deleted:
`DW <= -2`.
(D3) An equation with `<= 2` nonconstant terms has the form
`v = U + c` (`v` auxiliary — otherwise the `A`-row would be zero);
substitute everywhere, delete `v` and the equation (Section 8,
rerouting). As in D1, `DW <= -1`, and `<= -2` if no weight increases.
(D4) *Unit parameter deletion:* unpaired isolated `x` whose unique child
is an equation `E`. The `x`-column `N_x` of `N` is supported on `E`'s
row, so some auxiliary `Q` depends on `x` (`A^{-1} N_x != 0`); write
`Q = x + phi(u_{-x})`. Apply (AR) with this `Q`: now `x` is auxiliary,
still occurring only in `E`. Delete `x` and `E` (the column of `x` in
the new `A` has a single nonzero entry; determinant expansion as in
rerouting). The old weight `6` of `x` is removed, the new parameter `Q`
has weight `<= 7`, no other weight increases: `DW <= -2 - 6 + 7 = -1`,
with `Ds = -1`, `Dm = Dq = Da = 0`.

**Quadratic reparameterisation accounting (used in Case 10).** Suppose
the circuit contains an AND-type gate `G = (x + a)(y + b) + g` whose only
child is an XOR-type gate `D = G + z + c`, where `z` is an isolated
unpaired parameter whose only child is `D`, and `x, y` are unpaired with
`deg >= 2`. Apply (QR): `z <- w + (x + a)(y + b) + g + c`, adding pair
`{x, y}`. Then `D = w`; reroute `D`'s uses to `w`, delete `D`, then `G`
is unused and deleted by D2. `z` occurred nowhere else, so nothing else
changes. `Dm = 0`, `Ds = -2`, `Dq = +1`, `Da = -1`, and
`DW <= -4 (gates) + 3 (pair) + 7 (w) - 6 (z) - wt(x) - wt(y)
= -wt(x) - wt(y) <= -14`. ∎


## 9. The elimination algorithm and Theorem E

**Lemma E2.** If `a > d`, one iteration transforms the state so that
`m + s` strictly decreases, `Da in {0, -1, -2}`, and `DW <= 14 Da`,
in polynomial time.

Each iteration first moves all affine gates with no AND ancestor into
the system (`DW <= 0`, no degree increases), then performs the first
applicable of Cases 1–6 below; if none applies, Sections 9.2–9.3
restructure and eliminate a minimal AND (Cases 7–10). Preparatory moves,
(AR) exchanges and the reroutings of Section 9.2 preserve `m, s, a` and
do not increase `W`. In each case, earlier cases are assumed
inapplicable.

### 9.1 Preliminary reductions

**Case 1: a deletion is available.** Apply D1–D4 repeatedly until none
applies. Each single deletion has `Ds = -1`, `DW <= -1`, `Da = 0`, hence
already satisfies the iteration requirement `DW <= 14 Da` and strictly
decreases `m + s`; a maximal block of deletions is therefore a valid
iteration on its own, and after it we re-enter the case analysis.
When no deletion applies: every gate is XOR- or AND-type with two
distinct nonconstant inputs; an AND gate on two parameters/auxiliaries
has independent linear parts there; every non-output gate has a child;
every equation has exactly three nonconstant terms; no unpaired isolated
parameter occurs in the system (else D4).

**Case 2: an AND gate on parameters has an unpaired isolated input.**
Let `G` be AND-type on parameters `z, w` with `z` unpaired isolated. Fix
`w` to a constant making `G` constant; delete `G` by D1. If `w` is
unpaired, this is (UF): `Dm = -1, Ds = -1, Da = -1`, and `z, w` lose
weights `>= 6` each: `DW <= -6 - 6 - 2 = -14 = 14 Da`. If `w` is paired
with partner `r`, this is (PF): `Dm = -1, Ds = -1, Dq = -1, Da = 0`;
`r` becomes unpaired with weight `<= 7`: `DW <= -3 - 6 - 2 + 7 = -4
<= 0`.

*Consequence of Cases 1–2 failing.* Every isolated unpaired parameter
`z` is not the output (Lemma L3), does not occur in the system, so its
unique child is an acyclic gate `G`. The other input of `G` is not a
parameter: if `G` were XOR-type on two parameters it would have been
moved into the system; if AND-type, Case 2 would apply. Hence no gate or
equation has both an isolated unpaired parameter and another parameter
as parents.

**Case 3: a paired parameter has degree `0` or `>= 2`.** Let `p` be
paired, partner `r`. If `deg(p) = 0`: fix `r` to a constant by (PF);
`p` becomes unpaired of weight `0`: `Dm = -1, Dq = -1, Da = 0`,
`DW <= -3`. If `deg(p) >= 2`: fix `p` to a constant by (PF) and delete
two of its children by D1–D4; only the newly unpaired `r` can gain
weight (`<= 7`): `Dm = -1, Ds = -2, Dq = -1, Da = 0`,
`DW <= -3 - 4 + 7 = 0`.
*If Case 3 fails, every paired parameter is isolated.*

**Case 4: two partners share a child.** Partners `p, r` with common
child `H`. If `H` is an equation it reads `v + p + r = c`; if a gate, its
inputs are `p, r`. Impose `p = r + c'` by (PF), choosing `c'` so `H`
becomes constant (any `c'` for an equation or XOR gate; for AND-type
`H = (p + a)(r + b) + g` take `c' = a + b + 1`). Delete `H` (D1/D3);
`r` becomes unused: `Dm = -1, Ds = -1, Dq = -1, Da = 0`,
`DW <= -3 - 2 = -5`.

**Case 5: a paired parameter occurs in the system.** Paired `p` in the
system, partner `r`; by Cases 3–4, `p`'s unique child is an equation
`H_p`, `r`'s unique child `H_r != H_p`. Fix `r` to a constant by (PF)
and delete `H_r` (D1/D3); `H_r` has no isolated-unpaired parent (Cases
1–2 failing), so no initially-unpaired weight increases. Now `p` is
unpaired isolated in equation `H_p` with weight `6`: so far
`DW <= -3 - 2 + 6 = +1`. Apply D4 to `p, H_p`: `DW <= -1` more.
Totals: `Dm = -1, Ds = -2, Dq = -1, Da = 0`, `DW <= 0`.

**Case 6: an AND gate on parameters has a paired input.** `G` AND-type
on parameters `p, w`, `p` paired with partner `r`. By Cases 3–5: `G` is
`p`'s unique child; `r`'s unique child `H_r` is a gate, `H_r != G`;
`w != p, r`; `w` is not unpaired isolated (Case 2), so `w` is unpaired
of degree `>= 2` or paired. Fix `r` to a constant by (PF); fix `w` to a
constant making `G` constant ((UF) if unpaired, (PF) if paired); delete
`G` and `H_r` by D1. No initially-unpaired weight increases (`G` rerouted
to a constant; `H_r` has no isolated-unpaired parent).
If `w` was unpaired: `Dm = -2, Ds = -2, Dq = -1, Da = -1`, `w` loses
weight `7`: `DW <= -7 - 3 - 4 = -14 = 14 Da`.
If `w` was paired (partner `r'`): `Dm = -2, Ds = -2, Dq = -2, Da = 0`;
only `r'` newly unpaired, weight `<= 7`: `DW <= -6 - 4 + 7 = -3 <= 0`.

*If Cases 1–6 all fail:* every parameter in the system is unpaired of
degree `>= 2`, and every AND gate on two parameters has both inputs
unpaired of degree `>= 2`.

### 9.2 Restructuring the minimal AND

Assume `a > d` and no case of Section 9.1 applies. The acyclic part has a
gate (else the output would be affine in the parameters, against Lemma
L3). Let `G` be its topologically first gate; `G` is AND-type (affine
gates on system values were moved into the system), and its two inputs
are parameters or auxiliaries with independent linear parts (Case 1
fails). *Making the inputs parameters:* if an input is an auxiliary
`Q`, its solved affine form depends on some parameter `x`; `x` occurs in
the system (else its `N`-column, hence its `A^{-1} N`-column, would be
zero), so `x` is unpaired of degree `>= 2` (Section 9.1 summary). Apply
(AR) exchanging `Q` and `x`: the new parameter `w` has `deg'(w) =
deg(Q) >= 2` (`Q` occurs in an equation and is a parent of `G`), and
`wt(x) = wt'(w) = 7`, so `DW = 0`. Re-run the moves and Cases 1–6; if a
case applies, perform it and end the iteration. Repeat for the second
input if needed (linear independence is preserved by the invertible
exchange, so a parameter other than `w` exists). Now both inputs `x, y`
of `G` are unpaired parameters of degree `>= 2`. *Clearing children:*
if a child `J` of `G` also reads `x`, relabel `J` by the two-variable
operation `(x, y) |-> h(x, G(x, y))` and reroute its inputs to `x, y`:
the function is unchanged, `J` stops being a child of `G`, no weight
increases (`deg(y) >= 2` already, symmetrically for `y`). Re-run moves
and cases after each rerouting; each rerouting removes a child of `G`
and none is created, so this terminates. If it ends with no case
applicable, no child of `G` reads `x` or `y`, and `G` is not the output
and has `>= 1` child (Lemma L3, Case 1).

### 9.3 Eliminating the minimal AND

Order `deg(x) >= deg(y) >= 2`; let `k_c >= 1` be the number of children
of `G` (renamed to avoid clash with the field degree `k`; only here).

**Case 7: `deg(x) = deg(y) = 2` with a common second child `B != G`.**
`B` is an equation `v + x + y = c` (exactly three terms; `v` auxiliary)
or a gate on `x, y` with operation `h`. As a function on `F_2^2`, the
AND-type `G` is constant on three points; two of them share an `h`-value,
and two points of `F_2^2` lie on a line `a x + b y = g` on which both `G`
and `h` are constant. Restrict to that line by one (UF) (fix `x`, fix
`y`, or impose `x = y + g`). Delete `G` (D1) and `B` (D1, or D3 if an
equation, rerouting the auxiliary `v` to its constant). Both `x, y`
become unused and lose weight `7`: `Dm = -1, Ds = -2, Da = -1`,
`DW <= -7 - 7 - 4 = -18 <= -14`.
*Henceforth, if `deg(x) = deg(y) = 2`, the second children
`B_x != B_y` are distinct.*

**Case 8: the generic restriction.** Suppose `deg(x) >= 3`, or
`k_c >= 2`, or the unique child of `G` has no unpaired isolated other
parent. Let `k_1 <= k_c` be the number of children of `G` whose other
parent is an unpaired isolated parameter. Fix `x` to a constant making
`G` constant ((UF)); delete all `deg(x)` children of `x` (each by D1 if
a gate, D3 if an equation — all become deletable once `x` is fixed),
then all `k_c` children of `G` (gates of the acyclic part, deletable by
D1 once `G` is rerouted to its constant). These `deg(x) + k_c` values
are distinct (no child of `G` reads `x`). Parameter `x` loses weight
`7`; each deletion contributes `-2`. Weight increases: only an initially
isolated unpaired `z` can gain (`6 -> 7`), and only if its unique child
is deleted and rerouted to a literal of `z`; that child is not a child
of `x` (no gate/equation has both an isolated unpaired parameter and
another parameter as parents), so it must be one of the `k_1` children
of `G` counted above. Hence
`DW <= -7 - 2(deg(x) + k_c) + k_1`, with `Dm = -1`, `Da = -1` and
`Ds = -(deg(x) + k_c)` (the deleted gates/equations), so `m + s`
strictly decreases.
If `deg(x) = 2` (so `deg(y) = 2` too), there is one more unit: `y`'s
children are `G` and `B_y`; `y` is not a parent of `B_x` (distinctness)
and no child of `G` reads `y`, so after `G` is rerouted `y` becomes
isolated and its weight drops `7 -> <= 6`. Therefore

    DW <= -7 - 2 deg(x) - k_c <= -14        if deg(x) >= 3
                                             (since -2 k_c + k_1 <= -k_c <= -1),
    DW <= -12 - 2 k_c + k_1 <= -14         if deg(x) = 2, k_c >= 2,
    DW <= -14 + k_1                        if deg(x) = 2, k_c = 1.

So Case 8 satisfies `DW <= 14 Da` except in the *bottleneck
configuration* `deg(x) = deg(y) = 2`, `k_c = k_1 = 1`: `G` has a unique
child `D` whose other parent is an unpaired isolated parameter `z`.
This configuration is tested for before Case 8 is applied and is handled
by Cases 9–10.

**Case 9: bottleneck, `D` AND-type.** Write
`D = (G + d0)(z + h0) + t0`. Choose `x_0, y_0` with
`G(x_0, y_0) = d0` (possible: `G` is nonconstant). Fix `x = x_0`,
`y = y_0` (two (UF)s). Then `G` computes `d0` and `D` computes `t0`;
delete `G, D` by D1, and also `B_x, B_y` by D1/D3 (deletable once `x, y`
are fixed). The four are distinct (`B_x != B_y`; neither equals `D`,
which is a child of `G` and reads neither `x` nor `y`). `z` becomes
unused. No weight increases: `G, D` are rerouted to constants, and
`B_x, B_y` have no isolated-unpaired parent (Cases 1–2 fail). So
`x, y, z` lose `7, 7, 6` and four deletions give `-8`:
`Dm = -2, Ds = -4, Da = -2`, `DW <= -28 = 14 Da`.

**Case 10: bottleneck, `D` XOR-type.** `D = G + z + c`: exactly the
(QR) configuration of Section 8 (AND gate `G` on unpaired `x, y` of
degree `2`, unique child `D`, isolated `z` with unique child `D`).
Apply (QR) as analysed there: `Dm = 0, Ds = -2, Dq = +1, Da = -1`,
`DW <= -14 = 14 Da`.

This proves Lemma E2: every iteration does preparatory steps
(`DW <= 0`, `m, s, a` preserved) plus one case with
`Da in {0, -1, -2}` and `DW <= 14 Da`; deletions and (QR) decrease `s`
with `m` fixed, every other case decreases `m`, and no case increases
`s`, so `m + s` strictly decreases. Each iteration is polynomial-time:
at most `s` affine moves, at most two (AR) exchanges, at most one
rerouting per child of `G`, Gaussian elimination over `F_2` for solved
forms and dependence tests, and constant-size truth-table checks. ∎

**Theorem E (Mishra–Santhanam Theorem 4.1).** Every affine disperser
`f : F_2^n -> F_2` for dimension `d` satisfies `C_B2(f) >= (7/2) n - 7 d`.
*Proof.* Initialise as in Section 8: `m = a = n`, `q = 0`, `s = |C|`,
`W_0 <= 2|C| + 7n`. Apply Lemma E2 while `a > d`; `m + s` strictly
decreases from `n + |C|` and stays `>= 0`, so the process stops with
`a_T <= d` after finitely (polynomially) many iterations. Summing
`DW <= 14 Da` over all iterations: `W_T - W_0 <= 14 (a_T - n)`. With
`W_T >= 0` and `a_T <= d`:
`2|C| + 7n >= W_0 >= 14 (n - a_T) >= 14 (n - d)`, so
`|C| >= (7/2) n - 7 d`. ∎

*Remark (optimality of the weights, source Remark 4.3).* For potentials
`V = s + (3/2) q + (7/2) I - (1/2) l` (`I` = used unpaired parameters,
`l` = isolated unpaired parameters) — which is `W / 2` in our notation —
the two bottleneck-adjacent cases are simultaneously tight, and no
choice of weights of this form improves the coefficient `7/2` in this
analysis. We make no claim beyond `7/2`.

## 10. Assembly

Let `N >= 2^22` and fix an irreducible prefix `Q` (Lemma F2). Let `g` be
`F_N` with the `k` prefix bits fixed to that prefix. Any B2 circuit for
`F_N` becomes a B2 circuit for `g` of no larger size by substituting the
prefix constants (substitution cannot create gates; degenerate gates are
kept and still counted), so `C_B2(F_N) >= C_B2(g)`.

By Corollary D1, `g` is nonconstant on every affine subspace of
`F_2^m` of dimension `>= D`; in particular it is an affine disperser for
dimension `D` in the sense of Theorem E (nonconstant on every affine
subspace of dimension exactly `D`; any larger subspace contains a
`D`-dimensional affine subspace, so the Proposition E1 consequence used
throughout Sections 8–9 holds at every `a >= D`). If `D > m`, the bound
`(7/2) m - 7 D` about to be derived is negative, hence trivially at most
`C_B2(g)`; if `D <= m`, Theorem E applies with `n := m`, `d := D`.
Either way

    C_B2(g) >= (7/2) m - 7 D = (7/2)(N - k) - 7 D
             = (7/2) N - [ (7/2) k + 7 D ]
             >= (7/2) N - 120 N^{4/5}        (Lemma P(f)).

This holds for every integer `N >= 2^22`, with the one fixed family,
verifier and witness bound of Section 3. Since `120 N^{4/5} = o(N)`,
`liminf C_B2(F_N) / N >= 7/2`, and the epsilon form in the Theorem
statement follows with `N_eps = max(2^22, ceil((120/eps)^5))`: for
`N >= N_eps`, `120 N^{4/5} <= eps * N`. All lengths are covered — the
argument is not restricted to a sparse subsequence — and the coefficient
and exponent come from this single theorem for this single family. ∎

## Sources

* Mishra–Santhanam: ECCC TR26-260, 2026-10-19, Theorem 4.1 and Sections
  2–4 (Theorem 1.1's `O(log n)`-disperser corollary is *not* used; the
  class theorem is applied to the explicit disperser below).
  https://eccc.weizmann.ac.il/report/2026/260/download/
  SHA-256 `38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc`.
  Reproduced as Theorem E in Sections 8–9.
* Ben-Sasson–Kopparty: author preprint 2010-03-12, Theorem 2.2/7.3,
  Lemma 7.1, Lemma 7.2, Theorem 6.1, Claim 7.4, Section 3
  (Lemmas 3.7, 3.10, 3.13 and Claims 3.8, 3.9, 3.12).
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA-256 `922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf`.
  Reproduced as Theorem D in Sections 5–7.
* Li–Yang: ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4, 3.3, for
  the pinned B2 conventions only (their 31/10 LP is not used).
  SHA-256 `4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062`.
* Bertrand's postulate (prime in `[B, 2B]`), used only in Section 2/4.
