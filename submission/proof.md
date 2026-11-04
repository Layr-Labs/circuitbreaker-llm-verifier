# CircuitBreaker: C_B2(f_N) >= 3N - O(N^(4/5)) for one uniform family in P

**Declared result.** For the single family `f` of Section 1 (a total
language in P, hence in NP), at actual input length `N`,

    C_B2(f_N) >= 3N - 252 N^(4/5)   for every integer N >= 2^29,

in the full Li-Yang B2 model. Hence `liminf_N C_B2(f_N)/N >= 3`: the
declared pair is `(c, p) = (3, 1)` with no logarithmic factors.

The proof is self-contained. It has three parts:

- `proof.md` (this file): the family and its eligibility, the parameter bounds,
  the **elimination theorem** (Theorem E, an affine-restriction gate-elimination
  argument with labelled inputs in the style of Demenkov and Kulikov, MFCS 2011),
  and the final assembly with all quantifiers;
- `disperser.md`: **Theorem D**, a complete proof that the family is an affine
  disperser for dimension `O(N^(4/5))` (a full reproduction of the
  Ben-Sasson-Kopparty subspace-polynomial argument);
- `fields.md`: elementary facts (irreducibility test for prime degree,
  existence of irreducible polynomials, the field `K`, primes in `(X, 8X]`).

No other premise is used. Literature names are given for provenance only;
every step needed is proved in these files.

## 1. The family

Fix `N >= 1` (`f_0 = 0` on the empty string). Let `k = k_N` be the least
prime `k` with `k^5 >= 32 N^3` (i.e. `k >= 2 N^(3/5)`). Put `m = N - k`,
`r = ceil(m/k)` and `t = ceil(sqrt(r))` (least integer with `t^2 >= r`).
Call `N` **admissible** if `k < N`, `r >= 2` and `r + 2 <= k`. If `N` is not
admissible, `f_N` is the constant 0.

If `N` is admissible, write the input as `w = (q_0,...,q_{k-1}, x_1,...,x_m)`
and let `q(z) = z^k + sum_{j<k} q_j z^j in F_2[z]`.

- If the test `T(q)` of `fields.md` F2 fails (i.e. `z^(2^k) != z mod q`, or
  `q(0) = 0`, or `q(1) = 0`), then `f_N(w) = 0`.
- Otherwise `q` is irreducible (F2.3) and `K_q = F_2[z]/(q)` is a field with
  `2^k` elements (F3). Append `rk - m` zeros to `x` and, for `1 <= i <= r`,
  let `u_i = sum_{j=0}^{k-1} x_{(i-1)k+j+1} z^j in K_q` (polynomial basis).
  Let `y_i = u_i^(alpha_i)` with `alpha_i = 2^(r+2-i) - 1`, and

      f_N(w) = Tr( E_t(y_1,...,y_r) ),   E_t(y) = sum_{I subset [r], |I|=t} prod_{i in I} y_i,

  where `Tr(a) = sum_{j<k} a^(2^j)` lies in `{0,1}` (F3).

This defines a single language `L = {w : f_{|w|}(w) = 1}`. Every length
`N` is covered, including small and non-admissible ones. There is no advice:
the "field" is part of the input, the function is 0 on inputs whose first `k`
bits fail `T`, and the lower bound below fixes those bits to one irreducible
polynomial whose existence is proved in F2.4.

## 2. Eligibility: `L` is in P, hence in NP

The following deterministic algorithm `Alg` computes `f_{|w|}(w)`:

1. Scan `k = 2, 3, ..., N-1`, testing primality by trial division and the
   integer inequality `k^5 >= 32 N^3`. If no such prime is found below `N`,
   output 0. This takes at most `N` candidates of `O(sqrt N)` divisions each.
2. Compute `m, r, t` by integer arithmetic; if `r < 2` or `r + 2 > k`, output 0.
3. Evaluate `T(q)` by `k` squarings modulo `q`; if it fails output 0.
4. Compute `u_i`; `y_i = u_i^(alpha_i)` by repeated squaring (`alpha_i < 2^(r+1)`,
   so `O(r)` field multiplications each).
5. Compute `E_t` without division: start `e_0 = 1`, `e_j = 0` (`1 <= j <= t`);
   for `i = 1..r` and `j = min(i,t)` down to 1 set `e_j <- e_j + y_i e_{j-1}`.
   After step `i`, `e_j` is the sum over `j`-subsets of `{1..i}` (split by
   whether `i` belongs to the subset). This is `O(rt)` field operations.
6. Output `Tr(e_t)` (`k` squarings).

Each field operation costs `poly(k)` bit operations (F3), and `k, r, t < N`,
so `Alg` runs in time polynomial in `N`. Verifier: `V(w, c) = 1` iff `c` is
empty and `Alg(w) = 1`; witness bound `W(N) = 0`. If `w in L` the empty
witness is accepted; if `V(w, c) = 1` then `c` is empty and `Alg(w) = 1`, so
`w in L`. Witness bits are not circuit inputs.

## 3. Parameter bounds

**Lemma P.** For every integer `N >= 2^29`: (a) `2 N^(3/5) <= k <= 16 N^(3/5)`;
(b) `N` is admissible, with `m >= N/2` and `2 <= r <= N^(2/5)`;
(c) `sqrt(r) <= t < sqrt(r) + 1`;
(d) `Lambda := kr/t + kt + r(r+1) <= 50 N^(4/5)`.

*Proof.* (a) The lower bound is the definition. For `N >= 2^29`,
`X = 2N^(3/5) >= 2^18`, so by F4 there is a prime in `(X, 8X]`, and `k`, the
least prime `>= X`, is at most `8X = 16 N^(3/5)`.
(b) `k <= 16 N^(3/5) <= N/2` (as `N^(2/5) >= 2^11 >= 32`), so `k < N` and
`m >= N/2`. Then `r >= m/k >= N^(2/5)/32 >= 2`, and
`r <= m/k + 1 <= N/(2N^(3/5)) + 1 <= N^(2/5)`. Finally
`r + 2 <= N^(2/5) + 2 <= 2N^(3/5) <= k`.
(c) is the definition of `t`.
(d) `kr/t <= k sqrt(r) <= 16 N^(3/5) N^(1/5) = 16 N^(4/5)`;
`kt <= k(sqrt(r) + 1) <= 16 N^(4/5) + 16 N^(3/5) <= 32 N^(4/5)`;
`r(r+1) <= 2r^2 <= 2 N^(4/5)`. QED.

## 4. The disperser property

**Lemma 4.1.** Let `N` be admissible and let `q` satisfy `T(q)`. Define
`g_q : F_2^m -> {0,1}` by `g_q(x) = f_N(q, x)`. Then `g_q` is not constant on
any affine subspace of `F_2^m` of dimension `> Lambda`.

*Proof.* `k` is prime with `k >= r + 2 >= 4`, so `k >= 5`; `r >= 2`;
`1 <= t <= r`. The map `phi(x) = (u_1,...,u_r)` from `F_2^m` to `K_q^r` is
`F_2`-linear and injective (zero padding followed by coordinates), so it
maps an affine subspace of dimension `D` onto an affine subspace of dimension
`D`. With `e_i = r + 1 - i` we have `alpha_i = 2^(e_i+1) - 1`, so
`g_q = F o phi` for the function `F` of Theorem D (`disperser.md`), whose
hypotheses (`k` prime, `r >= 2`, `1 <= t <= r`, `r + 2 <= k`) hold. Theorem D
gives the claim. QED.

## 5. The elimination theorem

**Theorem E.** Let `m >= 1`, `G : F_2^m -> F_2`, and let `d*` be an integer
such that `G` is not constant on any affine subspace of `F_2^m` of dimension
`>= d*`. Then every B2 circuit computing `G`, possibly containing constant
nodes, has at least `3m - 4d*` gates.

The circuit model is that of `rules.md`: acyclic, one output, two (possibly
equal) pins per gate, any of the 16 binary tables, unrestricted fan-out and
depth, constant nodes and inputs free, every gate counted (including
degenerate ones). In this section all subspaces are affine subspaces of
`F_2^m`, and affine functions `lambda(x) = c_0 + sum c_i x_i` are over `F_2`.

### 5.1 Labelled circuits

Fix an affine subspace `S`. Say `lambda` is *constant on `S`* if its
restriction to `S` is constant, and `lambda, lambda'` are *`S`-equivalent*
if `lambda + lambda'` is constant on `S`. A *hyperplane* of `S` is an affine
subspace `S' subset S` with `dim S' = dim S - 1`; if `lambda` is not constant
on `S` then `S cap {lambda = c}` is a hyperplane of `S` for each `c`.

A *labelled circuit over `S`* is a finite acyclic circuit whose nodes are
constant nodes (value 0 or 1), *label nodes* `v` (each carrying an affine
function `lambda_v`) and gates (an ordered pair of pins from earlier nodes,
possibly equal, and any table `g : F_2^2 -> F_2`), with one output node. On
`x in S` every node gets a value (`c`, `lambda_v(x)`, or `g` of the pin
values). It *computes `G` on `S`* if the output value equals `G(x)` for all
`x in S`. The *successors* of a node are the gates having it as a pin. Let
`s` be the number of gates, `iota` the number of label nodes, and the
**measure** `mu = s + iota`. Let `mu*(S)` be the least measure of a labelled
circuit over `S` computing `G` on `S`. Label nodes are a bookkeeping device:
each costs exactly 1 in `mu`, like a gate; only the inequality
`s >= mu - m` for ordinary circuits (Section 5.5) is used at the end.
Lemmas E1 and E2 are this proof's complete replacement for Li-Yang
normalization: they treat constant nodes, unary and degenerate gates,
repeated pins and negated or constant outputs explicitly, and they never
add an uncharged computational element.

*Table types.* A table depending on both arguments is either **XOR-type**,
`g(a,b) = a + b + c` (the 2 tables with exactly two ones that depend on both
arguments), or **AND-type**, `g(a,b) = (a + al)(b + be) + ga` (the 8 tables
with exactly one or three ones); tables with zero or four ones are constant,
and the remaining 4 tables with two ones depend on one argument only. An AND-type gate is constant (`= ga`) whenever its
first argument equals `al`, and whenever its second equals `be`.

*Hypothesis H(S).* `G` is not constant on `S` nor on any hyperplane of `S`.
Under H(S), `G` is not of the form `lambda + c` on `S`: if `lambda` were
constant on `S`, `G` would be constant on `S`; otherwise `G` would be
constant on the hyperplane `S cap {lambda = 0}`.

### 5.2 Restriction lemma

**Lemma E1.** Let `Gamma` compute `G` on `S`, and let `S' subset S` be an
affine subspace on which `G` is not constant. Let `Z` be the set of gates
having at least one pin whose value is constant on `S'`, and `Lambda_0` the
set of label nodes whose function is constant on `S'`. Then there is a
labelled circuit `Gamma'` over `S'` computing `G` on `S'` with
`s(Gamma') = s(Gamma) - |Z|` and `iota(Gamma') = iota(Gamma) - |Lambda_0|`.
Moreover, if a label node `v` not in `Lambda_0` is not the output of `Gamma`,
and its only successor in `Gamma` is a gate `A in Z` whose value is constant
on `S'`, then `v` has no successor in `Gamma'` and is not its output.

*Proof.* In a topological order assign to each node `w` a *representative*:
either a constant `c` with `val_w = c` on `S'`, or a pair `(w', b)` with `w'`
a label node not in `Lambda_0` or a gate not in `Z`, `w'` equal to `w` or an
ancestor of `w`, and `val_w = val_{w'} + b` on `S'`:

- a constant node `c`: representative `c`; a label in `Lambda_0` with
  constant value `c` on `S'`: `c`; any other label `v`: `(v, 0)`; a gate not
  in `Z`: itself, `(Y, 0)`;
- a gate `Y in Z`: some pin has a value constant (`= c`) on `S'`, so on `S'`
  `val_Y = h(val_p)` where `p` is the other pin (or the same pin, if both
  pins coincide) and `h` is unary. If `val_Y` is constant on `S'` its
  representative is that constant. Otherwise `h` is the identity or the
  negation and `val_p` is not constant, so `p` has a representative
  `(p', b)`; give `Y` the representative `(p', b + [h = neg])`.

(If a node's representative is a constant, its value on `S'` is constant.)
Let `Gamma'` keep exactly the label nodes not in `Lambda_0` and the gates not
in `Z`, with the same label functions. A kept gate `Y` with pins `p_1, p_2`
and table `g` has no pin with constant value on `S'`, so each `p_i` has a
representative `(p'_i, b_i)`; in `Gamma'` give `Y` the pins `p'_1, p'_2` and
the table `g(a_1 + b_1, a_2 + b_2)`. The pins precede `Y`, so `Gamma'` is
acyclic, and by induction each kept node has on `S'` the same value as in
`Gamma`. The output `o` of `Gamma` has non-constant value on `S'`, so it has
a representative `(o', b)`. Make `o'` the output; if `b = 1`, complement
`o'` (complement its table if it is a gate, replace `lambda_{o'}` by
`lambda_{o'} + 1` if it is a label) and complement the corresponding argument
of every gate of `Gamma'` having `o'` as a pin, which keeps the values of
all other nodes. Then `Gamma'` computes `G` on `S'` with the stated counts.

For the last statement, suppose some node `w != v` had a representative
`(v, b)`, and take the first such `w`. It is not a constant node, a label,
or a gate outside `Z` (these are represented by constants or by
themselves), so `w in Z` and its representative comes from a pin `p` with
representative `(v, b')`; by minimality `p = v`, so `w` is a successor of
`v`, i.e. `w = A`. But `A` has constant value on `S'` and hence a constant
representative, a contradiction. Thus a gate of `Gamma'` has pin `v` only if
it had pin `v` in `Gamma`, i.e. only `A`, which is not kept; and the output
`o'` of `Gamma'` equals `v` only if `o = v`, which is excluded. QED.

### 5.3 Normalization lemma

**Lemma E2.** Assume H(S) and let `Gamma` compute `G` on `S`. Then some
`Gamma*` computes `G` on `S`, has `mu(Gamma*) <= mu(Gamma)`, and satisfies:
(N1) no constant nodes; (N2) every label function is non-constant on `S`
and no two distinct label nodes are `S`-equivalent; (N3) every gate has two
distinct pins and a table depending on both arguments; (N4) every node other
than the output has a successor; (N5) no XOR-type gate has two label nodes
as pins; (N6) the output is a gate.

*Proof.* Apply the following rules while one applies. A gate is
*degenerate* if a pin is a constant node, or its two pins coincide, or its
table does not depend on one argument; a degenerate gate has value
`h(val_w)` on `S` for one of its pins `w` and `h in {0, 1, id, neg}`.

- R1. A label node with function constant (`= c`) on `S` becomes a constant
  node `c`.
- R2. If distinct label nodes `u, v` have `lambda_u + lambda_v = c` on `S`,
  redirect every pin at `v` to `u`, complementing that argument of the table
  when `c = 1`, and delete `v`. (`v` is not the output: `G` is not affine on
  `S`.)
- R3. Remove a degenerate gate `Y`. If `h` is constant, `Y` becomes a
  constant node. If `h = id`, redirect all pins at `Y` (and the output, if
  `Y` is the output) to `w`. If `h = neg`, redirect all pins at `Y` to `w`
  with that argument complemented; if moreover `Y` is the output, then `w` is
  not a constant node (`G` would be constant) nor a label (`G` would be
  affine), so `w` is a gate: complement the table of `w`, complement the
  `w`-argument of every successor of `w` (including the pins just
  redirected), and make `w` the output.
- R4. Delete a gate, label node or constant node that has no successor and
  is not the output.
- R5. Replace an XOR-type gate `Y`, `g(a,b) = a + b + c`, whose pins are two
  distinct label nodes `u, v`, by a new label node with function
  `lambda_u + lambda_v + c` (pins at `Y` now point to it). `Y` is not the
  output, since `G` is not affine on `S`.

Each rule preserves the value of every remaining node on `S` (up to the
compensated complementations) and hence the computed function. R1, R2 and R4
do not increase `s` and decrease `iota` or the number of nodes; R3 decreases
`s` by 1 without changing `iota`; R5 decreases `s` by 1 and increases `iota`
by 1. So `mu` never increases, and the triple (`s`, `iota`, number of
constant nodes) decreases lexicographically at each step; the process stops.
At the end: a constant node would either have a successor (R3 applies) or
none (R4 applies; it is not the output as `G` is not constant on `S`), so
(N1) holds; (N2) by R1, R2; (N3) by R3; (N4) by R4; (N5) by R5; and the
output is neither a constant node nor a label (`G` is neither constant nor
affine on `S`), giving (N6). QED.

### 5.4 One elimination step

**Lemma E3.** If H(S) holds, then `S` has a hyperplane `S'` with
`mu*(S') <= mu*(S) - 4`.

*Proof.* Let `Gamma` have measure `mu*(S)` and let `Gamma*` be given by
Lemma E2; then `mu(Gamma*) = mu*(S)`. By (N6) there is a gate; let `A` be a
first gate in a topological order. Its pins are not gates, hence (N1) label
nodes `u != v` (N3). By (N3) and (N5) `A` is AND-type:
`A = (val_u + al)(val_v + be) + ga`. By (N2) `lambda_u` is not constant on
`S`, so `S' = S cap {lambda_u = al}` is a hyperplane of `S`, and by H(S) `G`
is not constant on `S'`. On `S'`, `u in Lambda_0` and `A in Z` has constant
value `ga`. `A` is not the output (else `G` would be constant on `S'`), so
by (N4) `A` has a successor.

*Case 1: `u` has a successor `B != A`.* Then `B in Z`.
(1a) If `A` has a successor `C != B`, then `C in Z` and `A, B, C` are
distinct. (1b) Otherwise `B` is the only successor of `A`; by (N3) the pins
of `B` are exactly `u` and `A`, both constant on `S'`, so `B` has constant
value on `S'`; `B` is not the output (else `G` would be constant on `S'`), so
by (N4) it has a successor `E`, which lies in `Z`, and `E != A` (`A` is a
pin of `B`, so it precedes `B`), `E != B`. In both subcases `|Z| >= 3` and
`|Lambda_0| >= 1`, so Lemma E1 gives `mu*(S') <= mu*(S) - 4`.

*Case 1': `v` has a successor other than `A`.* Symmetric, with
`S' = S cap {lambda_v = be}`.

*Case 2: `A` is the only successor of `u` and of `v`.* Use
`S' = S cap {lambda_u = al}`. A successor `C` of `A` lies in `Z`, so
`|Z| >= 2`, and `u in Lambda_0`. Also `v` is not in `Lambda_0`: if
`lambda_v` were constant on `S'`, the linear part of `lambda_v` on the
direction space `W` of `S` would vanish on the kernel of the (nonzero) linear
part of `lambda_u` on `W`, hence equal 0 or that linear part, so `lambda_v`
would be constant on `S` or `S`-equivalent to `lambda_u`, contradicting
(N2). `v` is not the output (N6), its only successor is `A in Z`, and `A`
has constant value on `S'`. By Lemma E1, `Gamma'` has
`s(Gamma') <= s - 2`, `iota(Gamma') <= iota - 1`, and `v` has no successor in
`Gamma'` and is not its output. Deleting `v` gives a labelled circuit over
`S'` computing `G` with measure at most `mu*(S) - 4`. QED.

### 5.5 Proof of Theorem E

*Claim.* For every affine subspace `S` with `dim S >= d*`,
`mu*(S) >= 4(dim S - d*)`. Induction on `D = dim S`. For `D = d*` this is
`mu* >= 0`. For `D > d*`, `S` and its hyperplanes have dimension
`>= D - 1 >= d*`, so H(S) holds; Lemma E3 gives a hyperplane `S'` with
`mu*(S) >= mu*(S') + 4 >= 4(D - 1 - d*) + 4`.

Now let `C` compute `G` with `s` gates. If `m < d*` then `3m - 4d* < 0` and
there is nothing to prove. Otherwise turn each input node `x_i` of `C` into a
label node with function `x_i` (constant nodes stay). This is a labelled
circuit over `S = F_2^m` computing `G`, with measure `s + m`. By the claim,
`s + m >= 4(m - d*)`, i.e. `s >= 3m - 4d*`. QED.

## 6. The lower bound and the quantifiers

**Theorem.** For every integer `N >= 2^29`, `C_B2(f_N) >= 3N - 252 N^(4/5)`.

*Proof.* `N` is admissible (Lemma P). By F2.4 fix a monic irreducible
`q^0` of degree `k`; `T(q^0)` holds (F2.3). Let `C` be any B2 circuit
computing `f_N` on its `N` actual input bits. Replace each of the first `k`
input nodes by the constant node `q^0_j`. No gate is added or removed, and
the result computes `g_{q^0}` on the remaining `m` inputs. Let
`d* = floor(Lambda) + 1`. By Lemma 4.1, `g_{q^0}` is not constant on any
affine subspace of dimension `>= d*`, so by Theorem E

    |C| >= 3m - 4d* >= 3(N - k) - 4 Lambda - 4
        >= 3N - 48 N^(3/5) - 200 N^(4/5) - 4 >= 3N - 252 N^(4/5),

using Lemma P(a),(d). QED.

*Quantifiers.* Given `eps > 0`, let `N_eps = max(2^29, ceil((252/eps)^5))`.
For every integer `N >= N_eps`, `252 N^(4/5) <= eps N`, so
`C_B2(f_N) >= (3 - eps) N`. Hence
`liminf_{N -> infinity} C_B2(f_N)/N >= 3` over **all** integer lengths, for
this one fixed family, with the remainder `252 N^(4/5) = o(N)` (no linear
loss). The full guarantee is `g(N) = 3 N^1 (log2 N)^0 (log2 log2 N)^0`, and
the finite projection `(3, 1)` is exactly that expression. No larger
coefficient or exponent is claimed.

*Model and accounting.* The circuits are arbitrary B2 circuits on the `N`
actual input bits `w = (q, x)`; `N` includes the `k = O(N^(3/5))` bits of
`q`, and the zero padding is inside the definition of `f_N`, not extra
inputs. Every gate is counted. Constant nodes, repeated pins, degenerate
and negation gates, and a negated or constant output are all handled by
Lemmas E1 and E2, which only remove charged gates or replace one charged
object by another of equal cost. Label nodes cost 1 each in the measure, so
no uncharged basis element is introduced: the final inequality concerns the
gate count of the original circuit `C`. The bound holds at every integer
length `N >= 2^29`, not on a subsequence.
