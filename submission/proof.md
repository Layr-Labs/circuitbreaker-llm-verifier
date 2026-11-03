# One total uniform family with C_B2(F_N) >= 4N - O(N^(4/5)): the declared theorem

## 0. The claim, in one place

Let `F` be the single total Boolean language defined in Section 1 below, on
actual input length `N` (so `N = |x|` is the number of input bits, not a
padded or shortened length). Then

    C_B2(F_N) >= 4N - 87 N^{4/5}                                   (0.1)

for every integer `N >= 1024`, where `C_B2(F_N)` is the minimum number of gates
in a full-`B2` circuit computing the length-`N` restriction of `F`. Since
`87 N^{4/5} = o(N)`,

    liminf over all integers N of  C_B2(F_N)/N  >=  4 .          (0.2)

The language `F` is decided by one fixed deterministic polynomial-time
algorithm and lies in `NP` with witness bound `0` (empty witnesses), so the
family is eligible under `rules.md`. The exponent is `p = 1`, and the
coefficient is `c = 4`. Both coordinates of `bound.json` belong to this one
theorem and this one language.

We assert no exponent above `1`, no endpoint, and no optimality of the
coefficient `4` over all families.

**Provenance.** The gate-elimination class theorem used here is proved in
full in `elimination.md`; the algebraic disperser is proved in full in
`disperser.md`. The construction is a *uniform* variant of the
Ben-Sasson-Kopparty affine disperser: instead of choosing an irreducible
polynomial by a construction algorithm, the irreducible polynomial is read
from a prefix of the input and *tested*, so that the whole family is decided
by one uniform deterministic algorithm with no advice. We prove everything we
use; the only external statements we rely on are the pinned model conventions
of `rules.md` (Li-Yang ECCC TR21-023) and the two papers named in Section 8.

## 1. The family: one evaluator on every bit string

The first bit of each field block is its constant coefficient. On an input
`x` of length `N`:

1. If `N = 0`, output `0`.
2. Let `B` be the smallest positive integer with `B^5 > 32 N^3`. Let `k` be
   the smallest prime at least `B`. Find `k` by deterministic trial division.
3. If `k >= N`, output `0`. Otherwise put `m = N - k`,
   `r = ceil(m/k)`, `t = ceil(sqrt(r))`. If `t >= r` or `k t < r (r+1)`,
   output `0`.
4. Read the first `k` bits `b_0,...,b_{k-1}` of `x` as the lower coefficients
   of a monic degree-`k` polynomial

       Q(Z) = Z^k + sum_{j=0}^{k-1} b_j Z^j   in F_2[Z].

   Check that `Z^{2^k} = Z` modulo `Q` and that `gcd(Q, Z^2 + Z) = 1`. If
   either check fails, output `0`.
5. Let `K = F_2[Z]/(Q)`. Append `r k - m` zero bits to the remaining `m`
   payload bits and interpret the resulting `r k` bits as `r` blocks
   `a_1,...,a_r` of `k` bits each. The first bit of a block is its constant
   coefficient. This is an injective `F_2`-linear map of the payload onto
   `K^r` (Section 5.2).
6. For `i = 1,...,r` set `w_i = r + 2 - i` and `y_i = a_i^{2^{w_i} - 1}`.
   Output the constant coefficient, in the basis `1, Z, ..., Z^{k-1}` of
   `K`, of the elementary symmetric polynomial of order `t` in `y_1,...,y_r`:

       F_N(x) = pi ( sum_{|I| = t} prod_{i in I} y_i ) ,

   where `pi : K -> F_2` is the constant-coefficient functional.

This specifies every choice and every input length. `Q` is part of `x`; it is
not advice, and the evaluator does not search for a hard prefix.

### 1.1 Evaluation in polynomial time

The elementary symmetric polynomial of order `t` is computed without
division, which matters because the characteristic is two. Initialize
`E_0 = 1` and `E_1 = ... = E_t = 0`. For `i = 1,...,r` update

    E_j <- E_j + y_i E_{j-1}     for j = min(i,t), t-1, ..., 1, in that order.

*Proof of the recurrence.* Let `E_j^{(i)}` be the sum of products over all
`j`-element subsets of `{1,...,i}`. Splitting those subsets according to
whether they contain `i` gives
`E_j^{(i)} = E_j^{(i-1)} + y_i E_{j-1}^{(i-1)}`, with `E_0^{(i)} = 1`.
The descending order of `j` ensures the right-hand side uses the values from
stage `i-1`. QED

Each `y_i` is computed by binary powering, using `O(w_i)` field
multiplications, since `w_i <= r+1` and the exponent `2^{w_i}-1` has `w_i`
binary digits. Field multiplication in the polynomial basis costs `O(k^2)`
coefficient operations by schoolbook multiplication. The symmetric recurrence
costs `O(r t)` field operations, so the total is `O(r^2 k^2 + r t k^2)` field
operations, each on `k`-bit values: polynomial in `N`.

## 2. The prefix test is exact, and hard prefixes exist

We must know (i) that the test in step 4 of Section 1 accepts exactly the
irreducible degree-`k` polynomials, and (ii) that at least one accepted prefix
exists for each prime `k`, so that the restriction used in Section 5 is a
restriction onto a nonconstant function.

**Fact A.** In an algebraic closure of `F_2`, the polynomial `Z^{2^k} - Z` has
`2^k` distinct roots, because its derivative is `1`. Its roots are closed
under addition, multiplication and inversion of nonzero elements, hence form
the field `K` of size `2^k`. The Frobenius map `z |-> z^2` permutes them and
its `k`-th iterate is the identity on them.

**Fact B.** A Frobenius orbit of length `d` gives an irreducible degree-`d`
polynomial over `F_2`: the product of the linear factors over the orbit is
Frobenius-fixed, while a proper factor over `F_2` would require a proper
invariant subset of the orbit. Conversely, an irreducible factor's roots form
one orbit. Thus the degrees of the irreducible factors of `Z^{2^k}-Z` divide
`k`.

**Fact C (existence).** Since `k` is prime, the orbits have size `1` or `k`.
The only `F_2`-elements fixed by Frobenius are `0` and `1`, and `2^k > 2`, so
there is an orbit of size `k`. Its product gives an irreducible degree-`k`
polynomial. Hence an irreducible degree-`k` polynomial exists for every prime
`k`. No algorithm selecting it is needed.

**Fact D (exactness of the test).** The first test says `Q` divides the
squarefree polynomial `Z^{2^k} - Z`. The second excludes the factors `Z` and
`Z + 1`, i.e. it says `gcd(Q, Z^2+Z) = 1`. By Fact B, every irreducible factor
of `Q` then has degree `k`; since `Q` has degree `k` and is monic, `Q` itself
is irreducible. Conversely, an irreducible degree-`k` polynomial `Q` divides
`Z^{2^k}-Z` (Fact B) and is coprime to `Z^2+Z` (its roots are not `0` or `1`),
so it passes both tests. The test is therefore exact.

**Fact E (polynomial-time computability of the test).** The first test is
decided by `k` squarings modulo `Q`, starting from `Z` and repeatedly
squaring: `Z^{2^k} = Z` is checked by computing `Z^{2^k}` in `O(k)` modular
squarings and comparing with `Z`. The second is a polynomial Euclid
computation with `gcd(Q, Z^2+Z)`. Both have polynomial bit cost in `k`.

## 3. Uniform P and NP membership with empty witnesses

**Polynomial time of the whole evaluator.** Bertrand's postulate gives
`k < 2B`, and `B^5 > 32 N^3` gives `B > (32 N^3)^{1/5} = 2 N^{3/5}`, so

    2 N^{3/5} < k < 4 N^{3/5} + 2 .                               (3.1)

The prime scan and the trial division that certifies primality are
polynomial in `N`. `B` is computed by integer binary search with fifth-power
comparisons, hence exactly and in polynomial time. All the integer,
polynomial and field operations of Section 1 have polynomial bit cost. A
coarse bound `O((N+1)^6)` suffices for the whole algorithm, including on
rejected inputs.

**NP membership, both directions.** Let `A` be the fixed evaluator of
Section 1, let `L = { x : A(x) = 1 }`, and let `W(N) = 0`. Define the
verifier

    V(x, w) = 1   iff  w is empty and A(x) = 1 .

Then `x in L` iff there is a `w` with `|w| <= W(|x|) = 0` and `V(x,w) = 1`:
if `x in L` then `A(x) = 1` and the empty witness is accepted; if a permitted
witness is accepted then `A(x) = 1`, so `x in L`. The verifier is one fixed
deterministic polynomial-time algorithm, and the witness bound is a fixed
polynomial (zero). This proves `L in NP` (in fact `L in P`).

The witness is *not* an input to the circuits that the lower bound concerns:
the circuits compute the length-`N` function `F_N` on the `N` bits of `x`
alone.

## 4. The guards hold at all sufficiently large integer lengths

We have from (3.1) `k = Theta(N^{3/5})`, so `m = N - k = Theta(N)`,
`r = ceil(m/k) = Theta(N^{2/5})` and `t = ceil(sqrt(r)) = Theta(N^{1/5})`.

* `k < N` and `t < r` hold eventually, since `k = O(N^{3/5}) = o(N)` and
  `t = O(N^{1/5}) = o(N) = o(r)`.
* `r < N/k + 1 < (1/2) N^{2/5} + 1` for large `N`, because `k > 2 N^{3/5}`
  and `m < N`.
* The guard `k t >= r (r+1)`: using `t >= sqrt(r)`,

      r (r+1)/t <= sqrt(r) (r+1)
                <= sqrt((1/2) N^{2/5} + 1) * ((1/2) N^{2/5} + 2) .

  Dividing the right-hand side by `N^{3/5}` and letting `N` grow gives the
  limit `1/(2 sqrt(2)) < 1`, whereas `k/N^{3/5} > 2`. Hence
  `k t >= r (r+1)` for all sufficiently large `N`.

All these estimates hold for **every** sufficiently large integer `N`, not
only for primes or a sparse subsequence. We choose `N_0` so that they hold for
all `N >= N_0`, and set `N_0 >= 1024`.

## 5. The hard restriction and its precise loss

### 5.1 Fixing a prefix

For a fixed sufficiently large `N`, fix **any** irreducible degree-`k`
polynomial `Q_0` over `F_2` (it exists, by Fact C) and fix the first `k` bits
of `x` to the coefficients of `Q_0`. The remaining function `G_{Q_0}` has
`m = N - k` inputs. Put

    d_N = floor(k r / t) + k t + r (r+1) + 1 .                   (5.1)

By `disperser.md` (D), with `pi` the constant-coefficient functional, the
function

    (a_1,...,a_r) |-> pi ( sum_{|I|=t} prod_{i in I} a_i^{2^{w_i}-1} ) ,  (5.2)

is nonconstant on every affine subspace of `K^r` of dimension strictly
greater than `k r / t + k t + r (r+1)`, provided `k` is prime, `1 <= t < r`
and `k t >= r (r+1)`. All three hold by Section 4 and by the choice of prime
`k`. The payload embedding of Section 1, step 5, is injective and
`F_2`-linear, and `pi` is nonzero because `pi(1) = 1`.

Since `k t` and `r (r+1)` are integers, the smallest *integer* dimension
strictly greater than `k r / t + k t + r (r+1)` is

    floor(k r / t) + k t + r (r+1) + 1 = d_N .                  (5.2')

Hence `G_{Q_0}` is an affine disperser for dimension `d_N` on its `m` actual
input bits: it is nonconstant on every affine subspace of dimension `d_N`, and
therefore on every affine subspace of any larger dimension (a subspace of
dimension `>= d_N` contains one of dimension exactly `d_N`).

**Off-by-one remark.** The threshold `d_N` is the smallest integer dimension
that BK's hypothesis covers; `d_N - 1` is *not* covered, because
`floor(k r / t) + k t + r (r+1) <= k r / t + k t + r (r+1)`. We use `d`
exactly, with no slack, so the bound below is `4 m - 10 d_N` and not
`4 m - 10 (d_N - 1)`.

The zero padding in step 5 preserves the dimension of every payload affine
subspace: padding is the injective `F_2`-linear map that appends zero
coordinates, which sends an affine subspace injectively to an affine subspace
of the same dimension. It adds no circuit inputs.

### 5.2 The embedding is injective and linear

The `r k` padded payload bits are cut into `r` blocks of `k` bits, and each
block is read as an element of `K` in the basis `1, Z, ..., Z^{k-1}`. Reading
bits as coefficients in a fixed basis is an `F_2`-linear isomorphism
`F_2^k -> K`, so the block map is `F_2`-linear; concatenating blocks is
`F_2`-linear; and the map is injective because the bases are bases.

### 5.3 The loss

Let `C` be any `B2` circuit computing `F_N`. Fixing the first `k` inputs of
`C` to the prefix of `Q_0` yields a circuit computing `G_{Q_0}`, with **no
additional gates**: substitution of constants into a circuit's input pins
cannot increase its gate count. By the class theorem of `elimination.md`,
applied to `G_{Q_0}` on its `m` actual input bits,

    |C| >= 4 m - 10 d_N .                                        (5.3)

Since `m = N - k`,

    C_B2(F_N) >= 4 (N - k) - 10 d_N .                            (5.4)

Now `d_N = floor(k r/t) + k t + r(r+1) + 1` with `k = O(N^{3/5})`,
`r = O(N^{2/5})`, `t = Theta(N^{1/5})`, so

    k r / t = O(N^{3/5} N^{2/5} N^{-1/5}) = O(N^{4/5}) ,
    k t     = O(N^{3/5} N^{1/5})      = O(N^{4/5}) ,
    r (r+1) = O(N^{4/5}) .

Hence `d_N = O(N^{4/5})` and `k = O(N^{3/5}) = O(N^{4/5})`, giving

    C_B2(F_N) >= 4N - 4k - 10 d_N = 4N - O(N^{4/5}).              (5.5)

### 5.4 An explicit remainder constant

We now bound `4k + 10 d_N` by an explicit multiple of `N^{4/5}`, so that the
displayed inequality is valid for every `N >= 1024` and not merely
asymptotically. Throughout, `N >= 1024`, so `N^{1/5} >= 4`, `N^{2/5} >= 16`,
`N^{3/5} >= 64` and `N^{4/5} >= 256`.

**(i) The field degree.** From (3.1), `k < 4 N^{3/5} + 2`. Since
`N^{3/5} >= 64 > 4`, we have `2 <= (1/2) N^{3/5}`, so

    k <= 4 N^{3/5} + (1/2) N^{3/5} = (9/2) N^{3/5}.                 (5.6)

**(ii) The block count.** `r = ceil(m/k) <= m/k + 1 < N/k + 1`, and
`k > 2 N^{3/5}` gives `N/k < (1/2) N^{2/5}`. Also
`1 <= N^{2/5}/16`, so

    r <= (1/2) N^{2/5} + 1 <= (9/16) N^{2/5}.                    (5.7)

**(iii) The symmetric order.** `t = ceil(sqrt r) <= sqrt r + 1`, and
`sqrt r <= (3/4) N^{1/5}` because `r <= (9/16) N^{2/5}` and
`sqrt(9/16) = 3/4`. Since `1 <= N^{1/5}/4`, we get

    t <= (3/4) N^{1/5} + 1 <= N^{1/5}.                            (5.8)

**(iv) The disperser dimension.** Now, term by term,

    floor(k r / t) <= k r / t <= k sqrt(r) <= (9/2)(3/4) N^{4/5}
                     = (27/8) N^{4/5} ,                          (5.9)

    k t <= (9/2) N^{3/5} * N^{1/5} = (9/2) N^{4/5} ,             (5.10)

    r (r+1) <= (9/16) N^{2/5} * ( (9/16) N^{2/5} + 1 )
             <= (9/16)(9/16 + 1/16) N^{4/5} = (45/128) N^{4/5} ,   (5.11)

    1 <= (1/256) N^{4/5} .                                       (5.12)

Adding (5.9)-(5.12),

    d_N <= (27/8 + 9/2 + 45/128 + 1/256) N^{4/5}
        = (2107/256) N^{4/5} .                                   (5.13)

**(v) The total loss.** From (5.6), `4k <= 18 N^{3/5}`, and
`N^{3/5} <= (1/4) N^{4/5}` because `N^{1/5} >= 4`. Hence `4k <= (9/2)
N^{4/5}`, and with (5.13),

    4 k + 10 d_N <= ( 9/2 + 10 * 2107/256 ) N^{4/5}
                 = ( 1152/256 + 21070/256 ) N^{4/5}
                 = (22222/256) N^{4/5}
                 = (11111/128) N^{4/5}
                 < 87 N^{4/5} .                                  (5.14)

Since `11111/128 = 86.8046875 < 87`, we conclude

    C_B2(F_N) >= 4N - 87 N^{4/5}     for every integer N >= 1024 .  (5.15)

This is the inequality (0.1). We make no claim that `87` is optimal; only that it is a valid explicit
constant for every `N >= 1024`, and that the remainder is `o(N)`.

### 5.5 The choice of `Q_0` is not advice

For each `N` we chose *one* irreducible `Q_0` to carry out the lower-bound
restriction. This is not advice to the evaluator: the evaluator reads the
prefix from the input and tests it, so it computes the same function for
*every* prefix, including prefixes that are rejected. A circuit for `F_N`
must be correct for all inputs, hence in particular for the prefix of `Q_0`;
restricting to that prefix gives a circuit for `G_{Q_0}` at no extra cost.
The language, the verifier and the witness bound are the same for every `N`.

## 6. The quantifiers, the model transfer, and the two coordinates

**The guarantee.** For every real `epsilon > 0` there is an `N_epsilon` such
that for every integer `N >= N_epsilon`,

    C_B2(F_N) >= (4 - epsilon) N .                               (6.1)

Indeed, take `N_epsilon` large enough that the guards of Section 4 hold and
that the remainder in (5.5) satisfies `O(N^{4/5}) <= epsilon N` for all
`N >= N_epsilon`; this is possible because `N^{4/5}/N = N^{-1/5} -> 0`. Then
(5.5) gives (6.1). Equivalently,

    liminf over all integers N of C_B2(F_N)/N >= 4 ,

which is the `polynomial_log` guarantee of `bound.json` with
`c = 4`, `p = 1`, `q = 0`, `r = 0`; and the declared `finite` projection is
the same pair `(4, 1)`.

**Model transfer.** The class theorem of `elimination.md` is stated and
proved in exactly the model that `rules.md` pins: Li-Yang's general `B2`
model, i.e. deterministic acyclic circuits, one output, two input pins per
gate, all sixteen binary truth tables, unrestricted depth and fan-out, exact
computation on every input, nonuniform circuits permitted. We spell out the
points that `rules.md` singles out, and where each is handled:

* *Inputs are `x` alone, not NP witnesses.* The circuits in (5.3) compute
  `G_{Q_0}` on the `m` payload bits; the witness bound is `0` and no witness
  bit is ever an input.
* *`n` is the actual input bit length.* `N = |x|`. The prefix of `k` bits is
  part of `x` and is charged in full: (5.4) uses `m = N - k` for the
  disperser argument and adds back `4k`, and `d_N` is a function of the
  payload parameters. No padding is counted as free input length; the
  zero-padding inside the payload is counted as part of the `r k` block
  structure and its dimension is preserved, not inflated.
* *All 16 truth tables; a retained negated projection costs one.* The
  sixteen tables are enumerated in `elimination.md` Section 1.1; constant and
  unary tables are handled by the deletion (D1), and a negated projection
  that is retained is one gate, as the counting of `s` records.
* *Constant nodes during substitutions cost zero and must be handled by
  Li-Yang normalization.* Rerouting to a constant is free (it deletes the
  object) and is exactly the normalization of `rules.md`; `elimination.md`
  Sections 3 and 5.3 state each deletion's effect on the potential.
* *Repeated pins and shared subgraphs are allowed.* Degree counts distinct
  children, not pins; a gate reading the same value twice is one child. The
  deletions of `elimination.md` count *distinct* deleted objects, so shared
  descendants are not double-counted, and an object removed early by
  propagation is counted once.
* *A normalized proof must cover constant, unary and terminal-output cases
  without adding an uncharged computational basis.* The terminal-output
  cases are handled by the explicit lemma `elimination.md` Section 4
  (Lemma 4.1), which constructs a `d`-dimensional affine subspace on which the
  circuit is constant; it does not presume a free extra output gate.
* *Formula, monotone, bounded-depth, XOR-free, multiple-output, randomized or
  average-case bounds require a proved transfer.* None of these is used: the
  bound is proved directly for the full `B2` model with one output, arbitrary
  depth, arbitrary fan-out, unrestricted gates and exact computation.
* *No direct-sum additivity under unrestricted sharing.* The only use of a
  decomposition into blocks is in `disperser.md` Section 4, which is ordinary
  linear Gaussian elimination on an affine subspace, and in the selection of
  `t` blocks in Section 5 of that file, which is a choice of indices with a
  dimension count. No circuit complexity is additively combined.
* *Do not identify a potential value with a gate count.* `elimination.md`
  Section 6 ends with an explicit paragraph on why `s` bounds the number of
  original gates: rows are introduced one-for-one from charged gates, roles
  are exchanged without adding computation, and restrictions never add a
  gate.

**The two coordinates belong to this one theorem.** The coefficient `4` and
the exponent `1` come from the single inequality (5.4)-(5.5) about the single
family `F` at actual input length `N`. They are not supplied by two
unrelated results.

## 7. What is not claimed

* We claim no exponent above `1`, no lane endpoint, no global endpoint, and
  no statement about `P` versus `NP`.
* We claim no coefficient above `4`. `optimality.md` proves, by exact
  rational arithmetic, that within the family of potentials
  `V = s + beta q + alpha I` and the row estimates of `elimination.md`
  Section 7, no choice of weights yields a leading coefficient above `4`.
  That is an optimisation bound for these estimates, not a barrier to any
  other proof.
* We do not claim that the constant `87` in (0.1) is sharp, nor that the
  disperser dimension `d_N` is optimal, nor that the construction is optimal
  among affine dispersers.
* The result is a written proof, not a formal proof-assistant certificate.
* No checker output, reviewer verdict, prior acceptance, or comparison hint
  is used as a premise of any step. Passing a checker certifies arithmetic
  only.

## 8. Sources, versions and how their hypotheses are met

* **Mishra, O. S., and Santhanam, R.** *An Explicit 4n - O(log(n)) Size Lower
  Bound against General Boolean Circuits.* ECCC TR26-260, **Revision 1**,
  3 November 2026. `https://eccc.weizmann.ac.il/report/2026/260/revision/1/download/`.
  SHA256 `a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921`.
  Used for: the gate-elimination class theorem `4n - 10d` for affine
  dispersers (its Theorem 4.1), the fair semi-circuit formalism, the
  parameter operations, the potential `W = s + 4q + 6I`, the
  one-iteration Lemma 4.2, and the terminal-output Lemma 2.2. Its
  hypotheses are met as follows: its Theorem 4.1 assumes `f` is an affine
  disperser for dimension `d` computable in polynomial time; our
  `G_{Q_0}` is an affine disperser for dimension `d_N - 1` by `disperser.md`
  (D), and it is computable in polynomial time because it is the
  prefix-restriction of the uniform evaluator of Section 1. The whole gate
  argument is reproduced in `elimination.md`; nothing here depends on an
  unquoted premise of that paper. Note that Revision 0 of the same report
  (19 October 2026) proves only `3.5n - o(n)`; Revision 1 is the version used.
* **Ben-Sasson, E., and Kopparty, S.** *Affine Dispersers from Subspace
  Polynomials.* Author preprint, 12 March 2010.
  `https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf`. SHA256
  `922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf`.
  Used for: the disperser construction and its dimension bound (their
  Theorem 7.3 and Theorem 2.2), and the image-polynomial algebra (their
  Lemmas 3.7, 3.10, 3.13, 6.1, 7.2 and Claim 7.4). Their hypotheses are met
  as follows: their Theorem 7.3 assumes `k` prime, `t < r` and
  `k t >= r(r+1)`; Section 4 of this file guarantees all three for all
  sufficiently large `N`, and `pi` is nonzero because `pi(1) = 1`. The
  argument is reproduced in `disperser.md`.
* **Li, J., and Yang, T.** *3.1n - o(n) Circuit Lower Bounds for Explicit
  Functions.* ECCC TR21-023, front-page date 20 February 2021.
  `https://eccc.weizmann.ac.il/report/2021/023/download/`. SHA256
  `4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062`.
  Used for: the pinned model conventions required by `rules.md` (Sections
  1.1, 2.1, 2.4 footnotes 3-4, and 3.3). Their circuit theorem and their
  linear program are *not* used; this submission does not invoke them.

**Not used as premises.** The starting repository seed (Li-Yang's `31/10`
coefficient) is not used: the coefficient here is `4`, from a different and
later argument. No claim, review verdict, acceptance, checker result or
comparison key from any other submission is relied upon.

## 9. Evidence index

* `elimination.md` — the class theorem `4n - 10d` for affine dispersers, in
  full, with its potential and all case estimates.
* `disperser.md` — the algebraic disperser, in full, in the binary
  prime-degree instance used here.
* `optimality.md` — an unscored lemma: the coefficient `4` is optimal for
  this potential family, by exact rational arithmetic.
* `proof.md` — this file: the family, the evaluator, the NP verifier, the
  guards, the restriction and its loss, the quantifiers, and the model
  transfer.
* `bound.json` — the canonical declaration of the joint theorem.
