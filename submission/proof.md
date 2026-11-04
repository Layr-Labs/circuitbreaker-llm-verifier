# Coefficient/exponent pair (1, 1) for the fixed total AND family

This proof is self-contained: it uses no external theorem, only the pinned
model conventions from `rules.md` and elementary counting. It claims a
finite score of [1,1], i.e. `C_B2(f_n) >= n - 1` at every length `n >= 1`.

## Fixed total family

For `n >= 1`, let `f_n : {0,1}^n -> {0,1}` be the AND of all `n` input bits.
For `n = 0`, set `f_0(empty) = 1` (the empty conjunction is vacuously true).
This defines one total family at every length, with no exceptions.

The language is `L = {x in {0,1}* : f_|x|(x) = 1}`.

## Uniform P and NP membership

One fixed deterministic algorithm decides `L`: on input `x` of length `n`,
scan the `n` bits once and accept iff every bit is 1 (accept the empty input).
This runs in `O(n)` time, so the family is uniform polynomial-time.

NP membership uses the fixed polynomial-time verifier `V(x, w)` that ignores
`w`, scans the bits of `x`, and outputs 1 iff every bit of `x` is 1, with the
fixed polynomial witness bound `W(n) = 0`. Both directions hold: if `x in L`
then `V(x, empty) = 1` with `|empty| = 0 <= W(|x|)`; if `V(x, w) = 1` for some
`w` with `|w| <= W(|x|)` then every bit of `x` is 1, so `x in L`. Per
`rules.md`, a fixed polynomial-time evaluator with an empty witness suffices.
The same language, verifier, and witness bound are used throughout.

## Circuit model: native B2, no transfer

The lower bound is proved directly in the pinned Li-Yang general B2 model of
`rules.md` (deterministic acyclic circuits, one output, two input pins per
gate, all 16 binary truth tables, unrestricted depth and fan-out, exact
computation on every input, every retained gate counted once). No transfer
from another model is needed. There are no blocks, field encodings, or
padding: the circuit input is exactly the `n` bits of `x`, so `n = |x|`
is the actual input bit length and there is no linear loss. Witness bits
are not circuit inputs.

## Lemma (variable dependence): a g-gate circuit depends on at most g+1 inputs

Say a signal in a circuit *depends on* input variable `x_i` if there exist
two assignments to the inputs, differing only at `x_i`, that give different
values to that signal. A bare-input feed depends on exactly its own variable
(one variable); a constant feed depends on none.

For a gate `G`, let `cone(G)` be the set of gates that feed `G`, including
`G` itself (all gates with a directed path to `G`). Claim: the output signal
of `G` depends on at most `|cone(G)| + 1` input variables.

Proof by strong induction on `|cone(G)| >= 1`. The two pins of `G` are each
fed by a signal that is either a bare input, a constant, or the output of
some gate `H` with `cone(H) subset of cone(G) \ {G}` (acyclicity: every path
into a pin-feed through another gate extends to `G`, and no path from `G`
reaches its own pin-feeds). The output of `G` is a fixed Boolean function of
its two pin values, so its dependency set is contained in the union of the
two pin-feed dependency sets.

Bound each pin-feed: a bare input contributes at most 1 variable; a constant
contributes 0; a gate `H` contributes at most `|cone(H)| + 1` by the induction
hypothesis, since `|cone(H)| <= |cone(G)| - 1`. In the worst case both pins
are fed by gates `H1, H2` with disjoint cones, giving at most
`(|cone(H1)| + 1) + (|cone(H2)| + 1) <= (|cone(G)| - 1) + 2 = |cone(G)| + 1`
variables, using `|cone(H1)| + |cone(H2)| <= |cone(G)| - 1` for disjoint
subsets of `cone(G) \ {G}`. Shared subgraphs (overlapping cones) and repeated
pins (both pins fed by the same signal) only shrink the union, and constant
pins contribute nothing, so all degenerate cases satisfy the same bound.
This closes the induction.

Consequences used below:

- Every gate in a circuit with `g` gates depends on at most `g + 1`
  variables, since `|cone(G)| <= g`.
- If the circuit output selects a bare input directly (in particular any
  `g = 0` circuit), it depends on exactly 1 variable, which is `<= g + 1`.
- Constant nodes, unary behavior (both pins carrying the same signal),
  negated projections, and shared subcircuits are all covered above; no
  substitution or normalization step is performed, and every retained gate
  is counted once in `g`. No direct-sum additivity is assumed: the argument
  upper-bounds a dependency set by a union, which holds under arbitrary
  sharing.

## The AND function depends on all n inputs

For `n >= 1` and each position `i`, the all-ones input gives `f_n = 1`
while flipping only bit `i` to 0 gives `f_n = 0`. Hence `f_n` depends on
every one of its `n` input variables.

## Main bound

Let a B2 circuit with `g` gates compute `f_n` exactly, `n >= 1`. Its output
signal equals `f_n`, so it depends on all `n` variables. By the Lemma it
depends on at most `g + 1`. Therefore `g + 1 >= n`, i.e. `g >= n - 1`.
Since the circuit was arbitrary:

    C_B2(f_n) >= n - 1 = (1)*n^1 - 1, for every integer n >= 1.

With `r(n) = 1`, `r(n)/n^1 -> 0`, so this is a `c*n^p - r(n)` bound with
`c = 1`, `p = 1`. It implies `liminf C_B2(f_n)/n >= 1`. Explicitly, for
every `eps > 0`, taking `N_eps = ceil(1/eps)`, every integer `n >= N_eps`
satisfies `n - 1 >= (1 - eps)*n` (equivalent to `eps*n >= 1`). All
sufficiently large integer lengths are covered, with the same family,
verifier, and witness bound throughout.

This establishes the declared full guarantee `polynomial_log` with
`c = 1/1`, `p = 1/1`, `q = 0`, `r = 0` and its finite projection `[1, 1]`.
