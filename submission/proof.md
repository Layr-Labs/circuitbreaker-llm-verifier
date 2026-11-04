# Coefficient 1 at exponent 1 for the OR family: a self-contained proof

This proof uses no external premise. Every statement needed for the scored
conclusion is proved below from the definitions in `rules.md`. No result
from any cited paper is used. The 31/10 literature seed from the baseline
is prior work; it is not claimed here, for the reasons in Section 6.

## 1. The family: every length defined

Let `n = |x|` be the actual input bit length. Define the total Boolean
family `f_n` on `{0,1}^n` by

- `f_0` (the empty input) `= 0`;
- for `n >= 1`, `f_n(x) = 1` iff at least one input bit `x_i` is `1`,
  i.e. the `n`-ary OR of the input bits.

There is no encoding, block, field, padding, or divisibility convention:
the input bits are exactly the function's arguments at every length. The
same definition is used at every length, including the small lengths
`n = 0` and `n = 1`, which are handled explicitly in Sections 2 and 5.

## 2. Uniform polynomial time, hence NP membership

One fixed deterministic evaluator `A` computes the whole family: on input
`x` of length `n`, scan the bits once from left to right, output `1` on
the first `1` seen and `0` if the scan finishes without seeing a `1`
(for `n = 0` the scan is empty and the output is `0`). The running time
is linear in `n`, hence polynomial. The same algorithm, with no
per-length advice and no length-dependent case split beyond the
definition above, works at all lengths.

Define the verifier `V(x, w) = 1` iff `w` is empty and `A(x) = 1`, and
`V(x, w) = 0` otherwise, with witness bound `W(n) = 0` for all `n`.
`V` runs in polynomial time because `A` does. Both directions hold:

- If `x` is in `L` (i.e. `f_n(x) = 1`), the empty witness `w` satisfies
  `V(x, w) = 1`, since `A(x) = f_n(x) = 1`.
- If some `w` with `|w| <= 0` satisfies `V(x, w) = 1`, then `w` is empty
  and `A(x) = 1`, so `f_n(x) = 1` and `x` is in `L`.

Thus `L` is in P, hence in NP, with one fixed verifier and witness bound.
Witness bits are not circuit inputs: the circuits below receive `x`
alone, exactly as `rules.md` requires.

## 3. The circuit model and what is counted

We use the Li-Yang general B2 model from `rules.md`: deterministic
acyclic circuits, one output, two input pins per gate, all 16 binary
Boolean truth tables, unrestricted depth and fan-out, exact computation
on every input. Circuits may be nonuniform; the argument below is made
separately at each fixed length `n`, so nonuniformity is allowed.
Inputs, wires, and output selection cost zero. Each gate present in
the circuit counts once, including trivial or degenerate gates.
Constant nodes cost zero. Negations can be absorbed into descendant
truth tables; a retained negated-projection gate costs one like any
other gate. None of these conventions is used to hide a gate: the
counting argument in Section 4 counts gates directly.

Fix `n >= 2` and a B2 circuit `C` that computes `f_n` exactly on every
input. Let `g_total` be its total number of gates. We prove
`g_total >= n - 1`. The proof covers the constant, unary, and
terminal-output cases explicitly:

- **Terminal-output case.** If the circuit output were a bare input
  wire `x_j`, the computed function would be the projection `x_j`,
  but `f_n(e_j) = 1` and `f_n(e_k) = 0` for `k != j` (where `e_i` is
  the `i`-th unit vector), while `f_n(0^n) = 0`; a projection cannot
  equal `f_n`. If the output were a constant node, the computed
  function would be constant, but `f_n` is nonconstant for `n >= 2`.
  Hence the output is produced by a gate, the output gate.
- **Trimming.** Let `G` be the set of gates that are ancestors of the
  output gate (including the output gate itself): a gate is in `G` if
  there is a directed path of wires from its output to the output
  gate. Removing every gate outside `G` does not change the computed
  function, and `|G| <= g_total`. It suffices to prove `|G| >= n - 1`.
  Write `g = |G|`.
- **Constants and unary gates.** Constant nodes that feed gates of `G`
  are vertices of the graph below; they are never charged as gates,
  which can only make the bound harder to prove, never easier. A
  unary-behaving gate (for example one whose output ignores a pin, or
  a negated projection) is still a gate with two pins in the model
  and is counted once in `g`; the edge bound below uses at most two
  incident edges per gate, which holds for every gate whatever its
  truth table.

## 4. Lower bound: every input is essential, and the graph is connected

**Essentiality.** For each `i` in `1..n`, `f_n(0^n) = 0` while
`f_n(e_i) = 1`. These two inputs differ only in bit `i`, so the
computed function genuinely depends on every input bit: no input bit
is dummy. Consequently every input node `x_i` feeds (through wires,
possibly via other gates) at least one gate of `G`: if some `x_i`
reached no gate of `G`, the output value could not depend on `x_i`,
contradicting essentiality. So all `n` input nodes occur in the
trimmed circuit.

**The graph.** Form the underlying undirected graph `H` whose vertices
are: the `n` input nodes, the `g` gates of `G`, and the `k` constant
nodes that feed gates of `G` (so `k >= 0`). Two vertices are joined by
an (undirected) edge when a wire runs between them, i.e. when one is a
gate of `G` and the other is the source of one of that gate's two
input pins (an input node, a constant node, or another gate's output).
Each distinct edge is charged to the unique gate whose input pin it
feeds: acyclicity forbids the same undirected edge arising from two
different gates' pins, and repeated pins to the same source contribute
one distinct edge. A gate's fan-out as a source creates no additional
edges under this charging, since each outgoing wire is charged to its
consumer gate's pin. Each gate has exactly two pins, so with `E` the
number of distinct edges of `H`,

    E <= 2g.

**Connectivity.** Every vertex of `H` is joined to the output gate by
an undirected path: each gate of `G` reaches the output gate by a
directed path by the definition of `G`; each input or constant node
feeds a gate of `G` (shown above for inputs; constants are included
only when they feed a gate of `G`). Forgetting edge directions, `H`
is therefore connected. A connected graph on `V` vertices has at
least `V - 1` edges. Here `V = n + g + k`, so

    n + g + k - 1 <= E <= 2g,

hence `g >= n + k - 1 >= n - 1`. Since `g_total >= g`,

    C_B2(f_n) >= n - 1     for every n >= 2.

No substitution, normalization, or case split on gate truth tables
was needed: the bound follows from essentiality, the two-pin limit,
and connectivity alone, and it holds under unrestricted sharing,
fan-out, and depth.

## 5. Upper bound and exact size; small lengths

For `n >= 2`, a binary OR tree computes `f_n`: pair the inputs, OR
each pair with one OR gate (OR is one of the 16 permitted truth
tables), and repeat on the results until one value remains. A full
binary tree with `n` leaves has exactly `n - 1` internal nodes, so
`C_B2(f_n) <= n - 1`. Combined with Section 4,

    C_B2(f_n) = n - 1     for every integer n >= 2.

Small lengths, for completeness: `f_0 = 0` is constant, so
`C_B2(f_0) = 0` (a constant output, no gate charged); `f_1(x_1) = x_1`
is the identity, so `C_B2(f_1) = 0` (a bare input wire). These values
do not affect any asymptotic statement, because the quantifiers below
use thresholds `N >= 2`.

## 6. Asymptotics: the declared guarantee

From `C_B2(f_n) = n - 1` for `n >= 2`,

    C_B2(f_n) / n = 1 - 1/n  ->  1,

so `liminf_{n -> infty} C_B2(f_n) / n >= 1`. In the equivalent
epsilon form: for every `epsilon > 0`, let
`N_epsilon = max(2, ceil(1/epsilon))`. For every integer
`n >= N_epsilon`, `epsilon * n >= 1`, hence
`C_B2(f_n) = n - 1 >= (1 - epsilon) * n`. This is exactly the
guarantee `g(n) = 1 * n^1` with `c = 1`, `p = 1`, log and loglog
exponents `0`, holding at every sufficiently large integer length
for the one fixed family of Section 1. The coefficient and the
exponent belong to this same theorem and this same family; no two
lemmas are combined.

## 7. Why the 31/10 seed is not claimed here (prior work, not a premise)

The baseline submission claimed `(31/10, 1)` for a BK/Shoup family by
applying Theorem 1.1 of Li-Yang (ECCC TR21-023) and Theorem 2.2 of a
Bourgain-Katz preprint, with a retained source discrepancy between
that paper's Corollary 3.8 and Theorem 3.9 remainders. That argument
is not reproduced in this submission, and no part of it is used as a
premise above: the proof in Sections 1-6 is complete without it.
Under the supplied-source evidence policy in `rules.md`, a cited
theorem whose proof is not reproduced in the supplied evidence cannot
carry the scored conclusion, and the baseline proof itself records
that it "does not re-prove the paper's long case analysis". Retuning
the Li-Yang linear program alone cannot improve on `31/10` either:
the appendix below verifies that `31/10` is the optimum of that
unchanged LP. No stronger unconditional bound is proved in this
session, and none is claimed. The coefficient `1` proved here is a
complete, unconditional, fully evidenced theorem; it is weaker than
the literature seed's claimed number and stronger than no evidence.

## Appendix (unscored): the unchanged LP arithmetic

The following is pure rational arithmetic about the Li-Yang Theorem
4.1 linear program at the published point. It is a relevant lemma
only; it establishes no circuit bound by itself, because the
construction and case analysis behind the LP are not reproduced here.

Let `alpha_I = 48/5`, `alpha_Q = 9/5`, `alpha_phi = 1/5`,
`delta = 64/5`, and write `phi = alpha_phi`, `Q = alpha_Q`,
`u = delta - alpha_I`. The six expressions bounding `u` are

    alpha_I/3 = 16/5,
    2 - 2*phi + Q = 17/5,
    4 - 4*phi = 16/5,
    3 + phi = 16/5,
    5 - Q = 16/5,
    (5 - 2*phi + Q)/2 = 16/5,

each verified by direct fraction arithmetic: `alpha_I/3 = 16/5`;
`2 - 2/5 + 9/5 = 17/5`; `4 - 4/5 = 16/5`; `3 + 1/5 = 16/5`;
`5 - 9/5 = 16/5`; `(5 - 2/5 + 9/5)/2 = (32/5)/2 = 16/5`.
Moreover `u = 64/5 - 48/5 = 16/5`, and the admissibility conditions
`0 <= Q <= alpha_I` and `0 <= phi <= 1/2` hold. Weighting the
fourth, fifth, and sixth expressions by `7/10`, `1/10`, and `1/5`
(which sum to `1`) gives

    u <= (7/10)*(3 + phi) + (1/10)*(5 - Q) + (1/5)*((5 - 2*phi + Q)/2)
      = 31/10 + phi/2,

since the constant parts sum to `(21 + 5 + 5)/10 = 31/10` and the
`phi`/`Q` parts cancel to `phi/2`. Hence with
`c = delta - alpha_I - phi/2`, the LP optimum is `c = 31/10` at this
feasible point, and no retuning of this LP alone exceeds it. The
optional `lp-certificate` checker inputs in the baseline recorded
exactly this arithmetic; they are not used as evidence for any
circuit claim in this submission.
