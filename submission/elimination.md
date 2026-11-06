# A simple parameter-potential gate bound for affine dispersers

This file proves the gate-elimination lemma used by `proof.md`, written from
scratch. It is the classical flat-weight special case of the parameter
method: no quadratic pairs, no weight split, a single potential rate.

**Lemma E.** Let g be an affine disperser of dimension d on m input bits:
g is nonconstant on every affine subspace of F_2^m of dimension at least d,
and hence on every larger one. Then every full-B2 circuit computing g has
size at least 3m - 6d.

All circuits in this file are full-B2: fan-in-two gates with arbitrary
binary truth tables, unbounded fan-out, arbitrary depth, one output;
repeated pins and constant pins are permitted and every gate is counted
once. Negations on edges are absorbed into the tables of the reading
gates, which does not change any gate count.

## 1. States and the invariant

A **state** consists of:

* a positive integer M and pairwise distinct symbols u_1,...,u_M, the
  **parameters**;
* an injective affine map Phi from F_2^M to F_2^m;
* a linear system A*v = N*u + c over F_2, where A is a square,
  invertible matrix called the **auxiliary system**; its variables v
  (the **auxiliary values**) and the parameters u together satisfy that
  each row of the system contains at most three distinct nonconstant
  terms among {u_1,...,u_M} union {v};
* an acyclic set of B2 gates over the node set {u} union {v} plus
  constant nodes, with a designated **output** that is a constant, a
  literal of one node, or (equivalently, by table absorption) one gate.

Write s = (number of rows of A) + (number of gates). The **invariant** of a
state is: on every assignment of the parameters, the designated output
evaluates to g(Phi(u)).

The **initial state** of a circuit C computing g on inputs x_1,...,x_m is:
M = m, parameters the x_i, Phi the identity, A the empty invertible
matrix, the gates of C in their given order. The invariant then reads
output = g(x), which holds because C computes g. Initially
s = |C|, and at most m parameters are read by anything.

**Terminal values.** While a := M > d, the output of a valid state can
neither be a constant, nor an affine function of the parameters, nor the
literal of an AND-type gate on two parameters:

* The image of Phi is an affine subspace of dimension exactly M
  (injectivity of an affine map), and M >= d+1 > d in the loop, so g is
  nonconstant on it.
  A constant or affine-in-u output composed with Phi would make
  g(Phi(u)) constant or affine. If the output were constant, g would be
  constant on that subspace, impossible. If the output were a nonconstant
  affine function of u, it would have some parameter u_i with nonzero
  coefficient; restricting u_i to the affine function of the OTHER
  parameters given by the same affine form (a legal M2 restriction, not
  necessarily a constant) makes the output identically constant, while
  the image of the new Phi still contains an affine subspace of
  dimension M-1 >= d, again impossible.
* If the output is a literal of an AND-type gate G on parameters x,y, fix
  a forcing value of x (every nondegenerate AND-type table has a constant
  forcing input); the output becomes constant while M-1 >= d, the same
  contradiction. If G's pins are not both parameters, this paragraph does
  not apply and G is treated by the reduction of Section 4 instead.

These impossibilities are the only places the disperser property is used.

## 2. Legal moves

Each move below preserves the invariant, and we track exactly how the
potential H of Section 3 changes.

**(M1) Absorb an affine gate.** An XOR-type gate (table u+v+c) none of
whose ancestors is an AND-type gate can be moved into the system: append
its value as a new auxiliary value and its defining equation as a new row.
The row has at most three nonconstant terms (its two parents, each a
parameter, an auxiliary value, or a constant). The new auxiliary matrix is
block triangular over the old one with last diagonal entry 1, so it stays
invertible. One gate becomes one row: s is unchanged, M is unchanged, and
no degree (Section 3) increases, since the gate's readers now read an
auxiliary value whose own readers are exactly the gate's old readers plus
possibly nothing new. Gates are absorbed in topological order, so each
gate moves at most once per phase.

**(M2) Fix a parameter.** Restrict u_i to an affine function of the other
parameters (a constant in all our uses). Substitute it in Phi, in N, in c,
and in every gate pin that reads u_i. Then delete u_i from the parameter
list: M drops by 1. Gates that read u_i acquire a fixed or constant pin;
rows that mention u_i lose one term. No new node is created.

**(M3) Affine pivot.** Let x be a parameter and Q an auxiliary value whose
expression (as an affine form in u, obtained by solving the system) has a
nonzero coefficient on x. Exchange the roles of x and Q: declare Q a
parameter and x the auxiliary value taking over Q's node. Concretely, in
A*v = N*u + c the column of v_Q is replaced by the column that expressed
Q through u, moved to the left side, while the old parameter x appears on
the right wherever it did before, and the value computed at Q's node is
redefined to x (whose old readers keep reading the same node). The new
matrix differs from A in one column; its determinant is the Q-coordinate
of A^(-1) times the new column, which is nonzero by the choice of Q, so
it stays invertible. Row term counts are preserved: each term is a symbol
and symbols only swap sides. Substitute the new parameterization
x = Q + (old expression minus x-part) in Phi only; Phi stays injective
because the swap is an invertible affine reparameterization. No gate or
row is added, s is unchanged.

**(M4) Delete.** The following deletions each remove at least one row or
gate (s drops) and create no new node:

* (D1) A gate whose table has one effective input (constant pin, repeated
  pin, or unary table) is rerouted: its readers are repointed to the
  surviving value, adjusting their tables by the absorbed unary function
  (a constant, a literal, or a negated literal; in the constant case the
  readers' tables are adjusted by constants). The gate is removed. If the
  surviving parent is a parameter, it was already read by this gate, so
  its weight class cannot worsen (Section 3).
* (D2) A gate with no readers that is not the output is removed.
* (D3) A row containing at most two distinct nonconstant terms must
  mention an auxiliary value (the parameter-only rows would make A
  singular: a row containing only parameters contradicts invertibility of
  A over the auxiliary variables). Writing it v = U + c with U a single
  other term or absent, substitute this value into every reader of v and
  delete the row and v. The remaining auxiliary matrix stays invertible:
  row-reducing by the deleted row and expanding the determinant on v's
  column. Concretely, adding the (at most three-term) row to every other
  row mentioning v eliminates v's column outside the deleted row, and
  cofactor expansion on that column reduces the determinant to that of
  the smaller invertible system. Term counts of remaining rows do not
  increase.
* (D4) A parameter x of degree exactly one (Section 3), whose unique
  reader is a row E. Since A^(-1)N has a nonzero x-column, some auxiliary
  value Q depends on x. After the D5 pass below, Q has a reader (or is
  the output, an affine literal, impossible while M > d), so a pivot
  exchanging x and Q by (M3) is available and H-neutral. After the pivot,
  the old node x appears only in E; it is not the output (the output
  would be affine, excluded while M > d). Delete the node x and the row E;
  invertibility again follows by cofactor expansion on x's column, whose
  entries outside E's row are the (now zero) references to x. The row
  disappears from s.
* (D5) An auxiliary value referenced by no row, no gate, and not the
  output, together with its row, is deleted. Its column in A is zero
  outside its own row, so cofactor expansion on that column proves the
  remaining matrix invertible. After a D5 pass to exhaustion, every
  auxiliary value is referenced (or is the output).

## 3. Potential

For a parameter z, its **degree** is the number of distinct rows and
gates (not pins, not the output reference, not Phi) that read z. Call z
**used** if its degree is at least one. Define

    H = 2s + 6 * (number of used parameters).

Moves M1, M3 preserve H or change it as follows: M1 unchanged; M3 swaps a
used parameter for a used parameter (the new parameter is read by Q's old
readers, which include at least the gate or row that forced the pivot),
so H is unchanged. D1, D2, D3, D4 each lower s by at least one and never
increase the number of used parameters (deletions only reroute to
constants or to values that were already read), so Delta H <= -2 for each.
M2 fixes one parameter: if it was used, its weight 6 disappears, and H
changes by the effect of the induced reroutings, which never add a used
parameter (M2 substitutes constants; fixed-pin gates become D1 cases).
In particular H never increases under any move, and initially

    H_0 = 2|C| + 6*(used input variables) <= 2|C| + 6m,

while H >= 0 always. The refinements specific to the two eliminating
cases appear in Section 4.

## 4. The two eliminating reductions

Work while M > d. First run M1 on all eligible gates in topological order,
then D1-D5 to exhaustion; each pass strictly lowers s or the gate count
and the passes terminate (each removes an object). Afterwards, if no
acyclic gate remains, the output is a constant, parameter literal, or
auxiliary literal: an affine function of u, impossible by Section 1. So
an AND-type gate exists. Take a topologically first AND-type gate G (one
no other AND-type gate precedes among the acyclic gates; its parents are
then values with no AND-type gate among their ancestors, hence
parameters or auxiliary values).

**Restart rule.** Any deletion D1-D5 that becomes applicable at any point
below is executed immediately. It lowers s by at least one with
Delta M = 0, hence Delta H <= -2, and Section 4 then restarts from the
top; this terminates because s strictly decreases.

**Preparation.** If a parent of G is an auxiliary value Q, express Q in
terms of the parameters through v = A^(-1)(N u + c). If that expression
is a constant, repoint the pin of G to that constant: G becomes
degenerate and D1 deletes it (restart rule). Otherwise pick a parameter
x in Q's support and pivot Q with x (M3). After the pivot, that parent is
the parameter Q. Re-run M1/D1-D5 checks. If the second parent is
auxiliary, express it likewise: if it is a constant, or equals x + c with
x the parameter on the other pin, repoint the pin accordingly (G becomes
degenerate, D1 deletes it, restart rule); otherwise its support contains
a parameter other than x and we pivot again. At most two pivots.
Now G reads two distinct parameters x, y. While some gate J other than G
reads both G and x, note that J computes a binary function of (x, y)
(since G = G(x,y)); give J that table with pins x and y. This rerouting
strictly reduces the number of gates reading G, so it terminates; it does
not change s or degrees of parameters other than possibly adding readers
of y (already used). Repeat with x and y exchanged. After this, no gate
other than G reads G together with x or y. Re-run D1-D5.

After preparation, either x or y has degree one — Case A — or both have
degree at least two — Case B.

**Case A (an input of G is an isolated parameter).** Say deg(x) = 1, so
G is the only reader of x. Fix the other input y to a forcing constant
(M2): G's table becomes constant, so G is removed by D1, y's weight 6
disappears with y (M drops by one), and x loses its only reader, so its
weight 6 disappears as well. While M > d, G is not the output (Section 1,
terminal values). A reader of y that is the designated output is fine:
D1's rerouting repoints the output reference to the surviving literal or
constant, which Section 1 permits. Thus

    Delta M = -1,   Delta H <= -2 - 6 - 6 = -14 <= -12.

**Case B (both inputs have degree at least two).** Write u = deg(x),
v = deg(y), both >= 2, and h = the number of readers of G. If h = 0, G is
not the output (Section 1) and D2 deletes it; so h >= 1. Fix x to a
forcing constant for G's table (M2). Then:

* every reader of G acquires a constant pin and is deleted by D1: that is
  h distinct gates;
* the readers of x — there are u of them, and G itself is one of them by
  the definition of degree; none of the other u-1 is a reader of G by
  the preparation, so the two groups are disjoint — each becomes a row
  with at most two nonconstant terms (deleted by D3) or a gate with a
  constant pin (deleted by D1); G's own deletion is counted here, in the
  u group, not separately.

In total at least u + h >= 2 + 1 = 3 objects leave s, the parameter x
disappears with weight 6, and no parameter gains a reader (deletions
only reroute to constants or to nodes already read, Section 3). If a
reader of x is the designated output, D1's rerouting repoints the output
reference to the surviving literal or constant, which Section 1 permits.
Hence

    Delta M = -1,   Delta H <= -2(u+h) - 6 <= -12.

In both cases Delta H <= -12 * Delta M with Delta M = -1. Objects that
become deletable as consequences (cascade) only lower H further; they are
deletions, never additions.

## 5. Termination and conclusion

Each eliminating case strictly decreases M + s, and the preparation steps
of a single iteration are finite (each absorption, pivot, or rerouting
either is used at most twice, or strictly reduces the acyclic gate set).
Hence the process runs at most m + |C| iterations and must stop with
M_T <= d, all the while maintaining the invariant. Summing the potential
inequality over the run: every unit of lost M paid at least 12, and the
non-eliminating moves only lowered H further, so

    H_0 >= H_T + 12*(M_0 - M_T) >= 12*(m - d)      (H_T >= 0, M_0 = m).

With H_0 <= 2|C| + 6m,

    2|C| + 6m >= 12(m - d),   i.e.   |C| >= 3m - 6d.

If m <= d the right-hand side is nonpositive and there is nothing to
prove. This completes Lemma E.

No pair structure, weight split, or quadratic pivot was used; the rate
12 = 2*(u+h) + 6 at the minimum configuration u = 2, h = 1 is what the
flat potential permits, and Section 7 of the parameter-method literature
explains why heavier machinery is needed to push past 3.

## 6. Coverage of the B2 model

Every binary table is XOR-type (affine: 0, 1, u, u+1, v, v+1, u+v,
u+v+1), AND-type ((u+b)(v+c)+e), or degenerate (one effective input);
degenerate tables are D1. So the case split of Sections 2-4 is exhaustive
over full B2. Charging is by retained objects only: rows enter one-for-one
from absorbed gates (M1) and never from elsewhere; pivots (M3) add no row
or gate; the composed-table rerouting of the preparation replaces one
table by another on the same gate. Therefore the potential argument bounds
the original |C|, under arbitrary sharing, fan-out and depth.
