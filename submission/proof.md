# MOD3: C_B2(f_n) >= 2n-3 at every length n >= 2, hence the pair (2, 1)

This proof is fully self-contained. It cites no external theorem: it is Schnorr's
gate-elimination argument (1974), written out with every step proved here,
including the handling of constants, repeated pins and degenerate gates.

## 1. The family and its uniform NP membership

For r in {0,1,2} and m >= 0 let M_r^m : {0,1}^m -> {0,1},
M_r^m(y) = 1 iff (y_1 + ... + y_m) = r (mod 3). Residues r are always taken mod 3.

The language is L = {x in {0,1}* : the number of 1s in x is divisible by 3}, so
f_n = M_0^n for every n >= 0. In particular f_0(empty) = 1.

Evaluator A: scan x, keep a counter mod 3, and output 1 iff it ends at 0. This is
one fixed deterministic O(n)-time algorithm for all lengths.

NP verifier: W(n) = 0 and V(x, w) = 1 iff w is empty and A(x) = 1.
- If x is in L, then w = empty (|w| = 0 <= W) gives V(x, w) = 1.
- If V(x, w) = 1 for some w with |w| <= 0, then A(x) = 1, so x is in L.

So L is in P and hence in NP. There is no advice and no per-length choice.
Witnesses are not circuit inputs.

Restriction rule (R). Fixing one variable of M_r^m to c in {0,1} gives M_(r-c)^(m-1)
on the remaining variables.

## 2. Circuit model and conventions

A circuit is a finite DAG with:
- input nodes (one per variable);
- constant nodes, carrying the fixed values 0 or 1;
- gates, each with an ordered pair of pins (u, v) wired from earlier nodes, where
  u = v is allowed, and a truth table T in the set of all 16 maps {0,1}^2 -> {0,1}.

Gate G computes T(val(u), val(v)). Fan-out and depth are unrestricted. One node is
selected as the output. The size |C| is the number of gates; inputs, constants,
wires and output selection cost zero. C_B2(h) is the minimum size of a circuit
computing h exactly on every input. It exists because every Boolean function has
some circuit. This is the full B2 model of rules.md. The proof below uses no other
property of the model.

## 3. Facts about the family

(F1) For m >= 2 and every r, M_r^m depends on every variable. Hence it is neither
constant nor a literal (y_j or NOT y_j), because it depends on at least two
variables.
Proof. Fix a variable y_j and set it to 0. The other m-1 >= 1 variables can have
sum t = 0 or t = 1. Set y_j = 0 and choose t:
- r = 0: t = 0 gives value 1, and flipping y_j (sum 1) gives 0;
- r = 1: t = 0 gives value 0, and flipping gives sum 1, value 1;
- r = 2: t = 1 gives value 0, and flipping gives sum 2, value 1.

(F2) For m >= 1, the three functions M_0^m, M_1^m, M_2^m are pairwise distinct.
Proof. At y = 0...0 the values are 1, 0, 0. At y = 10...0 they are 0, 1, 0.

(F3) For n >= 3, any two distinct variables x_i, x_j, and every r: the four
restrictions of M_r^n at (x_i, x_j) = (a, b) give at least three distinct
functions of the other n-2 >= 1 variables.
Proof. By (R) the restriction is M_(r-a-b)^(n-2). As (a, b) ranges over {0,1}^2,
a+b takes the values 0, 1, 2, so all three residues r, r-1, r-2 occur. By (F2)
these three functions are distinct.

## 4. Removing a degenerate gate

Call a gate G with pins (u, v) and table T degenerate if
(i) u = v, or (ii) u or v is a constant node, or (iii) T depends on at most one of
its two arguments.
In each case there is a node w in {u, v} and a map phi: {0,1} -> {0,1} with
val(G) = phi(val(w)) on every input:
- (i): w = u and phi(a) = T(a, a);
- (ii): if u is the constant c, then w = v and phi(b) = T(c, b), and symmetrically
  for v;
- (iii): w is the argument T depends on, or either one if T is constant.
phi is one of 0, 1, id, NOT.

Lemma D. Let C compute h, where h is neither constant nor a literal, and let G be a
degenerate gate of C. Then there is a circuit C' computing h with |C'| = |C| - 1.
Every gate of C' except G is a gate of C with the same pin sources, except that pins
that came from G are rewired, and possibly with a changed truth table. In
particular, every pin from a constant node in C that did not come from G is still
present in C'.

Proof. Take w and phi as above; w precedes G.

Case A: G is not the output. Delete G. For each pin of a gate H that was wired
from G:
- phi = c constant: wire that pin from a constant node c (adding one if needed;
  constants are free);
- phi = id: wire it from w;
- phi = NOT: wire it from w and replace H's table by the table with that argument
  negated, e.g. T'(a, b) = T(1-a, b) if it is the first pin.
w precedes G, which precedes H, so the graph stays acyclic. By induction in
topological order every remaining node keeps its value on every input. The output
is unchanged, and one gate is removed.

Case B: G is the output. Then h = phi(val(w)).
- If phi is constant, or w is a constant node, h is constant: excluded.
- If w is an input node, h is a literal: excluded.
So w is a gate W.
- If phi = id, rewire G's successors, if any, from W as in Case A, select W as the
  output, and delete G.
- If phi = NOT, replace W's table by its negation and, for every pin wired from W
  in a gate other than G, negate that argument in that gate's table, so all those gates
  keep their values. Then select W, which now computes NOT(old W) = h, as the
  output, delete G, and rewire G's own successors as in Case A, using phi = id
  from the new W.
In both cases the gate count drops by one.

Constant-node pins of gates other than G are never touched, which proves the last
sentence. []

Corollary M. In a minimum-size circuit for h, where h is neither constant nor a
literal, no gate is degenerate (otherwise Lemma D gives a smaller circuit for h), and
the output is a gate (otherwise h would be constant or a literal).

## 5. The induction

Let s(m) = min over r of C_B2(M_r^m).

Base, m = 2. By (F1), M_r^2 is neither constant nor a literal, so the output of any
circuit for it is a gate: s(2) >= 1 = 2*2 - 3.

Step, n >= 3. Fix r and a minimum circuit C for f = M_r^n. f is neither constant
nor a literal by (F1), so by Corollary M C has a gate, and none of its gates is
degenerate. Let G be the first gate in a topological order. Its pin sources are not
gates, so each is an input or constant node. Since G is not degenerate, neither
source is constant and they are distinct, so they are two distinct inputs x_i, x_j.

Claim: x_i or x_j has a pin into a gate other than G.
Suppose not. Then no node other than G reads x_i or x_j. By induction in topological
order, the value of every node other than x_i, x_j and G is a function of val(G) and
the other n-2 inputs. Hence f = F(T(x_i, x_j), rest) for some function F. So the
four restrictions at (x_i, x_j) = (a, b) are among F(0, .) and F(1, .), which is at
most two distinct functions. This contradicts (F3).

Say x_i has a pin into a gate H != G. H != G because G has exactly one pin from x_i,
its sources being the distinct x_i and x_j.

Substitute x_i := 0: replace the input node x_i by a constant node 0. The result C0
is a circuit with |C| gates on the remaining n-1 inputs computing f|_(x_i=0) =
M_r^(n-1) by (R). This function is neither constant nor a literal by (F1), since
n-1 >= 2. In C0, both G and H have a pin from a constant node, so both are
degenerate (ii).
- Apply Lemma D to G: this gives C1 with |C| - 1 gates computing M_r^(n-1).
- H is still a gate of C1 and still has its pin from the constant node, because
  that pin did not come from G. So H is still degenerate.
- Apply Lemma D again: this gives C2 with |C| - 2 gates computing M_r^(n-1).

Therefore C_B2(M_r^n) = |C| >= C_B2(M_r^(n-1)) + 2 >= s(n-1) + 2. Minimizing over r
gives s(n) >= s(n-1) + 2, and with the base case s(n) >= 2n - 3 for all n >= 2.
In particular

    C_B2(f_n) = C_B2(M_0^n) >= 2n - 3   for every integer n >= 2.        (*)

## 6. The declared asymptotic guarantee

The declaration is the polynomial_log guarantee g(n) = 2 * n^1 * (log2 n)^0 *
(log2 log2 n)^0 = 2n, with finite projection (c, p) = (2, 1).

By (*), C_B2(f_n)/g(n) >= 1 - 3/(2n) for every integer n >= 2. Given eps > 0, let
N_eps = max(2, ceil(3/(2 eps))). Then C_B2(f_n) >= (1 - eps) * 2n for EVERY integer
n >= N_eps, not just a subsequence. So liminf over all integers n of
C_B2(f_n)/g(n) >= 1, which is the declared full guarantee with c = 2 > 0 and p = 1.
Its finite projection, liminf C_B2(f_n)/n >= 2, follows for the same fixed family,
verifier and witness bound, so both coordinates come from this one theorem.

## 7. Scope

This is a classical, fully proved result. It is below the literature value 31/10
and claims no new technique.
