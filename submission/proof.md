# C_B2(MOD3_n) >= 2n - 4: a self-contained gate-elimination proof

This file proves the declared theorem completely. It uses no external premise.
The argument is the classical gate-elimination method of Schnorr (1974), written
out in full for the pinned Li–Yang B2 model. It claims no new result and is
weaker than the published 31/10 seed; it is submitted as a fully verified pair (2,1).

## 1. Family and NP eligibility

For an integer m >= 0, residue s in Z_3 and x in {0,1}^m, write |x| for the
Hamming weight and define

    M_{m,s}(x) = 1  iff  |x| = s (mod 3).

The scored family is f_n = M_{n,0} for every integer n >= 0 (so f_0 = 1, since
the empty sum is 0). It is total and defined at every length, with no exceptions.

**Uniform P and NP.** One fixed algorithm A reads x once, keeping a counter
modulo 3, and accepts iff the counter is 0 at the end. It runs in O(n) time at
every length. Let V(x,w) = 1 iff w is the empty string and A(x) accepts, and let
W(n) = 0. If x is in L then the empty w satisfies |w| <= W(|x|) and V(x,w)=1.
Conversely, if V(x,w)=1 for some w, then A(x) accepts, so x is in L. Hence L is
in P and so in NP. Witnesses are never circuit inputs; circuits for f_n have
exactly the n input bits of x.

## 2. Circuit model and conventions

A circuit is a finite directed acyclic graph. Its nodes are the input nodes
x_1..x_m, constant nodes 0 and 1 (free, any number), and gates. Each gate G has
an ordered pair of pins (u,v), each pin being a node earlier in a fixed
topological order (u = v allowed), and a truth table phi_G: {0,1}^2 -> {0,1},
any of the 16. One node is the output. Size |C| is the number of gates; inputs,
constants, wires and output selection cost nothing. The circuit computes its
output node's function, exactly, on every input. C(F) is the minimum size over
circuits computing F. This is the Li–Yang B2 model of rules.md: one output, two
pins per gate, all 16 tables, unbounded fan-out and depth, every retained gate
counted (including degenerate ones), constants free.

The **fan-out** of a node is the number of (gate, pin) pairs that point to it.

A gate G with pins (u,v) and table phi is **degenerate** if u = v, or u or v is
a constant node, or phi does not depend essentially on both of its arguments.
(phi depends essentially on its first argument if phi(0,b) != phi(1,b) for some b;
similarly for the second.) A circuit is **normalized** if it has no degenerate gate.

**Observation 2.1.** If G is degenerate then, as a function of the circuit inputs,
the value of G equals h(w) for a single node w and some unary h in
{0, 1, id, not}. Indeed: if u = v, G = phi(u,u); if u is a constant c,
G = phi(c,v), a unary function of v (symmetric for v); if phi ignores one argument,
G is a unary function of the other pin. If w is a constant node, h(w) is constant,
so we may take h constant.

A Boolean function F is a **literal** if F = x_k or F = not x_k for some k.

## 3. Removing degenerate gates

**Lemma 3.1.** Let C compute F, where F is neither constant nor a literal.
If C has a degenerate gate G, then some circuit C' with |C'| = |C| - 1 computes F.
Moreover C' is obtained from C by deleting G, changing only (i) pins that pointed
to G, (ii) truth tables of gates by negating one argument, and (iii) possibly
the output designation; every pin of C that pointed to an input or constant node
still points to the same node in C'.

*Proof.* Take h, w from Observation 2.1. In every case, redirected pins go to w,
which precedes G topologically, so acyclicity is kept. Negating an argument of a
table gives another of the 16 tables.

- h constant c. G is not the output (F is not constant). Point every pin that
  pointed to G at the constant node c; each reading gate sees the same value as
  before. Delete G.
- h = id. If G is the output, w is neither an input (F would be a literal) nor a
  constant (F would be constant), so w is a gate; make w the output. Point every
  pin that pointed to G at w. Delete G.
- h = not, with w not a constant. If G is not the output: point each pin that
  pointed to G at w and negate that argument in the reading gate's table, so its
  value is unchanged. Delete G. If G is the output: w is not an input (F would be
  a literal), so w is a gate. Replace phi_w by its complement, and negate the
  w-argument in the table of every gate reading w; all those gates keep their
  values, and w now computes the old value of G. Point the pins that pointed to G
  at w (no negation needed), make w the output and delete G.

In each case every node other than G computes the same function as before,
the output computes F, and exactly one gate is removed. Only items (i)–(iii)
change. []

**Corollary 3.2.** If F is neither constant nor a literal, some normalized
circuit of size at most C(F) computes F. (Apply Lemma 3.1 repeatedly to a
minimum circuit; the size strictly decreases, so this terminates.)

## 4. Properties of M_{m,s}

**Lemma 4.1 (restrictions).** Fixing x_k = a in M_{m,s} (m >= 1) gives
M_{m-1,s-a} on the remaining m-1 variables (by symmetry of M, after renaming).

**Lemma 4.2.** For m >= 2, M_{m,s} depends essentially on every variable; hence
it is neither constant nor a literal.

*Proof.* Fix k. Set every other variable to 0 except possibly one other variable
x_l (it exists since m >= 2) set to e in {0,1}. With x_k = 0 or 1 the weight is
e or e+1. Choosing e = 0 gives residues {0,1}; choosing e = 1 gives {1,2}. For
any s, some choice makes exactly one of the two weights congruent to s, so
flipping x_k changes the value. A function depending on m >= 2 variables is not
constant and not a literal. []

**Lemma 4.3 (three subfunctions).** For m >= 3 and distinct k, l, the four
restrictions of M_{m,s} given by (x_k,x_l) = (a,b) include at least three
distinct functions of the other m-2 >= 1 variables.

*Proof.* By Lemma 4.1 the restriction is M_{m-2,s-a-b}; a+b takes values 0,1,2,
giving the three residues s, s-1, s-2. For q = m-2 >= 1 the functions
M_{q,0}, M_{q,1}, M_{q,2} are pairwise distinct: on the all-zero input and on
the input with a single one they take the value pairs (1,0), (0,1), (0,0). []

## 5. The elimination step

**Lemma 5.1.** Let m >= 3 and let C be a normalized circuit computing M_{m,s}.
Then for some variable x_k and every a in {0,1} there is a circuit of size at
most |C| - 2 computing M_{m-1,s-a}.

*Proof.* C has a gate, since M_{m,s} is not a literal or constant (Lemma 4.2).
Let G be the first gate in topological order. Its pins are not gates, and they
are not constants and not equal (C is normalized), so they are two distinct
inputs x_k and x_l, and phi_G depends on both.

*Case B: x_k and x_l each have fan-out 1* (each is read only by its pin of G).
Then by induction along the topological order, the value of every node other
than x_k and x_l is a function of g = G(x_k,x_l) and the remaining variables y
(no gate other than G reads x_k or x_l, and the output is a gate); in particular
M_{m,s}(x) = Phi(G(x_k,x_l), y) for some Phi. For each (a,b), the restriction is
Phi(G(a,b), .), which is one of Phi(0,.), Phi(1,.). So there are at most two
distinct restrictions, contradicting Lemma 4.3. Case B is impossible.

*Case A: otherwise.* Some pin x in {x_k, x_l}, say x_k after renaming, is read by
a gate H with H != G (G reads x_k through one pin only, as its pins are distinct).
Fix a in {0,1}. Replace the input node x_k by the constant node a, so that every
pin that pointed to x_k now points to the constant a. The result C_a has |C|
gates and, as a circuit on the other m-1 variables, computes M_{m-1,s-a}
(Lemma 4.1). Since m-1 >= 2, this function is neither constant nor a literal
(Lemma 4.2), so Lemma 3.1 applies to C_a.

In C_a, both G and H have a pin pointing to a constant node, so both are
degenerate. Apply Lemma 3.1 to G, giving C_a' with |C| - 1 gates. By the
"moreover" part, H is still present and its pin that points to the constant a
still points to it, so H remains degenerate in C_a'. Apply Lemma 3.1 to H. The
resulting circuit computes M_{m-1,s-a} with |C| - 2 gates. []

## 6. The theorem

**Theorem 6.1.** For every integer m >= 2 and every s in Z_3,
C(M_{m,s}) >= 2m - 4. In particular C_B2(f_n) >= 2n - 4 for every n >= 2.

*Proof.* Induction on m. For m = 2 the bound is 0, which is trivial.
Let m >= 3 and let C be a minimum circuit for M_{m,s}. M_{m,s} is neither
constant nor a literal (Lemma 4.2), so by Corollary 3.2 there is a normalized
circuit C* computing it, with |C*| <= |C|. By Lemma 5.1 there are a and a circuit
of size at most |C*| - 2 computing M_{m-1,s-a}. By induction,
|C*| - 2 >= 2(m-1) - 4, so |C| >= |C*| >= 2m - 4. []

## 7. The declared guarantee and projection

Declared: g(n) = 2*n (coefficient 2, exponent 1, log and loglog exponents 0),
finite projection (2,1). For every epsilon > 0 let N_epsilon = max(2, ceil(4/epsilon)).
For every integer n >= N_epsilon, Theorem 6.1 gives
C_B2(f_n) >= 2n - 4 >= (2 - epsilon) n, since epsilon*n >= 4. Equivalently
liminf over all integers n of C_B2(f_n)/(2n) >= 1. The remainder is the constant 4,
which is o(n). This holds at every sufficiently large integer length, with the one
fixed family, verifier and witness bound of Section 1. Coefficient and exponent
come from the same theorem.

## 8. Scope and limitations

- This is the classical Schnorr-type bound. It is weaker than the 31/10 seed and
  it is not a new circuit lower bound.
- The theorem is unconditional and covers every length n >= 2; n = 0, 1 need no
  bound for the asymptotic claim.
- No external result is used. The only facts used are proved above: the removal
  lemma (Section 3), the restriction properties of M (Section 4) and the
  elimination step (Section 5).
- No finite computation or checker supports or replaces any step.
