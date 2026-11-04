# Coefficient 2 for MOD-3: a self-contained gate-elimination proof

This submission proves, with no external premise, that the MOD-3 family below
satisfies C_B2(f_n) >= 2n-3 for every integer n >= 2. Hence the liminf of
C_B2(f_n)/n is at least 2, and the declared pair is (2,1). The method is the
classical Schnorr (1974) gate-elimination argument. It is reproduced here in
full for the pinned model, so the proof does not depend on the citation. This
is NOT an improvement on the published 31/10 seed. It is a weaker theorem whose
every step is supplied and checkable from this file alone.

Notation: Z_3 = {0,1,2} with arithmetic mod 3. For a bit string a, wt(a) is the
number of ones in a.

## 1. The fixed total family

For every integer n >= 0 and every x in {0,1}^n, define

    f_n(x) = 1  if wt(x) = 0 (mod 3),   f_n(x) = 0 otherwise.

This covers every length, including n = 0: the empty string has weight 0, so
f_0 = 1. The language is L = { x in {0,1}^* : wt(x) = 0 (mod 3) }.

More generally, for a finite set Y of m >= 0 Boolean variables and r in Z_3, let
M^Y_r (written M^m_r when only m = |Y| matters) be the function on assignments
a in {0,1}^Y with M^Y_r(a) = 1 iff wt(a) = r (mod 3). Then f_n = M^n_0 on the
variables x_1,...,x_n. Renaming input nodes maps circuits to circuits of the
same size, so C*(M^Y_r), defined in Section 3, depends only on |Y|.

## 2. Uniform P and NP membership

Algorithm A reads x once and keeps a counter s in {0,1,2}, starting at s = 0.
For each input bit b it sets s := (s + b) mod 3. It accepts iff s = 0 at the
end. By induction on the prefix length, s = wt(prefix) mod 3 after each step,
so A(x) = f_{|x|}(x). A runs in O(n) time at every length n = |x| with no
advice. One algorithm serves every length.

The NP verifier is V(x,w) = 1 iff w is the empty string and A(x) = 1. The
witness bound is W(n) = 0. If x is in L, then w = empty gives V(x,w) = 1 with
|w| = 0 <= W(|x|). Conversely, if V(x,w) = 1 for some w with |w| <= W(|x|), then
A(x) = 1, so x is in L. Hence L is in P, and so in NP. Witness bits are never
circuit inputs. The circuits below have exactly the n input bits of x.

## 3. Circuit model

Fix a finite variable set Y. A circuit C over Y is a finite sequence of nodes:

- one input node for each variable of Y;
- any number of constant nodes, each labelled 0 or 1;
- gates G_1, ..., G_s in this order. Each gate has a truth table
  g: {0,1}^2 -> {0,1} (any of the 16 binary tables) and an ordered pair of pins.
  Each pin points to an input node, a constant node or an EARLIER gate. The two
  pins of one gate may point to the same node.

One node, of any kind, is the output. Fan-out and depth are unrestricted, and
subcircuits may be shared. On an assignment a in {0,1}^Y, input nodes take the
values a(y), constant nodes take their labels, and a gate takes
g(value of pin 1, value of pin 2). Gates are evaluated in order. We write
val_v(a) for the value of node v. C computes F if val_out(a) = F(a) for every a.
The size |C| = s counts EVERY gate, including degenerate, constant-valued and
negated-projection gates. Inputs, constants, wires and the choice of output cost
zero. Let C*(F) be the least size of a circuit computing F. This minimum exists
because every Boolean function has some circuit, for example a DNF built from
AND, OR and NOT tables (NOT is a gate with both pins equal).

Relation to the pinned model (Li–Yang, ECCC TR21-023, Sections 1.1, 2.1, 2.4
footnotes 3–4, 3.3). A Li–Yang B2 circuit has binary gates with all 16 tables,
unrestricted fan-out and depth, one output, and computes the function exactly on
every input. Its size counts all retained gates. Order its gates topologically,
and turn any constant inputs into constant nodes. The result is a circuit in the
sense above, of the same size. The model above may be more permissive:
constant nodes are free and the output may be any node. Any restriction the
pinned model adds can only increase the minimum, so

    C_B2(F) >= C*(F)  for every F.                                   (3.1)

All lower bounds below are proved for C*. They therefore hold in the pinned
model. Constant nodes created by substitution are free, as the rules require.
Lemma 5.2 does the job of Li–Yang normalization: it absorbs constant pins,
unary or repeated-pin gates and negations into truth tables. Negations go into
descendant tables, or into the predecessor table at the output. Inequality (3.1)
carries the resulting bound into the pinned model. No extra computational basis is introduced:
every circuit built in the proof has the binary gates defined above.

## 4. Facts about M^m_r

Lemma 4.1 (restriction). Let y be in Y and c in {0,1}. Substituting y = c into
M^Y_r gives M^{Y \ {y}}_{r-c}. Substituting y = a and y' = b, for distinct y and y',
gives M^{Y \ {y,y'}}_{r-a-b}.

Proof. The weight of the full assignment equals the weight of the remaining
variables plus c, or plus a+b in the second case. Compare residues mod 3. QED.

Lemma 4.2 (distinctness). For k >= 1, the functions M^k_0, M^k_1, M^k_2 are
pairwise distinct.

Proof. Take distinct residues s and s'. The weights 0 and 1 are attainable when
k >= 1, and they cover two different residues. So at least one of s and s', say
s, is the residue of an attainable weight w in {0,1}. On an assignment of weight
w, M^k_s = 1 and M^k_{s'} = 0. QED.

Lemma 4.3 (essential dependence). For m >= 2 and any r, M^m_r depends
essentially on every one of its m variables. In particular it is not constant,
and it is not equal to the value of a single variable or its negation.

Proof. Weights 0, 1 and 2 are attainable because m >= 2. So some assignment a has
wt(a) = r (mod 3), and M^m_r(a) = 1. Flipping any coordinate changes the weight
by +1 or -1, so the residue is no longer r and the value becomes 0. So every
variable is essential. A function with two or more essential variables is not
constant and not a literal. QED.

Lemma 4.4 (three subfunctions). Let m >= 3 and let y, y' be distinct variables.
The four restrictions M^m_r|_{y=a, y'=b}, for (a,b) in {0,1}^2, include at least
three pairwise distinct functions.

Proof. By Lemma 4.1 the restriction is M^{m-2}_{r-a-b}. The sums a+b = 0, 1, 2
give the three residues r, r-1 and r-2. These functions are pairwise distinct by
Lemma 4.2, because m-2 >= 1. QED.

## 5. Eliminating a reducible gate

Definition 5.1. In a circuit C, a gate G is reducible if there are a node w that
precedes G (w is not G) and a map h: {0,1} -> {0,1} such that
val_G(a) = h(val_w(a)) for every assignment a. The four possible h are
constant 0, constant 1, identity and negation.

Two sufficient conditions. Let G have pins (u,v) and table g.
(i) If u is a constant node labelled c, then val_G = g(c, val_v). Take w = v and
h(z) = g(c,z). The case where v is constant is symmetric.
(ii) If u = v, then val_G = g(val_u, val_u). Take w = u and h(z) = g(z,z).
In both cases w precedes G.

Lemma 5.2 (elimination). Let C compute F, where F depends essentially on at
least two variables. Let G be a reducible gate of C. Then some circuit C' over
the same variables computes F and satisfies:

- (a) the gates of C' are exactly the gates of C other than G, so
  |C'| = |C| - 1, though some truth tables may change;
- (b) every pin of a gate of C' that pointed in C to a node other than G still
  points to that same node. Only pins that pointed to G are redirected.

Input and constant nodes are kept. New constant nodes may be added; any gate may
reference a constant node, so where it sits in the sequence does not matter.

Proof. Fix w and h as in Definition 5.1. In every case, G is deleted after
redirecting the pins that point to it. Redirected pins point to w, to W = w, or
to a new constant node. All of these precede G, and G precedes every gate that
read it. So the new pins point to earlier nodes and C' is a valid circuit.

Case 1: h is constant, or w is a constant node. Then val_G is a constant c'.
G is not the output, since otherwise F would be constant, contradicting the
hypothesis that F has essential variables. Add a
constant node labelled c' and redirect to it every pin that pointed to G. Every
gate reads the same values as before, so every value, and the output, is
unchanged.

Case 2: h is the identity. Redirect to w every pin that pointed to G. If G was
the output, make w the output. Values are unchanged.

Case 3: h is negation and w is not a constant node.

- 3a: G is not the output. For each gate H that has one or both pins pointing to
  G, redirect those pins to w. Replace H's table by the table obtained by
  negating each argument whose pin was redirected. For example, if only pin 1
  was redirected, use g'_H(p,q) = g_H(1-p, q). H then sees val_w where it used to
  see val_G = 1 - val_w, so H computes the same value. All other values are
  unchanged.
- 3b: G is the output. Then F = 1 - val_w. If w is an input node, F is the
  negation of one variable, contradicting the assumption that F has two
  essential variables. So w is a gate W. Replace W's table g_W by 1 - g_W.
  W's own pins point to nodes before W, which are unaffected, so the new value
  of W is 1 - (old val_W) = val_G. For every gate H other than G and W, first redirect
  to W each pin that pointed to G, keeping that argument of H's table unchanged,
  since the new val_W equals the old val_G. Then, for each pin of H that pointed
  to W in C, keep it pointing to W and negate that argument in H's table, since
  the new val_W is the negation of the old val_W. Make W the output.

  By induction along the gate order, every gate other than W and G has the same
  value as in C. Its arguments are either unchanged nodes, or W with the
  matching adjustment. The output value is the new val_W = val_G = F.

In every case only G is deleted, which gives (a). Only pins that pointed to G
change target, which gives (b): in Case 3b, pins pointing to W still point to W.
QED.

Corollary 5.3. Let F depend essentially on at least two variables, and let C be
a circuit of size C*(F) computing F. Then C has no reducible gate. In
particular, by (i) and (ii), no gate of C has a pin to a constant node, and no
gate has both pins pointing to the same node.

Proof. Otherwise Lemma 5.2 yields a circuit of size C*(F) - 1 computing F,
contradicting minimality. QED.

## 6. Main lemma

Theorem 6.1. For every integer m >= 2 and every r in Z_3, C*(M^m_r) >= 2m - 3.

Proof, by induction on m.

Base m = 2. M^2_r depends essentially on both variables (Lemma 4.3). An input or
constant output node would compute a constant or a single variable. So the
output of any circuit computing M^2_r is a gate, and C*(M^2_r) >= 1 = 2*2 - 3.

Step m >= 3. Assume the claim for m-1 >= 2 and every residue. Let F = M^Y_r with
|Y| = m, and let C be a circuit of size C*(F) computing F. F has essential
variables (Lemma 4.3), so the output is a gate and C has at least one gate. Let
G be the first gate. Its pins point to earlier nodes, and no gate is earlier, so
both pins point to input or constant nodes. By Corollary 5.3, neither pin is
constant and the pins differ. So G reads two distinct input nodes y and y'.

Claim: some gate H other than G has a pin pointing to y or to y'.

Suppose not. The output is not y or y', since F is not a literal (Lemma 4.3).
Let C° be C with G replaced by a fresh input node z, and its pins dropped. In C,
the nodes y and y' influence no gate except through G. So C° computes a function
Phi(z, b) of z and the assignment b to Y \ {y,y'}, with

    F(a) = Phi( g(a(y), a(y')), b )  for every a,

where b is the restriction of a to Y \ {y,y'}. This includes the case where G is
the output, with Phi(z,b) = z. Therefore F|_{y=s, y'=t} = Phi(g(s,t), .) for each
(s,t). Since g(s,t) takes at most two values, at most two of the four
restrictions are distinct. This contradicts Lemma 4.4, because m >= 3. So the
claim holds.

Rename the two variables if necessary so that y is read by G and by a gate H
other than G. Let C_1 be C with the input node y relabelled as a constant node 0.
C_1 has size |C|, and it computes the restriction F|_{y=0} = M^{Y \ {y}}_r on
Y \ {y} (Lemma 4.1). By Lemma 4.3 this function has m-1 >= 2 essential
variables. In C_1, the gates G and H both have a pin to the constant node 0, so
both are reducible by (i).

Apply Lemma 5.2 to G in C_1. This gives C_2, computing M^{Y \ {y}}_r, of size
|C| - 1. By 5.2(a), H is still a gate of C_2. By 5.2(b), the pin of H that
pointed to the constant node 0 (a node other than G) still points to it. So H is
reducible in C_2 by (i). Apply Lemma 5.2 to H in C_2. This gives C_3, computing
M^{m-1}_r, of size |C| - 2. By the induction hypothesis,

    |C| - 2 >= C*(M^{m-1}_r) >= 2(m-1) - 3,

so C*(F) = |C| >= 2m - 3. QED.

Remarks on coverage:
- No upper bound on C is used, and no specific constant was needed: the
  substitution y = 1 works equally well.
- G and H are distinct gates of C. Lemma 5.2(a) keeps H when G is removed, so
  exactly two gates are charged.
- Degenerate gates, constant nodes, repeated pins, unary behaviour and a
  negated output are all covered by Lemma 5.2 and Corollary 5.3.
- Gates that are not ancestors of the output are still counted in |C|. The
  argument never needs to delete them.

## 7. The scored asymptotic statement

By (3.1) and Theorem 6.1 with Y = {x_1,...,x_n} and r = 0,

    C_B2(f_n) >= C*(f_n) >= 2n - 3   for EVERY integer n >= 2.        (7.1)

Here n is the actual input bit length. There is no padding, block encoding or
witness input. The remainder 3 is constant, so 3/n tends to 0. Explicitly, given
epsilon > 0, let N_epsilon = max(2, ceil(3/epsilon)). For every integer
n >= N_epsilon, (7.1) gives C_B2(f_n) >= 2n - 3 >= 2n - epsilon*n
= (2 - epsilon) n. So

    liminf_{n -> infinity} C_B2(f_n) / (2 n) >= 1.

This is the declared polynomial_log guarantee with c = 2, p = 1, q = 0 and r = 0,
that is g(n) = 2n. It holds at every sufficiently large integer length, not on a
subsequence. The declared finite projection (2,1) is the same statement:
liminf C_B2(f_n)/n^1 >= 2. The family, the verifier V and the witness bound
W = 0 are fixed before every quantifier.

## 8. Scope and limitations

- The coefficient 2 is below the published 31/10 seed (Li–Yang). This
  submission does not claim to improve the record. Its value is that the proof
  is complete and self-contained under the supplied-evidence policy.
- Nothing is claimed about the exact circuit complexity of MOD-3, and no upper
  bound is used or asserted.
- Provenance of the method: C. P. Schnorr, "Zwei lineare untere Schranken für
  die Komplexität Boolescher Funktionen", Computing 13 (1974) 155–171. It is
  also presented in standard texts, e.g. I. Wegener, The Complexity of Boolean
  Functions (1987), Chapter 5. Every step used here is re-proved above, so no
  citation is a premise.
