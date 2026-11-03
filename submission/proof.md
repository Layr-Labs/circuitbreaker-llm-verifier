# Coefficient 2 at exponent 1 for the total family MOD3: complete supplied proof

This file is the complete evidence for the declaration in `bound.json`. It defines one
total language L in P (hence in NP with empty witnesses), gives a single polynomial-time
evaluator for every length, and proves for one and the same family, at every sufficiently
large integer length and in the full B2 model fixed by `rules.md`,

    C_B2(f_n) >= 2n - 8   for every integer n >= 4,

so that liminf over all integers n of C_B2(f_n)/n >= 2, with the explicit threshold
N_eps = max(4, ceil(8/eps)) for every eps > 0. No endpoint, no exponent above one, no
class separation and no optimality of the constant is claimed.

---

## 1. Model and conventions

**1.1 Circuits.** The B2 model pinned by `rules.md`: deterministic acyclic circuits with
one output, gates of fan-in two with ordered pins, all sixteen binary Boolean operations
available, unrestricted depth and fan-out, repeated pins and shared subgraphs allowed, and
exact computation required on every input. Size is the number of retained gates; input
vertices, wires and output selection cost zero. Constant nodes introduced during
substitutions cost zero.

**1.2 Negations are free on wires and at the output.** A negation may be absorbed into the
truth table of the gate that reads it, and the output of a circuit may be a constant or a
*literal* of an internal value (`rules.md`, circuit model). Consequently, if a gate
computes a constant or a literal of one of its inputs, deleting it and rerouting all its
uses to that constant or literal is free of charge. We use only this standard
normalisation; no uncharged computational basis is introduced, since every retained gate
is still counted.

---

## 2. Part A: the fixed total family and eligibility

**2.1 Definition.** For every integer n >= 0 and every x in {0,1}^n, let

    f_n(x) = 1  if and only if  wt(x) is congruent to 0 modulo 3,

where wt(x) is the Hamming weight of x. Since n = |x| is the actual input bit length, this
is a total family on every length with no exceptions: f_0(empty) = 1 (the empty sum is 0),
f_1 is 1 exactly on the all-zero string, f_2 = 1 exactly on 00 and 11, and so on. Let
L = { x : f_|x|(x) = 1 }, so that f_n is the characteristic function of L on {0,1}^n.
The family is fixed once and for all: no per-length advice, no oracle, no nonconstructive
choice, and no dependence of the algorithm on n.

**2.2 One polynomial-time evaluator.**

*Lemma A1.* There is a single deterministic algorithm A that, on input x in {0,1}^N,
computes f_N(x) in O(N) time.

*Proof.* A scans x once from left to right, maintaining a counter c in {0,1,2} initialised
to 0; for each 1-bit read it replaces c by c + 1 mod 3. After the scan it outputs 1 if
c = 0 and 0 otherwise. The counter is a constant-size value, so the running time is linear
in N and the algorithm is deterministic. The output is 1 exactly when wt(x) = 0 mod 3,
which is the definition of f_N(x). QED

**2.3 NP membership, both directions.** Let A be the algorithm of Lemma A1 and define the
verifier V(x,w) = 1 iff w is the empty string and A(x) = 1, with polynomial witness bound
W(N) = 0. Then:
 (i) if x lies in L then the empty witness w = empty satisfies |w| <= W(|x|) = 0 and
     V(x, empty) = 1;
 (ii) if V(x,w) = 1 for some w with |w| <= W(|x|) = 0 then w is empty and A(x) = 1, so
     x lies in L.
Both directions hold, V runs in time O(|x|) (Lemma A1 plus a constant-size check), so L
lies in P and hence in NP. The witness is not an input to the target circuits: the
circuits of section 3 receive x alone. This discharges eligibility: one fixed
polynomial-time evaluator and one fixed polynomial witness bound work for every length,
and neither changes across the asymptotic quantifiers.

---

## 3. Part B: the lower bound

**3.1 The three-settings lemma.** For a setting s of two selected input coordinates of an
n-bit input, write f_n|_s for the function of the remaining m = n - 2 coordinates obtained
by fixing those two coordinates to s.

*Lemma B1.* Let n >= 4, let i != j be two coordinates, and let m = n - 2 >= 2. The three
restrictions of f_n to the settings (0,0), (0,1) and (1,1) of (x_i, x_j) are pairwise
distinct functions of the remaining m coordinates.

*Proof.* Write sigma(s) for the sum of the two fixed bits and y for the remaining
coordinates. By the definition of f_n,

    f_n|_s(y) = 1  iff  wt(y) + sigma(s) is congruent to 0 modulo 3,

i.e. f_n|_s is the indicator of the residue class -sigma(s) modulo 3 of wt(y). The three
settings (0,0), (0,1), (1,1) have sigma = 0, 1, 2, so the three restrictions are the
indicators of the residue classes 0, 2 and 1 respectively. Since m >= 2, the m-bit strings
realise the weights 0, 1 and 2 (take the all-zero string, a string with one 1, and a
string with two 1s; for m = 2 the last two coincide as weights 1 and 2 on distinct
strings), so each of the three residue classes is realised by some y and missed by
another. Hence the three indicators are pairwise distinct functions. QED

**3.2 Normal form of a minimum circuit.**

*Lemma B2.* Let n >= 1 and let C be a circuit of minimum size g computing f_n. Then no
gate of C is constant or unary: every gate reads two distinct internal values on its two
pins, and in particular no pin is a constant and no gate reads the same value twice. The
designated output is not a constant or unary gate either.

*Proof.* Suppose some gate G of C is constant or unary, i.e. computes a constant or a
literal of one (or both) of its pins. Then G's function is a constant or the value of a
single internal value t (possibly negated). Reroute every use of G to that constant or
literal: by section 1.2 this costs nothing, and the function computed by the circuit is
unchanged. If G is the designated output, the same operation makes the output a constant
or a literal, which the model allows at no cost, and if G is not the output, drop it. In
either case we obtain a circuit computing f_n with at most g - 1 gates (one gate fewer, no
gate added), contradicting the minimality of C. Hence no such gate exists. QED

*Remark.* Lemma B2 is the only place where the normalisation convention of section 1.2 is
used. It is exactly the normalisation required by `rules.md` (constants and unary effects
absorbed, with every retained gate still counted); in particular we never delete a gate
whose value is needed by two different pins of the same descendant without rerouting both.

**3.3 The elimination step.**

*Lemma B3.* Let n >= 4 and let C be a circuit of minimum size g computing f_n. Then
C_B2(f_(n-1)) <= g - 2.

*Proof.* Let G be the first gate of C in a topological order of its gates. Every pin of G
is an input vertex (a variable or a constant), because any earlier gate would precede G.
By Lemma B2, no pin of G is a constant and G's two pins carry distinct values, so they are
two distinct input variables; call them x_i and x_j.

Consider the three settings s_1 = (0,0), s_2 = (0,1), s_3 = (1,1) of (x_i, x_j). The
two-input Boolean function computed by G takes at most two values, so among the three
settings two, say s_p and s_q (p != q), give the same value of G. Suppose, for
contradiction, that neither x_i nor x_j is read by any gate of C other than G. Then the
only wires leaving {x_i, x_j} lead into G, so the value computed by C is obtained from the
value of G and the remaining n - 2 variables alone; equivalently, for each setting s of
(x_i, x_j) the circuit computes f_n|_s as the function of the remaining variables obtained
by substituting the constant G(s) for the output of G. Since G(s_p) = G(s_q), the two
resulting functions of the remaining variables coincide, i.e. f_n|_(s_p) = f_n|_(s_q),
contradicting Lemma B1. Hence at least one of x_i, x_j is read by a gate of C other than
G; that variable, call it x, is read by at least two gates of C (G and that other gate).

Fix x = 0 and consider the circuit C' obtained as follows: delete every gate of C that
reads x, and reroute each of their uses to the constant or literal that the deleted gate
computes (after the substitution x = 0 such a gate reads one remaining value or no value
at all, so its output is a constant or a literal of its other pin; rerouting is free by
section 1.2). At least two gates are deleted, namely G and the other gate reading x, and no
gate is added, so |C'| <= g - 2. Moreover C' reads only the n - 1 remaining input
variables, and for every assignment y of those variables it outputs the value that C
outputs on the corresponding input with x = 0 (deleting gates and rerouting their uses
preserves the computed function, since each deleted gate's output was replaced by an equal
value everywhere). Therefore, for every y,

    C'(y) = C(0, y) = f_n(0, y) = [wt(y) + 0 is congruent to 0 modulo 3] = f_(n-1)(y),

because fixing a coordinate to 0 does not change the weight. Hence C' is a circuit of size
at most g - 2 computing f_(n-1), which is the same family at the length n - 1. Therefore
C_B2(f_(n-1)) <= g - 2 = C_B2(f_n) - 2. QED

**3.4 The lower bound.**

*Theorem B4.* For every integer n >= 4, C_B2(f_n) >= 2n - 8.

*Proof.* Let g_n = C_B2(f_n), which is finite and attained (there is a circuit of size
O(n): compute the weight modulo 3 by a chain of constant-size gadgets, or evaluate the
linear-time algorithm A directly, which yields a circuit of size O(n)). By Lemma B3,
applied to a minimum circuit of size g_n for n >= 4,

    g_(n-1) <= g_n - 2, i.e. g_n >= g_(n-1) + 2.

Iterating from n down to 4 and using the trivial bound g_4 >= 0 gives
g_n >= g_4 + 2(n - 4) >= 2n - 8. QED

---

## 4. Part C: quantifiers and the explicit threshold

**4.1 The liminf.** By Theorem B4, for every integer n >= 4,

    C_B2(f_n) / n  >=  (2n - 8) / n  =  2 - 8/n.

Hence liminf over all integers n of C_B2(f_n)/n >= 2, with the same single family at every
length: the family, its evaluator and its verifier do not depend on n, and the bound is an
all-length statement, not a statement about a sparse set of lengths or about limsup.

**4.2 Explicit threshold.** Let eps > 0 and let N_eps = max(4, ceil(8/eps)). For every
integer n >= N_eps we have 8/n <= eps, so C_B2(f_n) >= (2 - eps) n. This is the
requirement "for every eps > 0 there is N_eps such that for every integer n >= N_eps,
C_B2(f_n) >= (2 - eps) n^1", i.e. the declared guarantee liminf C_B2(f_n)/n >= 2 with the
declared finite projection (coefficient 2, exponent 1) for that one and the same theorem.

**4.3 Why the remainder is not a proof obligation beyond this.** The declared bound
c*n^p - R(n) with R(n) = 8 and R(n)/n -> 0 is of the form admitted by `rules.md`
("a proved bound c*n^p - r(n) suffices when r(n)/n^p tends to zero"). No suitable-family
change occurs between the coefficient and the exponent: both coordinates describe
Theorem B4 for the family of section 2.1.

---

## 5. Part D: scope and limitations

 * **The coefficient is 2, and no more is claimed.** Nothing here proves 2n - o(1) with a
   smaller remainder, nothing proves coefficient above 2, nothing proves any exponent
   above 1, and nothing separates complexity classes. In particular this file makes no use
   of affine dispersers, of the Ben-Sasson-Kopparty construction, or of any external
   theorem about them; it is a self-contained elementary gate-elimination argument.
 * **The remainder 8 is an artefact of the induction base.** The proof gives
   C_B2(f_n) >= C_B2(f_4) + 2(n-4) and only C_B2(f_4) >= 0 is used, so the constant 8 could
   be improved by a better base case; no such improvement is claimed, and the remainder
   plays no role in the coefficient.
 * **The normalisation is the rules' normalisation.** Deleting a unary or constant gate and
   rerouting to a constant or literal is free precisely because negations may be absorbed
   into descendant truth tables and the output may be a literal (`rules.md`); every
   retained gate is counted, so no uncharged computational basis is introduced.
 * **Provenance.** The gate-elimination bound for MOD3 is classical (Schnorr, 1974); the
   argument as used here is stated in the introduction of Mishra and Santhanam, ECCC
   TR26-260 revision 1 (2026-11-03), whose abstract and introduction I read while working
   on this benchmark. No prior submission's archive, proof text or code was consulted for
   this entry: the package here was written from the definition of MOD3 and the
   conventions of `rules.md`, and every step used is proved in Parts A-C above. The family
   is the canonical MOD3 and therefore necessarily coincides with the MOD3 families of
   other entries; the proof is this file's own.
 * **No computational claim.** No experiment, checker or finite computation is needed or
   declared: the argument is symbolic at every length.
