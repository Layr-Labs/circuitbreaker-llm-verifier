# A complete proof of the (31/10, 1) circuit lower bound for the BK/Shoup family

**Theorem.** Let (f_n) be the single total BK/Shoup family defined in §1. Then, in the full Li–Yang B2 model (deterministic acyclic circuits, one output, two input pins per gate, all 16 binary Boolean truth tables, unrestricted depth and fan-out, exact computation on every input, nonuniform circuits allowed),

    liminf_{n→∞} C_B2(f_n) / n ≥ 31/10,

i.e. C_B2(f_n) ≥ (31/10)n − o(n) for every sufficiently large integer length n.

This document and the companion files `cases-a.md`, `cases-b.md`, `disperser-a.md` and `disperser-b.md` contain the complete proof: Li & Yang's gate-elimination proof (ECCC TR21-023) for their Theorem 1.1, reproduced in full, Section 4.1 case analysis included. The one premise the first version of this submission cited externally — the affine-disperser property of the §1 construction (Ben-Sasson–Kopparty) — is now proved here: §2 reduces it to Theorem 2.2 of Ben-Sasson–Kopparty with every hypothesis verified for our parameters, and `disperser-a.md` and `disperser-b.md` reproduce that theorem's proof (their Sections 3, 6 and 7) in full. The only facts used without proof are classical background not about this construction: Bertrand's postulate (a prime lies between x and 2x, used once to bound k in §2) and the basic existence and trace facts of finite fields (the trace fact, Proposition 3.3, is proved in `disperser-a.md`). Everything else is proved here, including the potential bound of Lemma 5.1 (FGHK16), whose short proof §5 reproduces. Parameter values: α_I = 48/5, α_Q = 9/5, α_φ = 1/5, δ = 64/5.

## 1. The fixed total family

Set f_0(empty) = 0. For n ≥ 1, let k be the least prime strictly above 2·n^{3/5}, r = ⌈n/k⌉, t = ⌈√r⌉. Fix one deterministic polynomial-time procedure (Shoup's deterministic irreducible-polynomial construction over F_2, all choices fixed by that procedure) producing a monic irreducible q_k ∈ F_2[Z] of degree k; let K = F_2[Z]/(q_k) with the polynomial basis 1, Z, …, Z^{k−1}. This fixes a single language, not per-length advice.

Append rk − n zeros to x ∈ {0,1}^n, split into r blocks of k bits, and read block i as an element φ_i(x) ∈ K with its first bit the constant coefficient; the combined map φ is an injective F_2-linear map. Let π: K → F_2 take the constant coefficient; π is F_2-linear and π(1) = 1. For 1 ≤ i ≤ r put y_i = φ_i(x)^{2^{i+1} − 1} (Frobenius powers in K). Define

    f_n(x) = π( Sym_r^t(y_1, …, y_r) ),

where Sym_r^t is the elementary symmetric polynomial of degree t: the sum of the products indexed by the t-element subsets of {1, …, r}.

## 2. Uniformity, membership in P and NP, and the disperser property

One deterministic algorithm computes f_n at every length: scan for the least prime k above the threshold (trial division is polynomial, as k = O(n^{3/5}) by Bertrand's postulate), run the fixed Shoup procedure (polynomial in k), then evaluate the symmetric polynomial, avoiding characteristic-two division, by the descending recurrence E_0 = 1, E_j = E_j + y_i·E_{j−1} after reading y_i, for j from min(i, t) down to 1; output π(E_t). Splitting the subsets by membership of i proves the recurrence. Cost: O(r·t) field operations of polynomial bit cost; the powers y_i need O(i+1) multiplications each. Hence one polynomial-time evaluator A computes the whole family. With V(x, w) = 1 iff w is empty and A(x) = 1, and witness bound W ≡ 0, both directions of the NP-membership equivalence are immediate; the family is in P, hence in NP, with empty witnesses. The witness is not an input to the target circuits.

The lower bound uses one property of this family: its affine-disperser property. We prove it by verifying that the §1 construction is exactly the construction of Theorem 2.2 (equation (4)) of Ben-Sasson and Kopparty, *Affine Dispersers from Subspace Polynomials* (cited as [BK12] in the source paper), whose proof is reproduced in `disperser-a.md` (their Section 3) and `disperser-b.md` (their Sections 6 and 7).

*Parameters.* Apply the theorem with characteristic p = 2 and input length m := n (the actual length of f_n, not the padded length rk; see the padding paragraph below). Its block length, the least prime strictly above 2·m^{3/5}, is exactly our k; its number of blocks, r = ⌈m/k⌉, exactly our r; its symmetric-polynomial degree, t = ⌈√r⌉, exactly our t.

*The map.* The theorem quantifies over an arbitrary injective F_2-linear map φ: F_2^m → (F_{2^k})^r and an arbitrary nontrivial F_2-linear functional π: F_{2^k} → F_2. The map of §1 (append rk − n zeros, split into blocks of k bits, read each block in the polynomial basis of K) is F_2-linear and injective (a coordinate embedding followed by fixed per-block basis identifications), is defined because n ≤ rk, and under the fixed identification K = F_2[Z]/(q_k) ≅ F_{2^k} is one of the theorem's admissible maps φ. No divisibility of n by k is assumed or needed: the short last block is absorbed by φ, its last component having its top rk − n coordinates constantly 0. The π of §1 (the constant coefficient) is F_2-linear with π(1) = 1, hence nontrivial and admissible; by Proposition 3.3 (`disperser-a.md`) it has the form π(x) = Tr(µx) for a unique µ ≠ 0, the form used in the proof.

*The exponents.* In equation (4) the i-th block is raised to the power 1 + p + … + p^i. At p = 2 this is 2^{i+1} − 1 — exactly the exponent of y_i in §1, for i = 1, …, r. (Section 7 of the source indexes the exponents in reverse order; `disperser-b.md` keeps that order and identifies the two forms in the deduction of Theorem 2.2 by renaming the blocks.)

*Conclusion.* Theorem 2.2 therefore gives: **for every sufficiently large integer length n, f_n is nonconstant on every affine subspace of F_2^n of dimension strictly greater than 6·n^{4/5}.** Its proof is reproduced in full in the two disperser files; the only steps using our parameters are the parameter deduction (Theorem 7.3 ⇒ Theorem 2.2), whose arithmetic is written out at the end of `disperser-b.md`: Bertrand's postulate gives k ≤ 4n^{3/5}; hence r(r+1)/t ≤ (r+1)^{3/2} ≤ (n^{2/5}/2 + 2)^{3/2} < 2n^{3/5} < k, and the threshold nr/t + nt + r(r+1) is at most about 4.25·n^{4/5} < 6n^{4/5} when r ≈ n/k is estimated jointly with k (not by maximizing the two independently). Restricting to large n costs nothing: for n ≤ 6^5 = 7776, 6n^{4/5} ≥ n, so no affine subspace exceeds the threshold and the assertion is vacuous; and the argument below is a liminf statement using the disperser property only at large lengths. Set d(n) = ⌊6·n^{4/5}⌋ + 1; eventually 1 ≤ d(n) ≤ n and d(n)/n → 0.

**Definition (affine disperser).** f: F_2^n → F_2 is an *affine disperser for dimension d* if f is nonconstant on every d-dimensional affine subspace. By the above, f_n is an affine disperser for dimension d(n) for all sufficiently large n, with d(n) = o(n); the rest of the proof uses only this property, at the actual input length n.

## 3. Restrictions: rdq-sources and substitutions

**Definition (rdq-source).** Let X = {x_1, …, x_n} and let (F, L, Q) be a disjoint partition of X into *free*, *linear* and *quadratic* variables. An rdq-source is the partition together with equations:
1. for each x_j ∈ Q, one quadratic equation x_j = ((x_i ⊕ c_1) ∧ (x_k ⊕ c_2)) ⊕ c_3 with x_i, x_k free and the right-hand-side variables of distinct quadratic equations disjoint;
2. for each x_j ∈ L, one affine equation x_j = (⊕_{x_i ∈ F_j} x_i) ⊕ (⊕_{x_k ∈ Q_j} x_k) ⊕ c_1 with F_j ⊆ F, Q_j ⊆ Q.
The *dimension* of the source is |F|; it defines a subset of F_2^n of size 2^{|F|}. A free variable occurring on the right of a quadratic equation is *protected*; the two protected variables of one quadratic equation are *coupled* (each is the other's *couple*); other free variables are *unprotected*.

A circuit C *computes f restricting to R* (C computes f | R) if the linear and quadratic variables of R are 0-variables in C (Definition in §4 below) and C agrees with f on all inputs satisfying R's equations.

**Proposition 3.1.** An rdq-source of dimension 2d contains an affine subspace of dimension d; hence an affine disperser for dimension d is nonconstant when restricted to any rdq-source of dimension > 2d. *Proof.* There are 2d free variables and at most d quadratic equations (each consumes two free variables). For each affine equation involving two free variables of one quadratic equation, fix one to 0; the quadratic equations become affine, leaving an affine subspace of dimension ≥ d on which f is nonconstant. □

**Substitutions.** A substitution further restricts R by adding or modifying equations; redundant gates of C are then deleted by the normalization rules of §5. Following §2.4 of the source paper, there are four types (a *g-substitution* if the substituted value is g of the other variables):

1. *Constant substitution:* x_j ← c for a free variable x_j and constant c.
2. *Affine substitution to couples:* x_j ← x_k ⊕ c, where x_j, x_k are coupled protected variables. The quadratic equation they serve becomes affine, so its variable moves from Q to L and x_j becomes a linear variable.
3. *General affine substitution:* x_j ← x_k ⊕ c, where x_j is an *unprotected* free variable and x_k is a free variable. Valid because x_j, being unprotected, occurs on the right-hand side of no quadratic equation; the substitution unfolds x_j as x_k ⊕ c at all its occurrences in the affine equations, and x_j becomes a linear (affine) variable with equation x_j = x_k ⊕ c.
4. *Quadratic substitution:* x_j ← ((x_i ⊕ c_1) ∧ (x_k ⊕ c_2)) ⊕ c_3, where x_j, x_i, x_k are unprotected free variables and x_j can be made a 0-variable by appropriate rewiring; one quadratic equation is added and x_j becomes a quadratic variable (the right-hand variables here are the existing feeders named in the case, not fresh variables).

Each type eliminates exactly one free variable; hence, by Proposition 3.1, an affine disperser for dimension d is not trivial after n − 2d substitutions. After a substitution the newly non-free variable is disconnected by rewiring, so every variable of non-zero out-degree is free. Constants introduced by substitution have in-degree 0 and are removed by the normalization rules of §5; size counts gates, not wires, inputs or constants.

A variable is a *k-variable* if its out-degree in C is exactly k (see §4). An input is *influential* if it is a 1+-variable or protected.

**Semicircuits and xor-reconstruction.** Case 8 (§6) substitutes constants to *gates*, not only to variables; this uses the generalized model of FGHK16 (§2.5–2.6 of the source paper), fixed now.

*Definition (cyclic xor-circuit, semicircuit, fair).* A *cyclic xor-circuit* is a directed graph of ⊕-type gates and input variables, possibly with cycles, with arbitrarily many output gates. A *semicircuit* is the composition of a cyclic xor-circuit C_1 and an (acyclic) circuit C_2, the outputs of C_1 feeding C_2. If the xor-circuit has gates G_1, …, G_m computing f_i ∈ {⊕, ≡} of their two inputs, it corresponds to a linear system A g = b over F_2 (inputs viewed as constants); it is *fair* if A has full rank, i.e. the system has a unique solution for every assignment of the input variables. A semicircuit is fair if its cyclic xor-circuit is fair.

*Proposition (substitutions in fair semicircuits; Prop. 2.6 of the source paper).* Let C be a fair semicircuit and R an rdq-source. A g-substitution (of any of the four types above) to x_j that does not introduce cycles is valid: it produces an rdq-source R′ restricting R to x_j = g(·), and the resulting C′ is still a fair semicircuit computing f | R′. (Stated without proof in the source paper; the content is that a substitution only rewrites R's equations and the occurrences of x_j — an invertible change of the xor-circuit's linear system, preserving full rank, i.e. fairness — and by hypothesis creates no cycle.)

*Proposition (path property; Prop. 2.7 of the source paper).* Let I be a gate in a fair cyclic xor-circuit computing an affine function that depends on x_i. Then there is a path from x_i to I in the circuit. *Proof.* Let X be the set of gates unreachable from x_i. The subcircuit on X has inputs other than x_i and is itself a fair semicircuit, so its gate outputs are uniquely determined without x_i; since I depends on x_i, I ∉ X. □

*Xor-reconstruction.* Let G be a topologically minimal ∧-type gate fed by I_1, I_2 with I_1 a gate. The subcircuit computing I_1 contains no ∧-type gate, hence is a cyclic xor-circuit. To substitute I_1 ← b: if I_1 depends on a variable x_i, the path property gives a path x_i = I_0 → I_1 → … → I_k = I_1 with I_j = f_j(I_{j−1}, T_j), f_j ∈ {⊕, ≡}; each f_j is invertible in its first argument, so the system rewrites as I_{j−1} = g_j(I_j, T_j). Replacing I_k by the constant b, x_i by a new ⊕-type gate Z, and rewiring along the rewritten system yields the reconstructed circuit. *Proposition (Prop. 2.8 of the source paper).* This xor-reconstruction produces a fair semicircuit computing f | (I_1 = b): the new circuit's linear system is equivalent to the old one, with I_1's output additionally fixed to b. □ In the case files, "constant substitution to a gate I" means this operation, applied when I depends on an unprotected variable x_i; the restriction makes x_i a linear variable, so it is one substitution in the count t.

*Model transfer.* Every acyclic B2 circuit is a fair semicircuit (cyclic part empty), with the same gate count. The elimination below is therefore carried out for fair semicircuits — as in the source paper, we say "circuit" for "fair semicircuit" — and every bound proved for them applies in particular to the acyclic circuits counted by C_B2.

## 4. Gate taxonomy, troubled gates, packing, potential

A gate (or variable) is a *k-gate* (resp. *k-variable*) if its out-degree in C is exactly k; k+ means out-degree at least k, and a *1+-variable* has out-degree ≥ 1. By functionality (§2.1 of the source paper), the 16 kinds of B2 gates fall into four classes: *trivial* (constant function: f(x, y) = c); *degenerate* (depends on exactly one input: f(x, y) = (x ∧ c_1) ⊕ (y ∧ ¬c_1) ⊕ c_2 — outputs one input or its negation); *⊕-type* (affine, depends on both inputs: f(x, y) = x ⊕ y ⊕ c); and *∧-type* (quadratic, depends on both inputs: f(x, y) = ((x ⊕ c_1) ∧ (y ⊕ c_2)) ⊕ c_3). A gate that *becomes* trivial or degenerate after a substitution (e.g. by being fed a constant) is called trivialized or degenerate in the same sense; Rules 2 and 3 of §5 eliminate such gates when fed by a constant.

**Definition (troubled gate).** A gate G is *troubled* if it is an ∧-type 1-gate fed by two 2-variables. Two troubled gates are *adjacent* if they share an input variable.

**Definition (troubled pair, packing).** A *troubled pair* is a pair of adjacent troubled gates; the variables feeding both are its *inner* variables, those feeding exactly one are *boundary* variables; a pair with two inner variables and no boundary variable is *compact*. A *packing* P is a set of troubled pairs pairwise disjoint in gates; a pair in P is a *pack*. A troubled gate lying in no pack of P is *truly troubled*.

**Definition (potential).** Φ(C, P) = (number of truly troubled gates) + |P| = (number of troubled gates) − |P|.

If a truly troubled gate were introduced by unpacking a pack (G, G′), its mate G′ would have to become non-troubled unless the inner variables no longer feed G′; Proposition 5.6 in §5 (Proposition 3.10 of the source paper) shows a normalization rule can never destruct a pack while leaving both mates troubled, so below the potential can grow only through directly introduced truly troubled gates or packs.

## 5. The measure, the potential bound, normalization, pair elimination

**Definition (measure).** For positive constants α_I, α_Q, α_φ fixed in §7,

    μ(C, P, R) = g + α_I·i + α_Q·q + α_φ·Φ(C, P),

where g is the number of gates, i the number of influential inputs (1+-variables or protected variables), and q the number of quadratic equations of R.

**Lemma 5.1 (potential bound).** If C computes an affine disperser for dimension d, then for every packing P, Φ(C, P) ≤ n/2 + 5d/2. *Proof.* Φ(C, P) = #troubled − |P| ≤ #troubled, so bound t, the number of troubled gates; this is Lemma 3 of FGHK16 (Lemma 3.7 of the source paper), whose proof we reproduce. Call two inputs *neighbors* if they feed a common troubled gate. A variable feeding a troubled gate is a 2-variable, so if v_i counts the variables feeding exactly i troubled gates, v_i = 0 for i > 2. Double-counting the wires from inputs to troubled gates gives 2t = v_1 + 2v_2, and v_1 + v_2 ≤ n. Suppose t ≥ n/2 + 5d/2; then n + 5d ≤ 2t = v_1 + 2v_2 ≤ n + v_2, so v_2 ≥ 5d. Let T be the set of inputs feeding two troubled gates, |T| = v_2 ≥ 5d. Construct disjoint X ⊆ T and Y ⊆ V = {inputs}, with |X| = d and |Y| consistent linear equations (on variables of Y) making C independent of the variables of X: repeat d times — pick any x ∈ T; x feeds two troubled gates, whose other feeders are y_1, y_2 (possibly y_1 = y_2); add x to X and y_1, y_2 to Y. Fixing y_1, y_2 to constants eliminates both troubled gates fed by x and leaves C independent of x (if y_1 = y_2, one of x = c, y_1 = c, or y_1 = x ⊕ c eliminates both gates for some constant c; in the x = c case, interchange the roles of x and y_1). Each of y_1, y_2 has at most one neighbor different from x; remove x, y_1, y_2 and the other neighbors of y_1, y_2 — at most five variables — from T. Since |T| ≥ 5d the routine runs d times, and since those neighbors are removed from T, the neighbors of every variable picked later lie outside the Y chosen earlier, so all the assignments stay mutually consistent. Now impose the |Y| equations and set every remaining variable in V ∖ (X ∪ Y) to a constant: C then evaluates to a constant, since it does not depend on X and everything else is fixed. The number of independent linear equations is |Y| + |V ∖ (X ∪ Y)| = n − d, i.e. f is constant on an affine subspace of dimension d, contradicting the disperser property. Hence t < n/2 + 5d/2, and Φ(C, P) ≤ n/2 + 5d/2. □

**Corollary 5.2 (size from measure).** If every circuit C computing f with every packing P satisfies μ(C, P, ∅) ≥ μ_0, then |C| ≥ μ_0 − α_I·n − α_φ·(n/2 + 5d/2); in particular if f is an affine disperser for sublinear dimension d = o(n) and μ_0 = δn − o(n), then |C| ≥ (δ − α_I − α_φ/2)·n − o(n).

**Normalization.** After substitutions, redundant gates are removed by the five rules of §3.3 of the source paper. A circuit to which no rule applies is *normalized*.

*Rule 1.* Let G be a 0-gate fed by I_1, I_2. Remove G; the out-degrees of I_1, I_2 drop by one. A new troubled gate can only be I_1, I_2, or a gate fed by one of them; troubled gates are fed by 2-variables, so I_1 (resp. I_2) feeds at most two troubled gates, and two fed by the same input are packed together. Hence each input of G produces at most a troubled gate or a pack: ΔΦ ≤ 2, and the step's measure decrement is β ≥ 1 − 2α_φ.

*Rule 2.* Let G be a trivialized gate fed by a constant b. Remove G; its descendants are now fed by a constant c. If the other input I_1 is not a constant, its out-degree drops by one and a new troubled gate is either I_1 or a gate fed by I_1, so as in Rule 1, ΔΦ ≤ 1 and β ≥ 1 − α_φ.

*Rule 3.* Let G be a degenerate gate fed by a constant b and another node I_1. Remove G and rewire so G's descendants are fed directly by I_1 (if G computes the negation of I_1, the descendants' functions are modified accordingly; if G is the output gate, I_1 becomes the output gate with functions adjusted). As in Rule 2, ΔΦ ≤ 1 and β ≥ 1 − α_φ.

*Rule 4.* A gate G is *useless* if it is a 1-gate with only descendant Q, such that another input of Q also feeds G. Eliminate G and rewire the nodes feeding it to feed Q. If one input I_2 of G is a gate, it is the only gate that may become troubled; if I_2 is a variable, new troubled gates must be fed by it, forming a troubled gate or a pack. Hence ΔΦ ≤ 1 and β ≥ 1 − α_φ.

*Rule 5.* Let G be a gate whose two inputs coincide. Make it degenerate (or trivial) by changing one (or both) of its inputs to a constant, and eliminate it via Rule 2 or Rule 3. Hence ΔΦ ≤ 1 and β ≥ 1 − α_φ.

**Proposition 5.6 (pack destruction; Proposition 3.10 of the source paper).** If a pack (G_1, G_2) is destructed after applying a normalization rule, at least one of G_1, G_2 becomes non-troubled; hence potential increases can only be caused by direct introduction of troubled gates or packs. *Proof.* Suppose a rule destructs (G_1, G_2) while both remain troubled. Destructing the pack requires the descendants of the pair's inner variable to change; but the normalization rules change only the descendants of the gates feeding the eliminated gate. Hence one of G_1, G_2 is itself the eliminated gate — a contradiction. □

**Lemma 5.3 (normalization lemma; Lemma 3.11 of the source paper).** Let β be the measure decrement of one normalization step. For Rule 1, β ≥ 1 − 2α_φ; for Rules 2–5, β ≥ 1 − α_φ. Moreover, if the eliminated gate G is fed by I_1, I_2, every potential increase in the step is a directly introduced troubled gate or pack (Proposition 5.6 rules out increases by pack destruction), and each introduced troubled gate or pack is I_1 or I_2 or is fed by one of them after the step, as verified rule by rule above. As α_φ < 1/2, β ≥ 0: normalization never hurts, so C may be assumed normalized whenever a case of §6 begins.

**Terminal-output convention.** The case analysis of §6 repeatedly eliminates gates at or near the output of the circuit, so we record once how the output gate is treated; the source paper uses the same convention (its Rule 3 states the exception in (T2), and its §4.1 preamble and Case 6.1 argue as in (T1) and (T3)).

(T1) *The output gate is never trivialized.* Throughout the elimination, the current circuit computes f | R on a source R of dimension > 2d, and by Proposition 3.1 the function f | R is not constant. Hence the output gate never computes a constant function, and no operation in Cases 0–8 trivializes it: whenever a case trivializes a gate, that gate is not the output gate. (This also legitimizes Case 0.3: a constant-output gate is never the output gate.)

(T2) *A degenerate output gate is eliminated by Rule 3 as usual.* If the output gate becomes degenerate, Rule 3 applies under the exception stated with it: its remaining input I_1 becomes the output gate, with the functions of I_1 and its descendants adjusted for a possible negation. This counts as one eliminated gate with ΔΦ ≤ 1, exactly as the case ledgers record for any Rule 3 elimination.

(T3) *Rule 1 never removes the output gate.* Rule 1 is applied in the case files only to a gate whose out-degree has dropped to 0 during the operation; the output gate has out-degree 0 from the start and is never removed by Rule 1. In every invocation in `cases-a.md` and `cases-b.md`, the gate removed by Rule 1 has a descendant named in the configuration at hand — hence is not the output gate — or does not reach the output gate (Case 0.2). Examples: in Case 5.4.2.2 the gate B removed by Rule 1 has the descendant E; in Case 5.4.1.4, step (2), the gate G removed by Rule 1 feeds D; in Case 7 the gate Q eliminated by Rule 1 feeds the minimal gate G. Case 7 is in any event the dedicated treatment of the output gate when its inputs are a 2-variable and a 1-gate.

**Lemma 5.4 (troubled pair elimination).** Let C compute f | R with a non-empty packing P, and suppose f is not trivial under two substitutions. Then substitutions can be performed whose measure decrement per substitution is at least δ_p = min{ (3/2)α_I, 3 − 3α_φ + α_I + α_Q }. *Proof.* Take a pack (G_1, G_2), inner variable y, boundary variables x, z. If the pack is not compact: if y is protected, substitute a constant to y trivializing G_1; G_2 and the only descendants of G_1 become degenerate; normalization removes three gates with ΔΦ ≤ 3, and the substitution kills one quadratic equation and one influential input, so Δμ ≥ 3 − 3α_φ + α_I + α_Q. If y is unprotected, substitute constants to x and z trivializing both gates; all three variables become non-influential, so Δμ ≥ 3α_I over two substitutions, i.e. (3/2)α_I each. If the pack is compact, fed by x and y: if one is protected, argue as in the first subcase; otherwise both are unprotected, and enumerating the two gates' functions shows a constant substitution to x or y, or the affine substitution y ← x ⊕ b, trivializes both, giving Δμ ≥ 2α_I > (3/2)α_I for one substitution. □

At the parameters of §7, δ_p = min{72/5, 69/5} = 69/5 > δ = 64/5, so pair elimination is never the bottleneck.

## 6. The elimination theorem

**Theorem 6.1.** Let C compute an affine disperser f for dimension d, let P be any packing, and set

    δ = α_I + min{ α_I/3, 2 − 2α_φ + α_Q, 4 − 4α_φ, 3 + α_φ, 5 − α_Q,
                   (5 − 2α_φ + α_Q)/2 }.

If the dimension of the ambient rdq-source is > 2d + 2, then either f becomes trivial after at most two substitutions (impossible while the dimension condition holds, by Proposition 3.1), or there is a sequence of 1 ≤ t ≤ 3 substitutions after which μ drops by at least δ per substitution: μ(C′, P′, R′) ≤ μ(C, P, R) − δ·t for the resulting circuit, source and a suitable packing P′.

The proof is the case analysis in `cases-a.md` (Cases 0–5) and `cases-b.md` (Cases 6–8), ordered: each case assumes all earlier configurations absent. Case 0 disposes of unnormalized circuits, constant-output gates and non-empty packings (via Lemma 5.4). Cases 1–4 dispose of protected variables with high out-degree, protected 0-variables, and high-out-degree ∧-type configurations. Case 5 and its subcases handle the troubled configuration itself (an ∧-type gate fed by two 2-variables), distinguishing whether and where new troubled gates appear. Case 6 handles substituting a constant to a protected variable, together with the configurations left by its exceptions. Case 7 handles the output gate when its inputs are a 2-variable and a 1-gate. Case 8 handles the rest: both inputs of the minimal ∧-type gate are gates — via the xor-reconstruction of §3 — or a 2-gate remains; the deepest subcases also use quadratic substitutions. Every leaf records t, the gates eliminated, and the bounds on ΔΦ, Δi, Δq, hence the per-substitution decrement; the minimum over all leaves gives exactly the six quantities in δ (the term 2 − 2α_φ + α_Q is attained by Cases 1 and 6.1 and is not the minimum at the chosen parameters).

Iterating Theorem 6.1 from (C, P, ∅) while the dimension exceeds 2d + 2: the dimension starts at n and each substitution lowers it by at most 1, so the process runs for at least n − 2d − 2 substitution steps before it can stop, and f never becomes trivial on a source of dimension > 2d (Proposition 3.1). Hence

    μ(C, P, ∅) ≥ δ·(n − 2d − 2).

## 7. Parameters and the final assembly

Maximize δ − α_I − α_φ/2 over α_I, α_Q, α_φ ≥ 0 with α_φ ≤ 1/2, where δ is as in Theorem 6.1. Writing u = δ − α_I, the constraints are u ≤ α_I/3, u ≤ 2 − 2α_φ + α_Q, u ≤ 4 − 4α_φ, u ≤ 3 + α_φ, u ≤ 5 − α_Q, u ≤ (5 − 2α_φ + α_Q)/2. Take

    α_I = 48/5,   α_Q = 9/5,   α_φ = 1/5.

Then the six quantities are 16/5, 17/5, 16/5, 16/5, 16/5, 16/5, so u = 16/5 and δ = 48/5 + 16/5 = 64/5. Optimality: for any feasible parameters, take the nonnegative combination of the three constraints u ≤ 3 + α_φ, u ≤ 5 − α_Q and u ≤ (5 − 2α_φ + α_Q)/2 with weights 7/10, 1/10 and 1/5 respectively:

    (7/10)(3 + α_φ) + (1/10)(5 − α_Q) + (1/5)·((5 − 2α_φ + α_Q)/2)
      = 21/10 + 1/2 + 1/2 + (7/10 − 1/5)α_φ + (−1/10 + 1/10)α_Q
      = 31/10 + (1/2)α_φ.

Hence u ≤ 31/10 + α_φ/2 for every feasible choice, so the objective u − α_φ/2 is at most 31/10, and the displayed parameters attain it.

Finally, Corollary 5.2 with μ_0 = δ(n − 2d − 2), d = d(n) = o(n):

    |C| ≥ δ(n − 2d − 2) − α_I·n − α_φ·(n/2 + 5d/2)
        = (64/5 − 48/5 − 1/10)·n − o(n)
        = (32/10 − 1/10)·n − o(n)
        = (31/10)·n − o(n).

Since C was an arbitrary B2 circuit computing f_n exactly, C_B2(f_n) ≥ (31/10)n − o(n) for all large n, and therefore liminf C_B2(f_n)/n ≥ 31/10. ∎
