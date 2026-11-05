# The (31/10, 1) lower bound for the BK/Shoup family: complete proof

**Theorem.** For the single total family (f_n) of §1, in the full Li–Yang B2 model (deterministic acyclic circuits, one output, two input pins per gate, all 16 binary Boolean truth tables, unrestricted depth and fan-out, exact computation on every input, nonuniform circuits allowed),

    liminf_{n→∞} C_B2(f_n) / n ≥ 31/10,   i.e. C_B2(f_n) ≥ (31/10)n − o(n)

for every sufficiently large integer n. Proof: Li–Yang's gate elimination (ECCC TR21-023, Theorem 1.1) reproduced in full, §4.1 case analysis in `elimination.md`. The one premise v1 cited externally — the affine-disperser property (Ben-Sasson–Kopparty) — is proved here: §2 reduces it to BK Theorem 2.2 with every hypothesis verified; `disperser.md` reproduces BK Sections 3, 6, 7 in full. Used without proof only: Bertrand's postulate (a prime between x and 2x; once, §2) and basic finite-field facts (the trace fact, Proposition 3.3, is proved in `disperser.md`). The one construction earlier versions cited externally — Shoup's deterministic irreducible-polynomial construction — is reproduced in full in §1A, specialized to p = 2 and prime degree; it fixes q_k in §1 and supplies the evaluator's runtime in §2. Lemma 5.1 (FGHK16) is proved in §5. Parameters: α_I = 48/5, α_Q = 9/5, α_φ = 1/5, δ = 64/5.

## 1. The fixed total family

f_0(empty) = 0. For n ≥ 1: k = least prime strictly above 2·n^{3/5}; r = ⌈n/k⌉; t = ⌈√r⌉. Fix the monic irreducible q_k in F_2[Z] of degree k output by the
deterministic construction of §1A (Shoup 1990, specialized to p = 2 and
prime degree k, reproduced and proved there); K = F_2[Z]/(q_k),
basis 1, Z, ..., Z^{k-1}.

**Construction S (Shoup, specialized).** The polynomial q_k is the
output of the deterministic algorithm of §1A (V. Shoup, *New algorithms
for finding irreducible polynomials over finite fields*, Mathematics of
Computation 54 (1990), 435--447, specialized to p = 2 and prime degree
k), which §1A reproduces and proves to run in time polynomial in k. The
algorithm is deterministic, so q_k is a single definite polynomial for
each k: no per-length advice and no randomness anywhere in the family.
Which monic irreducible is output is immaterial to everything else in
this proof: K is a field with 2^k elements for every such choice; the
identifications of §2(ii) (the block embedding phi, and pi = constant
coefficient with respect to the power basis) are F_2-linear and
admissible for every choice, and Theorem 2.2 of [BK12] quantifies over
all such phi and pi; the exponent and Frobenius structure is
field-theoretic. Hence the only content of Construction S used below is:
a definite q_k exists for every k, and it can be computed from k in time
polynomial in k — both proved in §1A.

(The code artifact `family.py` does not implement the construction of
§1A; for its finite evaluations it fixes q_k by the numerically smallest
monic irreducible, as a concrete stand-in. No audited quantity depends
on the choice of irreducible, and the polynomial-time premise for the
evaluator rests on §1A, not on `family.py`.) One language, no per-length advice. Pad x ∈ {0,1}^n with rk − n zeros, split into r blocks of k bits, read block i as φ_i(x) ∈ K, first bit = constant coefficient; φ = (φ_i) is injective F_2-linear. π: K → F_2 = constant coefficient: F_2-linear, π(1) = 1. y_i = φ_i(x)^{2^{i+1} − 1} (Frobenius powers in K), 1 ≤ i ≤ r.

    f_n(x) = π( Sym_r^t(y_1, …, y_r) ),

Sym_r^t = elementary symmetric polynomial of degree t (sum of products over t-subsets of {1, …, r}).

## 1A. The deterministic construction of q_k (Shoup, specialized)

We reproduce the case p = 2, k an odd prime, of Shoup's deterministic
construction (Mathematics of Computation 54 (1990), §§2–3), with the
factoring subroutine made explicit (Berlekamp's algorithm over F_2,
Lemma 1A.5). For k = 2, take q_2 = Z^2 + Z + 1; in §1, k is the least
prime above 2·n^{3/5}, hence odd for every n ≥ 1. Nothing in the argument
uses a bound on a least nonresidue or a least irreducible polynomial,
character sums, or any unproved hypothesis: the one nonresidue the
construction needs is manufactured as a root of unity and recognized by
Lagrange's theorem (Lemma 1A.2). This is the only known deterministic
route with that feature — the alternatives either search for a
nonresidue (deterministic only under ERH: Adleman–Lenstra 1986) or
prescribe coefficients of an irreducible (all known proofs use Weil
bounds) — which is why this construction is the one reproduced.

**The algorithm.** Let m = ord_k(2), the multiplicative order of 2
modulo k; m | k − 1 (Fermat), so m < k. Write 2^m − 1 = l·k^s with
k ∤ l; s ≥ 1 since k | 2^m − 1.

(A) Factor Φ_k = Z^{k−1} + ⋯ + 1 over F_2 (Φ_k = (Z^k − 1)/(Z − 1), k
prime) and let f⁽¹⁾ be any irreducible factor.

(B) For i = 2, …, s, factor f⁽ⁱ⁻¹⁾(Z^k) over F_2 and let f⁽ⁱ⁾ be any
irreducible factor. Put f = f⁽ˢ⁾, K = F_2[Z]/(f), θ = Z mod f.

(C) Form E = K[Y]/(Y^k − θ) and let α = Y mod (Y^k − θ).

(D) Let H be the unique subfield of E with 2^k elements, and put
β = Σ_{i=0}^{m−1} α^{2^{ik}} (the trace of α from E to H, see Lemma
1A.4). Output the minimal polynomial of β over F_2, computed by
F_2-linear algebra in the F_2-space E.

**Lemma 1A.1 (factor degrees).** Every irreducible factor of Φ_k over
F_2 has degree m. For 1 ≤ j ≤ s the roots of f⁽ʲ⁾ are primitive kʲ-th
roots of unity, f⁽ʲ⁾ has degree m, and for 2 ≤ j ≤ s the polynomial
f⁽ʲ⁻¹⁾(Z^k) is a product of distinct irreducible polynomials over F_2,
all of degree m. In particular θ is a primitive k^s-th root of unity
in K.

*Proof.* The roots of Φ_k are the primitive k-th roots of unity; for
such a root ζ, [F_2(ζ):F_2] is the least d with ζ ∈ GF(2^d), i.e. with
k | 2^d − 1, i.e. with ord_k(2) | d; so every irreducible factor of Φ_k
has degree m, giving the claim for j = 1. Let j ≥ 2 and suppose the
roots of f⁽ʲ⁻¹⁾ are primitive k^{j−1}-th roots of unity. A root ζ of
f⁽ʲ⁻¹⁾(Z^k) has ζ^k a primitive k^{j−1}-th root; writing o = ord(ζ),
ord(ζ^k) = o/gcd(o, k) = k^{j−1} forces o = kʲ (o = k^{j−1} would give
ord(ζ^k) = k^{j−2} ≠ k^{j−1}). So every root of f⁽ʲ⁻¹⁾(Z^k) is a
primitive kʲ-th root of unity. Its F_2-degree is ord_{kʲ}(2): since
j ≤ s, kʲ | 2^m − 1, so ord_{kʲ}(2) | m; and k | kʲ gives
ord_k(2) | ord_{kʲ}(2); hence ord_{kʲ}(2) = m. Finally
f⁽ʲ⁻¹⁾(Z^k) is separable: it divides Z^{kʲ} − 1 (every root has order
dividing kʲ), and Z^{kʲ} − 1 is separable because kʲ is odd (its
derivative Z^{kʲ−1} shares no root with it). So its irreducible
factors are distinct and all of degree m, and any one of them, f⁽ʲ⁾,
has primitive kʲ-th roots. ∎

**Lemma 1A.2 (the manufactured nonresidue).** θ is not a k-th power
in K.

*Proof.* θ has order exactly k^s in the cyclic group K*. If θ = γ^k
for some γ ∈ K*, then ord(γ^k) = ord(γ)/gcd(ord(γ), k) = k^s forces
k^{s+1} | ord(γ). By Lagrange's theorem ord(γ) divides
|K*| = 2^m − 1 = l·k^s, which k^{s+1} does not divide — contradiction.
∎

**Lemma 1A.3 (prime-case Kummer criterion).** Let K be a field of
characteristic ≠ k containing a primitive k-th root of unity ζ, k
prime, and let a ∈ K* not be a k-th power in K. Then Y^k − a is
irreducible over K.

*Proof.* Let α be a root in an algebraic closure; the roots of
Y^k − a are α·ζ^j, j = 0, …, k − 1, all distinct (k ≠ characteristic).
Suppose Y^k − a = g·h with g monic of degree d, 1 ≤ d ≤ k − 1. The
constant term of g is, up to sign, a product of d of the roots, i.e.
c·α^d with c a product of powers of ζ; c lies in the group of roots
of unity of K (ζ ∈ K), so α^d ∈ K. Since gcd(d, k) = 1, write
ud + vk = 1 with integers u, v; then
α = (α^d)^u·(α^k)^v = (α^d)^u·a^v ∈ K, so a = α^k is a k-th power in
K — contradiction. ∎

Apply Lemma 1A.3 with a = θ: K = GF(2^m) contains a primitive k-th
root of unity (k | 2^m − 1), and θ ∉ K^k by Lemma 1A.2; hence
Y^k − θ is irreducible over K, [E : K] = k and [E : F_2] = mk.

**Lemma 1A.4 (trace descent).** β ∈ H, and β has degree exactly k
over F_2. Hence the output of step (D) is a monic irreducible
polynomial of degree k over F_2.

*Proof.* Subfields of E = GF(2^{mk}) correspond to divisors of mk, so
H = GF(2^k) exists and is unique. The Galois group of E over H is
generated by x ↦ x^{2^k}, of order m, so β = Σ_{i=0}^{m−1} α^{2^{ik}}
is the trace of α from E to H and lies in H. Hence
[F_2(β):F_2] divides k, so it is 1 or k. Suppose β ∈ F_2 ⊆ K. Write
2^{ik} = x_i·k + y_i with 0 ≤ y_i < k; then
α^{2^{ik}} = (α^k)^{x_i}·α^{y_i} = θ^{x_i}·α^{y_i}, so
β = Σ_i θ^{x_i}·α^{y_i}. The residues y_i are pairwise distinct:
y_i = (2^k)^i mod k, and 2^k ≡ 2 (mod k) by Fermat's little theorem,
so y_i = 2^i mod k, pairwise distinct for i = 0, …, m − 1 because
ord_k(2) = m; in particular y_0 = 1 and its coefficient is
θ^{x_0} = θ^0 = 1. Thus α is a root of the polynomial
Σ_i θ^{x_i}·Y^{y_i} − β over K, which is nonzero (its Y^1 coefficient
is 1) and of degree at most k − 1 < k, contradicting Lemma 1A.3:
Y^k − θ is the minimal polynomial of α over K. So β ∉ K,
[F_2(β):F_2] = k, and the minimal polynomial of β is monic irreducible
of degree k. No separate representation of H is needed: β is computed
inside E, and its minimal polynomial is found by linear algebra over
F_2 in the mk-dimensional space E. ∎

**Lemma 1A.5 (Berlekamp factoring over F_2).** There is a
deterministic algorithm that, given a squarefree g ∈ F_2[Z] of degree
D, outputs the irreducible factorization of g in time polynomial in
D. (Squarefree inputs suffice: the construction factors only Φ_k,
which divides the separable polynomial Z^k − 1 with k odd, and the
polynomials f⁽ⁱ⁻¹⁾(Z^k), separable by Lemma 1A.1.)

*Proof.* Two classical steps. (i) *Distinct-degree factorization:*
proceed inductively on d = 1, 2, …: having removed from g all
irreducible factors of degree < d, gcd(g, Z^{2^d} − Z) is exactly the
product of the irreducible factors of g of degree d. Indeed an
irreducible h of degree e divides Z^{2^d} − Z iff every root of h
lies in GF(2^d), iff e | d; and after the removal, every remaining
factor has degree ≥ d, so e | d forces e = d. Compute Z^{2^d} mod g
by repeated squaring. (ii) *Equal-degree splitting:* let g be a
product of r ≥ 2 distinct irreducibles, all of the same degree.
Consider the algebra A = {h mod g : h² ≡ h (mod g)} of
Frobenius-fixed points of F_2[Z]/(g). By the Chinese remainder
theorem F_2[Z]/(g) ≅ Π_{j=1}^r GF(2^{d_j}), Frobenius acts
factorwise as the field Frobenius, whose fixed field in GF(2^{d_j})
is F_2; hence A ≅ F_2^r as an F_2-algebra, of dimension r. Compute a
basis of A as the kernel of Q − I, where Q is the matrix of the
Frobenius map h ↦ h² on F_2[Z]/(g) in the power basis. Any
h ∈ A ∖ {0, 1} corresponds to a {0,1}-tuple not all equal, so one of
gcd(g, h), gcd(g, h + 1) is a proper nontrivial factor of g; recurse
on the factors. A nonconstant h exists because dim A = r ≥ 2. Every
step is a gcd or linear algebra over F_2 in dimension at most D:
polynomial in D, and deterministic (the only choices are basis
vectors of a computed kernel). ∎

**Runtime tally.** The polynomials factored in steps (A)–(B) have
degree at most mk ≤ k(k − 1) (Φ_k has degree k − 1), and there are
s + 1 ≤ k factoring calls; by Lemma 1A.5 each costs polynomially many
F_2-operations in its degree — O(D³) by the linear algebra stated
there, so O(k⁷) in total with crude bounds. Arithmetic in K and E is
polynomial arithmetic modulo polynomials of degree at most mk; the
trace of step (D) is m ≤ k Frobenius powers, and the minimal
polynomial of β is one linear-algebra computation in dimension
mk ≤ k². Every step is effective and the total is polynomial in k.
This proves Construction S: q_k is definite and computable from k in
time polynomial in k. ∎

## 2. Uniformity, P and NP, the disperser property

*Evaluator A (all lengths):* trial-division scan for the least prime k above the threshold (polynomial: k = O(n^{3/5}) by Bertrand); q_k computed by the construction of §1A (deterministic; time polynomial in k; and k = O(n^{3/5}), so polynomial in n); field arithmetic in K at polynomial bit cost per operation (elements are k-bit vectors; multiplication is polynomial multiplication followed by reduction modulo q_k); Sym by the descending recurrence E_0 = 1, E_j ← E_j + y_i·E_{j−1} after reading y_i, j = min(i, t) down to 1 — no characteristic-two division; splitting subsets by membership of i proves the recurrence. Powers y_i: O(i+1) multiplications each; recurrence: O(r·t) field operations, polynomial bit cost; output π(E_t). Hence (f_n) ∈ P. V(x, w) = 1 iff w empty and A(x) = 1, W ≡ 0 ⇒ NP membership both directions immediate: in P, hence NP, empty witnesses. Witnesses are not inputs to the target circuits.

*Disperser property.* §1 is exactly the construction of Theorem 2.2, equation (4), of Ben-Sasson–Kopparty, *Affine Dispersers from Subspace Polynomials* ([BK12] in the source paper); proof reproduced in `disperser.md` (their Sections 3, 6, 7). Hypotheses: (i) *Parameters:* p = 2; input length m := n (actual length of f_n, not padded rk); block length = least prime > 2·m^{3/5} = our k; blocks r = ⌈m/k⌉ = our r; degree t = ⌈√r⌉ = our t. (ii) *Map:* Theorem 2.2 quantifies over arbitrary injective F_2-linear φ: F_2^m → (F_{2^k})^r and arbitrary nontrivial F_2-linear π: F_{2^k} → F_2. Ours: padding + block split + per-block basis identification = coordinate embedding then fixed identifications; F_2-linear, injective, defined since n ≤ rk; under K ≅ F_{2^k}, admissible. No divisibility of n by k assumed or needed: the short last block is absorbed by φ, its last component's top rk − n coordinates constantly 0. π = constant coefficient: F_2-linear, π(1) = 1, nontrivial, admissible; Proposition 3.3 (`disperser.md`): π(x) = Tr(µx) for a unique µ ≠ 0 — the form the proof uses. (iii) *Exponents:* equation (4) raises block i to 1 + p + … + p^i = 2^{i+1} − 1 at p = 2 = the exponent of y_i in §1, i = 1, …, r. (BK Section 7 indexes exponents in reverse; `disperser.md` keeps that order and identifies the forms in the deduction of Theorem 2.2 by renaming blocks.)

*Conclusion:* **for every sufficiently large integer n, f_n is nonconstant on every affine subspace of F_2^n of dimension strictly greater than 6·n^{4/5}.** Only the parameter deduction (Theorem 7.3 ⇒ Theorem 2.2; arithmetic at the end of `disperser.md`) uses our parameters: Bertrand ⇒ k ≤ 4n^{3/5}; r(r+1)/t ≤ (r+1)^{3/2} ≤ (n^{2/5}/2 + 2)^{3/2} < 2n^{3/5} < k; threshold nr/t + nt + r(r+1) ≤ about 4.25·n^{4/5} < 6n^{4/5}, estimating r ≈ n/k jointly with k, not maximized independently. Small n: n ≤ 6^5 = 7776 ⇒ 6n^{4/5} ≥ n, no subspace exceeds the threshold, vacuous; and the claim is a liminf at large lengths. Set d(n) = ⌊6·n^{4/5}⌋ + 1; eventually 1 ≤ d(n) ≤ n, d(n)/n → 0.

**Definition.** f: F_2^n → F_2 is an *affine disperser for dimension d* if nonconstant on every d-dimensional affine subspace. So f_n is an affine disperser for dimension d(n), all large n, d(n) = o(n). Below, only this property is used, at actual length n.

## 3. Restrictions: rdq-sources, substitutions, semicircuits

**Definition (rdq-source).** X = {x_1, …, x_n} partitioned disjointly into free F, linear L, quadratic Q variables, with equations: (1) each x_j ∈ Q has one quadratic equation x_j = ((x_i ⊕ c_1) ∧ (x_k ⊕ c_2)) ⊕ c_3, x_i, x_k free, RHS variables of distinct quadratic equations disjoint; (2) each x_j ∈ L has one affine equation x_j = (⊕_{x_i ∈ F_j} x_i) ⊕ (⊕_{x_k ∈ Q_j} x_k) ⊕ c_1, F_j ⊆ F, Q_j ⊆ Q. *Dimension* = |F|; the source defines 2^{|F|} inputs. A free variable on the RHS of a quadratic equation is *protected*; the two protected variables of one equation are *coupled* (each the other's *couple*); other free variables are *unprotected*. C *computes f | R* if R's L and Q variables are 0-variables in C (§4) and C agrees with f on all inputs satisfying R's equations.

**Proposition 3.1.** A dimension-2d rdq-source contains a d-dimensional affine subspace; hence an affine disperser for d is nonconstant on any rdq-source of dimension > 2d. *Proof.* 2d free variables, ≤ d quadratic equations (two free variables each). For each affine equation involving both variables of one quadratic equation, fix one to 0: the quadratic equations become affine, a subspace of dimension ≥ d remains, f nonconstant on it. □

**Substitutions** (source §2.4; *g-substitution* if the value is g of the other variables) further restrict R; redundant gates are then deleted by §5 normalization. Four types: (1) *Constant:* x_j ← c, x_j free. (2) *Affine to couples:* x_j ← x_k ⊕ c, x_j, x_k coupled protected; their quadratic equation becomes affine — its variable moves Q → L and x_j becomes linear. (3) *General affine:* x_j ← x_k ⊕ c, x_j unprotected free, x_k free; valid since unprotected x_j occurs on no quadratic RHS; unfold x_j = x_k ⊕ c at all occurrences in affine equations, x_j becomes a linear variable. (4) *Quadratic:* x_j ← ((x_i ⊕ c_1) ∧ (x_k ⊕ c_2)) ⊕ c_3, x_j, x_i, x_k unprotected free, x_j made a 0-variable by rewiring; one quadratic equation added, x_j becomes quadratic (RHS variables are the existing feeders named in the case, not fresh variables).

Each type eliminates exactly one free variable ⇒ by Proposition 3.1 an affine disperser for d is not trivial after n − 2d substitutions. The newly non-free variable is disconnected by rewiring, so every non-zero-out-degree variable is free. Substitution constants have in-degree 0, removed by §5 normalization; size counts gates only, not wires, inputs or constants. A *k-variable* has out-degree exactly k (§4); an input is *influential* if a 1+-variable or protected.

**Semicircuits** (FGHK16 model, source §§2.5–2.6; needed since Case 8 substitutes constants to gates). *Cyclic xor-circuit:* directed graph of ⊕-type gates and input variables, cycles allowed, arbitrarily many outputs. *Semicircuit:* a cyclic xor-circuit C_1 composed with an acyclic circuit C_2, C_1's outputs feeding C_2. Gates G_1, …, G_m of C_1 computing f_i ∈ {⊕, ≡} of two inputs ⇔ a linear system A g = b over F_2 (inputs as constants); *fair* iff A has full rank, i.e. a unique solution for every input assignment; a semicircuit is fair iff its xor-circuit is.

*Proposition (source 2.6).* C a fair semicircuit, R an rdq-source: a g-substitution (any of the four types) to x_j introducing no cycles is valid — it yields R′ = R restricted to x_j = g(·), and C′ is a fair semicircuit computing f | R′. (Stated without proof in the source; the content is that only R's equations and occurrences of x_j are rewritten — an invertible change of the linear system, preserving full rank, i.e. fairness — and no cycle arises by hypothesis.)

*Proposition (source 2.7, path property).* A gate I of a fair cyclic xor-circuit computing an affine function depending on x_i has a path from x_i to I. *Proof.* X = gates unreachable from x_i: the subcircuit on X has inputs other than x_i and is itself a fair semicircuit, so its outputs are uniquely determined without x_i; I depends on x_i ⇒ I ∉ X. □

*Xor-reconstruction.* G a topologically minimal ∧-type gate fed by I_1, I_2, I_1 a gate; I_1's subcircuit contains no ∧-type gate, hence is a cyclic xor-circuit. To substitute I_1 ← b: if I_1 depends on a variable x_i, the path property gives a path x_i = I_0 → I_1 → … → I_k = I_1 with I_j = f_j(I_{j−1}, T_j), f_j ∈ {⊕, ≡}, each invertible in its first argument; rewrite the system as I_{j−1} = g_j(I_j, T_j). Replace I_k by the constant b, x_i by a new ⊕-type gate Z, and rewire along the rewritten system. *Proposition (source 2.8).* The result is a fair semicircuit computing f | (I_1 = b): its linear system is equivalent to the old one, with I_1's output additionally fixed to b. □ In `elimination.md`, "constant substitution to a gate I" is this operation, applied when I depends on an unprotected variable x_i; the restriction makes x_i a linear variable, so it is one substitution in the count t.

*Model transfer.* Every acyclic B2 circuit is a fair semicircuit (cyclic part empty) with the same gate count. Elimination is carried out for fair semicircuits — "circuit" means fair semicircuit, as in the source paper — and every bound proved for them applies in particular to the acyclic circuits counted by C_B2.

## 4. Gates, troubled gates, packing, potential

A *k-gate* (resp. *k-variable*) has out-degree exactly k; k+ means out-degree at least k; a *1+-variable* has out-degree ≥ 1. By functionality (source §2.1), the 16 kinds of B2 gates fall into four classes: *trivial* (constant function: f(x, y) = c); *degenerate* (depends on exactly one input: f(x, y) = (x ∧ c_1) ⊕ (y ∧ ¬c_1) ⊕ c_2 — outputs one input or its negation); *⊕-type* (affine, depends on both inputs: f(x, y) = x ⊕ y ⊕ c); *∧-type* (quadratic, depends on both inputs: f(x, y) = ((x ⊕ c_1) ∧ (y ⊕ c_2)) ⊕ c_3). A gate that becomes trivial or degenerate after a substitution (e.g. by being fed a constant) is called trivialized or degenerate in the same sense; Rules 2 and 3 of §5 eliminate such gates when fed by a constant.

**Troubled gate:** an ∧-type 1-gate fed by two 2-variables. Two troubled gates are *adjacent* if they share an input variable. A *troubled pair* is a pair of adjacent troubled gates; the variables feeding both are its *inner* variables, those feeding exactly one are *boundary* variables; a pair with two inner variables and no boundary variable is *compact*. A *packing* P is a set of troubled pairs pairwise disjoint in gates; a pair in P is a *pack*; a troubled gate lying in no pack of P is *truly troubled*. **Potential:** Φ(C, P) = (number of truly troubled gates) + |P| = (number of troubled gates) − |P|. By Proposition 5.6 (§5; source Proposition 3.10), no normalization rule can destruct a pack while leaving both mates troubled, so Φ can grow only through directly introduced truly troubled gates or packs (a truly troubled gate arising from unpacking a pack (G, G′) would force G′ to become non-troubled, unless the inner variables no longer feed G′).

## 5. Measure, potential bound, normalization, pair elimination

**Measure.** For positive constants α_I, α_Q, α_φ fixed in §7: μ(C, P, R) = g + α_I·i + α_Q·q + α_φ·Φ(C, P), where g is the number of gates, i the number of influential inputs (1+-variables or protected variables), and q the number of quadratic equations of R.

**Lemma 5.1 (potential bound).** If C computes an affine disperser for dimension d, then for every packing P, Φ(C, P) ≤ n/2 + 5d/2. *Proof.* Φ(C, P) = #troubled − |P| ≤ #troubled, so bound t, the number of troubled gates; this is Lemma 3 of FGHK16 (Lemma 3.7 of the source paper), whose proof we reproduce. Call two inputs *neighbors* if they feed a common troubled gate. A variable feeding a troubled gate is a 2-variable, so if v_i counts the variables feeding exactly i troubled gates, v_i = 0 for i > 2. Double-counting wires from inputs to troubled gates: 2t = v_1 + 2v_2, and v_1 + v_2 ≤ n. Suppose t ≥ n/2 + 5d/2; then n + 5d ≤ 2t = v_1 + 2v_2 ≤ n + v_2, so v_2 ≥ 5d. Let T be the inputs feeding two troubled gates, |T| = v_2 ≥ 5d. Construct disjoint X ⊆ T and Y ⊆ V = {inputs}, with |X| = d and consistent linear equations on Y making C independent of the variables of X: repeat d times — pick any x ∈ T; x feeds two troubled gates, whose other feeders are y_1, y_2 (possibly y_1 = y_2); add x to X and y_1, y_2 to Y. Fixing y_1, y_2 to constants eliminates both troubled gates fed by x and leaves C independent of x (if y_1 = y_2, one of x = c, y_1 = c, or y_1 = x ⊕ c eliminates both gates for some constant c; in the x = c case, interchange the roles of x and y_1). Each of y_1, y_2 has at most one neighbor different from x; remove x, y_1, y_2 and the other neighbors of y_1, y_2 — at most five variables — from T. Since |T| ≥ 5d the routine runs d times, and since those neighbors are removed from T, the neighbors of every variable picked later lie outside the Y chosen earlier, so all the assignments stay mutually consistent. Now impose the |Y| equations and set every remaining variable in V ∖ (X ∪ Y) to a constant: C then evaluates to a constant, since it does not depend on X and everything else is fixed. The number of independent linear equations is |Y| + |V ∖ (X ∪ Y)| = n − d, i.e. f is constant on an affine subspace of dimension d, contradicting the disperser property. Hence t < n/2 + 5d/2, and Φ(C, P) ≤ n/2 + 5d/2. □

**Corollary 5.2 (size from measure).** If every circuit C computing f with every packing P satisfies μ(C, P, ∅) ≥ μ_0, then |C| ≥ μ_0 − α_I·n − α_φ·(n/2 + 5d/2); in particular if f is an affine disperser for sublinear dimension d = o(n) and μ_0 = δn − o(n), then |C| ≥ (δ − α_I − α_φ/2)·n − o(n).

**Normalization.** After substitutions, redundant gates are removed by the five rules of §3.3 of the source paper. A circuit to which no rule applies is *normalized*. β = the measure decrement of one normalization step. *Rule 1:* G a 0-gate fed by I_1, I_2; remove G, out-degrees of I_1, I_2 drop by one. A new troubled gate can only be I_1, I_2, or a gate fed by one of them; troubled gates are fed by 2-variables, so I_1 (resp. I_2) feeds at most two troubled gates, and two fed by the same input are packed together; hence each input of G produces at most a troubled gate or a pack: ΔΦ ≤ 2, β ≥ 1 − 2α_φ. *Rule 2:* G a trivialized gate fed by a constant b; remove G, its descendants are now fed by a constant c. If the other input I_1 is not a constant, its out-degree drops by one and a new troubled gate is either I_1 or a gate fed by I_1: ΔΦ ≤ 1, β ≥ 1 − α_φ. *Rule 3:* G a degenerate gate fed by a constant b and another node I_1; remove G and rewire so G's descendants are fed directly by I_1 (if G computes the negation of I_1, the descendants' functions are modified accordingly; if G is the output gate, the elimination is the one governed by (T2) below): ΔΦ ≤ 1, β ≥ 1 − α_φ. *Rule 4:* G is *useless* if it is a 1-gate with only descendant Q, such that another input of Q also feeds G; eliminate G and rewire the nodes feeding it to feed Q. If one input I_2 of G is a gate, it is the only gate that may become troubled; if I_2 is a variable, new troubled gates must be fed by it, forming a troubled gate or a pack: ΔΦ ≤ 1, β ≥ 1 − α_φ. *Rule 5:* G a gate whose two inputs coincide; make it degenerate (or trivial) by changing one (or both) of its inputs to a constant, and eliminate it via Rule 2 or Rule 3: ΔΦ ≤ 1, β ≥ 1 − α_φ.

**Proposition 5.6 (pack destruction; Proposition 3.10 of the source paper).** If a pack (G_1, G_2) is destructed after applying a normalization rule, at least one of G_1, G_2 becomes non-troubled; hence potential increases can only be caused by direct introduction of troubled gates or packs. *Proof.* Suppose a rule destructs (G_1, G_2) while both remain troubled. Destructing the pack requires the descendants of the pair's inner variable to change; but the normalization rules change only the descendants of the gates feeding the eliminated gate. Hence one of G_1, G_2 is itself the eliminated gate — a contradiction. □

**Lemma 5.3 (normalization lemma; Lemma 3.11 of the source paper).** Let β be the measure decrement of one normalization step. For Rule 1, β ≥ 1 − 2α_φ; for Rules 2–5, β ≥ 1 − α_φ. Moreover, if the eliminated gate G is fed by I_1, I_2, every potential increase in the step is a directly introduced troubled gate or pack (Proposition 5.6 rules out increases by pack destruction), and each introduced troubled gate or pack is I_1 or I_2 or is fed by one of them after the step, as verified rule by rule above. As α_φ < 1/2, β ≥ 0: normalization never hurts, so C may be assumed normalized whenever a case of §6 begins.

**Terminal-output convention.** The case analysis of §6 repeatedly eliminates gates at or near the output of the circuit, so we record once how the output gate is treated; the source paper uses the same convention (its Rule 3 states the exception in (T2), and its §4.1 preamble and Case 6.1 argue as in (T1) and (T3)).

(T1) *The output gate is never trivialized.* Throughout the elimination, the current circuit computes f | R on a source R of dimension at least 2d — the iteration's steps start above 2d, and a multi-substitution leaf can land exactly on 2d at its last substitution — and by Proposition 3.1 the function f | R is not constant (a dimension-2d rdq-source already contains a d-dimensional affine subspace on which the disperser is nonconstant). Hence the output gate never computes a constant function, and no operation in Cases 0–8 trivializes it: whenever a case trivializes a gate, that gate is not the output gate. (This also legitimizes Case 0.3: the output gate never has a constant truth table or two constant inputs, since either would make f | R constant.)

(T2) *A degenerate output gate is eliminated by Rule 3, with one corner excluded.* If the output gate becomes degenerate with remaining input I_1, Rule 3 applies under the exception stated with it, in the following forms. If I_1 is a gate, I_1 becomes the output gate with its own function adjusted for a possible negation — complementing a gate's function stays within the basis, and the descendants of I_1 are unaffected. If I_1 is a primary input x_j and the output gate passes x_j unchanged, x_j becomes the output, exactly. The remaining form — I_1 a primary input and the output gate computing ¬x_j — never occurs at a stage of the elimination: it would make the current circuit compute f | R = ¬x_j on a source R of dimension at least 2d; fixing the free variable x_j = 0 linearizes the at most one quadratic equation containing x_j (the equations of an rdq-source are disjoint), and linearizing the remaining equations exactly as in the proof of Proposition 3.1 fixes at most d variables in total, leaving an affine subspace of dimension at least d on which f is constant, contradicting the disperser property. A primary input has no truth table to adjust, and selecting ¬x_j as the output is not an operation of the model; it is the exclusion, not an adjustment, that licenses this convention. In every form that occurs, the elimination counts as one eliminated gate with ΔΦ ≤ 1, exactly as the case ledgers record for any Rule 3 elimination.

(T3) *Rule 1 never removes the output gate.* Rule 1 is applied in `elimination.md` only to a gate whose out-degree has dropped to 0 during the operation; the output gate has out-degree 0 from the start and is never removed by Rule 1. In every invocation in `elimination.md`, the gate removed by Rule 1 has a descendant named in the configuration at hand — hence is not the output gate — or does not reach the output gate (Case 0.2). Examples: in Case 5.4.2.2 the gate B removed by Rule 1 has the descendant E; in Case 5.4.1.4, step (2), the gate G removed by Rule 1 feeds D; in Case 7 the gate Q eliminated by Rule 1 feeds the minimal gate G. Case 7 is in any event the dedicated treatment of the output gate when its inputs are a 2-variable and a 1-gate.

**Lemma 5.4 (troubled pair elimination).** Let C compute f | R with a non-empty packing P, and suppose f is not trivial under two substitutions. Then substitutions can be performed whose measure decrement per substitution is at least δ_p = min{ (3/2)α_I, 3 − 3α_φ + α_I + α_Q }. *Proof.* Take a pack (G_1, G_2), inner variable y, boundary variables x, z. If the pack is not compact: if y is protected, substitute a constant to y trivializing G_1; G_2 and the only descendants of G_1 become degenerate; normalization removes three gates with ΔΦ ≤ 3, and the substitution kills one quadratic equation and one influential input, so Δμ ≥ 3 − 3α_φ + α_I + α_Q. If y is unprotected, substitute constants to x and z trivializing both gates; all three variables become non-influential, so Δμ ≥ 3α_I over two substitutions, i.e. (3/2)α_I each. If the pack is compact, fed by x and y: if one is protected, argue as in the first subcase; otherwise both are unprotected, and enumerating the two gates' functions shows a constant substitution to x or y, or the affine substitution y ← x ⊕ b, trivializes both, giving Δμ ≥ 2α_I > (3/2)α_I for one substitution. □ At the parameters of §7, δ_p = min{72/5, 69/5} = 69/5 > δ = 64/5, so pair elimination is never the bottleneck.

## 6. The elimination theorem

**Theorem 6.1.** Let C compute an affine disperser f for dimension d, let P be any packing, and set

    δ = α_I + min{ α_I/3, 2 − 2α_φ + α_Q, 4 − 4α_φ, 3 + α_φ, 5 − α_Q,
                   (5 − 2α_φ + α_Q)/2 }.

If the dimension of the ambient rdq-source is > 2d + 2, then either f becomes trivial after at most two substitutions (impossible while the dimension condition holds, by Proposition 3.1), or there is a sequence of 1 ≤ t ≤ 3 substitutions after which μ drops by at least δ per substitution: μ(C′, P′, R′) ≤ μ(C, P, R) − δ·t for the resulting circuit, source and a suitable packing P′.

The proof is the case analysis in `elimination.md` (Cases 0–8), ordered: each case assumes all earlier configurations absent. Case 0: unnormalized circuits, constant-output gates, non-empty packings (via Lemma 5.4). Cases 1–4: protected variables with high out-degree, protected 0-variables, high-out-degree ∧-type configurations. Case 5 and subcases: the troubled configuration itself (an ∧-type gate fed by two 2-variables), split by whether and where new troubled gates appear. Case 6: substituting a constant to a protected variable, with the configurations left by its exceptions. Case 7: the output gate when its inputs are a 2-variable and a 1-gate. Case 8: both inputs of the minimal ∧-type gate are gates — via the xor-reconstruction of §3 — or a 2-gate remains; the deepest subcases also use quadratic substitutions. Every leaf records t, the gates eliminated, and the bounds on ΔΦ, Δi, Δq, hence the per-substitution decrement; the minimum over all leaves gives exactly the six quantities in δ (the term 2 − 2α_φ + α_Q is attained by Cases 1 and 6.1 and is not the minimum at the chosen parameters).

Iterating Theorem 6.1 from (C, P, ∅) while the dimension exceeds 2d + 2: the dimension starts at n and each substitution lowers it by at most 1, so the process runs for at least n − 2d − 2 substitution steps before it can stop, and f never becomes trivial on a source of dimension > 2d (Proposition 3.1). Hence

    μ(C, P, ∅) ≥ δ·(n − 2d − 2).

## 7. Parameters and the final assembly

Maximize δ − α_I − α_φ/2 over α_I, α_Q, α_φ ≥ 0 with α_φ ≤ 1/2, where δ is as in Theorem 6.1. Writing u = δ − α_I, the constraints are u ≤ α_I/3, u ≤ 2 − 2α_φ + α_Q, u ≤ 4 − 4α_φ, u ≤ 3 + α_φ, u ≤ 5 − α_Q, u ≤ (5 − 2α_φ + α_Q)/2. Take

    α_I = 48/5,   α_Q = 9/5,   α_φ = 1/5.

Then the six quantities are 16/5, 17/5, 16/5, 16/5, 16/5, 16/5, so u = 16/5 and δ = 48/5 + 16/5 = 64/5. Optimality: for any feasible parameters, take the nonnegative combination of the three constraints u ≤ 3 + α_φ, u ≤ 5 − α_Q and u ≤ (5 − 2α_φ + α_Q)/2 with weights 7/10, 1/10 and 1/5 respectively:

    (7/10)(3 + α_φ) + (1/10)(5 − α_Q) + (1/5)·((5 − 2α_φ + α_Q)/2)
      = 21/10 + 1/2 + 1/2 + (7/10 − 1/5)α_φ + (−1/10 + 1/10)α_Q
      = 31/10 + (1/2)α_φ.

Hence u ≤ 31/10 + α_φ/2 for every feasible choice, so the objective u − α_φ/2 is at most 31/10, and the displayed parameters attain it. Finally, Corollary 5.2 with μ_0 = δ(n − 2d − 2), d = d(n) = o(n):

    |C| ≥ δ(n − 2d − 2) − α_I·n − α_φ·(n/2 + 5d/2)
        = (64/5 − 48/5 − 1/10)·n − o(n)
        = (32/10 − 1/10)·n − o(n)
        = (31/10)·n − o(n).

Since C was an arbitrary B2 circuit computing f_n exactly, C_B2(f_n) ≥ (31/10)n − o(n) for all large n, and therefore liminf C_B2(f_n)/n ≥ 31/10. ∎

Evidence: `disperser.md` (BK Sections 3, 6, 7 in full); `elimination.md` (Cases 0–8 in full).

**Supplementary exact finite checks.** `family.py` implements the §1 recipe (integer arithmetic only) for small lengths, under the stand-in rule stated in that file — the numerically smallest monic irreducible in place of Shoup's output; by §1 the choice is immaterial to every audited quantity. `audit.py` runs exact finite checks of the decidable arithmetic this proof relies on: the parameter recipe and its §2 consequences (least prime k and its minimality, the Bertrand consequence k < 4n^{3/5}, the chain link r(r+1)/t < k, the threshold d(n) and its vacuity range n ≤ 7776), the rational arithmetic of §7 (the six quantities of Theorem 6.1, the LP optimality certificate with weights 7/10, 1/10, 1/5, the final assembly (31/10)), Lemma 5.4's δ_p = 69/5, the Case 5.2 ledgers of the branch split in `elimination.md` (Δμ ≥ 72/5 at t = 1 when D is trivialized; over t = 2: Δμ ≥ 128/5 = 2δ in branch 5.2.β0, Δμ ≥ 129/5 in branch 5.2.γ and the remaining ⊕-type branches 5.2.β1–β3, Δμ ≥ 171/5 in the residual branch 5.2.β4, Δμ ≥ 172/5 in the A = C configuration of 5.2.γ, and Δμ ≥ 169/5, 173/5 and 178/5 in the ∧-type branch), the Case 5.4.1.4 coincidence ledger (Δμ ≥ 129/5 over t = 2), the Case 8.2.4.3.1 constant-R ledger (Δμ ≥ 64/5 = δ at t = 1), and the evaluator's agreement with the subset definition of Sym_r^t at small lengths. No finite run proves the asymptotic theorem or the disperser estimate; those rest on the prose above and in `disperser.md` alone. The audit certifies only that the finite arithmetic stated in the proof is correct as stated.
