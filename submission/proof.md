# (31/10, 1) for the BK/Shoup family

**Theorem.** For the total family (f_n) of §1, in the Li-Yang B2 model (acyclic circuits: one output, two pins/gate, all 16 truth tables, unrestricted depth/fan-out, exact, nonuniform), liminf C_B2(f_n)/n >= 31/10, i.e. C_B2(f_n) >= (31/10)n - o(n) for all large n. Proof: Li-Yang gate elimination (ECCC TR21-023 Thm 1.1); cases in `elimination.md`. The disperser premise is proved: §2 reduces it to BK Thm 2.2, hypotheses checked; `disperser.md` gives BK §§3,6,7. Unproved only: Bertrand (once, §2), basic field facts; Prop 3.3 is proved in `disperser.md`. Parameters aI = 48/5, aQ = 9/5, ap = 1/5, δ = 64/5 (aI,aQ,ap := α_I,α_Q,α_φ; DPhi,Dmu,Di,Dq := ΔΦ,Δμ,Δi,Δq). Operators: + = xor, & = and, ! = negation.

## 1. Family

f_0({})=0. For n>=1: k = least prime > 2n^{3/5}; r = ⌈n/k⌉; t = ⌈√r⌉. Fix Shoup's deterministic poly-time F_2 irreducible-polynomial procedure (choices fixed by it), giving monic irreducible q_k in F_2[Z], deg k; K = F_2[Z]/(q_k), basis 1,Z,..,Z^{k-1}: one language, no per-length advice. Pad x in {0,1}^n by rk-n zeros, split into r k-bit blocks; block i is φ_i(x) in K, first bit the constant coefficient; φ injective F_2-linear. π:K->F_2 = constant coefficient, F_2-linear, π(1)=1. y_i = φ_i(x)^{2^{i+1}-1}. f_n(x) = π(Sym_r^t(y_1,..,y_r)), Sym_r^t the degree-t elementary symmetric polynomial (sum over t-subsets).

## 2. Uniformity and disperser

Evaluator A: trial-division scan for k (poly, k = O(n^{3/5}) by Bertrand); fixed Shoup procedure (poly in k); Sym by E_0=1, E_j <- E_j + y_i*E_{j-1}, j = min(i,t) down to 1 (subset split on i; no char-2 division); O(r*t) field ops; output π(E_t), so (f_n) in P. V(x,w)=1 iff w empty and A(x)=1, W=0 gives NP membership both ways (witnesses are not circuit inputs).

§1 is BK Thm 2.2, eq. (4). Checks: (i) p=2; m := n (true length); block length least prime > 2m^{3/5} = k; r = ⌈m/k⌉; t = ⌈√r⌉. (ii) Thm 2.2 allows any injective F_2-linear φ: F_2^m -> (F_{2^k})^r and any nontrivial F_2-linear π. Ours (padding, block split, fixed basis identifications) is a coordinate embedding: F_2-linear, injective (n <= rk). No divisibility n/k assumed: the short last block is absorbed by φ, whose last component has top rk-n coordinates constantly 0. π is F_2-linear, π(1)=1, nontrivial; Prop 3.3 (`disperser.md`): π(x)=Tr(µx), unique µ!=0. (iii) Eq. (4)'s exponent for block i, 1+p+..+p^i, equals 2^{i+1}-1 at p=2, matching y_i (BK §7 reverses indexing; renamed in `disperser.md`). Hence f_n is nonconstant on every affine subspace of dim > 6n^{4/5}, all large n. The deduction Thm 7.3 => Thm 2.2 (in `disperser.md`) estimates jointly in r ~ n/k: Bertrand => k <= 4n^{3/5}, whence r(r+1)/t <= (r+1)^{3/2} <= (n^{2/5}/2 + 2)^{3/2} < 2n^{3/5} < k, and nr/t + nt + r(r+1) ≲ 4.25*n^{4/5} < 6n^{4/5}. For n <= 7776 = 6^5, 6n^{4/5} >= n: vacuous. Put d(n) = ⌊6n^{4/5}⌋+1; eventually 1 <= d(n) <= n, d(n)/n -> 0. f is an *affine disperser for d* if nonconstant on every d-dim affine subspace; f_n is one for d(n) = o(n).

## 3. rdq-sources, substitutions

rdq-source R: X = {x_1,..,x_n} splits disjointly into free F, linear L, quadratic Q. For x_j in Q, one equation x_j = ((x_i+c_1)&(x_k+c_2))+c_3, x_i,x_k free, RHS vars of distinct quadratic eqs disjoint. For x_j in L, x_j = (+_{x_i in F_j} x_i)+(+_{x_k in Q_j} x_k)+c_1, F_j⊆F, Q_j⊆Q. Dimension = |F| (2^{|F|} inputs). A free variable on a quadratic RHS is *protected*; the two in one equation are *coupled*; other free variables *unprotected*. C *computes f | R* if L,Q vars are 0-vars in C (§4) and C agrees with f on inputs satisfying R.

Prop 3.1. A dim-2d rdq-source contains a d-dim affine subspace; hence a disperser for d is nonconstant on any rdq-source of dim > 2d. *Proof.* 2d free vars, <= d quadratic equations; for each affine equation using both vars of one quadratic equation, fix one to 0: all equations become affine, dim >= d remains. □

Substitutions restrict R (a *g-substitution* sets x_j to g of the others); §5 deletes redundant gates. (1) Constant: x_j<-c, x_j free. (2) Affine to couples: x_j<-x_k+c, x_j,x_k coupled protected; their quadratic equation becomes affine, its variable moves Q->L, x_j linear. (3) General affine: x_j<-x_k+c, x_j unprotected free, x_k free; valid because unprotected x_j occurs on no quadratic RHS - unfold x_j = x_k+c in all affine equations, x_j becomes linear. (4) Quadratic: x_j<-((x_i+c_1)&(x_k+c_2))+c_3, x_j,x_i,x_k unprotected free, rewiring makes x_j a 0-var; one quadratic eq. added (RHS = existing feeders named in the case), x_j quadratic. Each type kills one free variable; by Prop 3.1 a disperser for d stays nontrivial after n-2d substitutions. After rewiring, nonzero-out-degree vars are free; subst. constants (in-degree 0) are removed by §5. A *k-variable* has out-degree k; an input is *influential* if 1+-variable or protected.

Semicircuits (for Case 8). A *cyclic xor-circuit* is a directed graph of xor-type gates and inputs, cycles allowed, many outputs. A *semicircuit* is C_1 composed with an acyclic C_2 fed by it. C_1's gates G_1,..,G_m (f_i in {+,=} of two inputs) form a linear system Ag = b over F_2 (inputs as constants); *fair* iff A has full rank; a semicircuit is fair iff C_1 is.

Prop 2.6. If C is fair and no cycle is created, any g-substitution to x_j is valid, giving R′ = R restricted by x_j = g(*) and a fair C′ computing f | R′ (rewriting R's equations and x_j's occurrences preserves the linear system's rank).

**Prop 2.7 (path).** If gate I in a fair cyclic xor-circuit computes an affine function depending on x_i, a path x_i -> I exists. *Proof.* Gates unreachable from x_i form a fair semicircuit omitting x_i as input, so their outputs are determined without x_i; I depends on x_i. □

Xor-reconstruction / Prop 2.8. Let G be a topologically minimal &-gate fed by I_1,I_2, I_1 a gate; I_1's subcircuit is &-free, a cyclic xor-circuit. To set I_1 <- b (I_1 depending on x_i): the path x_i = I_0->I_1->...->I_k = I_1 has I_j = f_j(I_{j-1},T_j), f_j in {+,=}, invertible in argument 1; rewrite I_{j-1} = g_j(I_j,T_j), replace I_k by b and x_i by a new xor-gate Z, rewire along it. Result: a fair semicircuit computing f | (I_1 = b) - the old linear system with I_1's output fixed to b. This is `elimination.md`'s "constant substitution to gate I"; x_i becomes linear, one subst. in t.

Transfer. Every acyclic B2 circuit is a fair semicircuit (empty cyclic part), same gate count; "circuit" = fair semicircuit below, and bounds transfer to C_B2.

## 4. Gates and potential

A *k-gate*/*k-variable* has out-degree k (k+: >=k). Gate kinds by function: *trivial* f(x,y)=c; *degenerate*, f(x,y)=(x&c_1)+(y&!c_1)+c_2 (one input or its negation); *xor-type*, f(x,y)=x+y+c; *&-type*, f(x,y)=((x+c_1)&(y+c_2))+c_3. "Trivialized/degenerate": becoming so after substitution; Rules 2-3 remove such gates when constant-fed.

A *troubled gate* is an &-type 1-gate fed by two 2-variables; TGs sharing an input are *adjacent*. A *troubled pair* is an adjacent TG pair; variables feeding both are *inner*, feeding exactly one *boundary*; a pair with two inner and no boundary variables is *compact*. A *packing* P is a set of TG pairs disjoint in gates (members: *packs*); a TG in no pack is *truly troubled*. Φ(C,P) = (# truly troubled) + |P| = (#TG) - |P|.

## 5. Measure and normalization

μ(C,P,R) = g + aI*i + aQ*q + ap*Φ(C,P): g gates, i influential inputs, q quadratic eqs of R.

Lemma 5.1. If C computes a disperser for d, Φ(C,P) <= n/2 + 5d/2 for all P. *Proof* (FGHK Lemma 3). Φ <= #TG =: t (TG = troubled gate, §4). Inputs feeding a common TG are *neighbors*; they are 2-variables. Let v_i = # feeding exactly i TGs (v_i = 0, i > 2). Double-counting wires: 2t = v_1 + 2v_2, v_1 + v_2 <= n. If t >= n/2 + 5d/2 then n + 5d <= 2t <= n + v_2, so v_2 >= 5d. Put T = inputs feeding two TGs, |T| = v_2 >= 5d. Build disjoint X⊆T, Y, |X| = d, with consistent linear eqs on Y making C independent of X: d times pick x in T, its gates' other feeders y_1,y_2 (possibly equal); x->X, y_1,y_2->Y. Fixing y_1,y_2 kills both TGs fed by x, so C ignores x (if y_1 = y_2, one of x = c, y_1 = c, y_1 = x+c kills both; if x = c, swap x,y_1). Each y_i has <=1 neighbor != x; delete x,y_1,y_2 and their other neighbors (<=5 vars) from T: |T| >= 5d allows d rounds, deletions keeping later picks consistent. Impose the |Y| equations, set vars outside X∪Y constant: C is constant under |Y| + (n - |X| - |Y|) = n - d equations - f constant on a d-dim affine subspace, contradiction. Hence t < n/2 + 5d/2. □

Cor 5.2. If μ(C,P,{}) >= μ_0 for all C,P, then |C| >= μ_0 - aI*n - ap*(n/2 + 5d/2); for d = o(n), μ_0 = δn - o(n): |C| >= (δ - aI - ap/2)*n - o(n).

Normalization (source §3.3); β = a step's μ-decrement; *normalized*: no rule applies.
R1: 0-gate G fed by I_1,I_2: remove G; I_1,I_2 lose out-edges. New TGs lie only at I_1,I_2 or gates they feed; <=2 per input, shared ones packed: DPhi <= 2, β >= 1 - 2ap.
R2: trivialized G fed by constant b: remove G; descendants get a constant; other input I_1 (if not constant) loses an out-edge; new TGs: I_1 or fed by it: DPhi <= 1, β >= 1 - ap.
R3: degenerate G fed by constant b and I_1: remove G, descendants fed directly by I_1 (adjust for negation; if G is the output, I_1 becomes it): DPhi <= 1, β >= 1 - ap.
R4: useless G (1-gate, only descendant Q, another input of Q also feeds G): remove G, its feeders feed Q. Only an input gate can become troubled; variable-fed new TGs share the input (gate or pack): DPhi <= 1, β >= 1 - ap.
R5: G with equal inputs: make one/both inputs constant so G is degenerate/trivial, then R2/R3: DPhi <= 1, β >= 1 - ap.

Prop 5.6 (source 3.10). If a rule destructs a pack (G_1,G_2), some G_i becomes non-troubled; hence Φ grows only by directly introduced gates/packs. *Idea:* destructing needs the inner variable's descendants to change; rules change only descendants of gates feeding the eliminated gate, so some G_i would itself be eliminated - contradiction.

Lemma 5.3. Per step, R1: β >= 1 - 2ap; R2-R5: β >= 1 - ap. Any Φ increase is a directly introduced TG/pack (Prop 5.6), each being I_1, I_2, or fed by one (per the rules). Since ap < 1/2, β >= 0: normalization never hurts; assume C normalized at each §6 case.

Terminal output. (T1) *Output never trivialized:* throughout, C computes f | R, dim R > 2d, and Prop 3.1 makes f | R nonconstant, so the output never computes a constant; no Case 0-8 operation trivializes it - any trivialized gate is not the output (legitimizing Case 0.3). (T2) *Degenerate output:* Rule 3 applies by its stated exception - I_1 becomes the output, functions adjusted for negation, counting one eliminated gate with DPhi <= 1. (T3) *Rule 1 never removes the output:* the output has out-degree 0 from the start, while R1 applies only to gates whose out-degree dropped to 0 during it; in every `elimination.md` invocation the removed gate has a descendant named in its configuration (hence is not the output) or does not reach the output (Case 0.2) - Case 5.4.2.2 (B; descendant E), Case 5.4.1.4 step (2) (G; feeds D), Case 7 (Q; feeds the minimal gate G).

Lemma 5.4 (pairs). If P != {} and f is nontrivial under two subst., some subst. give per-subst. Dmu >= δ_p = min{(3/2)aI, 3 - 3ap + aI + aQ}. *Sketch.* Pack (G_1,G_2), inner y, boundary x,z. Non-compact, y protected: constant to y trivializes G_1; G_2 and G_1's descendants degenerate; 3 gates removed, DPhi <= 3, one quadratic eq. and one influential input die: Dmu >= 3 - 3ap + aI + aQ. Non-compact, y unprotected: constants to x,z trivialize both; all three vars non-influential: Dmu >= 3aI over t=2, i.e. (3/2)aI each. Compact: a protected x/y reduces to the first case; else a constant to x or y, or y<-x+b, trivializes both: Dmu >= 2aI > (3/2)aI in one subst. □ At §7 values δ_p = min{72/5, 69/5} = 69/5 > δ = 64/5: pairs never bottleneck.

## 6. Elimination

Thm 6.1. For C computing a disperser for d, any P, let δ = aI + min{aI/3, 2 - 2ap + aQ, 4 - 4ap, 3 + ap, 5 - aQ, (5 - 2ap + aQ)/2}. If dim R > 2d + 2, either f becomes trivial after <= 2 subst. (impossible, Prop 3.1) or some 1 <= t <= 3 subst. give μ(C′,P′,R′) <= μ(C,P,R) - δ*t.

Proof: `elimination.md`, Cases 0-8 (each assumes earlier configs absent): Case 0 - unnormalized circuits, constant-output gates, nonempty packings (Lemma 5.4); Cases 1-4 - protected vars of high out-degree, protected 0-vars, high-out-degree &-configurations; Case 5 - the TG configuration (&-gate on two 2-vars), split by whether/where new TGs arise; Case 6 - constant subst. to a protected var, with exceptions; Case 7 - output gate fed by a 2-var and a 1-gate; Case 8 - both inputs of the minimal &-gate are gates (§3 reconstruction), or a 2-gate remains; deepest subcases use quadratic subst. Each leaf records t, gates eliminated, DPhi, Di, Dq, its decrement; the leaf minimum is δ's six quantities (attained at Cases 1, 6.1).

Iterating while dim > 2d + 2 (dim drops <=1 per subst., f never trivial) gives >= n - 2d - 2 steps, so μ(C,P,{}) >= δ(n - 2d - 2).

## 7. Parameters

Maximize δ - aI - ap/2 (aI,aQ,ap >= 0, ap <= 1/2). With u = δ - aI: u <= aI/3, 2 - 2ap + aQ, 4 - 4ap, 3 + ap, 5 - aQ, (5 - 2ap + aQ)/2. Take aI = 48/5, aQ = 9/5, ap = 1/5. The six quantities are 16/5, 17/5, 16/5, 16/5, 16/5, 16/5; hence u = 16/5 and δ = 48/5 + 16/5 = 64/5. Optimality certificate: weights 7/10, 1/10, 1/5 on the rows u <= 3 + ap, u <= 5 - aQ, u <= (5 - 2ap + aQ)/2 give

  (7/10)(3 + ap) + (1/10)(5 - aQ) + (1/5)*((5 - 2ap + aQ)/2)
    = 21/10 + 1/2 + 1/2 + (7/10 - 1/5)ap + (-1/10 + 1/10)aQ = 31/10 + (1/2)ap,

so u <= 31/10 + ap/2 always, and u - ap/2 <= 31/10 is attained. Cor 5.2 with μ_0 = δ(n - 2d - 2), d = d(n):

  |C| >= δ(n - 2d - 2) - aI*n - ap*(n/2 + 5d/2)
    = (64/5 - 48/5 - 1/10)*n - o(n) = (32/10 - 1/10)*n - o(n) = (31/10)*n - o(n).

Hence C_B2(f_n) >= (31/10)n - o(n) for all large n, i.e. liminf C_B2(f_n)/n >= 31/10. ∎

