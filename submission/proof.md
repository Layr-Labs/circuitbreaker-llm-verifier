# A full-B2 lower bound of 9413355971276602268255192551/1837861338540586116300000000 for a fixed norm-one predicate

The graph law has cap 12252408923603907442/10^20. Two exact replay phases check the graph and the entropy inequalities.

## A. Fixed family and eligibility

Use the norm-graph construction of Kollár, Rónyai and Szabó (Combinatorica 16, 399–406, 1996, DOI 10.1007/bf01261323) and the norm-one language of earlier entries.

### A.1 Uniform evaluator

For actual length N<64 set f_N=0. Otherwise put
```
L=floor(log2 N), B=ceil(log2(L+1)), a=L−B, p=2^a,
ell=floor(N/p), m=p*ell, d=2^floor(a/2), k=m/d, Q=2^k.
```
Here a≥3 and d divides m. The bounds 2^B≤ell<2^(B+1)<4(L+1), N/[4(L+1)]≤p≤N/(L+1), and N−p<m≤N give m/N→1 and k,d log2 d=O(sqrt(N log N))=o(N).

Use the first monic irreducible binary polynomial of degree ell in coefficient order. Test by division by all monic polynomials of degree at most floor(ell/2), in polynomial time since ell=O(log N). Existence follows from the field of roots of X^(2^ell)−X. Its proper subfields have dimensions at most ell/2 and their union has size at most Σ_{j≤ell/2}2^j<2^ell. Take the ordered power basis.

Double the field a times. In binary dimension r, the trace Σ_{j<r}z^(2^j) is a nonzero binary linear map, since its polynomial degree is less than 2^r. Take its first basis vector u of trace one. The image of z↦z²+z is the trace-zero space by telescoping and its two-element kernel. Adjoin β with β²+β+u=0. Append β times the old basis to that basis. At most five old multiplications suffice per new multiplication. Trace tests, basis scans and exponentiation are polynomial in N by a=O(log N).

In the resulting field F of dimension m interpret the first m bits as z and ignore the suffix. Accept iff z^e=1, with e=(2^m−1)/(Q−1)=Σ_{0≤j<d}Q^j. The exponent has at most m bits. This defines one total language in P, also in NP with the verifier that accepts precisely the empty witness on accepted inputs.

### A.2 Density and essential bits

The polynomial X^e−1 divides the split squarefree X^(2^m−1)−1. Thus the accepting set S has e points and density at least 2^(-k). The count e is odd, so no active bit is inessential, since flipping any such bit would pair the accepted words. Fixing the suffix to zero reduces a circuit on actual length N to one on all m essential bits without adding gates.

### A.3 Common neighbors

For a matrix a_ij over a field with d distinct entries down each column, the d equations ∏_j(z_j−a_ij)=1 have at most d^d solutions. Homogenize with L_ij=X_j−a_ijY_j and W=∏_jY_j. Work modulo ∏_jL_ij−W. In column j the forms Q_ij=∏_{r≠i}L_rj form a basis of degree d−1 by evaluation at (a_ij,1). Their monomial multiples span every larger degree.

Take a form of degree D≥d in each pair. Factor out Y_j successively in this quotient. For the terms without Y_j, expand the X_j^D coefficient in products of Q forms and monomials on the other d−1 columns. Their degrees there are D or D−1 if a Y was already removed. Each Q omits one row, so some row i occurs in all d−1 factors. Replace X_j^D by (L_ij+a_ijY_j)X_j^(D−1). The entire product of row i's L factors equals W in the quotient. Both terms acquire Y_j. Thus the original form reduces to W times one of degree D−1 in each pair.

Set every Y_j to one. Homogenize any affine polynomial to common degree D and iterate the reduction. The affine quotient is spanned by d^d monomials with all exponents below d. Evaluation onto any finite set of distinct solutions is surjective. At each point a product of normalized linear factors, one distinguishing it from each other point, evaluates to its delta function. This bounds the number of solutions by d^d.

### A.4 Sumsets

If A+B⊆S and |A|≥d, choose distinct u_i∈A. Every v∈B satisfies ∏_{j<d}(v^(Q^j)+u_i^(Q^j))=1 for all i. Frobenius is injective, so the common-neighbor lemma applies in F. Its solution tuple contains v, giving |B|≤d^d. Set K=d^d+1 so log2 K=o(N). On any partition of the active coordinates, each factor of an accepting rectangle maps injectively to its coordinate-supported field elements. One factor has size less than K by the sumset bound.

Write n=m for the active dimension and k=ceil(log2 K) in the rest of the proof. Then n/N→1. Use the uniform measure ν on S of density 2^(−D), where D=o(n). Thus H_ν(X)=n−D in bits. To pass from a full-cube probability bound to ν, multiply by at most 2^D.

Let c=12252408923603907442/10^20, L=1+1/(2c), d0=12323231/300000000. The graph certificate must supply a finite local matched-cubic-tree bit process, equivariant under automorphisms preserving its matching, fair at each vertex, whose expected crossing count per vertex (each edge counted once) is c_*<c. The two edge types have probabilities a_M,a_C with a_M/2+a_C<c. Finite shell/window parameters and deterministic activation rate are chosen using strict margins, before circuit size tends to infinity. This is the only numeric graph input. The following argument turns it into a weighted-cut law without importing a nontrivial perfect-matching theorem.

## 1. Finite-local matching and edgewise transfer

On the cubic tree start with no matched edges. Fix J=2l+1. At each round mark vertices independently fair low or high. Augment all length≤J alternating paths with unmatched endpoints that are whole low components. The paths are disjoint. Send half a unit from every vertex changed by an augmentation to each endpoint. Invariant expected mass sent and received agree because an automorphism exchanges any two vertices. If u_r is the unmatched density, the expected changed fraction is at most (J+1)(u_r−u_{r+1})/2. Summing shows that every vertex stabilizes almost surely. A limit augmenting path of length≤J would persist in late rounds, each with fresh selection probability at least 2^(−4(J+1)), which has probability zero.

In the limit let U be the unmatched set of density u, X_i the vertices reached from U by even alternating paths of length≤2i, and Y_i=N(X_i). For i≤l every Y_i vertex is matched to a vertex of X_(i+1) outside U. A neighbor on a witnessing path is its matched predecessor. An outside neighbor reached on an unmatched edge either completes a forbidden augmentation or extends the even path via its matched partner. Cubic incidence counting gives density(Y_i)≥density(X_i), and invariant transport along partners gives density(X_(i+1))≥density(Y_i)+u. Thus (l+2)u≤1. Choose l and then finitely many rounds to make u_r as small as needed.

Run the matched rule on fully matched balls of the needed radius, including partners, and use fresh fair bits on bad balls. Good balls and joined adjacent good balls extend to the perfectly matched tree by completing each unused branch. Perfect matchings on it are locally isomorphic by branching. Give each pair fresh independent uniform anchor labels using 1−(1−min(U,V))² on independent endpoint uniforms. Hence both good and bad balls have fair marginals. Each edge is matched with probability (1−u_r)/3 and has bad-pair probability O(u_r). The crossing probability is at most (a_M+2a_C)/3+O(u_r)≤2c/3 for small enough u_r. Each direction costs at most c/3 by symmetry and fairness.

Complete each missing stub of a finite subcubic multigraph with a separate K4 whose edge is subdivided at the attachment. A loop counts twice in degree. The completion has the cubic tree as universal cover, described by reduced incidence walks. Added gadgets give no shortcut between original vertices.

Begin with independent tree labels. For each finite-graph vertex fiber, condition on untreated labels outside the fiber and on previous replacement uniforms. Replace each covered row by its Bernoulli conditional quantile with one fresh common uniform for that fiber. This retains conditional marginals and uses only the projected original local support. Deck transformations preserve the conditional law, so after all fibers are treated equivalent rows coincide. The law descends to the finite graph.

For a submodular payoff on finitely many rows, order rows by decreasing conditional probabilities. If S_i is the first i rows, diminishing returns gives f(A)≥f(∅)+Σ_{i∈A}(f(S_i)−f(S_(i−1))). Its expected lower bound uses only marginals. The common quantile produces prefix sets and attains equality. This preserves upper directional star and arm costs. For a bounded sum of local payoffs on O(n) vertices, only O(n) pairs of supports overlap, so variance is O(n). Chebyshev and a finite union bound give simultaneous o(n) errors.

## 2. Bounded bags in the circuit graph

Include a dummy output leaf. Before binary copy expansion the connected graph has n+t+1 vertices and 2t+1 wires, counting repeated input pins separately. Its cycle rank is r=t−n+1. Copy expansion adds equally many vertices and edges. Contracting one nonloop edge also preserves rank, leaves parallel edges as loops, and reduces the sum of endpoint degrees by two. Each loop counts twice.

## 3. Raw-bundle exclusions and physical accounting

Prune dead gates. Substitute constant and signed projection gates into consumers until every gate depends on two distinct actual inputs. The remaining affine functions are XOR or XNOR. The remaining nonlinear functions are signed AND. An output projection or negation of an internal producer makes that producer the output, with its truth table complemented when needed. Constants and signed raw projections cannot use all essential bits. Every native gate reaches the dummy output, which has no gate charge.

Let u count degree-one raws. A gate that receives one is a private recipient. There are o(n) occurrences of each of these patterns:

- (i) Two private raws at one gate.
- (ii) A degree-two raw feeds two private recipients with distinct private raws.
- (iii) Two degree-two raws have distinct private recipients and the same other consumer, which receives those two raws.

Within one type the raw bundles are disjoint by native degrees. A bundle of r raw bits feeds at most r−1 listed outputs and has no other exterior use. On the full cube, at least half its assignments have output fibers of size at least two, since singleton fibers number at most 2^(r−1). With Ω(n) bundles, independence and the binomial bound give exponentially small probability of fewer than 2k nontrivial fibers. This probability is o(1) under ν. But an accepted input with 2k such bundles gives an accepted product cell by fixing exterior raws and all bundle outputs. Split the nontrivial fibers into two groups to contradict §A. This proves the counts.

Use binary copy trees. A gate with d native uses, including the dummy output, has 2d−1 carrier edges. For a raw with d≥2 merge the leaf, first fork and one child. For d≥3 choose a child fork to retain every consumer arm. For d=2 choose one consumer as owner and retain the other raw arm. The raw then has 2d−3 carrier edges. Degree-one raws keep their leaves. Bags have degree at most three and own at most two raws. All vertices except the u+1 leaves are trivalent and rank stays r. Each owned raw is known in its bag from the start.

If a degree-two raw has exactly one private recipient, absorb it into its other consumer. Both-private exceptions number o(n) by type (ii). Absorb other degree-two raws into the earlier consumer. At a nonlinear gate designate one retained arm by these priorities: private input (later in fixed order if two), nonraw arm, then raw arm of native degree at least three. Omitted gates thus have two raw inputs of native degree two. Their pairs form a multigraph of maximum degree two and admit a greedy three-color edge coloring. In a color class the pairs are disjoint, so subadditivity bounds their sum of 2−H(U,V) in a mode by n−H(X|σ). Signing operands preserves the pair entropy, and I(U;V)≤2−H(U,V). Thus only o(n) total designation credit is lost in each high-entropy mode.

## 4. Pinning and high-entropy innovations

Use the master list of O(n) binary variables under ν. It contains each true producer bit with 2d+2 copies for native use count d and needed raw or terminal copies, and the four signed operand-pair indicators of each nonlinear gate. There are s≤sqrt(n) fixed list indices whose revealed word σ satisfies E_σ Σ_{i,j}|Cov(W_i,W_j|σ)|=O(n² n^(-1/4)). To prove this, reveal fresh independent uniform indices and values one at a time. The average conditional entropy of a fresh indexed bit starts at most one. Its expected drop is average conditional mutual information of two fresh indexed bits. Choose a step before sqrt(n) with drop O(n^(-1/2)), then a fixed tuple of indices with no larger average drop. In bits, relative entropy is at least 2 TV²/ln 2 by log-sum grouping of the positive difference set and the binary second derivative bound. Covariance is bounded by TV from the product law. Cauchy–Schwarz yields the stated covariance sum. Conditional Chebyshev gives concentration of any fixed finite list of bounded-coefficient linear combinations about their mode means to within ξn for fixed ξ>0, except o(1) mass.

Since H_ν(X)=n−D, the chain rule bounds each mean deficit by D+s: n−H(X|σ), Σ_i(1−H(X_i|σ)), and |U|−H(X_U|X_complement,σ). The last uses H(X_complement|σ)≤n−|U|. The same bound holds for R−H(N|σ) for R independent raw linear forms N, by completing them to an invertible chart. These deficits are nonnegative. Markov with a vanishing tolerance gives o(n) deficits outside o(1) mode mass. Use conditional entropies within a positive mode below.

Write m for the nonlinear gate count. On each current affine cell, fixing one operand makes that signed AND affine. Partition the cube chronologically by one operand per nonlinear gate and a final affine output query. Every cell has dimension at least n−m−1. Any accepted cell of dimension ≥2k violates the sumset bound by splitting its basis into two groups of k directions. Thus m≥n−o(n).

Express each signed nonlinear operand as a fixed raw linear form plus an affine function of earlier nonlinear outputs treated as symbols. Let j_g be the rank increment of the pair of raw forms for g beyond preceding forms. Fixing all R=Σ_g j_g independent forms and the final raw output form fixes the output, so R≥n−o(n) by the same cell obstruction. Let N_g be the new independent forms and P_g earlier basis values. Earlier nonlinear outputs depend on P_g. The mean sum of j_g−H(N_g|P_g,σ) equals R−H(N|σ)≤D+s.

Compare each (P_g,N_g) law with the same P_g marginal and independent fair N_g. This relative entropy equals that gate's deficit. Data processing and Cauchy–Schwarz bound the sum of mean TV errors in the three-bit laws (U,V,Z=UV) by O(sqrt(m(D+s)))=o(n). Thus on good modes the total errors in finite-alphabet entropies, clipped credits and bounded scores are o(n) by uniform continuity. At rank zero the law stays exact. Rank two has fair independent U,V. Rank one either has one fresh fair U independent of old V with probability q of one, or U fair and V=U+T(P_g) in the same new coset. In the latter case q=Pr(T=0), both marginals are fair, I(U;V)=1−H_2(q) and Pr(Z=1)=q/2. These are comparison laws for the scores. The transcript still uses true inputs.

## 5. Flags and affine forcing

Write \(Z=UV\) for a nonlinear gate on its two signed operands. Its native output is \(Z\) or \(1+Z\). Flag an operand from a nonlinear producer with one native use. Classify the slot by its consumer: \(A\) means one flag and rank zero, \(B\) means one flag and rank one, and \(T\) means two flags. Rank is the raw rank increment \(j\) from §4. A flag has zero formal raw row, so these classes include every slot. Write \(b_g\) for the number of flags at gate \(g\), \(i_g\) for their common class, and \(o_g\) for the class of the flagged slot fed by \(g\). Use class zero when the corresponding flags are absent.

Put \(\phi(q)=\min(q,1-q)\). For each \(i\in\{A,B,T\}\), define
\[
 M_{g,i}=b_g[i_g=i]-[o_g=i],\qquad
 B_{g,i}=[i_g=i]\sum_{F\text{ flagged at }g}\phi(\Pr(F=1))
              -[o_g=i]\phi(\Pr(Z_g=1)).
\]
The two sums over gates vanish on the true law in each mode. This follows by matching each flag with its distinct producer, and \(\phi(q)=\phi(1-q)\). By the variation bound in §4, the ideal laws satisfy \(\sum_g B_{g,i}=o(n)\) on good modes. The \(M_{g,i}\) depend only on the circuit.

An outgoing \(B\) slot has code 4 instead of 2 if its consumer itself feeds a \(B\) slot. Both codes still mean \(B\) in the proof, with the same \(a_B\) and \(\eta_B\). Let \(d_g^{out}\) indicate code 4 and \(d_g^{in}=[i_g=B,o_g=B]\). Matching flagged edges gives zero sums for
\[
 D_g=d_g^{in}-d_g^{out},\qquad
 B_{g,D}=d_g^{in}\phi(\Pr(F_g=1))-d_g^{out}\phi(\Pr(Z_g=1)).
\]
On ideal laws the bias sum is \(o(n)\) by §4.

Fix an accepting input \(x\). Substitute the nonlinear outputs at \(x\) for their formal symbols in every signed operand and in the final output. Each expression is now affine in the raws and each flagged operand is constant. The original equations at a zero-flag gate require one zero operand to stay zero if \(Z(x)=0\), and both operands to stay one if \(Z(x)=1\). At a one-flag gate they require the other operand to stay one or zero as at \(x\) when its flag \(F(x)=1\). At a two-flag gate no equation is needed. The corresponding costs \(L_g(x)\) are \(1+Z_g(x)\), \(F(x)\), and zero. Include the final output equation.

Modify the equations as follows.

- At a class \(B\) gate \(g\), if the producer \(p(g)\) of its flag has fewer than two flags and the unflagged operand \(W_g(x)=0\), replace the old equations by \(W_g=0\) and skip \(p(g)\).
- At a class \(A\) gate \(g\) whose flagged producer \(p(g)\) has zero flags, do the same when \(W_g(x)=0\). Write \(u_g=1\) for this class of consumers and set it to zero elsewhere. Each such producer feeds exactly one such slot, so \(\sum_g(u_g-[b_g=0,o_g=A])=0\).
- At a class \(T\) gate \(g\), let \(t_g\) count flagged producers with zero flags of their own. If \(t_g=2\) and \(Z_g(x)=0\), keep one zero flag and skip the other producer. If \(t_g=1\), call the flag from the zero-flag producer \(U\). When the other flag \(V(x)=0\), keep \(V\) and skip \(U\)'s producer. The consumer needs no equation. Put \(R_g(x)=[Z_g(x)=0]\) for \(t_g=2\), \([V(x)=0]\) for \(t_g=1\), and zero otherwise. Set \(t_g=0\) elsewhere. Then \(\sum_g(t_g-[b_g=0,o_g=T])=0\).

On a flagged chain \(q\to p\to g\), also skip \(q\) if \(g\) is class \(B\), \(W_g(x)=0,W_p(x)=1\), and \(p\) is class \(B\) or has \(u_p=1\). The rule at \(g\) already skipped \(p\), since \(p\) has one flag. The rule at \(p\) has not skipped \(q\) because \(W_p(x)=1\). No other consumer can skip \(q\). These extra skips remove their entire \(L_q(x)\) allowances.

Impose equations only at non-skipped gates. They hold at \(x\). A skipped producer has one native use, at a masked or skipped consumer, and cannot affect any other substituted operand through an affine path. By topological induction each non-skipped gate has the prescribed output. A modified one-flag gate forces its product to zero by \(W_g=0\), whose true value equals the substituted form by induction. A modified two-flag gate keeps one zero flag whose producer cannot be skipped elsewhere and is not the skipped producer for this gate. At all other non-skipped gates the old equations suffice. The final affine equation uses no skipped symbol and forces acceptance.

For each modified one-flag gate the skip removes \(L_{p(g)}(x)\) equations and adds at most \(1-F_g(x)\), also if the consumer was itself skipped. Each class \(T\) skip removes at least one equation. The skipped producers are distinct. Thus the old total \(\sum_g L_g(x)\) minus these gains bounds the new count before the final equation. It must be at least \(n-2k\), since otherwise the solution space contains a \(2k\)-dimensional affine space, forbidden by §A.

For class \(B\), compare the full \((P_g,N_g)\) as in §4 so that the fair fresh \(W_g\) is independent of \(L_{p(g)}\) and \(F_g\). The aggregate error is \(o(n)\). If all such masks ran, the expected gain at \(g\) would be \((\ell_{p(g)}-(1-\Pr(F_g=1)))/2\), with \(\ell_p=E[L_p]\). If \(p(g)\) has two flags, \(\ell_p=0\). Omitting that mask reduces the equation bound by \((1-\Pr(F_g=1))/2\ge\phi(\Pr(Z_p=1))/2\). Charge this lower reduction at \(p\).

At a class \(A\) gate with \(u_g=1\), the gain when \(W_g=0\) is \(1+Z_{p(g)}-(1-F_g)\). If \(F_g=Z_{p(g)}\), its mean is \(2p_{10}\); if \(F_g=1-Z_{p(g)}\), its mean is \(p_{00}+p_{10}\), where \(p_{ij}=\Pr(F_g=i,W_g=j)\). Thus its mean is at least \(A_g=p_{10}+\min(p_{00},p_{10})\). Put \(A_g=0\) if \(u_g=0\).

For a chain with \(p,g\) of class \(B\), compare the full basis law first at \(g\), then at \(p\). The earlier \(L_q,W_p\) are functions of \(P_g\), and \(L_q\) is a function of \(P_p\). Each fresh comparison gives a factor \(1/2\). Thus the extra gain is \(\ell_q/4\), with total variation error \(o(n)\), because each gate appears at most once in each role. Charge this gain at \(q\) as \(d_q^{out}\ell_q/4\). If \(p\) is class \(A\) with \(u_p=1\), only \(g\) needs the fresh comparison. Since \(L_q=1+Z_q\) and \(Z_q\) is \(F_p\) or its complement, the gain is at least \(A'_p=(p_{01}+p_{11}+\min(p_{01},p_{11}))/2\) under the law of \((F_p,W_p)\). Use \(A'_p=0\) unless \(u_p=1,o_p=B\). The checker writes the expression as \(-1/4+3(p_{01}+p_{11})/4+\phi(p_{01}+(p_{00}+p_{10})/2)/2\).

With \(\sum_g j_g\le n\), the equation bound gives \(\sum_g C^\ast_g\ge-o(n)\), where on ideal laws
\[
 C^\ast_g=\ell_g-j_g-\tfrac12[o_g=B]\ell_g
  +\tfrac12[i_g=B](1-\Pr(F_g=1))
  -\tfrac12[b_g=2,o_g=B]\phi(\Pr(Z_g=1))-E[R_g]-A_g-A'_g-\tfrac14 d_g^{out}\ell_g .
\]
The aggregate variation estimate of §4 applies to every displayed local term.

## 6. Query capacity and physical entropy ledger

Put \(a_0=0,a_A=a_B=1,a_T=1/2\). At a zero-flag gate, choose \(r\in\{U,V,U+V\}\) minimizing \(H(r)+a_oH(Z\mid r)\) in the mode. At a one-flag gate query its flag, and at a two-flag gate choose one flag uniformly, independently at each gate. Fix the choices before seeing the input and extend the schedule to the full cube. The fiber of each successive query word is affine. Inductively earlier signals are affine on the current fiber. Fixing \(U\) or \(V\) makes \(UV\) affine. Fixing \(U+V\) gives \(UV=U\) at parity zero and \(UV=0\) at parity one. XOR and XNOR preserve affinity. Intersect the final fiber with output one. This is an accepting affine subspace of dimension below \(2k\) by §A, and its intersection with the mode has at most \(2^{2k}\) inputs. Thus the entropy of the query word is at least \(H(X\mid\sigma)-2k=n-o(n)\) in each fixed schedule on a good mode.

In the chain rule for the query word, bound each zero-flag query by its entropy. A flag query from \(p\) costs at most \(H(Z_p\mid Q_p)\) because \(Q_p\) is an earlier query. Each producer has at most one slot of this type, selected with probability \(a_{o_p}\) independently of \(Q_p\). Averaging the upper bounds proves \(\sum J_g\ge n-o(n)\), where \(J=\min_r(H(r)+a_oH(Z\mid r))\) for zero flags and \(J=a_o\,\operatorname{mean}_{F} H(Z\mid F)\) for positive flag count. Uniform continuity transfers this bound to the ideal laws.

Write \(D'(Y)=1-H(Y)\) and \(I'=I(U;V)\). Save the producer's \(D'\) on each physical carrier. Add \(I'\) on each retained designated final arm from §3. The total saved weight is \(S_0\). Each arm weight is a marginal or a conditional entropy, so each saving is in \([0,1]\). First assign zero input deficit to flagged arms. In the algebraic ledger use \(c_z=1-\eta_o\) on outputs, \(c_u,c_v=\eta_i\) on flagged inputs, and one on unflagged inputs. The signed numerators are \(0,6000,6000,52792516\) over \(10^8\). On each flagged edge the signed correction cancels exactly because its two bits agree up to complementation. The physical weights are unchanged. For \(d\ge2\) native uses, a gate has \(2d-1\ge d+1\) carriers available for its output deficit and consumer deficits. A one-use producer shares the deficit exactly when the consumer has a flagged nonlinear input. Affine producers have no output deficit to reserve. Raw \(d\ge3\) has \(2d-3\ge d\) carriers, raw \(d=1\) has one, and raw \(d=2\) misses at most one credit, totaling \(o(n)\) by §4. By §3 only \(o(n)\) of the mutual information is lost on missing designated arms. Therefore \(S_0\ge\sum_g E_g-o(n)\) where \(E_g=c_zD'(Z)+c_uD'(U)+c_vD'(V)+I'\). The same inequality holds with ideal local laws by §4.

At every nonlinear gate put \(\Delta_g=\min(I(U;V),H(U,V\mid Z))\). Section 7 earns this additional credit on the exposed gates and shows that dropping the others loses only \(o(n)\). The checker in replay.py proves
\[
 E_g+\Delta_g\ge \beta j_g+\theta_\ast C^\ast_g+\kappa J_g
  +\sum_{i\in\{A,B,T\}}(\xi_iM_{g,i}+\gamma_i B_{g,i})
  +\zeta_T(t_g-[b_g=0,o_g=T])+\zeta_A(u_g-[b_g=0,o_g=A])+\xi_D D_g+\gamma_D B_{g,D}+\sum_{h,i}\gamma_{h,i}B^{(h)}_{g,i}+\mu .
\]
The numerators over \(10^8\) for \(\beta,\theta_\ast,\kappa,\xi_A,\xi_B,\xi_T,\gamma_A,\gamma_B,\gamma_T,\zeta_T,\zeta_A,\xi_D,\gamma_D,\mu\) are \(6427690,3616495,5895441,-4469698,-9346728,27190256,34592891,35844536,-13423863,1048845,1727650,1364857,9489348,100\). One-sided multipliers are nonnegative. For \(h=9/32,21/64\), form \(B^{(h)}_{g,i}\) from the bias measures \(B_{g,i}\) by using \(\min(q,1-q,h)\) on both sides of the flag. Include \(i=A,B,T,D\). Complement invariance gives exact balance before ideal replacement, hence aggregate \(o(n)\) after it. The numerators for \(\gamma_{h,i}\) in that order are \((-50623992,-4675462,-43951127,-7299297)\) and \((28777192,-2351967,-17985221,-9098469)\) over \(10^8\).

For rank zero use \(U\) for the flag at a one-flag gate, and order corners \(00,10,01,11\). At zero flags bound \(J\) by the average of \((1-a_o)H(r)+a_oH(r,Z)\) over the two operand queries \(U,V\). Each joint law has masses \(p_{11},p_x,1-p_{11}-p_x\) at one of the other corners \(x\). At positive flag count, \(H(Z\mid F=0)=0\). At rank one, \(E=c_z(1-H_2(q/2))+w(1-H_2(q))\) with \(w=1\) if \(b=0\), and \(w=\eta_B\) if \(b=1\). For \(b=0\), \(\ell=1+q/2\) and \(J\le H_2(q)\) if \(o=0\) or \(J\le1+a_o H_2(q)/2\) otherwise, by querying the old form or a fresh operand respectively. For \(b=1\), \(\ell=q\) and \(J=a_o q\). At rank two, \(\ell=5/4\), \(J=1+a_o/2\) and \(E=c_z(3\log_2(3)/4-1)\).

The checker subtracts the proposed local ledger bound. For rank zero it checks both branches of the added minimum, using the four corner masses as simplex coordinates. Binary and joint entropies are sums of \(\phi(x)=-x\log_2 x\), with \(\phi(0)=0\), at affine masses. Conditional entropies are differences of joint and marginal entropies. Combine any shared same-sign coefficients on complementary masses into binary entropy. Its interval values add the two value bounds, and its derivative subtracts the complementary derivative. For rank zero with zero or two flags the checker starts with \(p_{10}\le p_{01}\), cutting at the midpoint of \(10,01\). The slack is symmetric except for \(t=1\), where larger \(V\) lowers \(E[R_g]\) and increases the right side. The other terms are symmetric. At rank one it drops the nonnegative added minimum and checks \(0\le q\le1\). Rank two is a single point.

On each simplex cell, an affine mass has its minimum and maximum at vertices. For positive coefficients, concavity bounds each entropy from below by the affine interpolation of its values at the vertices. For negative coefficients, use an upper tangent at an interior mass before multiplying. For the bias minima, including the clipped ones, use concavity at the vertices or an affine upper branch. Check the resulting affine lower bound at each vertex with directed rounding. If any comparison fails, bisect a longest edge and check both children. For rank zero the squared distance includes differences of the operand marginals and product mass. The two children cover the parent, and the depth limit keeps these dyadic coordinates exact.

For positive \(x\), write \(x=2^e u\) with \(1\le u<2\). For \(w=(u-1)/(u+1)\), enclose \(\ln u\) by the first twelve terms of \(2\sum_{j\ge0}w^{2j+1}/(2j+1)\) and the positive tail bound \(2w^{25}/[25(1-w^2)]\). The same enclosure at \(u=2\) gives \(\ln 2\). These intervals bound \(\log_2 x\) and \(\phi'(x)=-\log_2 x-1/\ln 2\). The strict comparisons over the complete partition permit a fixed increase in the reserve.

Summing the local inequalities and using the rank, query, forcing and balance bounds gives \(S_0+\sum_g\Delta_g\ge s_0n-o(n)\), with \(s_0=\beta+\kappa+\mu=12323231/100000000\). The ideal-law replacement affects \(\Delta_g\) by aggregate \(o(n)\), by its continuity. At least \(n-o(n)\) nonlinear gates contribute the reserve.

### Contract private-input leaves and output stems

A private raw has one native use. Mark each private recipient gate with exactly one private raw whose other true input was not absorbed into the same bag in Section 3. At most o(n) private raws lie outside marked recipients, by types (i) and (ii). For each marked g, let X_g be its private raw bit, Y_g its native output, and V_g its other true input. Its outgoing stem is the physical edge leaving g, including the edge to the dummy leaf when g is the final gate. Marking uses only the circuit graph.

At each accepting word, fewer than 2k marked nonlinear recipients have the signed value of V_g that makes the gate constant in X_g. Otherwise change any subset of 2k of those private raws. All native gate outputs remain unchanged by induction in topological order. Each changed raw occurs only at its own marked recipient, where the unchanged control makes the output constant in that raw. The resulting accepting coordinate cube splits into two sets of size 2^k whose sum lies in S, as in Section 5. Thus, in any positive mode, the sum of the probabilities of these constant cases is at most 2k. Affine gates have no such case.

Group marked recipients by the topological index of their earliest native consumer, with the dummy leaf last. For one group, let A_i be its set of private raws. There are at most two such raws since a consumer takes at most two pins. Let P_i consist of all raws outside the marked recipients and the private raws of earlier groups. The input V_g of any g in group i is determined by P_i. In fact any other marked output used to compute it is used by a gate strictly before the earliest consumer of g. The earliest consumer of that other marked output is therefore earlier than group i. A private raw itself can enter the calculation only at its own recipient.

The probability cost of a tag on g's outgoing stem can use the consumer's other operand. Such a tag requires a sole native consumer, since with fan-out greater than one the stem goes to a copy fork. When the stem of g is tagged, let H be that consumer. Section 3 gives a private raw priority over a producer stem, so H has no private raw input. No marked output from a later group can affect the other operand of H. Any marked output from group i that affects it must be the direct other operand, because any use strictly before H would put its first use earlier. The other operand cannot be Y_g itself after normalization. Consequently the other operand is determined by P_i and the other raws in A_i. The claim also covers signed operands, because complementation preserves the relevant entropies.

Write E_g for the conditioning on P_i and A_i without X_g, in addition to the fixed mode. The chain rule gives

    1-H(X_g | E_g) <= |A_i|-H(A_i | P_i,mode).

The sum of the left sides over marked recipients is at most twice the entropy deficit of all marked raws given the other raws and the mode. Its mean over modes is O(D+s) by Section 4. Thus the sum is o(n) on good modes. On a marked raw leaf the total saved weight is at most 1-H(X_g | E_g), since E_g determines the other input of g. On its outgoing stem the saved weight is at most 1-H(Y_g | E_g). For an untagged stem this follows by conditioning its marginal entropy. For a tagged stem it follows from the other-operand observation and D'+I'=1-H(Y_g | other operand,mode). Conditional on E_g, the gate is a signed identity in X_g except in the constant cases above. Therefore

    1-H(Y_g | E_g) <= 1-H(X_g | E_g) + Pr(g is constant in X_g | mode).

The marked leaves and stems together carry o(n) saved weight. Unmarked private leaves number o(n) and each carries bounded weight.

Absorb all private raw leaves into their recipients. Each marked bag now has degree two. Select its outgoing stem for a further contraction. A marked head has at most one selected incoming edge, and an unmarked head at most three. The selected edges cannot form a directed cycle because they join producers to consumers or to copy forks in their native order. The selected components are therefore chains, with at most three chains ending at one unmarked root. Make bounded bags by partitioning each chain into consecutive pieces. Allow the last piece of each incoming chain to share the root bag. For example, away from the root use pieces with original mass between B/8 and B/4, except for the last piece of each chain, when B is sufficiently large. This leaves O(n/B) selected edges between pieces, since every nonterminal piece has mass at least B/8. Every resulting bag has mass at most B. Each merger uses a degree-two end and preserves degree at most three and rank r. Off the contracted chains the only initial degree deficiencies are the dummy leaf, the marked bags and o(n) exceptional private recipients. Each contracted marked bag removes one unit of deficiency by the degree formula in Section 2. The final edge and savings totals satisfy

    E <= 3r + o(n) + O(n/B),       S >= S0 - o(n).

Here the first inequality uses the degree sum and connectedness. The contractions only change which gates and raw inputs belong to one party, and all evaluations still use the original native order.

## 7. Causal weighted transcripts and strict contradiction

### Exposed gates

Call a nonlinear gate exposed if its final bag is that gate alone, both actual final operand arms survive, and one arm is designated. Absorb optional degree-two raws into their earlier consumer. Keep the private-priority orientation from Section 3. We show that all other gates have total mutual information \(o(n)\). For a raw \(X\) absorbed into a later owner, its earlier recipient \(h\) has a single private raw \(Y\), except at the \(o(n)\) bundle exceptions. The recipient \(h\) is marked since its \(X\) arm survived. If \(h\) is nonlinear, it is constant in \(Y\) with probability at least \(1/4\) when both probabilities of \(X\) are at least \(1/4\). The good raw entropy bound and the private constant-case sum show that only \(o(n)\) such earlier recipients are nonlinear. For each remaining forced-late owner, change coordinates \((X,Y)\mapsto(X,H)\), where \(H=X+Y+\epsilon\) is the native affine output of \(h\). The pairs are disjoint by native use counts. This is an invertible raw chart preserving the mode entropy. In the analysis only, supply \(H\) directly at \(h\). Then \(X\) is first used at its later owner.

Group chart coordinates by first native use, introducing \(H\) at \(h\). Each group \(A_i\) has at most two coordinates. If \(X\) is an absorbed degree-two input first used at its owner, the other operand \(V\) is determined by the earlier groups \(P_i\) and \(A_i\setminus\{X\}\). Hence \(I(X;V)\le1-H(X\mid P_i,A_i\setminus\{X\})\le |A_i|-H(A_i\mid P_i)\). The sum is at most twice the chart deficit. An absorbed private raw has mutual information \(o(n)\) in aggregate by the same estimate with the private groups and \(E_p\) from Section 6. A selected stem into a copy fork removes no native operand arm. If it removes an arm from a marked producer \(p\) into its sole native consumer \(g\) without a private raw, the other operand of \(g\) is determined by \(E_p\) by the same earliest-consumer argument, whether or not the arm was designated. Thus this mutual information is at most \(1-H(Y_p\mid E_p)\), summing to \(o(n)\) by the private-identity estimate. A consumer with a private raw was already covered. These are all ways that a nonlinear gate acquires other bag contents or loses an operand arm. Section 3 also bounds the mutual information at gates with no designation by o(n) in total.

### Capacity per direction

At an exposed gate use \(j=H(U,V)\), \(a=H(U\mid V)\) for its designated \(U\), \(b=H(V\mid U)\), \(d=H(Z)\), \(D=\max(a,b,d)\), and \(B=\max(b,D-a)\). For either direction let \(A_U,A_V\) indicate incoming crossings. Replace that gate's two input charges by

\[
 F=(D-B)A_U+(D-a)A_V+(a+B-D)\mathbf 1[A_U\lor A_V].
\]

All three coefficients are nonnegative. With receiver bits one the last indicator is \(z(1-xy)\), where \(x,y\) are operand-supplier colors and \(z\) is the consumer color. Its three mixed differences are nonpositive, so it is submodular. Simultaneous complementation and identifying parallel suppliers preserve submodularity. The four capacities are \(0,a,B,D\). The old capacities are \(0,a,H(V),j\). The saving is at least \(\Delta_g\) whenever \(A_V=1\), because \(H(V)-B=\min(I,j-D)=\Delta_g\) and \(j-D\ge\Delta_g\). Here \(j-D=\min(H(U),H(V),H(U,V\mid Z))\). Before fiber folding every tree arm has the same directional crossing probability \(p\le c/3\), by the unmarked fair law. Expected cost is therefore at most \(p(E-S-\sum_{\rm exposed}\Delta_g)\). Fold both submodular directional sums and take the same balanced cut by concentration. This gives the bound with \(c/3\) in place of \(p\) and an \(o(n)\) error.

### The exact simulation

At a nonexposed gate, use the original carrier plan. For each crossed designated arm, replace the send with a private incoming guess and retain the would-be supplied bit. If exactly one operand of an exposed gate crosses, guess it after the other operand is available. If both cross, evaluate the gate on the common supplier side and send its output to the original gate side at its native step. Omit both input sends. Route all later uses from the original-side output slot, even at another moved gate. The nonempty cases charge at most \(a,B,D\) respectively, since the moved output costs \(d\le D\). Signs preserve these entropies.

After simulation exchange both shares of the guess-equality vector and require equality. The schedule and equality give the true output slots by topological induction. Each true raw pair has one passing guess pair. Announce at most two producer values per pin indicator to check the mode. Check the final output and the two following cost tests.

For each receiver sum the negative binary log probabilities of its ordinary receives and of its guesses given the other operand. Use the mode laws, mixed with a fixed small \(\rho>0\) of the fair law, or use the fair law alone on a null conditioning event. On true inputs the mean cost per receive exceeds its entropy by at most \(-\log_2(1-\rho)\). Each cost is a bounded function of a listed producer bit or of the signed-pair indicators in §4. The extra output occurrence per moved gate fits the producer-copy budget. For the total capacity \(W_{\rm dir}\) in either direction, §4 gives cost at most \(W_{\rm dir}+\epsilon n\) in both directions on all but \(o(1)\) of the accepting mass, after choosing \(\rho\) sufficiently small compared to fixed \(\epsilon>0\).

A passing transcript is a Cartesian leaf in the extended raw-plus-guess domains, by local consistency at each message. Each projection from a passing leaf to the side's raw words is injective: otherwise fix the other side's extended word and contradict uniqueness of the passing guesses. Each projected leaf is an accepting rectangle, with one side of size less than \(K\).

Fix the receiver's raw word. The receiver can reconstruct its trajectory and outgoing messages from its preceding received values and guesses as they are used. These choices form a causal binary tree, also on invalid histories, with the normalized smoothed laws at each choice. A passing causal path has probability at least \(2^{-W_{\rm dir}-\epsilon n}\). Its receiver-calculated equality share fixes the sender's equality share. Pins and tests add at most \(2s+O(1)\) bits. Thus there are at most \(2^{W_{\rm dir}+\epsilon n+2s+O(1)}\) passing leaves for this word.

Suppose a fixed \(\tau>0\) and subsequence has \(t\le(L+d_0-\tau)n\), where \(d_0=s_0/3\). The contracted ledger and extra credits give \(S+\sum_{\rm exposed}\Delta_g\ge s_0n-o(n)\). The two directional bounds, \(E\le3(t-n+1)+o(n)+O(n/B)\), and \(n/2=c(L-1)n\) give \(q_{\rm dir}-W_{\rm dir}\ge\gamma n\) for fixed \(\gamma>0\) when the bag mass cap \(B\) is sufficiently large. Here \(q_{\rm dir}=n/2+o(n)\) counts sender raws. The mode cut is fixed before the word is seen. With \(\epsilon<\gamma/2\), count each leaf in a direction with a thin sender. It contains fewer than \(K\) sender words for each fixed opposite word, and there are at most \(2^s\) modes. The covered full-cube mass is \(2^{-\Omega(n)}\) because \(s+\log_2 K=o(n)\). Dividing by the accepting density \(2^{-D}\) gives \(o(1)\) accepting mass, contrary to passing mass \(1-o(1)\). This proves the limit inferior bound.

The actual input coefficient is
    L+d0 = 1+1/(2*(12252408923603907442/10^20))+12323231/300000000
         = 9413355971276602268255192551/1837861338540586116300000000.
Since n/N tends to one, this is also the coefficient for actual length N. The prefix certifies the first continuation endpoint. The terminal checks both entropy partitions, 46,081 cells with (b,o)=(1,3) and 94,706 complementary cells. The strict comparisons permit a fixed increase in μ. Applying the proof with that larger reserve gives a strictly larger limit inferior. Therefore the declared lower bound holds for every sufficiently large actual length, with zero remainder.

## Provenance

The family follows classical norm graphs and Arthur Ostapenko 808da377. The immediate source is Yukon c716d171-ab2f-424d-8369-b7df52e0fa61, evaluated at 76e8c5b8e2e7cb25ee94e95ac7cd66f8167a093c. It supplies the graph certificate, linked queries, transcript, and carrier methods. This version adds private-stem contraction, two-step masks, movable gates, the continuation certificate, and the strengthened Gaussian moment and repair checks.

## Graph addendum. A strict additional repair on a single C edge

For a C edge with raw pins 11, order its paired perpendicular rays with independent labels. Test the six prefixes in gain.py. The 18 rows are disjoint and the checker tests every exterior trit. A matching pin stays one exactly when more endpoints have C degree two than zero. The qualifying edge is positive, with C degrees one and two in order and both matching bits zero. Only one orientation qualifies. The first endpoint charges its unique positive edge. Double the sum of the 18 directed Gaussian lower bounds for the two C incidences per vertex. This gives gain at least 17120035/10^12. The combined count below allows overlapping deletions.

For a marked path let w_i=t for a straight internal passage and w_i=a for a turn. Set b_i=a² or at respectively. Expanding the incoming independent Gaussian modes gives Cov(R_u,R_v)=(a+t/2) prod_i w_i + (1/2) sum_i b_i prod_{j!=i} w_j. The endpoint modes contribute the first term. The two off-path modes at each internal vertex have combined coefficient product 2b_i before the two score factors 1/2. The diagonal is one. The checker uses the exact precision matrix and the directed halfnormal bound of graph.md §3.

The one-edge repair and its marked-path covariance first appeared in yukonautominer’s Yukon submission 086221d5-56f1-4068-9241-b118a66143fc, source commit 56d0efa296a3af3b8d0c478c5ef96165be237882. This version uses a larger disjoint cylinder sum and a direct bound on the majority loss.


## Combined repairs and longer tails
For each tail length 8, 10, 12 or 14, the two opposite signed patterns are disjoint. Apply lower_signed to their marked-line covariance and signs; take the maximum of the summed lower bound and the Jensen bound.

Let d_v be positive degree after matching repair. Deleting any positive edge set D saves sum_{uv in D}(4-d_u-d_v) plus the number of adjacent edge pairs in D, by expanding the split indicator (3d-d²)/2. The singleton repair above has individual gain one. Also select both C edges centered at each vertex satisfying one of the 79 adjacent_floor cylinders. Order rays by independent continuous labels. The checker tests contradictory pins for disjointness, every possible exterior trit, positive edges, and nonnegative individual gains in 5,976 cases. Each center supplies a distinct extra pair charge. All tests use the original pre-motif state, even when selections overlap.

With Gaussian covariance C, K=C^-1, nonnegative lower triangular L of positive diagonal, nonnegative b, U=diag(sign)L and z=T*1+diag(sign)b, the map x=z+U y for independent standard halfnormals y injects into the cylinder. Jensen bounds its mass below by
    det(L)/(2^n sqrt(det C)) exp(n/2-tr/2+(tr-q)/pi-lin sqrt(2/pi)-z^T K z/2),
where tr=trace(U^T K U), q=1^T U^T K U 1 and lin=z^T K U 1. This uses E y_i²=1 and E y_i=sqrt(2/pi). The code constructs C from the independent Gaussian innovations of graph.md, checks the rational inverse, bounds pi and square roots outward, and uses an odd Taylor exponential minorant with downward squaring. The disjoint sum is at least 26772897894/10^15.

Activate each qualifying singleton or center independently with probability p and delete their union. Each touches at most three vertices and pays a distinct charge; hence saving is at least one third of majority loss. Full activation saves at least the certified sum B of the two gains. The pre-motif marginal is verified above .503 and B<.001. If full repair takes the marginal to at most half, continuity gives a fair p with saving greater than B. Otherwise apply all repairs and vary the existing terminal patch activation to reach half; this never increases energy. No positive bits or matching pins are created, so the same terminal relaxation applies. The exact rational guard in energy.py checks that the unrounded split-energy upper bound is strictly below c. The printed 20-digit ceiling can equal c. Choose finite approximations that preserve this strict inequality, the early-stop margin, and the relaxed terminal majority below one half. Then choose the activation parameter by continuity on finite neighborhoods.

Attribution: linked-query and terminal code from xlib, Yukon 7596d878-f88a-497e-a98c-277cd59d7045, source c3ae7d9af7d80fdc1407c1a5279b3f48b8d64233. Longer-tail tilts from Holindauer, Yukon 5a954ca2-0b1d-48d4-834c-cc3fb897cd10, source 54f02a09674b425383aea12d494dfd0491b35b86 (public PR 937; no prior verdict assumed). Adjacent repair from the lab's directed cone certificate. All numerical premises are recomputed, not imported judgments.
