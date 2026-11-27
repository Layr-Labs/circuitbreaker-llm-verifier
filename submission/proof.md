# A full-B2 bound for the norm-one predicate


## A. Fixed family and eligibility

The family uses the norm-graph idea of Kollár, Rónyai and Szabó, Combinatorica 16 (1996), 399–406, DOI 10.1007/bf01261323, and the norm-one languages used in earlier CircuitBreaker entries.

### A.1 A fixed uniform evaluator at every actual input length

For n<64 define f_n=0, including n=0. For every integer n>=64 put

```
L=floor(log2 n), B=ceil(log2(L+1)), a=L-B, p=2^a,
ell=floor(n/p), m=p*ell, d=2^floor(a/2), k=m/d, Q=2^k.
```

Here $a\ge3$, $d\ge2$, $d\mid m$. The dyadic bounds give $2^B\le\ell<2^{B+1}<4(L+1)$, $0\le n-m<p\le n/(L+1)$ and $p\ge n/[4(L+1)]$. Thus $m/n\to1$, $d=\Theta(\sqrt{n/\log n})$, and $k,d\log_2d=O(\sqrt{n\log n})=o(n)$.

Build the degree-$\ell$ field from the first monic binary polynomial with no divisor of degree $1,\ldots,\lfloor\ell/2\rfloor$. The roots of $X^{2^\ell}-X$ form a field whose proper subfields cover at most $\sum_{j\le\ell/2}2^j<2^\ell$ roots. An irreducible exists and the division search is polynomial since $\ell=O(\log n)$. Start with its power basis.

At each of $a$ steps select the first basis element $u$ of trace one. The linear trace is nonzero since its polynomial on the degree-$r$ field has degree $2^{r-1}$. The map $z\mapsto z^2+z$ has kernel two and image equal to the trace-zero space. Thus $X^2+X+u$ is irreducible. Append $\beta$ times the old basis for its root $\beta$. Reducing $\beta^2=\beta+u$ uses at most five old multiplications per new multiplication. There are $O(\log n)$ levels, so arithmetic and trace scans are polynomial.

Read the first $m$ bits as $z$ in this basis. Ignore the rest. Accept precisely when $z^e=1$, with $e=(2^m-1)/(Q-1)=\sum_{j<d}Q^j$. Its $m$-bit exponent permits polynomial repeated squaring. The fixed NP verifier applies this total evaluator and requires the witness empty, with witness bound zero.

### A.2 Accepting density and essential input count

Since $X^e-1$ divides the split squarefree polynomial $X^{2^m-1}-1$, exactly $e$ inputs accept on the active bits and density is at least $2^{-k}$. Odd $e$ makes every active bit essential, because otherwise flipping it pairs the accepting inputs. Fix ignored bits to zero and normalize by §3, without increasing gates.

### A.3 Elementary common-neighbor lemma

For $a_{ij}\in K$ distinct down each of $d$ columns, the equations $\prod_j(z_j-a_{ij})=1$ have at most $d^d$ common solutions. Set $L_{ij}=X_j-a_{ij}Y_j$, $W=\prod_jY_j$ and work modulo $\prod_jL_{ij}-W$. The $d$ forms $Q_{ij}=\prod_{r\ne i}L_{rj}$ span degree $d-1$ in column $j$ by evaluation at $(a_{ij},1)$. Their monomial multiples span every higher degree.

For a form of degree $D\ge d$ in each pair, extract $Y_j$ successively. After earlier factors are extracted, terms missing $Y_j$ contain $X_j^D$. Expand each other column in the $Q$ forms times monomials. Since $d-1$ omissions leave a common row $i$, every product contains $L_{ih}$ for $h\ne j$. Write $X_j^D=(L_{ij}+a_{ij}Y_j)X_j^{D-1}$. The $L$ product becomes $W$, giving $Y_j$ in both terms. All columns reduce to $W$ times degree $D-1$.

Set all $Y_j=1$. Pad degrees by homogenizing with $Y_j$ and repeat, so the affine quotient spans at most $d^d$ monomials of exponents below $d$. For any finite set of $s$ solutions, evaluation surjects onto $K^s$: at a chosen point multiply normalized linear factors vanishing at the other individual points, using a differing coordinate per factor. Thus $s\le d^d$.

### A.4 Additive sumset exclusion on the entire field

For $A+B$ accepting and distinct $u_i\in A$, $1\le i\le d$, Frobenius gives $\prod_{j<d}(v^{Q^j}+u_i^{Q^j})=1$ for $v\in B$. The column entries $-u_i^{Q^j}$ are distinct. Apply §A.3, with first coordinate $v$, to get $|B|\le d^d$. Put $K=d^d+1$, $k=\lceil\log_2 K\rceil$. Hereafter $n$ is active length and $N$ actual length, $n/N\to1$ and $k=o(n)$. The active set $S$ has density $2^{-D}$ for $D=o(n)$. It excludes sumsets with both factors of size at least $K$, hence accepting rectangles have a side smaller than $K$ by embedding in complementary coordinate subspaces. Use uniform $\nu$ on $S$, with $H_\nu(X)=n-D$ in bits.

Let c=12233148667589113194/10^20, L=1+1/(2c), d0=13528991/300000000. The only numeric graph input is the finite-local matching-equivariant fair tree law with edge probabilities a_M,a_C satisfying a_M/2+a_C<c. Fix its shell, window and activation parameters before sending circuit size to infinity.

## 1. Elementary finite-local almost-perfect matching and edgewise transfer


Start unmatched on the cubic tree. Fix $J=2l+1$. Each round augment all alternating paths of length at most $J$ with unmatched ends that form whole low components of fresh fair marks. They are disjoint. Send half a unit from each changed vertex to each end. Expected transport bounds total changes per vertex by $(J+1)/2$, giving stabilization. A surviving path has a uniformly positive selection chance every late round, a null event. Let $U$ be unmatched at the limit, density $u$. For even alternating reachability $X_i$ within $2i$, let $Y_i=N(X_i)$. For $i\le l$, each $Y_i$ matches into $X_{i+1}\setminus U$ by predecessor or extension, else a short augmentation exists. Densities satisfy $Y_i\ge X_i$ by cubic incidence and $X_{i+1}\ge Y_i+u$ by matching. Thus $(l+2)u\le1$. Choose $l$ then finite rounds for small unmatched density $u_r$.

Use good balls with their incident matching data. A connected union of good balls extends to a perfectly matched tree by keeping prescribed partners and pairing unmatched vertices outward on attached branches. At an anchor use uniform $1-(1-\min(U,V))^2$ from independent endpoint uniforms. Use fair fallback on bad balls. The edge matching probability is $(1-u_r)/3$, bad probability $O(u_r)$. Total crossing is at most $(a_M+2a_C)/3+O(u_r)\le2c/3$ by strictness. Fair marginals give at most $c/3$ per direction.

Complete each missing half-edge by attaching a subdivided $K_4$ at its subdivision vertex. Loops count twice. New costs and raw weights are zero. Each component has the cubic tree of reduced incidence walks as cover. Its deck group acts transitively on each vertex fiber, also for loops and parallel edges.

Lift the equivariant radius-$R$ outputs with independent input labels on cover vertices or incidences. Enumerate their base fibers $F_i$. Give each a fresh independent uniform $U_i\in[0,1)$. At step $i$, fix earlier uniforms and remaining original labels outside $F_i$. For every output, including previously changed outputs, integrate its current bit over labels in $F_i$ to get $p_v$, then replace it by $[U_i<p_v]$. The integral uses only original support coordinates. Inductively each output uses surviving labels in its original support and uniforms of fibers meeting it. If a fiber misses that support, its integral is the existing bit and the update leaves it fixed.

For a local submodular $f$, sort conditional marginals in decreasing order with prefixes $S_j$. Diminishing returns gives $f(A)\ge f(\varnothing)+\sum_{j\in A}[f(S_j)-f(S_{j-1})]$. The common uniform attains this bound with unchanged marginals, so every specified payoff has nonincreasing mean. Deck maps permute each $F_i$ and preserve product measure. The finite integral recursion commutes pointwise with this permutation of the remaining labels. Only the $U_i$ remain at the end and deck maps fix them, so outputs agree pointwise on each vertex fiber. No intermediate equality is needed.

Each finite vertex uses uniforms only in its projected radius-$R$ support. Receiver stars and wires use at most radius $R+1$. Lift each term through its own incidence star. Bounded degree and fixed $R,B$ give $O(n)$ intersecting support pairs. Disjoint supports, also across components, are independent. Thus each directional sum and the raw mass have variance $O(n)$. Raw weight per bag is at most $B$, and fair means give raw mass mean $n/2$. Chebyshev at $n^{2/3}$ and a union bound give one cut meeting all three expectations within $o(n)$.

## 2. Bounded bags in the circuit graph

Represent fan-out by binary copy trees at no gate charge, counting repeated pins separately. Include a dummy output leaf. Before expansion a pruned circuit with $n$ raws and $t$ gates has $n+t+1$ vertices and $2t+1$ wires. Copy expansion preserves cycle rank $t-n+1$. Contracting a nonloop edge preserves rank and adds the two degrees minus two. Other parallel edges become loops, each counted twice in degree.

## 3. Raw-bundle exclusions and physical accounting

Prune dead gates and absorb constants, repeated pins, projections and negations. Retained gates use both inputs, hence are signed AND or XOR. For a signed terminal projection keep its producer with the chosen sign. Essentiality excludes a constant or raw output. All gates reach the output and $t$ cannot increase.

Let $u$ count private bits, meaning degree-one raws. Consider disjoint bundles within each type: two private bits at one gate, a degree-two raw at two recipients with distinct private bits, or two degree-two raws with distinct private recipients and one common further consumer. Only their listed gate outputs leave a bundle. With $r$ raws and at most $r-1$ outputs, at most half of cube assignments have singleton fibers. With $\Omega(n)$ bundles, fewer than $2k$ nonsingleton fibers has cube mass $2^{-\Omega(n)}$ by the independent binomial bound, hence accepting mass $o(1)$. Otherwise fix outputs and outside raws. Split $2k$ varying bundles into two $k$-groups to contradict §A.4. Each type occurs $o(n)$ times.

A gate with $d$ uses, final use included, has $2d-1$ carriers. For raw degree $d\ge2$ merge the leaf, first fork and a child. At $d\ge3$ choose a child fork to preserve arms, and at $d=2$ absorb into one consumer. This leaves $2d-3$ carriers. Degree one has one. Bags know their at most two raws initially. They are trivalent except the $u$ private leaves and dummy. Rank is unchanged.

Absorb a degree-two raw into its other consumer if it has one private recipient, otherwise into its earlier consumer. Both-private cases are $o(n)$ by §3. Designate an arm in order private, nonraw, then degree-at-least-three raw, with fixed ties. Omitted nonlinear pairs have two degree-two raws. Greedy coloring of their degree-two multigraph gives three matchings. In each, subadditivity bounds $\sum I(U;V)\le\sum(2-H(U,V))\le n-H(X\mid\sigma)$, unchanged by signs. Thus omissions cost $o(n)$ on good modes.

## 4. Pinning and high-entropy innovations

List $O(n)$ bits: $2d+2$ producer copies at use count $d$, raw and terminal entries, and four indicators per nonlinear gate. Reveal up to $\sqrt n$ uniform entries with replacement. Mean fresh conditional entropy drops by mean conditional mutual information with the next entry, so some prefix drop is $O(n^{-1/2})$. Fix its indices and let $\sigma$ be their values. Log-sum and the binary second derivative give $D_{\rm KL}\ge2\,\mathrm{TV}^2/\ln2$ in bits. Absolute covariance is at most product TV. Cauchy–Schwarz bounds $E_\sigma\sum_{i,j}|\operatorname{Cov}(W_i,W_j\mid\sigma)|=O(n^{7/4})$. Chebyshev gives $o(n)$ deviations on good mass for fixed sums with bounded mode-dependent weights.

The mean deficit $n-H(X\mid\sigma)$ is $\le D+s$ by the chain rule. It bounds the sum of raw marginal deficits, $|U|-H(X_U\mid X_{U^c},\sigma)$, and $R-H(N\mid\sigma)$ for $R$ independent basis forms $N$ completed to a full chart. Markov makes each loss $o(n)$ outside $o(1)$ mode mass. Below use fixed positive modes.

Let $m$ count nonlinear gates. Query one affine operand chronologically at each, then the output. Every cell has dimension at least $n-m-1$ and an accepting one cannot have dimension $\ge2k$ by splitting a basis. Thus $m\ge n-o(n)$.

Treat earlier nonlinear outputs as formal symbols. Each signed operand has a linear raw form and an affine function of these symbols. The raw rows at $g$ add rank $j_g$ beyond preceding rows. Fixing a basis of $R=\sum j_g$ forms and the final raw row fixes output, so $R\ge n-o(n)$.

Let $N_g$ be the new forms, $P_g$ the earlier basis. Previous nonlinear outputs depend only on $P_g$. The sum of mean conditional deficits $j_g-H(N_g\mid P_g,\sigma)$ is $R-H(N\mid\sigma)\le D+s$ (averaged over modes). Each is relative entropy to the preserved $P_g$ marginal times fair innovations. Data processing and Cauchy–Schwarz give aggregate TV $O(\sqrt{m(D+s)})=o(n)$ on good modes. Bounded continuous scores transfer. Rank zero is unchanged and rank two gives fair independent operands. At rank one $U$ is fair and either independent of old $V$ or $V=U+T(P_g)$. In the latter case $q=\Pr(T=0)$ gives $I=1-H_2(q)$ and $\Pr(Z=1)=q/2$.

Tag a flag whose producer has rank two. Its signed bias is $1/4+\epsilon/2$ on the ideal law of §4. Correct each receiving marginal to this value by moving mass across that bit alone, preserving the other marginal. The TV cost equals the marginal change. Tags have unique consumers, so total TV is $o(n)$. At rank one keep the independent fair innovation. Keep schedule choices fixed during the bounded-cost transfer. For class $i$ let $R_{g,i}$ count incoming tags minus outgoing tags. This balances exactly.

For a sole tag put it first, also respecting the zero-flag-parent order at $T$. At rank zero with one fixed bias $q$, use the convex hull of $(1-q,q,0,0),(1-q,0,0,q),(0,q,1-q,0),(0,0,1-q,q)$. Divide the rectangle into two triangles. With two tags the earlier product is in the later basis past. Replacing later innovations makes its product independent of that past. Correct the earlier marginal last, at aggregate TV $o(n)$ by unique uses. Use $(1-q-r+qr,q-qr,r-qr,qr)$. Test only matching signs. Rank one uses the fixed point $q$. The zero-tag domain is unchanged.

## 5. Flags and affine forcing

Write $Z=UV$ for a nonlinear gate on its two signed operands. Its native output is $Z$ or $1+Z$. Flag an operand from a nonlinear producer with one native use. Classify the slot by its consumer: $A$ means one flag and rank zero, $B$ means one flag and rank one, and $T$ means two flags. Rank is the raw rank increment $j$ from §4. A flag has zero formal raw row, so these classes include every slot. Write $b_g$ for the number of flags at gate $g$, $i_g$ for their common class, and $o_g$ for the class of the flagged slot fed by $g$. Use class zero when the corresponding flags are absent.

Put $\psi(q)=\min(q,1-q)$. For each $i\in\{A,B,T\}$, define
\[
 M_{g,i}=b_g[i_g=i]-[o_g=i],\qquad
 B_{g,i}=[i_g=i]\sum_{F\text{ flagged at }g}\psi(\Pr(F=1))
              -[o_g=i]\psi(\Pr(Z_g=1)).
\]
Matching slots cancels both counts and biases on true laws by $\psi(q)=\psi(1-q)$. Section 4 gives ideal bias sums $o(n)$ on good modes.

For flag sign $\epsilon$, each class balances incoming $(\epsilon,\Pr(F=1),\Pr(F+\epsilon=1))$ against outgoing $(\epsilon,\epsilon+(1-2\epsilon)z,z)$ at $z=\Pr(Z=1)$. Denote the differences by $U_{g,i}$. Matching slots cancels them on the true law. The ideal errors sum to $o(n)$.

### Affine certificates on the flagged forest

Link each flagged producer to its sole successor. Nonflagged affine operands and final output depend on raws and nonlinear roots only. At an accepting $x$ substitute its root values in these forms. For a zero-flag gate fix both one operands, or any zero operand. For one flag $F$ from $p$ with other operand $W$, fix $W=0$ if zero. Otherwise certify $p$, fixing $W=1$ only if $F=1$. For two flags certify both ones or the cheaper zero. Let $N_g$ count equations. Use every root and the final accepting equation. Native induction gives correct roots and requested flags, since omitted parents have sole use at skipped or masked gates. This accepting affine fiber has dimension below $2k$, so $\sum_{\rm roots}N_g\ge n-2k-1$.

Set $L_g=1+Z_g,F_g,0$ for zero, one and two flags respectively, and $G_g=L_g+\sum_{p\text{ flagged at }g}N_p-N_g$. Counting each flagged producer once gives $\sum_{\rm roots}N_g=\sum_g(L_g-G_g)$. Induction yields $N_g\ge1+Z_g$. At a one-flag gate
\[
 N_g=1+Z_g+W_g(N_p-1),\qquad G_g=(1-W_g)(N_p+F_g-1).
\]
Write $F=Z_p+\epsilon$. At $A$ use $(G_g,N_g)\ge(2U(1-V),1+2UV)$ for sign zero and $(1-V,1+UV+(1-U)V)$ for sign one. At $T$, signs $\epsilon,\delta$ give the lower vectors
\[
 (1+\max(\epsilon,\delta),2-\epsilon,2-\delta,0),\qquad
 (1+\min(\epsilon,\delta),1+\delta,1+\epsilon,4-\epsilon-\delta)
\]
on $00,10,01,11$ respectively. Use $N_p\ge1+Z_p$. On $00$ omit the larger certificate, on a mixed state omit the one flag, and on $11$ keep both. Let $A_g,T_g$ and $K_g$ be the expected gain and size bounds at these classes.

### Three-step class-B propagation

Let $h_g$ count consecutive $B$ consumers on the outgoing flagged path, capped at three. Codes $2,4,5$ denote horizons $1,2,3$ and retain class $B$. A $B$ consumer has incoming horizon $h_g^{\rm in}=\min(1+h_g,3)$, and all other gates have incoming horizon zero. For $d=2,3$ put
\[
 D_{g,d}=[h_g^{\rm in}=d]-[h_g=d],\qquad
 B_{g,d}=[h_g^{\rm in}=d]\psi(\Pr(F_g=1))-[h_g=d]\psi(\Pr(Z_g=1)).
\]
These balances vanish on true laws. Index $D$ adds the clipped balance for horizons at least two. Counts $t_g$ and $u_g$ match zero-flag parents into $T$ and $A$ with the corresponding outgoing slots.

At $B$ substitute $N_g=1+W_g(F_g-1)+W_gN_{p(g)}$ into $G_v$ for three predecessor steps, dropping the nonnegative tail. Reverse §4 replacement averages factors $W,1-W$ to $1/2$ and the base to $(1+F)/2$, using earlier-basis measurability. Unique successors and bounded depth give error $o(n)$. The base weight is $1-2^{-h_g}$, with lower mean $K_g$ as above at $A,T$, $1+\Pr(Z_g=1)$ at zero flags and $(1+\Pr(F_g=1))/2$ at $B$. Each $B$ consumer contributes the remaining $-\tfrac12(1-\Pr(F=1))$.

With $\ell_g=E[L_g]$ and $\sum j_g\le n$, the equation bound now gives $\sum_g C_g^*\ge-o(n)$ for
\[
 C_g^*=\ell_g-j_g-(1-2^{-h_g})K_g+
       \tfrac12[i_g=B](1-\Pr(F_g=1))-[i_g=A]A_g-[i_g=T]T_g .
\]

## 6. Query capacity and physical entropy ledger

Fix the mode and selectors. Put $a_0=0,a_A=a_B=1,a_T=1/2$. At one flag request it, and at two choose one uniformly. Draw the requested flag conditional on its producer's earlier query $Q_p$, costing $H(Z_p\mid Q_p)$. The fixed flag makes the product affine. Its later conditional charge is $a_o\operatorname{mean}_F H(Z\mid F)$.

At zero flags and $o=0$ draw a resolving pair. The incidence polytope is $v_i\in[0,1],\sum v_i=2$, since its vertices are integral. Use $v_i=2p_i$ unless $r=p_i>1/2$, then set $v_i=1$ and $v_j=p_j/(1-r)$ otherwise, arbitrary at null residual. This costs $H(U,V)-1$ or $H(U,V)-H_2(r)$.

For $o=A,B$, draw $11$ with probability $z=p_{11}$ or a pair from the zero corners using the relative incidence rule. The cost is $H(U,V)-(1-z)$, or $H(U,V)-h(r,1-z-r)$ when a zero mass $r>(1-z)/2$. Here $h(a,b)=\phi(a)+\phi(b)-\phi(a+b)$, $\phi(t)=-t\log_2t$, $\phi(0)=0$. The inferred product supplies the requested flag.

Mark every zero-flag parent of $T$ and delay it to its sole use. Query a uniform parent, using a constant certificate when marked. Charge half that cost per marked parent. When the chosen flag is one, resolve an unchosen marked parent. An unmarked parent is already affine. On zero the product vanishes. The saved query remains the ordinary query for subsequent groups.

For a resolving draw use a uniform pair containing $11$ with probability $a=i/8$, $1\le i\le8$, or a uniform zero pair. Incidences are $a$ and $v=(2-a)/3$. Relative to one bit the expected gain is $q\log_2(2v)+w\log_2(a/v)$ for $q=\Pr(F_{\rm chosen}=1)$ and $w=\Pr(F_{\rm chosen}=1,Z_{\rm unchosen}=1)$. Use mrates' directed affine floors and average query choices. Put the marked parent first for $t_g=1$, so the baseline resolution charge at $t_g>0$ is $(q_V+[t_g=2]q_U)/2$ for flag-one probabilities $q_U,q_V$. The normalized $Z$ in $w$ accounts for the unchosen flag sign.

Group a zero-flag producer with up to three initial $B$ gates. Delay sole-use internal gates to the last. Other operands use raws and earlier roots. Reverse §4 replacement gives aggregate TV $o(n)$. Set $Z_0=UV$, $F_i=Z_{i-1}+\epsilon_i$, $Z_i=F_iW_i$, $z=p_{11}$, $q_i=\Pr(F_i=1)$. Debit $\sum_{i<r}q_i+a q_r$ for $r$ successors and last request $a$, retaining their old charges. A requested certificate fixes the product and returns the query. An unused result needs affinity.

With positive signs before the last, the product on $d=r+2$ bits is a singleton indicator plus an affine term. For no request, draw a uniform line through the exception with probability $1/(2^d-2)$ or a uniform avoiding hyperplane. For a positive request use the point with probability $2^{-d}$ or those hyperplanes. Nonexceptional incidence is $1/2$ and exception mass $z/2^r$. For complemented request at $r\ge2$, query the last bit, then fix the preceding product on one. For complement at $r=1$ draw $W=0$, a uniform line from $111$ to it, or a uniform pair of other $W=1$ corners with weights $11/24,1/6,3/8$. Respective incidences are $1/2,1/6,1/4$.

At $r=2$ with complemented first sign, for no request draw a uniform line in $W=X=1$ with probability $1/4$, or a uniform $(W,X)$ hyperplane avoiding $(1,1)$. Affine incidences are $1/8,1/2$ inside and outside. If requested with complemented second sign, query $X$ then use the complemented $r=1$ certificate on one. For positive second sign use $u=(U,V)+(1,1),v=(W,X)+(1,1)$ and product $[u\ne0,v=0]$. Uniform hyperplanes $\ell(v)=1$, linear graphs $u=Mv$, or nonzero pairs at $v=0$ get weights $57/80,1/10,3/16$. Incidences are $1/2,1/8,1/10$ at $v\ne0$, $u\ne0,v=0$, and $(0,0)$. At $r=3$ stop with the request at two when either of the first two signs is complemented.

On $z\in[1/8,3/8]$, the costs after debits are bounded by $1+4z/5,1-19z/64,1-z$ at horizons $1,2,\ge3$. For a complemented first sign use $1/2+z$ at horizon two and $1/5+21z/20$ at higher horizons. The respective two-successor request bounds are $19z/12,1/4+13z/20$. The logarithm bounds follow from $6^5<2^{13},14^{16}<2^{61},10^3<2^{10},30<32$. Use the stated band at rank zero, $z=1/4$ at rank two, and $z=q/2$ at rank one with $q\in[1/2,5/8]$ at horizon one and $[3/8,3/4]$ otherwise. Use the old charge outside. Select by ideal mode law and keep choices fixed through the bounded-cost transfer.

For flagged rank zero at horizon at least two, keep incoming query $F$ and group with first $B$, whose output is requested. Costs agree on $F=0$. On $F=1$ the new bit is fair and the signed operand has bias $r$. The old cost is $H_2(r)+r$. The zero-pair cost is $H_2(r)+r/2$ for $r\le2/3$, else $H_2(r)+1-h(r/2,1-r)$. The entropy chord gives rebate $\min(r/2,1-r)$ or more. Weight by $\Pr(F=1)$ for $\min(z/2,p_{10})$ on positive first sign or $\min(p_{10}/2,z)$ on complement. At $T$ average the $p_{01}$ version.

For complemented first sign and horizon at least three, instead keep ordinary queries and pair the first two $B$ gates. Given incoming flag one, fix the fair product for $5/4$. Given zero query the final fresh bit only if the internal sign $e=1$. For $q=1-z$, new cost $5q/4+e(1-q)$ saves $q/4$ from old $q+q/2$ or $q+1-q/2$. Keep the next inferred query. Paths have distinct non-$B$ starts and use §4's reverse replacement.

Call the charge $\widetilde J$. For $J^R$ keep groups and average both parities of adjacent pairs on each unused suffix. Each pair has a following $B$, giving $q/8$ per eligible start. Debit removed starts at their non-$B$ origin at horizon at least three. At flagged starts debit $q_1/8$ for positive sign. For complement group the first two only on ordinary query $F_1=1$. Costs agree on zero and $F_2$ is fair on one for either sign, so debit $3q_1/16$. Choose parity independently.

For a zero-flag start with positive first sign the total is at most $1-z+7z/32$. A complemented second sign stops at two, costing $1/2+(1+z/2)/8$. Three positive signs lose at most $7z/4$ before dividing by eight. A complemented third sign with request costs $1/2-7z/8+(1+5z/4)/8$. Without request the path ends there and only the first $B$ can start a pair, with bias $z$.

For a complemented first sign, the two next signs give net charges $19z/12,1/4+13z/20$ and lost sums at most $3(1-z)/2,1+(1-z)/2$. Add one eighth of each loss. The band endpoints certify the upper affine bound $1/5+21z/20+27/128-z/4$.

For $J^G$ omit complemented flagged prefixes. Scan each unused $B$ suffix for starts followed by two $B$ gates. Group three if the ordinary incoming flag is one and no previous group uses the start. Absorb the first internal sign into the fair first bit. The next two flag biases are $1/2$ and $1/4$ or $3/4$. On request the three-bit constant certificates cost $5/4$ or $5/4+\log_2(6)/8$, saving at least $1/2$ or $43/64$ from the old three charges. Without request the singleton-affine cover costs $7/8+\log_2(6)/8$, saving at least $19/64$. Only the last eligible start can lack a request.

Each eligible one-flag start is chosen or blocked by one earlier group. A choice pays for itself and at most $1/2+1/4$ or $1/2+3/4$ expected later ones, or only itself at the last start. Its saving is at least $2/7$ times that count. The choice depends on the earlier basis, and reverse §4 replacement transfers costs and counts with bounded reuse. Credit $2q/7$ per eligible start, and return any inferred request.

Debit $2z/7$ at positive flagged prefixes and undo $q/4$ at complemented ones. At grouped zero-flag prefixes multiply the same missed sums by $2/7$. Their affine costs at both band ends give added debits $17z/28-3/224$ or $67/140-19z/35$ by first sign. The three ledgers hold separately.

Fixing the sample leaves an affine accepting fiber of size at most $2^{2k}$. On matching prefixes all queries are true. Multiply answer and incidence probabilities, including the marked resolving incidence. Thus $\sum_{x\in S_\sigma}p(x)\le2^{2k}$. Jensen gives $E[-\log_2p(X)]\ge H(X\mid\sigma)-2k$. Average choices and apply §4 to the continuous charges.

The balances $Q_{g,i}$ count $t_g$ at $T$ and minus one at each marked parent, also weighting these counts by the matched flags' $\psi$. Complement invariance gives true sum zero and ideal sum $o(n)$.

Let $S_0$ sum physical credits $D'(Y)=1-H(Y)$ and retained $I'=I(U;V)$, in $[0,1]$. Transfer $\eta_iD'$ from producer to class-$i$ flag, canceling by complement invariance. Then $c_z=1-\eta_o$ and $c_u,c_v=\eta_i$ on flags, one otherwise. Gate degree $d\ge2$ has $2d-1\ge d+1$ credits for inputs and output. One use shares with its flag. Affine outputs need no reservation. Raw degree one or $d\ge3$ covers its uses, with one possible shortfall at two. By §§3–4 the losses are $o(n)$, so $S_0\ge\sum_g E_g-o(n)$ for $E_g=c_zD'(Z)+c_uD'(U)+c_vD'(V)+I'$.

Put $\Delta_g=\min(I(U;V),H(U,V\mid Z))$. Section 7 earns it on exposed gates, with loss $o(n)$ elsewhere. The checker proves
\[
 E_g+\Delta_g\ge \beta j_g+\theta_\ast C^\ast_g+(\kappa-\rho-\rho_G)\widetilde J_g+\rho J^R_g+\rho_G J^G_g
  +\sum_{i\in\{A,B,T\}}(\xi_iM_{g,i}+\gamma_i B_{g,i})
  +\zeta_T(t_g-[b_g=0,o_g=T])+\zeta_A(u_g-[b_g=0,o_g=A])+\sum_{d=2,3}(\xi_d D_{g,d}+\gamma_d B_{g,d})+\sum_{h,i}\gamma_{h,i}B^{(h)}_{g,i}+\sum_i(\alpha_i Q_{g,i}+\nu_i\cdot U_{g,i}+\chi_i R_{g,i})+\mu .
\]
Weights in gain.py use denominator $10^8$. $(\beta,\theta_*,\kappa,\lambda,\mu)=(5223940,2014996,8305050,8305050,1)/10^8$. $(\rho,\rho_G)=(0,8305050)/10^8$. The one-sided multipliers are nonnegative. For $h=6/64,9/64,\ldots,30/64$, form $B^{(h)}_{g,i}$ by using $\min(q,1-q,h)$ on both sides of the flag for $i=A,B,T,D$. Complement invariance gives exact balance before ideal replacement, hence aggregate $o(n)$ afterward.

Use corners $00,10,01,11$ and the first operand for a sole flag. For rank-one zero-flag cases put $J_c=1+H_2(q)-h(q/2,1-q)$, from a zero corner of mass $q/2$ with incidence one. Use $H_2(q),J_c,J_c/2$ at $0,A/B,T$, with the stated group bounds. At rank one with a flag use $J=a_oq$. At rank two the unmarked bounds and paired rebates above apply, and a marked parent uses $5/8$ in the marked schedule.

The checker tests both branches of $\Delta$ at rank zero, drops $\Delta$ at rank one, and checks the fair point at rank two. Use affine corner masses for $\phi$, combining complementary terms of the same sign into $H_2$. Symmetry halves unpartitioned symmetric cases. The resolving domains use four corner tetrahedra and a triangulated octahedron. Constant-product domains cone four zero-face triangles to $11$. At each chain cut $z=1/8,3/8$ triangulate the slab between corresponding ordered faces as $(l_0,\ldots,l_i,u_i,\ldots,u_2)$. Use the stated rank-one intervals and all $T$ parent counts. For $t_g=1$ retain the full simplex. Choose one marked affine gain per sign and branch using the simplex center, then apply it to every vertex. Cache choices by dyadic cell using exact ordered line intersections.

Bound positive entropy terms by chords, negative terms by upper tangents, and clipped biases by the corresponding bounds. The vertex minimum over signed corrections and concave rebates gives a common lower chord. Check vertices with directed rounding, bisecting a longest edge on failure. Distance includes rank-zero marginals and product. Depth 70 preserves exact dyadic masks.

For $x=2^e u$, $1\le u<2$, set $w=(u-1)/(u+1)$. Use twelve terms of $2\sum_{j\ge0}w^{2j+1}/(2j+1)$ and tail $2w^{25}/[25(1-w^2)]$, also at $u=2$. Bound $\log_2 x$ and $\phi'(x)$ outward. Strict finite comparisons allow a fixed reserve increase.

Summing the local inequalities and using the rank, query, forcing and balance bounds gives $S_0+\sum_g\Delta_g\ge s_0n-o(n)$, with $s_0=\beta+\kappa+\mu=13528991/100000000$. The ideal-law replacement affects $\Delta_g$ by aggregate $o(n)$, by its continuity. At least $n-o(n)$ nonlinear gates contribute the reserve.

### Contract private-input leaves and output stems

Mark each single-private recipient with retained other input. By the first two bundle exclusions in §3, $o(n)$ private raws remain unmarked. For a marked $g$, use private bit $X_g$, other operand $V_g$, output $Y_g$, and outgoing stem including final use. At an accepting word there are fewer than $2k$ masked nonlinear recipients. Varying their private raws keeps all outputs fixed by native induction, and a split of $2k$ such directions violates §A.4. Thus their constant-case probabilities sum to at most $2k$ in every mode.

Group recipients by earliest consumer, dummy last. Each private group $A_i$ has at most two bits. Let $P_i$ contain unmarked raws and preceding groups. It determines $V_g$, since any marked output reaching $V_g$ has an earlier use. A designated marked stem has a sole consumer with no private raw by priority. For any such consumer, its other operand uses only earlier marked groups or the same group's direct second input, since later groups have no earlier uses and pins differ. Thus it is known given $E_g=(P_i,A_i\setminus\{X_g\},\sigma)$. This also holds for an undesignated sole consumer without private raws.

The chain rule gives $1-H(X_g\mid E_g)\le |A_i|-H(A_i\mid P_i,\sigma)$. Sum at most twice the marked-raw deficit given unmarked raws, hence $o(n)$ by §4. Leaf credit is at most $1-H(X_g\mid E_g)$. Stem credit is at most $1-H(Y_g\mid E_g)$, using $D'+I'=1-H(Y_g\mid\text{other operand})$ for designated stems. A signed identity except on masked cases gives
\[
 1-H(Y_g\mid E_g)\le1-H(X_g\mid E_g)+\Pr(g\text{ constant in }X_g\mid\sigma).
\]
Thus leaf and stem losses are $o(n)$, including bounded unmarked leaves.

Absorb leaves and select each marked degree-two bag's outgoing stem. Native order gives chains, at most three per unmarked root and one incoming selected edge per marked bag. Break into mass-$[B/8,B/4]$ pieces with the final residual merged into its root. Leave $O(n/B)$ stems. Each degree-two end merge preserves rank and reduces deficiency by one. The dummy, marked bags, and $o(n)$ exceptions give the initial deficiency. Degree sum and connectedness imply $E\le3r+o(n)+O(n/B)$ and the credit bound gives $S\ge S_0-o(n)$. Native evaluation order holds within contracted bags.

## 7. Causal weighted transcripts and strict contradiction

### Exposed gates

An exposed gate stands alone with two retained arms and a designation. A later owner of a degree-two raw $X$ has an earlier marked recipient $h$ with private raw $Y$, up to $o(n)$ exceptions. Except for $o(n)$ raws, mode bias lies in $[1/4,3/4]$. Nonlinear $h$ is then constant in $Y$ with probability at least $1/4$, so §6 allows $o(n)$ such gates. At affine $h$, use the invertible raw chart $(X,Y)\mapsto(X,H)$, $H=X+Y+\epsilon$, on disjoint pairs. Supply $H$ at $h$, moving first use of $X$ to its later owner.

Group chart coordinates by first use into sets $A_i$ of at most two, with earlier coordinates $P_i$. An absorbed $X$'s other operand is known given $P_i,A_i\setminus\{X\}$. Thus $I(X;V)\le1-H(X\mid P_i,A_i\setminus\{X\})\le |A_i|-H(A_i\mid P_i)$. Sum at most twice the chart deficit. For private raws and stems into sole consumers without private raws, §6 bounds lost information by the $E_p$ deficits even without designation. Stems into forks lose no operand arm. Missing designations cost $o(n)$ by §3.

### Capacity per direction

At an exposed gate use $j=H(U,V)$, $a=H(U\mid V)$ for its designated $U$, $b=H(V\mid U)$, $d=H(Z)$, $D=\max(a,b,d)$, and $B=\max(b,D-a)$. For either direction let $A_U,A_V$ indicate incoming crossings. Replace that gate's two input charges by

\[
 F=(D-B)A_U+(D-a)A_V+(a+B-D)\mathbf 1[A_U\lor A_V].
\]

Coefficients are nonnegative. For the direction with receiver bit one, $A_U=z(1-x)$, $A_V=z(1-y)$ and their union is $z(1-xy)$ for supplier bits $x,y$, receiver bit $z$. Their mixed differences are nonpositive. Complement all bits for the other direction. Identifying parallel suppliers preserves submodularity. Costs $0,a,B,D$ replace $0,a,H(V),j$. Here $H(V)-B=\min(I,j-D)\ge\Delta_g$ by $j-D=\min(H(U),H(V),H(U,V\mid Z))$.

Each oriented tree crossing has probability at most $p=c/3$. At an exposed star, $Pr(A_U\lor A_V)\le2p$ gives cost at most $p(a+B)$. Every other arm uses its nonnegative mode charge with factor at most $p$. Thus each lifted directional sum has mean at most $p(E-S-\sum_{\rm exposed}\Delta_g)$. Apply §1 to both sums and the raw mass with the mode charges fixed. The same cut gives both mean bounds plus $o(n)$ and both raw counts $n/2+o(n)$.

### The exact simulation

At nonexposed gates use the carrier plan with designated guesses and true bits for equality. At an exposed gate guess a single crossing after obtaining the local operand. For two crossings compute and send the product at its native step on the supplier side. Costs are at most $a,B,D$, since $d\le D$. Route uses from the original output slot.

Exchange and compare both guess-equality shares. Native induction forces truth, hence at most one passing guess pair per raw pair. Check the pins with at most two producer values per indicator. Check the output and both costs.

Charge each receive or guess its negative binary log probability in the mode, conditioning guesses on the other operand. Mix with fixed fair weight $\omega>0$, using fair laws on null events. The mean excess is at most $-\log_2(1-\omega)$ per bit. Costs are bounded functions of §4's list, including moved outputs. Choose $\omega$ small for fixed $\epsilon>0$. Section 4 bounds realized costs by $W_{\rm dir}+\epsilon n$ outside $o(1)$ accepting mass.

Passing transcripts are Cartesian in raws and guesses. Raw projection on either side is injective by uniqueness after fixing the other side. Each projected accepting rectangle has a side smaller than $K$.

Fix receiver raws. Its receives and guesses specify its trajectory and sends, in a tree of normalized smoothed laws even on invalid histories. Passing paths have probability at least $2^{-W_{\rm dir}-\epsilon n}$. Local equality fixes the sender share. Pins and tests add at most $2s+O(1)$ bits, so there are at most $2^{W_{\rm dir}+\epsilon n+2s+O(1)}$ passing leaves per raw word.

Suppose $t\le(L+d_0-\tau)n$ with fixed $\tau>0$ on a subsequence. Use $S+\sum_{\rm exposed}\Delta_g\ge s_0n-o(n)$, $E\le3(t-n+1)+o(n)+O(n/B)$ and $n/2=c(L-1)n$. For large fixed $B$, both $q_{\rm dir}=n/2+o(n)$ exceed $W_{\rm dir}$ by fixed $\gamma n>0$. Take $\epsilon<\gamma/2$. Count each leaf towards its thin sender, with fewer than $K$ words per opposite word. The at most $2^s$ modes cover cube mass $2^{-\Omega(n)}$ since $s+\log_2K=o(n)$. Dividing by $2^{-D}$ contradicts accepting pass mass $1-o(1)$.

The actual input coefficient is
    L+d0 = 1+1/(2*(12233148667589113194/10^20))+13528991/300000000
         = 3139241126417034843716601209/611657433379455659700000000.
As n/N tends to one, a fixed reserve increase from the strict 42,470 and 180,348 cell checks absorbs ignored inputs. Both phases check the shared prefix digest. Terminal adds the second path and directed energy.

## Provenance


Declared sources include V31 overlap and V30 graph bounds.

## One-edge repair and marked covariances

Thirty disjoint six-prefix rows at 11 force a positive C edge of ordered degrees (1,2) and zero matching bits for all exterior trits. M stays one iff more ends have C degree two than zero. Twice the oriented mass gives ≥145636395/10^12 per vertex.

On a marked path (w_i,b_i)=(t,a²) straight or (a,at) turning. Covariance is (a+t/2) prod_i w_i + (1/2) sum_i b_i prod_{j!=i} w_j by the independent off-path modes. Diagonal one. Source: Yukon 086221d5-56f1-4068-9241-b118a66143fc.

## Adjacent-center repairs and longer tails

For degrees d_v, the saving for deletion set D is the sum of 4−d_u−d_v over D plus one per adjacent pair in D, by expanding (3d−d²)/2. Count each singleton and each certified two-C-edge center with nonnegative flanks once.

The 184 disjoint center cylinders check 11,813 trit cases, covariances and triangular charts. Directed floors are 33885419485/10^15 for the first 48 and 14098240991/10^15 for the next 136.

Delete a certified run of k=1,2,3 positive weight −1 C edges with both positive nonnegative flanks. The k negative costs plus k+1 pair terms yield one. Charge sets are disjoint even with shared flanks. At a fixed start, lengths are exclusive and 17 cylinders pass same-length disjointness. Recognize either orientation once. Terminal verifies the cold C-line run floor of 9479849696849877/2^72 per vertex.

For a positive M edge with degrees 3 and 2, require its degree-2 end's unique positive C neighbor to have degree 2. At the degree-3 end certify each positive C neighbor of degree 1. Flanks have weight zero. Delete the M edge and these flanks together. The low pair pays weight −1; each high pair gains one. Group by M edge, so simultaneous high-side witnesses earn distinct pairs. These M-containing charges are disjoint from C charges. Six disjoint cylinders per oriented choice give degrees (1,3,2,2). Eight end/direction choices and M density 1/2 multiply mass by four, giving 3100375520549915/2^69 per vertex at the cold threshold.

Activate repairs at common p and delete their union. Each operation gains one or more and touches at most six vertices, so the gain covers one sixth of the majority loss. The broad motif saving and all four cold floors are added exactly in the strict cap computation. Their sum B is below .001. The directed initial majority exceeds .504 even after the extra motif debit. If activation crosses half, use the pre-repair margin above .568, so the gain exceeds .068/6>B. Otherwise finish with terminal patches. Deletions preserve terminal domination. Finite approximations preserve the strict margins.

The simplex-root overlap bound follows Yukon b2cc7bbd-6e3e-4054-ad13-9fa30a1bb5ae. The fractional, first-pass and singleton terms follow yukonautominer 171fcf07-f20e-4f0b-858e-c400c4de21c3.
