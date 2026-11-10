# Finite support graph,linear sampler,and conditional breaking
Fix Q_S=213938373412081/100000000000000. With L=floor(log2(actual input length)),S=floor(L/Q_S)
is an integer;S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer rounding.
All logarithms are binary,addition is XOR unless a field is specified,and
statistical distance is TV = half the L1 distance. A k-source is a finite
random variable with every atom at most 2^(-k);k may be real. Graphs have
labelled edges and may have loops or repeated edges. Every constant is fixed
independently of the eventual outer input length N.
The support-growth Lemma S1' (§SP.6.2),the affine-only
three-condenser/kernel base of linear Section7f and the affine breaker yield the finite
graph/sampler and conditional-breaking reductions below. proof.md supplies
the tagged affine-permutation outer assembly and finite guards.
## 1. Supplied contracts and elementary probability facts
### 1.1 Selected seeded-linear contracts
The supplied general-source maps are PXG(v) from linear7g and
PXTM(j,v) from linear7h. Each fixed ordered rate/arity vector defines
an exact integer threshold K_0(u,z,e),complete seed length D(u,z,e),
and deterministic map Ext:{0,1}^u x {0,1}^D -> {0,1}^z,source-linear
for every fixed seed. On every legal call,for every source X with
H_inf(X)>=K_0 and an independent uniform complete seed Y,
TV((Ext(X,Y),Y),(U_z,Y))<=2^-e. (1)
All ignored suffixes are retained;a shorter output takes a prefix.
The literal algorithms,threshold and seed slices include all ceilings.
Before field setup require K_0<=u. PXTM hasj_i>=2,R_i>j_i-1 and
w>=max(2^20,e+lg(product j_i)+6+2lg(4u)). Section7h proves its
source bound K_0<=Lam*z+g*w. BOTH PXG and PXTM use the clean seed
bound D<=z/B+[4+sum_i(R_i+1)]w only under the additional output-cap,
terminal-rate/arity/kap and height guards of7g. Closed slots and B9
check these hypotheses for every selected call.
Ext0 and Ext2 use their respective PXTM vectors;inner and outer PXG
maps remain binary. Canonical condenser fields are supplied in7a;leaf matrices in7g.
Enumerated base degrees are O(log actual N) at selected calls,so exhaustive base
field setup costs N^O(1);larger relative extensions use polynomial-cost
towers. Linear5/proof.md prove uniform polynomial time in actual N,
without asserting polynomial exact-degree search for arbitrary u.
The selected BASE instead uses linear7f's affine-only three-condenser/
kernel algorithm with its complete seed and integer threshold. That
affine contract alone is consumed in Sections3--5;it is not substituted
for the general-source contracts needed inside the breakers.
### 1.2 All-seeds-retained affine-breaker contract
Let F:{0,1}^n x {0,1}^d x {0,1}^a -> {0,1}^m be one fixed
deterministic map,with a>=1. Its supplied degree-r contract at source
entropy K and error eta means the following universal statement. For
arbitrary correlated tuples
Left=(A,A_1,...,A_r),
Right=(B,B_1,...,B_r,Y,Y_1,...,Y_r),
assume Left and Right are independent,H_inf(A)>=K,Y is uniform on d
bits,and fixed advice strings alpha_i all differ from alpha. Set
V=F(A+B,Y,alpha),V_i=F(A_i+B_i,Y_i,alpha_i).
Then,with U_m independent of all displayed side information,
TV((V,V_1,...,V_r,Y,Y_1,...,Y_r),
(U_m,V_1,...,V_r,Y,Y_1,...,Y_r)) <= eta. (4)
Only A needs entropy and only Y uniformity. No A_i need be a function
of A. All within-side correlations,repetitions,and seed equalities are allowed. Tampered advice
strings need not differ from one another. For a smaller degree,pad with
repeated permissible copies and discard their outputs/seeds. This
branch uses the all-seeds-retained affine contract.
### 1.3 Probability identities
For any common marginal Z,grouping the TV sum by z gives
TV((V,Z),(V',Z))=sum_z P(Z=z) TV(V|z,V'|z). (5)
A common deterministic map or random conditional kernel cannot increase TV:
sum differences over each fiber and use the triangle inequality. In all
ideal laws below a displayed U is fresh uniform independent of the entire
retained side information.
Define average conditional min-entropy by
Havg_inf(A|Z)=-log sum_z P(Z=z) max_x P(A=x|Z=z).
If Z has at most 2^ell possible values,then
sum_z max_x P(A=x,Z=z) <= 2^ell max_x P(A=x),
Havg_inf(A|Z) >= H_inf(A)-ell. (6)
If Havg_inf(A|Z)>=K+log(1/delta),with 0<delta<=1,the total probability
of z with H_inf(A|Z=z)<K is at most delta: on that set the maximum
conditional atom exceeds 2^-K,while its expectation is at most
delta*2^-K. These facts require no unstated independence.
## 2. An explicit small-bias Cayley graph
Fix a field dimension r_F>=1 and put c=2r_F;Section3 chooses r_F=S,s=2ceil((50/7-2^-100)*S/2).
When r_F is a power of two,one admissible choice is the canonical tower field of binary dimension r_F
(§3 does not use it). Here is its complete
finite choice. Start with F_2 and trace-one element a_F=1. Repeatedly
adjoin w satisfying w^2+w+a_F=0 and replace a_F by a_F*w. The image of
x -> x^2+x is exactly the trace-zero subspace: its kernel has size two,
each image has trace zero,and the trace is nonzero because the current
a_F has trace one. Thus the quadratic has no root and is irreducible.
The two conjugates of w are w,w+1,so the relative trace of a_F*w is
a_F,and trace transitivity proves absolute trace one at the next level.
Trace transitivity follows by grouping the Frobenius-power sum. This
constructs every power-of-two dimension with a fixed ordered binary basis,
without an irreducible-polynomial search.
In one tower step write x=x0+x1*w,y=y0+y1*w. With
p=x0*y0,q=x1*y1,r=(x0+x1)(y0+y1),t=a_F*q,their product is
(p+t)+(r+p)*w. Four recursive products and linear extra work give O(r_F^2)
bit operations for multiplication in the final field.
Section3 uses r_F=S,s=2ceil((50/7-2^-100)*S/2) and,for every r_F (also powers of two),the canonical field of
weak-seed §4 of linear-and-bilinear-primitives.md (the first irreducible binary polynomial of degree
r_F in lexicographic order,found in 2^(r_F)poly(r_F) bit operations) with its
polynomial basis;schoolbook multiplication costs O(r_F^2) bit operations.
Everything below uses only field arithmetic and a fixed binary basis.
For each m>=1 and each label (u,w) in F^2,define an m-bit increment
g(u,w)_j=<w,u^j>,j=0,...,m-1,(8)
where brackets are the ordinary coordinate dot product and u^0=1 also
when u=0. On vertices F_2^m use the 2^c labelled edges x -> x+g(u,w).
A fixed label is an involution,so the transition matrix P is symmetric
and doubly stochastic,with the uniform stationary distribution.
For any nonzero a in F_2^m,the character (-1)^<a,x> has eigenvalue
E_(u,w)(-1)^<w,p_a(u)>,where p_a(T)=sum_j a_j*T^j. The polynomial is
nonzero of degree at most m-1. Averaging over w annihilates each nonzero
coordinate vector p_a(u),leaving its fraction of zeros. A nonzero degree-j
polynomial has at most j roots,by repeated division by T-alpha. Hence
0<=eigenvalue<=(m-1)/2^r_F. (9)
The binary characters form an orthogonal basis,so ||P-J||<=lambda=(m-1)/2^r_F,
where J projects onto constants in the uniform-normalized Euclidean inner
product. A labelled step is evaluated by successive field powers and dot
products in O(m*r_F^2) bit operations.
## 3. The affine-base walk
Fix Q_S as fixed above,beta=50/7-2^-100,q_walk=17,C=4(q_walk+1)=72 and
p=144>3(beta+2q_walk),a retained upper budget. For actual n let L=floor(log2 n),
S=floor(L/Q_S)>=2^160,s=2ceil(beta*S/2)<=8S,ellc=lg(s+4),
R(s)=0,
r_F=S and m=s+2q_walk*r_F<=48S+320ellc.
The integers s,m are even. Use Section2's graph on F_2^s with field
degree r_F. An m-bit string specifies an s-bit start and q_walk complete,
independent 2r_F-bit labels. Its positions are W(u,j),j=0,...,q_walk.
The transition matrix is symmetric and ||P-J||<=lambda=(s-1)2^-r_F.
Lemma W_A. If |U|>2^(9m/16-8) and N subset F_2^s has
|N|<=5*2^(-r_F)2^s,some walk position of some u in U lies outside N.
Proof. Put rho_N=|N|/2^s and let Pi project coordinates onto N.
Then ||Pi P Pi||<=rho_N+lambda. A uniform walk stays in N with
probability at most rho_N*(rho_N+lambda)^q_walk. Thus the number
of all-bad m-bit strings is at most
5*2^(m-r_F)*((s+4)2^-r_F)^q_walk.
Its logarithm minus9m/16-8 equals
(7/16)s+log2(5)+8-r_F-q_walk*(r_F/8-log2(s+4))
<=-(7/16)*2^-100*S+100+20log2 S<0.
Here s<=(50/7-2^-100)S+2 and s+4<=16S.
For 2^j<=S<2^(j+1),j>=160,the negative term in magnitude is
at least (7/16)*2^(j-100),while the positive part is at most100+20(j+1).
At j=160 the first exceeds3320;it doubles at each step while the
second increases by20 and is positive. Induction proves strict negativity
at EVERY allowed S. U cannot consist only
of all-bad strings. This proves W_A at every allowed integer S.
## 4. The support graph and property (G_A)
Set s_tag=3m<=144S+960ellc,A_rows=2^s_tag,A0=2^s,
Lambda_S=2^(3m/2-17),epsilon0=2^-r_F. A_rows is the number of rows;
it is unrelated to the source coefficient A_B in linear Section7f.
Use the SAME four AB+C row maps of §SP.6.2 and define
Gamma(alpha,(i,j))=W(row_i(alpha),j),i=1,...,4,j=0,...,17.
The index order is lexicographic in (i,j),giving degree exactly C=72.
Theorem G_A. For every set Q0 of at least Lambda_S left vertices,
|Gamma(Q0,[C])|>5*epsilon0*A0.
Proof. Since m is even Lambda_S is an integer,and
m>=41*2^160>134217728. The complete direct incidence support lemma S1'
gives a union U of row images with |U|>2^(9m/16-8). If its neighborhood
had size <=5*2^(-r_F)2^s,W_A would give a row image and a walk
position outside that neighborhood,a contradiction. Canonical prime selection,field construction,AB+C and
the increments are the supplied algorithms and cost2^O(s_tag).
Proof.md outer Sections2,6 prove s_tag<=2pS<=(2p/Q_S)L,so the
complete row enumeration is polynomial in actual n.
## 5. The evaluated source-linear sampler
Let E be the actual three-condenser/kernel base of linear-and-bilinear-
primitives.md Section7f,output d,seed padded to s bits. On legal calls it
has strong error<=epsilon0=2^-r_F for every uniform affine source of
dimension>=K_B. Define Samp(x,alpha,j)=E(x,Gamma(alpha,j)),using Section4's
actual graph. This map is source-linear at every alpha,j. AF1--AF4 supply
its complete seed,padding,entropy and runtime bounds. The outer proof uses
this strong base error and G_A directly for fixed rank-good neighbors.
## 6. Conditional affine breaking
### 6.1 Conditional/average-entropy lift with arbitrary copies
Suppose(4)holds at degree r,entropy K,error eta. Put
L=(A,A_1,...,A_r),R=(B,B_1,...,B_r,Y,Y_1,...,Y_r).
Assume L--Z--R,(Y,Z)=(U_d,Z),fixed legal advice,and
Havg_inf(A|Z)>=K+log(1/delta). By(6),outside Z-mass at most delta,
Hinf(A|Z=z)>=K. Conditional independence of the full tuples and uniform Y
then meet every hypothesis of(4). Averaging by(5),joint TV retaining Z is
at most eta*(1-P[bad])+P[bad]<=eta+delta. The arbitrary source-copy
quantifiers persist. Taking eta=delta=e_F/2 needs entropy K+log(2/e_F).
This is conditional averaging,not a new evaluated algorithm.
## 8. Selected outer use
Proof.md composes this graph/sampler and conditional breaker with tags,
fixed rank-good neighbors,whole-right profile postselection and its four-moment argument.
Attribution: extractor-to-sampler conversion is due to Zuckerman (1997);
the low-rate-condense/high-rate-walk architecture to Zuckerman,Theory of
Computing 3 (2007),§§2--3;the one AB+C support step to
Barak--Impagliazzo--Wigderson (credited in Part SP);the row-XOR/majority architecture to
Chattopadhyay--Liao,arXiv:2110.12652v1,§5. The affine interface
is compared to Li,ECCC TR23-023 revision 1.
## 0. Atomic proof map and selected construction
Legacy names denote included proofs: linear-and-bilinear-primitives.md
Sections1--5,7a--h,8 cover seeded-linear,weak-seed-hashing and bilinear tools;
restricted-nonmalleability-and-breaking.md Sections1--8 cover prime fields,
encoding,local conditional tools and antichain8.2;wrappers Sections2--4 cover
affine transformation,advice,CZ/E6 and NB,with restricted8.4's F# interface.
Part SP.1--5 proves rings,SP.7 proves M/I,and SP.6.2 proves S1'.
No missing GUV/CGL theorem is a premise. The dependency order is
linear primitives -> local core -> F#;rings/CZ/E6 -> advice condenser;
condenser/primitives/F# -> NB -> affine transformation;
SP.7 -> S1' -> graph/walk -> proof.md's tagged rank-good four-moment assembly.
## Paired read-once zero-period lemma
A proper read-once full-B2 formula on 2q distinct leaves,q>=1,with any
pairing of its leaves,has a nonzero zero-period on some product of q
allowed lines in the paired three-corner domain D. A zero-period v means
F(u+v)=F(u) everywhere on that flat. Flip input bits so every pair forbids
00;this preserves properness and the assertion. A safe axis line fixes
one endpoint to 1;the other allowed line is the diagonal 01,10.
First suppose all gates are AND or XNOR (1+x+y). Both have neutral input
1. Induct on q. If a root child contains m>=1 complete pairs internally,
fix its leaves with partners outside that child to 1. Those fixings use
at most one endpoint of each crossing pair. Neutrality reduces the child
to a formula on its 2m internally paired leaves,where m<q. Induction gives
a zero-period supported only on those leaves. The other root child is
unchanged. Extend the chosen internal lines by the safe crossing lines
and arbitrary allowed lines on the other pairs. This proves this case.
Otherwise every pair crosses the root,and each child has q leaves.
If one child contains an XNOR gate,choose a leaf from each of that gate's
two child subtrees,say u,v. Fix every other leaf of the chosen root child
to 1. Leave u,v free by fixing their partners in the other root child to
1. These are consistent one-endpoint-per-pair fixings. Neutrality reduces
the chosen child to XNOR(u,v). Flipping u and v preserves it and leaves
the other child unchanged,giving the desired period. If neither child
contains XNOR,both are pure ANDs. For an AND root choose one pair's
diagonal;the whole output is identically zero. For an XNOR root choose
all pair diagonals: flipping all q free coordinates swaps the two AND
monomials,preserving their XNOR. This includes q=1 and completes induction.
Now take arbitrary proper B2 gates and suppose no such zero-period exists.
For any product gate,consider a child whose leaf set contains at most
one endpoint of every pair. Fix all those leaves to safe 1. Its output
cannot be a controlling value at that gate's pin: otherwise the entire
other child is killed. Choose a leaf z in that other child,fix its mate
to 1,and complete to a product flat leaving z free. This is consistent
with the earlier fixings. Flipping z would preserve the whole formula,
a contradiction. Thus this child's all-ones output is noncontrolling.
Complement internal signals as needed so every signal takes value 1 at
the global all-ones input,absorbing complements in parent truth tables.
Leaves stay unchanged;complementing the root preserves zero-periods.
Every affine proper gate becomes XNOR. Every product gate whose two child
leaf sets each contain at most one endpoint per pair becomes AND: both
pins have noncontrolling value 1,and its output at 11 is 1.
Choose a minimal subtree V containing a complete pair. Every proper
subtree of V has at most one endpoint per pair. Consequently all gates
within V are AND or XNOR. Fix V's leaves with mates outside V to 1;
neutrality leaves a nonempty formula on its internally paired leaves.
Apply the auxiliary induction. Its period changes only those internal
leaves,leaving the rest of the full formula unchanged. Extend the flat
by the safe crossing lines and arbitrary remaining pair lines. This
contradicts the supposition and proves the lemma. The requirement that
ALL leaves be paired is essential to this proof.
## Tight incidence tree and terminal exclusion
Use the representation and incidence graph of proof.md §7. Suppose
at its terminal state a=k,q=k,M=2k,s=M-1. All M parameters belong to the
output component. The component bound s_0>=M-1 forces s_0=s. There are M+s
variable vertices and s object vertices,with at most 3s incidences.
Connectedness and M=s+1 force exactly 3s incidences and a bipartite tree;
every object has three distinct variable neighbors.
The output is a gate output. Indeed an affine output on M>=2 parameters
has a nonzero kernel direction v. For any direction v the domain D has
a safe a-flat parallel to v: in each pair choose the allowed line parallel
to its nonzero projection (three corners supply all three directions),
or any allowed line if the projection is zero;keep free coordinates.
An affine output would thus have zero derivative on such a flat,
contradicting the separation hypothesis of §7.
Every auxiliary variable has degree at least two. If it had degree one,
its sole incidence would be a row,since it is not the output. Delete
that row and its uniquely supported auxiliary column. Invertibility of
the remaining block follows by expansion on that column. The represented
function,parameters,Phi and D are unchanged,while s drops to M-2,
contradicting the same component bound. Every nonoutput gate variable has
degree at least two (its defining gate and a consumer);dead gates have
already been removed. The output gate variable has degree one. Parameters
have degree at least one. Their total degree is thus at least
M+2(s-1)+1=3s,so equality forces every parameter to have degree one and
every nonoutput gate/auxiliary variable to have degree two.
Root this incidence tree at the output variable. Each object has one
rootward and two outward variables. For a row the rootward variable is
auxiliary: parameters are leaves and gate variables never occur in rows.
Solve the row for that auxiliary,obtaining an XOR or XNOR of its two
outward variables. For a gate the rootward variable is its designated
output: the native directed path through consumers to the final output
exists by dead-gate removal and must be the unique incidence path.
Thus no product gate is inverted. Orienting the tree computes a proper
full-B2 read-once formula from all M parameters,each used exactly once.
It agrees with the represented function (even off D). The paired lemma
supplies a nonzero zero-period on a safe k-flat in D. Injective affine Phi
maps it to a nonzero original direction and a k-flat,contradicting
separation. The alleged tight terminal state is impossible.
# SP. Finite rings and direct incidence support
The following proofs supply E_p(((2-nu)/(14+(h_R+5/2)/D_E))) from restricted-nonmalleability-and-breaking.md
for every prime p >= 2^(2^90),without CZ Lemma 7.3,Bourgain's filling
lemma or an unspecified constant.
Attribution: the ring/stabilizer framework builds on the
self-improving approximate-ring construction of Terence Tao,*The
sum-product phenomenon in arbitrary rings*,Contributions to Discrete
Mathematics 4(2),2009,pp. 59--82,especially Lemma 1.4,Proposition 3.1,
and Theorem 1.2. Official source:
https://cdm.ucalgary.ca/article/view/61994
The sumset tools below are Ruzsa's inequalities and Petridis's minimal-
growth argument. We work only in commutative rings
with identity. No inclusion AB+AC subset A(B+C) is used: that inclusion is
false in general. Scalar distributions and every product of sets below
are distinguished explicitly.
Notation: nA denotes an n-fold sumset,with 0A={0};AA={ab:a,b in A};
aA is multiplication by the single element a. A unit is an invertible
ring element. All sets used in cardinality estimates are finite and
nonempty unless indicated otherwise. Products in a product of fields are coordinatewise. The selected ring below
is Kp x Kp,using QUADRATIC-FIELD multiplication inside each coordinate;
its four F_p basis coordinates do not multiply coordinatewise.
## SP.1. Fully quantified additive tools
### SP.1.1 Ruzsa triangle and covering
For finite nonempty subsets X,Y,Z of an abelian group,
|X-Z| |Y| <= |X-Y| |Y-Z|. (SP-1)
For each w in X-Z choose one representation w=x_w-z_w. The map
(w,y) -> (x_w-y,y-z_w) is injective: summing the two coordinates recovers
w and the chosen x_w then recovers y. This proves (SP-1).
If |X+Y| <= M|Y|,there is T subset X with |T|<=M and
X subset Y-Y+T. (SP-2)
Choose a maximal collection of disjoint translates t+Y,t in X. Their
number is at most M. Every other x+Y intersects one of them,which gives
x in t+Y-Y. The same conclusion follows from |X-Y|<=M|Y| by using -Y.
### SP.1.2 Explicit iterated sumset estimate
If |U+V|<=K|U| for finite nonempty U,V and K>=1,then for all nonnegative
integers m,n with m+n>0,
|mV-nV| <= K^(m+n) |U|. (SP-3)
Here is a proof rather than an invocation of Plunnecke's theorem. Choose
nonempty X subset U minimizing k=|X+V|/|X|. Thus k<=K and |Y+V|>=k|Y|
for every Y subset X,including Y empty. For any finite set C,order its
elements and add the translates X+c successively. For the next c,let
Y={x in X:x+c belongs to the previous union X+C_old}. Then Y+V+c is
contained in the previous union X+V+C_old. Consequently the newly added
part of X+V+c has size at most |X+V|-|Y+V|<=k(|X|-|Y|). The latter
parenthesis is exactly the newly added size of X+c. Summing proves
|X+V+C| <= k |X+C|.
Induction gives |X+jV|<=k^j|X| for every j>=0. Apply (SP-1) to mV,-X,nV
to obtain |mV-nV|<=k^(m+n)|X|<=K^(m+n)|U|. This includes m=0 or n=0.
In particular,|A+A|<=K|A| implies
|mA-nA| <= K^(m+n)|A|. (SP-4)
If P contains a set C of size N and |P-P|<=L N,then N<=|P|<=L N,
and (SP-3) with U=P,V=-P yields
|mP-nP| <= L^(m+n+1)N. (SP-5)
The next sharper bound proves(SP-5),since L>=1.
### SP.1.3 A sharper form of (SP-5)
If C subset P,|C|=N and |P-P|<=LN,then,for k=m+n>=1,
|mP-nP| <= L^(m+n) N. (SP-5s)
Indeed put Q=|P|>=N and k0=|P-P|/Q<=LN/Q. Applying (SP-3) with
U=P,V=-P,and interchanging the signs if necessary,gives
|mP-nP|<=k0^k Q<=(LN)^k/Q^(k-1)<=L^k N. All sets are nonempty;
k0>=1 since a translate of P lies in P-P. The weaker (SP-5) follows.
## SP.2. Dyadic overlaps and weighted unit actions
Let C0 be M units in a commutative ring,K>=2,with
|C0-C0|,|C0/C0|<=KM. Put H=1+ceil(log2(2K)). Assume either
Z=|(C0-C0) outside units|<M/(2^20*K^14*H^2),
or the ring is a product of fields,each coordinate fiber of C0 has size
at most F_j,and W=sum_j F_j<M/(2^20*K^14*H^2). (SP-Mi0)
Then some b0 in C0,B0 subset C0 and finite unital ring T satisfy
|B0|>=M/(4KH),B0 subset b0T,|T|<=2^17*K^12*M. (SP-Mi1)
Product energy equals quotient energy: a/b=c/d iff ad=bc. Cauchy on C0/C0 gives sum_(a,b)|aC0 intersect bC0|>=M^3/K.
Normalize intersections by M. Values below1/(2K) contribute at most
M^2/(2K). Partition the rest into H ranges[mu,2mu),mu=2^-j,
j=0,...,ceil(log2(2K));the top includes1. One range contributes
at least M^2/(2KH),hence has at least M^2/(4KHmu) pairs. Choose b0 with
at least M/(4KHmu) partners B0,adding b0 if absent. Normalize
C=b0^-1 C0,G=b0^-1 B0. Put kappa=4KHmu,L=K^2/mu. Thus
1 in G subset C,|G|>=M/kappa,1/(4K)<=mu<=1,
|gC intersect C|>=mu M for every g in G.
Triangle with middle W=gC intersect C gives |gC-C|<=|gC-W||W-C|/|W|<=K^2*M/mu=LM. (SP-Mi2)
Let S={x:|xC-C|<=2048K^2L^4 M}. Cauchy counts at least
|G|^4/(2048K^2L^4 M) ordered quadruples satisfying xg1-g2=xg1'-g2'.
Put d=g1-g1',so xd=g2-g2'. For each nonunit d at most |G|^2
quadruples occur,giving bad count<=Z|G|^2 in the first case.
In a product of fields,a nonunit d has some zero coordinate j,hence
g1_j=g1prime_j. There are at most |G|F_j such ordered pairs.
For each pair choose g2prime in G;the collision equation forces the
ENTIRE element g2=x*(g1-g1prime)+g2prime uniquely. Its membership in G
can only reduce the count. Thus this coordinate contributes at most
|G|^2 F_j collisions,and the bad count is at most |G|^2 W,where
W=sum_j F_j. This also holds for nonunit x. Unit normalization and
taking subsets preserve the fiber caps.
In either case write V=Z or W. Since
4096kappa^2 K^2L^4=2^16*K^12*H^2/mu^2<=2^20*K^14*H^2 and |G|>=M/kappa,
V<|G|^2/(4096K^2L^4 M).
At least half the quadruples have unit d. Choose one and write
d=g1-g1',xd=g2-g2'. For D=union_{g in {g1,g1',g2,g2'}}gC,
|D-C|<=4LM;(SP-3) gives |2D-2D|<=(4L)^4 M.
Both xdC-dC and xdC+dC lie in this set. Unit dilation yields
|xC-C|,|xC+C|<=256L^4 M. (SP-Mi3)
Repeat Cauchy at this sharper bound. The bad count bound remains
|G|^2 V<|G|^4/(4096K^2L^4M)<=|G|^4/(512L^4M),since K>=2.
Hence each x in S supplies at least |G|^4/(512L^4M) unit-d quadruples.
Every such quadruple determines at most one x. There are at most |G|^4
quadruples in total,proving S finite and
|S|<=512L^4M<=2^17*K^12*M.
For every finite nonempty A the extension estimate is
|A+C-C|<=K^2|A+C|. (SP-Mi4)
Choose nonempty P subset C minimizing k=|P-C|/|P|<=K. The Petridis
extension in SP.1 gives |A+P-C|<=k|A+P|. Triangle with
X=A-C,Y=-P,Z=-C gives |A+C-C|<=k^2|A+P|<=K^2|A+C|.
Thus (x+1)C-C subset xC+C-C and the plus bound in(SP-Mi3),K^2<=L imply translation
by1 preserves S. For W=gC intersect C,triangle gives
|gxC-C|<=|gxC-W||W-C|/|W|<=256(K/mu)L^4M<=1024K^2L^4M,
since |gxC-W|<=|gxC-gC| and |W-C|<=KM. Thus every g in G preserves S. These bijections permute finite S;their
inverses preserve S. Conjugation gives translation by every G-monomial
and its negative. Since0 belongs to S,the additive group T generated by
these monomials lies in S. It contains1,G and is closed under products.
B0=b0G subset b0T proves (SP-Mi1).
## SP.3. Coordinate-fiber bound on bad differences
Let C be a nonempty subset of a product F^d of a finite field F,|C+C|<=K|C|,and let pi_j denote the
j-th coordinate projection. Put
Z_j=(C-C) intersect ker(pi_j).
Choose one representative c from C for each value of pi_j(C). Their
translates c+Z_j are disjoint,so
|pi_j(C)| |Z_j| <= |C+Z_j| <= |2C-C| <= K^3|C|. (SP-17)
A vector in a product of fields is a nonunit exactly when at least one coordinate is
zero. Taking a union bound in (SP-17) proves
|(C-C) outside the units|
<= K^3 |C| sum_(j=1)^d 1/|pi_j(C)|. (SP-18)
This bounds the number of distinct bad differences,not just the number
of pairs yielding them. That distinction is needed in §SP.2.
## SP.4. Quadratic-field pair hypotheses
Put nu=2^-80,D_E=3080960812491/181138871323+5nu,
h_R=543816613969/181138871323+2nu,gamma=(2-nu)/(14+(h_R+5/2)/D_E),
eta0=gamma/D_E,delta_R=gamma+nu,q_R=14gamma-eta0+nu.
These are fixed rationals;the selected eta=gamma/D_E. Exact algebra gives
2-(h_R+7/2)*eta0=q_R,q_R>1,delta_R<1/7,h_R<7,gamma<1/7,q_R-delta_R>1.
Let p>3 be prime and Kp its canonical quadratic field from restricted1.
For nonempty A subset R=Kp x Kp require:
(a) p^(2-h_R*eta0)<=|A|<=p^2;
(b) C0=A intersect units has size>=|A|/2;
(c0) both field projections of C0 have size>=p^q_R and every field-coordinate fiber of C0 has size<=4p^((5/2)eta0);
(c1) every D subset A with |D|>=p^-delta_R*|A| has both projections>p;
(d) no such D is contained in either graph(z,c*z) or(z,c*conj(z)).
## SP.5. Quadratic-field expansion by the minus stabilizer
Under SP.4,for p>=2^(2^90),
max(|A-A|,|C0/C0|)>p^gamma|A|. (SP-19m)
Assume failure. Set N=|A|,M=|C0|>=N/2,K=2p^gamma and
H=1+ceil(log2(2K)). Write B=log2p>=2^90. Since gamma<1/7,
H<=4+gamma*B<=B. By(c0),W=sum_j F_j<=8p^((5/2)eta0).
Using M>=N/2>=p^(2-h_R*eta0)/2 gives
2^20*K^14*H^2*W/M<=2^41*H^2*p^-nu<1. (SP-21m)
Indeed14gamma+(h_R+5/2)eta0=2-nu. On2^j<=B<2^(j+1),j>=90,
41+2log2B<43+2j,while nu*B>=2^(j-80). Base223<1024 and
doubling dominate slope2. SP.2's fiber case thus gives B0,T with
|B0|>=M/(4KH)>=N/(16H*p^gamma),
|T|<=2^17*K^12*M<=2^29*p^(12gamma)*N<p^4.
Here N<=p^2,12gamma<12/7 and29/2^90<2/7 prove the last inequality.
B0 qualifies for(c1),(d),since4+log2H<nu*B anddelta_R=gamma+nu.
Unit scaling preserves projection cardinalities.
The unital subring T is an F_p vector subspace: addition has characteristic
p,so scalar multiplication is repeated addition. Its two projections
are F_p subspaces of Kp and contain more than p elements by(c1),hence both
are all Kp. We now classify this SUBDIRECT product without a black box.
Let J={b:(0,b) in T}. It is an ideal of Kp: for every c in Kp the second
surjection supplies(a,c) in T,so(0,b)(a,c)=(0,bc) is in T.
If J has a nonzero element,J=Kp. Subtracting(0,c) from each(a,c) shows
T=Kp x Kp,contrary to |T|<p^4. Thus J=0;the first projection is
bijective and T is the graph of a unital F_p-field automorphism phi of Kp.
The canonical quadratic-field proof gives phi=identity or conjugation.
Writing b0=(b1,b2),its unit scaling is consequently a graph
(z,c*z) or(z,c*conj(z)),c=b2/phi(b1). Yet B0 qualifies for(d) and lies
in that graph. This contradiction proves(SP-19m),with all subset quantifiers.
The selected family uses gamma=((2-nu)/(14+(h_R+5/2)/D_E)),nu=2^-80,
eta=((2-nu)/(14D_E+h_R+5/2)),theta=eta/4=((2-nu)/(4*(14D_E+h_R+5/2))).
Its cap eta<=1/120 holds by cross multiplication. CZ uses unit deficit
181538871323/181138871323+nu,relative energy exponent3080960812491/181138871323+4nu and hereditary cutoff
p^(2-h_R*eta). The finite guard eta*logp>=64/nu pays2^64.
Every actual NB call has eta*(r+a)>=4*p0*w_N>=64*2^80.
Since eta<=1/120,r+a>=7680*2^80>2^90,independently paying the ring guard.
The image cutoff,unit deletion,clipping and fiber partition costs are bounded explicitly in restricted Sections6--7.
The supplied conditional proof gives45p^-eta,the matched marginal90p^-eta,
and the four-block advice lemma14*2^-sig,with all actual algorithm guards.
# SP.6. Direct point-line incidence support growth
§SP.7's full incidence proof gives the AB+C four-row support bound.
The evaluator uses fixed integer arithmetic,canonical primes,AB+C and
bit truncation;no scalar hitting/peeling lemma is used.
## SP.6.0. One support-growth step
### SP.6.1 Exact finite endpoint
Require m>=134217728;this conservative guard implies m>=4. The direct
incidence count below proves the finite endpoint |U|>2^(9m/16-8),keeping
its explicit additive constant rather than replacing it by a fixed exponent loss.
### SP.6.2 The one-step map and its support lemma
For an integerm>=1,letp_m be the least prime strictly above2^m. The supplied
prime-fields Sections1--2 inrestricted-nonmalleability-and-breaking.md prove
2^m<p_m<2^(m+1) and deterministic trialdivision bitcost
O(2^(3m/2+2)(m+1)^2). Read an m-bitstring as its binaryinteger in[0,2^m),
which injects intoF_(p_m). Define
E_m(x1,x2,x3)=(x1*x2+x3 modp_m) mod2^m. (SP-S2)
For a3m-bitstringalpha=(alpha1,alpha2,alpha3),thefourROWS arealpha1,
alpha2,alpha3,E_m(alpha). This is the same actual deterministic row map.
**Lemma S1' (direct incidence support).** For everyintegerm>=134217728 and
Q subset{0,1}^(3m) with|Q|>=2^(3m/2-17),the unionU of itsfourrowimages
satisfies|U|>2^(9m/16-8).
Proof. SupposeN=|U|<=2^(9m/16-8). N>=1 sinceQ is nonempty. RegardU
as an injective subset ofF_p,p=p_m,and put
W={y inF_p:the lowm bits of its canonical representative belongtoU}.
Sincep<2^(m+1),each member ofU has atmosttwo preimages,so|W|<=2N.
TakeP=U xW andL={y=a*x+c:a,c inU}. There areatmost2N^2 points and
EXACTLYN^2 distinctnonvertical lines: slope andintercept determinea,c,
includinga=0. Every triple inQ determines its line(a=alpha1,c=alpha3)
and its point(x=alpha2,y=alpha1*alpha2+alpha3 modp). The fourth-row
condition puts y inW. This map into point-line incidences is INJECTIVE:
the line determinesalpha1,alpha3 and the point's first coordinate
determinesalpha2. Thus|Q|<=I(P,L). Zero slopes require no deletion.
Choose d=ceil(N^(4/3)) as an ANALYSIS parameter,not an evaluator operation.
SinceN<=2^(9m/16),N^(4/3)<=2^(3m/4),andm>=4 gives
1<=d<=2^(3m/4)+1<2^m<p.
Indeed2^m-2^(3m/4)=2^(3m/4)(2^(m/4)-1)>=8>1 at m>=4.
The complete incidence theorem(SP-I) below thus applies. Set
v=floor(2N^2/(d+2))<=2N^(2/3).
As d-1<N^(4/3) andd<=N^(4/3)+1,it gives
|Q|<=I(P,L)<=2N^2(v+1)+d(d-1)+N^2
<5N^(8/3)+3N^2+N^(4/3)<=9N^(8/3)<16N^(8/3).
The middle inequality usesN>=1. Taking logarithms with the lower bound
onQ gives log2N>9m/16-63/8>9m/16-8,a contradiction.
All counts concern arbitraryQ and its actual row images;no independence
is asserted. Complete polynomial/multiplicity theorem (SP-I) is proved below.
At the selected call s=2ceil((50/7-2^-100)*S/2),r_F=S,
m=s+34S>=41S>=41*2^160>134217728.
## SP.7. Complete polynomial incidence proofs
The complete multiplicity and polynomial-incidence argument below is adapted
from published1d75672c-5080-4151-89d2-cd776a09d06b,actualimmutable009e40109ba442d132016de118d39c2fdf92f206,
and Mark Lewko,arXiv2609.27023v3,Lemma2/Section4(CCBY4.0). The full proof
is supplied;neither source acceptance,the paper nor its Lean formalization
is an external premise. The direct support count in SP.6.2 is this candidate's
application. Write t_+=max(t,0),andp=infinity incharacteristiczero.
Algebra used below. Polynomial rings over fields are unique factorization
domains. One elementary justification is induction on the variables:
over a UFD,the product of primitive polynomials is primitive (reduce
modulo each prime divisor of the coefficients;the resulting polynomial
ring over a domain has no zero divisors). Clearing denominators then
transfers factorization and its uniqueness from the fraction field's
Euclidean univariate polynomial ring. This is Gauss's lemma,and also
shows that primitive irreducibles remain irreducible on passing to the
fraction field. Degrees,including degree in a chosen group of variables,
add under multiplication of nonzero polynomials,by their leading forms.
A nonzero univariate polynomial of degree j has at most j distinct roots:
divide successively by the distinct linear factors.
### SP.7.1 Lemma M: multiplicities on parallel lines
Let K be a field,0<=d<p an integer,and F in K[x,y] nonzero and squarefree
(no repeated irreducible factors),with total degree <=d. Let lambda_i be
finitely many distinct parallel lines,with q_i on lambda_i,such that
F restricts to a nonzero polynomial on every lambda_i. For a nonzero common
direction v let mu_i be the order of u in F(q_i+uv). Then
sum_i (mu_i-1)_+ <= d(d-1). (SP-M)
Proof. Empty collections and constant F are immediate. Interchange x,y
if the lines are vertical,then shear X=x,Y=y-a*x,where a is their common
slope. The lines become Y=b_i with distinct b_i,and q_i=(xi_i,b_i).
The invertible substitution preserves degree and squarefreeness. Factor
the transformed F as A(Y)G(X,Y),removing all factors independent of X.
Since no restriction is zero,A(b_i)!=0;multiplicities are unchanged.
If G is constant,all mu_i=0. Otherwise each irreducible factor f of G has
1<=deg_X f<=d<p,so f_X is nonzero of smaller X-degree. Write G=fQ.
Squarefreeness gives f not dividing Q;from G_X=f_X Q+f Q_X,primality
shows f cannot divide G_X. Thus G_X!=0 and gcd(G,G_X)=1.
Put h_i=(mu_i-1)_+,d0=deg G and d1=deg G_X<=d0-1. Both
G(xi_i+u,b_i) and G_X(xi_i+u,b_i) are divisible by u^h_i,since
(u^mu g)'=mu*u^(mu-1)g+u^mu g'. Let Pi_T be the space of polynomials of
total degree <=T,of dimension binom(T+2,2). For an integer T sufficiently
large,in particular T>=d0+d1,let W_T be the common kernel of the
sum h_i coefficient functionals [u^j]H(xi_i+u,b_i),0<=j<h_i.
They are independent: for each (i,j) the polynomial
(X-xi_i)^j product_(v!=i) (Y-b_v)/(b_i-b_v)
has just that prescribed coefficient equal to one and all others zero.
Take T also at least all their degrees. Thus codim W_T=sum h_i.
Every multiple of G or G_X lies in this kernel,so
J_T=G Pi_(T-d0)+G_X Pi_(T-d1) is contained in W_T.
Coprimality and degree additivity give the intersection of these summands
as G G_X Pi_(T-d0-d1). Multiplication is injective. Hence
codim J_T=binom(T+2,2)-binom(T-d0+2,2)-binom(T-d1+2,2)
+binom(T-d0-d1+2,2)=d0*d1.
Thus sum h_i=codim W_T<=codim J_T=d0*d1<=d(d-1). This includes d1=0.
The restriction d<p is essential: x^p-y has zero x-derivative and order p
on every horizontal line through (xi,xi^p).
### SP.7.2 Lemma I: the d-form incidence bound
Let P be m points of k^2 and L be n distinct nonvertical lines over a field
k of characteristic p. For every integer 1<=d<p put v=floor(2n/(d+2)). Then
I(P,L)<=m(v+1)+d(d-1)+n<=m+n+2mn/(d+2)+d(d-1). (SP-I)
Proof. Empty P or L is immediate. Seek S=sum_(j=0)^v z^j S_j(x,y),with
deg S_j<=d,satisfying S(x,ax+b,a)=0 as a POLYNOMIAL IDENTITY in k[x]
for each y=ax+b in L. There are (v+1)binom(d+2,2)>n(d+1) coefficients and
at most n(d+1) homogeneous linear conditions,so a nonzero S exists.
Replace S by the product of its distinct irreducible factors in k[x,y,z].
Neither degree increases. Each old zero restriction forces a factor's
restriction to vanish,since k[x] is a domain,and that factor survives.
Thus S is squarefree,deg_(x,y) S<=d and deg_z S<=v,with all identities.
For q=(q1,q2) in P let r(q) count its incident lines and write
R_q(u,z)=S(q1+u,q2+zu,z)=sum_j c_(q,j)(z)u^j.
A monomial x^alpha y^beta z^gamma contributes to u^j with z-degree at
most gamma+j,so deg c_(q,j)<=v+j. Every incident slope a is a root of
every c_(q,j),by the imposed identity. Distinct lines through q have
distinct slopes. Root counting gives c_(q,j)=0 when j<r(q)-v. Thus
u^((r(q)-v)_+) divides R_q in k[z,u].
Work over K=k(z),keeping z indeterminate,and set
lambda_q={(q1+u,q2+zu):u in K}. These parallel lines are distinct:
equality of their intercepts q2-z*q1 forces equality of both coordinates.
S stays squarefree in K[x,y]. Indeed each irreducible factor involving
x or y is primitive over k[z],hence irreducible over k(z) by Gauss's
lemma. Distinct such factors f,g cannot become associates: clearing
denominators would give af=bg for nonzero a,b in k[z];primality of f and
f not dividing b would force f dividing g. Factors in z alone become units.
Let E={q:R_q=0},h=|E|. On the remaining lines apply (SP-M),using
mu_q>=(r(q)-v)_+. Their incidences total at most
(m-h)(v+1)+sum_(q notin E)(r(q)-v-1)_+
<=(m-h)(v+1)+d(d-1).
For q in E,substitute u=x-q1 in R_q=0. Division in y by the monic
linear polynomial l_q=y-q2-z(x-q1) shows l_q divides S in k[x,y,z].
These factors are irreducible (their quotient ring is k[x,z]),pairwise
non-associate,and of z-degree one. Their product divides S,so h<=v.
A given line containing s points of E contributes s<=1+binom(s,2);
each pair of distinct points lies on at most one line. Consequently
I(E,L)<=n+binom(h,2)<=n+h(v+1).
Adding gives (SP-I). All zero restrictions above are polynomial identities,
not merely pointwise zeros on a finite field.
## SP.8. Sources and selected scope
The multiplicity and d-form incidence proofs are adapted with attribution
from Mark Lewko,*A Szemeredi-Trotter Theorem in Arbitrary Fields*,
arXiv2609.27023v3(7Oct2026),https://arxiv.org/html/2609.27023v3,
CCBY4.0,and the complete transcription in publishedYukon1d75672c,
commit009e40109ba442d132016de118d39c2fdf92f206. All needed algebra,counting,polynomial identities and characteristic
restrictions are proved here.
The direct Q-to-incidence injection,zero-slope handling and finite xi
calibration above are this candidate's changes.
The AB+C support/sampler architecture has its historical hitting/nearhalf
origins in b15a467b and f7447a6f/04e91ff4/7d6289bb through the credited
published precursors.
The failed/unpublishedc1b72 RG4/9/4 route and ownsubmittedfinite-rg-antichain
were material research precursors,credited in the public note;the current
support theorem is the independent direct count just proved.
