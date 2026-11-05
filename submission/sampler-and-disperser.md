# Finite support graph, linear sampler, and conditional breaking

All logarithms are binary, addition is XOR unless a field is specified, and
statistical distance is TV = half the L1 distance. A k-source is a finite
random variable with every atom at most 2^(-k); k may be real. Graphs have
labelled edges and may have loops or repeated edges. Every constant is fixed
independently of the eventual outer input length N.

The supplied ordinary condenser, seeded-linear extractor and affine breaker
yield the finite graph/sampler and conditional-breaking reductions below.
proof.md supplies the tagged affine-permutation outer assembly and finite guards.

## 1. Supplied contracts and elementary probability facts

### 1.1 Finite seeded-linear contract

The supplied seeded-linear module constructs a deterministic strong extractor
E: {0,1}^n x {0,1}^D -> {0,1}^(2^h ell), linear in its first argument
at every fixed seed. Its parameters are positive integers n,h,ell,e, with

  a_E=e+h+4, U=4n, rho_E=a_E+2 ceil(log U)+2,
  R=4h, s0=least power of two >=2rho_E,
  S=least power of two >=(R+1)rho_E,
  K_h=ell+2a_E,
  K_(i-1)=ceil[(2R/(R-1))(K_i+a_E+S/2+1)]  (i=h,...,1).

The required guard is K_0<=n. Put
m_h=ceil[S ceil(K_(h-1)/(S-rho_E))/2] and let W be the least power of
two >=m_h. The exact complete seed length is D=s0+hS+W. The contract is

  TV((E(X,U_D),U_D),(U_(2^h ell),U_D)) <= 2^(-e),
  D <= 6ell+(h+7)S+18a_E+12,
  K_0 <= 2*2^h(ell+4a_E+S+4),                               (1)

for every k-source X with k>=K_0 and an independent uniform seed. Its
canonical field/setup/evaluation algorithms have a supplied O(n^8) bit-time
bound on admissible parameters. A target shorter than 2^h ell is obtained
by taking a prefix. A longer prescribed seed is implemented by ignoring its
suffix. Neither operation harms linearity or strongness. Its construction and proof are supplied by the seeded-linear module.

### 1.2 Ordinary-condenser contract used by the graph

The supplied ordinary-condenser theorem gives, for fixed rational
0<beta<tau<1 and 0<nu<tau, fixed explicit integers ell_O>=2,T>=1,m_*>=1
and rational gamma>0. For input length r satisfying its guards, its row
map has exactly D0=(ell_O+1)^T rows of v(r)=ceil(r/ell_O^T) bits.
Every source of entropy beta*r is mapped within

  T*2^(-gamma*v/4)+D0*2^(-nu*v)                              (2)

of a convex mixture of row matrices, each component having one fixed row
of min-entropy at least (tau-nu)v. The selected good row can vary between
components; neither the selection nor the mixture is algorithm input. Its
uniform row evaluator and parameter algorithm cost 2^O(r) bits.

Here beta=1/2. First select J and ell_O below, then put
R_G=176128000ell_O, nu=1/R_G, tau=1/2+2/R_G. The ordinary module's explicit
choices may be recorded as

  J=least integer with (1+1/21500)^J*(beta/4)>=9/10,
  ell_O=8*3^J,
  gamma=min(beta/344000,beta/100,(1-tau)/100,1/100),
  g_O=gamma/(64ell_O), T=ceil((tau-beta)/g_O),
  M=ell_O^T, D0=(ell_O+1)^T, kappa=1/M.

Here gamma=1/688000 and g_O=1/(44032000ell_O). Hence
T=ceil((2/R_G)*44032000ell_O)=ceil(1/2)=1.
Thus M=ell_O and D0=ell_O+1.

The primitive supplies the exact m_*. Integer log ceilings below use bit length:

  delta_O=min(gamma/8,nu/2),
  Q*=max(m_*,ceil((8/gamma)*ceil(log(2T))),
                ceil((2/nu)*ceil(log(2D0))),ceil(1/delta_O),1),
  t0=M*Q*.                                                  (3)

For every integer t>=t0 the primitive's guards hold: t>=M*m_* and
v=ceil(t/M)>=Q* meet its two logarithmic error guards. In particular,
T<=2^(gamma*v/8)/2 and D0<=2^(nu*v/2)/2, so (2) is at most
2^(-delta_O*v)<=1/2. Also kappa*t<=v<=t. Thus the resulting exact
contract, denoted O below, is a D0-row map T_O on t bits, at input entropy
beta*t, within 1/2 of a mixture with one row of entropy at least real (1/2+1/R_G)v.
There is no ceiling on that real entropy threshold. This invokes the supplied ordinary theorem.

### 1.3 All-seeds-retained affine-breaker contract

Let F:{0,1}^n x {0,1}^d x {0,1}^a -> {0,1}^m be one fixed
deterministic map, with a>=1. Its supplied degree-r contract at source
entropy K and error eta means the following universal statement. For
arbitrary correlated tuples

  Left=(A,A_1,...,A_r),
  Right=(B,B_1,...,B_r,Y,Y_1,...,Y_r),

assume Left and Right are independent, H_inf(A)>=K, Y is uniform on d
bits, and fixed advice strings alpha_i all differ from alpha. Set

  V=F(A+B,Y,alpha), V_i=F(A_i+B_i,Y_i,alpha_i).

Then, with U_m independent of all displayed side information,

  TV((V,V_1,...,V_r,Y,Y_1,...,Y_r),
     (U_m,V_1,...,V_r,Y,Y_1,...,Y_r)) <= eta.                 (4)

Only A needs entropy and only Y uniformity. No A_i need be a function
of A. All within-side correlations, repetitions, and seed equalities are allowed. Tampered advice
strings need not differ from one another. For a smaller degree, pad with
repeated permissible copies and discard their outputs/seeds. This chosen
branch directly uses the all-seeds-retained contract of the actual affine
construction; it does not need a dummy-tampering adapter for a weaker
contract that retains only tampered seeds.

### 1.4 Probability identities

For any common marginal Z, grouping the TV sum by z gives

  TV((V,Z),(V',Z))=sum_z P(Z=z) TV(V|z,V'|z).                 (5)

A common deterministic map or random conditional kernel cannot increase TV:
sum differences over each fiber and use the triangle inequality. In all
ideal laws below a displayed U is fresh uniform independent of the entire
retained side information.

Define average conditional min-entropy by

  Havg_inf(A|Z)=-log sum_z P(Z=z) max_x P(A=x|Z=z).

If Z has at most 2^ell possible values, then

  sum_z max_x P(A=x,Z=z) <= 2^ell max_x P(A=x),
  Havg_inf(A|Z) >= H_inf(A)-ell.                             (6)

If Havg_inf(A|Z)>=K+log(1/delta), with 0<delta<=1, the total probability
of z with H_inf(A|Z=z)<K is at most delta: on that set the maximum
conditional atom exceeds 2^-K, while its expectation is at most
 delta*2^-K. These facts require no unstated independence.

## 2. An explicit small-bias Cayley graph

Fix integers v,z with z>=200R_G+2 and100R_G*ceil(log v)<=z.
Set r_F=z+ceil(log v) EXACTLY and c=2r_F. Then

  2z<=c=2(z+ceil(log v))<=2(1+1/(100R_G))z.              (7)

Use the first irreducible monic binary polynomial of degree r_F in the
fixed coefficient ordering, with the quotient's standard binary basis.
The finite existence, irreducibility test and 2^r_F*poly(r_F) search bound
are proved in linear-and-bilinear-primitives.md, Section4 'Deterministic
field setup at the actual outer length'. This search is part of the
evaluator. Binary quotient multiplication uses O(r_F^2) bit operations.
In the outer schedule r_F=L+lg(v)=O(log N), so setup is polynomial in
the actual N. No table, field advice or randomness is supplied.

For each 1<=m<=v and each label (u,w) in F^2, define an m-bit increment

  g(u,w)_j=<w,u^j>, j=0,...,m-1,                             (8)

where brackets are the ordinary coordinate dot product and u^0=1 also
when u=0. On vertices F_2^m use the 2^c labelled edges x -> x+g(u,w).
A fixed label is an involution, so the transition matrix P is symmetric
and doubly stochastic, with the uniform stationary distribution.

For any nonzero a in F_2^m, the character (-1)^<a,x> has eigenvalue
E_(u,w)(-1)^<w,p_a(u)>, where p_a(T)=sum_j a_j*T^j. The polynomial is
nonzero of degree at most m-1. Averaging over w annihilates each nonzero
coordinate vector p_a(u), leaving its fraction of zeros. A nonzero degree-j
polynomial has at most j roots, by repeated division by T-alpha. Hence

  0<=eigenvalue<=(m-1)/2^r_F<=2^-z.                           (9)

The binary characters form an orthogonal basis, so ||P-J||<=lambda=2^-z,
where J projects onto constants in the uniform-normalized Euclidean inner
product. A labelled step is evaluated by successive field powers and dot
products in O(m*r_F^2) bit operations. The label count 2^c is not the
sampler degree, and the algorithm does not enumerate those labels.

## 3. Exact high-rate walk and strict support expansion

Set

  q=floor((1/2-1/R_G)*v/(z-2)), m=v-qc.                                   (10)

An input of exactly v bits parses into an m-bit start x0 followed by q
c-bit labels. Follow (8) to x1,...,xq and output W(x,j)=xj for
j=0,...,q. This fixes a deterministic degree-(q+1) map, with

  m>=v/R_G, q+1<=1+(1/2-1/R_G)*v/(z-2).                               (11)

For every Q subset {0,1}^v it has the strict implication

  |Q|>=2^((1/2+1/R_G)v-1)  =>  |W(Q,[q+1])|>2^(1-z)*2^m.            (12)

Proof. Let R be any right set of density rho<=2^(1-z) and Pi its coordinate
projection. Since ||Pi J Pi||=rho, (9) gives
||Pi P Pi||<=rho+lambda. A uniform v-bit encoding is exactly a uniform-start
walk with independent uniform labels. Thus

  Pr[every xj belongs to R]
    =<1_R,(Pi P Pi)^q 1_R> <=rho*(rho+lambda)^q.

This includes q=0, when both probability and bound equal rho. Since
rho+lambda<=3*2^-z<=2^(2-z), the number of all-bad encodings is at most

  2^[v+1-z+q(2-z)]=2^[v-1-(q+1)(z-2)].                      (13)

By (7) and (10), qc/v<1-1/R_G: it is at most
(1-2/R_G)(1+1/(100R_G))^2<1-1/R_G.
Indeed z/(z-2)<=1+1/(100R_G), and (1+a)^2<=1+3a for 0<=a<=1.
This proves the output-length claim, also immediate for q=0.
Moreover (q+1)(z-2)>(1/2-1/R_G)v, so (13) is
strictly smaller than 2^((1/2+1/R_G)v-1). If (12) failed, taking R=W(Q,[q+1])
would make every point of Q all-bad, a contradiction. This includes q=0.
A real size threshold means its ceiling for integer |Q|; strict slack
covers every floor. Only uniform encodings are independent-label walks.

## 4. Ordinary condensation to the exact finite support graph

For Q subset {0,1}^t of size at least Lambda=2^(t/2-20), take X uniform
on Q and independent uniform Z on40bits. Then Hinf(X,Z)>=t/2+20=(t+40)/2.
Apply O at length t+40, with v=ceil((t+40)/M). Let U be the union of all
D0 row supports over EVERY x in Q and z in{0,1}^40.
The actual row matrix lies in U^D0 surely. A repaired mixture Y within
TV 1/2 therefore satisfies P[Y in U^D0]>=1/2. In each component select
its promised good row. The event that all rows lie in U implies that row
lies in U, whose probability is at most |U|*2^(-(1/2+1/R_G)v). Averaging gives

  |U|>=2^((1/2+1/R_G)v-1).                   (14)

There is no log D0 entropy loss, no assumption that the same row is good
in all components, and no rounding-up of min-entropy. Apply (12) to U to obtain strict
support density greater than 2^(1-z) after walking all condenser rows.

For integers s>=1,z>=2 put A0=2^s, epsilon0=2^-z. For the supplied
ordinary condenser set

  b=R_G*M, t=bs, A=2^t, Dg=2^40*D0, kappa_G=R_G/2-1.                                      (15)

Here M>40: J>=2, ell_O>=72 and T>=1. Thus v=R_G*s+1<=t. Require

  t>=t0, v>=64, z>=200R_G+2, 100R_G*ceil(log v)<=z.           (16)

The walk supplies m>=v/R_G>=s, enough for projection.
For label (zpad,i,j), zpad in{0,1}^40, i in[D0], j=0,...,q, define
Gamma0(x,(zpad,i,j)) by walking row i of T_O(x||zpad). Define Gamma as its first s bits.
For every H subset {0,1}^m, each prefix has at most 2^(m-s) preimages,
so projection preserves any strict lower bound on support density. Equations
(12)--(16) therefore prove exactly

  |Gamma(Q,[C])|>2*epsilon0*A0 for every |Q|>=Lambda,        (G)
  C=Dg(q+1)<=Dg(1+floor((kappa_G+1)s/z))<=g_G(1+s/z), g_G=(kappa_G+1)Dg.         (17)

Here q<=(kappa_G*s+1/2-1/R_G)/(z-2)<=(kappa_G+1)s/z:
(s-1/2+1/R_G)z>=(s/2)z>=2(kappa_G+1)s.
Lambda is integral on the outer schedule (t is even and t>=40). Repeated edges pad the degree without changing support.

For z>=s+2 a constant map already proves (G), since its singleton
neighborhood exceeds 2^(1-z)*2^s<=1/2. Select this branch before constructing
a field depending on z. For z=0,1, (G) is impossible. Finitely many
excluded pairs in a bounded s interval with z>=2 can be handled by listing
all A0 neighbors and absorbing their degrees into a fixed constant. This
is not a general all-(s,z) construction: fixed z excludes infinitely many
s through (16), and finite fallback does not repair that quantifier. The
intended schedule below satisfies (16) eventually at every integer N.

### 4.1 Intended graph schedule and its actual cost

For fixed h>=11 put u=h+1,B0=1200u^2. This is a coarser schedule
than proof.md's final specialization. Set L=ceil(log(N+1)) and

  s=B0L, z=L, p=bB0, A0=2^(B0L), A=2^(pL),
  Cstar=Dg*(1+(kappa_G+1)B0).                            (18)

Check the actual integer conditions (16). They hold eventually at every
integer N: v=R_G*B0L+1<=pL, and the explicit fixed-R_G logarithmic
growth argument in proof.md Section4 applies to every fixed p.
Pad the genuine degree to Cstar by repetitions. Evaluation uses the
ordinary map, field increments, walk and projection in 2^O(pL) time.
The 2^40 padding factor is constant. All A=2^(pL) left evaluations still
cost 2^O(pL), polynomial in actual N since 2^L<=2(N+1).
No random tables or graph certificates are inputs.

## 5. Extractor-to-sampler conversion and graph amplification

### 5.1 Basic conversion with the correct entropy direction

Let E:{0,1}^n x [A0] -> {0,1}^d extract with error epsilon0 from every
source of min-entropy at least K_E, with an independent uniform index.
Strongness is not needed here. For 0<delta<1 it is an
(epsilon0,delta)-sampler at any threshold

  k_S>=K_E+ceil(log(1/delta)):                                (20)

for every such source X and every test set B of density <=epsilon0,

  P_X[ |{j:E(X,j) in B}|/A0>2epsilon0 ]<=delta.              (21)

If the bad-input event had probability p_bad>delta, conditioning X on it
would leave entropy at least H_inf(X)+log p_bad>=K_E, because each
conditional atom is at most its old atom divided by p_bad. Every point
of that conditional source has acceptance probability >2epsilon0. But
extraction bounds the average acceptance by density(B)+epsilon0<=2epsilon0,
a contradiction. Thus sampler entropy increases by the positive overhead
in (20). No inequality with the opposite sign is being invoked. Source
linearity survives because the evaluated map is unchanged.

### 5.2 Amplification, including integer thresholds

Suppose S0 obeys (21) and a fixed map Gamma:[A] x [C] -> [A0] has

  |Gamma(Q,[C])|>2epsilon0*A0 whenever |Q|>=Lambda.

For 0<epsilon1<=epsilon0 and an integer Lambda>=1 require

  Lambda-1<=2epsilon1*A.                                    (22)

Define Samp(x,alpha,j)=S0(x,Gamma(alpha,j)). For every k_S-source X
and every B of density at most epsilon1,

  P_X[ |{alpha: Samp(X,alpha,j) in B for every j}|>
                                      2epsilon1*A ]<=delta. (23)

Proof. Outside the exceptional delta mass of source points, the old bad
indices R_x={u:S0(x,u) in B} have size <=2epsilon0*A0. If Lambda or
more new indices had their entire neighborhoods in R_x, the expansion
hypothesis would give a neighborhood larger than R_x. Hence there are at
most Lambda-1 such new indices, proving (23). Gamma has no randomness in
this argument. Its construction and uniform evaluator, not just its
existence, are required. As it only changes the old seed index, source
linearity is preserved regardless of whether Gamma or its condenser is
linear.

For (G) take Lambda=2^-20*A^(1/2) and epsilon1=2^-20*A^(-1/2).
Then Lambda is integral and Lambda-1<Lambda<2epsilon1*A. On (18), epsilon1<=2^-L
because p/2>=1 (in fact b>=4 and p>=4800). The all-bad new-index count
is strictly less than Lambda, stronger than the allowance 2Lambda.

### 5.3 Concrete base sampler from the finite primitive

At actual input length n=N choose the seeded-linear parameters

  e=L, ell=128uL, d=128*2^h*uL, L>=h+6, K_0<=N.             (24)

For L=ceil(log(N+1)), ceil(log(4N))<=L+2. Thus a_E<=2L,
rho_E<=5L, S<=50uL, and (1) gives

  D<=1200u^2L=s, K_0<=3d.                                   (25)

For example the seed bound follows from
6ell+(h+7)S+18a_E+12<=768uL+50u(u+6)L+36L+12<=1200u^2L
for h>=11 and L>=h+6; the entropy bound follows from
2*2^h(178uL+8L+4)<=384*2^h*uL=3d. Use the exact recurrence to
check K_0<=N, not merely its upper bound. Pad the seed to s bits by
ignoring the suffix. Set delta=1/100; since ceil(log 100)=7, (20) yields

  k_S=K_0+7<=3d+7, epsilon0=2^-L, A0=2^s.                   (26)

Combine with (18) and (23), setting C=Cstar. This is a uniform linear
somewhere sampler with output length d, error-density epsilon1=2^-20*A^(-1/2),
failure probability 1/100, and no source-entropy change during amplification.
The length d must still satisfy the actual affine breaker's seed guards;
(24) does not prove simultaneous compatibility.

## 6. Conditional affine breaking and the simultaneous bad-seed lemma

### 6.1 Conditional/average-entropy lift with arbitrary copies

Suppose (4) holds at degree r, entropy K, error eta. Let Z be finite,
and replace unconditional independence by the vector Markov condition

  (A,A_1,...,A_r) -- Z -- (B,B_1,...,B_r,Y,Y_1,...,Y_r).

Assume (Y,Z)=(U_d,Z), meaning Y is uniform independently of Z, and
Havg_inf(A|Z)>=K+log(1/delta). With the same fixed advice, the joint
distance in (4), now retaining Z as well, is at most eta+delta.

Proof. By (6), except on z of total mass <=delta, the conditional law of
A has min-entropy >=K. For each such z, conditional independence of the
complete source tuples, uniformity of Y, and every additive identity and
advice inequality give a legal instance of (4). Apply it with error eta.
Use the trivial distance one on exceptional z and average by (5), obtaining
eta*(1-P[bad])+P[bad]<=eta+delta. The arbitrary source-copy quantifiers
are unchanged. Choosing eta=delta=e_F/2 yields joint error e_F at entropy
K+log(2/e_F). This step uses no new evaluated algorithm.

### 6.2 Bad seeds and witness consistency

Fix a source A on n bits; fixed differing advice alpha,alpha_1,...,alpha_r;
and deterministic maps L_0,...,L_r:{0,1}^n -> {0,1}^d. Suppose their
joint leakage Z=(L_0(A),...,L_r(A)) has range size at most 2^ell and

  H_inf(A)>=K+ell+log(1/delta).                              (27)

For fixed b,w,w_1,...,w_r, define random outputs over A only by

  R_0(b,w)=F(A+b,w+L_0(A),alpha),
  R_i(b,w_i)=F(A+b,w_i+L_i(A),alpha_i).

For 0<theta<=1, call w bad if some b,w_1,...,w_r give

  TV((R_0(b,w),R_1(b,w_1),...,R_r(b,w_r)),
     (U_m,R_1(b,w_1),...,R_r(b,w_r)))>theta.                 (28)

The density of bad seeds among all d-bit strings is at most

  (eta+delta)/theta.                                        (29)

Proof. Every witness space is finite. For each bad w select one witness
(g(w),f_1(w),...,f_r(w)) by a fixed ordering; define these functions
arbitrarily on other w. Take W uniform independent of A and put

  B=B_i=g(W), A_i=A,
  Y=W+L_0(A), Y_i=f_i(W)+L_i(A),
  X=A+B, X_i=A_i+B_i, Z=(L_0(A),...,L_r(A)).

The use of f_i(W), not W, in Y_i is essential: these are exactly the
chosen tampered unshifted seeds of (28). Since Z is a function of A,
conditional on Z the complete left tuple is independent of W, and the
complete right tuple is a function of W and fixed leakage values. Thus
the vector Markov condition holds. Given Z, Y is a fixed translate of
uniform W, so (Y,Z)=(U_d,Z). Equations (6) and (27) supply the requisite
average conditional entropy.

Apply Section 6.1. The retained pair (Y,Z) determines W=Y+L_0(A), because
L_0(A) is a coordinate of Z. Project the guarantee to the outputs and W;
this is legitimate on both real and ideal laws because W is determined
from retained side information. Its TV is at most eta+delta. Given W=w,
A retains its original unconditional law, because A and W were independent.
Consequently (5) identifies the projected TV with the average, over w,
of exactly the local distances selected in (28). Every bad w contributes
more than theta, proving (29). The witnesses are proof devices, not evaluator inputs or subroutines. Taking
ell=(r+1)d always suffices; an actual smaller leakage range may be used.

## 8. Selected outer use

The source-linear graph/sampler and conditional breaking above are used
by proof.md's tagged affine-permutation outer assembly. That file supplies the
modified advice, simultaneous bad-seed accounting and four-moment proof.
The unused ordinary four-moment outer evaluator and affine consequence
have been omitted; no reference to them is a mathematical dependency.

Attribution: extractor-to-sampler conversion is due to Zuckerman (1997);
the low-rate-condense/high-rate-walk architecture to Zuckerman, Theory of
Computing 3 (2007), Sections 2--3; the row-XOR/majority architecture to
Chattopadhyay--Liao, arXiv:2110.12652v1, Section 5. The affine interface
is compared to Li, ECCC TR23-023 revision 1. The finite proofs here are
independently written, not a claim of authorship of those constructions.

## 0. Atomic proof map and selected construction

This map describes proof.md and its companions. Legacy component names refer
to their included bodies, not external files. Source URLs, scratch prefixes
and component hashes give provenance, not mathematical premises.

- linear-and-bilinear-primitives.md supplies seeded-linear.md Sections 1--6,8,
  near-full-seeded.md1--4, weak-seed-hashing.md1--4 and bilinear-extraction.md.
  The finite weaker linear primitive is selected, without a GUV/CGL premise.
- sum-product-and-ordinary-condensation.md supplies finite-approximate-ring.md
  Sections1--5,7 and ordinary-condenser.md, whose6--7 prove scalar sum-product.
- restricted-nonmalleability-and-breaking.md supplies prime-fields-and-encoding.md
  1--3, cz-conditional-entropy.md1--8 and constant-advice-breaker.md1--8.
  The ring/CZ guard is gamma=1/10000, eta=1/10000000. CZ3 is elementary,
  so ordinary condensation's use of it is noncircular.
- affine-and-standard-wrappers.md supplies compact-affine-wrappers.md. Its2
  gives the used m<=r affine transformation,3/5 give the advice/equal-length
  reductions,3.4 gives the CZ/IP reduction,6 gives the general-input proof.
  Only weak-seed hashing is selected. Constant-advice breaking is used only
  inside the fixed-degree inner affine breaker.
- sampler-and-disperser.md retains the finite real-entropy graph, sampler,
  seed-retention and conditional breaking. proof.md now gives its tagged
  four-moment affine-permutation assembly. The original ordinary-majority
  outer alternative is omitted.
- proof.md retains minimal-circuit-bridge.md and breaker-parameter-assembly-v2.md,
  whose exact degree-nine contract(F) supplies the directional outer assembly.

The selected acyclic chain is finite linear extraction -> fixed-advice
breaker -> inner affine breaker; finite ring/CZ + ordinary condensation +
weak hash -> equal-length breaker; near-full extraction -> standard breaker;
finite affine transformation -> outer degree-nine breaker; deterministic
graph/sampler + tags and four moments -> one fixed affine-permutation family; native B2 bridge
-> the eventual full 4n-O(log n) inequality. All parameter and evaluator
steps for that exact chain are specified in the final two components.

Precise legacy equation aliases: condenser-integration.md Section 5.10's P1
bilinear-prefix calculation is affine-and-standard-wrappers.md Section 5,
equation (17). Its Section 6.2 and interface (N1) are that file's Section 3.4,
equations (11)--(12), with the joint smoothed conditional-entropy law spelled
out there and in the retained CZ Sections 2,8. They are not generic ordinary
extractor assumptions. References to original preparation/audit files in the
sampler source list are provenance only; all their mathematical interfaces
are reproved in the retained sampler sections and tagged outer assembly. The broader original
circuit-reduction.md is provenance; the complete used bridge follows below.

The direct circuit proof source is GordoAR, submission 6622943f-7b2d-4f8f-b73c-979c5232484a:
https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e/submission
Attribution record:
https://api.yukon.org/api/submissions/6622943f-7b2d-4f8f-b73c-979c5232484a/attribution
Its upstream affine-system framework is Mishra--Santhanam:
https://eccc.weizmann.ac.il/report/2026/260/revision/1/download .
The circuit proof is independently written and proved below, not adopted
through the source's status or an assumed license to republish its prose.




## Paired read-once zero-period lemma

A proper read-once full-B2 formula on 2q distinct leaves, q>=1, with any
pairing of its leaves, has a nonzero zero-period on some product of q
allowed lines in the paired three-corner domain D. A zero-period v means
F(u+v)=F(u) everywhere on that flat. Flip input bits so every pair forbids
00; this preserves properness and the assertion. A safe axis line fixes
one endpoint to 1; the other allowed line is the diagonal 01,10.

First suppose all gates are AND or XNOR (1+x+y). Both have neutral input
1. Induct on q. If a root child contains m>=1 complete pairs internally,
fix its leaves with partners outside that child to 1. Those fixings use
at most one endpoint of each crossing pair. Neutrality reduces the child
to a formula on its 2m internally paired leaves, where m<q. Induction gives
a zero-period supported only on those leaves. The other root child is
unchanged. Extend the chosen internal lines by the safe crossing lines
and arbitrary allowed lines on the other pairs. This proves this case.

Otherwise every pair crosses the root, and each child has q leaves.
If one child contains an XNOR gate, choose a leaf from each of that gate's
two child subtrees, say u,v. Fix every other leaf of the chosen root child
to 1. Leave u,v free by fixing their partners in the other root child to
1. These are consistent one-endpoint-per-pair fixings. Neutrality reduces
the chosen child to XNOR(u,v). Flipping u and v preserves it and leaves
the other child unchanged, giving the desired period. If neither child
contains XNOR, both are pure ANDs. For an AND root choose one pair's
diagonal; the whole output is identically zero. For an XNOR root choose
all pair diagonals: flipping all q free coordinates swaps the two AND
monomials, preserving their XNOR. This includes q=1 and completes induction.

Now take arbitrary proper B2 gates and suppose no such zero-period exists.
For any product gate, consider a child whose leaf set contains at most
one endpoint of every pair. Fix all those leaves to safe 1. Its output
cannot be a controlling value at that gate's pin: otherwise the entire
other child is killed. Choose a leaf z in that other child, fix its mate
to 1, and complete to a product flat leaving z free. This is consistent
with the earlier fixings. Flipping z would preserve the whole formula,
a contradiction. Thus this child's all-ones output is noncontrolling.

Complement internal signals as needed so every signal takes value 1 at
the global all-ones input, absorbing complements in parent truth tables.
Leaves stay unchanged; complementing the root preserves zero-periods.
Every affine proper gate becomes XNOR. Every product gate whose two child
leaf sets each contain at most one endpoint per pair becomes AND: both
pins have noncontrolling value 1, and its output at 11 is 1.
Choose a minimal subtree V containing a complete pair. Every proper
subtree of V has at most one endpoint per pair. Consequently all gates
within V are AND or XNOR. Fix V's leaves with mates outside V to 1;
neutrality leaves a nonempty formula on its internally paired leaves.
Apply the auxiliary induction. Its period changes only those internal
leaves, leaving the rest of the full formula unchanged. Extend the flat
by the safe crossing lines and arbitrary remaining pair lines. This
contradicts the supposition and proves the lemma. The requirement that
ALL leaves be paired is essential to this proof.

## Tight incidence tree and terminal exclusion

Use the representation and incidence graph of proof.md Section 7. Suppose
at its terminal state a=k, q=k, M=2k, s=M-1. All M parameters belong to the
output component. The component bound s_0>=M-1 forces s_0=s. There are M+s
variable vertices and s object vertices, with at most 3s incidences.
Connectedness and M=s+1 force exactly 3s incidences and a bipartite tree;
every object has three distinct variable neighbors.

The output is a gate output. Indeed an affine output on M>=2 parameters
has a nonzero kernel direction v. For any direction v the domain D has
a safe a-flat parallel to v: in each pair choose the allowed line parallel
to its nonzero projection (three corners supply all three directions),
or any allowed line if the projection is zero; keep free coordinates.
An affine output would therefore have zero derivative on such a flat,
contradicting the separation hypothesis of Section 7.

Every auxiliary variable has degree at least two. If it had degree one,
its sole incidence would be a row, since it is not the output. Delete
that row and its uniquely supported auxiliary column. Invertibility of
the remaining block follows by expansion on that column. The represented
function, parameters, Phi and D are unchanged, while s drops to M-2,
contradicting the same component bound. Every nonoutput gate variable has
degree at least two (its defining gate and a consumer); dead gates have
already been removed. The output gate variable has degree one. Parameters
have degree at least one. Their total degree is therefore at least
M+2(s-1)+1=3s, so equality forces every parameter to have degree one and
every nonoutput gate/auxiliary variable to have degree two.

Root this incidence tree at the output variable. Each object has one
rootward and two outward variables. For a row the rootward variable is
auxiliary: parameters are leaves and gate variables never occur in rows.
Solve the row for that auxiliary, obtaining an XOR or XNOR of its two
outward variables. For a gate the rootward variable is its designated
output: the native directed path through consumers to the final output
exists by dead-gate removal and must be the unique incidence path.
Thus no product gate is inverted. Orienting the tree computes a proper
full-B2 read-once formula from all M parameters, each used exactly once.
It agrees with the represented function (even off D). The paired lemma
supplies a nonzero zero-period on a safe k-flat in D. Injective affine Phi
maps it to a nonzero original direction and a k-flat, contradicting
separation. The alleged tight terminal state is impossible.
