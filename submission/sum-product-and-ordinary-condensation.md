# sum-product-and-ordinary-condensation.md

Component references resolve through proof.md Section 0.

# Explicit finite approximate rings and four-coordinate expansion

Status: complete original quantitative proofs of the finite lemmas below,
including E_p(256/20481) from cz-conditional-entropy.md for every prime
p >= 2^3200000. These proofs do not use CZ Lemma 7.3, Bourgain's filling
lemma, or an unspecified absolute constant. The surrounding condenser arguments
require their separate audits; this file claims no external review.

Attribution: the ring/stabilizer framework builds on the
self-improving approximate-ring construction of Terence Tao, *The
sum-product phenomenon in arbitrary rings*, Contributions to Discrete
Mathematics 4(2), 2009, pp. 59--82, especially Lemma 1.4, Proposition 3.1,
and Theorem 1.2. Official source:
https://cdm.ucalgary.ca/article/view/61994
The sumset tools below are Ruzsa's inequalities and Petridis's minimal-
growth argument. The explicit powers, normalization, proofs and application
below are independently supplied here. We work only in commutative rings
with identity. No inclusion AB+AC subset A(B+C) is used: that inclusion is
false in general. Scalar distributions and every product of sets below
are distinguished explicitly.

Notation: nA denotes an n-fold sumset, with 0A={0}; AA={ab:a,b in A};
aA is multiplication by the single element a. A unit is an invertible
ring element. All sets used in cardinality estimates are finite and
nonempty unless indicated otherwise. Products in F_p^d are coordinatewise.

## 1. Fully quantified additive tools

### 1.1 Ruzsa triangle and covering

For finite nonempty subsets X,Y,Z of an abelian group,

  |X-Z| |Y| <= |X-Y| |Y-Z|.                              (1)

For each w in X-Z choose one representation w=x_w-z_w. The map
(w,y) -> (x_w-y,y-z_w) is injective: summing the two coordinates recovers
w and the chosen x_w then recovers y. This proves (1).

If |X+Y| <= M|Y|, there is T subset X with |T|<=M and

  X subset Y-Y+T.                                        (2)

Choose a maximal collection of disjoint translates t+Y, t in X. Their
number is at most M. Every other x+Y intersects one of them, which gives
x in t+Y-Y. The same conclusion follows from |X-Y|<=M|Y| by using -Y.

### 1.2 Explicit iterated sumset estimate

If |U+V|<=K|U| for finite nonempty U,V and K>=1, then for all nonnegative
integers m,n with m+n>0,

  |mV-nV| <= K^(m+n) |U|.                                (3)

Here is a proof rather than an invocation of Plunnecke's theorem. Choose
nonempty X subset U minimizing k=|X+V|/|X|. Thus k<=K and |Y+V|>=k|Y|
for every Y subset X, including Y empty. For any finite set C, order its
elements and add the translates X+c successively. For the next c, let
Y={x in X:x+c belongs to the previous union X+C_old}. Then Y+V+c is
contained in the previous union X+V+C_old. Consequently the newly added
part of X+V+c has size at most |X+V|-|Y+V|<=k(|X|-|Y|). The latter
parenthesis is exactly the newly added size of X+c. Summing proves

  |X+V+C| <= k |X+C|.

Induction gives |X+jV|<=k^j|X| for every j>=0. Apply (1) to mV,-X,nV
to obtain |mV-nV|<=k^(m+n)|X|<=K^(m+n)|U|. This includes m=0 or n=0.

In particular, |A+A|<=K|A| implies

  |mA-nA| <= K^(m+n)|A|.                                 (4)

If P contains a set C of size N and |P-P|<=L N, then N<=|P|<=L N,
and (3) with U=P,V=-P yields

  |mP-nP| <= L^(m+n+1)N.                                 (5)

The hypotheses of (5) include C subset P. The bound |P|<=|P-P| follows
by fixing one element of P. Interchanging m and n in (3) gives the stated
signs. The harmless extra power in (5) avoids any unrecorded constants.

### 1.3 A sharper form of (5)

If C subset P, |C|=N and |P-P|<=LN, then, for k=m+n>=1,

    |mP-nP| <= L^(m+n) N.                                (5s)

Indeed put Q=|P|>=N and k0=|P-P|/Q<=LN/Q. Applying (3) with
U=P,V=-P, and interchanging the signs if necessary, gives
|mP-nP|<=k0^k Q<=(LN)^k/Q^(k-1)<=L^k N. All sets are nonempty;
k0>=1 since a translate of P lies in P-P. This improves (5) without
changing its hypotheses. The older bound (5) remains a valid corollary.

## 3g. A minus stabilizer and direct Ruzsa unit actions

Let C0 be a nonempty set of units,|C0|=M,K>=2,with
 |C0+C0|,|C0C0|<=KM,
 |(C0-C0) outside units|<M/(32768K^17).                 (Mi0)
Then there areb0 in C0,B0 subsetC0 and a finiteunitalringT such that
 |B0|>=M/(2K), B0 subsetb0T, |T|<=32768K^18M.           (Mi1)
Multiplicativeenergy andCauchy--Schwarz give
sum_(a,b in C0)|aC0 intersectbC0|>=M^3/K. Every summand is at mostM.
Pairs with intersection at leastM/(2K) number at leastM^2/(2K).
Chooseb0 with at leastM/(2K) suchpartnersB0,includingb0.
Normalize C=b0^(-1)C0,G=b0^(-1)B0. All hypotheses are preserved;
1 in G subsetC,|G|>=M/(2K),and|gC intersectC|>=M/(2K) foreachg inG.
ForH=gC intersectC,triangle(1)with middle-H gives
 |gC-C|<=|gC+H||C+H|/|H|<=2K^3M.
Covering(2) in BOTH directions therefore gives
 gC subsetC-C+X_g, C subsetgC-gC+Y_g,
 |X_g|,|Y_g|<=2K^3.                                   (Mi2)
Here|gC|=M sinceg is a unit;|C-gC|=|gC-C| supplies the reversecover.

DefineS={x:|xC-C|<=8192K^15M}. Collisionsofxg1-g2 withgi inG
number at least|G|^4/(8192K^15M). For fixedd=g1-g1',at most|G|^2
representations contribute. Thus at leastM/(32768K^17) distinctd
inG-G satisfyxd inG-G. By(Mi0),one is a UNITd.
Write d=g1-g2,xd=g3-g4,and putD=g1C uniong2C uniong3C uniong4C.
By(Mi2),|(-C)+D|<=8K^3M. Apply the FULL two-set estimate(3)
withU=-C,V=D,m=n=2 to obtain
 |2D-2D|<=(8K^3)^4M=4096K^12M.
Matched productsxdC-dC are contained ing3C-g4C-g1C+g2C,which is
contained in2D-2D. Fixed-unit dilation consequently proves
 |xC-C|=|xdC-dC|<=4096K^12M.                           (Mi3)
No product-of-sets distributivity is used;D need not have small doubling.
RepeatingG-collisioncount at(Mi3) gives at leastM/(16384K^14)
differences. Fewer thanM/(32768K^17)<=M/(32768K^14) are nonunits,
leaving more thanM/(32768K^14) unitchoices. Each pair(d,xd)withd a
unit determinesx. The pairset has size at most|C-C|^2<=K^4M^2,
soS is finite; doublecounting gives|S|<=32768K^18M.
Foranyfinite nonemptyA,the additive extension inequality is
 |A+C-C|<=K^2|A-C|.                                   (Mi4)
To proveit,choose nonemptyP subset-C minimizingk=|P-C|/|P|<=K.
ThePetridisextension argument in theproof of(3),withU=-C,V=-C,
gives|A+P-C|<=k|A+P|. Triangle(1)atX=A-C,Y=-P,Z=-C gives
|A+C-C|<=|A-C+P||P-C|/|P|<=k^2|A+P|<=K^2|A-C|.
HereP subset-C impliesA+P subsetA-C. No cardinality ofxC isassumed.

Translationby1 preservesS: (x+1)C-C subsetxC+C-C,so(Mi4)withA=xC
and(Mi3) bounditby4096K^14M<=8192K^15M. Forunitg inG,
triangle(1)atX=gxC,Y=gC,Z=C gives
 |gxC-C|<=|gxC-gC||gC-C|/M
          <=4096K^12M*2K^3M/M=8192K^15M.
Onlyg,notx,isrequiredtobeaunit. HencegS subsetS.
Multiplicationbyg andtranslationby1 are ambientbijections restricting
to permutations of finiteS; their inverses preserveS. EveryG-monomialm,
includingemptyproduct1,is a unit andits multiplicationM_m preservesS
in both directions. Conjugating translationby1 givestranslationbym.
ThusS is invariant under monomialtranslations,their negatives andfinite
integerlinearcombinations. Since0 inS,the additivegroupT generatedby
G-monomials is contained inS. It contains1 andG andis closed under
products,so is a finiteunital subring. B0=b0G subsetb0T proves(Mi1).
The ring claim concernsT,not arbitraryclosure ofS.

## 4. Coordinate-fiber bound on bad differences

Let C subset F_p^d be nonempty, |C+C|<=K|C|, and let pi_j denote the
j-th coordinate projection. Put

  Z_j=(C-C) intersect ker(pi_j).

Choose one representative c from C for each value of pi_j(C). Their
translates c+Z_j are disjoint, so

  |pi_j(C)| |Z_j| <= |C+Z_j| <= |2C-C| <= K^3|C|.         (17)

A vector in F_p^d is a nonunit exactly when at least one coordinate is
zero. Taking a union bound in (17) proves

  |(C-C) outside the units|
      <= K^3 |C| sum_(j=1)^d 1/|pi_j(C)|.                (18)

This bounds the number of distinct bad differences, not just the number
of pairs yielding them. That distinction is needed in Section3g.

## 5. Four-coordinate hypotheses

Letp beprimeandA a nonempty subsetofF_p^4. The hypotheses are:

(a) |A|<p^(7/3-1/16);
(b) at least half of A consists of units;
(c) every D subset A with |D|>=p^(-1/16)|A| has both coordinate-pair
    projections pi_12(D), pi_34(D) of size at least p^(5/4);
(d) no such D is contained in either kind of plane

    {(x,y,c1*x,c2*y):x,y in F_p},
    {(x,y,c2*y,c1*x):x,y in F_p},

    for any c1,c2 in F_p.

## 5g. Four-coordinate expansion from the minus stabilizer

Under exactly Section5's hypotheses(a)--(d),forp>=2^3200000,
 max(|A+A|,|AA|)>p^(256/20481)|A|.                     (19m)
Putgamma=256/20481 and assume failure. WriteN=|A|,C0=A intersectunits,
M=|C0|,K=2p^gamma. ThenM>=N/2,|C0+C0|,|C0C0|<=KM.
Sincep^(1/16)>2,C0 qualifies for(c),so everyindividualcoordinate
projection has at leastp^(1/4)values. Equation(18) gives
 Z=|(C0-C0) outsideunits|<=4K^3M p^(-1/4).
Now1/4-20gamma=1/81924 and
 131072K^20p^(-1/4)=2^37p^(-1/81924)<1,               (21m)
since37*81924=3031188<3200000. ThusZ<M/(32768K^17).
Section3g givesb0,B0,T with
 |B0|>=M/(2K)>=N/(8p^gamma),B0 subsetb0T,|T|<=32768K^18M.
The subsetB0 qualifies for(c)--(d),sincep^(1/16-gamma)>8.
Also32768K^18M<=2^33p^(18gamma)N<p^(1/4)N,
because18*256/20481+33/3200000<1/4. Consequently
 |T|<p^(7/3-1/16+1/4)=p^(121/48)<p^3.
An additive subgroupofF_p^4 is anF_p-vector subspace,so dimT<=2.
Unit scaling preserves pairedprojectioncardinalities. By(c),each
projection ofT has size at leastp^(5/4),hence dimensionat leasttwo.
ThusdimT=2 and both projections bijectontoF_p^2. T is thegraph of
a unitalF_p-linear ringisomorphismphi:F_p^2->F_p^2. Imagesof(1,0),(0,1)
are nonzero complementaryidempotents,so are thatpair in eitherorder.
F_p-linearity makesphi identityorcoordinateswap. Forb0=(b1,b2,b3,b4),
allbi nonzero,b0T is therefore oneof
 {(x,y,(b3/b1)x,(b4/b2)y):x,y inF_p},
 {(x,y,(b3/b2)y,(b4/b1)x):x,y inF_p}.
These are the forbiddenplanes(d),yet qualifyingB0 liesinb0T.
This contradiction proves(19m)withallsubsetquantifiers.

The selected family uses eta=256000000/2109543000001, gamma=256/20481.
The CZ cap eta<=1/8000 holds; 103eta<gamma follows from
2109543*1000000<2109543000001. Theta=eta/8=32000000/2109543000001.
The separate explicit guard u>=3200000
in the advice map ensures actual p>2^u meets this theorem's threshold.
The weaker eta*u>=64 guard is not used to infer that prime threshold.
AllremainingfiniteCZ/uniform-sourceconditions retain their suppliedproofs.
The new threshold is proved directly,not inferredfrom a weaker theorem.

# Scalar sum-product and one support-growth step

Section 6 proves the scalar sum-product estimate (30') used by Lemmas RG
and S1' in Section 5. No somewhere-condenser or scalar entropy-growth
theorem is assumed. All
logarithms are binary, TV is half-L1, and H_inf is min-entropy. Statements
with real entropy thresholds mean literal atom bounds, so integer ceilings are
not silently omitted.

## 5. One support-growth step

This section replaces the original reduction to scalar entropy growth. It
uses the two branches of (30') in Section 6.3, with weighted cost
16*2+18*2=68; the sumset inequality (6) and the
three-edge-walk lemma (7) of restricted-nonmalleability-and-breaking.md
Sections 3.1 and 3.2, cited below as (R6) and (R7) (not this file's (6) and
(7)), used here in (F_p,+) and (F_p^*,x); and prime-fields Sections 1-2,
recalled in Section 5.2. U_D is uniform on D.

Constants: e=eta_RG=10/711, xi=1/427, r=(1-2xi)/(1+2xi)=425/429.
The two separate branch margins and the support margin are
 1-71e=1/711, r-34e=52105/101673, e/2-3xi=2/303597.      (RG0)
In particular 800000/711>252, 800000*52105/101673>122, and
300000000*(2/303597)>17+60/711. These are rational comparisons.

### 5.1 Lemma RG

**Lemma RG.** Let p be prime and M an integer with log M>=800000 and
p>M^(2/(1+2xi)). Let A1,A2,C0 be subsets of F_p^* of size exactly M and
X1,X2,X3 independent and uniform on them. Then for every W subset F_p with
|W|<=4M,

    P[X1*X2+X3 in W] <= 2M^(-e).                                      (RG)

Proof. Call a nonempty D subset F_p^* plus-good if P[U_D+U_V in W]<=M^(-e)
for every V subset F_p^* with |V|=M and every W subset F_p with |W|<=4M,
and times-good if P[U_D*U_V in W]<=M^(-e) for all such V,W. A disjoint
union of times-good sets is times-good: for fixed V,W its hitting
probability is the size-weighted average of the pieces'. Start with R=A1
and repeat: stop if R is plus-good or |R|<=M^(1-e); otherwise find a
times-good D subset R (below) and move D out of R. Each round removes a
point, so the procedure stops; it is a proof device, never evaluated.

Finding D. R is not plus-good: some V,W as above have
P[U_R+U_V in W]>alpha:=M^(-e), so the bipartite graph on R x V with edges
{(x,v):x+v in W} has more than alpha|R|M edges. (R7) gives D subset R and
C' subset V with |D|>=alpha|R|/8, |C'|>=alpha M/8, each (x,v) in D x C'
joined by at least alpha^5|R|M/2^18 walks x,v1,x1,v. A walk gives
x+v=s1-s2+s3 with s1=x+v1, s2=x1+v1, s3=x1+v in W; given (x,v) the triple
fixes v1=s1-x and x1=s2-v1, and the triple fixes the sum. So distinct sums
use disjoint sets of triples, and |D+C'|alpha^5|R|M/2^18<=|W|^3<=64M^3.

Suppose D is not times-good: some V',W' have P[U_D*U_V' in W']>alpha. The
edges {(x,v):xv in W'} have products in W' minus 0, a subset of F_p^*, and
the same counting in (F_p^*,x) (xv=s1*s2^(-1)*s3) gives D' subset D and
B' subset V' with |D'|>=alpha|D|/8, |B'|>=alpha M/8 and
|D'B'|alpha^5|D|M/2^18<=64M^3. As |R|>M^(1-e), |D|>M^(1-2e)/8 and
|D'|>M^(1-3e)/64, (R6) additively (with D' subset D) and multiplicatively
give, for K+=|D'+D'|/|D'| and Kx=|D'D'|/|D'|,

    K+<=|D+C'|^2/(|C'||D'|)<=2^54 M^(3+12e)/(|R|^2|D|)<2^57 M^(16e),
    Kx<=|D'B'|^2/(|B'||D'|)<=2^54 M^(3+12e)/|D|^3<2^63 M^(18e).

N'=|D'|>=2: (1-3e)log M>=800000*(681/711)>7, and
N'>M^(1-3e)/64. Also N'<=M and p/N'>=p/M>M^r.
Put X=K+*Kx<2^120 M^(34e) and apply the two branches of Section 6.3.
If p>=3XN', then N'<=33X^2<33*2^240 M^(68e), so
 M^(1-71e)<64*33*2^240<2^252,
contradicting (RG0). If p<3XN', then
 M^r<p/N'<3X<3*2^120 M^(34e)<2^122 M^(34e),
contradicting the second margin in (RG0). This treats both branches
without weakening the lower bound on N' by a factor r in its exponent.
Therefore D is times-good.

Distribution. Let G be the times-good union and R the final remainder.
Given X1 in G and X3=c0, P[X1X2 in W-c0]<=M^(-e) (|W-c0|=|W|); average over
c0. Given X1 in R with R plus-good, fix X2=a2: X1a2+X3=a2(X1+a2^(-1)X3),
where a2^(-1)C0 is a size-M subset of F_p^* and |a2^(-1)W|=|W|, so the bound
is again M^(-e); average over a2. If instead |R|<=M^(1-e), the R branch
has probability at most M^(-e). Empty branches contribute zero. The
G contribution is <=M^(-e); the remainder contributes <=M^(-e),
proving (RG).

### 5.2 The one-step map and its support lemma

For an integer m>=1 let p_m be the least prime above 2^m; prime-fields
Sections 1-2 (restricted-nonmalleability-and-breaking.md, first part) prove
2^m<p_m<2^(m+1) and find p_m by trial division in O(2^(3m/2+2)(m+1)^2) bit
operations. Read an m-bit string as a binary integer in [0,2^m), hence as an
element of F_(p_m) (injective, since 2^m<p_m). Define

    E_m(x1,x2,x3) = (x1*x2+x3 mod p_m) mod 2^m,                        (S2)

the low m bits of the canonical representative in {0,...,p_m-1}. For alpha
in {0,1}^(3m) write alpha=(alpha1,alpha2,alpha3) with m-bit blocks; its four
ROWS are alpha1, alpha2, alpha3 and E_m(alpha).

**Lemma S1' (support growth).** Let m>=300000000. For every Q subset
{0,1}^(3m) with |Q|>=2^(3m/2-15), the union U of the four row images of Q
(a subset of {0,1}^m) satisfies |U|>2^((1/2+xi)m).

Proof. Suppose N=|U|<=2^((1/2+xi)m). As Q is a subset of U^3,
N>=2^(m/2-5), so M:=|U\{0}|>=N/2>=2^(m/2-6) and log M>=m/2-6>=149999994
>=800000. Also M<=N<p_m^(1/2+xi) (as 2^m<p_m), so p_m>M^(2/(1+2xi)). Let
X1,X2,X3 be independent and uniform on U, viewed in F_(p_m). Lower bound: on
Q the fourth row lies in U, so P[E_m(X) in U]>=|Q|/N^3>=2^(-15-3xi*m).
Upper bound: E_m(X) lies in U exactly when Y=X1X2+X3 lies in
U*={y<p_m: y mod 2^m in U}, and |U*|<=2N since p_m<2^(m+1). P[some Xi=0]
<=3/N, and conditioned on no zero the Xi are independent and uniform on
U\{0}, of size M. The hitting-probability Lemma RG with |W|=|U*|<=2N<=4M gives
(the conservative upper bound below weakens its factor2):

    P[E_m(X) in U] <= 3/N+2M^(-eta_RG) <= 3M^(-eta_RG) <= 4M^(-eta_RG),

since 3/N<=3/M<=M^(-eta_RG) (M^(1-eta_RG)>=3). So 4M^(-eta_RG)>=2^(-15-3xi*m), that
is, eta_RG*log M<=17+3xi*m, and with log M>=m/2-6: (eta_RG/2-3xi)m<=17+6eta_RG. But
(eta_RG/2-3xi)m=2m/303597>=600000000/303597>18>17+6eta_RG.
Contradiction.

Only set counting and independent uniform variables on U appear; no
independence of the rows of an actual source is used. sampler-and-
disperser.md Sections 3-4 use Lemma S1' with the explicit guard
m=P3*B0*L>=300000000, automatic at the outer L>=2^40.

## 6. Scalar sum-product from polynomial incidences

We prove every incidence estimate used here. The argument is adapted from
M. Lewko, *A Szemeredi-Trotter Theorem in Arbitrary Fields*,
arXiv:2609.27023v3, Lemma 2 and Section 4, equation (14), and Elekes's
point-line construction. No theorem from either source is a premise.
Write t_+=max(t,0); characteristic zero is denoted p=infinity.

Algebra used below. Polynomial rings over fields are unique factorization
domains. One elementary justification is induction on the variables:
over a UFD, the product of primitive polynomials is primitive (reduce
modulo each prime divisor of the coefficients; the resulting polynomial
ring over a domain has no zero divisors). Clearing denominators then
transfers factorization and its uniqueness from the fraction field's
Euclidean univariate polynomial ring. This is Gauss's lemma, and also
shows that primitive irreducibles remain irreducible on passing to the
fraction field. Degrees, including degree in a chosen group of variables,
add under multiplication of nonzero polynomials, by their leading forms.
A nonzero univariate polynomial of degree j has at most j distinct roots:
divide successively by the distinct linear factors.

### 6.1 Lemma M: multiplicities on parallel lines

Let K be a field, 0<=d<p an integer, and F in K[x,y] nonzero and squarefree
(no repeated irreducible factors), with total degree <=d. Let lambda_i be
finitely many distinct parallel lines, with q_i on lambda_i, such that
F restricts to a nonzero polynomial on every lambda_i. For a nonzero common
direction v let mu_i be the order of u in F(q_i+uv). Then

    sum_i (mu_i-1)_+ <= d(d-1).                              (M)

Proof. Empty collections and constant F are immediate. Interchange x,y
if the lines are vertical, then shear X=x,Y=y-a*x, where a is their common
slope. The lines become Y=b_i with distinct b_i, and q_i=(xi_i,b_i).
The invertible substitution preserves degree and squarefreeness. Factor
the transformed F as A(Y)G(X,Y), removing all factors independent of X.
Since no restriction is zero, A(b_i)!=0; multiplicities are unchanged.
If G is constant, all mu_i=0. Otherwise each irreducible factor f of G has
1<=deg_X f<=d<p, so f_X is nonzero of smaller X-degree. Write G=fQ.
Squarefreeness gives f not dividing Q; from G_X=f_X Q+f Q_X, primality
shows f cannot divide G_X. Thus G_X!=0 and gcd(G,G_X)=1.

Put h_i=(mu_i-1)_+, d0=deg G and d1=deg G_X<=d0-1. Both
G(xi_i+u,b_i) and G_X(xi_i+u,b_i) are divisible by u^h_i, since
(u^mu g)'=mu*u^(mu-1)g+u^mu g'. Let Pi_T be the space of polynomials of
total degree <=T, of dimension binom(T+2,2). For an integer T sufficiently
large, in particular T>=d0+d1, let W_T be the common kernel of the
sum h_i coefficient functionals [u^j]H(xi_i+u,b_i), 0<=j<h_i.
They are independent: for each (i,j) the polynomial

    (X-xi_i)^j product_(v!=i) (Y-b_v)/(b_i-b_v)

has just that prescribed coefficient equal to one and all others zero.
Take T also at least all their degrees. Therefore codim W_T=sum h_i.
Every multiple of G or G_X lies in this kernel, so
 J_T=G Pi_(T-d0)+G_X Pi_(T-d1) is contained in W_T.
Coprimality and degree additivity give the intersection of these summands
as G G_X Pi_(T-d0-d1). Multiplication is injective. Hence

 codim J_T=binom(T+2,2)-binom(T-d0+2,2)-binom(T-d1+2,2)
             +binom(T-d0-d1+2,2)=d0*d1.

Thus sum h_i=codim W_T<=codim J_T=d0*d1<=d(d-1). This includes d1=0.
The restriction d<p is essential: x^p-y has zero x-derivative and order p
on every horizontal line through (xi,xi^p).

### 6.2 Lemma I: the d-form incidence bound

Let P be m points of k^2 and L be n distinct nonvertical lines over a field
k of characteristic p. For every integer 1<=d<p put v=floor(2n/(d+2)). Then

 I(P,L)<=m(v+1)+d(d-1)+n<=m+n+2mn/(d+2)+d(d-1).            (I)

Proof. Empty P or L is immediate. Seek S=sum_(j=0)^v z^j S_j(x,y), with
deg S_j<=d, satisfying S(x,ax+b,a)=0 as a POLYNOMIAL IDENTITY in k[x]
for each y=ax+b in L. There are (v+1)binom(d+2,2)>n(d+1) coefficients and
at most n(d+1) homogeneous linear conditions, so a nonzero S exists.
Replace S by the product of its distinct irreducible factors in k[x,y,z].
Neither degree increases. Each old zero restriction forces a factor's
restriction to vanish, since k[x] is a domain, and that factor survives.
Thus S is squarefree, deg_(x,y) S<=d and deg_z S<=v, with all identities.

For q=(q1,q2) in P let r(q) count its incident lines and write
 R_q(u,z)=S(q1+u,q2+zu,z)=sum_j c_(q,j)(z)u^j.
A monomial x^alpha y^beta z^gamma contributes to u^j with z-degree at
most gamma+j, so deg c_(q,j)<=v+j. Every incident slope a is a root of
every c_(q,j), by the imposed identity. Distinct lines through q have
distinct slopes. Root counting gives c_(q,j)=0 when j<r(q)-v. Therefore
u^((r(q)-v)_+) divides R_q in k[z,u].

Work over K=k(z), keeping z indeterminate, and set
lambda_q={(q1+u,q2+zu):u in K}. These parallel lines are distinct:
equality of their intercepts q2-z*q1 forces equality of both coordinates.
S stays squarefree in K[x,y]. Indeed each irreducible factor involving
x or y is primitive over k[z], hence irreducible over k(z) by Gauss's
lemma. Distinct such factors f,g cannot become associates: clearing
denominators would give af=bg for nonzero a,b in k[z]; primality of f and
f not dividing b would force f dividing g. Factors in z alone become units.

Let E={q:R_q=0}, h=|E|. On the remaining lines apply (M), using
mu_q>=(r(q)-v)_+. Their incidences total at most
 (m-h)(v+1)+sum_(q notin E)(r(q)-v-1)_+
 <=(m-h)(v+1)+d(d-1).
For q in E, substitute u=x-q1 in R_q=0. Division in y by the monic
linear polynomial l_q=y-q2-z(x-q1) shows l_q divides S in k[x,y,z].
These factors are irreducible (their quotient ring is k[x,z]), pairwise
non-associate, and of z-degree one. Their product divides S, so h<=v.
A given line containing s points of E contributes s<=1+binom(s,2);
each pair of distinct points lies on at most one line. Consequently
 I(E,L)<=n+binom(h,2)<=n+h(v+1).
Adding gives (I). All zero restrictions above are polynomial identities,
not merely pointwise zeros on a finite field.

### 6.3 Elekes theorem, with both characteristic branches

For A subset F_p^*, N=|A|>=1, let Kplus=|A+A|/N, Ktimes=|AA|/N and
X=Kplus*Ktimes>=1. Then

 p>=3XN implies N<=33X^2; p<3XN implies p/N<3X;
 min(N,p/N)<=33 Kplus^2 Ktimes^2.                         (30')

Proof. Take P=(A+A) times (AA), of size m=XN^2, and the N^2 lines
L={y=a(x-b):a,b in A}. They are distinct: their slopes determine a and,
since a!=0, their intercepts determine b. For each a,b,c in A the point
(b+c,ac) lies on its line, and distinct c give distinct points on that
line. Thus I(P,L)>=N^3. If p>=3XN, choose d=ceil(3XN)-2. Because XN>=1
and p is an integer, 1<=d<=p-2 and d+2>=3XN, while d<3XN.
Applying (I) gives
 N^3<=(X+1)N^2+(2/3)N^3+9X^2N^2,
so N<=3X+3+27X^2<=33X^2, using X>=1. Otherwise p/N<3X<=33X^2,
proving both branches and (30'). No deletion of zero is hidden: A is
assumed nonzero; Lemma RG applies it to D' subset F_p^*.

## 8. Sources

Barak, Impagliazzo, Wigderson, *Extracting Randomness Using Few Independent
Sources*, Sections 3.1--3.4 and Appendix A:
https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/BIW04/biw.pdf
(the AB+C growth step and the additive/multiplicative peeling architecture);
The finite BSG and sumset tools are proved in the companion module, with
attribution there. The hitting-probability goodness and additive-first
peeling architecture adapt published b15a467b-c75f-4152-befe-2b634e34a474
and the exact-size support route of mpjunior92's f7447a6f-a996-404e-bb7f-df63d9488837.
Section 6 adapts Lewko's multiplicity and incidence proofs (v3, 7 Oct 2026):
https://arxiv.org/html/2609.27023v3 ; source license CC BY 4.0.
The final construction is due to G. Elekes, *On the number of sums and
products*, Acta Arith. 81 (1997), 365-367. The cost-68 branch calibration
and the explicit constants here are derived in this package. The preprint's
Lean formalization is not assumed. The earlier ring/CZ growth argument
is unchanged; eta_RG is distinct from its parameter eta_CZ.
