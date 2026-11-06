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

The selectedfamily useseta=1/53683,gamma=256/20481.
The directlyprovedCZcap iseta<=1/10000. Its strictgamma guard holds:
 671*20481=13742751<13742848=256*53683.
Theta=eta/9=1/483147. The CZguardeta*u>=64 forces
 u>=3435712>3200000,so actualp>2^u meets this theorem's threshold.
AllremainingfiniteCZ/uniform-sourceconditions retain their suppliedproofs.
The new threshold is proved directly,not inferredfrom a weaker theorem.

# Scalar sum-product and one support-growth step

From the original ordinary-condenser module only the scalar sum-product
estimate (30) of Section 6 is used here, by the new Lemmas RG and S1' of
Section 5; its somewhere condenser, the entropy-growth statement G and G's
proof are not part of this package. All
logarithms are binary, TV is half-L1, and H_inf is min-entropy. Statements
with real entropy thresholds mean literal atom bounds, so integer ceilings are
not silently omitted.

## 5. One support-growth step

This section replaces the original reduction to scalar entropy growth. It
uses only (30) of Section6.2: K=max(|A+A|,|AA|)/N satisfies
K^16>=2^-15 min(N,p/N); the sumset inequality (6) and the
three-edge-walk lemma (7) of restricted-nonmalleability-and-breaking.md
Sections 3.1 and 3.2, cited below as (R6) and (R7) (not this file's (6) and
(7)), used here in (F_p,+) and (F_p^*,x); and prime-fields Sections 1-2,
recalled in Section 5.2. U_D is uniform on D.

Constants: e=eta_RG=1/292,xi=1/1753,r=(1-2xi)/(1+2xi)=1751/1755.
 r(1-3e)-288e=599/512460>0, e/2-3xi=1/1023752>0.        (RG0)
These identities follow by integer multiplication; 900000*599/512460>1029.

### 5.1 Lemma RG

**Lemma RG.** Let p be prime and M an integer with log M>=900000 and
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

N'=|D'|>=2 (as (1-3e)logM>=7),D' subset F_p^*,N'<=M,
and p/N'>=p/M>M^r. Thus min(N',p/N')>=M^(r(1-3e))/64.
Both displayed growth ratios are <2^63 M^(18e), since M>=1.
By(30), 2^-21 M^(r(1-3e))<=K(D')^16<2^1008 M^(288e).
Consequently [r(1-3e)-288e]logM<1029, contradicting(RG0)
and logM>=900000. Thus D is times-good.

Distribution. Let G be the times-good union and R the final remainder.
Given X1 in G and X3=c0, P[X1X2 in W-c0]<=M^(-e) (|W-c0|=|W|); average over
c0. Given X1 in R with R plus-good, fix X2=a2: X1a2+X3=a2(X1+a2^(-1)X3),
where a2^(-1)C0 is a size-M subset of F_p^* and |a2^(-1)W|=|W|, so the bound
is again M^(-e); average over a2. If instead |R|<=M^(1-e), the R branch
has probability at most M^(-e). This proves (RG).

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
>=900000. Also M<=N<p_m^(1/2+xi) (as 2^m<p_m), so p_m>M^(2/(1+2xi)). Let
X1,X2,X3 be independent and uniform on U, viewed in F_(p_m). Lower bound: on
Q the fourth row lies in U, so P[E_m(X) in U]>=|Q|/N^3>=2^(-15-3xi*m).
Upper bound: E_m(X) lies in U exactly when Y=X1X2+X3 lies in
U*={y<p_m: y mod 2^m in U}, and |U*|<=2N since p_m<2^(m+1). P[some Xi=0]
<=3/N, and conditioned on no zero the Xi are independent and uniform on
U\{0}, of size M. The hitting-probability Lemma RG with |W|=|U*|<=2N<=4M gives
(the conservative upper bound below weakens its factor2):

    P[E_m(X) in U] <= 3/N+M^(-eta_RG)+2M^(-eta_RG) <= 4M^(-eta_RG),

since 3/N<=3/M<=M^(-eta_RG) (M^(1-eta_RG)>=3). So 4M^(-eta_RG)>=2^(-15-3xi*m), that
is, eta_RG*log M<=17+3xi*m, and with log M>=m/2-6: (eta_RG/2-3xi)m<=17+6eta_RG. But
(eta_RG/2-3xi)m=m/1023752>=300000000/1023752>25>18. Contradiction.

Only set counting and independent uniform variables on U appear; no
independence of the rows of an actual source is used. sampler-and-
disperser.md Sections 3-4 use Lemma S1' at m=P3*B0*L>=878*500*1024.

## 6. An elementary finite scalar sum-product theorem

The argument uses Garaev's prime-field pivot (R58, Lemma 3), specialized to
one set. We give its proof and a deliberately weak explicit bound. It does not
invoke Garaev's full unequal-size theorem or an unspecified implied constant.
The sumset tool used here is proved in cz-conditional-entropy.md Section 3.1:
if finite nonempty X,C in an abelian group satisfy |X+C|<=K|X|, there is a
nonempty X' subset of X with

    |X'+jC|<=K^j |X'|, hence |jC|<=K^j|X|, for every j>=1.  (23)

Indeed choose X' minimizing |Y+C|/|Y| over nonempty Y subset of X. Its
minimal-growth translate argument applies to every additional set, so induction
gives the displayed bound. This is the exact elementary lemma, not a call to
an unproved asymptotic Plunnecke estimate.

We also use the elementary difference triangle inequality

    |U-V||Z| <= |U+Z||V+Z|.                                (24)

For each d in U-V fix a representation u_d-v_d. The map
(d,z) -> (u_d+z,v_d+z) is injective, proving (24).

### 6.1 A prime-field pivot without an extra product-energy loss

LetB subsetF_p^* haveb=|B|>=2. Define
 R=(B-B)/(B-B),usingonlynonzerodenominators.
IfR isproper,thenR-1!=R: invarianceunder1wouldimplyinvarianceunder
allintegertranslationmultiples,whichfillthePRIMEfield. Chooser inR
withr-1 outsideR;write r=(x1-x2)/(y1-y2),xi,yi inB,y1!=y2.
Themap(b',c)->b'+(r-1)c onB xB isinjective. A collisionwithc!=c'
wouldputr-1inR (R=-R);equalc forcesb'=b''. Itsb^2imagepoints liein
B-B+rB. Scalingbythenonzeroy1-y2 givesasubsetof
 (x1-x2)B+(y1-y2)B-(y1-y2)B,                           (25)
ofsize at leastb^2. Every displayedset multiplication isbysinglefixed
scalars;theinclusionusesonly(r-1)B subsetrB-B.

IfR=F_p,averageovert inF_pthenumberofcollisionsofb'+tc onB xB.
Equalc contributeb^2;differentc haveatmostonetforeachorderedpair
ofb'values. Theaverageis atmostb^2+b^4/p. ThusCauchy--Schwarz
givesonetwithimagesizeatleastb^2/(1+b^2/p)>=(1/2)min(b^2,p).
Representt=(x1-x2)/(y1-y2)usingR,withy1!=y2. Scalingembedsitinside
 (x1-x2)B+(y1-y2)B.                                   (26)
InbothcasesfourcoefficientsfromBproducefourorsixsigneddilatesofB
whose sumhascardinalityatleast
 (1/2)min(b^2,p).                                     (27)
ThisdirectproofusesC=B,s=1,ratherthananunnecessarysmallerintersection
set. Itassertsnoextensionoftheprime-fieldtranslationargumentto
arbitraryextensionfields.

### 6.2 Popular overlaps and an explicit bound

Let A subset of F_p^*, N=|A|>=2, and set

    K=max(|A+A|,|AA|)/N >=1.

Multiplicative energy gives sum_(a,b in A)|aA intersect bA|>=N^3/K.
Choose b0 in A with sum_(b in A)|bA intersect b0A|>=N^2/K, and let

    B={b in A : |bA intersect b0A| >= N/(2K)}.

The complementary contribution is at most N^2/(2K), while any one summand is
at most N. Consequently |B|>=N/(2K). For each b in B put I=bA intersect b0A
and C_b=bA union b0A. Since I lies in both dilates,

    |I+C_b| <= 2|A+A| <=2KN.

Apply (23) with j=2 and base I; this gives

    |bA+b0A| <= |2C_b| <= (2KN)^2/|I| <=8K^3 N.

Equation (24) gives |bA-b0A|<=2K^3N, so both signs are bounded by 8K^3N.
We use an additive extension with no ring or unit-density hypotheses:
if finite nonempty C in an abelian group satisfies |C+C|<=K|C|, then
 |U+C-C|<=K^2|U-C| for every finite nonempty U.          (28a)
Choose nonempty P subset -C minimizing k=|P-C|/|P|<=K.
The minimal-growth extension argument of (3), at base -C and added -C,
gives |U+P-C|<=k|U+P|. Triangle(1) at X=U-C,Y=-P,Z=-C gives
|U+C-C|<=|U-C+P||P-C|/|P|<=k^2|U+P|<=K^2|U-C|.
This is the general additive proof underlying(Mi4); it uses only |C+C|.

If |B|>=2 apply6.1. Write n=x1-x2,d=y1-y2!=0, with all four
coefficients in B. Put D=x1A union x2A union y1A union y2A.
The popular minus bounds give |D-b0A|<=8K^3N. Lemma(3), at
U=-b0A,V=D,m=n=2, therefore yields
 |2D-2D|<=4096K^12N.                                  (28b)
In particular |nA-dA|<=4096K^12N. Since C=dA has N elements and
|C+C|<=KN, (28a) at U=nA proves
 |nA+dA-dA|<=4096K^14N.                               (29)
The proper ratio-set case embeds |B|^2 points in this set; the full
ratio-set case embeds at least half min(|B|^2,p) points in nA+dA,
which lies in2D-2D. Thus in either case
 min(N^2/(8K^2),p/2)<=4096K^14N.
Either K^16>=N/2^15 or K^14>=p/(2^13N). Since K>=1,
 K^16>=2^(-15)min(N,p/N).                             (30)
If |B|<2, the nonempty popular B has size1, so K>=N/2.
Then K^16>=(N/2)^16>=N/2 for N>=2, which also proves(30).
When N<=p^(9/10), p/N>=N^(1/9), so
 max(|A+A|,|AA|)>=(1/2)N^(1+1/144).                    (31)
Here 2^(-15/16)>1/2. No product-energy intersection is lost in6.1.

## 8. Sources

Barak, Impagliazzo, Wigderson, *Extracting Randomness Using Few Independent
Sources*, Sections 3.1--3.4 and Appendix A:
https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/BIW04/biw.pdf
(the AB+C growth step and the additive/multiplicative peeling architecture);
M. Z. Garaev, *A Quantified Version of Bourgain's Sum-Product Estimate in F_p
for Subsets of Incomparable Sizes*, Electronic Journal of Combinatorics 15
(2008), R58, Lemma 3 (the prime-field pivot, of which only the pivot is used).
The finite BSG and minimal-growth tools are proved in the companion module,
with attribution there. The exact-size witness form of (R8) used in Lemma RG
(|V|=M exactly, alpha=1/(2M^e)) and the use of (30) near rate 1/2 follow
mpjunior92's submission f7447a6f-a996-404e-bb7f-df63d9488837 (its
sum-product Section 7.3 and (39)), kept in its revisions
04e91ff4-9abd-474d-810a-1e48d1c4db6c and 7d6289bb-6624-4b01-ae0c-cb7dd033f0e0.
The coarse constants, Lemma RG and Lemma S1' here are independently written;
no new extractor theorem is claimed. All peeling partitions are proof devices;
the evaluator uses only fixed integer arithmetic, AB+C and bit truncation.


Attribution: hitting-probability goodness, the additive-first three-edge graph peeling and the small sampler-unit idea are adapted from published b15a467b-c75f-4152-befe-2b634e34a474. The minus pivot/general additive extension and scalar exponent16, new RG1/292 and the coupled calibration are this submission's new work.
