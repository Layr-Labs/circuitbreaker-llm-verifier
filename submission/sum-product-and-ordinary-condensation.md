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

The selectedfamily useseta=1213/10000000,gamma=256/20481.
The directlyprovedCZcap iseta<=1/8000. Its strictgamma guard holds:
 gamma-103eta=1124341/204810000000>0.
Theta=eta/8=1213/80000000. The separate explicit guard u>=3200000
in the advice map ensures actual p>2^u meets this theorem's threshold.
The weaker eta*u>=64 guard is not used to infer that prime threshold.
AllremainingfiniteCZ/uniform-sourceconditions retain their suppliedproofs.
The new threshold is proved directly,not inferredfrom a weaker theorem.

# Scalar sum-product and one support-growth step

From the original ordinary-condenser module only the scalar sum-product
estimate (30') of Section 6 is used here, by the new Lemmas RG and S1' of
Section 5; its somewhere condenser, the entropy-growth statement G and G's
proof are not part of this package. All
logarithms are binary, TV is half-L1, and H_inf is min-entropy. Statements
with real entropy thresholds mean literal atom bounds, so integer ceilings are
not silently omitted.

## 5. One support-growth step

This section uses (30') below, in its difference/quotient form, and an
elementary four-edge extraction lemma proved here. U_D is uniform on D.
Fix e=eta_RG=4027/500000, xi=1/745, r=(1-2xi)/(1+2xi)=743/747.
The exact calibration is
 r-(1-3e)=7024507/373500000>0,
 1-124e=163/125000>0, 800000(1-124e)=5216/5>398,
 (e/2-3xi)*300000000=6900/149>17+6e,
 (1-3e)*800000>=7.                                      (RG0)
These are exact rational comparisons. In the lower bound below we keep
min(1-3e,r)=1-3e, avoiding the weaker product r(1-3e).

### 5.0 Four-edge extraction

Let a bipartite graph on sets A,V of sizes m,n have >=alpha*m*n edges,
0<alpha<=1. There is D subset A, |D|>=alpha*m/4, such that every pair
x,y in D has >=alpha^5*m*n^2/2^16 four-edge walks x,v1,z,v2,y.

Proof. Put beta=alpha^2/128. Call an ordered left pair bad if its
common neighborhood has size <beta*n. Sum over right vertices v the
number bad(v) of bad ordered pairs in N(v); the sum is <=beta*m^2*n.
Sample v with probability proportional to deg(v). The probability that
deg(v)<alpha*m/2 is <=1/2. The expectation of bad(v)/deg(v) is
<=beta*m/alpha; hence its probability of exceeding 4beta*m/alpha
is <=1/4, by summing the nonnegative random variable over that event.
Choose v0 satisfying both bounds, and let A0=N(v0). Delete vertices
with >8beta*m/alpha bad partners in A0; at most half are deleted.
The remaining D has size >=alpha*m/4. For x,y in D, at least
 |A0|-16beta*m/alpha>=alpha*m/4
vertices z in A0 are good partners of both. For each z there are at
least (beta*n)^2 choices v1,v2. Their product is the claimed walk count.
Repeated vertices are allowed throughout.

If graph edges carry sums x+v in W in an abelian group, each walk gives
 x-y=s1-s2+s3-s4.
For fixed endpoints, its four labels recover v1,z,v2 uniquely. Distinct
differences use disjoint quadruples. Thus this same D satisfies

 |D-D|<=2^16 |W|^4/(alpha^5*m*n^2).                  (RG4)

The identical proof in a multiplicative abelian group gives |D/D| with
label s1/s2*s3/s4. No passage from difference to sum, or quotient to
product, is needed.

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
P[U_R+U_V in W]>alpha:=M^(-e). Apply (RG4) to its sum graph, with
m=|R|, n=M, |W|<=4M. It gives D subset R with
 |D|>=alpha|R|/4, |D-D|<=2^24 M^2/(alpha^5 |R|).
If D is not times-good, apply (RG4) multiplicatively to its witness
V',W' (all edge products are nonzero). Obtain D' subset D with
 |D'|>=alpha|D|/4, |D'/D'|<=2^24 M^2/(alpha^5 |D|).
Because |R|>M^(1-e), we have |D|>M^(1-2e)/4 and
N'=|D'|>M^(1-3e)/16. Hence the TWO ratios satisfy

 Kminus=|D'-D'|/N'<=2^26 M^2/(alpha^6 |R||D|)
       <=2^28 M^2/(alpha^7 |R|^2)<2^28 M^(9e),
 Kdivide=|D'/D'|/N'<=2^26 M^2/(alpha^6 |D|^2)
       <=2^30 M^2/(alpha^8 |R|^2)<2^30 M^(10e).

Here N'>=2 by (1-3e)logM>=7, and N'<=M. Also p/N'>=p/M>M^r.
Since r>=1-3e>0 by (RG0),
 min(N',p/N')>=M^(1-3e)/16.
The difference/quotient form of (30') gives
 M^(1-3e)/16<=2^22(2^28 M^(9e))^9(2^30 M^(10e))^4
                 =2^394 M^(121e).
Thus [1-124e]logM<=398, contradicting (RG0). D is times-good.

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
>=800000. Also M<=N<p_m^(1/2+xi) (as 2^m<p_m), so p_m>M^(2/(1+2xi)). Let
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
(eta_RG/2-3xi)m>=6900/149>18>17+6eta_RG. Contradiction.

Only set counting and independent uniform variables on U appear; no
independence of the rows of an actual source is used. sampler-and-
disperser.md Sections 3-4 use Lemma S1' at m=P3*B0*L>=374*415*2^40>300000000.

## 6. An elementary finite scalar sum-product theorem

The argument uses Garaev's prime-field pivot (R58, Lemma 3), specialized to
one set. We give its proof and an explicit full-range bound. It does not
invoke Garaev's full unequal-size theorem or an unspecified implied constant.
The sumset tool used here is proved in Section 1.2 of this file:
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

### 6.2 A mixed-sum refinement with explicit constants

For finite nonempty X,U,V in an abelian group there is X' subset X,
|X'|>=|X|/2, such that

 |X'+U+V|<=16|X+U||X+V|/|X|.                         (28)

Here is a proof using only the minimal-growth argument of (3). Put
x=|X|, a=|X+U|, b=|X+V|. In the direct product with the free abelian
group on b+a basis vectors, let T consist of b separately tagged copies
of U and a separately tagged copies of V. For every Y subset X,
| (Y x {0})+T |<=2ab. While the remaining Y has size >x/2, choose a
nonempty Z subset Y minimizing |(Z x {0})+T|/|Z|. Its ratio is <=4ab/x.
The extension proof in (3) gives
 |(Z x {0})+2T|<=(4ab/x)^2 |Z|.
Remove Z and continue until the union X' of removed pieces has size
>=x/2. This terminates and |X'|<=x. The union bound therefore gives
|(X' x {0})+2T|<=16a^2b^2/x. There are ab distinct mixed tags,
each carrying exactly X'+U+V; no unmixed tag coincides with a mixed tag.
Dividing by ab proves (28). Tags are a proof device only.

### 6.3 A full-range scalar estimate for either sign and either product

Let A subset F_p^*, N=|A|>=2, sigma in {1,-1}, and put
 K=|A+sigma*A|/N, L=|AA|/N or L=|A/A|/N.
Then, for EVERY prime p and every such N,

 min(N,p/N)<=2^22 K^9 L^4.                            (30')

In particular this gives exponents (9,4) for Kplus,Ktimes, and also
for Kminus=|A-A|/N,Kdivide=|A/A|/N. Both K,L are >=1.

Multiplicative energy E counts aa'=bb'. Equivalently it counts
 a/b=b'/a', so Cauchy--Schwarz gives E>=N^3/L for either choice of L.
Also E=sum_(a,b in A)|aA intersect bA|. Fix b0 with
 sum_(b in A)|bA intersect b0A|>=N^2/L.
List these N intersection sizes as h1>=...>=hN>=0 and choose an index
j maximizing Q=j^3*h_j^4. Since h_i<=Q^(1/4)*i^(-3/4) and
 sum_(i=1)^N i^(-3/4)<=4N^(1/4),
we have Q>=N^7/(256L^4). The sum inequality follows by summing
 i^(-3/4)<=4(i^(1/4)-(i-1)^(1/4)); factor the difference of fourth
powers to verify each term. Put h=h_j>0 and B={b:|bA intersect b0A|>=h},
b=|B|. Thus

 b^3 h^4>=N^7/(256L^4),  h^4>=N^4/(256L^4).          (30a)

For each z in B, I=zA intersect b0A has size >=h. Ruzsa triangle
with middle -I if sigma=1, or I if sigma=-1, gives
 |zA-b0A|<=K^2 N^2/h.
If b=1, (30a) and h<=N imply N^3<=256L^4, which proves (30').
Hence assume b>=2 and put R=(B-B)/(B-B), excluding zero denominators.
For four coefficients x1,x2,y1,y2 in B set n=x1-x2,d=y1-y2 !=0,
and D=x1A union x2A union y1A union y2A. Then
 |D-b0A|<=4K^2 N^2/h.
Apply (3) at base -b0A and added D, with m=n=2, to obtain

 |nB+dB|<=|2D-2D|<=256K^8 N^5/h^4.                  (30b)

If R is proper, translation by sigma cannot preserve R: repeated
translation by sigma reaches every element of the prime field. Choose
r in R with r+sigma outside R and represent r=n/d by coefficients in B.
Apply (28) with X=dB,U=sigma*dB,V=nB. It gives B' subset B, |B'|>=b/2,
and, using |B+sigma*B|<=KN and (30b),
 |dB'+sigma*dB+nB|<=4096K^9 N^6/(b h^4).
The map (u,v)->d*u+(n+sigma*d)*v is injective on B' x B': a collision
with unequal v would put r+sigma in R, and equal v forces equal u.
Its image is contained in the preceding sumset. Therefore
 b^3 h^4<=16384K^9 N^6.
Combine with (30a) to get N<=2^22 K^9 L^4, proving (30') in this case.

If R=F_p, the averaging proof in 6.1 supplies n,d with
 |nB+dB|>=(1/2)min(b^2,p). Together with (30b),
 min(b^2 h^4,p h^4)<=512K^8 N^5.
Since b<=N, (30a) also gives b^2 h^4>=N^6/(256L^4).
Consequently min(N,p/N)<=2^17 K^8 L^4<=2^22 K^9 L^4.
This proves the full range, without a size restriction relative to sqrt(p).
No incidence theorem, asymptotic constant or dyadic logarithmic loss is used.
The mixed refinement follows the Katz--Shen strategy; the weighted overlap
selection, explicit constants and the two sign/quotient versions are proved here.

## 8. Sources

Barak, Impagliazzo, Wigderson, *Extracting Randomness Using Few Independent
Sources*, Sections 3.1--3.4 and Appendix A:
https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/BIW04/biw.pdf
(the AB+C growth step and the additive/multiplicative peeling architecture);
N. H. Katz and C.-Y. Shen, *A Slight Improvement to Garaev's Sum Product
Estimate*, https://arxiv.org/abs/math/0703614 (mixed-sum refinement).
M. Z. Garaev, *A Quantified Version of Bourgain's Sum-Product Estimate in F_p
for Subsets of Incomparable Sizes*, Electronic Journal of Combinatorics 15
(2008), R58, Lemma 3 (the prime-field pivot, of which only the pivot is used).
The finite BSG and minimal-growth tools are proved in the companion module,
with attribution there. The earlier exact-size witness and near-rate-1/2
peeling architecture follow
mpjunior92's submission f7447a6f-a996-404e-bb7f-df63d9488837 (its
sum-product Section 7.3 and (39)), kept in its revisions
04e91ff4-9abd-474d-810a-1e48d1c4db6c and 7d6289bb-6624-4b01-ae0c-cb7dd033f0e0.
The coarse constants, Lemma RG and Lemma S1' here are independently written;
no new extractor theorem is claimed. All peeling partitions are proof devices;
the evaluator uses only fixed integer arithmetic, AB+C and bit truncation.


Attribution: hitting-probability goodness, the additive-first peeling and
small sampler-unit architecture adapt published b15a467b-c75f-4152-befe-2b634e34a474.
The present scalar proof replaces the old 10/6 estimate by 9/4, and
four-edge extraction replaces the two three-edge/cross-sumset conversions.
The separate four-coordinate CZ proof and its gamma are unchanged; eta is
retuned within its proved cap and strict 103eta<gamma guard.
