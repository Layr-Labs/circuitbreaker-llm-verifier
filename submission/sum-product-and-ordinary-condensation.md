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

The selectedfamily usesgamma=256/20481,nu=2^-20,
 eta=gamma/(80+7nu)=268435456/1718070947847,
 theta=eta/7=268435456/12026496634929.
Itscapeta<=1/6000holdsbycrossmultiplication. ThecompleteCZproofuses
unitdeficit5+nu,relativeenergyexponent80+6nu,andthefiniteguard
eta*logp>=64/nu toabsorb2^64intotheremainingnu*eta exponent.
EveryactualNBcallhaseta*r>=7*pi0*w_N>=112*2^20>64*2^20;
r>=7*pi0*2^20/eta>3200000,p>2^(r+a) supplytheindependentringguard.
Theimagecutoffp^(2-5eta)costs30p^-eta;unitdeletionandclippingcost
5p^-eta,stillbelow36p^-eta. RestrictedSections4--8andwrappers3--4
provethefullnewcap,45smoothing,91N7,14U7andallactualalgorithmguards.

# Direct point-line incidence support growth

The scalar sum-product/hitting-peeling detour is not a selected premise.
The actual AB+C four-row map is retained. We prove its support property
directly by counting point-line incidences,using the COMPLETE polynomial
multiplicity/incidence proof inSection6. The evaluator is still only fixed
integer arithmetic,canonical prime selection,AB+C and bit truncation.

## 5. One support-growth step

### 5.1 Exact finite endpoint

Require m>=134217728; this conservative guard implies m>=4. The direct
incidence count below proves the finite endpoint |U|>2^(9m/16-8), keeping
its explicit additive constant rather than replacing it by a fixed exponent loss.

### 5.2 The one-step map and its support lemma

For an integerm>=1,letp_m be the least prime strictly above2^m. The supplied
prime-fields Sections1--2 inrestricted-nonmalleability-and-breaking.md prove
2^m<p_m<2^(m+1) and deterministic trialdivision bitcost
O(2^(3m/2+2)(m+1)^2). Read an m-bitstring as its binaryinteger in[0,2^m),
which injects intoF_(p_m). Define
 E_m(x1,x2,x3)=(x1*x2+x3 modp_m) mod2^m.                (S2)
For a3m-bitstringalpha=(alpha1,alpha2,alpha3),thefourROWS arealpha1,
alpha2,alpha3,E_m(alpha). This is the same actual deterministic row map.

**Lemma S1' (direct incidence support).** For everyintegerm>=134217728 and
Q subset{0,1}^(3m) with|Q|>=2^(3m/2-15),the unionU of itsfourrowimages
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
The complete incidence theorem(I) below therefore applies. Set
 v=floor(2N^2/(d+2))<=2N^(2/3).
As d-1<N^(4/3) andd<=N^(4/3)+1,it gives
 |Q|<=I(P,L)<=2N^2(v+1)+d(d-1)+N^2
 <5N^(8/3)+3N^2+N^(4/3)<=9N^(8/3)<16N^(8/3).
The middle inequality usesN>=1. Taking logarithms with the lower bound
onQ gives log2N>9m/16-57/8>9m/16-8, a contradiction.
All counts concern arbitraryQ and its actual row images; no independence
is asserted. Complete polynomial/multiplicity theoremI is proved below.
At the selected call s=8S, r_F=S+8lg(s+255),
m=s+40r_F>=48S>=48*2^40>134217728.

## 6. Complete polynomial incidence proofs

The complete multiplicity and polynomial-incidence argument below is adapted
from published1d75672c-5080-4151-89d2-cd776a09d06b,actualimmutable009e40109ba442d132016de118d39c2fdf92f206,
and Mark Lewko,arXiv2609.27023v3,Lemma2/Section4(CCBY4.0). The full proof
is supplied;neither source acceptance,the paper nor its Lean formalization
is an external premise. The direct support count in5.2 is this candidate's
application. Write t_+=max(t,0),andp=infinity incharacteristiczero.

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


## 8. Sources and selected scope

The multiplicity and d-form incidence proofs are adapted with attribution
from Mark Lewko,*A Szemeredi-Trotter Theorem in Arbitrary Fields*,
arXiv2609.27023v3(7Oct2026),https://arxiv.org/html/2609.27023v3,
CCBY4.0,and the complete transcription in publishedYukon1d75672c,
commit009e40109ba442d132016de118d39c2fdf92f206. All needed algebra,counting,
polynomial identities and characteristic restrictions are supplied here.
No paper retrieval or formalization/verdict is a mathematical premise.
The direct Q-to-incidence injection,zero-slope handling and finite xi
calibration above are this candidate's changes.

The AB+C support/sampler architecture has its historical hitting/nearhalf
origins in b15a467b and f7447a6f/04e91ff4/7d6289bb through the credited
published precursors. The earlier exact-M witness,scalar10/6,9/4,R7/RG4
peeling and their constants are NOT invoked by the current direct proof.
The failed/unpublishedc1b72 RG4/9/4 route and ownsubmittedfinite-rg-antichain
were material research precursors,credited in the public note;the current
support theorem is the independent direct count just proved. The separate
finite ring/CZ proof remains above,with its own eta_CZ,not a scalarRG exponent.
