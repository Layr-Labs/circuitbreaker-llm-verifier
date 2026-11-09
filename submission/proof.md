# Provenance
Uses the proofs and architecture of 64a9e40a-e9cc-4568-8982-fc9d9bfcb56e. NB uses a compressed seed record. No verdict is a premise.
# A full-B2 circuit-to-affine-flat bridge
All algebra is over F_2. A circuit is a finite acyclic,single-output DAG;
each gate is any one of the sixteen two-input Boolean tables and costs one.
Wires,reuse,depth,sharing and fan-out are unrestricted. Free constants
are allowed,a stronger model for the lower bound.
Attribution: the upstream source named in the supplied proof is GordoAR,
submission 6622943f-7b2d-4f8f-b73c-979c5232484a,version
226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e,proof.md,Chapters II–III;
its public antecedents are 5e9fcf7f-8338-4729-9faa-a18994661af6 and
bf313f35-d97b-4836-b1e8-30afd1157ac2. The supplied proof attributes the
affine-system framework to Mishra–Santhanam,ECCC TR26-260 (2026-10-19).
This reconstruction reproduces that circuit-reduction argument below.
Fix an integer d>=0 and suppose f is nonconstant on every affine flat
of dimension>d. The following reduction supplies the operations used in Section7.
### Representation and potential
During the proof we permit signed references to signals. Signs entering a
gate are absorbed in its truth table. A constant,unary or repeated-parent
gate is replaced by a constant or signed earlier signal. These proof
representations do not assert free native input negation: we start with the
native count and need no conversion of the terminal representation back
to a circuit. The output can be a constant or signed parameter,auxiliary
or gate output.
The six inessential tables are two constants and four unary tables. The
remaining ten are XOR,XNOR and the eight tables
c + (x+b)(y+e),b,e,c in F_2. (2)
Call the latter product gates. These exhaust the basis.
Use M parameter coordinates u,of which 2q are in disjoint pairs. Each pair
ranges over three of its four corners,and every unpaired (free) coordinate
ranges over F_2. Write D for this product domain. Maintain an affine map
Phi:F_2^M -> F_2^n injective on the entire cube,and a representation of
f(Phi(u)) for every u in D. Equality outside D is not required.
Auxiliaries v satisfy a square system
A v = N u + c,A invertible,(3)
with at most three distinct nonconstant variables in each row,after
cancellation. A row has an auxiliary,since A is invertible. Remaining
binary gates read parameters,auxiliaries or earlier gates. Gate outputs
never enter (3). Let s be rows plus gates. A parameter's children are the
distinct rows and gates containing it,not pins and not the output reference.
Set a=M-q. (4)
Initially M=n,q=0,Phi is the identity,the system is empty and s=|C|.
Thus a=n,including circuits omitting inputs.
For each pair,fixing one endpoint opposite its forbidden coordinate
leaves the other endpoint unrestricted. Such safe fixings give a flat
inside D of dimension a; Phi maps it to an a-flat in the original inputs.
Consequently a constant represented output is impossible if a>d. If
a>=d+2,an affine output is also impossible: on a safe a-flat its constant
fiber has dimension at least a-1. A product gate with affine parents cannot
be the output either: on that flat a nonconstant affine factor has a zero
fiber of dimension a-1,and constant factors give a constant output. This covers
constants,parameter/auxiliary literals and terminal gates with affine parents.
### Algebra and cost-free preparations
All coordinate changes below compose Phi with an affine injection. Replacing
a free x by w=x+ell(other parameters),with ell affine,is a bijection of
the whole coordinate cube and of the domain's free fibers; ell may involve
paired coordinates. Fixing a free coordinate lowers a by one. A safe fixing
of a paired endpoint,or restriction of its pair to any two allowed points
parametrized as a line,lowers M,q by one and preserves a. Restricting two
free coordinates to three corners increases q by one and lowers a by one.
Move an affine gate with no product-gate ancestor into (3),processing its
affine predecessors first. Its output becomes a new auxiliary with its
one defining row. The new auxiliary matrix is block triangular with new
diagonal entry one. The row has width at most three; s and parameter use
do not change. Degenerate gates are deleted instead.
The following simplifications decrease s by one,preserve a:
* Delete a constant or literal gate,substituting its output. This includes
constant pins,equal/complementary pins,unary tables,and a gate with
affine parents whose linear parts are dependent. In that last case a
parent is constant or the parents are equal up to a bit,so its result is
a constant or parent literal.
* Delete a gate with no child unless it supplies the output.
* For a row of width at most two,solve for one auxiliary v as a constant
or signed other variable,substitute,and delete that row and column.
Other widths do not increase. Add this row to every other row containing
v,then expand the determinant along v's column: the remaining auxiliary
matrix is invertible.
We need one further elimination. An auxiliary Q with coefficient one on a
free x has the form Q=x+ell(other parameters). Exchange their roles: Q is
free and x auxiliary,without expanding any equation. If b=A^(-1)N_x,
then b_Q=1. Replacing A's Q-column with N_x gives determinant
det(A)b_Q=1,preserving invertibility and row widths. Update Phi by
x=Q+ell. The exchange preserves s,M,q.
If free x has exactly one child,a row E,N_x is nonzero. Some Q thus
has b_Q=1. Exchange x,Q. Now auxiliary x occurs only in E and,during the
loop a>=d+2,cannot be the output by the output exclusion. Delete
E and x,which E determines but nothing else uses. Expanding along
x's one-entry column proves invertibility of the remaining system. This deletion loses one row charge and preserves a.
### Clean all pairs
Perform available affine transfers and the preceding eliminations. Then
apply the first applicable pair rule below,propagate constant/literal
substitutions and short rows,and restart. A fixed parameter deletes each
of its old gate children,since a binary table with one fixed pin is unary,
and each old row child,since that row shortens to width at most two.
Count each removed object once. All four rules preserve a.
1. If endpoint p is unused,safely fix its partner h. The new free p is
unused and q falls by one.
2. If an endpoint has at least two children,safely fix it. At least two
distinct objects disappear;its partner becomes free and q falls by one.
3. Otherwise each endpoint has exactly one child. If partners share that
child B,its gate output or sole auxiliary is a function of the pair
(a row B has just the endpoints and that auxiliary). Two of the three
allowed points give the same bit. Restrict the pair to their affine line.
Delete B after substituting that constant bit. No other object used the
pair,so the new line coordinate is unused.
4. Otherwise suppose p's sole child is a row E. Its partner h has a distinct
sole child B. Safely fix h and delete B. If propagation deletes E too,
two distinct charges have been consumed. If not,now-free p still occurs
only in E: its only original occurrence cannot spread without E's
deletion. Apply the one-child-row elimination,consuming E's distinct
charge. Thus at least two objects disappear and q falls by one.
When these fail,every paired endpoint has exactly one gate child; partners
have different children,and paired parameters do not occur in (3).
In particular every auxiliary depends only on free parameters.
### Prepare a product source
Continue only while a>=d+2. If all gates transfer away,the output is affine,
a contradiction. Choose a product gate G with no product ancestor; its
predecessors are affine and have been transferred. After the simplifications
its parents have nonzero independent linear parts. Turn any auxiliary
parent Q into a free parameter by a role exchange. Such a free pivot exists
because Q is nonconstant and all paired columns of N vanish. If the other
parent is free,independence supplies a pivot different from it; if paired,
Q is independent of that coordinate. With two auxiliary parents,exchange
the first and then use this same argument for the second.
After each exchange,check for simplifications and pair cleanup. If any
is possible,perform a completed action and restart. Otherwise keep this
same G: at most two exchanges make its parents distinct parameters. G is
not the output,by the product-output exclusion,and has a gate child by
dead-gate cleanup. It cannot have a row child.
### Restrictions and the sharing cases
Each restriction lowers a by one. The new safe dimension
is at least d+1,so the restricted output is nonconstant. Freeze old gates
and incidences immediately before the restriction. Propagate topologically,
recording old outputs as constants or signed aliases. Every old gate hit
by a constant pin is removed. If an old output becomes constant,all its
old children are removed too,even if it first became an alias. Charges
refer to distinct removed old objects,not to propagation events.
**Two free parents x,y.** Restrict to the three corners on which G has its
majority value,creating a pair. G and its k>=1 distinct children disappear,
even with shared descendants or children that also read x or y. Thus
Delta q=1,Delta M=0,Delta s<=-(1+k),Delta a=-1.
**Paired p and free x.** Safely fix p's partner h and fix x to kill one
factor of G. Its constant output deletes G. Newly free p loses its sole
use and is absent from the output reference,so is unused. The net
changes are Delta M=-2,Delta q=-1,Delta a=-1.
**Paired p and paired r.** Their pairs differ by pair rule3. Safely fix
their partners h,k,then fix now-free p to kill G's product. The newly
free r is unused,since its sole child G is constant and paired coordinates
were absent from all affine rows. Here Delta M=-3,Delta q=-2,Delta a=-1.
These recipes exhaust the parent types;all substituted coordinates are
restricted through the affine injections specified above.
### Termination and conclusion
Every completed action preserves representation on the current D and
strictly decreases the nonnegative integer M+s: eliminations decrease s,
pair rules decrease M,and each parent-type action decreases M+s. Neither M nor s increases.
Preparatory transfers decrease the number of gates; there are finitely
many. A fixed-source preparation uses at most two role exchanges before
an action. Thus preparation cannot cycle. The geometric exclusions prevent
termination inside a>=d+2,so the process reaches a_T<=d+1.
# Exact finite direct-inner and affine-outer assembly
The dependency order is noncircular:
1. PXT/PXG7g and PXTM7h and exact fields feed restricted Sections1--3 and8.2.
2. Restricted8.4's shortened four-round8-copy corollary feeds wrappers4.2's F#.
3. CZ/E6,canonical LinExt/W2,Ext0,Ext2,F# feed wrappers4's NB.
4. NB with B9's complete seed/threshold feeds outer3's A# at the sampler degree.
Constants use exact integer/rational operations or proved terminating searches.
For positive rational x,lg(x) is the least nonnegative integer j with2^j>=x.
## 1. Inputs and their precise scope
All inputs are included. Linear Sections2--4,7a,7g--7h prove the arbitrary-source
PXT/PXG/PXTM vectors and exact fields; the same file supplies bilinear/weak-seed/AGHP.
Only BASE uses affine Section7f.
Sampler Part SP proves finite rings and complete incidence: gamma=((2-nu)/(14+(h_R+5/2)/D_E)),
prime guard2^(2^90). Restricted supplies prime/encoding,CZ smoothing,
conditional tools and the four-round antichain8.2 theorem and its five-label shortened corollary8.4 retaining all9
seeds. Wrappers4.2 instantiates the zero-shift8-copy F# in NB4.6 step7.
Wrappers2,4 prove the all-seeds affine transformation and NB. Outer3 builds
A# at t287,B9 pays its full NB call. Sampler3--6 prove W_A/G_A/AS2 and the
general conditional bad-seed lemmas. Sampler0 resolves legacy
component names to these included bodies. Every component seed is retained.
## 2. Primitive convention
Selected PXG calls use PT2 in7g; their heights
are at most13. Establish threshold guards by PT1 before applying PT2.
PXT source coefficients Lam,g remain unchanged.
Selected slots use the proved PXT/PXG7g and PXTM7h primitives;
explicitly named AGHP calls use the separately proved matrix extractor. Each has
its explicit source,seed and dimension guards.
## 3. The affine assemblies and the standard-breaker interface
### 3.1 The standard interface required here
A standard one-tampering breaker C is used in the outer assembly through
its actual seed lengthd_C and thresholdk_C: arbitrary correlated source
pair independent of an arbitrary correlated seed pair,honestseed uniform,
distinct advice,bothcomplete seeds retained,error<=2^-e. HereC=NB has
explicitguardsk_C<=u,d2<=u,(PG),(AG1) s>=a+40,(AG2) eta*(r+a)>=(905)(a+sig+2) of wrappers4.3
and works at every advicelengtha>=1 passing them; failed guards give 0^z.
NB's innerF# instead has the four-round antichain8-copy tuple contract of wrappers4.2
and restricted8.4 at five fixed labels; NB4.6 step7 uses exactly that
zero-shift contract. The outer assembly retains arbitrary additive shifts.
### 3.2 The selected outer affine assembly
The outerA# is wrappers2.1 withC=NB,the affine-left guard of2.5,
and the explicit outerSection3 slots.
F# is the map of wrappers4.2. Failed public guards return the all-zero
string; every guard holds at every sufficiently large actual length.
## 4. Fix the uniform-source advice condenser first
Fixgamma=((2-nu)/(14+(h_R+5/2)/D_E)),eta=((2-nu)/(14D_E+h_R+5/2)),theta=eta/4=((2-nu)/(4*(14D_E+h_R+5/2))).
Every actual NB call haseta*(r+a)>=64*2^80,as proved
in wrappers(I3),so this guard and the independent ring-prime guard hold.
NB's repaired honest Ext0 output is exactly uniform,so it uses the four-block Lemma U
of affine-and-standard-wrappers.md Section3.5.
Put r0=r+a (the pre-subtraction length from wrappers4.3),sig=Ee+4 and s=floor((eta*(r+a)-11sig-8)/4). NB's finite parameter
algorithm proves s>=m0+2e_L; with its explicit guards (AG1) s>=a+40 and
(AG2) eta*(r+a)>=(905)(a+sig+2) it proves all Lemma U and E6 guards (wrappers(I1),(I2)).
Its advice row error is at most14*2^-sig<2^-Ee.
## 5. The standard breaker NB and Lemma B9
Theorem NB (affine-and-standard-wrappers.md §4) is a standard breaker
in the sense of §3.1,with guards k_C<=n,d2<=n,(PG),(AG1) and (AG2),at every a,z,E>=1.
It uses F#,the Section4 advice condenser,canonical LinExt/W2,Ext0 and Ext2; the outer
assembly uses NB only through A#,and Lemma B9 gives its seed and threshold.
**Lemma B9 (NB at the ACTUAL outer call).** Source n'=d_O,output
z=(43389108270635/1099511627776)Om,advice a=s_tag+2+lg C,exponent e,M=2^40,Om=Mceil((e+L+54)/M). Put
eps_N=10^-18,e_B=1728/7-6*2^-100,alpha=e_B/(e_B+Q_S)+eps_N.
The exact slots below have d_O=floor((Dslot)Om)<2^30Om,and O2a gives
Om<=c_O L<1024S: c_O*Q_S*(1+2^-80)<1024 and S>=2^80 (5'').
Hence 4n'<2^32Om<2^42S and lg(4n')<=42+lg S. With Ee=e+4,the actual NB unit satisfies
w_N=max(2^80,Ee+44+2lg(4n'))<=alpha*Om.
Indeed Om>=E+L+67 and E<=e_B*S+47,L>=Q_S*S imply
alpha*Om-E >= eps_N*(e_B+Q_S)S+67alpha-47(1-alpha).
The constant is positive. The nonfloor part of w_N is at most
E+146+2lg S. The dyadic certificate
5770*80+29028=490628<eps_N*(e_B+Q_S)*2^80
dominates 146+2lg S at every S>=2^80,hence under the new guard.
Since s>=7S and Om=Mceil((6s+204S+L+89)/M),Om>=(246+Q_S)S. Exact arithmetic gives
alpha*(246+Q_S)>1. Thus alpha*Om>=S
pays the floor2^80.
The inner unit w_F=max(2^80,Ee+2lg(4n')+32,ceil(z/128))<=alpha*Om,
because its nonfloor terms are bounded by w_N,while
2^80<=S<=alpha*Om pays the new fixed floor,and
ceil(z/128)<=z/128+1<=alpha*Om follows from 50/128+eps_N<=alpha.
All seven inner leaf-log budgets are<=8; NB
leaf-log budgets14,5 obey Ee+depth+6<=Ee+44. All generic component proofs apply.
The new NB guards are discharged without using Ext0's seed bound first.
Since Om<=c_O*L and L<Q_S*(S+1),the exact scalar check
c_O*Q_S*(1+2^-80)<440 proves Om<=440S. Also w_N>=Ee>=246S.
Thus DX=43389108270635/1099511627776 and DX*440<=128*246 imply ceil(z/128)<=w_N;
the other inner terms are <=w_N,so w_F<=w_N. First
m1=d_F<=(d_inner+2^-80)*w_N (paying its final ceiling),
and d_inner+2^-80<=256*BY proves the Ext2 cap. Its clean PTM.2 bound gives
m0<=(d_inner/BY+pi2+2^-80/BY)w_N. Put
den=theta-bc0,Runit=(d_inner/BY+pi2+p0+6)/den.
The positive extra (1/4-2^-80/BY)w_N pays56+den:
check (1/4-2^-80/BY)2^80>=56+den,then monotonicity. Thus the literal r0
ceiling obeys r0<=Runit*w_N. Since r=r0-a<=r0,the bound v<=b*r+1 and the exact check
b*Runit+2^-80<=256*BM proves Ext0's cap BEFORE invoking its PTM.2.
Terminal rates,kap,arities and heights obey7g for both calls.
The F# guards hold: h>=100w_F>=100(w_N-12)>=24600S-1200>z.
For the final call,Ee>=e1*S+39 with e1=1728/7-6*2^-100,
so Ee>=(1728/7-2/10^6)S. The exact tag formula,L<Q_S*(S+1),Q_S<3 give
Om<=(1728/7+Q_S)S+M+116<=(1728/7+Q_S+1/10^6)S,
since(M+116)*10^6<2^160. Exact multiplication gives
(43389108270635/1099511627776)*(1728/7+Q_S+1/10^6)<(40)*(1728/7-2/10^6).
Hence z=(43389108270635/1099511627776)Om<=40w_F.
These bounds precede all seed envelopes.
For each integer slot ceil(A*w_F),use ceil(A*w_F)<=A*alpha*Om+1.
Thus m1<=D_F*Om and k_F<=K_F*Om,with
D_F=alpha*d_inner+eps_N,K_F=alpha*k_inner+eps_N.
Let (c_P,pi_P,Lam_P,g_P) be the exact PTM.1/PTM.2 for both arity vectors coefficients at vNB0,vNB2.
Put b=6-eta/2. Ext0 outputs v=6r-floor(eta*r/2)<=b*r+1,so
the finite PTM.2 bound is d0<=bc0*r+p0*w_N+D0+c0,
p0=200846087453311/1000000000000,D0=0. Define
M0=c2*D_F+pi2*alpha+eps_N,
K2bar=Lam2*D_F+g2*alpha+eps_N,
Amin=(12342857/100000)/(1728/7+Q_S+1/10^6),
Rnew=(M0+p0*alpha+(23/4)*alpha+(56+D0)eps_N)/(theta-bc0)+eps_N-theta*Amin/(theta-bc0),
d0bar=bc0*Rnew+p0*alpha+(D0+2)eps_N,
d1bar=M0+2d0bar+6alpha+eps_N,
DNB=2d0bar+2d1bar+D_F+K2bar+alpha.
All denominators are positive and p0>=16. Exact substitution gives
ceil(M0)=67, ceil(K2bar)=1735,
 ceil(Rnew)=164067, ceil(d0bar)=265,
 ceil(d1bar)=602, ceil(DNB)=4258.
Ee<w_N<=alpha*Om and 56+D0<=(56+D0)eps_N*Om give
(23/4)Ee+56+D0<=(23/4)*alpha*Om+(56+D0)eps_N*Om; the r0 ceiling includes -bc0*a in its numerator. The raw length is r=r0-a.
Since s>=beta*S,a=s_tag+9=3(s+34S)+9>=(3beta+102)S>=(12342857/100000)S.
The last exact comparison uses beta=50/7-2^-100. The Om upper bound
proved above gives a>=Amin*Om. The two discounts sum to (1+bc0/(theta-bc0))*a=theta*a/(theta-bc0),
paying the displayed theta*Amin/(theta-bc0) before invoking Ext0. Its output cap still follows from r<=r0.
Advice guards. By (O2),e=2s_tag+35,so2a=2s_tag+18=Ee-21 EXACTLY;
hence a<Ee/2. With sig=Ee+4,w_N>=Ee+44 and p0>=16:
(AG1) wrappers(I1) gives s>=m0+bc0*r+p0*w_N+D0+c0+3Ee+41>=3Ee+40>6a+40>=a+40,
termwise for every S; at S=2^80,s/(a+40)>500.
(AG2) Ee>=S>=2^160 gives w_N=Ee+44+2lg(4n'). Ext2 has B=30 and rho_2>=w_N-33,so its two slices total at least (1-2^-40)*sum_i(R_i+1)*w_N. Put a_px=Ee+lg(B)+4. Then w_N-a_px<=119+2lg S<=S/1099511627776<=w_N/1099511627776: on2^j<=S<2^(j+1),j>=160,the left bound is121+2j,starting441<2^160/1099511627776 and growing by2 while the right doubles. Also w_F>=w_N-12>=(1-2^-40)w_N,hence m1>=(1-2^-40)*d_inner*w_N.
The new leaf seed is ell+4a_px+lg(4(m_h+1)),with ell>=m1/B.
Including it gives m0>=SD*w_N,
SD=(1-2^-40)*(sum_i(R_i+1)+d_inner/B+4).
Here33<=2^80/1099511627776 pays both slice and w_F rounding.
For this actual untrimmed caller r_F=S. Set W=1728/7+10^-6,
G=2715/2 and Gplus=G+10^-6.
The exact tag formula gives E<=(1728/7-6*2^-100)S+47.
With Ee=E+17 and lg(4n')<=42+lg S,
w_N=max(2^80,Ee+44+2lg(4n'))<=(1728/7)S+256+2lg S<=W*S.
For2^j<=S<2^(j+1),j>=80,the last overhead is at most258+2j;
its base418<2^80/10^6 and increment2 are dominated by doubling.
The maximum's floor is also at most W*S. Since Ee>=246S,
the floor is inactive and w_N-Ee<=128+2lg S<=S/1099511627776<=w_N/1099511627776.
Indeed its dyadic bound130+2j starts290<2^80/1099511627776 and grows by2.
Thus Ee>=(1-2^-40)w_N,without using any extractor guarantee.
From2a=Ee-21 and sig=Ee+4 we have the exact identity
2(a+sig+2)=3Ee-9. As Ee<=w_N and S>=w_N/W,
(905)(a+sig+2)<=G*w_N-(8145/2).
The r0 formula subtracts bc0*a but retains its positive(23/4)Ee term.
Since a<Ee/2<=w_N/2,
eta*(r+a)=eta*r0>=eta*(p0+SD+(23/4)*(1-2^-40)-bc0/2)*w_N/(theta-bc0)>=Gplus*w_N.
The last rational inequality is checked exactly. Hence(AG2) holds
for every legal S;10^-6*2^80>6*(905) also pays the old finite constant.
At S=2^80 the actual ratio eta*(r+a)/((905)(a+sig+2))>1.0000001.
The exact inequality (15393/905+5*2^-80)*eta<=gamma pays the new energy budget. Both guards hold at the only NB call (outer A#),so Theorem NB applies there.
The r0 ceiling costs <=eps_N*Om and a>=Amin*Om,so r_NB=r0-a<=Rnew*Om. Both source branch coefficients are at most KC-eps_N,
KC=1233515101920841359870607/500000000000000000. These exact inequalities are recomputed from the displayed recurrences.
For LinExt put J=lg(4(Rrow+1)). Since Rrow>=3,m0>=1,
4(Rrow+1)<=(Rrow*m0+1)^2,hence J<=2lg(Rrow*m0+1).
Now 2e_L<=d0+3Ee and ceil(2*(Rnew+2)*M0)+1<2^240 give
2lg(Rrow*m0+1)<=480+4lg Om<=eps_N*Om. The k80 base is
804<eps_N*2^80,with slope4<=804; doubling covers every Om>=2^80.
Thus d1<=m0+2d0+6Ee+eps_N*Om<=d1bar*Om and d_C<=DNB*Om.
Define U18(x)=ceil(10^18*x)/10^18 and Slot(x)=U18(x)+2eps_N.
Put SC=Slot(DNB)=4257477769389152692011/1000000000000000000,KC=Slot(max(K00bar,Kotherbar))=1233515101920841359870607/500000000000000000,
where K00bar=b*Lam0*Rnew+g0*alpha+eps_N and
Kotherbar=K_F+10M0+43389108270635/1099511627776+2alpha+eps_N. The extraeps_N pays
the three-bit selected-row label in the compressed10m0-bit seed record and the output rounding: exact checks give
eps_N*Om>=max(4,Lam0),sinceOm>=2^80.
Since eps_N*Om>=1,floor(SC*Om)>=d_C and floor(KC*Om)>=k_C.
These are the actual padded seed and entropy slots used below.
The outer d3c contract ensures k_Cslot<=n'=d_O and d2<=d_C<=s_C<=n'.
The NB algorithm still has eta*(r_NB+a)>=4p0*w_N>=64*2^80 and (AG1),(AG2) above,
and r_NB+a=r0>=4p0*w_N/eta supplies the separate ring-prime guard u>=2^90.
All base/hash degrees remain O(log actual n); large relative fields use towers.
## 6. Uniformity,exactness,and remaining obligations
For fixed t and m,a,E=O(log n),slot and NB/inner dimensions are O(log n).
Closed slots and LemmaB9 give F#,Ext0 and Ext2 logarithmic enumerated
field degrees,so setup costs2^O(logactualn). CZ/AGHP/bilinear setup is
2^O(local dimension),making the composition polynomial in actual n.
Its exponent may depend on the fixed t=4C-1; no joint bound in t,n is claimed.
For ordinary seeded primitives or standard breakers with the complete
source tuple independent of the complete seed tuple,unused suffixes may be
retained by their conditional kernel given all used prefixes. The outputs
depend only on the source tuple and those prefixes. This applies to the
primitive pads of F# and the outer A# (restricted-nonmalleability-and-breaking.md
§§8.2--8.4,affine-and-standard-wrappers.md §4.2,outer §3) and to the
zero-additive-shift independent-side use of NB's inner breaker
(affine-and-standard-wrappers.md §4.6,step 7). It does not assert a generic
kernel principle when arbitrary additive B-tuples correlate with suffixes.
# Tagged affine-permutation logarithmic-loss family: outer assembly
The affine-breaker contract (F') is proved below from the outer A# (§3) and Lemma
B9.
Credits: sumset framework,Chattopadhyay--Liao,
https://arxiv.org/abs/2110.12652v1 ; advice-condenser architecture,Xin Li,
https://eccc.weizmann.ac.il/report/2023/023/revision/1/ ; native full-B2 bridge,
minimal-circuit-bridge.md.
All logarithms are binary. Write lg(x)=ceil(log2 x) for positive x,computed
by integer bit length and a power-of-two comparison.
Foractualintegern>=1,flg(n)=floor(log2n),computedasbitlength(n)-1.
ONLYtheactualouterlengthusesflg.Everylgforfielddegrees,seeds,graphconstants
remainsCEILINGlog2.
## 1. The exact inputs to this assembly
Fix Q_S=2139/1000. With L=floor(log2(actual input length)),S=floor(L/Q_S)
is an integer; S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer rounding.
(L_A) Put L=floor(log2 n),S=floor(L/Q_S)>=2^160,s=2ceil(beta*S/2),
ellc=lg(s+4) and r_F=S. Linear Section7f supplies the literal
three-condenser/kernel base E_n,source-linear at each seed and strong
for every affine source of dimension at least
K_B=ceil(A_B*d+B_B*rho),rho=r_F+225+2lg S,
d=floor((Dslot)Om),with Om defined below and d<=753316552L.
Its complete seed length D<=s,error<2^-r_F and public dimensions are
AF1. The full AF2--AF4 proof gives K_B<=217527117634L<=HL.
(G_A) Sampler Sections3--5 use q_walk=17,C=72,m_walk=s+34r_F,
s_tag=3m_walk<=144S+960ellc,A=A_rows=2^s_tag and
Lambda_S=A^(1/2)/2^17. The full four-row AB+C incidence proof followed
by W_A expands every Lambda_S-row set to more than5*2^-r_F*2^s
neighbors. Direct affine averaging gives AS2: for every uniform affine
source of dimension>=K_B and fixed test of density<=2^-r_F,fewer than
Lambda_S rows are wholly bad except on probability<=2/5.
The evaluated sampler is source-linear and costs2^O(s_tag).
(F') For every n passing the guards (5'') of outer §4,the outer A#
(outer §3) with C=NB,at the request (O2) and the unit Om,returns an
exact seed length d_F,source entropy requirement k_F and a one-bit
generalized affine breaker,where
d_F=d_O<=delta_O*L,k_F=k_O<=kappa_O*L. (F')
For arbitrary correlated tuples (A0,...,At) independent of
(B0,...,Bt,Y0,...,Yt),with H_inf(A0)>=k_F+e,Y0 uniform on d_F bits,and
alpha0 different from every alphaj,its honest output on Xj=Aj+Bj is
2^-E-close to uniform jointly with every tampered output and EVERY seed
Y0,...,Yt;here e is the component exponent in(O2). For the actual affine-left tuple A_i=T_i Z+c_i with uniform affine Z
and invertible T_i,threshold k_F alone suffices by wrappers2.5;(F_aff').
The field is F_2,so addition is XOR. Neither source copies nor
seed copies are required to be deterministic functions of their honest copy.
For the fixed t and a,E=O(log n),its complete evaluator is polynomial in
the actual n (the evaluator of affine-and-standard-wrappers.md §2.1
with NB,whose cost is stated in §4.7 there).
## D1. Algorithm and finite hypotheses
Write L=flg(n),S=floor(L/Q_S),s=2ceil(beta*S/2),ellc=lg(s+4),r_F=S,m_walk=s+34r_F,s_tag=3m_walk=log2 A,w=2. Require S>=2^160,L>=2048 and s_tag>=2L; the latter follows from the fixed constants below.
Build the canonical field F_(2^L),and split x into L-bit coefficients of
P_x,padding its final block with zeros. For a full s_tag-bit row alpha,
let beta,gamma be its first and second L-bit blocks. Define
H_alpha(x)=Prefix_2(gamma*P_x(beta)),a two-bit linear map in x.
Each (beta,gamma) has exactly A/2^(2L) row preimages. For Delta!=0,
P_Delta is nonzero of degree ceil(n/L)-1<n/L,so its root fraction is
at most n/(L*2^L)<2/L<=1/1024. Heren<2^(L+1),L>=2048;
thedegreeisalso<2^L,sotherootboundapplies. Off those roots,multiplication by P_Delta(beta)
permutes gamma and exactly 1/4 of gamma values have zero two-bit
prefix. Thus at most A*(2/L+1/4)<=257A/1024 rows have H_alpha(Delta)=0.
The root bound follows by factoring X-beta at distinct roots. The supplied
finite-field proof makes the canonical search halt in 2^O(L) time.
Full alpha remains in the advice label despite repeated (beta,gamma).
No independence of hash values across different rows is asserted.
Use the source-linear AFFINE sampler Samp(x,alpha,j),j in [C],seed
length d,with AS2 at dimension k_S=K_B: for every uniform affine source
and fixed test of density<=epsilon0=2^-r_F,fewer than Lambda_S rows
are wholly bad except on probability<=2/5. Only uniform affine sampler
sources are used below. The enumerated field index alpha replaces [A].
For t=4C-1 and a=s_tag+w+lg(C),request the padded all-seeds-retained affine
breaker F of output length one,affine-left threshold K and ordinary
threshold K_gen=K+e,and error at most
epsilon = 1/(2^22 A^2).
The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1,and w-bit H_alpha(x). Define
R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),advice(x,alpha,j)),
D(x)=1 iff 2 sum_alpha R_alpha(x)>A.
For the affine application below,assume the source entropy satisfies
k >= max(k_S,K+4w+4Cd). (D1)
All concatenations,zero padding and majority ties are specified. There
are no seed or advice inputs to D. The tags are internal computations on
the actual n-bit x.
## D2. Mean error with arbitrary correlated left copies
Fix an arbitrary correlated tuple(A_0,...,A_r),fixed advice differing
from honest advice,and deterministic seed maps L_i of A_i. Assume
Hinf(A_0)>=K_gen+ell+lg(1/epsilon),with joint seed leakage at most2^ell values.
For each seed y let f(y) be the MAXIMUM TV,over all finite choices of
b_0,...,b_r,y_1,...,y_r,of the honest bit in
F(A_0+b_0,y+L_0(A_0),advice_0),
(F(A_i+b_i,y_i+L_i(A_i),advice_i))_(i=1)^r
from a fresh uniform bit jointly with the other outputs. Then
E_(uniform y) f(y)<=2epsilon. (D2mean)
Indeed choose a maximizing witness tuple for EVERY y; finite maxima exist.
Take Y uniform independent of the left tuple and let all b_i,y_i depend on
Y through those witnesses. Conditional on the complete L_i leakage,left
and right tuples are independent and honest seed Y+L_0(A_0) is uniform.
Bounded-range conditional entropy loses at most ell average bits. Its tail
at probability epsilon and the supplied conditional lift of(F') give TV
at most2epsilon,retaining all seeds and leakage. The honest seed and
leakage determine Y. Drop other retained data,then leakage: contraction
bounds TV retaining Y and tampered outputs. Since Y is independent of the
original left tuple,that TV is exactly the average of f(y). This proves
(D2mean). Markov also gives density<=2epsilon/theta for{f>theta}.
Witnesses and this function are finite proof objects,never evaluator searches.
Affine specialization: let A_i=T_i Z+c_i for one uniform affine Z and
invertible linear T_i,and let every L_i be linear. Conditioning the joint
seed leakage costs at most ell affine codimensions. Every nonempty fiber
is uniform affine,and every honest A_i has entropy>=Hinf(Z)-ell.
Thus Hinf(Z)>=K+ell lets the same maximizing-witness construction apply
(F_aff') on EVERY fiber without an entropy bad set. The fresh W remains
uniform independent of the left tuple,and W+L_0(A_0) is uniform on each
fiber. Averaging retaining W,and then projecting,gives E f(W)<=epsilon,
hence the displayed2epsilon bound without an entropy-tail logarithm.
This remains true after fixed linear-tag conditioning.

## D3. Fixed rank-good neighbors and actual parities
Fix independent uniform affine sources X_1,X_2 of dimension at least k.
Put A_0=X_1,A_1=Lambda X_1,B_0=X_2,B_1=Lambda X_2+b,with Lambda
invertible and tau(x)=Lambda x+b fixed-point-free on the sum support.
For each sigma define Bad_sigma as base seeds gamma for which the linear
map E_n(.,gamma) has rank<d on the direction of B_sigma. Its output on
that affine source is uniform on an affine subspace. If rank deficient,
its distance from U_d is at least1/2. The strong affine base error epsilon0
therefore gives |Bad_sigma|<=2epsilon0*2^s.
The graph G_A expands every Lambda_S-row set to more than5epsilon0*2^s
neighbors. Consequently fewer than Lambda_S rows have every neighbor in
Bad_sigma. Intersect the two good row sets to get a FIXED G,omitting
q<=2(Lambda_S-1)<2^-16*A^(1/2) rows. Thus g=|G|>=255A/256.
For every row-copy sigma,alpha in {0,1}xG choose a fixed neighbor j whose
base map has full rank on B_sigma. This choice depends on the source law,
not its sampled point. Its offset Samp(B_sigma,alpha,j) is exactly uniform
on d bits. All expectations below are over the joint X_1,X_2 law.
For a two- or four-set Q of distinct row-copy indices,expand its at most
4C calls. Choose an honest row-copy as in the singleton/fully-doubled
rule below and use its fixed rank-good neighbor;degree t=4C-1 suffices.
Condition the selected left tags,at most4w linear bits,and then the at
most4Cd complete left seed-leakage bits. Every nonempty fiber is uniform
affine and each honest A_sigma has entropy>=k-4w-4Cd>=K. The right law
X_2 is unconditioned and independent of this left fiber. Source-linearity
splits every seed into its fixed left leakage and right offset;the honest
seed remains uniform by the fixed full-rank choice.
Now enumerate the at most2^(4w) selected RIGHT tag profiles. For each fixed
profile with honest advice distinct from all tampered advice,apply(F_aff')
and its WHOLE-right retention to this conditional left law and R=X_2.
This bounds the absolute parity expectation restricted to the right-profile
event by2epsilon: that indicator is a function of retained R,and the ideal
honest bit is independent uniform. We do NOT condition the right law on
that event or assert uniform seed after such conditioning. Sum these
restricted expectations,then average the left fibers. Thus each legal
parity has bias<=b_Q=2^(4w+1)epsilon. If an alpha occurs only once in Q,
choose such an index as honest;all profiles are legal because alpha,j
identify that call. Otherwise Q is fully doubled;choose its least index.
Only profiles with a selected collision can then be illegal. Charge their
probability separately,giving bias<=b_Q+Pr[at least one selected collision].
The tag root count bounds collisions at EVERY(X_1,X_2) point by
N_col=257A/1024. No independence of collision events is asserted.
Use a_Q=1 for four-sets;for two-sets use A^2 for matched copies of one row
and A otherwise. The same counting gives sum_Q a_Q<=(2/3+3/256)A^4.
Hence Bsum=sum_(Q on G) a_Q*b_Q
<=(2/3+3/256)*2^9*epsilon*A^4
<=A^2/8192<3gA/35,by epsilon=1/(2^22*A^2),w=2.
The graph,rank argument and error budget require no pointwise bad-seed
set and no exceptional right-source probability. The two actual row-copy
arrays and every collision term are retained for the quartic identities.
## D4. A robust quartic sign certificate
Let X_alpha,Y_alpha be actual signs on G and S_0=sum X_alpha,S_1=sum Y_alpha.
Put Z_alpha=X_alpha Y_alpha,Z=sum Z_alpha,zeta=E Z. Only matched pairs occur in zeta,each with weight A^2. §D3 gives
|zeta|<=N_col+Bsum/A^2<N_col+3g/(35A).
Since g>=255A/256 and A>=2^256,
r=|zeta|/g<=257/1020+3/(35A)<13/50.
For Q_+(x,y)=xy(x-y)^2 and Q_-(x,y)=-xy(x+y)^2,expand the ACTUAL
ordered quartic monomials. All-even terms contribute-2g^2. A matched-pair
parity has coefficient plus or minus2(3g-2),and fully doubled two-row
parities contribute-8 sum_(alpha<beta) E Z_alpha Z_beta.
Every remaining parity has a singleton row and is bounded by b_Q.
A size-four parity has coefficient6,6,or-8 according to its two-copy
counts(3,1),(1,3),(2,2). For size two,a mixed nonmatched pair has absolute
coefficient6g-4; a same-copy pair has coefficient-4g. For a nonmatched pair its coefficient is at most6A and a_Q=A; for a four-set it is at most8 and a_Q=1. Thus every remaining absolute coefficient is at most8a_Q.
Thus the exact identities give
E Q_+(S)=2(3g-2)zeta-2g^2+4g-4E Z^2+u_+,
E Q_-(S)=-2(3g-2)zeta-2g^2+4g-4E Z^2+u_-,|u_+|,|u_-|<=8Bsum.
All doubled collision terms are retained exactly,not charged separately.
The variance identity E Z^2>=zeta^2 now bounds BOTH expectations by
2(3g-2)|zeta|-2g^2+4g-4zeta^2+8Bsum
<=[6r-4r^2-2+4/g+24A/(35g)]g^2
<=[-444/625+1/1000+6144/8925]g^2<-g^2/50. (D4')
Here6r-4r^2 increases up to13/50,A/g<=256/255,and
4/g<=2^(3-2L)<1/1000 forL>=2048.
Marginal fourth moments have all-even part3g^2-2g. Remaining parity
coefficients are24 for size four and12g-16 for size two,each bounded by
24a_Q. Thus
E S_sigma^4<=3g^2+24Bsum
<[3+72*256/(35*255)]g^2<16g^2.
Their L4 norms are consequently at most2sqrt(g). This weaker norm bound
is sufficient for the deterministic removed-row perturbation below.
Let Y_sigma be the ACTUAL sum over ALL rows in that copy. Pointwise
|Y_sigma-S_sigma|<=delta=q<2^(-16)A^(1/2). Since g>=A/2,
delta/sqrt(g)<2^(-16)A^(1/2)/(A/2)^(1/2)=2^(-31/2)<1/40000.
The L4 norms of S_sigma are at most2sqrt(g),and those of Y_sigma at
most3sqrt(g). The norm triangle inequality follows from Holder applied
to (|X|+|Z|)|X+Z|^3 and division by the cube of the L4 norm (zero norm
is trivial). For four factors,Holder itself follows by Cauchy-Schwarz
twice: E|abcd|<=(E a^4 E b^4 E c^4 E d^4)^(1/4).
Telescoping the four factors in any degree-four monomial shows its
expectation changes by at most4delta(3sqrt(g))^3. Each Q has sum of
absolute coefficients four. Thus
|E Q_+(Y)-E Q_+(S)|,|E Q_-(Y)-E Q_-(S)|
<=432delta g^(3/2)<(432/40000)g^2<g^2/80. (D5')
The expectation margins are -g^2/50 and the perturbations are less than g^2/80,so both actual expectations are strictly negative. If the actual majority
XOR were identically zero,Y_0Y_1>=0 pointwise and Q_+(Y)>=0,impossible.
If that XOR were identically one,Y_0Y_1<=0 and Q_-(Y)>=0,also impossible.
The weak product inequalities include all ties: the evaluator outputs one
exactly for a negative sign sum,and zero for a nonnegative sum.
Thus both actual majority-XOR values occur. No sixth or eighth moments,
small-tail inference,or independence of the full sign vector is used.
All sums and collision events above are actual deterministic quantities.
For an affine S=a+V of dimension at least k,take X_1 uniform on V and
X_2 independently uniform on a+V. Their sum is uniform on S,and D(x)+D(tau(x)) depends only on that sum.
The joint-law contradiction therefore proves both values occur on S. Thus the property holds
for every invertible affine tau without fixed points on S,and includes
every nonzero translation.
## 2. Fix every constant before seeing n
Fix M=2^40,Q_S as fixed above,beta=50/7-2^-100,q_walk=17,C=72,p=144,
ACTUAL outer degree t=4C-1=287 and j=lg(t+1)=lg(t)=9.
For actual n>=1 compute L=floor(log2 n),S=floor(L/Q_S),
s=2ceil(beta*S/2)<=8S,0<=delta_v=s-beta*S<2,
R(s)=0,
ellc=lg(s+4),r_F=S,
m_walk=s+2q_walk*r_F<=48S+320ellc,
s_tag=3m_walk<=144S+960ellc,A=A_rows=2^s_tag,
delta_O=753316552; target d=floor((Dslot)Om) after Om below.
Here A_rows (or A) is a row count,A_B=1/G_B is the BASE source
coefficient,gamma_B is the exact degree coefficient defined in linear7f,
and G_B is its rank slack. None denotes CZ's gamma or the graph Gamma.
Request the outer breaker at
t=287,output m=1,a=s_tag+2+lg C,w=2,
E=2s_tag+22=6s+204S+22<=252S+22. (O2)
We retain the larger coarse coefficient433 below.
Because lg C=7 and13*2^9-8=6648<2^13,its internal exponent is
e=E+13<=433S+5768ellc+54,
Om=Mceil((e+L+54)/M)<=433S+5768ellc+L+M+107.
Thus 2^-E=epsilon=1/(2^22*A^2). This error is independent of the
base rank-test error epsilon0=2^-r_F.
Tags and seeds are internal. R=0<=8ellc,so all retained coarse upper calibrations remain valid. Sampler Section3 proves the walk count under the new S>=2^160 guard.
(O2a) For 2^k<=S<2^(k+1),k>=40,s+4<9S implies ellc<=k+5.
Use the common bound J=max(111+lg C+lg(26t),146+2lg C)=160.
At k40 the exact inequality is
2*5768*45+160=519280
<2^-20*(2139/1000)*2^40.
The left side grows by11536 when k grows by1; its positive intercept
makes this increase at most its value at k40,while the right side doubles.
Induction proves
2*5768*ellc+J<2^-20*L
at EVERY allowed integer S,since L>=Q_S*S. This auxiliary comparison is retained. For Om use
the separate tighter bound below. Here s<=8S;the unit ceiling costs at most M-1.
For Om alone start the dyadic induction at k=160. Its exact base is
5768*165+M+107 < 10^-15*(2139/1000)*2^160.
The left side increases by5768 per step,at most its positive base;the
right side doubles. Thus for EVERY integer S>=2^160,
5768*ellc+M+107 < 10^-15*L.
This proves c_O below.
The separate k80 bases 1067<10^-18*2^80 and817+Q_S<10^-18*2^80 prove
rho-S<eps_B*S and rho0-(Q_S+1)S<eps_B*S in linear7f,eps_B=10^-18.
The row-stretch base 48q_walk*45+6=36726<p*2^40 and the same induction
give 48q_walk*ellc+6<=pS. Hence
s_tag<=2pS<=(2p/Q_S)L. (O2b)
This proves the runtime bound.
Constants precede n; each dimension uses its specified integer rounding.
## 3. The outer affine breaker A# and (F')
By(O2a),Om=Mceil((e+L+54)/M)<=c_O*L,where
c_O=(3p+1)/Q_S+1+10^-15=435139000000000002139/2139000000000000000.
Use wrappers2.5 for its affine-left source guard below;ordinary left
tuples instead use k_gen=k_O+e,restoring every guard of2.1.
Construct wrappers2.1 DIRECTLY at t=287,j=9,C=NB,exponent e,output1
and advice a=s_tag+2+lg C. Use U18 and Slot from B9. All coefficients
are fixed before n. Om is a multiple ofM=2^40; thus dy=(1954763634203/137438953472)Om and r=(4522683721114867/1099511627776)Om are integers,as are dx=(43389108270635/1099511627776)Om and all floor slots. The unit ceiling adds at most M-1 to e+L+54. In the following definitions c_P,pi_P,Lam_P,g_P are the exact7h
coefficients for L0,Lr from7h-alpha at xi=alpha,gamma_1=1,gamma_i=alpha(i>=2); for E use7h at w_N with g_E=alpha*g_generic,pi_E=alpha*pi_generic (c_E,Lam_E unchanged). Alpha is B9's bound:
SC=4257477769389152692011/1000000000000000000,KC=1233515101920841359870607/500000000000000000,S0=34275175387631717613/1000000000000000000,Dslot=1851528021673451056293079/500000000000000000,Kslot=2811521039718751218794657/1000000000000000000.
Equivalently SC=Slot(DNB),KC=Slot(max(K00bar,Kotherbar)),
S0=Slot(c_L0*SC+pi_L0),
TRbar=288*(S0+SC+18*(1954763634203/137438953472 +43389108270635/1099511627776)),TLbar=288*(SC+18*(1954763634203/137438953472)),
Dslot=Slot(max(288*(S0+SC)+KC+alpha,TRbar+Lam_E*(43389108270635/1099511627776)+g_E+287*(43389108270635/1099511627776)+alpha,
S0,Lam_E*(43389108270635/1099511627776)+g_E,KC)),
Kslot=Slot(max(TLbar+Lam_Lr*(4522683721114867/1099511627776)+g_Lr+287*(4522683721114867/1099511627776),Lam_L0*SC+g_L0)).
Actual integer dimensions are
dy=(1954763634203/137438953472)Om, dx=(43389108270635/1099511627776)Om, r=(4522683721114867/1099511627776)Om,
s_C=floor((4257477769389152692011/1000000000000000000)Om),s0=floor((S0)Om),
k_Cslot=floor((1233515101920841359870607/500000000000000000)Om),d_O=floor((Dslot)Om), k_O=floor((Kslot)Om).
Use PXTMloc(JL0=(3,4,4,4,4),vL0) for L0,PXTMloc(JLR=(3,3,4,4,4),vLr) for Lr,PXTM(JE=(3,4),vE) for E and the supplied
uniform-seed AGHP extractor for Lm(input r,output dy,exponent e),padded
and retained to seed dx. B9 supplies NB's complete seed and both source branches.
This gives
delta_O=ceil(Dslotc_O)=753316552,
kappa_O=ceil(Kslotc_O)=571950657.
Every source and input guard is checked separately. Since H>=delta_O
by(O3),d_O<=753316552L<=HL<=n. Also r<=k_O<=kappa_O L<=HL<=n
(k_O>=r by k3b); wrappers2.1(3) has no k_O>=d_O requirement.
E's actual input is d_O,and dKE gives K_E<=d_O. NB's source is n'=d_O; Lr,L0 read n.
L0/Lr use PXTMloc(n,z,e,d_O),a_px=e+14. Since lg(4n)<=L+3,
rho_1<=e+L+49<=Om-5; later rho_i=e+16+2lg(4d_O)<=w_N-32<=alpha*Om.
Public d_O precedes every call. Pure envelopes give K_0<=k_O<=min(d_O,n),
K_0<2^24Om from Kslot<min(Dslot,2^24). As a_px>=246S,Om<=440S,
Om<=2rho_i. All selected R_i>=j_i,so b_i<36(2^24+1)+2<2^30.
Thus every pre-field guard passes before using condensation.
For outer L0,Lr use7h-alpha at xi=alpha,w=Om,gamma_1=1,gamma_i=alpha(i>=2). Here a_px<=e+14,h<=5.
B9 gives alpha*Om-E>=eps_N*(e_B+Q_S)S+67alpha-47(1-alpha).
Since e=E+13,e_B+Q_S>1 and67alpha-47(1-alpha)>13,
alpha*Om-a_px>=eps_N*S-14. On2^j<=S<2^(j+1),j>=160,
lg S<=j+1; base262<eps_N*2^160 and doubling with slope1 prove
eps_N*S>=100+lg S. Hence the gap is>=86+lg S>=6.
Since Om<=440S,lg Om<=9+lg S,so4gap>=20+lg Om+h+1.
Exact0<alpha<1 completes every extra hypothesis before the refined
source/seed envelopes; the new local ambient calls have the same strong error.
E reads d_O,exactly NB's input n'. Its generic7h unit is w_N:
e+lg(B_E)+6+2lg(4d_O)<=e+48+2lg(4d_O)<=w_N,
and w_N>=2^80. Also dx=(43389108270635/1099511627776)Om<=(43389108270635/1099511627776)*440S<=256*12*246S<=256*12w_N.
Its terminal/height/rate guards are unchanged. Apply7h before normalizing:
K_E<=Lam_E*dx+g_generic*w_N,Seed_E<=c_E*dx+pi_generic*w_N.
B9's independently proved w_N<=alpha*Om gives the E coefficients above.
This uses a legal unit and cap before either envelope; no seed guard is circular.
The7g/7h guards require R_h<=64,kap<=505 and output<=256B*Om: every selected
terminal rate and the outputs h,a,s (F#) and s_C,r,dx (outer) obey it.
For Lm,r*dy+1<2^33*Om^2 by (4522683721114867/1099511627776)*(1954763634203/137438953472)+1<2^33.
Hence lg(r*dy+1)<=33+2lg Om. The exact base390<10^-18*2^80,
slope4<=390 and doubling prove66+4lg Om<10^-18*Om for EVERY Om>=2^80.
Its exact seed2[lg(r*dy+1)+dy+2e] fits
(2*(1954763634203/137438953472)+4+10^-18)Om<=(43389108270635/1099511627776)Om,
and its threshold dy+2e<=(1954763634203/137438953472+2)Om<=r. Every original seed and
padded suffix remains retained by the AGHP strongness proof.
R0=3270198837531/1099511627776 is fixed; RATE=289DY+R0. B9 gives w_N>=Ee=e+4 andw_N<=alpha*Om,so e<=alpha*Om-4. AGHP has q=e andKm=dy+2e. Thus(t+1)dy+Km+q=289dy+3e<=(289DY+3alpha)Om<=r,byR0>=3alpha. This proves3a before extraction;the stronger stated Lm seed bound remains valid.
The causal same-side ledgers of wrappers2.3 are
TR=288(s0+s_C+18(dy+dx))<=99242972577449094712841664507/65536000000000000000000Om,
TL=288(s_C+18dy)<=10648653622126077507570706131/8192000000000000000000Om;
the full transcript still retains the final (t+1)dx bits; no seed is dropped.
Let(c_P,pi_P,Lam_P,g_P) be these normalized outer coefficients:7h-alpha for L0,Lr,and7h at w_N for E.
The finite row right sides are listed below. For a fractional slot A,
U18(rhs)+eps_N<=A implies floor(A*Om)>=rhs*Om because eps_N*Om>=1.
Exact-ratio slots use rhs<=A directly; their bit lengths are integers becauseM dividesOm. B9 gives q=e<=alpha*Om.
dy_E: rhs <= 1954763634203/137438953472. Here rhs=c_E*(43389108270635/1099511627776)+pi_E.
dx_Lr: rhs <= 43389108270635/1099511627776. Here rhs=c_Lr*(4522683721114867/1099511627776)+pi_Lr.
dx_Lm: rhs <= 43389108270635/1099511627776. Here rhs=2*(1954763634203/137438953472)+4+10^-18.
s0_L0: U18(rhs)+eps_N <= S0. Here rhs=c_L0*SC+pi_L0.
r_3a: rhs <= 4522683721114867/1099511627776. Here rhs=289*(1954763634203/137438953472)+3alpha.
r_dy: rhs <= 4522683721114867/1099511627776. Here rhs=1954763634203/137438953472.
r_Km: rhs <= 4522683721114867/1099511627776. Here rhs=1954763634203/137438953472+2.
d_3c: U18(rhs)+eps_N <= Dslot. Here rhs=288*(S0+SC)+KC+alpha.
d_3d: U18(rhs)+eps_N <= Dslot. Here rhs=TRbar+Lam_E*(43389108270635/1099511627776)+g_E+287*(43389108270635/1099511627776)+alpha.
d_s0: U18(rhs)+eps_N <= Dslot. Here rhs=S0.
d_KE: U18(rhs)+eps_N <= Dslot. Here rhs=Lam_E*(43389108270635/1099511627776)+g_E.
d_kC: U18(rhs)+eps_N <= Dslot. Here rhs=KC.
k_3b: U18(rhs)+eps_N <= Kslot. Here rhs=TLbar+Lam_Lr*(4522683721114867/1099511627776)+g_Lr+287*(4522683721114867/1099511627776).
k_K0: U18(rhs)+eps_N <= Kslot. Here rhs=Lam_L0*SC+g_L0.
The exact outer coefficients use7h-alpha for JL0,vL0,JLR,vLr and normalized7h at w_N forJE,vE.
verify.py also checks upward 10^-40 coefficient enclosures.
All floor slots include the explicit rounding buffer; all dimensions are fixed before source data,and d3c includes q=e. Together with B9 these cover every source,
seed,output,public cap,relative-field and padding guard. The proof of
wrappers2 gives error(13*2^9-8)2^-e=6648*2^-e<2^13*2^-e=2^-E,
with arbitrary within-side correlations and arbitrary additive right shifts.
All original d-bit seeds remain retained by the following consequence of
the relevant(F') or(F_aff') universal right-side quantifiers.
Fix the admissible left tuple A. For an honest raw seed y and a deterministic
profile r=(B0,...,Bt,Y1,...,Yt),let Delta(y,r) be the TV between the output
tuple over A and an independent uniform honest output with the
unchanged tampered-output tuple. For each y choose a maximizing profile r_y
from this finite domain. Uniform y with profile r_y is a legal instance
of the relevant(F') or(F_aff') contract,
since the whole right side is independent of A and every retained seed is a
function of y. Retention of the honest seed gives
2^(-d_F)*sum_y max_r Delta(y,r) <= 2^-E.
For an actual full right-side variable R independent of A jointly with the
uniform honest raw seed y,conditioning on (y,R) fixes a profile. The joint
output distance retaining (y,R) is the average of its Delta(y,r) values,
and is thus bounded by the displayed average of maxima. Thus the
WHOLE right side,including all padded seed suffixes,is retained free.
## 4. One fixed all-length guard and entropy coefficient
Put eps=2^-18,eps_B=10^-18,eps_R=2^-80 and fix
T=1,
H_raw=4C*delta_O+kappa_O+T=217527117634,
H_base=ceil(A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R)=195869460589,
H=max(H_raw,H_base)=217527117634. (O3)
Here H>=delta_O,which supplies the separate d_O<=n guard.
Here4w=8<=L=T*L. Thus k_F+4w+4Cd
<=[kappa_O+4C*delta_O+T]L<=H*L.
AF4 separately proves k_S=K_B<=H_base*L<=H*L. The final eps_R in
H_base pays the K_B ceiling. This proves(D1) at k=HL,including all4C
linear seed leakages;the affine-fiber argument pays no entropy tails.
For integer floors S>=L/Q_S-1 and 123>3Q_S,L>=p+128 imply
123S>=123L/Q_S-123>3L-123>=2L.
Since linear7f proves s>=7S,s_tag=3(s+34S)>=123S>2L,and every ordered pair of L-bit tag blocks has
exactly A/2^(2L) preimages among the COMPLETE row indices. The tag
collision proof uses A itself,with no stronger lower power of n.
Even s,m_walk make Lambda_S=2^(3m_walk/2-17) an integer. Also
m_walk>=41*2^160>134217728,paying the incidence characteristic guard.
The affine-only sampler is used only on uniform affine B_0 and B_1.
The total evaluator tests these literal finite guards:
n>=1; L=floor(log2 n)>=max(2048,p+128);
S=floor(L/Q_S)>=2^160; s_tag>=2L; HL<=n. (5'')
It also computes/tests every actual AF1 and primitive length/seed/entropy
parameter guard,returning zero on any failed public check before setup.
AF1--AF4,B9 and all fourteen outer requirements prove that those guards
pass under (5'').
Since S tends to infinity at EVERY integer actual length,all fixed
lower thresholds eventually hold. Finally n>=2^L and2^L/L tends to
infinity,its consecutive ratio being2L/(L+1)>=4/3 for L>=2.
At L=128,H*128<2^128 by exact arithmetic. The same increasing ratio
proves HL<=2^L<=n at EVERY L>=128. Thus S>=2^80 already supplies
this native guard,and gives ONE common fixed threshold N for ALL
integer n>=N. All exceptional lengths,including n=0,use the zero branch.
## 5. The fixed evaluator
Precompute the constants in §2 by their finite algorithms. On an
input x of actual length n:
1. If n=0,output0. Compute L=floor(log2 n),S=floor(L/Q_S).
If L<max(2048,p+128) or S<2^160,output0 before any logarithm of S.
Otherwise compute the remaining public parameters. If s_tag<2L,
HL>n or any exact supplied parameter guard fails,output0.
2. Build the exact base extractor E_n from (L_A),seed-padded to s bits.
Its data argument is x and its output length is d.
3. Build the (G_A) graph on left s_tag-bit strings: the four rows of the AB+C
step,the short-start Cayley walk and the s-bit prefix (sampler-and-
disperser.md §4). Its degree is exactly C.
4. Build the canonical field of degree L and the tags in §D1.
For each alpha in lexicographic {0,1}^{s_tag} and j=1,...,C compute
y_(alpha,j)=E_n(x,Gamma(alpha,j)),
b_(alpha,j)=F_n(x,Prefix(y_(alpha,j),d_F),
alpha || encode(j-1) || H_alpha(x)),
using the ONE breaker parameter setting (O2). Let R_alpha be the XOR
of these C bits.
5. Output 1 iff strictly more than half of the A row bits are one;
otherwise output zero. This is exactly §D1's evaluator.
F_n is the raw d_F-seed breaker,identical to its padded form.
The affine base sampler tests density epsilon0=2^-r_F. By AS2,for every
uniform affine source of dimension>=K_B and every fixed test of that
density,fewer than Lambda_S rows are wholly bad except on mass<=2/5.
Source-linearity survives graph index selection. The degree and advice lengths match §D1,and (D1) was checked
above. §§D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.
## 6. Uniform polynomial time in actual input length
One-time constant computation terminates independently of n; parameter
arithmetic adds O(log n) bits. There are
A*C<=C*n^ceil(2p/Q_S)=C*n^135
breaker calls. Affine base degrees b0,b=O(log n) permit2^O(log n) canonical
searches; towers build its initial O(n)-dimensional extension. NB and inner dimensions are O(log n). Outer L0/Lr start on n bits;
their O(n) ambient extensions use7a towers,and their matrices and
unit-vector evaluation have polynomial cost. The graph prime p_m has m+1 bits,
m=s_tag/3,and is found in2^O(s_tag); evaluating all A rows and the walk
thus costs polynomial time in actual n.
By construction,each breaker call is polynomial in n for fixed t=4C-1,a,E=O(log n).
Tags cost2^O(L); all A evaluations are polynomial in n.
For every inner and outer PXT/PXG/PXTM leaf,the explicit131072w bound
makes its N,output ell and accuracy a all O(log n). This includes outer
L0 on the original n-bit input; its terminal dimensions depend on the
logarithmic output and unit. The canonical matrix family takes polynomial
time in2^N,2^(ell+4a+lg(4(N+1))),N,ell,a,hence fixed polynomial time in n.
Its integer weights have polynomial bitlength. B9 also gives LinExt
Rrow,m0,e_L=O(log n),so the same generator bound applies. There are polynomially many
calls,so their composition remains polynomial; n^135 bounds row count
only,and is not asserted as the entire evaluator's time exponent.
Majority uses an O(log A)-bit counter. The total evaluator is
polynomial in ACTUAL n.
## 7. Native circuit bound,NP eligibility and expression
The assembly is nonconstant under every fixed-point-free affine permutation
on a source flat of dimension at least k=HL. In particular,for every
nonzero v and every affine k-flat S its translation derivative is
nonconstant. We prove a finite bridge under the WEAKER hypothesis
for every nonzero v and affine k-flat S,some x in S has f(x)!=f(x+v).
For integers k>=1 this implies C_fullB2(f)>=4n-4k. For n<=k this is trivial.
Otherwise use the retained affine-map,three-corner-pair bridge with d=k-1,
stopping at a=M-q=k. Separation implies nonconstancy on every flat of
dimension at least k,by choosing a nonzero tangent translation of any
hypothetically constant flat.
At every state with a>=k an unused parameter is impossible: for a free
parameter keep its coordinate direction in a safe a-flat; for a paired
parameter fix its partner safely. Independence gives zero derivative
there,and injective Phi preserves its nonzero direction. Absence includes
the output reference. Pair rule1 is thus impossible; rule3 creates an
unused free line coordinate at the same a and is impossible as well.
Mixed and paired/paired productive actions create the unused coordinate
identified in the retained ledgers at new a>=k,so only FF actions occur.
Put W=s-2q. Each FF action removes at least two objects and adds one pair,
dropping W by at least 4 as a drops by one. Basic deletions,transfers,
role exchanges and one-child-row elimination do not increase W. Each
remaining pair rule2 or4 removes at least two objects and one pair,also
not increasing W. Termination thus gives exactly n-k FF actions and
W_T<=W_0-4(n-k)<=|C|-4(n-k).
All parameters lie in the output's undirected incidence component;
otherwise the invertible auxiliary equations split by components and
the output is independent of a parameter,already excluded. In each
component the row/auxiliary block is square: full row and column rank
of invertible A force both inequalities for that block. Each gate-output
variable belongs with its defining object. Thus the output component has
M+s_0 variable vertices and s_0 hyperedges of size at most three. A connected
hypergraph of s_0 such edges has at most 2s_0+1 vertices,since adjoining
an intersecting edge adds at most two. Hence s>=s_0>=M-1. At a=k,
W_T=s-2q>=k-q-1>=-1.
The only integer case with W_T<0 is q=k,M=2k,s=M-1. The tight incidence
tree proof and paired read-once zero-period lemma,proved in
sampler-and-disperser.md,exclude exactly this case under our weaker
separation hypothesis. Thus W_T>=0,proving the finite bound. Its uniform
additive constant zero is optimal for this weaker hypothesis: n=k=1,
f(x)=x satisfies separation and has zero gates.
For all sufficiently large actual integers n,the bridge gives the coarse
bound4n-4H*flg(n); K=4H=870108470536. The refined remainder is proved below.
All constants and exceptional zero branches define one fixed total P
evaluator,with f_0=0. The NP verifier accepts (x,empty) exactly when it
outputs1,and rejects nonempty witnesses (witness bound zero). Tags,seeds
and fields are internal. The coarse remainder is exactlyKlog2n; log2n/n tends to
zero,giving both joint projections(4,1).
The actual stretched row length obeys(O2b):
48q_walk*ellc+6<=pS and s_tag<=2pS<=(2p/Q_S)L. The exact k40
comparison36726<144*2^40 and doubling prove this at every integer S.
Thus every2^O(s_tag) enumeration is polynomial in actual n,including
all row,prime and field work.
# Explicit logarithmic slack pays fixed unit rounding
The evaluator uses beta=50/7-2^-100,s=2ceil(beta*S/2),r_F=S,M=2^40,
with S>=2^160. All constants and rational slots precede actual n.
Put L=floor(log2 n),h=log2 L,a=A_B,b=B_B,Q=Q_S,
u=Dslot,v=Kslot,J=288u+v,E1=1728/7-6*2^-100.
Since s<beta*S+2,E<=E1*S+47,Om<=E1*S+L+M+113.
Also rho<=L/Q+2h+227. Define Aideal=J*(E1/Q+1),
A_tight=Aideal+10^-18,C=(M+113)(J+1). Let c1=497895276875 and
A=c1/4. Exact arithmetic gives4A_tight<=c1<4A_tight+1,so
A-Aideal>=10^-18. Here A is the loss constant,not the row count. Thus
288d+k_O+8<=Aideal*L+C.
Exact arithmetic gives10^-18*Q*2^160>C+2. Since L>=Q*S>=Q*2^160,
this is<=AL-2. The declared coefficient rounds4A_tight upward by less
than one; it includes the paid logarithmic slack,with no additive loss.
For BASE put A0=a*u*(E1/Q+1)+b/Q,B0=2b,C0=(M+113)a*u+227b.
Then a*d+b*rho<=A0*L+B0*h+C0. Exact arithmetic gives Delta=A-A0>0.
On2^j<=L<2^(j+1),h/L<=(j+1)/2^j tends to zero,with successive ratio
at most3/4 forj>=1. Thus Delta*L>=B0*h+C0+3 eventually,and
K_B=ceil(a*d+b*rho)<=AL-2 at EVERY sufficiently large integer L.
Set proof-only k=ceil(AL)-1. It is>=AL-1,so meets both entropy bounds,
and k<=AL. For all large L,1<=k<=AL<=HL<=n. The evaluator does not
use k. The native bridge gives C_B2(f_n)>=4n-4k>=4n-4A*log2(n).
The declared integer log-loss coefficient is c1=497895276875. No constant loss is declared.
A single eventual threshold absorbs all finite guards for every larger
integer actual n,not merely a sparse subsequence. M is a fixed41-bit
integer;all evaluated logarithmic dimensions and uniformP claims persist.
