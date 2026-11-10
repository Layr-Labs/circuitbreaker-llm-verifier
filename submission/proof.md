# Provenance
Ports: DX=6145558683576305895489/1180591620717411303424,DY=5050055860307202828353/590295810358705651712,RATE=2933071716900772837762655/1180591620717411303424.

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
BASE uses7f;L0/Lr/Lm use7f-local.
Sampler Part SP proves finite rings and complete incidence: gamma=((2-nu)/(14+(h_R+5/2)/D_E)),
prime guard2^(2^90). Restricted supplies prime/encoding,CZ smoothing,
conditional tools and the four-round antichain8.2 theorem and its five-label shortened corollary8.4 retaining all9
seeds. Wrappers4.2 instantiates the zero-shift8-copy F# in NB4.6 step7.
Wrappers2,4 prove the all-seeds affine transformation and NB. Outer3 builds
A# at t287,B9 pays its full NB call. Sampler3--6 prove W_A/G_A,the source-linear sampler and
the general conditional lift. Sampler0 resolves legacy
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
explicitguardsk_C<=u,d2<=u,(PG),(AG1) s>=a+40,(AG2) eta*(r+a)>=(905694466927/1000000000)(a+sig+2) of wrappers4.3
and works at every advicelengtha>=1 passing them; failed guards give 0^z.
NB's innerF# instead has the four-round antichain8-copy tuple contract of wrappers4.2
and restricted8.4 at five fixed labels; NB4.6 step7 uses exactly that
zero-shift contract. The outer assembly retains arbitrary additive shifts.
### 3.2 The selected outer affine assembly
The outerA# is the evaluator of wrappers2.1 withC=NB,the rank guards and error of2.5,
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
(AG2) eta*(r+a)>=(905694466927/1000000000)(a+sig+2) it proves all Lemma U and E6 guards (wrappers(I1),(I2)).
Its advice row error is at most14*2^-sig<2^-Ee.
## 5. The standard breaker NB and Lemma B9
Theorem NB (affine-and-standard-wrappers.md §4) is a standard breaker
in the sense of §3.1,with guards k_C<=n,d2<=n,(PG),(AG1) and (AG2),at every a,z,E>=1.
It uses F#,the Section4 advice condenser,canonical LinExt/W2,Ext0 and Ext2; the outer
assembly uses NB only through A#,and Lemma B9 gives its seed and threshold.
**Lemma B9 (NB at the ACTUAL outer call).** Source n'=d_O,output
z=ceil((DX)Om),advice a=s_tag+2+lg C,exponent e,M=2^70,Om=e+L-43. Put
eps_N=10^-18,eps_alpha=10^-40,e_B=1728/7-6*2^-79,alpha=e_B/(e_B+Q_S)+eps_alpha.
The exact slots below have d_O=floor((Dslot)Om)<2^30Om,and O2a gives
Om<=c_O L<1024S: c_O*Q_S*(1+2^-80)<1024 and S>=2^80 (5'').
Hence 4n'<2^32Om<2^42S and lg(4n')<=42+lg S. With Ee=e+4,the actual NB unit satisfies
w_N=max(2^80,Ee+44+2lg(4n'))<=alpha*Om.
Indeed Om=E+L-25 and E<=e_B*S+47,L>=Q_S*S imply
alpha*Om-E >= eps_alpha*(e_B+Q_S)S+-25alpha-47(1-alpha).
The constant is greater than-47;the positive slope reserve pays it. The nonfloor part of w_N is at most
E+151+2lg S. The dyadic certificate
200+2*160<eps_alpha*(e_B+Q_S)*2^160
dominates 198+2lg S (ceil-log included) for every S>=2^160:the left side doubles
on each dyadic step while the positive right bound increases by2.
Since s>=7S,e>=246S+30 and L>=2048,we have Om=e+L-43>=246S.
Exact246alpha>1 gives alpha*Om>=S,paying the floor2^80.
The inner unit w_F=max(2^80,Ee+2lg(4n')+32,ceil(z/128))<=alpha*Om,
because its nonfloor terms are bounded by w_N,while
2^80<=S<=alpha*Om pays the new fixed floor,and
ceil(z/128)<=z/128+1<=alpha*Om follows from 50/128+eps_N<=alpha.
All seven inner leaf-log budgets are<=8; NB
leaf-log budgets14,5 obey Ee+depth+6<=Ee+44. All generic component proofs apply.
The new NB guards are discharged without using Ext0's seed bound first.
Since Om<=c_O*L and L<Q_S*(S+1),the exact scalar check
c_O*Q_S*(1+2^-80)<440 proves Om<=440S. Also w_N>=Ee>=246S.
Thus DX=DX and DX*440<=128*246 imply ceil(z/128)<=w_N;
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
For the final call,Ee>=e1*S+29 with e1=1728/7-6*2^-79,
so Ee>=(1728/7-2/10^12)S. The exact tag formula,L<Q_S*(S+1),Q_S<3 give
Om<=(1728/7+Q_S)S+M+116<=(1728/7+Q_S+1/10^12)S,
since(M+116)*10^12<2^160. Exact multiplication gives
(DX)*(1728/7+Q_S+1/10^12)<(1771192802333/62500000000)*(1728/7-2/10^12).
The difference Delta_Z between the right and left sides of the last
strict scalar inequality satisfies Delta_Z*2^160>1 exactly. Thus
z=ceil((DX)Om)<=1771192802333w_F/62500000000,
including its extra ceiling unit.
These bounds precede all seed envelopes.
For each integer slot ceil(A*w_F),use ceil(A*w_F)<=A*alpha*Om+1.
Thus m1<=D_F*Om and k_F<=K_F*Om,with
D_F=alpha*d_inner+eps_N,K_F=alpha*k_inner+eps_N.
Let (c_P,pi_P,Lam_P,g_P) be the exact PTM.1/PTM.2 for both arity vectors coefficients at vNB0,vNB2.
Put b=6-eta/2. Ext0 outputs v=6r-floor(eta*r/2)<=b*r+1,so
the finite PTM.2 bound is d0<=bc0*r+p0*w_N+D0+c0,
p0=502191152047999/2500000000000,D0=0. Define
M0=c2*D_F+pi2*alpha+eps_N,
K2bar=Lam2*D_F+g2*alpha+eps_N,
h=864/7-3*2^-79,Amin=h/(2h+Q_S),
Rnew=(M0+p0*alpha+(23/4)*alpha+(56+D0)eps_N)/(theta-bc0)+eps_N-theta*Amin/(theta-bc0),
d0bar=bc0*Rnew+p0*alpha+(D0+2)eps_N,
d1bar=M0+2d0bar+6alpha+eps_N,
DNB=2d0bar+2d1bar+D_F+K2bar+alpha.
All denominators are positive and p0>=16. Exact substitution gives
ceil(M0)=67, ceil(K2bar)=1734,
 ceil(Rnew)=164066, ceil(d0bar)=265,
 ceil(d1bar)=602, ceil(DNB)=4257.
Ee<w_N<=alpha*Om and 56+D0<=(56+D0)eps_N*Om give
(23/4)Ee+56+D0<=(23/4)*alpha*Om+(56+D0)eps_N*Om; the r0 ceiling includes -bc0*a in its numerator. The raw length is r=r0-a.
Since a>=hS+9,Om=2a+L-31,L<Q_S(S+1),L>31,
a/Om>=(hS+9)/((2h+Q_S)S+Q_S-31)>=Amin. The second
cross-product difference is h(49-Q_S)+9Q_S>0; denominators are positive.
The discounts sum to theta*a/(theta-bc0),
paying the displayed theta*Amin/(theta-bc0) before invoking Ext0. Its output cap still follows from r<=r0.
Advice guards. By (O2),e=2s_tag+30,so2a=2s_tag+18=Ee-16 EXACTLY;
hence a<Ee/2. With sig=Ee+4,w_N>=Ee+44 and p0>=16:
(AG1) wrappers(I1) gives s>=m0+bc0*r+p0*w_N+D0+c0+3Ee+41>=3Ee+40>6a+40>=a+40,
termwise for every S; at S=2^80,s/(a+40)>500.
(AG2) Ee>=S>=2^160 gives w_N=Ee+44+2lg(4n'). Ext2 has B=30 and rho_2>=w_N-33,so its two slices total at least (1-2^-40)*sum_i(R_i+1)*w_N. Put a_px=Ee+lg(B)+4. Then w_N-a_px<=119+2lg S<=S/1099511627776<=w_N/1099511627776: on2^j<=S<2^(j+1),j>=160,the left bound is121+2j,starting441<2^160/1099511627776 and growing by2 while the right doubles. Also w_F>=w_N-12>=(1-2^-40)w_N,hence m1>=(1-2^-40)*d_inner*w_N.
The new leaf seed is ell+4a_px+lg(4(m_h+1)),with ell>=m1/B.
Including it gives m0>=SD*w_N,
SD=(1-2^-40)*(sum_i(R_i+1)+d_inner/B+4).
Here33<=2^80/1099511627776 pays both slice and w_F rounding.
For this actual untrimmed caller r_F=S. Set W=1728/7+10^-6,
G=2717083400781/2000000000 and Gplus=G+10^-12.
The exact tag formula gives E<=(1728/7-6*2^-79)S+47.
With Ee=E+22 and lg(4n')<=42+lg S,
w_N=max(2^80,Ee+44+2lg(4n'))<=(1728/7)S+256+2lg S<=W*S.
For2^j<=S<2^(j+1),j>=80,the last overhead is at most258+2j;
its base418<2^80/10^6 and increment2 are dominated by doubling.
The maximum's floor is also at most W*S. Since Ee>=246S,
the floor is inactive and w_N-Ee<=128+2lg S<=S/1099511627776<=w_N/1099511627776.
Indeed its dyadic bound130+2j starts290<2^80/1099511627776 and grows by2.
Thus Ee>=(1-2^-40)w_N,without using any extractor guarantee.
From2a=Ee-16 and sig=Ee+4 we have the exact identity
2(a+sig+2)=3Ee-4. Since w_N>=Ee+44 and S>=w_N/W,
the following conservative upper bound remains valid:
(905694466927/1000000000)(a+sig+2)<=G*w_N-(8151250202343/2000000000).
The r0 formula subtracts bc0*a but retains its positive(23/4)Ee term.
Since a<Ee/2<=w_N/2,
eta*(r+a)=eta*r0>=eta*(p0+SD+(23/4)*(1-2^-40)-bc0/2)*w_N/(theta-bc0)>=Gplus*w_N.
The last rational inequality is checked exactly. Hence(AG2) holds
for every legal S;10^-12*2^80>6*(905694466927/1000000000) also pays the old finite constant.
At S=2^80 the actual ratio eta*(r+a)/((905694466927/1000000000)(a+sig+2))>1.000000000001.
The exact inequality (15404805937759/905694466927+5*2^-80)*eta<=gamma pays the new energy budget. Both guards hold at the only NB call (outer A#),so Theorem NB applies there.
The r0 ceiling costs <=eps_N*Om and a>=Amin*Om,so r_NB=r0-a<=Rnew*Om. Both source branch coefficients are at most KC-eps_N,
KC=1230045614569197952756989/500000000000000000. These exact inequalities are recomputed from the displayed recurrences.
For LinExt put J=lg(4(Rrow+1)). Since Rrow>=3,m0>=1,
4(Rrow+1)<=(Rrow*m0+1)^2,hence J<=2lg(Rrow*m0+1).
Now 2e_L<=d0+3Ee and ceil(2*(Rnew+2)*M0)+1<2^240 give
2lg(Rrow*m0+1)<=480+4lg Om<=eps_N*Om. The k80 base is
804<eps_N*2^80,with slope4<=804; doubling covers every Om>=2^80.
Thus d1<=m0+2d0+6Ee+eps_N*Om<=d1bar*Om and d_C<=DNB*Om.
Define U18(x)=ceil(10^18*x)/10^18 and Slot(x)=U18(x)+2eps_N.
Put SC=Slot(DNB)=1064150499197349597843/250000000000000000,KC=Slot(max(K00bar,Kotherbar))=1230045614569197952756989/500000000000000000,
where K00bar=b*Lam0*Rnew+g0*alpha+eps_N and
Kotherbar=K_F+10M0+DX+2alpha+eps_N. The extraeps_N pays
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
The affine-breaker contract (F_aff') is proved below from the outer A# (§3) and Lemma
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
Fix Q_S=213938239825351/100000000000000. With L=floor(log2(actual input length)),S=floor(L/Q_S)
is an integer; S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer rounding.
(L_A) Put L=floor(log2 n),S=floor(L/Q_S)>=2^160,s=2ceil(beta*S/2),
ellc=lg(s+4) and r_F=S. Linear Section7f supplies the literal
three-condenser/kernel base E_n,source-linear at each seed and strong
for every affine source of dimension at least
K_B=ceil(A_B*d+B_B*rho),rho=r_F+225+2lg S,
d=floor((Dslot)Om),with Om defined below and d<=750003150L.
Its complete seed length D<=s,error<2^-r_F and public dimensions are
AF1. The full AF2--AF4 proof gives K_B<=216406774195L<=HL.
(G_A) Sampler Sections3--5 use q_walk=17,C=72,m_walk=s+34r_F,
s_tag=3m_walk<=144S+960ellc,A=A_rows=2^s_tag and
Lambda_S=A^(1/2)/2^17. The full four-row AB+C incidence proof followed
by W_A expands every Lambda_S-row set to more than5*2^-r_F*2^s
neighbors.
The evaluated sampler is source-linear and costs2^O(s_tag).
(F_aff') For every n passing the public guards,the selected outer A#
uses wrappers2.5 at the slots below,with d_F=d_O<=delta_O L and
k_F=k_O<=kappa_O L. For A_i=T_i Z+c_i with uniform affine Z,invertible
T_i,independent of an arbitrary correlated tuple(B,Y),Hinf(Z)>=k_F,
honest Y0 uniform d_F-bit,and honest advice distinct from every tampered
advice,the honest one-bit output on A_i+B_i is2^-E-close to uniform
jointly with all tampered outputs and EVERY seed. No functional-copy
assumption is made on(B,Y). Addition is XOR. For fixed t and advice
length O(log n),the complete evaluator is polynomial in actual n.
No ordinary-source outer theorem is asserted at these slots.
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
Use the source-linear sampler Samp(x,alpha,j),j in [C],of sampler Section5,
seed length d,affine base entropy requirement k_S=K_B and error2^-r_F.
Only uniform affine sources are used. The field index alpha replaces[A].
For t=4C-1 and a=s_tag+w+lg(C),request the padded all-seeds-retained affine
breaker F of output length one,affine-left threshold K,and error at most
epsilon = 1/(2^12 A^2).
The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1,and w-bit H_alpha(x). Define
R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),advice(x,alpha,j)),
D(x)=1 iff 2 sum_alpha R_alpha(x)>A.
For the affine application below,assume the source entropy satisfies
k >= max(k_S,K+4w+4Cd). (D1)
All concatenations,zero padding and majority ties are specified. There
are no seed or advice inputs to D. The tags are internal computations on
the actual n-bit x.
## D2. Interface
Linear left leakage and tags leave affine fibers. D3 uses(F_aff') on
each fiber and enumerates right profiles without conditioning the right law.
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
<=521A^2/6144<3gA/35,by epsilon=1/(2^12*A^2),w=2 andg>=255A/256.
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
Fix M=2^70,Q_S as fixed above,beta=50/7-2^-79,q_walk=17,C=72,p=144,
ACTUAL outer degree t=4C-1=287 and j=lg(t+1)=lg(t)=9.
For actual n>=1 compute L=floor(log2 n),S=floor(L/Q_S),
s=2ceil(beta*S/2)<=8S,0<=delta_v=s-beta*S<2,
R(s)=0,
ellc=lg(s+4),r_F=S,
m_walk=s+2q_walk*r_F<=48S+320ellc,
s_tag=3m_walk<=144S+960ellc,A=A_rows=2^s_tag,
delta_O=750003150; target d=floor((Dslot)Om) after Om below.
Here A_rows (or A) is a row count,A_B=1/G_B is the BASE source
coefficient,gamma_B is the exact degree coefficient defined in linear7f,
and G_B is its rank slack. None denotes CZ's gamma or the graph Gamma.
Request the outer breaker at
t=287,output m=1,a=s_tag+2+lg C,w=2,
E=2s_tag+12=6s+204S+12<=252S+12. (O2)
We retain the larger coarse coefficient433 below.
Because lg C=7 and9*3^9-4=177143<2^18,its internal exponent is
e=E+18<=433S+5768ellc+59,
Om=e+L-43<=433S+5768ellc+L+M+116.
Thus 2^-E=epsilon=1/(2^12*A^2). This error is independent of the
base rank-test error epsilon0=2^-r_F.
Tags and seeds are internal. R=0<=8ellc,so all retained coarse upper calibrations remain valid. Sampler Section3 proves the walk count under the new S>=2^160 guard.
(O2a) For 2^k<=S<2^(k+1),k>=40,s+4<9S implies ellc<=k+5.
Use the common bound J=max(111+lg C+lg(26t),146+2lg C)=160.
At k40 the exact inequality is
2*5768*45+160=519280
<2^-20*(213938239825351/100000000000000)*2^40.
The left side grows by11536 when k grows by1; its positive intercept
makes this increase at most its value at k40,while the right side doubles.
Induction proves
2*5768*ellc+J<2^-20*L
at EVERY allowed integer S,since L>=Q_S*S. This auxiliary comparison is retained. For Om use
the separate tighter bound below. Here s<=8S;the old M-1 allowance is retained only as a conservative upper estimate.
For Om alone start the dyadic induction at k=160. Its exact base is
5768*165+M+116 < 10^-15*(213938239825351/100000000000000)*2^160.
The left side increases by5768 per step,at most its positive base;the
right side doubles. Thus for EVERY integer S>=2^160,
5768*ellc+M+116 < 10^-15*L.
This proves c_O below.
The separate k80 bases 1067<10^-18*2^80 and817+Q_S<10^-18*2^80 prove
rho-S<eps_B*S and rho0-(Q_S+1)S<eps_B*S in linear7f,eps_B=10^-18.
The row-stretch base 48q_walk*45+6=36726<p*2^40 and the same induction
give 48q_walk*ellc+6<=pS. Hence
s_tag<=2pS<=(2p/Q_S)L. (O2b)
This proves the runtime bound.
Constants precede n; each dimension uses its specified integer rounding.
Fixed gamma_0=500267319014381/250000000000000,gamma_R=530242185449/250000000000,xi_0=gamma_0/(gamma_0-1),G_0=1-(xi_0-1)^2,
A_0=1/G_0,B_0=xi_0*gamma_0/G_0. ModeR has xi_R=gamma_R/(gamma_R-1),G_R=1-(xi_R-1)^2,A_R=1/G_R,B_R=xi_R*gamma_R/G_R.
## 3. The outer affine breaker A# and (F_aff')
By(O2a),Om=e+L-43<=c_O*L,where
c_O=(3p+1)/Q_S+1+10^-15=43513938239825351213938239825351/213938239825351000000000000000.
Use wrappers2.5 for the affine-left rank guard and error below.
Only the supplied affine-left wrapper is selected;no arbitrary-left guarantee is asserted.
Construct wrappers2.1 DIRECTLY at t=287,j=9,C=NB,exponent e,output1
and advice a=s_tag+2+lg C. Use U18 and Slot from B9. All coefficients
are fixed before n. Om=e+L-43 is an integer working unit. M=2^70 is only the fixed
port coefficient denominator. Each port below is rounded up individually;
all floor slots retain their original rule. In the following definitions c_P,pi_P,Lam_P,g_P are the exact7h
coefficients for L0 are(0,1+2gamma_0*alpha,A_0,B_0*alpha),for Lr(0,1+2gamma_R*alpha,A_R,B_R*alpha);for E use7h at w_N with g_E=alpha*g_generic,pi_E=alpha*pi_generic (c_E,Lam_E unchanged). Alpha is B9's bound:
SC=1064150499197349597843/250000000000000000,KC=1230045614569197952756989/500000000000000000,S0=4967752109312344033/1000000000000000000,Dslot=230464019264325469731889/62500000000000000,Kslot=1995463379848692560200133/1000000000000000000.
Equivalently SC=Slot(DNB),KC=Slot(max(K00bar,Kotherbar)),
S0=Slot(1+2gamma_0*alpha),S_R=Slot(1+2gamma_R*alpha),
TRbar=288*(S0+SC+18*(DY +DX)),TLbar=288*(SC+18*(DY)),
Dslot=Slot(max(288*(S0+SC)+KC+alpha,TRbar+Lam_E*(DX)+g_E+287*(DX)+alpha,
S0,Lam_E*(DX)+g_E,KC)),
Kslot=Slot(max(TLbar+A_R*(RATE)+B_R*alpha+287*(RATE),A_0*SC+B_0*alpha)).
Actual integer dimensions are
dy=ceil((DY)Om), dx=ceil((DX)Om), r=ceil((RATE)Om),
s_C=floor((1064150499197349597843/250000000000000000)Om),s0=floor((S0)Om),
k_Cslot=floor((1230045614569197952756989/500000000000000000)Om),d_O=floor((Dslot)Om), k_O=floor((Kslot)Om).
**Individually rounded ports.** The normalized table is unchanged.
Exact rational arithmetic gives DY-c_E*DX-pi_E>16/M,
DX-c_Lr*RATE-pi_Lr>16/M,and RSHIFT-alpha>16/M.
For Om>=246*2^160,their products with Om exceed2^20,which pays
all extra seed ceilings and the288 dy ceilings in r3a. In particular
r>= (t+1)dy+floor(KM*Om),dy>=Seed_E and dx>=Seed_Lr whenever the pure
primitive envelopes apply. Every lambda for L0,Lr,Lm,E is below2^20.
Changing each port from its normalized real value to its ceiling changes
TR,TL and the d3c,d3d,k3b arithmetic right sides by less than2^30:
t=287,j=9,so the transcript additions are at most4j(t+1)+t,
and each changed source-output term adds at most its lambda. Each Slot
has a surplus of at least eps_N*Om before its final floor. Exact
eps_N*246*2^160>2^31 pays these additions and that floor. S0 uses
the same Slot margin. NB's k_C expression gains at most one from z;
its existing eps_N term pays the finite constant4 instead of3.
These are pure length comparisons;no extractor guarantee is invoked here.
The primitive guards below are checked before invoking those guarantees.
Thus every displayed caller guard
holds for these actual integer dimensions. The sampler seed d_O and
source threshold k_O are still floor slots;their combined bound remains
4C*d_O+k_O+8<=(288Dslot+Kslot)Om+8.
Use AInit of linear7f-local for L0/Lr,PXTM(JE=(2,),vE) for E and AInit7f-local for Lm(input r,output dy,exponent e). All three AInit calls receive the same outer(n,L,Om);their actual data ambients are n,n,r. Padded dx-bit seeds are retained. B9 supplies NB's complete seed and both source branches.
This gives
delta_O=ceil(Dslotc_O)=750003150,
kappa_O=ceil(Kslotc_O)=405866994.
Every source and input guard is checked separately. Since H>=delta_O
by(O3),d_O<=750003150L<=HL<=n. Also r<=k_O<=kappa_O L<=HL<=n
(k_O>=r by k3b); wrappers2.1(3) has no k_O>=d_O requirement.
E's actual input is d_O,and dKE gives K_E<=d_O. NB's source is n'=d_O; Lr,L0 read n.
All L0/Lr/Lm use AInit7f-local on uniform affine fibers,at modes gamma_0,gamma_R,gamma_M.
Its complete pre-field proof gives K0<2^24Om for all three output requests.
Each K0 includes the fixed xi/G additive allowance;Slot surplus pays it.
Initial padding is S0=Slot(1+2gamma_0*alpha);Lr uses S_R=Slot(1+2gamma_R*alpha)
and DX>S_R,so it ignores its suffix beyond s_R,not s0. Initial L0 entropy is
paid by k_K0;Lr's actual source threshold is paid by k_3b after all affine
transcript and tampered-output constraints. Its strong affine error2^-e
therefore supplies the exact-fiber wrapper proof without an ordinary-left
assumption. All three selected L maps use AInit;the unselected local library is omitted.
E reads d_O,exactly NB's input n'. Its generic7h unit is w_N:
e+lg(B_E)+6+2lg(4d_O)<=e+48+2lg(4d_O)<=w_N,
and w_N>=2^80. Since DX*Om<=DX*440S<=256*2*246S<=256*2w_N,
and the last bound is an integer,dx=ceil(DX*Om)<=256*2w_N.
Its terminal/height/rate guards are unchanged. Apply7h before normalizing:
K_E<=Lam_E*dx+g_generic*w_N,Seed_E<=c_E*dx+pi_generic*w_N.
B9's independently proved w_N<=alpha*Om gives the E coefficients above.
This uses a legal unit and cap before either envelope; no seed guard is circular.
The7g/7h guards require R_h<=64,kap<=505 and output<=256B*Om: every selected
terminal rate and generic outputs h,a,s (F#) and dx (E) obey it.
Lm=AInit(N_amb=r,z=dy,e;N_outer=n,L_outer=L,Om) of linear7f-local.
gamma_M=525060396321/200000000000;SM=S_M.
Set S_M=Slot(2gamma_M*alpha),KM=Slot(A_M*(DY+eps_N)+B_M*alpha),s_M=floor(S_M*Om),k_M=floor(KM*Om). Lm ignores its dx-bit seed after s_M;the ENTIRE seed is retained.
The three-mode proof pays the dy ceiling via A_M*eps_N*Om>=A_M,and xi/G,K0 ceiling/floor via Slot surplus. Hence literal K0<=k_M<=r before fields;ambient r<=n was proved above. For Lm the initial condenser is identity since r<=2^12Om<U1;only the two final independent seeds are used,and2b<=2gamma_M*alpha*Om+2 fits s_M<dx by Slot surplus before extraction.
Exact (RATE-288DY-KM)>16/M and Om>=246*2^160 pay all288 dy ceilings: r>=288dy+k_M. This is the rank guard(3aff) on every affine fiber,not an entropy tail. Fixed-seed source-linearity and full padded-seed strongness are7f-local's.
This proves3a before extraction.
The causal same-side ledgers of wrappers2.3 are
TR=288(s0+s_C+18(dy+dx))<=18277115034237853680561685512565759167/14073748835532800000000000000000Om+2^30,
TL=288(s_C+18dy)<=44692989301668757446011251251574476777/35184372088832000000000000000000Om+2^30;
the full transcript still retains the final (t+1)dx bits; no seed is dropped.
Let(c_P,pi_P,Lam_P,g_P) be the normalized coefficients:for L0 use(0,1+2gamma_0*alpha,A_0,B_0*alpha),for Lr(0,1+2gamma_R*alpha,A_R,B_R*alpha),and for E7h at w_N.
The table lists normalized rows. U18(rhs)+eps_N<=A implies
floor(A*Om)>=rhs*Om+eps_N*Om-1. The rounded-port lemma bounds each
actual-row increase by2^30,paid by eps_N*Om>2^31. Port lengths are
integers by their ceilings,not divisibility. B9 gives q=e<=alpha*Om.
dy_E: rhs <= DY. Here rhs=c_E*(DX)+pi_E.
dx_Lr: rhs <= DX. Here rhs=0*(RATE)+(1+2gamma_R*alpha).
dx_Lm: rhs <= DX. Here rhs=SM.
s0_L0: U18(rhs)+eps_N <= S0. Here rhs=1+2gamma_0*alpha.
r_3a: rhs <= RATE. Here rhs=288*(DY)+KM.
r_dy: rhs <= RATE. Here rhs=DY.
r_Km: rhs <= RATE. Here rhs=KM.
d_3c: U18(rhs)+eps_N <= Dslot. Here rhs=288*(S0+SC)+KC+alpha.
d_3d: U18(rhs)+eps_N <= Dslot. Here rhs=TRbar+Lam_E*(DX)+g_E+287*(DX)+alpha.
d_s0: U18(rhs)+eps_N <= Dslot. Here rhs=S0.
d_KE: U18(rhs)+eps_N <= Dslot. Here rhs=Lam_E*(DX)+g_E.
d_kC: U18(rhs)+eps_N <= Dslot. Here rhs=KC.
k_3b: U18(rhs)+eps_N <= Kslot. Here rhs=TLbar+A_R*(RATE)+B_R*alpha+287*(RATE).
k_K0: U18(rhs)+eps_N <= Kslot. Here rhs=A_0*SC+B_0*alpha.
The exact selected coefficients use7f-local for L0/Lr/Lm and normalized7h at w_N forJE,vE.
verify.py also checks upward 10^-40 coefficient enclosures.
All floor slots include the explicit rounding buffer; all dimensions are fixed before source data,and d3c includes q=e. Together with B9 these cover every source,
seed,output,public cap,relative-field and padding guard. The proof of
wrappers2.5 gives error(9*3^9-4)2^-e=177143*2^-e<2^18*2^-e=2^-E,
with arbitrary within-side correlations and arbitrary additive right shifts.
Whole-right retention follows directly from wrappers2.5:at the final
rectangular history the full right tuple is independent of the remaining
left law. Append it before translating the right contributions and taking
prefixes. Hence the bound retains all original seeds and every additional
right-side variable R;the rank-good history event depends only on the
retained history and its error is averaged with the actual marginal.
## 4. One fixed all-length guard and entropy coefficient
Put eps=2^-18,eps_B=10^-18,eps_R=2^-80 and fix
T=1,
H_raw=4C*delta_O+kappa_O+T=216406774195,
H_base=ceil(A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R)=216406773785,
H=max(H_raw,H_base)=216406774195. (O3)
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
Source-linearity survives graph index selection. The degree and advice lengths match §D1,and (D1) was checked
above. §§D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.
## 6. Uniform polynomial time in actual input length
All fixed rational/constant searches terminate before n;parameter arithmetic
has O(log n) bits. There are A*C<=C*n^ceil(2p/Q_S)=C*n^135 calls.
BASE and AInit binary field degrees are O(log n),so their canonical
irreducible searches,condenser evaluation and kernel row reduction cost
polynomial time in actual n. Both outer left maps use the explicitly analyzed AInit.
For every used PXT/PXG/PXTM leaf,131072w bounds its terminal ambient N,
output ell and accuracy a by O(log n),Remaining matrix calls have the stated logarithmic terminal bounds.
The canonical matrix algorithm is polynomial in2^N,
2^(ell+4a+lg(4(N+1))),N,ell,a,hence polynomial in n;weights have
polynomial bitlength. NB's LinExt also has Rrow,m0,e_L=O(log n).
Graph primes have O(s_tag)=O(log n) bits and searches take2^O(s_tag);
all row/walk evaluations,tags and O(log A)-bit majority counting are
polynomial. Polynomially many fixed-depth calls preserve polynomial time.
The n^135 bound counts rows,not the entire evaluator time exponent.
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
bound4n-4H*flg(n); K=4H=865627096780. The refined remainder is proved below.
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
## 8. Exact logarithmic loss without a fixed loss
Keep the single fixed evaluator above,with the new unit Om=e+L-43.
Put Q=Q_S,E1=1728/7-6*2^-79,J=288Dslot+Kslot and
A=J*(E1/Q+1). The exact declared logarithmic coefficient is
c1=684182439186457328558933673490230439058128908659049362379903523/1381262385173351218365413725936025600000000000000000=4A. There is no constant monomial in bound.json.
Since s<beta*S+2 and e=6s+204S+30,
Om=6s+204S+L-13<E1*S+L-1<=(E1/Q+1)L-1.
Both entropy floor slots therefore satisfy
4C*d_O+k_O+8<=J*Om+8<A*L-J+8<A*L-1,
as the exact fixed J>9. This pays the raw source requirement,including
all left tags/seed leakage,without a fixed positive loss.
For BASE let a=A_B,b=B_B,u=Dslot. Its entropy requirement is bounded by
A0*L+B0*log2(L)+C0,where
A0=a*u*(E1/Q+1)+b/Q,B0=2b,C0=(M+113)a*u+227b.
This remains an upper bound with the smaller unit. The exact rational
margin A-A0 is positive. Thus there is one integer j0 such that for all
j>=j0,(A-A0)2^j>=B0(j+1)+C0+2:choose a base large enough;the
left doubles and the right has fixed increment B0. For2^j<=L<2^(j+1),
this gives the BASE requirement<=A*L-2. Hence it holds for EVERY
sufficiently large integer L,not only the dyadic sequence.
Choose proof-only k=ceil(A*L)-1. Then k>=A*L-1 pays both source
budgets,and k<=A*L. Also A<H-1,so eventually1<=k<=H*L<=n.
D3 applies whenever those two explicit source budgets hold,so its full
affine-permutation property holds at this k. The evaluator never receives k.
The native bridge gives C_B2(f_n)>=4n-4k>=4n-4A*L
>=4n-c1*log2(n) for every sufficiently large ACTUAL integer n.
All public guards and exceptional zero branches are those proved above;
the small-length branch returns zero before the remaining parameters.
The fixed evaluator uses the stated unit and primitive modes. Its exact
logarithmic coefficient is explicitly declared above;no fixed loss term
occurs. No comparison routine or checker verdict is a
mathematical premise. Since log2(n)=o(n),the joint projection is(4,1).

Material source:771814cc-1e2c-4f72-8e53-b9e1a4125e14;historical credits remain above. Selected parameters and the seed-margin proof are fixed above. Verdicts are not premises.
