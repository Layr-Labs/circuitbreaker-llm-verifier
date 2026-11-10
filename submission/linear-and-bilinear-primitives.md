# linear-and-bilinear-primitives.md
Fix Q_S=213938242726479/100000000000000. With L=floor(log2(actual input length)),S=floor(L/Q_S)
is an integer; S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer rounding.
Component references resolve through the map in sampler-and-disperser.md.
# Finite seeded and linear extraction: an elementary logarithmic-regime replacement
All logarithms are binary; distance means half-L1.
## 1. Selected interfaces
Sections3--4 supply the strong linear condenser,exact block-source repair
and shared hash. Sections7g--7h give complete selected recurrences and
finite seed/entropy envelopes;7a supplies canonical exact fields,and7f
supplies the affine base and7f-local the selected L0/Lr. 
## 2. Deterministic fields,with no irreducible-polynomial search
Construct nested fields F_i of binary dimension 2^i. Start with F_0=F_2 and
a_0=1. Given F_i and an element a_i of absolute trace 1,adjoin z_i satisfying
z_i^2+z_i+a_i=0,
a_(i+1)=a_i z_i.
Here absolute trace in a field of binary dimension v means
sum_(j=0)^(v-1) x^(2^j). Finite-field powering x^(2^v)=x follows from
Lagrange's counting argument on the nonzero multiplicative group. Relative
trace in a quadratic extension of a field of size Q is x+x^Q; grouping
the Frobenius sum proves trace transitivity.
This really is a quadratic field extension. The map x -> x^2+x is F_2-linear,
has kernel {0,1},and its image has half the field's elements. Every image
has absolute trace zero,because Tr(x^2)=Tr(x). Trace is nonzero by the
induction hypothesis Tr(a_i)=1; thus its zero kernel also has half the
field's elements. The two sets coincide,and a_i is not an image.
The polynomial has no root and is irreducible. Its two roots are z_i,z_i+1,
so the relative trace of a_i z_i is a_i. Hence a_(i+1) has absolute trace 1.
This proves the induction.
Use ordered tower coordinates,with constants embedded by appending zero
coordinates. Addition is XOR. For x=x0+x1 z_i and y=y0+y1 z_i compute
p=x0*y0,q=x1*y1,r=(x0+x1)*(y0+y1),t=a_i*q in the preceding field.
The product is (p+t)+(r+p)z_i,because z_i^2=z_i+a_i. These four recursive
products give T(v)<=4T(v/2)+O(v),hence O(v^2) bit operations. Inverses can be computed
by exponentiation to 2^v-2 in a field of binary dimension v,using O(v)
squarings/multiplications. Zero inversion is never requested.
We also require a polynomial-basis representation of an extension of degree
l,where both l and the base-field binary dimension s are powers of 2 and
l>=2. Build the tower to binary dimension sl. Let g be the last adjoined
root. It is not in the immediately preceding subfield. Its orbit under
q-th powering,q=2^s,has length exactly l. Indeed the orbit length divides
l; a proper divisor of the power of 2 l is at most l/2. Then g would satisfy
x^(q^(l/2))=x. The preceding subfield has exactly q^(l/2) elements and fills
the roots of that polynomial,contradicting the choice of g.
Consequently
E(T)=product_(j=0)^(l-1) (T-g^(q^j))
has coefficients fixed by q-th powering,whose roots are exactly F_q.
Every F_q-polynomial vanishing at g also vanishes at all l orbit points,
so no nonzero polynomial of degree below l can vanish there. Thus E is
irreducible of degree l over F_q. Compute that product in the
tower and read its coefficients in the embedded base field. This constructs
E in polynomial time in sl and gives a deterministic canonical choice.
No factoring,primitive-element search,or enumeration of 2^(sl) field
elements is involved. All choices,including coordinate orders,are fixed.
## 3. A finite strong linear lossless condenser
This section fixes integers U>=1,1<=k<=U,a>=1 and R>=1. Let
rho=a+2 ceil(log U)+2,
s=the least power of 2 at least (R+1)rho,
u=s-rho,q=2^s,H=2^u,
l=max(2,the least power of 2 at least ceil(U/s)),
b=ceil(k/u),v=bs.
Use §2 for F_q and an irreducible E(T) of degree l. Pad the U-bit
input with zeros to sl bits and view it as the coefficient vector of
f(T),degree <l,over F_q. A seed is y in F_q. Output the b field elements
Con(f,y)=(f(y),(f^H mod E)(y),...,
(f^(H^(b-1)) mod E)(y)). (4)
The output has v bits. Each map f -> f^(H^j) mod E is F_2-linear,because
H is a power of 2. Evaluation and padding are linear too. Moreover
s < 2(R+1)rho,u>=Rs/(R+1),
v <= (1+1/R)k+s. (5)
We prove that Con is a strong lossless condenser: for every source X of
min-entropy k,(Y,Con(X,Y)) is within 2^(-a) of a distribution (Y,Z) with
Y uniform and H_infty(Z | Y=y)>=k for every y.
### Expansion proof
Write A=q-(l-1)(H-1)b. The graph with neighbors (y,Con(f,y)) expands every
set of J<=H^b inputs to at least AJ neighbors. To prove this,suppose a
set T of fewer than AJ neighbors contains every neighbor of some inputs.
For j=0,...,J-1,write j in base H with b digits,and let M_j(Z_0,...,Z_(b-1))
be the corresponding monomial. There are AJ coefficients in a polynomial
Q(Y,Z_0,...,Z_(b-1))=sum_(j<J) p_j(Y) M_j(Z_0,...,Z_(b-1)),
deg p_j < A.
The fewer than AJ conditions of vanishing on T admit a nonzero solution.
Choose one with least degree in Y. Not all p_j are divisible by E: otherwise
dividing by E preserves vanishing on T (E has no F_q root because l>=2) and
lowers that degree.
For any input f whose neighbors all lie in T,substitute the polynomials
f_j=f^(H^j) mod E. The resulting polynomial in Y has degree at most
(A-1)+(l-1)(H-1)b=q-1 and vanishes at every y in F_q. It is zero. Reducing
modulo E thus says that f,as an element of F_q[T]/E,is a root of
sum_(j<J) (p_j mod E) Z^j.
This is a nonzero polynomial of degree <=J-1,so there are at most J-1 such
inputs. In particular a J-input set cannot have fewer than AJ neighbors.
The same statement holds after restricting the source to the padded inputs.
Here l<=2U and b<=k<=U. Also q/H=2^rho>=4 U^2 2^a. Thus
((l-1)(H-1)b)/q <= 2U^2/2^rho <= 2^(-a-1).
Thus A>0 and A>=(1-2^(-a))q. Since H^b>=2^k,the required size-2^k
sets have at least (1-2^(-a))2^k q distinct neighbors.
### Strongness,not only joint min-entropy
For a source uniform on a size-K set,K=2^k,consider each y separately.
Keep one preimage for each distinct value Con(x,y). To obtain a K-element
set of output values at that y,replace the duplicates by unused values.
There is room,because v>=k (s>=u and b>=k/u). The number replaced,
summed over all y,is at most 2^(-a)Kq by expansion. The repaired conditional
output at each y is uniform on K values,and the seed marginal is unchanged.
The joint distance is at most 2^(-a). An arbitrary k-source is a convex
combination of uniform size-K sources: iteratively express a probability
vector with coordinates <=1/K as a mixture of vertices of that capped
probability simplex; its vertices have exactly K entries 1/K. Mixing the
repairs keeps each conditional atom <=1/K and preserves the distance bound.
This proves the asserted strong condenser,including the seed marginal.
## 4. Three exact block operations
A block source (X_1,...,X_j) has block entropy k if every conditional law
X_i | X_<i=x_<i has min-entropy >=k on its support. All statements below
also hold conditioned on any fixed external variable.
**Apply a strong condenser with one shared seed.** If Con has error eta,
apply Con(X_i,Y) to all j blocks with the same uniform independent Y.
The output tuple,retaining Y,is within j eta of a block source with the
condenser's guaranteed conditional entropy. A fully safe coupling argument
is as follows. Retain the original source blocks as latent variables. For
each original prefix and seed,choose the condenser's repaired conditional
output law and a coupling to the actual output. Generate these repairs
causally,with independent fresh coupling randomness. Given an original
prefix and the seed,previous repairs carry no additional information about
the next original block. Each repaired next block has the desired entropy;
mixing over latent prefixes preserves the atom bound. Each repair disagrees
with its original output with probability at most eta on average,because
the original prefix is independent of the uniform seed. A union bound proves
the distance. No coupling is part of the evaluator.
**Split one high-entropy block.** Pad a v-bit k-source by one zero if needed
and split into two m-bit halves,m=ceil(v/2). If k>=m+k'+a,the pair is within
2^(-a) of a block source with entropy k' in each half. The first half already
has entropy >=k-m>=k'+a. The average maximum probability of the second half
given the first is at most 2^(m-k). Markov bounds by 2^(-a) the mass of first
halves whose conditional second-half entropy is below k'. Replace the
second-half law by uniform only on those prefixes. The first-half law is
unchanged. Applied to j input blocks with the causal latent-prefix coupling
above,the total distance is at most j2^(-a). This proof also explains why
arbitrary conditioning on a rare prefix is not silently allowed.
**Hash every block with one shared seed.** For m-bit blocks,put W equal to
the least power of 2 at least m. Use the tower field of binary dimension W.
Pad each block to W bits. For a uniform z in that field define Hash_z(x)
as the first ell coordinates of zx. For distinct x,x',z(x-x') is uniform,
so collision probability is exactly 2^(-ell).
For a single source of entropy >=ell+2a,its conditional output collision
probability,averaged over z,is at most 2^(-ell)+2^(-ell-2a). Cauchy-Schwarz
and Jensen give strong TV distance <=(1/2)2^(-a). For a j-block source this
yields joint jell output bits with error <=(j/2)2^(-a),retaining z. Replace
outputs from the LAST block backwards: the extractor bound can retain the
entire previous source prefix,which is independent of z. Project that
prefix to its prior hashed outputs,then repeat. This avoids conditioning
on already hashed outputs as though the shared seed were still independent.
This hash is linear for each z and its seed length is exactly W<=2m.
## 5. Arithmetic cost and selected outer runtime
For the exact selected recurrences,m_i>K_i implies K_(i-1)>j_i*K_i,
so B<=K_0<=u. Also b_i>=1 and K_(i-1)>((j_i-1)/j_i)S_i give
S_i<2u; every later input and leaf hash has dimension<=u. These guards
precede field setup. For the Section3 condenser,with U=4u,the least
relative degree l=max(2,pow2ceil(ceil(U/S_i))) satisfies
S_i*l<=2U+2S_i=O(u). The orbit polynomial uses
O(l*S_i) squarings and O(l^2) products in that tower. Its coefficients
occupy the initial base-field coordinates,so no expensive basis search
is required. Polynomial multiplication/reduction,repeated squaring and
Horner evaluation cost polynomially many bit operations once the base
fields are constructed. There are at most B nodes per level and a fixed
number of levels at every selected call.
Exact base/hash field searches need 2^O(v)poly(v) time at degree v;
no polynomial bound in arbitrary component input length is asserted.
At the actual outer calls these degrees are O(log actual N),by the
literal envelopes and closed slots. Their searches thus cost
N^O(1). Larger relative extensions use the explicit polynomial-cost
towers,including the affine base's initial O(N)-dimensional extension.
All other arithmetic,matrices,kernels and output evaluation are
polynomial in actual N. Sections7a/7f and proof.md supply the complete
all-length guards and the fixed total evaluator's exceptional branches.
## 7. Exact multiway block-source repair
For a positive integer j,pad a v-bit k-source to jm bits,m=ceil(v/j),
and split into Z_1,...,Z_j. If k>=(j-1)m+k'+a,the tuple is within
(j-1)2^-a of a block source with entropy k' per block. Indeed the prefix
of length i has entropy>=k-(j-i)m. Averaging its largest next-block atom
across earlier prefixes gives at most2^((j-1)m-k). Markov bounds prefixes
of conditional entropy<k' by2^-a. Generate a repaired tuple sequentially:
use the original given good prefixes,otherwise uniform on m bits.
Its conditional entropy is>=k',since k<=jm implies k'<=m; the union bound
over the j-1 possibly bad steps gives the error. The first block is already
good since its marginal entropy is>=k-(j-1)m>=k'+a. Conditional external
variables and causal repairs obey Section4's same argument. This lemma
is used at each selected j_i by7h,and at j=2 by7g.
## 7a. Exact binary field degrees for the base and both NB primitives
Condenser/AGHP fields use exact binary degrees.
For each positive degree v
choose the first monic irreducible binary polynomial by lexicographic
enumeration and the modular-squaring/gcd test proved in weak-seed Section4
later in this file. That proof supplies existence and correctness at every
v. Quotient fields use polynomial coordinates and schoolbook arithmetic.
The trace is a nonzero binary linear map: sum X^(2^j),j=0,...,v-1,
is a nonzero polynomial of degree2^(v-1)<2^v,so cannot vanish everywhere.
Test the ordered basis units and choose the first of trace1. Starting
from that base field,use Section2's Artin--Schreier tower. Its trace
transitivity and kernel/image proof give irreducibility and propagate
trace1 at every quadratic step. Only the relative extension degree is a
power of two. The last root's full Frobenius orbit and its orbit product
give the required condenser modulus. Larger extensions are thus
built by towers,not enumeration in their full dimension. Section4's optional field hash
uses canonical degree W and the same collision proof.
Section3's expansion requires only the stated field arithmetic and degree
inequalities; Sections7g/7h give each selected exact recurrence,padding,
error and complete-seed proof directly.
F#/NB inputs/outputs/error and selected output/error/base/leaf degrees
are O(log actual N) by B9; outer L0/Lr inputs can be N bits.
BASE7f has b0,b=O(log N),U1=2^(110+lg S)=O(log N); towers build its initial
O(N)-dimensional relative extension. BASE's later ambient/matrix/kernel dimensions
are O(log N). Canonical searches cost2^O(log N); tower/matrix operations
are polynomial in actual N,without field or source advice. No joint bound
for arbitrary component dimensions is asserted.
## 7f. A three-seed affine base from linear condenser duality
Only this BASE is specialized to uniform affine sources. All PXT,NB and
internal breaker source/right-tuple contracts remain general-source.
An affine k-source is uniform on a translate of a binary k-dimensional
vector subspace. Write Con[U,k,rho,b] for Section3's actual polynomial
map at integer field degree b>rho: Hpow=2^(b-rho),relative degree
l=max(2,pow2ceil(ceil(U/b))),output length b*ceil(k/(b-rho)). Pad shorter
inputs by zeros. Section7a supplies the canonical exact-degree field and
relative power-of-two tower. The polynomial maps f -> f^(Hpow^i) mod E
are binary linear,as are padding and evaluation. Thus this is source-linear.
Section3's expansion and seed-preserving repair proof uses these integers
and the defect (l-1)(Hpow-1)*ceil(k/(b-rho))/2^b,not an integer-rate
hypothesis. For k<=U,rho=a+2lg(U)+2 and b-rho>=1,l<=2U and
ceil(k/(b-rho))<=k<=U bound the defect by 2^(-a-1). The same proof
thus gives strong lossless condensation with error <=2^-a.
Rank preservation. For an affine source of dimension k,a source-linear
condenser's output at seed y is uniform on a translate of a rank-r_y
subspace,r_y<=k. If r_y<k,every repaired conditional law with atom cap
2^-k gives this support mass <=2^(r_y-k)<=1/2. Its conditional TV distance
is thus >=1/2. Since the repair retains the uniform seed marginal,
joint TV is the average of these distances. At most a 2eta fraction of
seeds can fail injectivity on any fixed k-dimensional direction space
when the strong condenser error is eta.
Kernel duality. Let V<=F_2^v have dimension k and W=V-perp under the
binary coordinate dot product. Suppose G:F_2^v -> F_2^w is injective
on W. Choose any d independent vectors from ker G as matrix rows M.
If M restricted to V did not have rank d,a nonzero row combination
would annihilate V and hence belong to W intersect ker G,a contradiction.
A rank-d linear image of a uniform affine V-source is exactly uniform
on d bits,since all fibers are cosets of the same kernel. This argument
allows G to have dependent rows; it assumes no surjectivity of G.
For actual n>=1 fix Q_S as fixed above,beta=50/7-2^-79,delta_O=750003139,
eps_B=10^-18,eps_R=2^-80,N0=2^60,a0=1+2^-60,
The guard S>=2^80 pays each unit rounding by eps_R*S>=1.
gamma_B=(50/7-a0*(Q_S+1+eps_B)-3eps_R)/(2(1+eps_B)),
xi_B=gamma_B/(gamma_B-1),G_B=1-(xi_B-1)^2,
A_B=1/G_B,B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)=1155.
First compute L=floor(log2 n),S=floor(L/Q_S),d=0.
If S<2^160 return 0^d BEFORE evaluating any lg S. Otherwise set
these ACTUAL integer parameters in the following order:
s=2ceil(beta*S/2)<=8S,
R(s)=0,
ellc=lg(s+4),r_F=S,a_B=r_F+3,
U=4n,U1=2^(110+lg S),rho=a_B+2lg(U1)+2=r_F+225+2lg S,
M=2^40,Om=6s+204S+L-13,d=floor((Dslot)Om),K_B=ceil(A_B*d+B_B*rho),
rho0=a_B+lg U+lg K_B+2,
b0=ceil((N0+1)rho0/N0),u0=b0-rho0,
b=ceil(gamma_B*rho),u=b-rho,
v0=b0*ceil(K_B/u0),v1=b*ceil(K_B/u),
c=v1-K_B,v2=b*ceil(c/u),D=b0+2b. (AF1)
O2a gives d<=delta_O*L; all calibrations below use this upper bound.
First require S>=2^160. Test K_B<=n,rho0>=N0,u0>=1,u>=1,
v0<=U1,v1<=U1,1<=c<=U1,v1-v2>=d and D<=s. A failed
public dimension guard returns 0^d,before field construction.
Parse the first D seed bits as independent uniform y0,y1,y2 of lengths
b0,b,b. Compute z0=Con[U,K_B,rho0,b0](x,y0),then
z1=Con[U1,K_B,rho,b](z0,y1). Define G_y2 by restricting
Con[U1,c,rho,b] to its first v1 input coordinates through zero padding.
Form its v2-by-v1 binary matrix by evaluating it on the ordered v1
coordinate units. Row-reduce it over F_2,processing columns increasingly,
choosing the first remaining row with a pivot1,and eliminating that column
from every other row. For each free column in increasing order,set that
free coordinate1,the other free coordinates0,and solve pivot coordinates.
These vectors are the canonical ordered kernel basis. Use the first d
vectors as rows M_y2 and return M_y2*z1. Ignore the final s-D seed bits.
Only the PUBLIC matrix kernel is computed.
Initial condensation is legal even at its smaller degree. Since
u0>=rho0/N0>=1,ceil(K_B/u0)<=K_B,and b0/u0<=N0+1,
its output has v0<=(N0+1)K_B+b0 bits. Its defect is bounded by
2U*K_B/2^rho0<=2^(-a_B-1).
The two remaining condensers have thresholds K_B and c,both <=U1,
ambient U1 and rho=a_B+2lg(U1)+2,so their defects have the same bound.
Dimension,including BOTH ceilings. The exact constants satisfy
gamma_B>2,G_B>0,G_B*A_B=1,
G_B*B_B-xi_B*gamma_B>=xi_B/2^40.
Because b>=gamma_B*rho,b/(b-rho)<=xi_B. Consequently
v1<=xi_B*K_B+b,c<=(xi_B-1)K_B+b,v2<=xi_B*c+b.
Using the EXACT identity v1=K_B+c,not a false lower bound on v1,gives
v1-v2>=K_B-(xi_B-1)c-b>=G_B*K_B-xi_B*b
>=G_B*A_B*d+(G_B*B_B-xi_B*gamma_B)rho-xi_B
>=d+xi_B*(rho/2^40-1)>=d. (AF3)
Here ceil(K_B/u)<=K_B/u+1 and b<=gamma_B*rho+1 pay the two
roundings; the outer K_B ceiling helps this lower bound. Also c>0 since
b>u and v1=b*ceil(K_B/u)>K_B. Thus for EVERY y2,rank(G_y2)<=v2
and dim ker G_y2>=v1-v2>=d. The evaluated rows always exist.
For an affine K_B-source,the first condenser preserves its direction
dimension except on fraction <=2*2^-a_B of y0. For every good y0,the
second preserves that dimension except on fraction <=2*2^-a_B of fresh
y1. For a good pair let V be its K_B-dimensional direction space in
F_2^v1. Its perpendicular space has dimension c. Fresh y2 makes G_y2
injective on that space except on fraction <=2*2^-a_B. Kernel duality
then makes the actual output exactly uniform on every good seed triple.
A union bound gives TV <=6*2^-a_B=(3/4)2^-r_F<2^-r_F,retaining the entire seed
triple. The independent ignored suffix can also be retained. For a larger
affine source choose any K_B-dimensional direction subspace only in this
analysis: surjectivity on it implies surjectivity on the full source.
Complete all-length calibration. For 2^k<=S<2^(k+1),k>=80,
s+4<9S implies ellc<=k+5. Thus
rho-S<=10k+267<eps_B*S.
The exact coefficient bound is
K_B/S<=A_B*delta_O*Q_S*(1+2^-40)+B_B*(1+eps_B)+eps_R
<2^48.
Here L<Q_S(S+1),rho<=(1+eps_B)S,and eps_R*S>=1 pay K_B's ceiling.
Thus lg K_B<=48+lg S. Since lg(4n)<=L+3,
rho0-(Q_S+1)S<=9k+Q_S+97<eps_B*S.
The k80 bases are 1067<10^-18*2^80 and 817+Q_S<10^-18*2^80. Their affine left sides
increase by10 and9,at most their bases,while eps_B*2^k doubles.
This proves both bounds for EVERY integer S>=2^80. Each field-degree
ceiling costs <=1<=eps_R*S. By the definition of gamma_B,
D/S<=a0*(Q_S+1+eps_B)+2gamma_B*(1+eps_B)+3eps_R=50/7. (AF2)
Retain gamma_B,which uses budget50/7. Before paying the three ceilings by
eps_R*S,the calculation gives D<=(50/7-3eps_R)S+3. The positive margin
marg=3eps_R-2^-79 obeys marg*2^160>3,so for all S>=2^160,
D<=(50/7-2^-79)S<=s. Also beta>7 gives s>=7S,and s<=8S.
Thus the COMPLETE three-slice seed and ignored suffix fit s.
Put b0coef=a0*(Q_S+1+eps_B)+eps_R. The exact inequalities
(2^60+1)2^48+b0coef<2^110,
xi_B*2^48+gamma_B*(1+eps_B)+eps_B<2^110
imply v0,v1<2^110*S<=U1,and 0<c<v1<=U1. Also
K_B<=[A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R]L
<=H_base*L,H_base=216408572747. (AF4)
The last allowance pays the K_B ceiling. H*L<=n and H>=H_base give
K_B<=n. Finally rho0>=(Q_S+1)S>2S>=2^60=N0 and rho>=S>=2^40 prove every
public AF1 dimension guard and the sharpened final inequality in AF3.
Each enumerated exact field degree b0,b is O(log actual n). The initial
relative extension has dimension O(n) and uses the explicit towers of
Sections2/7a,never enumeration in that large dimension. U1,v0,v1,
v2 and the evaluated public matrix dimensions are O(log n),with enormous
fixed coefficients. Matrix evaluation and binary Gaussian elimination are
polynomial in these dimensions. Thus this is one deterministic polynomial
evaluator in actual n,source-linear at every full seed.
The classical condenser-duality principle is credited to Mahdi Cheraghchi
Bashi Astaneh,Applications of Derandomization Theory in Coding,EPFL thesis
4767 (July2010),Corollary5.11 and following discussion,p.148,
https://eccc.weizmann.ac.il/resources/pdf/astaneh.pdf . The full rank and
kernel proofs needed here were supplied above; the citation is attribution.
## 7f-local. Three affine initializers
Fixed gamma_0=250133067710547/125000000000000,gamma_R=87419532791/40000000000. Mode0 is L0 with output z=s_C and gamma=gamma_0;
modeR is Lr with output z=r and gamma=gamma_R. ModeM is Lm with z=dy,gamma_M=87419532791/40000000000. In each mode set
xi=gamma/(gamma-1),G=1-(xi-1)^2,A=1/G,B=xi*gamma/G,N0=2^60.
Pass public outer N_outer,L_outer=floor(log2 N_outer),Om=e+L_outer-43 and actual data ambient N_amb<=N_outer. All three modes use THIS outer Om;for Lm,N_amb=r,not N_outer. Invalid positive-length/unit guards return0^z before logs. Set L=L_outer,
U=4N_amb,U1=2^(110+lg Om),rho=e+225+2lg Om,
K0=ceil(A*z+B*rho+xi/G),rho0=e+5+lg(4N_amb)+lg K0,
b0=ceil((N0+1)rho0/N0),u0=b0-rho0,b=ceil(gamma*rho),u=b-rho,
v0=b0*ceil(K0/u0),v1=b*ceil(K0/u),c=v1-K0,v2=b*ceil(c/u).
Run the7f canonical condensers,kernel row reduction and first z kernel
rows. Set S_gamma=Slot(1+2gamma*alpha),with
Slot(x)=ceil(10^18*x)/10^18+2*10^-18. L0 seed length is s0=floor(S_0 Om).
Lr takes the first s_R=floor(S_R Om) bits of its dx-bit seed;DX>S_R is
checked and every ignored suffix is retained. Lm takes the first s_M=floor(S_M Om) bits of its dx-bit seed,with S_M=Slot(1+2gamma_M*alpha);DX>S_M. Its threshold is floor(KM*Om),KM=Slot(A_M*(DY+10^-18)+B_M*alpha). All ignored suffixes are retained. Failed public guards return0^z before fields;fixed-seed maps are linear.

Here are all caller guards BEFORE condensation. All three fixed gamma values
are in(2,16],so1<xi<2,0<G<1. The unit bound gives alpha*Om-e>=10^-40*S-66. With
lg Om<=9+lg S,the base313+2*160<10^-40*2^160 and doubling
pay225+2lg Om,so rho<=alpha*Om.
The source rows A_0*SC+B_0*alpha and A_R*RATE+B_R*alpha have their Slot margins. ModeM uses A_M*(DY+10^-18)+B_M*alpha:dy<=DY*Om+1 and10^-18*Om>=1 pay its rounded output. Slot surplus pays xi/G+1<2^12+1 before the floor,as10^-18*Om/4>2^12+2.
The separate seed slot is unchanged.
 The source cap is
K0<2^24Om in all three modes . Thus
lg K0<=24+lg Om and rho0<=e+L+32+lg Om<=Om+75+lg Om<=2Om;
84+lg S<=246S and Om>=246S pay the last inequality.
Since e>=246S,L>2S,S>=2^160,we have rho0>=2^60,rho>=2^40,u0,u>=1.
Since N_amb<=N_outer,lg(4N_amb)<=L+3 gives the same rho0 bound. L0/Lr callers pay K0<=N_outer. For Lm,the literal Km bound above and r>=floor(KM*Om) give K0<=N_amb=r before condensation. Also b0<=2rho0<=4Om and
b<=16rho+1<=16Om+1,so
v0<=(N0+1)K0+b0<2^86Om<U1,
v1<=xi*K0+b<2^26Om<U1,0<c<v1.

The full seed obeys,uniformly over all three modes,
b0+2b<=(1+2gamma)e+L+8000+65lg Om+rho0/2^60.
Indeed rho0's bound contributes32+lg Om,the two later degrees contribute
2gamma*(225+2lg Om)+2,and b0's ceiling contributes at most1;
gamma<=16 makes this smaller than the stated constants.
Since e+L=Om+43 and e<=alpha*Om-4,this is at most
(1+2gamma*alpha)Om+8100+65lg Om+(2/2^60)Om.
Exact2/2^60+10^-18/4<2*10^-18 leaves a quarter of10^-18*Om to pay
8101+65lg Om,including the final floor. For Om>=246S,Om<=440S,
lg Om<=9+lg S,the base at S=2^160 and doubling prove this bound.
Hence the complete seed fits its respective padded length in all three modes.

For kernel dimension,b/(b-rho)<=xi and b<=gamma*rho+1 yield
v1-v2>=G*K0-xi*b>=z+(G*B-xi*gamma)rho+xi-xi=z.
The fixed xi/G in K0 explicitly pays the later degree ceiling;
G*B=xi*gamma exactly in all three modes.
The7f direction-preservation proof now applies:each of three independent
condenser seeds fails on mass<=2*2^(-e-3),the last on the c-dimensional
perpendicular space. Kernel duality gives uniform output on every good
triple. Union error<2^-e retains the full padded seed. Larger affine
sources use a K0-dimensional direction subspace only in analysis.
This guarantee requires an affine source.
The pre-field bounds use no condenser guarantee.
Canonical fields/row reduction are7f's;degrees and condensed dimensions
are O(log actual N_outer),so every call is uniformly polynomial-time in N_outer.

## 7g. Exact finite rational-rate primitives PXT
PXT is exactly PXTM of7h with every j_i=2 and R_i>1.
PXT[h,Ro,Rl] sets R_1=...=R_(h-1)=Ro,R_h=Rl; PXT[h,R] sets all R_i=R.
Write PT0,PT0a,PT1 for this specialization of PTM.0,its capacity/split
inequalities,and PTM.1. Thus B=2^h,a=e+h+4,lambda_i=2R_i/(R_i-1),
and e_i=0. The complete algorithm,guards,strong error and source proof
are7h's; its terminal hash is defined next independently of the recurrence.
**Canonical general-source linear leaf family.** For positive integers N,m,a set
q=2^-m,eps=2^-2a,b=1+eps/2,d=m+4a+lg(4(N+1)),T=2^d.
Construct T binary m-by-N matrices H_s with
Z_x=|{s:H_s x=0}|<Tq(1+eps) for every nonzero x.
Initially take independent uniform matrices and potential F=sum_(x!=0)b^Z_x.
Its expectation is (2^N-1)(1+q*eps/2)^T. For t=eps/2,
ln(1+t)>=t-t^2/2 and ln(1+qt)<=qt give
q(1+eps)ln b-ln(1+qt)>=q*eps^2/4.
Since Tq*eps^2>=4(N+1),the initial expectation is strictly below
b^(Tq(1+eps)): indeed ln(2^N-1)<N ln2<N+1.
Fix entries in seed/row/column order by minimum conditional expectation,
ties choosing0. This is an exact canonical finite algorithm. For each x,
a partially fixed row is either a determined parity or uniform. Thus
p_x=Pr[H_s x=0|prefix] is0 if a determined row is1,otherwise2^-r,
where r counts undetermined rows. Future matrices have the common factor
1+q*eps/2. Set Den=2^(2a+1),Num=Den+1,and maintain integer
W_x=Num^(completed zero count)*Den^(completed nonzero count).
Compare the trial-bit integers sum_x W_x*(Den*2^m+[p_x>0]*2^(m-r)).
Omitted denominators and future factors are common and positive.
At matrix completion multiply W_x by Num or Den according as H_s x=0
or not. Conditional expectation never increases,proving the claimed Z_x.
Weights have O(Ta)bits. The TmN bit choices each enumerate2^N-1 vectors
and inspect mN entries; all arithmetic and table generation have time
polynomial in2^N,T,N,m,a. There is no source enumeration or advice.
Linearity gives collision at most q+q*eps for every distinct pair.
For any source of minentropy k>=m+2a,the mean output collision is at most
2^-k+q+q*eps. Cauchy--Schwarz bounds retained-seed TV from uniform by
.5 sqrt(2^m*(2^-k+q*eps))<=2^-a/sqrt(2)<2^-a.
The same calculation gives mean_seed TV^2<=(1/2)2^-2a.
This is strong extraction for every source,with a source-linear map for
fixed seeds. The table depends only on dimensions.
Apply this family to N=m_h,m=ell with a FRESH shared d-bit seed.
The causal repairs precede this seed Y and give each leaf entropy ell+2a
given the retained condenser seeds and its full repaired prefix.
Reverse-hybrid replace leaves B,...,1,retaining Y,earlier source leaves
and independent later uniform outputs. Apply the single-block average-Y
bound,then project earlier source leaves to their hashes. Each step costs
<2^-a even with Y retained; unused independent suffixes remain retained.
Together with 2(B-1)2^-a this gives total<(3B-2)2^-a<2^-e.
PXG is now this SAME terminal variant as PXT; all thresholds stay PT0/PT1.
**Unified legal seed envelope (PT2).** Require h<=13,z<=256B*w,
terminal arity2<=j<=18,R_h<=64 (binary:253/252<=Rh<=64) and kap=(R_h+1)/(R_h-j+1)<=505,
in addition to K0<=u and w>=max(2^20,e+lg(B)+6+2lg(4u)).
For binary j=2. The terminal ceiling calculation gives
m_h<kap*(ell+3a+j-1)+S_h/j+1
<=kap*z/B+[3kap+(R_h+1)/j]w+1+1/j+(j-18)kap.
Since kap<=505,j<=18 and R_h<=64,505*259+65/2=261655/2<130828<131072. The residual is<=1+1/j<=3/2;w>=2^20 pays it,so m_h<=131072w under the output cap. The binary endpoint253/252 and mixed j=18 endpoint477/28 both give kap=505. No hash dimension or logarithmic seed allowance grows.
Put J=lg(4(m_h+1))<=20+lgw. The actual seed length is
D=sum_iS_i+ell+4a+J<=z/B+sum_i(R_i+1)w+4a+J+h+1.
If u<w/2,PT0a gives K0>rho>=a,so a<w/2 and J+h+1<2w.
Otherwise lg(4u)>=lgw+1,whence w-a>=2lgw+4.
Then4(w-a)-(J+h+1)>=7lgw-h-5>0 since lgw>=20,h<=13.
Both cases prove
D<=c*z+pi*w,c=1/B,pi=4+sum_i(R_i+1). (PT2)
Ignored padding bits remain retained. Outside these extra cap/rate guards
the hash theorem and literal algorithm still hold,but PT2 is not asserted.
Every selected actual call satisfies them,as checked below and in B9.
At every selected call m_h,ell,a=O(log actual n). Even when u=actual n,
the terminal bound is131072w with w=O(log n). Thus2^N,T and the integer
weight bitlength are fixed polynomials in actual n. Construction is uniform
polynomial in actual n,without a joint-polynomial claim in unrelated
component dimensions. Condenser fields are still canonical fields of7a.
## 7h. Mixed-arity exact primitive PXTM
Fix h>=1,integers j_i>=2,and rational R_i>j_i-1.
For integers u,z,e>=1 set B=product_i j_i,ell=ceil(z/B),
a=e+lg(B)+4,U=4u,rho=a+2lg(4u)+2,
S_i=ceil((R_i+1)rho),K_h=ell+2a,
b_i=ceil(j_i(K_i+a+j_i-1)/(S_i-j_i*rho)),m_i=ceil(b_i*S_i/j_i),
K_(i-1)=(j_i-1)m_i+K_i+a,i=h,...,1. (PTM.0)
If K_0>u return0^z before field setup. Otherwise begin on the original
u-bit source. Level i applies Section3's condenser with ambient4u,
threshold K_(i-1),degree S_i and b_i outputs,using ONE fresh S_i-bit seed
shared across that level. Pad each result to j_i*m_i bits and split into
j_i blocks. Hash all B leaves to ell bits with the canonical family of7g,
using one fresh shared seed of length ell+4a+lg(4(m_h+1)),then truncate to z. Retain every seed and ignored suffix.
All fields are canonical exact fields/towers from7a. The leaf count B need
not be a power of two; only relative field extensions use quadratic towers.
Here S_i>j_i*rho and b_i(S_i/j_i-rho)>=K_i+a+j_i-1,so
K_(i-1)<(j_i-1)b_i*S_i/j_i+j_i-1+K_i+a<=b_i(S_i-rho),
b_i<=m_i<K_(i-1),K_(i-1)=(j_i-1)m_i+K_i+a.
Section3's expansion guards hold. Each later input has length
m_i<K_(i-1)<=K_0<=u,so all padding is legal. Section7's exact multiway
repair applies. With M_i=product_(l<i)j_l,the node counts satisfy
sum_i M_i(j_i-1)=B-1 and sum_i M_i<=B-1.
Condenser plus split errors cost at most2(B-1)2^-a,and the shared hash
costs less than B2^-a. Thus total<3B2^-a<2^-e,since B<=2^lg(B).
The causal repair proof of4 retains the complete seed and every latent
prefix. Every operation is source-linear for fixed seeds. This proves
strong linear extraction for arbitrary K_0-sources,not just affine sources.
For integer w>=max(2^20,e+lg(B)+6+2lg(4u)),rho<=w,a<=w-6.
Put lambda_i=j_i*R_i/(R_i-j_i+1),A_h=2,B_h=1/B and backwards
e_i=max(0,(j_i-7)lambda_i+(j_i-1)/j_i),
A_(i-1)=lambda_i(A_i+1)+((j_i-1)/j_i)(R_i+1)+e_i/2^20,
B_(i-1)=lambda_i*B_i. Define Lam=B_0,g=A_0.
The actual two ceilings give
K_(i-1)<[j_i(S_i-rho)/(S_i-j_i*rho)](K_i+a+j_i-1)
+((j_i-1)/j_i)S_i
<=lambda_i(K_i+a+j_i-1)+((j_i-1)/j_i)S_i.
Since a+j_i-1<=w+j_i-7 and S_i<=((R_i+1)w+1),the residual
(j_i-7)lambda_i+(j_i-1)/j_i is at most e_i. Since w>=2^20,
this contributes at most (e_i/2^20)w,exactly the added recurrence term.
For 2<=j_i<=6,e_i=0 and the old envelope is recovered.
With K_h<=z/B+2w-11,induction gives
K_0<=Lam*z+g*w,Lam=product_i R_i/(R_i-j_i+1). (PTM.1)
Let j=j_h and kap=(R_h+1)/(R_h-j+1). The terminal block obeys
m_h<kap*(ell+3a+j-1)+S_h/j+1
<=kap*z/B+[3kap+(R_h+1)/j]w+1+1/j+(j-18)kap.
Under the additional legal seed-envelope guards of7g,its proof applies
verbatim to these ceilings,giving c=1/B,pi=4+sum_i(R_i+1) and
D<=c*z+pi*w. (PTM.2)
The source envelope PTM.1 is unchanged. Generic NB tests BOTH output caps
and all the terminal rate/arity/height guards before claiming PTM.2.
Ext0 selects JS=(4,5,5,5,5,6),h=6,B=15000,lg(B)=14,j=6.
Ext2 selects JY=(5,6),h=2,B=30,lg(B)=5,j=6.
Their complete ordered rate vectors are supplied separately.
The actual unit w_N>=Ee+44+2lg(4n) satisfies the above unit guard.
Output6r and all threshold/input/output/retained-seed requirements are
charged by B9 and wrappersI1--I3. All field degrees and inputs are
O(log actual N); the fixed15000 and30 leaf counts and canonical searches/towers
have polynomial bit time in actual N.
## 10. Attribution and source convention
The PV-code expansion/block method is due to Guruswami,Umans and Vadhan:
https://people.seas.harvard.edu/~salil/research/PV-condenser-June08.pdf .
The recursive linear condensation/shared-seed strategy follows
Chattopadhyay,Goodman and Liao,ECCC TR21-075,Lemmas 4.7--4.10 and
Theorem 5.6: https://eccc.weizmann.ac.il/report/2021/075/download/ .
Chattopadhyay--Liao use a sharper variable-parameter form in Lemma 3.13:
https://arxiv.org/pdf/2110.12652v1 . Xin Li's later application is
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/ .
These finite constructions are proved here. Strong output excludes retained seed bits.
# An elementary weak-seed replacement for the local Raz call
This local interface uses the small-bias construction of Alon,Goldreich,
Hastad and Peralta (1992) and the leftover-hash collision argument.
Sources inspected:
https://www.cs.tau.ac.il/~nogaa/PDFS/aghp4.pdf
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/
(Algorithm 3 and Lemma 5.1).
All logarithms are base two,TV is half the L1 norm,and vector/matrix
arithmetic is over F_2. This proof is finite and elementary. A finite-field existence argument,
deterministic field algorithm and outer-length complexity are supplied in
§4.
## 1. A finite strong linear seeded extractor
Fix input length d>=1,output length m>=1 and integer e>=1. Put
N=d*m,
r=ceil(log(N+1))+m+2e.
Fix one presentation of F_{2^r} as an r-dimensional binary vector space.
For a uniform seed (u,v) in F_{2^r}^2,define N bits
b_i=<bin(u^i),bin(v)> for i=0,...,N-1,
using the ordinary nondegenerate binary coordinate dot product. Place these
bits,in fixed row-major order,into an m-by-d binary matrix M(u,v). Define
E(x,(u,v))=M(u,v)*x.
For every fixed seed this map is binary linear in x. If X is independent of
the uniform seed and H_infinity(X)>=m+2e,then
TV((E(X,U_{2r}),U_{2r}),(U_m,U_{2r})) <= 2^-e. (1)
Proof of small bias. A nonempty parity of the b_i equals the dot product of
bin(P(u)) with bin(v),where P is a nonzero polynomial with binary
coefficients and degree at most N-1. If P(u) is nonzero,the dot product is
balanced over uniform v; if P(u)=0,the sign equals one. A nonzero polynomial
has at most its degree many roots,by repeated division by X-a. Thus the
absolute bias is at most
delta=(N-1)/2^r < 2^(-m-2e).
For every nonzero w in F_2^d and every nonzero a in F_2^m,the parity
<a,Mw> uses a nonempty collection of matrix entries: choose any nonzero
coordinate of a and of w. Consequently the same bias bound applies.
Fourier inversion of the zero event gives
Pr[Mw=0]=2^-m*sum_a E[(-1)^(<a,Mw>)] <= 2^-m+delta. (2)
Let X' be an independent copy of X. For each seed s let p_s be the law of
E(X,s). Averaging its collision probability and separating X=X' yields
E_s sum_y p_s(y)^2
<= Pr[X=X']+(1-Pr[X=X'])*(2^-m+delta)
<= 2^-m+2^(-m-2e)+delta.
Cauchy--Schwarz gives TV(p_s,U_m)<=sqrt(2^m*sum_y(p_s(y)-2^-m)^2)/2.
Jensen's inequality for the square root now bounds the average TV by
(1/2)*sqrt(2^m*(2^(-m-2e)+delta)) < 2^-e.
That average is exactly the joint TV in (1),proving the strong guarantee.
Selected affine variant AGHPaff. Use the SAME matrix formula,but at
field degree R=ceil(log(N+1))+m+e+1 and2R seed bits. Every nonempty
matrix parity has bias delta=(N-1)/2^R<2^(-m-e-1),by the same root proof.
For X uniform on an affine K-dimensional space with direction V,output
is exactly uniform if M restricted to V has rank m. For nonzero a in F_2^m,
Pr[a^T M annihilates V]=2^-K sum_(v in V) E(-1)^(a^T Mv)
<=2^-K+delta. The zero vector contributes1;each nonzero v gives a
nonempty matrix parity. Union over a!=0 bounds failure by
(2^m-1)(2^-K+delta)<2^-e if K>=m+e+1. Affine shifts do not affect
rank;bad seeds cost at most1. This proves strong affine threshold m+e+1,
full2R-bit seed retained,and fixed-seed source-linearity. This AGHPaff is an unselected library map. Ordinary AGHP and W2 below retain their original field
degree and arbitrary-source threshold m+2e;they are different instances.
## 2. A weak seed whose entropy deficit is paid explicitly
**Lemma W2 (the selected weak-seed bound).** For the AGHP matrix
of Section1 or the canonical family7g,assume Hinf(X)>=m+2e and an independent full q-bit
seed S has entropy>=q-Delta; unused prescribed suffixes may be retained.
Then TV((E(X,S),S),(U_m,S))<=2^((Delta-1)/2)*2^-e.
For each seed let f be its conditional TV. Cauchy--Schwarz pointwise gives
f^2<=2^m/4 times the output collision excess. Section1 and7g each bound
its uniform-seed mean by(1/2)2^-2e. Weighting
by the seed law multiplies this mean by at most2^Delta; Jensen bounds the
mean f by its root. The entire seed is retained. This bound uses the supplied collision calculation.
## 4. Deterministic field setup at the actual outer length
For the local application,d,q,r=O(log n),where n is the actual input length
of the desired final unseeded family. Enumerate monic binary polynomials of
degree r in lexicographic coefficient order. Test irreducibility by the
standard finite test: X^(2^r)=X modulo P,and for each prime divisor p of r,
gcd(P,X^(2^(r/p))-X)=1. These conditions exclude every proper factor degree
and guarantee that the quotient is the field with 2^r elements.
Here are the algebraic details to make the setup premise explicit. Every
polynomial over a field has a finite splitting field: factor off an
irreducible polynomial,form the field quotient by that factor to adjoin a
root,and repeat for the remaining factors. The quotient is a field because
Euclid's identity supplies inverses modulo an irreducible polynomial; an
irreducible factor exists by descending degree. Each extension is finite
when the base field is finite. In a splitting field of X^(2^r)-X,its roots
form a subfield: Frobenius preserves sums/products,and inverses of nonzero
roots are roots. Its derivative is one,so there are exactly 2^r distinct
roots. We have thus constructed a field F with exactly 2^r elements.
Frobenius x -> x^2 is a permutation of F whose r-th power is identity.
For a in F its orbit length j divides r. Its minimal polynomial degree is
exactly j: the orbit product product_(i=0)^(j-1)(X-a^(2^i)) has coefficients
fixed by Frobenius,hence in F_2 (the fixed equation x^2=x has just roots
0,1),while any binary polynomial vanishing at a vanishes at every orbit
point. No proper-degree polynomial can vanish at j distinct points. Every
element of degree a proper divisor j of r is a root of X^(2^j)-X. If r>=2,
the union of all these proper-degree sets has size at most
sum_(j=1)^floor(r/2) 2^j = 2^(floor(r/2)+1)-2 < 2^r.
Thus some a has degree r,and its minimal polynomial is irreducible of
that degree. The r=1 case is immediate. This proves the enumeration will
terminate without a probabilistic choice of a per-length circuit.
For an arbitrary irreducible binary polynomial Q of degree j,its quotient
has 2^j elements. Each nonzero element has (2^j-1)-st power one by the
elementary finite-group counting proof of Lagrange's theorem,so Frobenius
to power j is identity. For the generating residue of X,its orbit has
length exactly j by the preceding minimal-polynomial argument. Thus
Q divides X^(2^h)-X exactly when j divides h. The displayed congruence test
makes every factor degree divide r and excludes repeated factors,since
X^(2^r)-X has derivative one. If a factor degree j were a proper divisor
of r,choose a prime p dividing r/j; then j divides r/p,and the relevant
gcd test would detect Q. The tests thus force a single degree-r
irreducible factor,proving their correctness.
There are 2^r candidates,and each test uses polynomially many r-bit
polynomial operations. Trial division factors r in polynomial time in r.
The search is thus 2^r*poly(r)=n^O(1). This argument specifically relies
on r=O(log n). It does NOT make an r=Theta(n) field search polynomial.
The matrix can be streamed using N=d*m field powers and row sums,so no
exponential table of seed values is stored or searched.
## Deterministic field presentation and exact asymptotic bit cost
For any positive integer h_b,use the canonical exact-degree-h_b field of §7a
above,which proves existence and correctness of its lexicographic
modular-squaring/gcd test for every positive degree,without advice.
Schoolbook multiplication and reduction of degree-below-h_b binary polynomials
cost O(h_b^2) bit operations. A modular Frobenius test uses at most h_b
squarings; the gcd tests use polynomially many operations on h_b-bit
polynomials. A deliberately loose upper bound per candidate is O(h_b^5) bit
operations,including factoring h_b by trial division. There are 2^h_b
candidates,hence setup costs at most O(2^h_b h_b^5) bit operations and
polynomial working space. Once P is fixed,field multiplication costs O(h_b^2) bit
operations and O(h_b) temporary storage.
For the intended local use h_b<=C log_2 n with one fixed C,where n is the
actual outer input length. Thus setup is at most n^C times a fixed
polylogarithmic factor and evaluation is polynomial in n. All finite lengths remain
well-defined by the same search,with no uncomputable exceptional branch.
Fixed rational rates in forward split order:
vMI=(2089386889935331/1000000000000000,1044689989485183/500000000000000,2088853616865983/1000000000000000).
vMF=(26117786802641/12500000000000,417882710485773/200000000000000,1044441427848347/500000000000000).
vMG=(1355357759318483/500000000000000).
vMJ=(2668419535131661/1000000000000000).
vNB0=(1422748394509221/50000000000000,8283980301785549/250000000000000,2069136158127037/62500000000000,8234610855322131/250000000000000,4003390277748289/125000000000000,7803456887141389/250000000000000).
vNB2=(1679072555706599/100000000000000,8770010831323477/500000000000000).
vL0=(775407558545041/250000000000000,263992752467657/62500000000000,4223621378833313/1000000000000000,844066643697397/200000000000000,1043775070687987/250000000000000).
vLr=(3968992880314711/1000000000000000,497527544116513/125000000000000,5249934268502773/1000000000000000,327288677273241/62500000000000,5114402132460489/1000000000000000).
vE=(7843501984127/7812500000000).
I/F usePXTM(JI=(3,3,3),vMI/vMF);G usesPXTM(JG=(3,),vMG);J usesbinaryvMJ. Fz usesPXTM(JFZ=(3,4),vMFZ);Gz/Jz usebinaryvMGZ/vMJZ. Each has independent rates and directly paid thresholds.
Ext0/2 useJS=(4,5,5,5,5,6)/JY=(5,6); L0 usesAInit7f-local; Lr usesAInit7f-local; E usesPXTM(JE=(2,),vE).
Lm usesAInit7f-local. B9 checks every7g/7h guard; vectors are fixed before n.

vMFZ=(2286924234837741/1000000000000000,3251779870093499/1000000000000000).

vMGZ=(998945090248823/500000000000000).

vMJZ=(212182348966477/100000000000000).

# Relocated sampler appendices (SP.1,SP.3--SP.4 and SP.7--SP.8)
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
Put nu=2^-80,D_E=15404805937759/905694466927+5nu,
h_R=2719083400781/905694466927+2nu,gamma=(2-nu)/(14+(h_R+5/2)/D_E),
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
