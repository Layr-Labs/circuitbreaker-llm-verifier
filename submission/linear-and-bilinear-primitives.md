# linear-and-bilinear-primitives.md
Fix Q_S=10696922043222096617696249/5000000000000000000000000. With L=floor(log2(actual input length)), S=floor(L/Q_S)
is an integer; S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer rounding.
Component references resolve through the map in sampler-and-disperser.md.
# Finite seeded and linear extraction: an elementary logarithmic-regime replacement
All logarithms are binary; distance means half-L1.
## 1. Selected interfaces
Sections3--4 supply the strong linear condenser, exact block-source repair
and shared hash. Sections7g--7h give complete selected recurrences and
finite seed/entropy envelopes;7a supplies canonical exact fields, and7f
supplies the affine-only base. No unselected extractor theorem is needed.
## 2. Deterministic fields, with no irreducible-polynomial search
Construct nested fields F_i of binary dimension 2^i. Start with F_0=F_2 and
a_0=1. Given F_i and an element a_i of absolute trace 1, adjoin z_i satisfying
z_i^2+z_i+a_i=0,
a_(i+1)=a_i z_i.
Here absolute trace in a field of binary dimension v means
sum_(j=0)^(v-1) x^(2^j). Finite-field powering x^(2^v)=x follows from
Lagrange's counting argument on the nonzero multiplicative group. Relative
trace in a quadratic extension of a field of size Q is x+x^Q; grouping
the Frobenius sum proves trace transitivity.
This really is a quadratic field extension. The map x -> x^2+x is F_2-linear,
has kernel {0,1}, and its image has half the field's elements. Every image
has absolute trace zero, because Tr(x^2)=Tr(x). Trace is nonzero by the
induction hypothesis Tr(a_i)=1; thus its zero kernel also has half the
field's elements. The two sets coincide, and a_i is not an image.
The polynomial has no root and is irreducible. Its two roots are z_i,z_i+1,
so the relative trace of a_i z_i is a_i. Hence a_(i+1) has absolute trace 1.
This proves the induction.
Use ordered tower coordinates, with constants embedded by appending zero
coordinates. Addition is XOR. For x=x0+x1 z_i and y=y0+y1 z_i compute
p=x0*y0, q=x1*y1, r=(x0+x1)*(y0+y1), t=a_i*q in the preceding field.
The product is (p+t)+(r+p)z_i, because z_i^2=z_i+a_i. These four recursive
products give T(v)<=4T(v/2)+O(v), hence O(v^2) bit operations. Inverses can be computed
by exponentiation to 2^v-2 in a field of binary dimension v, using O(v)
squarings/multiplications. Zero inversion is never requested.
We also require a polynomial-basis representation of an extension of degree
l, where both l and the base-field binary dimension s are powers of 2 and
l>=2. Build the tower to binary dimension sl. Let g be the last adjoined
root. It is not in the immediately preceding subfield. Its orbit under
q-th powering, q=2^s, has length exactly l. Indeed the orbit length divides
l; a proper divisor of the power of 2 l is at most l/2. Then g would satisfy
x^(q^(l/2))=x. The preceding subfield has exactly q^(l/2) elements and fills
the roots of that polynomial, contradicting the choice of g.
Consequently
E(T)=product_(j=0)^(l-1) (T-g^(q^j))
has coefficients fixed by q-th powering, whose roots are exactly F_q.
Every F_q-polynomial vanishing at g also vanishes at all l orbit points,
so no nonzero polynomial of degree below l can vanish there. Thus E is
irreducible of degree l over F_q. Compute that product in the
tower and read its coefficients in the embedded base field. This constructs
E in polynomial time in sl and gives a deterministic canonical choice.
No factoring, primitive-element search, or enumeration of 2^(sl) field
elements is involved. All choices, including coordinate orders, are fixed.
## 3. A finite strong linear lossless condenser
This section fixes integers U>=1, 1<=k<=U, a>=1 and R>=1. Let
rho=a+2 ceil(log U)+2,
s=the least power of 2 at least (R+1)rho,
u=s-rho, q=2^s, H=2^u,
l=max(2, the least power of 2 at least ceil(U/s)),
b=ceil(k/u), v=bs.
Use Section 2 for F_q and an irreducible E(T) of degree l. Pad the U-bit
input with zeros to sl bits and view it as the coefficient vector of
f(T), degree <l, over F_q. A seed is y in F_q. Output the b field elements
Con(f,y)=(f(y), (f^H mod E)(y), ...,
(f^(H^(b-1)) mod E)(y)). (4)
The output has v bits. Each map f -> f^(H^j) mod E is F_2-linear, because
H is a power of 2. Evaluation and padding are linear too. Moreover
s < 2(R+1)rho, u>=Rs/(R+1),
v <= (1+1/R)k+s. (5)
We prove that Con is a strong lossless condenser: for every source X of
min-entropy k, (Y,Con(X,Y)) is within 2^(-a) of a distribution (Y,Z) with
Y uniform and H_infty(Z | Y=y)>=k for every y.
### Expansion proof
Write A=q-(l-1)(H-1)b. The graph with neighbors (y,Con(f,y)) expands every
set of J<=H^b inputs to at least AJ neighbors. To prove this, suppose a
set T of fewer than AJ neighbors contains every neighbor of some inputs.
For j=0,...,J-1, write j in base H with b digits, and let M_j(Z_0,...,Z_(b-1))
be the corresponding monomial. There are AJ coefficients in a polynomial
Q(Y,Z_0,...,Z_(b-1))=sum_(j<J) p_j(Y) M_j(Z_0,...,Z_(b-1)),
deg p_j < A.
The fewer than AJ conditions of vanishing on T admit a nonzero solution.
Choose one with least degree in Y. Not all p_j are divisible by E: otherwise
dividing by E preserves vanishing on T (E has no F_q root because l>=2) and
lowers that degree.
For any input f whose neighbors all lie in T, substitute the polynomials
f_j=f^(H^j) mod E. The resulting polynomial in Y has degree at most
(A-1)+(l-1)(H-1)b=q-1 and vanishes at every y in F_q. It is zero. Reducing
modulo E therefore says that f, as an element of F_q[T]/E, is a root of
sum_(j<J) (p_j mod E) Z^j.
This is a nonzero polynomial of degree <=J-1, so there are at most J-1 such
inputs. In particular a J-input set cannot have fewer than AJ neighbors.
The same statement holds after restricting the source to the padded inputs.
Here l<=2U and b<=k<=U. Also q/H=2^rho>=4 U^2 2^a. Therefore
((l-1)(H-1)b)/q <= 2U^2/2^rho <= 2^(-a-1).
Thus A>0 and A>=(1-2^(-a))q. Since H^b>=2^k, the required size-2^k
sets have at least (1-2^(-a))2^k q distinct neighbors.
### Strongness, not only joint min-entropy
For a source uniform on a size-K set, K=2^k, consider each y separately.
Keep one preimage for each distinct value Con(x,y). To obtain a K-element
set of output values at that y, replace the duplicates by unused values.
There is room, because v>=k (s>=u and b>=k/u). The number replaced,
summed over all y, is at most 2^(-a)Kq by expansion. The repaired conditional
output at each y is uniform on K values, and the seed marginal is unchanged.
The joint distance is at most 2^(-a). An arbitrary k-source is a convex
combination of uniform size-K sources: iteratively express a probability
vector with coordinates <=1/K as a mixture of vertices of that capped
probability simplex; its vertices have exactly K entries 1/K. Mixing the
repairs keeps each conditional atom <=1/K and preserves the distance bound.
This proves the asserted strong condenser, including the seed marginal.
## 4. Three exact block operations
A block source (X_1,...,X_j) has block entropy k if every conditional law
X_i | X_<i=x_<i has min-entropy >=k on its support. All statements below
also hold conditioned on any fixed external variable.
**Apply a strong condenser with one shared seed.** If Con has error eta,
apply Con(X_i,Y) to all j blocks with the same uniform independent Y.
The output tuple, retaining Y, is within j eta of a block source with the
condenser's guaranteed conditional entropy. A fully safe coupling argument
is as follows. Retain the original source blocks as latent variables. For
each original prefix and seed, choose the condenser's repaired conditional
output law and a coupling to the actual output. Generate these repairs
causally, with independent fresh coupling randomness. Given an original
prefix and the seed, previous repairs carry no additional information about
the next original block. Each repaired next block has the desired entropy;
mixing over latent prefixes preserves the atom bound. Each repair disagrees
with its original output with probability at most eta on average, because
the original prefix is independent of the uniform seed. A union bound proves
the distance. No coupling is part of the evaluator.
**Split one high-entropy block.** Pad a v-bit k-source by one zero if needed
and split into two m-bit halves, m=ceil(v/2). If k>=m+k'+a, the pair is within
2^(-a) of a block source with entropy k' in each half. The first half already
has entropy >=k-m>=k'+a. The average maximum probability of the second half
given the first is at most 2^(m-k). Markov bounds by 2^(-a) the mass of first
halves whose conditional second-half entropy is below k'. Replace the
second-half law by uniform only on those prefixes. The first-half law is
unchanged. Applied to j input blocks with the causal latent-prefix coupling
above, the total distance is at most j2^(-a). This proof also explains why
arbitrary conditioning on a rare prefix is not silently allowed.
**Hash every block with one shared seed.** For m-bit blocks, put W equal to
the least power of 2 at least m. Use the tower field of binary dimension W.
Pad each block to W bits. For a uniform z in that field define Hash_z(x)
as the first ell coordinates of zx. For distinct x,x', z(x-x') is uniform,
so collision probability is exactly 2^(-ell).
For a single source of entropy >=ell+2a, its conditional output collision
probability, averaged over z, is at most 2^(-ell)+2^(-ell-2a). Cauchy-Schwarz
and Jensen give strong TV distance <=(1/2)2^(-a). For a j-block source this
yields joint jell output bits with error <=(j/2)2^(-a), retaining z. Replace
outputs from the LAST block backwards: the extractor bound can retain the
entire previous source prefix, which is independent of z. Project that
prefix to its prior hashed outputs, then repeat. This avoids conditioning
on already hashed outputs as though the shared seed were still independent.
This hash is linear for each z and its seed length is exactly W<=2m.
## 5. Arithmetic cost and selected outer runtime
For the exact selected recurrences, m_i>K_i implies K_(i-1)>j_i*K_i,
so B<=K_0<=u. Also b_i>=1 and K_(i-1)>((j_i-1)/j_i)S_i give
S_i<2u; every later input and leaf hash has dimension<=u. These guards
precede field setup. For the Section3 condenser, with U=4u, the least
relative degree l=max(2,pow2ceil(ceil(U/S_i))) satisfies
S_i*l<=2U+2S_i=O(u). The orbit polynomial uses
O(l*S_i) squarings and O(l^2) products in that tower. Its coefficients
occupy the initial base-field coordinates, so no expensive basis search
is required. Polynomial multiplication/reduction, repeated squaring and
Horner evaluation cost polynomially many bit operations once the base
fields are constructed. There are at most B nodes per level and a fixed
number of levels at every selected call.
Exact base/hash field searches need 2^O(v)poly(v) time at degree v;
no polynomial bound in arbitrary component input length is asserted.
At the actual outer calls these degrees are O(log actual N), by the
literal envelopes and closed slots. Their searches therefore cost
N^O(1). Larger relative extensions use the explicit polynomial-cost
towers, including the affine base's initial O(N)-dimensional extension.
All other arithmetic, matrices, kernels and output evaluation are
polynomial in actual N. Sections7a/7f and proof.md supply the complete
all-length guards and the fixed total evaluator's exceptional branches.
## 7. Exact multiway block-source repair
For a positive integer j, pad a v-bit k-source to jm bits, m=ceil(v/j),
and split into Z_1,...,Z_j. If k>=(j-1)m+k'+a, the tuple is within
(j-1)2^-a of a block source with entropy k' per block. Indeed the prefix
of length i has entropy>=k-(j-i)m. Averaging its largest next-block atom
across earlier prefixes gives at most2^((j-1)m-k). Markov bounds prefixes
of conditional entropy<k' by2^-a. Generate a repaired tuple sequentially:
use the original conditional on good prefixes, otherwise uniform on m bits.
Its conditional entropy is>=k', since k<=jm implies k'<=m; the union bound
over the j-1 possibly bad steps gives the error. The first block is already
good since its marginal entropy is>=k-(j-1)m>=k'+a. Conditional external
variables and causal repairs obey Section4's same argument. This lemma
is used at each selected j_i by7h, and at j=2 by7g.
## 7a. Exact binary field degrees for the base and both NB primitives
Condenser/AGHP fields use exact binary degrees.
For each positive degree v
choose the first monic irreducible binary polynomial by lexicographic
enumeration and the modular-squaring/gcd test proved in weak-seed Section4
later in this file. That proof supplies existence and correctness at every
v. Quotient fields use polynomial coordinates and schoolbook arithmetic.
The trace is a nonzero binary linear map: sum X^(2^j),j=0,...,v-1,
is a nonzero polynomial of degree2^(v-1)<2^v, so cannot vanish everywhere.
Test the ordered basis units and choose the first of trace1. Starting
from that base field, use Section2's Artin--Schreier tower. Its trace
transitivity and kernel/image proof give irreducibility and propagate
trace1 at every quadratic step. Only the relative extension degree is a
power of two. The last root's full Frobenius orbit and its orbit product
give the required condenser modulus. Larger extensions are therefore
built by towers, not enumeration in their full dimension. Section4's optional field hash
uses canonical degree W and the same collision proof.
Section3's expansion requires only the stated field arithmetic and degree
inequalities; Sections7g/7h give each selected exact recurrence, padding,
error and complete-seed proof directly.
F#/NB inputs/outputs/error and selected output/error/base/leaf degrees
are O(log actual N) by B9; outer L0/Lr inputs can be N bits.
BASE7f has b0,b=O(log N), U1=2^(110+lg S)=O(log N); towers build its initial
O(N)-dimensional relative extension. BASE's later ambient/matrix/kernel dimensions
are O(log N). Canonical searches cost2^O(log N); tower/matrix operations
are polynomial in actual N, without field or source advice. No joint bound
for arbitrary component dimensions is asserted.
## 7f. A three-seed affine base from linear condenser duality
Only this BASE is specialized to uniform affine sources. All PXT, NB and
internal breaker source/right-tuple contracts remain general-source.
An affine k-source is uniform on a translate of a binary k-dimensional
vector subspace. Write Con[U,k,rho,b] for Section3's actual polynomial
map at integer field degree b>rho: Hpow=2^(b-rho), relative degree
l=max(2,pow2ceil(ceil(U/b))), output length b*ceil(k/(b-rho)). Pad shorter
inputs by zeros. Section7a supplies the canonical exact-degree field and
relative power-of-two tower. The polynomial maps f -> f^(Hpow^i) mod E
are binary linear, as are padding and evaluation. Thus this is source-linear.
Section3's expansion and seed-preserving repair proof uses these integers
and the defect (l-1)(Hpow-1)*ceil(k/(b-rho))/2^b, not an integer-rate
hypothesis. For k<=U, rho=a+2lg(U)+2 and b-rho>=1, l<=2U and
ceil(k/(b-rho))<=k<=U bound the defect by 2^(-a-1). The same proof
therefore gives strong lossless condensation with error <=2^-a.
Rank preservation. For an affine source of dimension k, a source-linear
condenser's output at seed y is uniform on a translate of a rank-r_y
subspace, r_y<=k. If r_y<k, every repaired conditional law with atom cap
2^-k gives this support mass <=2^(r_y-k)<=1/2. Its conditional TV distance
is therefore >=1/2. Since the repair retains the uniform seed marginal,
joint TV is the average of these distances. At most a 2eta fraction of
seeds can fail injectivity on any fixed k-dimensional direction space
when the strong condenser error is eta.
Kernel duality. Let V<=F_2^v have dimension k and W=V-perp under the
binary coordinate dot product. Suppose G:F_2^v -> F_2^w is injective
on W. Choose any d independent vectors from ker G as matrix rows M.
If M restricted to V did not have rank d, a nonzero row combination
would annihilate V and hence belong to W intersect ker G, a contradiction.
A rank-d linear image of a uniform affine V-source is exactly uniform
on d bits, since all fibers are cosets of the same kernel. This argument
allows G to have dependent rows; it assumes no surjectivity of G.
For actual n>=1 fix Q_S=10696922043222096617696249/5000000000000000000000000, beta=50/7-2^-100, delta_O=760128341,
eps_B=10^-18, eps_R=2^-80, N0=2^60, a0=1+2^-60,
The guard S>=2^80 pays each unit rounding by eps_R*S>=1.
gamma_B=(50/7-a0*(Q_S+1+eps_B)-3eps_R)/(2(1+eps_B)),
xi_B=gamma_B/(gamma_B-1), G_B=1-(xi_B-1)^2,
A_B=1/G_B, B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)=1155.
First compute L=floor(log2 n), S=floor(L/Q_S), d=0.
If S<2^160 return 0^d BEFORE evaluating any lg S. Otherwise set
these ACTUAL integer parameters in the following order:
ell_w=lg(S+4), s=2ceil((beta*S-ell_w)/2)<=8S,
R(s)=0,
ellc=lg(s+4), r_F=S, a_B=r_F+3,
U=4n, U1=2^(110+lg S), rho=a_B+2lg(U1)+2=r_F+225+2lg S,
Om=9s+307S+2L+84, d=floor((1859462242691845490204131/500000000000000000)Om), K_B=ceil(A_B*d+B_B*rho),
rho0=a_B+lg U+lg K_B+2,
b0=ceil((N0+1)rho0/N0), u0=b0-rho0,
b=ceil(gamma_B*rho), u=b-rho,
v0=b0*ceil(K_B/u0), v1=b*ceil(K_B/u),
c=v1-K_B, v2=b*ceil(c/u), D=b0+2b. (AF1)
O2a gives d<=delta_O*L; all calibrations below use this upper bound.
First require S>=2^160. Test K_B<=n, rho0>=N0, u0>=1, u>=1,
v0<=U1, v1<=U1, 1<=c<=U1, v1-v2>=d and D<=s. A failed
public dimension guard returns 0^d, before field construction.
Parse the first D seed bits as independent uniform y0,y1,y2 of lengths
b0,b,b. Compute z0=Con[U,K_B,rho0,b0](x,y0), then
z1=Con[U1,K_B,rho,b](z0,y1). Define G_y2 by restricting
Con[U1,c,rho,b] to its first v1 input coordinates through zero padding.
Form its v2-by-v1 binary matrix by evaluating it on the ordered v1
coordinate units. Row-reduce it over F_2, processing columns increasingly,
choosing the first remaining row with a pivot1, and eliminating that column
from every other row. For each free column in increasing order, set that
free coordinate1, the other free coordinates0, and solve pivot coordinates.
These vectors are the canonical ordered kernel basis. Use the first d
vectors as rows M_y2 and return M_y2*z1. Ignore the final s-D seed bits.
Only the PUBLIC matrix kernel is computed.
Initial condensation is legal even at its smaller degree. Since
u0>=rho0/N0>=1, ceil(K_B/u0)<=K_B, and b0/u0<=N0+1,
its output has v0<=(N0+1)K_B+b0 bits. Its defect is bounded by
2U*K_B/2^rho0<=2^(-a_B-1).
The two remaining condensers have thresholds K_B and c, both <=U1,
ambient U1 and rho=a_B+2lg(U1)+2, so their defects have the same bound.
Dimension, including BOTH ceilings. The exact constants satisfy
gamma_B=161549143293335727659155991612931407279525089/80704505322479288400704505322479288320000000>2, xi_B=161549143293335727659155991612931407279525089/80844637970856439258451486290452118959525089,
G_B=22638309292754974231255878057592598165219246943507051605987025368933246665400420457921/6535855488638842764486344832123456581535745617036971673634641620742621669065400420457921>0, G_B*A_B=1,
G_B*B_B-xi_B*gamma_B>=xi_B/2^40.
Because b>=gamma_B*rho, b/(b-rho)<=xi_B. Consequently
v1<=xi_B*K_B+b, c<=(xi_B-1)K_B+b, v2<=xi_B*c+b.
Using the EXACT identity v1=K_B+c, not a false lower bound on v1, gives
v1-v2>=K_B-(xi_B-1)c-b>=G_B*K_B-xi_B*b
>=G_B*A_B*d+(G_B*B_B-xi_B*gamma_B)rho-xi_B
>=d+xi_B*(rho/2^40-1)>=d. (AF3)
Here ceil(K_B/u)<=K_B/u+1 and b<=gamma_B*rho+1 pay the two
roundings; the outer K_B ceiling helps this lower bound. Also c>0 since
b>u and v1=b*ceil(K_B/u)>K_B. Thus for EVERY y2, rank(G_y2)<=v2
and dim ker G_y2>=v1-v2>=d. The evaluated rows always exist.
For an affine K_B-source, the first condenser preserves its direction
dimension except on fraction <=2*2^-a_B of y0. For every good y0, the
second preserves that dimension except on fraction <=2*2^-a_B of fresh
y1. For a good pair let V be its K_B-dimensional direction space in
F_2^v1. Its perpendicular space has dimension c. Fresh y2 makes G_y2
injective on that space except on fraction <=2*2^-a_B. Kernel duality
then makes the actual output exactly uniform on every good seed triple.
A union bound gives TV <=6*2^-a_B=(3/4)2^-r_F<2^-r_F, retaining the entire seed
triple. The independent ignored suffix can also be retained. For a larger
affine source choose any K_B-dimensional direction subspace only in this
analysis: surjectivity on it implies surjectivity on the full source.
Complete all-length calibration. For 2^k<=S<2^(k+1), k>=80,
s+4<9S implies ellc<=k+5. Therefore
rho-S<=10k+267<eps_B*S.
The exact coefficient bound is
K_B/S<=A_B*delta_O*Q_S*(1+2^-40)+B_B*(1+eps_B)+eps_R
=3829375509427054392217624301414698808695802031482475436008044830040930066628828606859775073632253843274862847262146780155220951964920623333/8156310503611947103200570487561190992933448972625625424013960560814326081702769562918430317875591905280000000000000000000000000<2^48.
Here L<Q_S(S+1), rho<=(1+eps_B)S, and eps_R*S>=1 pay K_B's ceiling.
Thus lg K_B<=48+lg S. Since lg(4n)<=L+3,
rho0-(Q_S+1)S<=9k+Q_S+97<eps_B*S.
The k80 bases are 1067<10^-18*2^80 and 817+Q_S<10^-18*2^80. Their affine left sides
increase by10 and9, at most their bases, while eps_B*2^k doubles.
This proves both bounds for EVERY integer S>=2^80. Each field-degree
ceiling costs <=1<=eps_R*S. By the definition of gamma_B,
D/S<=a0*(Q_S+1+eps_B)+2gamma_B*(1+eps_B)+3eps_R=50/7. (AF2)
For the reduced-width variant retain this gamma_B, which uses the budget 50/7.
Before charging the three rounding units by eps_R*S, the same calculation gives
D<=(50/7-3eps_R)S+3. For S>=2^160,
put M=3eps_R-2^-100>0. On 2^j<=S<2^(j+1), j>=160,
ell_w<=j+2. The exact base M*2^160>=165 and doubling prove
M*S>=ell_w+3, since the right dyadic bound j+5 increases by1.
Thus D<=(50/7-2^-100)S-ell_w<=s.
Also (beta-7)*2^160>=162; the same doubling argument proves
(beta-7)S>=ell_w, hence s>=7S. All previous upper bounds s<=8S remain valid.
Thus the COMPLETE three-slice seed and ignored suffix fit the new s.
Put b0coef=a0*(Q_S+1+eps_B)+eps_R. The exact inequalities
(2^60+1)2^48+b0coef<2^110,
xi_B*2^48+gamma_B*(1+eps_B)+eps_B<2^110
imply v0,v1<2^110*S<=U1, and 0<c<v1<=U1. Also
K_B<=[A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R]L
<=H_base*L, H_base=219454948582. (AF4)
The last allowance pays the K_B ceiling. H*L<=n and H>=H_base give
K_B<=n. Finally rho0>=(Q_S+1)S>2S>=2^60=N0 and rho>=S>=2^40 prove every
public AF1 dimension guard and the sharpened final inequality in AF3.
Each enumerated exact field degree b0,b is O(log actual n). The initial
relative extension has dimension O(n) and uses the explicit towers of
Sections2/7a, never enumeration in that large dimension. U1, v0, v1,
v2 and the evaluated public matrix dimensions are O(log n), with enormous
fixed coefficients. Matrix evaluation and binary Gaussian elimination are
polynomial in these dimensions. Thus this is one deterministic polynomial
evaluator in actual n, source-linear at every full seed.
The classical condenser-duality principle is credited to Mahdi Cheraghchi
Bashi Astaneh, Applications of Derandomization Theory in Coding, EPFL thesis
4767 (July2010), Corollary5.11 and following discussion,p.148,
https://eccc.weizmann.ac.il/resources/pdf/astaneh.pdf . The full rank and
kernel proofs needed here were supplied above; the citation is attribution.
## 7g. Exact finite rational-rate primitives PXT
Fix h>=1 and a vector of fixed rational rates R_i>1, i=1,...,h.
PXT[h,Ro,Rl] abbreviates R_1=...=R_(h-1)=Ro,R_h=Rl;
PXT[h,R] means every rate is R. For integers u,z,e>=1, put
B=2^h,ell=ceil(z/B),a=e+h+4,U=U1=4u,
rho=a+2lg(4u)+2,S_i=ceil((R_i+1)rho),
lambda_i=2R_i/(R_i-1),K_h=ell+2a,
b_i=ceil(2(K_i+a+1)/(S_i-2rho)),m_i=ceil(b_i*S_i/2),
K_(i-1)=m_i+K_i+a,i=h,...,1. (PT0)
Check K_0<=u before field setup; otherwise return0^z. Start directly with
the original u-bit source. Every condenser uses Section3 with common
ambient bound U=4u, padding its current input with zeros. Section7a supplies
exact-degree fields and relative power-of-two extensions. At level i use
degree S_i and b_i field outputs for a binary condenser/split, one shared
fresh seed across that level. Hash all B
leaves with the canonical linear family below, using one fresh shared seed, then truncate
to z. Retain every seed and every padded suffix. This defines the map.
Since R_i>1, S_i>2rho. The block-count choice gives
b_i*(S_i/2-rho)>=K_i+a+1,
K_(i-1)<b_i*S_i/2+1+K_i+a<=b_i*(S_i-rho),
b_i<=m_i<K_(i-1),m_i>K_i. (PT0a)
Section3's expansion proof allows any b_i<=K_(i-1) with
b_i*(S_i-rho)>=K_(i-1): these imply its b<=U and H^b>=2^k
conditions; 2^rho>=4U^2*2^a is unchanged. The split condition is equality
in PT0. Thus the causal repairs of Section4 apply. The first input has
length u<=U; every later input has length m_i<K_(i-1)<=K_0<=u.
This proves all padding guards without an initial condenser. Such a stage
with seed 2rho was redundant because the ambient bound never shrinks here.
The condenser/split errors sum to 2(B-1)2^-a; the PXT hash costs
less than B2^-a, giving less than (3B-2)2^-a<3B2^-a<2^-e. Strongness and
source-linearity follow as in Section4, retaining every actual seed.
Here are sharper FINITE envelopes. For integer
w>=max(2^20,e+h+6+2lg(4u)),rho<=w,a<=w-6.
Define A_h=2,B_h=1/B and, backwards,
A_(i-1)=lambda_i*(A_i+1)+(R_i+1)/2,B_(i-1)=lambda_i*B_i.
Put Lam=B_0=product_i R_i/(R_i-1),g=A_0. Then
K_0<=Lam*z+g*w. (PT1)
This is the j_i=2 specialization of the PTM.1 ceiling induction below;
its residual e_i=0, with the same K_h<=z/B+2w-11 initial bound.
**Canonical general-source linear leaf family.** For positive integers N,m,a set
q=2^-m, eps=2^-2a, b=1+eps/2, d=m+4a+lg(4(N+1)), T=2^d.
Construct T binary m-by-N matrices H_s with
Z_x=|{s:H_s x=0}|<Tq(1+eps) for every nonzero x.
Initially take independent uniform matrices and potential F=sum_(x!=0)b^Z_x.
Its expectation is (2^N-1)(1+q*eps/2)^T. For t=eps/2,
ln(1+t)>=t-t^2/2 and ln(1+qt)<=qt give
q(1+eps)ln b-ln(1+qt)>=q*eps^2/4.
Since Tq*eps^2>=4(N+1), the initial expectation is strictly below
b^(Tq(1+eps)): indeed ln(2^N-1)<N ln2<N+1.
Fix entries in seed/row/column order by minimum conditional expectation,
ties choosing0. This is an exact canonical finite algorithm. For each x,
a partially fixed row is either a determined parity or uniform. Thus
p_x=Pr[H_s x=0|prefix] is0 if a determined row is1, otherwise2^-r,
where r counts undetermined rows. Future matrices have the common factor
1+q*eps/2. Set Den=2^(2a+1), Num=Den+1, and maintain integer
W_x=Num^(completed zero count)*Den^(completed nonzero count).
Compare the trial-bit integers sum_x W_x*(Den*2^m+[p_x>0]*2^(m-r)).
Omitted denominators and future factors are common and positive.
At matrix completion multiply W_x by Num or Den according as H_s x=0
or not. Conditional expectation never increases, proving the claimed Z_x.
Weights have O(Ta)bits. The TmN bit choices each enumerate2^N-1 vectors
and inspect mN entries; all arithmetic and table generation have time
polynomial in2^N,T,N,m,a. There is no source enumeration or advice.
Linearity gives collision at most q+q*eps for every distinct pair.
For any source of minentropy k>=m+2a, the mean output collision is at most
2^-k+q+q*eps. Cauchy--Schwarz bounds retained-seed TV from uniform by
.5 sqrt(2^m*(2^-k+q*eps))<=2^-a/sqrt(2)<2^-a.
The same calculation gives mean_seed TV^2<=(1/2)2^-2a.
This is strong extraction for every source, with a source-linear map for
fixed seeds. The table depends only on dimensions.

Apply this family to N=m_h,m=ell with a FRESH shared d-bit seed.
The causal repairs precede this seed Y and give each leaf entropy ell+2a
conditional on the retained condenser seeds and its full repaired prefix.
Reverse-hybrid replace leaves B,...,1, retaining Y, earlier source leaves
and independent later uniform outputs. Apply the single-block average-Y
bound, then project earlier source leaves to their hashes. Each step costs
<2^-a even with Y retained; unused independent suffixes remain retained.
Together with 2(B-1)2^-a this gives total<(3B-2)2^-a<2^-e.
PXG is now this SAME terminal variant as PXT; all thresholds stay PT0/PT1.

**Unified legal seed envelope (PT2).** Require h<=13, z<=256B*w,
terminal arity2<=j<=18, R_h<=64 (binary:16/15<=Rh<=64) and kap=(R_h+1)/(R_h-j+1)<=31,
in addition to K0<=u and w>=max(2^20,e+lg(B)+6+2lg(4u)).
For binary j=2. The terminal ceiling calculation gives
m_h<kap*(ell+3a+j-1)+S_h/j+1
<=kap*z/B+[3kap+(R_h+1)/j]w+1+1/j+(j-18)kap.
Since j<=18 and R_h<=64, this is <=8192w under the output cap.
Put J=lg(4(m_h+1))<=16+lgw. The actual seed length is
D=sum_iS_i+ell+4a+J<=z/B+sum_i(R_i+1)w+4a+J+h+1.
If u<w/2, PT0a gives K0>rho>=a, so a<w/2 and J+h+1<2w.
Otherwise lg(4u)>=lgw+1, whence w-a>=2lgw+4.
Then4(w-a)-(J+h+1)>=7lgw-h-1>0 since lgw>=20,h<=13.
Both cases prove
D<=c*z+pi*w, c=1/B, pi=4+sum_i(R_i+1). (PT2)
Ignored padding bits remain retained. Outside these extra cap/rate guards
the hash theorem and literal algorithm still hold, but PT2 is not asserted.
Every selected actual call satisfies them, as checked below and in B9.
At every selected call m_h,ell,a=O(log actual n). Even when u=actual n,
the terminal bound is8192w with w=O(log n). Thus2^N,T and the integer
weight bitlength are fixed polynomials in actual n. Construction is uniform
polynomial in actual n, without a joint-polynomial claim in unrelated
component dimensions. Condenser fields are still canonical fields of7a.
## 7h. Mixed-arity exact primitive PXTM
Both Ext0 and Ext2 use this variant; I/F also use JI=(3,3,3); G uses JG=(3,); J remains binary7g; outer L0/Lr/E use JL0=(3,4,4,4,4),JLR=(4,4,4,5),JE=(2,3,3).
Fix h>=1, integers j_i>=2, and rational R_i>j_i-1.
For integers u,z,e>=1 set B=product_i j_i,ell=ceil(z/B),
a=e+lg(B)+4,U=4u,rho=a+2lg(4u)+2,
S_i=ceil((R_i+1)rho),K_h=ell+2a,
b_i=ceil(j_i(K_i+a+j_i-1)/(S_i-j_i*rho)),m_i=ceil(b_i*S_i/j_i),
K_(i-1)=(j_i-1)m_i+K_i+a, i=h,...,1. (PTM.0)
If K_0>u return0^z before field setup. Otherwise begin on the original
u-bit source. Level i applies Section3's condenser with ambient4u,
threshold K_(i-1), degree S_i and b_i outputs, using ONE fresh S_i-bit seed
shared across that level. Pad each result to j_i*m_i bits and split into
j_i blocks. Hash all B leaves to ell bits with the canonical family of7g,
using one fresh shared seed of length ell+4a+lg(4(m_h+1)), then truncate to z. Retain every seed and ignored suffix.
All fields are canonical exact fields/towers from7a. The leaf count B need
not be a power of two; only relative field extensions use quadratic towers.
Here S_i>j_i*rho and b_i(S_i/j_i-rho)>=K_i+a+j_i-1, so
K_(i-1)<(j_i-1)b_i*S_i/j_i+j_i-1+K_i+a<=b_i(S_i-rho),
b_i<=m_i<K_(i-1), K_(i-1)=(j_i-1)m_i+K_i+a.
Section3's expansion guards hold. Each later input has length
m_i<K_(i-1)<=K_0<=u, so all padding is legal. Section7's exact multiway
repair applies. With M_i=product_(l<i)j_l, the node counts satisfy
sum_i M_i(j_i-1)=B-1 and sum_i M_i<=B-1.
Condenser plus split errors cost at most2(B-1)2^-a, and the shared hash
costs less than B2^-a. Thus total<3B2^-a<2^-e, since B<=2^lg(B).
The causal repair proof of4 retains the complete seed and every latent
prefix. Every operation is source-linear for fixed seeds. This proves
strong linear extraction for arbitrary K_0-sources, not just affine sources.
For integer w>=max(2^20,e+lg(B)+6+2lg(4u)), rho<=w,a<=w-6.
Put lambda_i=j_i*R_i/(R_i-j_i+1),A_h=2,B_h=1/B and backwards
e_i=max(0,(j_i-7)lambda_i+(j_i-1)/j_i),
A_(i-1)=lambda_i(A_i+1)+((j_i-1)/j_i)(R_i+1)+e_i/2^20,
B_(i-1)=lambda_i*B_i. Define Lam=B_0,g=A_0.
The actual two ceilings give
K_(i-1)<[j_i(S_i-rho)/(S_i-j_i*rho)](K_i+a+j_i-1)
+((j_i-1)/j_i)S_i
<=lambda_i(K_i+a+j_i-1)+((j_i-1)/j_i)S_i.
Since a+j_i-1<=w+j_i-7 and S_i<=((R_i+1)w+1), the residual
(j_i-7)lambda_i+(j_i-1)/j_i is at most e_i. Since w>=2^20,
this contributes at most (e_i/2^20)w, exactly the added recurrence term.
For 2<=j_i<=6, e_i=0 and the old envelope is recovered.
With K_h<=z/B+2w-11, induction gives
K_0<=Lam*z+g*w, Lam=product_i R_i/(R_i-j_i+1). (PTM.1)
Let j=j_h and kap=(R_h+1)/(R_h-j+1). The terminal block obeys
m_h<kap*(ell+3a+j-1)+S_h/j+1
<=kap*z/B+[3kap+(R_h+1)/j]w+1+1/j+(j-18)kap.
Under the additional legal seed-envelope guards of7g, its proof applies
verbatim to these ceilings, giving c=1/B,pi=4+sum_i(R_i+1) and
D<=c*z+pi*w. (PTM.2)
The source envelope PTM.1 is unchanged. Generic NB tests BOTH output caps
and all the terminal rate/arity/height guards before claiming PTM.2.
Ext0 selects JS=(4,5,5,5,5,6), h=6,B=15000,lg(B)=14,j=6.
Ext2 selects JY=(5,6), h=2,B=30,lg(B)=5,j=6.
Their complete ordered rate vectors are supplied separately.
The actual unit w_N>=Ee+44+2lg(4n) satisfies the above unit guard.
Output6r and all threshold/input/output/retained-seed requirements are
charged by B9 and wrappersI1--I3. All field degrees and inputs are
O(log actual N); the fixed15000 and30 leaf counts and canonical searches/towers
have polynomial bit time in actual N.

## 10. Attribution and source convention
The PV-code expansion/block method is due to Guruswami, Umans and Vadhan:
https://people.seas.harvard.edu/~salil/research/PV-condenser-June08.pdf .
The recursive linear condensation/shared-seed strategy follows
Chattopadhyay, Goodman and Liao, ECCC TR21-075, Lemmas 4.7--4.10 and
Theorem 5.6: https://eccc.weizmann.ac.il/report/2021/075/download/ .
Chattopadhyay--Liao use a sharper variable-parameter form in Lemma 3.13:
https://arxiv.org/pdf/2110.12652v1 . Xin Li's later application is
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/ .
These finite constructions are proved here. Strong output excludes retained seed bits.
# An elementary weak-seed replacement for the local Raz call
This local interface uses the small-bias construction of Alon, Goldreich,
Hastad and Peralta (1992) and the leftover-hash collision argument.
Sources inspected:
https://www.cs.tau.ac.il/~nogaa/PDFS/aghp4.pdf
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/
(Algorithm 3 and Lemma 5.1).
All logarithms are base two, TV is half the L1 norm, and vector/matrix
arithmetic is over F_2. This proof is finite and elementary. A finite-field existence argument,
deterministic field algorithm and outer-length complexity are supplied in
Section 4.
## 1. A finite strong linear seeded extractor
Fix input length d>=1, output length m>=1 and integer e>=1. Put
N=d*m,
r=ceil(log(N+1))+m+2e.
Fix one presentation of F_{2^r} as an r-dimensional binary vector space.
For a uniform seed (u,v) in F_{2^r}^2, define N bits
b_i=<bin(u^i),bin(v)> for i=0,...,N-1,
using the ordinary nondegenerate binary coordinate dot product. Place these
bits, in fixed row-major order, into an m-by-d binary matrix M(u,v). Define
E(x,(u,v))=M(u,v)*x.
For every fixed seed this map is binary linear in x. If X is independent of
the uniform seed and H_infinity(X)>=m+2e, then
TV((E(X,U_{2r}),U_{2r}),(U_m,U_{2r})) <= 2^-e. (1)
Proof of small bias. A nonempty parity of the b_i equals the dot product of
bin(P(u)) with bin(v), where P is a nonzero polynomial with binary
coefficients and degree at most N-1. If P(u) is nonzero, the dot product is
balanced over uniform v; if P(u)=0, the sign equals one. A nonzero polynomial
has at most its degree many roots, by repeated division by X-a. Thus the
absolute bias is at most
delta=(N-1)/2^r < 2^(-m-2e).
For every nonzero w in F_2^d and every nonzero a in F_2^m, the parity
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
That average is exactly the joint TV in (1), proving the strong guarantee.
## 2. A weak seed whose entropy deficit is paid explicitly
Suppose a seeded map uses q seed bits, is strong with error eta for uniform
seed, and ignores any prescribed unused seed suffix. Let S be any independent
q-bit seed with H_infinity(S)>=q-Delta. Then
TV((E(X,S),S),(U_m,S)) <= 2^Delta*eta. (3)
For each s the conditional distance is a number f(s) in [0,1]. Its mean
under uniform seed is at most eta. Since Pr[S=s]<=2^Delta/2^q,
its mean under S is at most 2^Delta*eta. This also proves the claim with the
ENTIRE q-bit seed retained. It does not claim that the result remains random
when both complete inputs are exposed.
**Lemma W2 (the selected weak-seed bound).** For the AGHP matrix
of Section1 or the canonical family7g, assume Hinf(X)>=m+2e and an independent full q-bit
seed S has entropy>=q-Delta; unused prescribed suffixes may be retained.
Then TV((E(X,S),S),(U_m,S))<=2^((Delta-1)/2)*2^-e.
For each seed let f be its conditional TV. Cauchy--Schwarz pointwise gives
f^2<=2^m/4 times the output collision excess. Section1 and7g each bound
its uniform-seed mean by(1/2)2^-2e. Weighting
by the seed law multiplies this mean by at most2^Delta; Jensen bounds the
mean f by its root. The entire seed is retained. This is special to the
supplied collision calculation, rather than a strengthening of(3) for an
arbitrary seeded map.
## 3. Omitted
The original Section 3 (a local call inside Li's equal-length wrapper) is not
used by this package and is omitted; Section 4 keeps its number.
## 4. Deterministic field setup at the actual outer length
For the local application, d,q,r=O(log n), where n is the actual input length
of the desired final unseeded family. Enumerate monic binary polynomials of
degree r in lexicographic coefficient order. Test irreducibility by the
standard finite test: X^(2^r)=X modulo P, and for each prime divisor p of r,
gcd(P,X^(2^(r/p))-X)=1. These conditions exclude every proper factor degree
and guarantee that the quotient is the field with 2^r elements.
Here are the algebraic details to make the setup premise explicit. Every
polynomial over a field has a finite splitting field: factor off an
irreducible polynomial, form the field quotient by that factor to adjoin a
root, and repeat for the remaining factors. The quotient is a field because
Euclid's identity supplies inverses modulo an irreducible polynomial; an
irreducible factor exists by descending degree. Each extension is finite
when the base field is finite. In a splitting field of X^(2^r)-X, its roots
form a subfield: Frobenius preserves sums/products, and inverses of nonzero
roots are roots. Its derivative is one, so there are exactly 2^r distinct
roots. We have therefore constructed a field F with exactly 2^r elements.
Frobenius x -> x^2 is a permutation of F whose r-th power is identity.
For a in F its orbit length j divides r. Its minimal polynomial degree is
exactly j: the orbit product product_(i=0)^(j-1)(X-a^(2^i)) has coefficients
fixed by Frobenius, hence in F_2 (the fixed equation x^2=x has just roots
0,1), while any binary polynomial vanishing at a vanishes at every orbit
point. No proper-degree polynomial can vanish at j distinct points. Every
element of degree a proper divisor j of r is a root of X^(2^j)-X. If r>=2,
the union of all these proper-degree sets has size at most
sum_(j=1)^floor(r/2) 2^j = 2^(floor(r/2)+1)-2 < 2^r.
Thus some a has degree r, and its minimal polynomial is irreducible of
that degree. The r=1 case is immediate. This proves the enumeration will
terminate without a probabilistic choice of a per-length circuit.
For an arbitrary irreducible binary polynomial Q of degree j, its quotient
has 2^j elements. Each nonzero element has (2^j-1)-st power one by the
elementary finite-group counting proof of Lagrange's theorem, so Frobenius
to power j is identity. For the generating residue of X, its orbit has
length exactly j by the preceding minimal-polynomial argument. Therefore
Q divides X^(2^h)-X exactly when j divides h. The displayed congruence test
makes every factor degree divide r and excludes repeated factors, since
X^(2^r)-X has derivative one. If a factor degree j were a proper divisor
of r, choose a prime p dividing r/j; then j divides r/p, and the relevant
gcd test would detect Q. The tests therefore force a single degree-r
irreducible factor, proving their correctness.
There are 2^r candidates, and each test uses polynomially many r-bit
polynomial operations. Trial division factors r in polynomial time in r.
The search is therefore 2^r*poly(r)=n^O(1). This argument specifically relies
on r=O(log n). It does NOT make an r=Theta(n) field search polynomial.
The matrix can be streamed using N=d*m field powers and row sums, so no
exponential table of seed values is stored or searched.
# Elementary bilinear extraction with both separate strongness guarantees
Here is a finite bilinear implementation of the ordinary strong two-source
primitive used in Section 6.2 of condenser-integration.md
(affine-and-standard-wrappers.md Section 3.4). The construction is standard field
inner-product hashing; the proof is included rather than imported from Li's
Theorem 2.20. Choose integers 1<=m_b<=h_b and one fixed binary polynomial basis
of F_(2^h_b). For field elements x,y define
IP_b(x,y)=Prefix_(m_b)(bin(x*y)).
If independent h_b-bit X,Y have entropies at least k_1,k_2, respectively, then
TV((IP_b(X,Y),Y),(U_(m_b),Y))
<= (1/2) sqrt((2^m_b-1)2^(h_b-k_1-k_2))
<= (1/2)2^((h_b+m_b-k_1-k_2)/2). (IP1)
The symmetric statement retaining X holds with the same error. These are
separate retained-input guarantees, not uniformity after retaining both.
Proof. For each nonzero output character a in F_2^m_b, write
ell_a(z)=a dot Prefix_(m_b)(bin(z)), a nonzero binary linear functional on the
field. For a fixed y put b_a(y)=E_X (-1)^(ell_a(X*y)). If Xhat is an independent
copy of X, orthogonality under uniform field y gives
E_(uniform y) b_a(y)^2
= E_(X,Xhat) E_(uniform y) (-1)^(ell_a((X+Xhat)*y))
= Pr[X=Xhat] = cp(X).
The last equality holds because multiplication by a nonzero X+Xhat is a
permutation of the field and a nonzero binary linear functional is balanced.
Since every y under Y has mass at most 2^(-k_2), comparison with the uniform
field law gives E_Y b_a(Y)^2<=2^(h_b-k_2)cp(X)<=2^(h_b-k_1-k_2).
For the conditional output probability vector p_y, character orthogonality
(Parseval, obtained by expanding the finite sums) gives
sum_z(p_y(z)-2^(-m_b))^2=2^(-m_b) sum_(a!=0)b_a(y)^2.
Cauchy--Schwarz bounds TV(p_y,U) by one half the square root of the rightmost
character sum. Jensen's inequality and the preceding bound prove (IP1).
Field multiplication is commutative, so exchanging X,Y proves the other
strong statement.
A polynomial-basis multiplication algorithm evaluates this map once an
irreducible degree-h_b binary polynomial is fixed. For the local logarithmic
regime h_b=O(log outer n), choose the first irreducible monic polynomial by
exhaustive lexicographic search and the finite modular-squaring/gcd tests
described in weak-seed-hashing.md Section 4. This costs 2^h_b poly(h_b), hence
polynomial in the actual outer n. Small lengths use the same finite search.
For the restricted CZ final extractor set
h_b=N and m_b=h in Section 6.2 of condenser-integration.md, and use (IP1) with k_1=kappa,k_2=k_T. Thus the
ordinary IP extraction analysis itself is supplied here; only the surrounding
deep conditional-entropy and condenser premises remain imports.
## Deterministic field presentation and exact asymptotic bit cost
Use the canonical exact-degree-h_b field constructed in weak-seed Section4
above, which proves existence and correctness of its lexicographic
modular-squaring/gcd test for every positive degree, without advice.
Schoolbook multiplication and reduction of degree-below-h_b binary polynomials
cost O(h_b^2) bit operations. A modular Frobenius test uses at most h_b
squarings; the gcd tests use polynomially many operations on h_b-bit
polynomials. A deliberately loose upper bound per candidate is O(h_b^5) bit
operations, including factoring h_b by trial division. There are 2^h_b
candidates, hence setup costs at most O(2^h_b h_b^5) bit operations and
polynomial working space. Once P is fixed, a call to IP_b costs O(h_b^2) bit
operations and O(h_b) temporary storage; prefixing the result costs O(m_b).
For the intended local use h_b<=C log_2 n with one fixed C, where n is the
actual outer input length. Therefore setup is at most n^C times a fixed
polylogarithmic factor and evaluation is polynomial in n. All finite lengths remain
well-defined by the same search, with no uncomputable exceptional branch.
## Attribution, application, and limits
Field bilinear extraction and its character/collision analysis are standard
inner-product hashing. The needed contract is the one cited by Xin Li,
*Two Source Extractors for Asymptotically Optimal Entropy, and (Many) More*,
ECCC TR23-023 revision 1, Theorem 2.20:
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/
The proved constant strengthens the error bound in condenser-integration.md.
The hypotheses include independence of the two honest sources, integer
1<=m_b<=h_b, and entropy lower bounds 0<=k_1,k_2<=h_b. No claim retains both complete inputs jointly: the output is a
deterministic function of those inputs.
The remaining ordinary/advice condenser and CZ proofs are supplied in the
companion modules; the final assembly is in proof.md.
Fixed rational rates in forward split order:
vMI=(416487/200000,2082427/1000000,21/10).
vMF=(1041243/500000,1041239/500000,21/10).
vMG=(3991127/1000000).
vMJ=(1188521/250000).
vNB0=(28521239/1000000,4151549/125000,16591293/500000,16507271/500000,8025451/250000,3911323/125000).
vNB2=(53086/3125,17748083/1000000).
vL0=(1566157/500000,2125683/500000,4251099/1000000,66369/15625,2100321/500000).
vLr=(1707757/250000,3413761/500000,6803097/1000000,1596337/200000).
vE=(1055443/1000000,25719/12500,32828/15625).
I/F usePXTM(JI=(3,3,3),vMI/vMF); G usesPXTM(JG=(3,),vMG); J usesbinaryvMJ.
Ext0/2 useJS=(4,5,5,5,5,6)/JY=(5,6); L0 usesPXTM(JL0=(3,4,4,4,4),vL0); Lr usesPXTM(JLR=(4,4,4,5),vLr); E usesPXTM(JE=(2,3,3),vE).
Lm usesAGHP. B9 checks every7g/7h guard; vectors are fixed before n.
