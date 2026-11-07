# linear-and-bilinear-primitives.md

Fix Q_S=749254341853/250000000000. With L=floor(log2(actual input length)), S=floor(L/Q_S)
is an integer; S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer ceilings.


Component references resolve through proof.md Section 0.

# Finite seeded and linear extraction: an elementary logarithmic-regime replacement

All logarithms are binary; distance means half-L1.

## 1. Main result and what it replaces

The condenser/block-hash construction supplies the required output/seed ratio.

Fix positive INTEGERS n,h,ell,e,
put B=2^h, and request error 2^(-e). The input length is n and leaf output
length is ell. Put

    a=e+h+4, U=4n, rho=a+2 ceil(log U)+2,
    s0=the least power of 2 at least 2 rho,
    R=4h,
    S=the least power of 2 at least (R+1)rho.

Define integers backwards by

    K_h=ell+2a,
    K_(i-1)=ceil[(2R/(R-1)) (K_i+a+S/2+1)], i=h,...,1. (1)

If K_0<=n, there is one fully specified deterministic strong LINEAR extractor

    H: {0,1}^n x {0,1}^D -> {0,1}^(B ell)

for min-entropy K_0 and error at most 2^(-e). Its exact seed length is
s0+hS+W, where W is defined in Section 5, and the useful finite bounds are

    D <= 6ell+(h+7)S+18a+12, (2)
    K_0 <= 2B(ell+4a+S+4). (3)

For an arbitrary target m, take ell=ceil(m/B) and discard extra output bits.
Every operation is specified using binary arithmetic and quadratic finite-field
towers. Its bit-time is polynomial in n, h, e and B; the outer application fixes
h and has e=O(log outer N), so its evaluator is polynomial in outer N. All field and probability prerequisites are proved below.

This finite GUV/CGL reconstruction supplies arbitrarily large fixed output/seed
ratios with logarithmic outer entropy. It does not assert CL Lemma 3.13;
applications must check the finite threshold (1).

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

## 5. Complete construction and finite ledger

Use the parameters in Section 1. Every condenser input is padded to U=4n
bits, so all condenser seeds and field choices are predetermined.

1. Apply Section 3 with entropy K_0, error 2^(-a), R=1, seed length s0.
   Call the output length m0=s0 ceil(K_0/(s0-rho)).
2. For i=1,...,h, apply Section 3 with entropy K_(i-1), error 2^(-a),
   R=4h, simultaneously to every current block using ONE fresh S-bit seed.
   Its output length is v_i=S ceil(K_(i-1)/(S-rho)). Pad each output if
   necessary and split it into two blocks of length m_i=ceil(v_i/2).
3. Hash all B final blocks to ell bits each with the shared field seed from
   Section 4. Here W is the least power of 2 at least m_h. Concatenate.

The source starts independent of all seeds, which are separate uniform slices
of the extractor's complete seed. The good conditional sources exist only in
the proof. All evaluator steps above are literal deterministic functions.

By (5), m_i <=(1+1/R)K_(i-1)/2+S/2+1. Recurrence (1) implies
K_(i-1)>=m_i+K_i+a, exactly the splitting condition. The threshold decreases
with i, and K_0>=S. The initial output m0<=2K_0+s0<=3K_0. Each later block
has length at most 3K_0 as well. Thus all can indeed be padded to U=4n when
K_0<=n, without deleting input bits.

The initial error is 2^(-a). At level i there are 2^(i-1) blocks, and the
condenser and split each contribute at most 2^(i-1)2^(-a). The final hash
contributes at most (B/2)2^(-a). The sum is

    [1+2(B-1)+B/2]2^(-a) <=3B2^(-a) <2^(-e).

All seeds are retained throughout, so the result is STRONG. Every stage is
linear in the original input when seeds are fixed, so the result is LINEAR.

For completeness, let lambda=2R/(R-1). R>=4 gives lambda<=8/3. From (1),

    K_(i-1) <= lambda(K_i+a+S/2+2).

Also lambda^h=2^h(1+1/(4h-1))^h<=2B. Summing the geometric recurrence gives
K_0<=2B(K_h+2a+S+4), which is (3). At the final split, the same inequalities
and (R+1)/(R-1)<=5/3 give

    m_h <=3(K_h+a+S+2)=3(ell+3a+S+2).

Thus W<=6ell+18a+6S+12. Since s0<=S, adding s0+hS proves (2).

### Explicit bit-time bound

Check K_0<=n before constructing any fields. On admissible parameters the
recurrence implies B<=n, h<=log n, a<=n and S<=n. Every condenser tower has
binary dimension sl<=8n; the final hash has binary dimension W<=6n.
The orbit polynomial uses O(l*s) tower squarings and O(l^2) tower products.
Its coefficients lie in the initial s coordinate positions, so no expensive
basis search is needed. Polynomial multiplication and reduction modulo the
monic E use classical schoolbook arithmetic. The consecutive PV powers need
b*u<=k+u<=2n modular squarings. Evaluation is Horner's rule. There are at
most B condenser calls, B hashes, and polynomially many integer operations.
With the four-product tower arithmetic above, a conservative O(n^8) BIT-TIME
bound covers all parameter calculation, field setup and evaluation on every
admissible input. Invalid
parameter inputs can be rejected immediately after the integer guard.

## 6. Selected scaling

The selected parameters are the exactPXT/PXG vectors of7g and affine-only
three-condenser/kernel base7f. Proof.md supplies all closed slots,feedback guards,sampler/error
accounting and all-length calibration.

## 7. A generalized finite recursive primitive (Theorem D3)

This section extends Sections 2-5; the theorem of Section 1 is its case j=2,
R=4h, U1=4n with powers of two. TheselectedbaseusesSection7f;its three-seed public-matrix/kernel construction is proved there. SelectedNBExt0andExt2usePXTofSection7g,
atj2. CurrentouterprimitivesuseexactPXG/uniformAGHP;F#usesPXG,NBusesPXT. lg(x) is the least integer i>=0 with 2^i>=x.

Parameters: positive integers n, e, ell; an integer j>=2 that is a power of
two; integers R>=max(4,j) and h'>=1; an integer U1>=1; a dimension rule,
either "powers of two" or "admissible" (integers r0*2^i with r0 in
{1,3,5,7}, i>=0). adm(x) is the least allowed dimension >=x. Put

    h=h' log2 j, B=2^h=j^h', a=e+h+4, U=4n,
    rho=a+2lg(U)+2, s0=adm(2rho), rho1=a+2lg(U1)+2, S=adm((R+1)rho1),
    K_h'=ell+2a, K_(i-1)=ceil[(jR/(R-j+1))(K_i+a+(j-1)(S/j+1))], i=h',...,1.

Evaluator. (1) Section 3's condenser on the n-bit source (padded to U),
entropy K_0, error 2^-a, R=1, field dimension s0; output
m0=s0*ceil(K_0/(s0-rho)) bits. (2) For i=1..h': zero-pad every current block
to U1 bits, apply Section 3's condenser (entropy K_(i-1), error 2^-a,
parameter R, ambient U1, field dimension S) to every block with ONE fresh
S-bit seed, zero-pad each output to a multiple of j and split it into j
blocks of m_i=ceil(v_i/j) bits. (3) Hash all B blocks to ell bits each with
one shared seed in the field of dimension W=adm(m_h') (Section 4's hash).
Seed length D=s0+h'S+W.

**Theorem D3.** If K_0<=n and 2K_0+s0<=U1, this is a strong linear extractor
for min-entropy K_0 with error <2^-e and output B*ell. With
Lambda_j=(R/(R-j+1))^h', kappa_j=(R+1)/(R-j+1), gamma=2 (powers of two) or
5/4 (admissible) and f_j=jR/((j-1)(R+1))<=2:

    D < s0+h'S+gamma*m_h', m_h'<=kappa_j(ell+3a+(j-1)(S/j+1)+1)+S/j+1, (2')
    K_0 <= B*Lambda_j*(ell+(2+f_j)a+f_j(j-1)S/j+f_j*j). (3')

Proof (only the changes to Sections 2-5).
(a) Admissible fields. Run Section 2's tower from F_0=F_2[x]/(f) with
f=x^3+x+1, x^5+x^2+1 or x^7+x+1 and a_0=1. Each f is irreducible: none has
a root in F_2 (f(0)=f(1)=1); x^5+x^2+1 and x^7+x+1 leave residue 1 modulo
x^2+x+1 (there x^3=1, so x^5=x^2 and x^7=x); x^7+x+1 leaves residue x
modulo x^3+x+1 and modulo x^3+x^2+1 (x^7-1=(x+1)(x^3+x+1)(x^3+x^2+1), so
x^7=1 modulo either cubic); a reducible polynomial of degree <=7 has an
irreducible factor of degree <=3, and x, x+1, x^2+x+1, x^3+x+1, x^3+x^2+1
are all of them. The absolute trace of 1 in F_(2^r0) is r0 mod 2=1, so
a_0=1 starts the induction, which is otherwise verbatim and gives fields of
binary dimension r0*2^i. The polynomial-basis extension E(T) of degree l (a
power of two) over F_q, q=2^s with s=r0*2^i, is built verbatim: the
immediately preceding tower level is the subfield with q^(l/2) elements.
Consecutive admissible integers a<b with a>=4 satisfy b<=5a/4 (pattern 4,5,6,7,8 times
powers of two), so adm(x)<5x/4 for x>4; adm(x)<=pow2ceil(x) always. The
hash of Section 4 works in any field (z(x-x') is uniform for x!=x').
(b) Section 3 at ambient U1 and dimension s: its proof uses k<=U, input
length <=U, rho=a+2lg(U)+2 and s>=(R+1)rho only (H=2^u stays a power of
two); (5) holds with s<gamma(R+1)rho in place of s<2(R+1)rho.
(c) j-way split. Pad a v-bit k-source to jm bits, m=ceil(v/j), blocks
Z_1..Z_j. If k>=(j-1)m+k'+a, the tuple is within (j-1)2^-a of a block
source with entropy k' per block. Indeed (Z_1..Z_i) has min-entropy
>=k-(j-i)m, so the average over z_<i of max_z P[Z_i=z|Z_<i=z_<i] is
<=2^((i-1)m-k+(j-i)m), and by Markov the prefixes with conditional
min-entropy <k' have mass <=2^-a. Generate Z' block by block: equal to Z's
conditional law on good prefixes, uniform on bad ones (and on prefixes Z
cannot reach); coupled with Z until the first bad prefix, which has
probability <=(j-1)2^-a. Z_1 needs no repair (entropy >=k-(j-1)m>=k'+a), and
uniform blocks have entropy m>=k' because k<=v<=jm. Conditioning on
external variables is as in Section 4.
(d) Ledger. By (5), v_i<=(1+1/R)K_(i-1)+S and m_i<=v_i/j+1, so the
recurrence gives K_(i-1)(R-j+1)/(jR)>=K_i+a+(j-1)(S/j+1), i.e.
K_(i-1)>=(j-1)m_i+K_i+a: the split condition. Hence m_i<=K_(i-1)<=K_0<=U1
for i>=1, m0<=2K_0+s0<=U1, and each level-i entropy K_(i-1)<=U1: every
padding is legal. Errors: 2^-a (level 0), j^(i-1) condensers and j^(i-1)
splits at level i (<=j^i 2^-a), hash (B/2)2^-a; total
<=(1+jB/(j-1)+B/2)2^-a<=3B 2^-a<2^-e.
(e) (3'): with lambda=jR/(R-j+1)>=j>=2 and c'=a+(j-1)(S/j+1)+1
=a+(j-1)S/j+j, K_(i-1)<=lambda(K_i+c'), so K_0<=lambda^h'K_h'
+c'(lambda+...+lambda^h')<=lambda^h'(K_h'+f_j*c') with
f_j=lambda/(lambda-1)=jR/((j-1)(R+1)) (as lambda-1=(j-1)(R+1)/(R-j+1));
this is (3'). (2'): W=adm(m_h')<gamma*m_h', and
m_h'<=(R+1)K_(h'-1)/(jR)+S/j+1 with K_(h'-1)<=lambda(ell+3a+(j-1)(S/j+1)+1).
Strongness and linearity are unchanged (every step is linear for fixed
seeds; all seeds are retained). Bit time is polynomial: every field
dimension is at most gamma*max(2rho,(R+1)rho1,m_h'), and each extension tower
of Section 3 has dimension s*l<=2max(U,U1)+2s.

**Corollary D3'** (the envelope behind Lemma P). Take j=2, U1=U=4u,
height h, R=R'>=4, admissible dimensions, leaf ell=ceil(z/B),
b=e+h+4, rho_u=b+2lg(4u)+2, S'_u=adm((R'+1)rho_u), Lambda'=(R'/(R'-1))^h,
kappa_2=(R'+1)/(R'-1), gamma=5/4. (2K_0+s0<=U1 holds when K_0<=u, since
s0<=S'_u<=K_0.) Then, on input length u and output z,

    K(u,z) <= Lambda'z+G'_u, G'_u=Lambda'B(4b+S'_u+5),
    D(u,z) < gamma*kappa_2*z/B+P''_u,
    P''_u=(h+1)S'_u+gamma*kappa_2(3b+S'_u/2+3)+gamma(S'_u/2+1).

(K: from (3') with f_2<=2 and B*ell<=z+B-1. D: s0<=S'_u, W<gamma*m_h,
m_h<=kappa_2(ell+3b+S'_u/2+2)+S'_u/2+1 and ell<=z/B+1.) Since
P''_u<=(h+4)(S'_u+b+1) for h>=11 (gamma*kappa_2<=25/12) and
G'_u>=Lambda'B(S'_u+b+1), 128P''_u<=G'_u whenever Lambda'2^h>=128(h+4).
Both G'_u and P''_u are nondecreasing in u.

**Lemma P** (one primitive in a unit w). Let P[h,R] be the primitive of
Theorem D3 at j=2, height h>=1 (h'=h), parameter R>=4, admissible
dimensions, U1=U=4u, leaf ceil(z/2^h), error exponent e, output truncated to
z bits: input u bits, exact threshold K(u,z) and exact seed length D(u,z)
from its algorithm. Put gamma=5/4 and

    Lam=(R/(R-1))^h, kap=(R+1)/(R-1), c=gamma*kap/2^h, g=Lam*2^h*(4+gamma(R+1)),
    pi=(h+1)gamma(R+1)+gamma*kap*(3+gamma(R+1)/2)+gamma*(gamma(R+1)/2+1/2).

For every integer w>=max(2, e+h+6+2lg(4u)):

    K(u,z) <= Lam*z+g*w, D(u,z) < c*z+pi*w. (P)

If K(u,z)<=u, P[h,R] is a strong linear (K(u,z),2^-e) extractor; a longer
seed slot ignores (and retains) its suffix.

Proof. Corollary D3' gives K<=Lam*z+G'_u and D<c*z+P''_u, with b=e+h+4,
rho_u=b+2lg(4u)+2=e+h+6+2lg(4u)<=w, S'_u=adm((R+1)rho_u)<gamma(R+1)rho_u
<=gamma(R+1)w ((R+1)rho_u>4 is an integer and adm(x)<5x/4 for x>4),
G'_u=Lam*2^h(4b+S'_u+5) and P''_u=(h+1)S'_u+gamma*kap(3b+S'_u/2+3)
+gamma(S'_u/2+1). Its derivation uses no lower bound on h (only its final
remark 128P''<=G' needs h>=11; it is not used here). Since w>=b+6 and w>=2:
4b+5<=4w, 3b+3<=3w and 1<=w/2; substitute. Legality: Theorem D3 needs
K_0<=u and 2K_0+s0<=U1=4u, and the second follows from the first because
s0<=S'_u<=K_0 (the recurrence gives K_0>=(2R/(R-1))(S'_u/2+1)>S'_u).

## 7a. Exact binary field degrees for the base and both NB primitives

Define D3_exact by the same parameters, entropy recurrence, condensers,
splitting, guards and concatenation as Theorem D3, but take exact binary
field degrees s0=2rho,S=(R+1)rho1,W=m_h'. For each positive degree v
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
built by towers, not enumeration in their full dimension. The final hash
uses the canonical field of exact degree W and the same collision proof.

D3's expansion, splitting, conditional coupling, error and linearity
proofs require s0>=2rho,S>=(R+1)rho1 and field arithmetic, so they still
apply. Its recurrence and bound(3') hold; its final split bound on m_h'
in(2') holds. The exact seed is now D=2rho+h'S+m_h', so the seed bound
uses gamma=1 and <=, without the admissible rounding factor5/4. All seeds
are retained. Both K0<=u and 2K0+s0<=U1 are checked by the algorithm.

For each binary exact PX primitive put j2,h'=h,R>=4,U1=4u,ell=ceil(z/2^h),
then truncate to z. Call this PX[h,R]. For integer unit
w>=max(2,e+h+6+2lg(4u)), define
 Lam=(R/(R-1))^h,kap=(R+1)/(R-1),B=2^h,
 c=kap/B,g=Lam*B*(R+5),
 pi=2+h(R+1)+kap*(3+(R+1)/2)+(R+1)/2+1/2.
Then

    K(u,z)<=Lam*z+g*w, D(u,z)<=c*z+pi*w. (PX)

Proof. rho=e+h+6+2lg(4u)<=w,b=e+h+4,w>=b+6,S=(R+1)rho.
Bound(3') gives K<=Lam*z+Lam*B*(4b+S+5);4b+5<=4w gives the result.
For the seed m_h<=kap(ell+3b+S/2+2)+S/2+1. Substitute
ell<=z/B+1,3b+3<=3w,1<=w/2 into D=2rho+hS+m_h.
The same recurrence gives K0>=S>=s0, so K0<=u implies
2K0+s0<=3u<4u.

F#/NB Ext0,Ext2 and outer PXT/PXG/AGHP have local input/output/error and
enumerated degrees O(log actual N) by their recurrences, closed slots and B9.
BASE7f has b0,b=O(log N), U1=2^(90+lg S)=O(log N); towers build its initial
O(N)-dimensional relative extension. Later ambient/matrix/kernel dimensions
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

For actual n>=1 fix Q_S=749254341853/250000000000, beta=8, delta_O=186359055259,
eps_B=157/2^38, eps_R=2^-40, N0=2^41, a0=1+2^-41,
gamma_B=(8-a0*(Q_S+1+eps_B)-3eps_R)/(2(1+eps_B)),
xi_B=gamma_B/(gamma_B-1), G_B=1-(xi_B-1)^2,
A_B=1/G_B, B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)=1345.
First compute L=floor(log2 n), S=floor(L/Q_S), d=delta_O*L.
If S<2^40 return 0^d BEFORE evaluating any lg S. Otherwise set
these ACTUAL integer parameters in the following order:
 s=2ceil(beta*S/2)=8S,
 ellc=lg(s+255), r_F=S+8ellc, a_B=r_F+3,
 U=4n, U1=2^(90+lg S), rho=a_B+2lg(U1)+2=r_F+185+2lg S,
 d=delta_O*L, K_B=ceil(A_B*d+B_B*rho),
 rho0=a_B+lg U+lg K_B+2,
 b0=ceil((N0+1)rho0/N0), u0=b0-rho0,
 b=ceil(gamma_B*rho), u=b-rho,
 v0=b0*ceil(K_B/u0), v1=b*ceil(K_B/u),
 c=v1-K_B, v2=b*ceil(c/u), D=b0+2b. (AF1)
First require S>=2^40. Test K_B<=n, rho0>=N0, u0>=1, u>=1,
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
 gamma_B=590735969153879965910620377279571/295147905347930292224000000000000>2, xi_B=590735969153879965910620377279571/295588063805949673686620377279571,
 G_B=260017433279356697616979501287034414431296556106419914693944041/87372303464550175490575249188365069280607296556106419914693944041>0, G_B*A_B=1,
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

Complete all-length calibration. For 2^k<=S<2^(k+1), k>=40,
s+255<9S implies ellc<=k+5. Therefore
 rho-S<=10k+227<eps_B*S.
The exact coefficient bound is
 K_B/S<=A_B*delta_O*Q_S*(1+2^-40)+B_B*(1+eps_B)+eps_R
 =13413847991861855304161360965324103717628046435565872889594198770784082604854719791329941109767374439/71473047828780739083743237909132974003771465430688168272230368342341320704000000000000<2^48.
Here L<Q_S(S+1), rho<=(1+eps_B)S, and eps_R*S>=1 pay K_B's ceiling.
Thus lg K_B<=48+lg S. Since lg(4n)<=L+3,
 rho0-(Q_S+1)S<=9k+Q_S+97<eps_B*S.
The k40 bases are 627<628 and 457+Q_S<628. Their affine left sides
increase by10 and9, at most their bases, while eps_B*2^k doubles.
This proves both bounds for EVERY integer S>=2^40. Each field-degree
ceiling costs <=1<=eps_R*S. By the definition of gamma_B,
 D/S<=a0*(Q_S+1+eps_B)+2gamma_B*(1+eps_B)+3eps_R=8. (AF2)
Thus the COMPLETE three-slice seed and ignored suffix fit s=8S.
Put b0coef=a0*(Q_S+1+eps_B)+eps_R. The exact inequalities
 (2^41+1)2^48+b0coef<2^90,
 xi_B*2^48+gamma_B*(1+eps_B)+eps_B<2^90
imply v0,v1<2^90*S<=U1, and 0<c<v1<=U1. Also
 K_B<=[A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R]L
 <=H_base*L, H_base=62621262444658. (AF4)
The last allowance pays the K_B ceiling. H*L<=n and H>=H_base give
K_B<=n. Finally rho0>=(Q_S+1)S>2S>=2^41=N0 and rho>=S>=2^40 prove every
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

## 7e. A directly proved exact primitive at rational rates

Define PXR[h,R] for integerh>=1 and fixed rationalR>1,inputu,outputz,errore>=1.
Put B=2^h,ell=ceil(z/B),a=e+h+4,U=4u,U1=4u,
rho=a+2lg(4u)+2,s0=2rho,S=ceil((R+1)rho),lambda=2R/(R-1),
Kh=ell+2a,Ki-1=ceil[lambda*(Ki+a+S/2+1)],i=h,...,1.
CheckK0<=u,else outputzero. Use7a's exact canonical fields/towers.
Initially condense at entropyK0,degree s0; performh binary condenser/split
levels of degreeS with one shared freshS-bit seed per level; hash allB
leaves toell bits with one shared exact degree-mh field seed,then truncate
toz bits. All seeds are retained. This defines a new map directly; no
integerR>=4 premise of PX or D3 is applied outside its stated domain.

Here is its full finite guarantee. Section3's integer condenser expansion
uses onlys>rho and q/H=2^rho>=4U^2*2^a; its ratio atS satisfies
S/(S-rho)<=(R+1)/R. Ifvi is its output,mi=ceil(vi/2)
<=(R+1)Ki-1/(2R)+S/2+1. Since1-(R+1)/(2R)=1/lambda,
the reverse recurrence givesKi-1>=mi+Ki+a. Thus4's binary split and
causal shared-seed repair apply. Initialoutput<=2K0+s0;K0>=S>=s0,
soK0<=u implies2K0+s0<=3u<U1. Every split block has length<=K0.
The same shared-seed hashing collision proof and level union bound give
error<3B*2^-a<2^-e,strongness and source-linearity,retaining all seeds.

For integerunitw>=max(2^20,e+h+6+2lg(4u)),put
 Lam=(R/(R-1))^h,kap=(R+1)/(R-1),c=kap/B,g=Lam*B*(R+5+2^-20),
 pi=2+h(R+1)+kap*(3+(R+1)/2)+(R+1)/2+1/2
    +(h+kap/2+1/2)*2^-20.
ThenK<=Lam*z+g*w andD<=c*z+pi*w.
To prove these bounds,geometric summation withf=lambda/(lambda-1)<=2
andlambda^h=B*Lam gives
 K0<=B*Lam*[ell+(2+f)a+f*S/2+2f]
 <=Lam*z+B*Lam*(4a+S+5)<=Lam*z+g*w.
The last inequality usesrho<=w,a<=w-6,ell<=z/B+1 and
S<= (R+1)w+1; the additional Lam*B is paid by Lam*B*w/2^20.
The final split givesmh<=kap*(ell+3a+S/2+2)+S/2+1;
Substitute intoD=2rho+hS+mh,using3a+3<=3w and1<=w/2.
The integerdegreeceiling adds at most h+kap/2+1/2 seedbits,
paid by the displayed additional coefficient timesw/2^20.
All exact seed/threshold/padding guards hold. At the selected logarithmic
component dimensions, canonical searches cost2^O(logactualN) and relative
extensions use towers, giving polynomial time in actualN.

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
leaves to ell bits with one shared exact degree-m_h seed, then truncate
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
(B/2)2^-a, giving (5B/2-2)2^-a<3B2^-a<2^-e. Strongness and
source-linearity follow as in Section4, retaining every actual seed.

Here are sharper FINITE envelopes. For integer
 w>=max(2^20,e+h+6+2lg(4u)),rho<=w,a<=w-6.
Define A_h=2,B_h=1/B and, backwards,
 A_(i-1)=lambda_i*(A_i+1)+(R_i+1)/2,B_(i-1)=lambda_i*B_i.
Put Lam=B_0=product_i R_i/(R_i-1),g=A_0. Then
 K_0<=Lam*z+g*w. (PT1)
Indeed K_h<=z/B+2w-11. The two integer ceilings in PT0 give
 K_(i-1)<[2(S_i-rho)/(S_i-2rho)]*(K_i+a+1)+S_i/2
           <=lambda_i*(K_i+a+1)+S_i/2. (PT1a)
Here 2(S_i-rho)/(S_i-2rho)<=lambda_i follows from
S_i>=(R_i+1)rho by cross-multiplication. Since a+1<=w-5 and
S_i/2<=((R_i+1)w+1)/2, induction proves PT1: the residual
-5lambda_i+1/2 is negative. Choosing b_i first avoids multiplying
the entire field-output rounding allowance S_i/2 by lambda_i.

For kappa=(R_h+1)/(R_h-1),put
 c=kappa/B,
 C_h=3/2-16kappa,
 pi=sum_i(R_i+1)+3kappa+(R_h+1)/2
      +max(0,h+C_h)/2^20.
The exact seed length satisfies
 D<=c*z+pi*w. (PT2)
For the terminal block, PT0 gives
 m_h<kappa*(ell+3a+1)+S_h/2+1
     <=kappa*z/B+[3kappa+(R_h+1)/2]w+C_h. (PT2a)
The last inequality uses ell<=z/B+1 and a<=w-6, including both
ceiling allowances. Now D=sum_iS_i+m_h and
sum_iS_i<=sum_i(R_i+1)w+h.
Any positive h+C_h is paid by its coefficient timesw/2^20; a negative
constant is simply discarded. This proves PT2 without a hidden output
or seed ceiling. Padding to the envelope is legal; ignored bits remain retained.

Section7a supplies canonical fields: selected base/hash degrees are
O(logactualN), and O(actualN) relative extensions use quadratic towers.
Thus these fixed-rate calls are polynomial in actualN, with no joint-polynomial
bound for unrelated component dimensions.

**PXG terminal variant.** PXG(v) keeps PXT's vector, thresholds and
condensers; replace the leaf hash by weak-seed Section1 of this file's AGHP matrix on
mh bits, output ell, exponent a, fresh shared seed of length
2rH=2(lg(mh*ell+1)+ell+2a). Section4's reverse hybrid requires,
conditional on retained prior seeds T, a block source independent of the
fresh uniform hash seed Y, and a strong single-block bound averaged over Y.
The causal repairs precede Y and give every leaf conditional entropy
K_h=ell+2a, including after its full repaired prefix. AGHP Section1 below
therefore applies to each such law with d=mh,m=ell,e=a and error<2^-a,
retaining Y. Replace leaves B,...,1, retain T,Y and earlier source leaves,
then project to earlier hashes and append independent uniform later outputs.
Summing these B average-seed bounds costs <B2^-a; independent unused
suffixes remain retained. Earlier costs are
[2(B-1)]2^-a, so total <(3B-2)2^-a<2^-e. Fixed-seed linearity composes;
K0,Lam,g and all padding/entropy guards are unchanged. Retain every new
seed bit and suffix: this defines a new algorithm, not seed deletion.
For the seed envelope require 5/3<=Rh<=64 and z<=256*2^h*w. Then
ell<=256w, kap<=4, (Rh+1)/(2Rh)<=4/5, Sh<=65w+1; the mh bound gives
mh<=1199w+9<=2048w. Thus lg(mh*ell+1)<=20+2lgw<=w/2^14:
on [2^k,2^(k+1)), k>=20, use 22+2k<=2^(k-14), true at20 and
preserved by doubling. Adding all seed slices gives
 D=sum Si+2rH<=c*z+pi*w,
 c=2/2^h, pi=4+sum_i(Ri+1)+2^-13+(h+2)eps.
Here a,rho<=w; ell<=z/2^h+1 pays two ceiling bits and Si pay h.
AGHP's proof and canonical fields are below. At selected calls
rH=O(log n), hence setup/evaluation are uniform polynomial time.

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
extractor of Section1, assume Hinf(X)>=m+2e and an independent full q-bit
seed S has entropy>=q-Delta; unused prescribed suffixes may be retained.
Then TV((E(X,S),S),(U_m,S))<=2^((Delta-1)/2)*2^-e.
For each seed let f be its conditional TV. Cauchy--Schwarz pointwise gives
f^2<=2^m/4 times the output collision excess. The preceding uniform-seed
collision calculation bounds its mean by less than(1/2)2^-2e. Weighting
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

The field setup can be specified without a random choice or per-length advice.
Enumerate monic binary degree-h_b polynomials P in lexicographic coefficient
order, choosing the first with both properties:

1. Z^(2^h_b)=Z modulo P.
2. For every prime divisor ell of h_b,
   gcd(P,Z^(2^(h_b/ell))-Z)=1.

For h_b=1 the prime-divisor list is empty; use the same lexicographic rule
and tests. Ordinary trial division of the integer h_b lists its prime divisors.
Repeated squaring modulo P and Euclid's algorithm perform the tests. The first
test makes P divide a polynomial with derivative one, so its factors are
squarefree and each irreducible factor has degree dividing h_b. If a factor
had proper degree j dividing h_b, some prime ell dividing h_b/j would make
j divide h_b/ell; the second test excludes that factor. Thus a passing P is
irreducible. Conversely every irreducible degree-h_b P passes those tests. Existence of an irreducible polynomial of every positive degree is the
standard finite-field existence theorem, not a nonuniform choice of a table.

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

Selected fixed ordered rational rate vectors,in FORWARD splitting order:
PXT(v) and PXG(v) use their entire Section7g vector;these constants are fixed before n.
 vMI=(5941137/5000000,11875281/10000000,5940167/5000000,11872563/10000000,1483507/1250000,11862179/10000000,16666667/10000000).
 vMF=(12826263/10000000,12802221/10000000,12808823/10000000,400159/312500,12829973/10000000,12802799/10000000,6378311/5000000,4166667/2500000).
 vMG=(29222817/10000000,13812957/5000000).
 vMJ=(37132243/10000000,34685317/10000000).
 vNB0=(27646149/1000000,276460263/10000000,276471543/10000000,276438699/10000000,276544001/10000000,34539999/1250000,1727279/62500,276412459/10000000,276444889/10000000,138180279/5000000,138175961/5000000,69033073/2500000,69004041/2500000,4303047/156250,274213091/10000000,271670427/10000000,266886261/10000000,257438723/10000000,120257419/5000000,42668011/2000000,93621427/5000000).
 vNB2=(41483973/10000000,41513673/10000000,41470411/10000000,20670899/5000000,20592331/5000000,40733253/10000000,3938283/1000000,54667481/10000000).
 vL0=(26988893/10000000,47654/15625,16041411/5000000,2917227/1000000,33677157/10000000,28681949/10000000,15661699/5000000,7010659/2500000,29099643/10000000,32259319/10000000,2848829/1000000,13993429/5000000,31631193/10000000).
 vLr=(6599701/1000000,2053963/312500,66425807/10000000,66014289/10000000,65724087/10000000,65120379/10000000,65502849/10000000,1603469/250000,64106441/10000000,14885741/2500000).
 vE=(2604709/2500000,2195549/2000000,2678239/2500000,5362737/5000000,438451/400000,2096013/1250000).
The inner usesvMI/vMF/vMG/vMJ at heights7/8/2/2;NB Ext0/2 usevNB0/vNB2 at21/8;outerL0/Lr/E usevL0/vLr/vE at13/10/6. OuterLm is the separately provedAGHP. Every rate exceeds1;PT1/PT2 use each actual ordered rate. No affine-only primitive is used inside NB.
