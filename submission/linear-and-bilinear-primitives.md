# linear-and-bilinear-primitives.md

Component references resolve through proof.md Section 0.

# Finite seeded and linear extraction: an elementary logarithmic-regime replacement

All logarithms are binary; distance means half-L1.

## 1. Main result and what it replaces

Small-bias hashing alone has seed longer than output, obstructing alternating
extraction and sampler feedback. The condenser/block-hash construction resolves this.

Here is the finite replacement theorem. Fix positive INTEGERS n,h,ell,e,
put B=2^h, and request error 2^(-e). The input length is n and leaf output
length is ell. Put

    a=e+h+4,  U=4n,  rho=a+2 ceil(log U)+2,
    s0=the least power of 2 at least 2 rho,
    R=4h,
    S=the least power of 2 at least (R+1)rho.

Define integers backwards by

    K_h=ell+2a,
    K_(i-1)=ceil[(2R/(R-1)) (K_i+a+S/2+1)],  i=h,...,1.       (1)

If K_0<=n, there is one specified deterministic strong LINEAR extractor

    H: {0,1}^n x {0,1}^D -> {0,1}^(B ell)

for min-entropy K_0 and error at most 2^(-e). Its exact seed length is
s0+hS+W, where W is defined in Section 5, and the useful finite bounds are

    D <= 6ell+(h+7)S+18a+12,                                  (2)
    K_0 <= 2B(ell+4a+S+4).                                   (3)

For target m use ell=ceil(m/B) and truncate. Binary arithmetic and quadratic
field towers give polynomial bit-time in n,h,e,B; fixed h and e=O(log outer N)
give polynomial time in outer N. The algorithm uses no truth-table search,
randomized setup, advice or unproved extractor.
This weaker finite GUV/CGL-style construction does not assert CL Lemma 3.13's
universal entropy constant or variable-t tradeoff. It gives arbitrarily large
fixed output/seed ratios at logarithmic outer entropy; applications must
check threshold (1), not merely big-O bounds.

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

The F_2-linear map x->x^2+x has kernel {0,1}, hence half-field image,
all of trace zero since Tr(x^2)=Tr(x). By Tr(a_i)=1, trace is nonzero,
so its kernel also has half-field size. Thus a_i is outside the image:
the polynomial is irreducible. Roots z_i,z_i+1 give relative trace a_i
for a_i*z_i, hence absolute trace 1 for a_(i+1), completing induction.

Use ordered tower coordinates, with constants embedded by appending zero
coordinates. Addition is XOR. For x=x0+x1 z_i and y=y0+y1 z_i compute
p=x0*y0, q=x1*y1, r=(x0+x1)*(y0+y1), t=a_i*q in the preceding field.
The product is (p+t)+(r+p)z_i, because z_i^2=z_i+a_i. These four recursive
products give T(v)<=4T(v/2)+O(v), hence O(v^2) bit operations. Inverses can be computed
by exponentiation to 2^v-2 in a field of binary dimension v, using O(v)
squarings/multiplications. Zero inversion is never requested.

For polynomial-basis extension degree l>=2, with l and base dimension s
powers of 2, build the tower to sl and take its last root g.
It is not in the immediately preceding subfield.
Its q-Frobenius orbit, q=2^s, has length l: a shorter orbit divides l,
hence divides l/2, forcing g into the preceding q^(l/2)-element subfield,
which fills all roots of x^(q^(l/2))=x, a contradiction.

Thus

    E(T)=product_(j=0)^(l-1) (T-g^(q^j))

has coefficients in F_q, the q-Frobenius fixed field. Any F_q-polynomial
vanishing at g vanishes on its l orbit points; none of degree below l can.
Thus E is irreducible of degree l. Compute the product in tower coordinates
and read the embedded base-field coefficients: canonical polynomial time
in sl, without factoring, primitive-element search or enumerating 2^(sl)
elements. All coordinate orders and choices are fixed.

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
                       (f^(H^(b-1)) mod E)(y)).             (4)

The output has v bits. Each map f -> f^(H^j) mod E is F_2-linear, because
H is a power of 2. Evaluation and padding are linear too. Also

    s < 2(R+1)rho,  u>=Rs/(R+1),
    v <= (1+1/R)k+s.                                        (5)

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
modulo E thus says that f, as an element of F_q[T]/E, is a root of

    sum_(j<J) (p_j mod E) Z^j.

This is a nonzero polynomial of degree <=J-1, so there are at most J-1 such
inputs. In particular a J-input set cannot have fewer than AJ neighbors.
The same statement holds after restricting the source to the padded inputs.

Here l<=2U and b<=k<=U. Also q/H=2^rho>=4 U^2 2^a. Thus

    ((l-1)(H-1)b)/q <= 2U^2/2^rho <= 2^(-a-1).

Thus A>0 and A>=(1-2^(-a))q. Since H^b>=2^k, the required size-2^k
sets have at least (1-2^(-a))2^k q distinct neighbors.

### Strongness, not only joint min-entropy

For a source uniform on a size-K set, K=2^k, consider each y separately.
Keep one preimage for each distinct value Con(x,y). To obtain a K-element
set of output values at that y, replace the duplicates by unused values.
There is room, because v>=k (s>=u and b>=k/u). The number replaced,
summed over all y, is <= 2^(-a)Kq by expansion. The repaired conditional
output at each y is uniform on K values, and the seed marginal is unchanged.
The joint distance is <= 2^(-a). An arbitrary k-source is a convex
combination of uniform size-K sources: iteratively express a probability
vector with coordinates <=1/K as a mixture of vertices of that capped
probability simplex; its vertices have exactly K entries 1/K. Mixing the
repairs keeps each conditional atom <=1/K and preserves the distance bound.
This proves the asserted strong condenser, including the seed marginal.

## 4. Three exact block operations

A block source (X_1,...,X_j) has block entropy k if every conditional law
X_i | X_<i=x_<i has min-entropy >=k on its support. All statements below
also hold conditioned on any fixed external variable.

**Apply a strong condenser with one shared seed.** Apply error-eta
Con(X_i,Y) to all j blocks with one independent uniform Y. Retaining Y,
the tuple is within j eta of a block source with the guaranteed entropy.
For each original prefix and seed, couple the actual output to a repaired
conditional law, using fresh independent coupling coins. Keep original
blocks latent. Given their prefix and seed, earlier repairs reveal nothing
about the next original block. Each repaired block has the required atom
bound, preserved on mixing latent prefixes. Each repair disagrees with
probability <=eta on average, since the original prefix is independent of Y.
Union bound gives j eta. Couplings are not evaluated.

**Split one high-entropy block.** Pad a v-bit k-source by one zero if needed
and split into two m-bit halves, m=ceil(v/2). If k>=m+k'+a, the pair is within
2^(-a) of a block source with entropy k' in each half. The first half already
has entropy >=k-m>=k'+a. The average maximum probability of the second half
given the first is <= 2^(m-k). Markov bounds by 2^(-a) the mass of first
halves whose conditional second-half entropy is below k'. Replace the
second-half law by uniform only on those prefixes. The first-half law is
unchanged. Applied to j input blocks with the causal latent-prefix coupling
above, the total distance is <= j2^(-a). This proof also explains why
arbitrary conditioning on a rare prefix is not silently allowed.

**Hash every block with one shared seed.** For m-bit blocks, put W equal to
the least power of 2 at least m. Use the tower field of binary dimension W.
Pad each block to W bits. For a uniform z in that field define Hash_z(x)
as the first ell coordinates of zx. For distinct x,x', z(x-x') is uniform,
so collision probability is exactly 2^(-ell).

For a single source of entropy >=ell+2a, its conditional output collision
probability, averaged over z, is <= 2^(-ell)+2^(-ell-2a). Cauchy-Schwarz
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

Check K_0<=n before field setup. The recurrence gives B<=n,h<=log n,a<=n,S<=n.
Condenser towers have dimension sl<=8n; the hash has W<=6n.
The orbit polynomial needs O(l*s) squarings and O(l^2) products; coefficients
are in the initial s coordinates, requiring no basis search.
Schoolbook multiplication/reduction modulo monic E, b*u<=k+u<=2n
PV squarings and Horner evaluation suffice. At most B condenser calls,
B hashes and polynomially many integer operations, with four-product tower
arithmetic, give O(n^8) BIT-TIME for parameters, setup and evaluation.
No field-element or truth-table enumeration occurs. Reject invalid parameters
immediately after the guard.

## 6. Concrete logarithmic scaling and the two feedback checks

This illustrative schedule is unselected. The family uses proof.md (L),
Section7f's entropy-aware two-rate base with leafceil(alpha*S),
S=floor(floor(log2 actualN)/144), alpha as in7f, and proof.md D1-D4.

Let the outer length be N, let L=ceil(log N), and assume L>=h+6. For input
length n=N, set e=L and ell=128(h+1)L. Then a<=2L, rho<=5L,
S<=50(h+1)L, and (2)-(3) imply the deliberately loose bounds

    output d=128·2^h(h+1)L,
    seed D<=1200(h+1)^2 L,
    entropy K_0<=3d.                                        (6)

For sufficiently large N the admissibility K_0<=N holds. The error is <=2^-L.
The output/seed ratio is >=

    (8/75)·2^h/(h+1),

which tends to infinity. Unlike small-bias hashing, this supplies unbounded
output/seed ratio; K_0 is absorbed in the logarithmic entropy coefficient.

**Alternating-extraction cycle.** Applying the small-bias primitive of Section
8 to both E:d x dy -> dx and Lm:r x dx ->dy requires dy>2dx and dx>2dy,
an immediate contradiction. Formula (2), in contrast, has leading seed
coefficient 6/B on its output: for output m after rounding,

    D<=6m/B+18+(h+7)S+18a.

Choose fixed B above the tampering degree, then entropies and lengths covering
the additive terms. This resolves the cycle but does not supply d,dx,dy,r,K_C
or the ordinary-breaker contract; all affine inequalities still require checking.

**Sampler/error feedback.** A sampler made by enumerating a seed D has at
least 2^D base rows. For the small-bias construction, D>2d, so after any
row expansion A>=2^D the final target e_F=1/(800 C A^9) has
log(1/e_F)>18d. That cannot be absorbed in an ordinary breaker whose
available d-bit inputs require error 2^(-c d) with a small fixed c. It is
not legitimate to 'increase d' while keeping that sampler relation fixed.
In contrast, (6) makes d/log(base row count) arbitrarily large.

If a proved amplifier has log A<=pD and degree C polynomial in D/L,
with p absolute or polynomial in h, then pD/L,C are polynomial in h
under (6), while d/L grows exponentially in h. Thus O(poly(C)*log A)
seed demand fits for large fixed h if its constants are h-independent;
(log A)^2 would not fit d=Theta(L).
This only reduces compatibility to finite graph/breaker contracts,
a fixed h and the final integer schedule; it supplies no final error bound.

## 7. A generalized finite recursive primitive (Theorem D3)

This section extends Sections 2-5; the theorem of Section 1 is its case j=2,
R=4h, U1=4n with powers of two. TheselectedbaseusesSection7f;itsentropy-awareinitialdegreeandtwofixedratesareprovedthere. SelectedNBExt0andExt2usePXTofSection7g,
atj2. Selected primitives use PXT/PXG of7g; LemmaP is unselected. lg(x) is the least integer i>=0 with 2^i>=x.

Parameters: positive integers n, e, ell; an integer j>=2 that is a power of
two; integers R>=max(4,j) and h'>=1; an integer U1>=1; a dimension rule,
either "powers of two" or "admissible" (integers r0*2^i with r0 in
{1,3,5,7}, i>=0). adm(x) is the least allowed dimension >=x. Put

    h=h' log2 j,  B=2^h=j^h',  a=e+h+4,  U=4n,
    rho=a+2lg(U)+2,  s0=adm(2rho),  rho1=a+2lg(U1)+2,  S=adm((R+1)rho1),
    K_h'=ell+2a,  K_(i-1)=ceil[(jR/(R-j+1))(K_i+a+(j-1)(S/j+1))],  i=h',...,1.

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

    D < s0+h'S+gamma*m_h',  m_h'<=kappa_j(ell+3a+(j-1)(S/j+1)+1)+S/j+1,   (2')
    K_0 <= B*Lambda_j*(ell+(2+f_j)a+f_j(j-1)S/j+f_j*j).                  (3')

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
dimension is <= gamma*max(2rho,(R+1)rho1,m_h'), and each extension tower
of Section 3 has dimension s*l<=2max(U,U1)+2s.

**Corollary D3'** (the envelope behind Lemma P). Take j=2, U1=U=4u,
height h, R=R'>=4, admissible dimensions, leaf ell=ceil(z/B),
b=e+h+4, rho_u=b+2lg(4u)+2, S'_u=adm((R'+1)rho_u), Lambda'=(R'/(R'-1))^h,
kappa_2=(R'+1)/(R'-1), gamma=5/4. (2K_0+s0<=U1 holds when K_0<=u, since
s0<=S'_u<=K_0.) Then, on input length u and output z,

    K(u,z) <= Lambda'z+G'_u,           G'_u=Lambda'B(4b+S'_u+5),
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

    Lam=(R/(R-1))^h,  kap=(R+1)/(R-1),  c=gamma*kap/2^h,  g=Lam*2^h*(4+gamma(R+1)),
    pi=(h+1)gamma(R+1)+gamma*kap*(3+gamma(R+1)/2)+gamma*(gamma(R+1)/2+1/2).

For every integer w>=max(2, e+h+6+2lg(4u)):

    K(u,z) <= Lam*z+g*w,     D(u,z) < c*z+pi*w.                        (P)

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

The selected outer primitives are PXT/PXG, with the four fixed ordered vectors
and complete rational ledger in proof.md outer Section3. Their constants
are those of Section7g below, without LemmaP's admissible rounding factor.

Legacy slot convention (restricted-nonmalleability-and-breaking.md Section 8.1; the
selected 14-copy F# certificate and the outer ledger state their own exact
rows and do not rely on this convention). Every
slot is an integer X=A_X*w+ceil(B_X*z), A_X a positive integer, B_X a
nonnegative multiple of 2^-20, w the construction's unit and z its requested
output. An inequality X>=sum_i c_i*Y_i+x_a*w+x_b*z+x_0 (all c_i,x>=0) follows
from B_X>=sum c_i*B_i+x_b with (A_X-sum c_i*A_i-x_a)w>=sum c_i+x_0,
because each Y_i<=A_i*w+B_i*z+1. A seed requirement X>=D(u,Y) is the case
c=c(h,R), x_a=pi(h,R) of (P); a threshold K(u,Y) is bounded by Lam*Y+g*w. Each
ledger line lists [sum c_i*A_i+x_a]w+[sum c_i*B_i+x_b]z+[sum c_i+x_0] (decimals
rounded up, and never above the slot's coefficients) and the slot; the stated
"valid for w>=" is an integer at least the largest [sum c_i+x_0]/(A_X-[sum
c_i*A_i+x_a]) over the lines, computed from the printed brackets, so every line
holds as printed for every unit at least that integer. The requirement q=e is charged as 1*w
(e<=w), and C's slots enter as Y_i.

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
give the required condenser modulus. Larger extensions are thus
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

    K(u,z)<=Lam*z+g*w,   D(u,z)<=c*z+pi*w.                (PX)

Proof. rho=e+h+6+2lg(4u)<=w,b=e+h+4,w>=b+6,S=(R+1)rho.
Bound(3') gives K<=Lam*z+Lam*B*(4b+S+5);4b+5<=4w gives the result.
For the seed m_h<=kap(ell+3b+S/2+2)+S/2+1. Substitute
ell<=z/B+1,3b+3<=3w,1<=w/2 into D=2rho+hS+m_h.
The same recurrence gives K0>=S>=s0, so K0<=u implies
2K0+s0<=3u<4u. No hidden extra guard is introduced.

Exact field enumeration is not generally polynomial in component dimensions.
For CA#/F#/NBExt0,Ext2, local input/output/error are O(log actualN) by
closed slots/B9, and (2') bounds base/hash degrees by O(log N).
For BASE, input is N, e=floor(floor(log2 N)/144), ell=ceil(alpha*e)
with7f's alpha, U1=2^(128+lg(e))=O(log N), and rho,rho1,S and final
block dimension are O(log N). Its O(N) initial relative extension uses
Section2's quadratic towers, never 2^O(N) enumeration. Outer PXT/PXG
also read N-bit sources but their fixed logarithmic output slots give
O(log N) rho and final-block bounds. Thus every selected enumeration
costs2^O(log N), with larger extensions built by polynomial-time towers.
No field advice or source/property search occurs.

## 7f. Rounded-leaf entropy-aware base with a smaller sampler unit

For actual N, put L=floor(log2 N), S=floor(L/144)>=2^40 and fix
 levels=26, alpha=384290345147091/93824902758400, Ro=50981815897738182537983470459/6296488643820466995200000000, Rl=6386477/1000000, initial rate R0=1099511627776.
Put ell=ceil(alpha*S), a=S+64>=S+52+4, U=4N,
U1=2^(128+lgS), rho1=a+2lgU1+2=S+322+2lgS,
So=ceil((Ro+1)rho1), Sl=ceil((Rl+1)rho1),
lo=4Ro/(Ro-3), ll=4Rl/(Rl-3).
Set K26=ell+2a, K25=ceil[ll*(K26+a+3(Sl/4+1))],
Ki-1=ceil[lo*(Ki+a+3(So/4+1))] for i=25,...,1.
After the thresholds set rho0=a+lgU+lgK0+2,
s0=ceil((R0+1)rho0/R0), u0=s0-rho0. Check K0<=N, R0<=rho0,
(R0+1)K0+s0<=U1; otherwise return zero. Condense at K0, degree s0, using7a's field/tower. Then u0>=rho0/R0>=1,
b=ceil(K0/u0)<=K0 and s0/u0<=R0+1. Expansion count <=2UK0*2^u0, while 2^rho0>=4UK0*2^a, so defect<=2^-a/2 and
output<=(R0+1)K0+s0. Section3 proves strongness/linearity using only s0>rho0.

Use 25 four-way levels at So and one at Sl. For final block length
m26=ceil(ceil(K25/(Sl-rho1))*Sl/4), put rH=lg(m26*ell+1)+ell+2a.
Extract ell bits from each of B=4^26 leaves with the weak-seed Section1
AGHP matrix and ONE fresh shared 2rH-bit seed. A repaired conditional leaf has
entropy K26=ell+2a, so the map is strong with error <2^-a and is linear.
Condition on prior seeds; replace outputs last to first. The repaired
previous-leaf prefix is independent of the fresh seed; retain it in the
single-leaf bound, then project to its outputs. Each replacement costs
<2^-a. Retaining every prior seed and both AGHP field elements costs
<B2^-a. Adding initial, condenser and split errors gives
[1+4(B-1)/3+B]2^-a<3B2^-a<2^-S, the required strong guarantee.
The integer ratios are <=(Ro+1)/Ro,(Rl+1)/Rl, and
1-3(Ro+1)/(4Ro)=1/lo, 1-3(Rl+1)/(4Rl)=1/ll. Reverse thresholds
prove each split guard. Later blocks have length <=K0, so the initial
padding suffices. Section4 supplies the causal split repairs used above. Define eps=2^-20,
spo=(Ro+1)(1+eps)+eps, spl=(Rl+1)(1+eps)+eps,
kl=ll[alpha+eps+3(1+eps)+(3/4)spl+eps]+eps, f=lo/(lo-1).
Then K25<=kl*S and Sl-rho1>=Rl*S. Exact rational substitution gives
 15Rl>=kl, hence b_last=ceil(K25/(Sl-rho1))<= 15.
The full seed is D=s0+25So+Sl+2rH. The fixed constants satisfy
alpha+eps<=16 and (b_last/4)spl+eps<=64, so ell<=16S and m26<=64S.
Thus lg(m26*ell+1)<=11+2lgS<=eps*S: on [2^k,2^(k+1)), k>=40,
13+2k<=2^(k-20), by the endpoint and doubling. Thus
2rH<=(2alpha+4+8eps)S. Reverse thresholds and output length are unchanged. Since rH=O(log N),
weak-seed Section4's field search and existing towers are polynomial in N.

All leaf/reverse ceilings are paid by eps. Summing the geometric series
in the reverse recurrence gives
 K0<=kfac*S, kfac=lo^25[kl+f((1+eps)+(3/4)spo+eps)+eps],
 ceilkfac=12484765614829834828034<2^75.
Indeed a,rho1<=(1+eps)S; the last recurrence gives kl, and each preceding
step has additive cost lo*((1+eps)+(3/4)spo+eps)*S+1. Their geometric
sum, divided by lo^25, is bounded by the f term plus eps
because 1/((lo-1)S)<=eps. Thresholds depend on rho1, not rho0.
All base/hash degrees are O(log actual N); larger relative extensions use
explicit towers. Sampler5.3 proves seed/output/padding/invocation bounds.

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
paid by the additional coefficient timesw/2^20.
These bounds prove all seed/threshold/padding guards. F#/NB dimensions are
O(logactualN); outer L0,Lr read actualN bits but fixed logarithmic slots
keep rho and final blocks O(logactualN). Base/hash enumeration is 2^O(logN),
with towers for larger extensions: polynomial in actualN, not jointly in
arbitrary component dimensions.

## 7g. Exact finite ordered-rate envelopes (PXT)

Fix h,u,z,e>=1 and a rational vector (R1,...,Rh)>1 before input length.
PXT uses a=e+h+4, rho=a+2lg(4u)+2, sinit=2rho, ell=ceil(z/2^h),
Si=ceil((Ri+1)rho), lambda_i=2Ri/(Ri-1), Kh=ell+2a and
 K_(i-1)=ceil(lambda_i*(K_i+a+S_i/2+1)), i=h,...,1.
Return zero if K0>u. Otherwise condense at entropy K0, degree sinit,
ambient 4u, then use one fresh shared Si-bit seed at binary level i:
v_i=ceil(K_(i-1)/(Si-rho))*Si, m_i=ceil(v_i/2). Hash all 2^h final leaves
with one shared exact-degree mh-bit field seed; truncate to z. Retain all
seeds and padded suffixes. PXT[h,R,b] means h-1 rates 1+b(R-1), then R
(just R if h=1); h,R,b>0 are fixed.
Section3's expansion needs only Si>rho and 2^rho>=4(4u)^2*2^a. It gives
mi<=(Ri+1)Ki-1/(2Ri)+Si/2+1, hence Ki-1>=mi+Ki+a by the recurrence.
Section4's binary split and causal shared-seed repair thus apply at each level. Also Ki-1>Si>=sinit and Ki-1>Ki: initial output
<=2K0+sinit<=3u<4u; later blocks have length <=K0, proving every padding
and entropy guard. The leaf collision proof and union bound give strong
error <3*2^h*2^-a<2^-e, retaining all seeds; fixed-seed linearity composes.

For integer w>=max(2^20,e+h+6+2lg(4u)), let eps=2^-20 and define
 Ah=2, Bh=2^-h, Ch=1;
 A_(i-1)=lambda_i*(A_i + 1 + (R_i+1)/2), B_(i-1)=lambda_i*B_i,
 C_(i-1)=lambda_i*(C_i+3/2)+1.
Since a,rho<=w and Si<=(Ri+1)w+1, descending induction in the exact
recurrence gives Ki<=Ai*w+Bi*z+Ci. Hence with Lam=B0 and g=A0+eps*C0,
 K0<=Lam*z+g*w. These A0,B0 are envelope coefficients, not sampler B0.
For the seed put R=Rh, kap=(R+1)/(R-1), tau=(R+1)/(2R),
 c=kap/2^h,
 pi=2+sum_i(Ri+1)+kap*(3+(R+1)/2)+(R+1)/2
       +eps*(h+5kap/2+tau+3/2).
Indeed the last actual block satisfies
 mh<=tau*Kh-1+Sh/2+1
   <=kap*(ell+3a+Sh/2+1)+tau+Sh/2+1.
Substituting ell<=z/2^h+1 and Sh<=(R+1)w+1 bounds mh by c*z plus
[ kap*(3+(R+1)/2)+(R+1)/2 ]w+5kap/2+tau+3/2.
Adding sinit+sum Si proves D<=c*z+pi*w, charging every internal ceiling.
Sections7a,3,4,7e supply canonical fields, towers and retention. Selected
u<=n, z,e=O(log n) imply O(log n) base/hash degrees and 2^O(log n)
setup; larger extensions use polynomial-time towers. No field/vector advice.

**PXG terminal variant.** PXG[h,R,b] keeps PXT's vector, thresholds and
condensers; replace the leaf hash by weak-seed Section1's AGHP matrix on
mh bits, output ell, exponent a, fresh shared seed of length
2rH=2(lg(mh*ell+1)+ell+2a). The7f reverse hybrid with B=2^h retains all
prior seeds and repaired leaf prefixes, costing <B2^-a. Earlier costs are
[1+2(B-1)]2^-a, so total <(3B-1)2^-a<2^-e. Fixed-seed linearity composes;
K0,Lam,g and all padding/entropy guards are unchanged. Retain every new
seed bit and suffix: this defines a new algorithm, not seed deletion.
For the seed envelope require 5/3<=Rh<=64 and z<=256*2^h*w. Then
ell<=256w, kap<=4, tau<=4/5, Sh<=65w+1; the mh bound gives
mh<=1199w+9<=2048w. Thus lg(mh*ell+1)<=20+2lgw<=w/2^14:
on [2^k,2^(k+1)), k>=20, use 22+2k<=2^(k-14), true at20 and
preserved by doubling. Adding all seed slices gives
 D=2rho+sum Si+2rH<=c*z+pi*w,
 c=2/2^h, pi=6+sum_i(Ri+1)+2^-13+(h+2)eps.
Here a,rho<=w; ell<=z/2^h+1 pays two ceiling bits and Si pay h.
AGHP's proof and canonical fields are below. At selected calls
rH=O(log n), hence setup/evaluation are uniform polynomial time.
PXT remains available with its original hash.

Design credit: mpjunior92, submission 8ef2ef8 (PEER-NOTE-8ef2ef8.md):
vector envelopes, capped unit and charged slots in wrappers4.2.
This material dependency is reproved here; its bound/status is no premise.

## 10. Attribution and source convention

The PV-code expansion/block method is due to Guruswami, Umans and Vadhan:
https://people.seas.harvard.edu/~salil/research/PV-condenser-June08.pdf .
The recursive linear condensation/shared-seed strategy follows
Chattopadhyay, Goodman and Liao, ECCC TR21-075, Lemmas 4.7--4.10 and
Theorem 5.6: https://eccc.weizmann.ac.il/report/2021/075/download/ .
Chattopadhyay--Liao use a sharper variable-parameter form in Lemma 3.13:
https://arxiv.org/pdf/2110.12652v1 . Xin Li's later application is
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/ .

This independently proved weaker construction excludes retained seed bits
from its strong output. It imports no sharper or unrestricted-error theorem.

# An elementary weak-seed replacement for the local Raz call

This local interface uses Alon--Goldreich--Hastad--Peralta (1992) small bias
and the leftover-hash collision calculation. Its high-rate-row/low-rate-source
use is reconstructed here, without asserting the full advice condenser,
breaker, affine extractor or circuit theorem.

Sources inspected:
https://www.cs.tau.ac.il/~nogaa/PDFS/aghp4.pdf
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/
(Algorithm 3 and Lemma 5.1).

Logs are binary, TV half L1, vector/matrix arithmetic over F_2.
Section 4 proves finite-field existence, deterministic setup and outer-length cost.

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

  TV((E(X,U_{2r}),U_{2r}),(U_m,U_{2r})) <= 2^-e.              (1)

Proof of small bias. A nonempty parity of the b_i equals the dot product of
bin(P(u)) with bin(v), where P is a nonzero polynomial with binary
coefficients and degree at most N-1. If P(u) is nonzero, the dot product is
balanced over uniform v; if P(u)=0, the sign equals one. A nonzero polynomial
has at most its degree many roots, by repeated division by X-a. Thus the
absolute bias is <=

  delta=(N-1)/2^r < 2^(-m-2e).

For every nonzero w in F_2^d and every nonzero a in F_2^m, the parity
<a,Mw> uses a nonempty collection of matrix entries: choose any nonzero
coordinate of a and of w. Thus the same bias bound applies.
Fourier inversion of the zero event gives

  Pr[Mw=0]=2^-m*sum_a E[(-1)^(<a,Mw>)] <= 2^-m+delta.        (2)

Let X' be an independent copy of X. For each seed s let p_s be the law of
E(X,s). Averaging its collision probability and separating X=X' yields

  E_s sum_y p_s(y)^2
  <= Pr[X=X']+(1-Pr[X=X'])*(2^-m+delta)
  <= 2^-m+2^(-m-2e)+delta.

Cauchy--Schwarz gives TV(p_s,U_m)<=sqrt(2^m*sum_y(p_s(y)-2^-m)^2)/2.
Jensen's inequality for the square root now bounds the average TV by

  (1/2)*sqrt(2^m*(2^(-m-2e)+delta)) < 2^-e.

That average is exactly the joint TV in (1), proving the strong guarantee.
No statement about an unexposed-seed extractor is being substituted for it.

## 2. A weak seed whose entropy deficit is paid explicitly

Suppose a seeded map uses q seed bits, is strong with error eta for uniform
seed, and ignores any prescribed unused seed suffix. Let S be any independent
q-bit seed with H_infinity(S)>=q-Delta. Then

  TV((E(X,S),S),(U_m,S)) <= 2^Delta*eta.                     (3)

For each s the conditional distance is a number f(s) in [0,1]. Its mean
under uniform seed is <= eta. Since Pr[S=s]<=2^Delta/2^q,
its mean under S is <= 2^Delta*eta. This also proves the claim with the
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
roots. We have thus constructed a field F with exactly 2^r elements.

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
length exactly j by the preceding minimal-polynomial argument. Thus
Q divides X^(2^h)-X exactly when j divides h. The displayed congruence test
makes every factor degree divide r and excludes repeated factors, since
X^(2^r)-X has derivative one. If a factor degree j were a proper divisor
of r, choose a prime p dividing r/j; then j divides r/p, and the relevant
gcd test would detect Q. The tests thus force a single degree-r
irreducible factor, proving their correctness.

There are 2^r candidates, and each test uses polynomially many r-bit
polynomial operations. Trial division factors r in polynomial time in r.
The search is thus 2^r*poly(r)=n^O(1). This argument specifically relies
on r=O(log n). It does NOT make an r=Theta(n) field search polynomial.
The matrix can be streamed using N=d*m field powers and row sums, so no
exponential table of seed values is stored or searched.

The separate seeded-extractor worker is deriving an elementary tower setup
for its larger-field recursive construction. That work is not a premise of
this local logarithmic-degree implementation.



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
    <= (1/2)2^((h_b+m_b-k_1-k_2)/2).                        (IP1)

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
polynomial in the actual outer n. It is not a claim that exhaustive setup is
polynomial in h_b for arbitrary h_b. Small lengths use the same finite search.

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
irreducible. Conversely every irreducible degree-h_b P passes those tests.
This is the usual elementary Frobenius characterization of binary finite
fields. Existence of an irreducible polynomial of every positive degree is the
standard finite-field existence theorem, not a nonuniform choice of a table.

Schoolbook multiplication and reduction of degree-below-h_b binary polynomials
cost O(h_b^2) bit operations. A modular Frobenius test uses at most h_b
squarings; the gcd tests use polynomially many operations on h_b-bit
polynomials. A deliberately loose upper bound per candidate is O(h_b^5) bit
operations, including factoring h_b by trial division. There are 2^h_b
candidates, hence setup costs at most O(2^h_b h_b^5) bit operations and
polynomial working space. Once P is fixed, a call to IP_b costs O(h_b^2) bit
operations and O(h_b) temporary storage; prefixing the result costs O(m_b).
Caching the polynomial is optional and is not part of the mathematical family.

For the intended local use h_b<=C log_2 n with one fixed C, where n is the
actual outer input length. Thus setup is <= n^C times a fixed
polylogarithmic factor and evaluation is polynomial in n. If h_b is not
O(log n), this exhaustive setup bound cannot be called polynomial in n or
in h_b. A different deterministic field-construction algorithm would be
needed for a general polynomial-in-h_b claim. All finite lengths remain
well-defined by the same search, with no uncomputable exceptional branch.

## Attribution, application, and limits

Field bilinear extraction and its character/collision analysis are standard
inner-product hashing. The needed contract is the one cited by Xin Li,
*Two Source Extractors for Asymptotically Optimal Entropy, and (Many) More*,
ECCC TR23-023 revision 1, Theorem 2.20:
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/
The derivation improves the constant used in condenser-integration.md;
it claims no new extractor theorem.

The hypotheses include independence of the two honest sources, integer
1<=m_b<=h_b, and entropy lower bounds 0<=k_1,k_2<=h_b. When the error
exceeds one it can be capped at one, but it provides no useful extraction in
that range. No claim retains both complete inputs jointly: the output is a
deterministic function of those inputs.

The remaining ordinary/advice condenser and CZ proofs are supplied in the
companion modules; the final assembly is in proof.md.
