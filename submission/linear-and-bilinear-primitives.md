# linear-and-bilinear-primitives.md

Component references resolve through proof.md Section 0.

# Finite seeded and linear extraction: an elementary logarithmic-regime replacement

All logarithms are binary; distance means half-L1.

## 1. Main result and what it replaces

An elementary small-bias hash family alone is **not** enough: its seed exceeds
its output, which creates an impossible cycle in alternating extraction and
also fails the sampler/error feedback. A balanced condenser-and-block-hash
construction below removes that particular obstruction.

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

If K_0<=n, there is one fully specified deterministic strong LINEAR extractor

    H: {0,1}^n x {0,1}^D -> {0,1}^(B ell)

for min-entropy K_0 and error at most 2^(-e). Its exact seed length is
s0+hS+W, where W is defined in Section 5, and the useful finite bounds are

    D <= 6ell+(h+7)S+18a+12,                                  (2)
    K_0 <= 2B(ell+4a+S+4).                                   (3)

For an arbitrary target m, take ell=ceil(m/B) and discard extra output bits.
Every operation is specified using binary arithmetic and quadratic finite-field
towers. Its bit-time is polynomial in n, h, e and B; the outer application fixes
h and has e=O(log outer N), so its evaluator is polynomial in outer N. There is
no truth-table search, randomized setup, advice table, or unproved extractor
primitive in this theorem. Elementary field and probability facts used below
are proved to the extent needed by the algorithm.

The construction is a deliberately weaker, finite-constant reconstruction of
the linear-condenser/block-source approach in GUV and CGL. It does NOT assert
CL Lemma 3.13 with that lemma's precise universal entropy constant or variable-t
tradeoff. It is enough to obtain arbitrarily large fixed output/seed ratios
with logarithmic outer entropy. Whether a given application can use the weaker
finite threshold (1) must be checked, not inferred from big-O notation.

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
                       (f^(H^(b-1)) mod E)(y)).             (4)

The output has v bits. Each map f -> f^(H^j) mod E is F_2-linear, because
H is a power of 2. Evaluation and padding are linear too. Moreover

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
admissible input. It is an intentionally loose single polynomial bound.
No step enumerates a field's elements or extractor truth tables. Invalid
parameter inputs can be rejected immediately after the integer guard.

## 6. Concrete logarithmic scaling and the two feedback checks

Let the outer length be N, let L=ceil(log N), and assume L>=h+6. For input
length n=N, set e=L and ell=128(h+1)L. Then a<=2L, rho<=5L,
S<=50(h+1)L, and (2)-(3) imply the deliberately loose bounds

    output d=128·2^h(h+1)L,
    seed D<=1200(h+1)^2 L,
    entropy K_0<=3d.                                        (6)

For sufficiently large N the admissibility K_0<=N holds. The error is <=2^-L.
The output/seed ratio is at least

    (8/75)·2^h/(h+1),

which tends to infinity. This is the precise property missing from the
simple small-bias construction. The constants here have not been optimized.
The output need not be near the full source entropy for a sampler; K_0 can
be absorbed in the logarithmic entropy coefficient.

**Alternating-extraction cycle.** Applying the small-bias primitive of Section
8 to both E:d x dy -> dx and Lm:r x dx ->dy requires dy>2dx and dx>2dy,
an immediate contradiction. Formula (2), in contrast, has leading seed
coefficient 6/B on its output: for output m after rounding,

    D<=6m/B+18+(h+7)S+18a.

Choose fixed B sufficiently larger than the fixed tampering degree and then
choose the available source entropies and lengths to dominate the explicit
additive terms. The former mandatory expansion cycle disappears. This is
NOT yet a full choice of d,dx,dy,r,K_C: the ordinary breaker and its entropy
contract remain to be supplied, and every inequality in the affine module
must still be tested together. In particular a bound O(log N) by itself
never licenses any of these calls.

**Sampler/error feedback.** A sampler made by enumerating a seed D has at
least 2^D base rows. For the small-bias construction, D>2d, so after any
row expansion A>=2^D the final target e_F=1/(800 C A^9) has
log(1/e_F)>18d. That cannot be absorbed in an ordinary breaker whose
available d-bit inputs require error 2^(-c d) with a small fixed c. It is
not legitimate to 'increase d' while keeping that sampler relation fixed.
In contrast, (6) makes d/log(base row count) arbitrarily large.

If the separately proved graph amplifier has log A<=pD and degree
C polynomial in D/L, with p absolute or polynomial in h, then both pD/L
and C are polynomial in h under (6), while d/L grows exponentially in h.
Hence a breaker seed requirement O(poly(C) * log A) fits for a sufficiently
large FIXED h, provided its constants are independent of h. An arbitrary
polynomial in log A would not suffice: (log A)^2 introduces L^2, which no
fixed h can absorb into d=Theta(L).
This is a useful compatibility reduction, not a proof that the missing
graph/breaker contracts actually have the required quantified bounds.
One must supply those finite functions, then choose h and verify the final
integer schedule. The promised final error is not supplied by this module.

## 7. A generalized finite recursive primitive (Theorem D3)

This section extends Sections 2-5; the theorem of Section 1 is its case j=2,
R=4h, U1=4n with powers of two, and the near-full, NB and inner-breaker uses
keep that case. lg(x) is the least integer i>=0 with 2^i>=x.

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
Consecutive admissible integers a<b satisfy b<=5a/4 (pattern 4,5,6,7,8 times
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

**Corollary D3'** (the envelope of proof.md Lemma A'). Take j=2, U1=U=4u,
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

## 7a. Exact binary degrees in the base sampler only

Define D3_exact with the same parameters, recurrence, condenser formulas,
j-way splits, exactK0 guard andU1 padding guard asTheoremD3, but choose
s0=2rho,S=(R+1)rho1,W=m_h' exactly. This variant is used ONLY for the base
sampler. Every breaker,NB andconstant-advice call retainsoriginalH, and
CorollaryD3' is not substituted into those calls.

For every positive binarydegreev choose the first monic irreducible binary
polynomial by lexicographic enumeration and the finite modular-squaring/gcd
test proved inweak-seedSection4 later inthisfile. Its independentalgebraic
proof shows existence andcorrectness at everydegree. The quotientfield
uses polynomialbasis coordinates andschoolbookarithmetic, withsetup2^v poly(v).
Absolute trace is a nonzero binarylinear map: its polynomialsumX^(2^j),
j=0..v-1,is nonzero ofdegree2^(v-1)<2^v, so it cannot vanish everywhere.
Test the orderedbasisunits tofind the first trace1unit; linearity ensures
one exists. Startingover this arbitrarydegree-s basefield, build Section2's
Artin--Schreiertower. The same trace-transitivity andimage/kernel proof
makes eachquadraticstepirreducible andpropagates trace1. Only the relative
extensiondegree l needs tobe a poweroftwo. The lastroot's q-Frobeniusorbit
has exactlyl elements by the same proper-divisorargument, so its orbitproduct
constructs the required irreducible condensermodulus overF_(2^s). Embedded
basecoordinates remain the firsts positions. This does not enumerate2^(sl)
objects. The finalhash uses the quotientfieldofdegreeW=m_h', withidentical
collision andsharedseedproofs, withoutpaddingblocks to a dyadicfielddegree.

D3's expansion, splitting, conditionalcoupling, error andlinearity arguments
use s0>=2rho,S>=(R+1)rho1 andfieldoperations, so remainvalid. Its same
recurrence implies bound(3'), and its finalsplitboundm_h' in(2') remainsvalid.
The exactseed is nowD=s0+h'S+m_h', so the seedupperbound(2') is read as
D<=s0+h'S+m_h' withgamma=1, not a strictinequality. In the application below,
strict slack beforethe next evenintegerB0 provesD<B0L.

Runtime: no polynomial inarbitrarycomponentdimensions isasserted for this
variant. InitsONLY useh',R,j arefixed,ell=17L,e=L,U1=2^(132+lgL),rho=O(L),
rho1=O(L),andm_h'=O(L) by(2'). All enumeratedbinarydegrees areO(L),
L=lg(N+1) foractualouterN. Their setupcost2^O(L) is a fixedpolynomial inN.
Initialcondenser extensions mayhavedimensionO(N); these are built by the
quadratic tower inpolynomial time, never enumerated atcost2^O(N). Remaining
arithmetic,evaluation,numberofblocks andallconstantparameters givefixed
polynomial time inactualN. Thereisnosource/property/truth-table search,
randomsetup oradvice. Invalidcalls stop at the exactguards.

## 8. The simpler small-bias hash: exact proof and limitations

For source length n, output m, and e>=1, let L=nm and
r=ceil(log max(1,L))+m+2e. Construct F_(2^r), for this section either by
exhaustive search over monic degree-r binary polynomials with polynomial-time
irreducibility testing, or another independently supplied deterministic
algorithm. The exhaustive search costs 2^r poly(r), polynomial in outer N
when n<=N^O(1), m,e=O(log N). It is not polynomial in log N alone.
For seed (u,v) uniform in the square of that field, define L binary entries

    b_j=<v,u^j>, j=0,...,L-1,

where the inner product is in the fixed binary coefficient basis. For any
nonzero binary c, sum c_j u^j is a nonzero polynomial of degree at most L-1.
Except at at most L-1 field elements u, averaging over v annihilates its
parity. Thus every nonzero parity of the entries has bias at most
(L-1)/2^r<=2^(-m-2e). Arrange entries into an m by n matrix M and output Mx.
For z!=0, Fourier expansion of the event Mz=0 gives

    Pr[Mz=0]<=2^-m+2^(-m-2e).

For a source of entropy k>=m+2e, average output collision probability is
at most 2^-m+2^-k+2^(-m-2e). Cauchy-Schwarz/Jensen gives strong distance
at most (1/2)sqrt(2^(m-k)+2^-2e)<=2^-e/sqrt(2). The seed has 2r bits.
This fully elementary construction remains useful for calls without a
short-seed/output expansion requirement. It does not close the feedback
loops discussed above.

## 10. Attribution and source convention

The PV-code expansion/block method is due to Guruswami, Umans and Vadhan:
https://people.seas.harvard.edu/~salil/research/PV-condenser-June08.pdf .
The recursive linear condensation/shared-seed strategy follows
Chattopadhyay, Goodman and Liao, ECCC TR21-075, Lemmas 4.7--4.10 and
Theorem 5.6: https://eccc.weizmann.ac.il/report/2021/075/download/ .
Chattopadhyay--Liao use a sharper variable-parameter form in Lemma 3.13:
https://arxiv.org/pdf/2110.12652v1 . Xin Li's later application is
https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/ .

This file proves its own finite weaker construction, not those sharper
statements. Its modern strong output excludes retained seed bits. No
unrestricted-error quotation of a near-full modern strong extractor is
used.
The constants, finite tower setup, recurrence, guards and bounds here are
independently derived, with no imported primitive theorem.


# An elementary weak-seed replacement for the local Raz call

This is original proof support for one local extractor interface. It does not
prove Li's advice condenser, the entire correlation breaker, the final affine
extractor, or a logarithmic circuit candidate. The small-bias construction is
credited to Alon, Goldreich, Hastad and Peralta (1992); the collision argument
is the standard leftover-hash calculation. The proposed use at this particular
high-rate-row/low-rate-source interface is an independent reconstruction.

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

  TV((E(X,U_{2r}),U_{2r}),(U_m,U_{2r})) <= 2^-e.              (1)

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
under uniform seed is at most eta. Since Pr[S=s]<=2^Delta/2^q,
its mean under S is at most 2^Delta*eta. This also proves the claim with the
ENTIRE q-bit seed retained. It does not claim that the result remains random
when both complete inputs are exposed.

Take now independent sources Y on d bits and S on q bits. Define

  m=floor(q/64), e=ceil(7q/32),
  r=ceil(log(d*m+1))+m+2e.

Assume the following exact finite guards:

  m>=1, 2r<=q, m+2e<=floor(d/10),
  H_infinity(Y)>=floor(d/10), H_infinity(S)>=4q/5.            (4)

Use the first 2r bits of S as (u,v), ignoring the remaining seed bits, and
put W=E(Y,S). Equations (1) and (3), with Delta=q/5, give

  TV((W,S),(U_m,S)) <= 2^(q/5-e) <= 2^(-3q/160).            (5)

The exponent identity is 7/32-1/5=3/160. Since m=floor(q/64), this is a
linear output length and exponentially small error as q grows.

For convenient sufficient parameter checks, put B=ceil(log(d*m+1)). Since
2m+4e <= (29/32)q+4, the condition

  q>=32B+64

implies 2r<=q. The entropy guard m+2e<=floor(d/10) should be checked directly;
for example q<=d/8 and sufficiently large d leave a strict linear margin.
If q is bounded above and below by fixed positive multiples of d, with the
upper multiple small enough, both guards hold at every sufficiently large d,
not merely a subsequence. Their eventual validity follows from log d=o(d).

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
actual outer input length. Therefore setup is at most n^C times a fixed
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
The derivation above is independently written and proves a slightly stronger
constant than the error bound used in condenser-integration.md. It is not a
claim of a new extractor theorem.

The hypotheses include independence of the two honest sources, integer
1<=m_b<=h_b, and entropy lower bounds 0<=k_1,k_2<=h_b. When the displayed error
exceeds one it can be capped at one, but it provides no useful extraction in
that range. No claim retains both complete inputs jointly: the output is a
deterministic function of those inputs.

The remaining ordinary/advice condenser and CZ proofs are supplied in the
companion modules; the final assembly is in proof.md.
