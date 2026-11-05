# Supplied extraction and robust-disperser arguments

The complete A1 strong linear seeded-extractor proof, including field and matrix existence, is supplied in seeded.md. Every reference to A1 below uses that supplied argument.

# B1. Deterministic linear rank condensation

## Fields without per-length advice

For every power-of-two integer s construct F_(2^s) by an Artin--Schreier
tower. At degree 1 the field is F_2. Given the degree-h field, adjoin T
with T^2+T+a=0, where a is its top basis vector, the product of all
previous tower generators (and a=1 when h=1). Its absolute trace is one.
Indeed Tr_(2h)/h(u+vT)=v; applying this identity successively proves the
assertion for the top basis vector. The image of x->x^2+x is exactly the
trace-zero hyperplane: its kernel has size two, its image has codimension
one, and trace annihilates it. Thus the quadratic is irreducible at each
step. The embedded degree-h field occupies the low h coordinates.
Multiplication is

 (u+vT)(r+wT)=(ur+a*vw)+(uw+vr+vw)T.

The implementation evaluates its five smaller-field products by an explicit
postorder stack. There is no unbounded interpreter recursion. The recurrence
uses polynomial time in the field degree; shifts and XORs have polynomial
bit cost.

Let s,D be powers of two with D>=2. In the degree-sD tower let theta be
the last quadratic generator. It lies outside the degree-sD/2 subfield.
Every proper intermediate field over F_(2^s) is contained in that subfield,
since its degree divides the power of two D. Hence theta has relative
degree D. Put q=2^s and

 E(Y)=product_(j=0)^(D-1) (Y+theta^(q^j)).

Frobenius permutes its roots, so its coefficients belong to the embedded
F_q. It is the irreducible minimal polynomial of theta. Computing the
product and the q-th powers requires polynomial time in sD and s. The
routine constructs E from its integer dimensions; no input prefix of
length comparable with the ambient input is hidden here.

## The polynomial graph and its expansion

Fix a field F_q, an irreducible E of degree D>=2, and h=2^a. For a
polynomial f of degree below D, define f_i=f^(h^i) mod E. For a seed y
in F_q output the m symbols f_0(y),...,f_(m-1)(y). Each fixed-seed map
is F_2-linear. View the seed and output together as a graph neighbor
Gamma(f,y)=(y,f_0(y),...,f_(m-1)(y)). Let

 A=q-(D-1)(h-1)m.

For every integer J<=h^m and every set S of J polynomials, its union of
neighbors has size at least A*J. Here and below a vacuous nonpositive A
causes no difficulty, although our parameter choices make A positive.

To prove this, suppose a set T of fewer than A*J neighbors contains all
neighbors of at least J polynomials. Interpolate a nonzero polynomial

 Q(Y,Z_1,...,Z_m)=sum_(j=0)^(J-1) p_j(Y) M_j(Z_1,...,Z_m),

where deg p_j<A and M_j has the base-h digits of j as exponents. There
are A*J unknown coefficients and fewer constraints. Among solutions choose
one of minimal Y-degree. Not all p_j are divisible by E: otherwise divide
Q by E; E has no F_q root, so the quotient still vanishes on T and has
smaller Y-degree. For a polynomial f whose neighbors all lie in T,
Q(Y,f_0(Y),...,f_(m-1)(Y)) has degree at most
(A-1)+(D-1)(h-1)m<q and vanishes at all q field points. It is identically
zero. Reducing modulo E, f is a root over the extension field of

 Q*(Z)=sum_(j=0)^(J-1) (p_j mod E)*Z^j.

This is nonzero and has degree at most J-1. It cannot have J distinct
roots. This contradiction proves expansion. Interpolation is used only in
the proof, not as an algorithm that enumerates source points.

## Exact parameters and rank preservation

For n input bits, an integer entropy threshold K>=1, and error log e,
let

 a=ceil(log2 K)+e+10,
 s=the least power of two >=ceil(log2 n)+ceil(log2 K)+a+e+3,
 D=the least power of two >=max(2,ceil(n/s)),
 m=ceil(K/a), h=2^a.

Pack n input bits injectively into D field coefficients. The deterministic tower and relative E are unchanged. Since D-1<=4n and m<=2K for positive integers,

 (D-1)(h-1)m<=8nK*2^a
 <=2^(ceil(log2 n)+ceil(log2 K)+a+3)<=2^(-e)*q.

Thus A>=(1-2^-e)q exactly as required. Also s>=a, so J=2^K<=h^m remains available and the output has enough bits for rank K. Every parameter is computed by integer arithmetic. The scalar audit separately tests408 positive plans against this inequality.

The new seed bound is s<2ceil(log2 n)+4a, because the lower rounding target is ceil(log2 n)+2a-7. Thus exact output nprime=m*s satisfies

 nprime<=4K+2K*ceil(log2 n)/a+2ceil(log2 n)+4a.          (B1.1)

This is NOT claimed O(K+log n+e) for arbitrary parameters. In each actual application here e>=c*log n eventually for a fixed absolute c>0: initial precision has c=1/28 up to constant slack, and internal precision e grows faster than log n. Hence ceil(log2 n)/a is bounded by a fixed constant, making nprime=O(K+log n+e) in these uses. On all other admitted calls the displayed general finite-parameter bound remains polynomial. No ambient-length term is omitted from source or time accounting.

For an affine input of dimension exactly K, each seed output is an affine
image of rank r_y<=K. Its graph-neighbor count is sum_y 2^r_y. Expansion
therefore gives E_y[2^(r_y-K)]>=1-2^-e. Whenever r_y<K this summand is at
most 1/2, so the fraction of rank-deficient seeds is at most 2^(1-e).
On every other seed the output is an affine source of dimension K. For
larger input dimension decompose it into parallel K-dimensional flats and
average; alternatively the image rank bound follows from any such flat.
This is the precise affine rank-condensation statement used below.

# B2. Supplied strong linear extraction routines

An affine distribution has a vector-space complement description: if Z is
an affine function of X, then X=A+B, with A an affine function of Z and
B an independent uniform linear source. Its dimension is H(X|Z).
Conditioning on a fixed linear observation of t bits reduces any remaining
affine rank by at most t. An observation determined by Z does not reduce
H(X|Z). These statements follow by choosing a basis of the projection
kernel and extending it to a basis of the direction space.

For any of the initial or internal maps, consider an affine extractor plan for source length n, entropy K, output
M<=K/2 and precision 2^-e. First use B1 with precision e+10. On every
good first seed its output is affine of entropy at least K. For its exact ambient length nprime from(B1.1), use A1 with the plan's stated dimension D and
per-call error 2^-(e+ceil(log2 K)+10). Choose the largest power-of-two
output m0<=min(M,K/4) satisfying A1's condition with threshold floor(K/2).
The executable plan halves m0 until it succeeds. Failed finite parameter
cases are rejected by the outer ready flag.

Repeat A1 on the SAME condensed input with fresh independent seeds until
M bits have been produced, truncating the last call. Earlier fixed seeds
and outputs leave an affine source of entropy at least K-M>=K/2. A
sequential coupling loses at most the sum of call errors; there are at
most K calls. Adding the B1 bad-seed probability gives total strong error
at most 2^-e. Each fixed combined seed is a linear map of the original
input. Padding a longer stipulated seed with ignored bits preserves this.

We prove a budget retaining the dimension explicitly. Put
H=D*ceil((n')^(1/D)), R=16D+16. A1's field degree satisfies
s=O(e+log K+log H+log R). Substituting a sufficiently small fixed
multiple of K/[R*H*(e+log K+log H+log R)] for m0 in A1's condition
shows that the largest admissible power of two is at least a fixed
multiple of that value, unless M is smaller and one call suffices.
The additive logarithmic term is o(K) in all our uses. Thus the combined
seed length is at most

 C(log n+e+log K)
 +C(D+1)*R*H*(e+log K+log H+log R)^2.                    (B2.1)

Here C is independent of the displayed dimensions. This proves eventual
success of the halving loop as well as its seed cost. Only actually
computed scalar dimensions enter prefix size or executable seed slicing.

Initial extraction keeps D=11 FIXED, source k, output Y and E0=ceil(L/28)+10. B7 now uses k=Theta_S(L^9 z^5),Y=Theta_S(L^7 z^5). Its advice Y*k^(1/11)*L has L exponent7+9/11+1<9. Therefore it is o(k) and one call eventually suffices. The generic A1 argument applies at this fixed dimension.

For D=11, r=192 and r-D=181, E0=ceil(L/28)+10. Exact A1 rounding gives s0<=(192/1267+o(1))L. The adaptive B1 block exponent is a0=(1/28+o(1))L and its field-rounding target is(15/14+o(1))L. Therefore its seed is at most(15/7+o(1))L and its output is O(k). Combining12s0 with this seed yields

 d0<=(2304/1267+15/7+o(1))L<4L.

The explicit5L cap passes eventually. Initial A1 and all B5 data are now canonical by C1. The header and four internal packages contribute L+O_S(L^2 z^2) actual bits. Thus the COMPLETE prefix P is O_S(L^2 z^2), not O(L), and remains o(k). Each asymptotic estimate is at fixed scale and no numerical threshold is claimed; C1 separately proves canonical-search runtime at every admitted header.

## B2.2 Variable-dimension internal maps with charged data

Put z=ceil(log2 L),D_int=max(2,z),C=2^30 and E=e. Each of the four internal maps is exactly B2.1 AffineSeedPlan at dimension D_int, source threshold K=C*E*z^2*M, output M and stipulated seed length ell. W/T output ell, E2 outputs q0 and E3 outputs d. Each reads its OWN C5 package of t-a(t) bits,t=D_int*s. C5 decodes a monic constant-one binary polynomial or failure; C3 reconstructs a field of size2^s and relative companion action; arbitrary data use its bounded guarded ring algorithm. Good primitive binary data exist by A1 and C3; no field/matrix search or uncharged advice is used. Arbitrary admitted data still give bounded binary-linear ring/matrix computations. Concatenation, truncation and ignored padding preserve fixed-seed linearity.

The dimension now grows; the following bounds retain it explicitly. At each fixed S, B7 gives largest K=O_S(L^6*z^8). B1 precision E+10 dominates every actual log ambient length, including log N=L for W. Its exact condensed nprime is O_S(K). As N tends to infinity, D_int=z and log2(nprime)/z has limsup at most6, because log2 L/z tends to1 and log z/z tends to0. Thus ceil(nprime^(1/D_int))<=128 eventually, uniformly across the four calls for this fixed S. Consequently H=D_int*ceil(nprime^(1/D_int))<=128z and R=16D_int+16<=32z.

The exact A1 field formula, with per-call precision E+ceil(log2 K)+10, gives s/E<=5 eventually: R/(R-D_int)<=16/15 and every other log term is o(E). The advice for the first power-of-two output candidate m0<=M is at most16*M*R*H*s<=327680*M*E*z^2. This is less than K/2 with fixed C=2^30; the additive logarithmic condition is negligible with that strict margin. Also K/4>=M. Hence the first candidate passes eventually, with at most TWO calls (one when M is a power of two). The source threshold in each call is floor(K/2), exactly as in the supplied B2.1 proof.

B1 seed is <5E eventually. The direct calls use at most2*(D_int+1)*s<=20zE bits for z>=1 in the eventual regime. Total actual seed is <25zE. Since ell/E>=S^2*z/2, any fixed S>=8 pays this absolute seed bound eventually. The source/error guarantee follows B2.1 with its full seed retained; no fixed-dimension constant is hidden in this estimate.

All four packages total O(D_int*E)=O_S(L^2*z^2) charged bits. Uniform evaluation remains polynomial even though D_int grows: on admitted plans D_int=O(loglog N), nprime,H and field degree are polynomial in N. The monomial generator visits ONLY its first nprime coefficient tuples. Each tuple has D_int entries, updates and exponentiation have polynomial cost in D_int,log H and field degree, and the outer degree loop is bounded by H. No full grid or nprime^D_int list is constructed. Matrix arithmetic costs polynomial(D_int,s). Thus the fixed polynomial-time algorithm does not acquire an input-dependent polynomial exponent.

The scalar A1 field computation uses ceil_log(2^j*R*H)=j+ceil_log(R*H), j=32+4*errorlog+4*log2(m0), for power-of-two m0. This is exactly the previous field degree and avoids allocating an integer of2^errorlog magnitude. Actual fields are bound only after all ready guards pass.

## B2.3 Multi-bit polynomial hashing for short extraction

Let n,K,M,E be positive integers with K<=n, M<=K/2 and K>=M+2E+10. Put s=P2(8*(M+E+ceil(log2 n)+10)). Use the deterministic degree-s tower field of B1. Pack the n input bits injectively into coefficients of a polynomial P_X over that field, padding the last coefficient with zeros. For a uniform seed(y,z) in F_(2^s)^2 output the first M tower-basis coordinates of z*P_X(y). Since M<=s, this coordinate projection is surjective. The exact seed length is2s, input data length is zero, and every fixed-seed map is binary linear. Unused stipulated seed bits are ignored.

Here is the full strong-extraction proof for ANY source with min-entropy at least K. Independent source copies coincide with probability at most2^-K. For distinct source points the polynomial difference is nonzero and has degree below ceil(n/s), hence at most n roots. At every other y, multiplication by its value is a bijection on the field, so uniform z makes the M-bit difference exactly uniform. Averaging over the ENTIRE seed(y,z), the conditional output collision probability is at most2^-M+2^-K+n/2^s. Cauchy--Schwarz therefore bounds statistical distance between(seed,output) and(seed,uniform M bits) by

 (1/2)*sqrt(2^M*(2^-K+n/2^s)).

The entropy margin K>=M+2E+10 bounds the first summand after multiplication by2^M by2^(-2E-10). The chosen field degree is at least M+2E+ceil(log2 n)+10, so the second summand has the same bound. The total strong error is at most(1/2)*sqrt(2^(-2E-9))<2^-E. This argument includes all seed bits, does not condition on future outputs, and requires no rank-condenser or external extractor theorem. Affine sources are included, so B3--B4's fixed-seed rank and weak-seed estimates apply.

The evaluator uses tower multiplication and Horner's rule on ceil(n/s) coefficients and then projects the first M bits. It has polynomial bit complexity in n,s,M on every admitted call. For short E1, M=ell, E=Eprime, and ambient n=N; Eprime dominates ell and log N at every fixed scale. Thus exact seed2s=O(Eprime), fitting the stipulated d once the fixed scale is large enough. K=10d meets the entropy margin by the explicit scale comparisons in B7. No field, matrix, seed table or good construction is read as input data for this map.

# B3. Parallel affine look-ahead extraction

Let at most h streams Q_i be affine functions of Z, each padded to a common
ambient bound. Let a<=h be the look-ahead time limit, separate from the
number of streams. Put z=ceil(log2 L),C=2^30,E=e, and let q0=(2h*a+C*E*z^2+10)ell. One Q_v is exactly uniform on q0 bits. Require q0-2h*(a-1)ell>=C*E*z^2*ell. Let H(X|Z) be at least h*a*ell+C*E*z^2*ell. The extractors W on X and T
on Q_i have output and seed length ell, source entropy thresholdC*E*z^2*ell,
and strong error eta. Set

 S_i,1=first ell bits of Q_i,
 R_i,t=W(X,S_i,t),
 S_i,t+1=T(Q_i,R_i,t).

For t<=a, the designated R_v,t is uniform conditional on all Q_i, all
current and earlier S values, and all earlier R values, except on total
history mass at most 4a*eta. More precisely, after the preceding messages
and the current seeds are fixed, split X=A+B by Z. Outside that error
mass the map from B to R_v,t has full output rank. It therefore remains
uniform under any further linear fixing determined by Z. This stronger
rank form is what is used in B4.

Here is a direct proof, including the conditioning order. Initially the
designated S_v,1 is uniform. At a step fix all current S_i,t. These are
linear observations of Z once earlier R values are fixed, so the B kernel
is unchanged. Its bad-seed set for W is fixed, and the marginal uniformity
of S_v,t bounds the bad history mass by 2eta. On a good seed the B-part
of R_v,t is exactly uniform and independent of A. This holds even though
the other current seeds may be arbitrarily correlated with it.

To advance, fix the A-parts of ALL current R_i,t, including the designated
one. Each is a linear observation of Z because its seed is already fixed.
The designated R_v,t is now a constant translate of a uniform B-part,
independent of A. Apply T to Q_v with that seed. Fix R_v,t and then the
other B-parts. The output S_v,t+1 is linear in A and exactly uniform except
for another bad-seed mass 2eta; the other B-fixings cannot affect it.
All next seeds are now linear functions of Z. This proves the induction.

Each completed step fixes at most h*ell bits from the seeds and h*ell
A-parts, hence reduces Q_v's rank by at most 2h*ell. The B dimension falls
by at most h*ell. Before the last needed transition, q0-2h*(a-1)ell>=C*E*z^2*ell and
the remaining B dimension is at leastC*E*z^2*ell. The initial uniform seed and the
sequential rank conclusions give the claimed cumulative error. All
conditioned components are affine: seeds are fixed before a seeded map's
output is conditioned on, and every subsequent observation is then linear.
A good-history flag is determined by a fixed transcript value and does
not impose a further nonlinear condition inside that affine component.

# B4. Many almost independent affine rows

The following explicit variant of index-based alternating extraction is
supplied in full because its conditioning order matters. It follows the
architecture of Li 2015, but refreshes the look-ahead source on EVERY
index round, including the last. No theorem from that paper is used as
an unexplained extractor premise.

Let k,h,ell,d,z be the parameters of B7, put C=2^30,E=e, and let Y_i be the initial B2
extractor output on seed i, of length Y. There are M0=2^d0 seeds.
The index alphabet is a=2 and b=d0. Write i as its d0-bit binary word,
and add one to each bit to obtain t_i,1,...,t_i,b in {1,2}. This is
injective and preserves numerical order. The independence order h is
not the index alphabet; it still bounds the number of parallel streams.

For each i start D_i,1=first d bits of Y_i. On every round j compute

 R0_i,j=E1(X,D_i,j),
 Q_i,j=E2(Y_i,R0_i,j),
 (R_i,j,1,...,R_i,j,a)=look-ahead(X,Q_i,j),
 V_i,j=R_i,j,t_i,j.

For j<b put D_i,j+1=E3(Y_i,V_i,j). The row W_i is V_i,b. It is enough
for the evaluator to compute each look-ahead stream only up to its selected
digit. The full prefix messages used in the proof are its identical
mathematical continuation and do not change the selected output.

The source/seed/output specifications are:

* E1: ambient N, entropy 10d, output ell, seed d, error eta'=2^-e'.
* E2: ambient Y, entropyC*E*z^2*q0, output q0, seed ell,
  error eta=2^-e.
* E3: ambient Y, entropyC*E*z^2*d, output d, seed ell, error eta.
* W: ambient N, entropyC*E*z^2*ell, output ell, seed ell, error eta.
* T: ambient d, entropyC*E*z^2*ell, output ell, seed ell, error eta.

Short Q_i are padded at the END before T is evaluated. Construction B2
provides all five maps with these seed lengths, as verified in B7.
Assume the following conservative bounds, also checked eventually there:

 Y>=1000*b*h*d,
 k-h*Y-100*b*h*a*ell>=10d,
 d>=(2h*a+10)ell,
 e'>=100h*a*ell+2e,
 Y-100bhd>=C*E*z^2*max(d,q0),
 q0-2h*(a-1)ell>=C*E*z^2*ell,
 k-hY-100bh*a*ell>=10d+C*E*z^2*ell.                         (B4.1)

For every affine source of entropy at least k, all but at most
2M0/N^(1/28) initial rows Y_i are exactly uniform, eventually. For ANY set of
at most h of those good indices, its W_i rows are jointly within
100*b*h^2*eta of independent uniform ell-bit rows.

## Proof by ordered prefixes

Fix such a set I, and fix a designated v in it. Let Z concatenate its
original Y_i; these are linear functions of X, and H(X|Z)>=k-h*Y.
Let E_j be predecessors i<v whose first j index digits are strictly
smaller than v's, and H_j the predecessors with the same j-digit prefix.
At j=0, E_0 is empty. At j=b, H_b is empty because the full index words
are distinct. A same-prefix predecessor's next digit is never larger
than v's digit.

We maintain an adaptive transcript of fixed seed and linear-output values.
At the start of round j, all previous chosen R values are fixed, so all
current D_i,j are linear functions of Z. The seeds D_i,j for i in E_(j-1)
are fixed as part of the transcript. Outside accumulated bad history mass,
D_v,j is exactly uniform conditional on that transcript. Moreover

 H(Y_v)>=Y-100(j-1)h*d,
 H(X|Z)>=k-h*Y-100(j-1)h*a*ell.                     (B4.2)

Initially this follows from uniformity of Y_v. The additional transcript
values described below make the induction stronger, not weaker: whenever
an A-part is fixed it is a linear function of Z, and whenever a seeded
output is fixed its seed has already been fixed.

First process the already-diverged rows E_(j-1). Their current seeds are
fixed. Fix R0_i,j, after which Q_i,j is a LINEAR function of Z, but do NOT
fix the whole Q_i,j. Instead fix its first seed S_i,1, then R_i,1, then
S_i,2, then R_i,2 if needed. In general each next S is linear in Z once
the preceding R is fixed, and each R is linear in X once its seed is fixed.
This conditions only the actual seeds and outputs up to at most a steps.
The full intermediate row Q_i,j need not be revealed.

For at most h already-diverged streams, this costs no more than
(2a+1)h*ell<=100h*a*ell bits from D_v,j and O(h*a*ell) from H(X|Z).
Thus D_v,j retains entropy at least d-100h*a*ell. Every chosen R of
these predecessors is fixed, which is precisely what makes their next
seed linear in Z. No same-prefix seed of length d has yet been fixed.

Split X=A+B by Z in the resulting affine component. The target seed is a
function of A and independent of B. Its weak-seed probability of choosing
an E1 map of deficient rank on B is at most
2^(100h*a*ell+1)*eta'<=eta, by (B4.1). Fix D_v,j first. On a good value
the B-part of R0_v,j is uniform independently of A. Now fix the current
seeds of H_(j-1); these are A-only observations and cannot spoil that
rank. Fix the A-contributions to R0_i,j for EVERY remaining stream,
including v itself. This last inclusion ensures the target R0 is purely
a translate of a B-variable, rather than a seed varying with Y_v.

The uniform R0_v,j is now independent of all Y_i. E2 therefore produces
a uniform Q_v,j except for bad-seed mass 2eta. Fix R0_v,j and then the
other B-contributions to R0_i,j. All Q_i,j are now linear functions of Z;
Q_v,j remains uniform. Its dimension is q0, and the source B entropy is
more than sufficient for B3 by (B4.1)--(B4.2).

Apply B3 to the target and the same-prefix streams up to the target's
current digit t. Outside mass at most 4a*eta, its final B-part is uniform
and stays so under further Z-linear fixings. A predecessor with digit
less than t has its chosen R already fixed in this prefix. Its next seed
E3(Y_i,chosen R) is therefore linear in Z. Fix these next seeds and those
of the already-diverged predecessors. Also fix the Q list and the A-parts
of the current R messages for all streams, including the target. These
are all A-only fixings after the current look-ahead seeds are fixed;
therefore the target final B-rank survives them. The remaining Y_v entropy
is at leastC*E*z^2*max(d,q0) by (B4.1). This meets both E2 and E3 source thresholds. The residual B entropy meets the larger W thresholdC*E*z^2*ell; the designated Q_v rank meets the T threshold by the displayed q0 guard. These are the actual charged A1-map thresholds from B2.2.

For j<b apply E3 to Y_v using the uniform target chosen R as an independent
seed. Fix this seed, and then fix the other current B-parts. Outside another
bad mass 2eta the resulting D_v,j+1 is linear in A and exactly uniform.
The other B-fixings leave this uniformity intact. Do NOT fix the next
seeds of still-same-prefix predecessors: those may be correlated with it.
They are fixed only on the next round, AFTER the target seed is fixed.
This gives precisely the induction assertion conditional on E_j next seeds.

All original Y_v rank losses in one round are bounded by 100h*d: the
seed observations and next-seed observations cost O(hd). The Q lists cost at most h*q0<=h*d, and look-ahead A-parts cost O(h^2*a*ell)<=O(hd), since d>=(2h*a+10)ell. Only O(h*a*ell) output bits
per round can reduce H(X|Z); observations determined by Z do not. These
bounds include the target's own previous seeds, Q values and A-parts.
Thus (B4.2) is preserved with the displayed conservative constants. The
failure mass per round is at most 10h*eta.

On the LAST round the same E1/E2 refresh is used. Every predecessor is
already diverged or has digit less than t, since H_b is empty. Hence all
its FINAL chosen W_i values are fixed before the target R_v,b,t. The
target row is uniform conditional on them outside the accumulated error;
no final E3 call is needed. Applying this conclusion to the ordered indices
and coupling one row at a time proves joint error at most
100b*h^2*eta. This proof permits arbitrary correlations among original
Y_i and among all nondesignated streams.

# B5. A logarithmic-seed strong generator by backward extraction

Put w=ceil(log2 r),delta=1/125,m=floor(r/125),s=1000w,c=250w. The generator has a c-bit polynomial-hash prefix, followed by repeated strong base extraction. Its full external seed length is O(w), proved below. The suffix is generally nonlinear in v; the robust argument B6 requires strong general extraction and shifted agreement, not a fixed-seed rank property.

The C1 canonical irreducible modulus of degree s defines the prefix field. For v in F_2^r pack its bits into polynomial coefficients P_v over this field. Let

 C(v;y0,z0)=first c coordinates of z0*P_v(y0).

A nonzero difference has at most r roots. Thus for distinct v,v' and any c-bit offset, prefix agreement has probability at most2^-c+r/2^s<r^-200. The full-seed collision calculation in B2.3 gives prefix strong error below r^-300 on every source of entropy at least2m under the finite guard m>=4(c+150w). The field size and c are the same as before; the remaining extractor is now constructed explicitly below.

## A general-source weak hash extractor

For the base construction use a sampled input of length ell<=r, source entropy kappa and output M1=floor(kappa/8). The field remains degree s=1000w, blocks have B=750w bits, and the weak design has y-row size a=s, universe T=5s^2+s and overlap weight5/4. There is one shared multiplier z. The exact seed length is

 D1=s+T=5s^2+2s.

There is NO output prefix in this weak map. The RecycledHash routine is invoked with code_c=0 and its unused prefix-seed coordinates fixed to zero; only its independent z and U are counted in D1. It packs the ell sampled bits into field coefficients. Every fixed seed gives a binary linear hash map. Its full design construction and polynomial-time CE proof follow.

## Exact conditional-expectation design

Construct J rows A_i, each an a-subset of [T], satisfying

 sum_(j<i)2^|A_i intersection A_j|<=(5/4)*i.            (B5.3)

The implementation compares integers4*sum<=5i. For a uniform a-subset A and a fixed previous row D0, expand2^|A intersection D0| as the number of its subsets. For a fixed j-subset of D0, its probability of lying in A is(a)_j/(T)_j<=(a/(T-a+1))^j. Consequently

 E[2^|A intersection D0|]
 <=(1+a/(T-a+1))^a<=exp(a^2/(T-a+1))<=exp(1/5)
 <=1+1/5+(1/5)^2/(2*(1-1/5))=49/40<5/4.

The series inequality uses1/n!<=1/2 for n>=2 and the geometric series; no numerical exponential approximation is assumed. Summing over i predecessors gives the required initial expectation.

For completeness, the code's exact conditional expectation is as follows. Let P be the selected positions, R the undecided positions and k the remaining cardinality, with t0=|R|. For each prior row D0 put h0=|P intersection D0| and l=|R intersection D0|. Its numerator contribution is

 2^h0 * sum_(j=0)^min(k,l) binomial(l,j)*binomial(t0-j,k-j),

with common denominator binomial(t0,k). This follows by expanding the overlap weight into subsets of R intersection D0 and counting completions containing each subset. It is the conditional average over every k-completion. Inclusion and exclusion of the next coordinate have the old expectation as their probability-weighted average, so one is no larger. Choose it by integer cross multiplication; handle k=0 or k=t0 directly. At the end the sum is at most its initial expectation, proving(B5.3), even when repeated rows occur. Each row depends only on predecessors, so the output is a deterministic prefix for fixed a,T. Binomial values and overlap weights have polynomial bit length; the number of decisions and terms is polynomial in J,T,a. Since T=O(log^2 r),J<=r, this is polynomial in r, not a hidden search or supplied advice table.

Here J=ceil(M1/B). The construction fixes a,T for all rows and is deterministic polynomial in its finite dimensions.

For ANY kappa-min-entropy source, compare hybrids with real earlier blocks and independent uniform later blocks. Fix design bits outside the current y-row and later uniform blocks, leaving the shared z unfixed. Store raw P_X(y_j) for every assignment on each earlier overlap. The table is independent of the current (y,z), determines all preceding outputs given them, and has length less than

 (5/4)*(s/B)*M1=5M1/3.

Table values below2^(-table_length-150w) have total mass at most2^-150w. The exact plan checks

 kappa-5M1/3-150w>=s.                                  (B5.4)

All other values leave general conditional entropy at least s. B2.3's collision/root calculation with ell<=r and at most750w output coordinates bounds the current strong block error by

 (1/2)*sqrt(2^(750w)*(2^(-1000w)+r/2^(1000w)))<r^-120

for r>=2. There are J<=M1<=r blocks, so the full strong error is at most2r^-119. This includes the entire weak seed, even though all rows share z. For clarity, outside current-row U bits are fixed first. Each earlier evaluation argument is then determined by its intersection assignment with the current row. The table depends only on X and those fixed outside bits, not the current y or shared z. Given that table and current(y,z), multiplication/projection reconstruct every earlier block. Conditioning on a table value therefore leaves(y,z) independent uniform, so the conditional strong hash guarantee transfers by data processing to the full hybrid with every seed bit retained. No affine-source promise is used.

## Strong two-block conversion

Let K0=floor(m/2), n=P2(r)=2^w and use A2 with sigma=1/2000,b=17,t=256,1024 coefficients per sampler. Its sample length and entropy are

 ell=floor(n/(64*2000*17)), kappa=floor(ell/16000).

The exact scalar sampler guards in A2 are required. They give strong sampler error at most(3/4)r^-150. The plan also checks

 K0>=n/1000, K0-ell-150w>=n/2000.                      (B5.5)

Take two independent sampler seeds A,B and output Y1=sample(V,A),Y2=sample(V,B). For every source V of entropy at least K0, their joint distribution, INCLUDING both seeds, is at distance at most(5/2)r^-150 from a conditional two-block source satisfying

 Hinf(Y1|A=a,B=b)>=kappa,
 Hinf(Y2|A=a,B=b,Y1=y)>=kappa

at every supported conditioning value. Here is the full repair argument. For each fixed first seed, values y of Y1 of probability below2^(-ell-150w) have total mass at most2^-150w. At every other y the conditional original V has entropy at least K0-ell-150w, enough for A2 by(B5.5). Apply A2 to the second sampler, repairing Y2 per second seed. On rare y replace Y2 by a uniform word. This retains the Y1 marginal and both independent uniform seed marginals, and makes every Y2 conditional kernel a kappa-source, at cost at most sampler error plus r^-150. Now repair the Y1 marginal using A2 for its first sampler; it is independent of the second seed. Retain the already-good Y2 conditional kernels. This costs at most another sampler error and preserves the seed marginals. It gives the asserted conditional block source, with total cost at most2*(3/4)r^-150+r^-150. Repairs are comparison distributions only; the evaluator computes only the two sampled bit strings.

## A short external seed for the base extractor

On Y2 use DIRECT A1, with dimension4 (reconstruction order80), source threshold kappa, output M2=P2(D1) and error log160w+10. There is no preceding affine condenser in this call; A1's supplied theorem is for general min-entropy sources. Its integer field parameters and all A1 hypotheses are checked by SeedParams.valid. C1 constructs its modulus and4-by4 cyclic matrix canonically; they consume zero input data bits. Here s2 is the helper field degree. Take its first D1 output bits as seed of the weak hash extractor on Y1, and return that weak map's M1 output bits.

For each fixed converter seed, the comparison pair is a block source. A1 applied to Y2 conditional on Y1 makes its D1-bit output uniform independently of Y1 and of the retained external A1 seed, up to error2^(-160w-10). The weak extractor then applies to the kappa-source Y1 with that independent uniform seed. Retain both sampler seeds and the helper A1 seed; drop the internal derived seed. Data processing and the preceding block-source comparison give base strong error at most

 (5/2)r^-150+r^-160+2r^-119<=4r^-119.

This is a strong statement in the EXTERNAL seed of length2*1024*w+5s2. No seed derived from Y2 is assumed independent before this coupling. The same fixed good helper matrix and moduli work for every conditional source and every round, by the universal supplied A1 theorem.

All helper guards hold eventually. Here ell and kappa are fixed positive fractions of r up to flooring. D1=O(w^2),M2=O(w^2), and A1's H at fixed dimension4 is O(r^(1/4)), while s2=O(w) at error log160w+10. Hence its advice cost O(r^(1/4)*w^3) is o(kappa), and every exact A1 condition passes. The sampler guards hold eventually by A2; (B5.4)--(B5.5) likewise have linear-in-r positive margins over logarithmic terms. These conclusions hold at every sufficiently large integer r, not selected powers of two, using r<=n<2r.

## Recovering the required output length

Repeat this strong base map on the SAME input V with fresh independent external seeds, concatenate M1 bits per round and truncate only the final round to m-c bits. Let

 rounds=ceil((m-c)/M1).

After earlier seeds and outputs are fixed, prior output length is at most m. Values of probability below2^(-prior_length-150w) have total mass at most r^-150. At every other history a7m/4-source retains entropy at least3m/4-150w. The plan checks4K0+600w<=3m, so this is at least K0. Fresh round seeds remain independent of this conditional source. Compare hybrids with real earlier outputs and uniform later outputs; every step costs at most4r^-119+r^-150, retaining all external seeds. Since rounds<=r, total error is at most5r^-118.

For the prefix, values of probability below2^(-c-150w) have total mass at most r^-150. All others leave entropy at least2m-c-150w>=7m/4 under m>=4(c+150w). Apply the repeated-base result, then replace the prefix by uniform using its strong error below r^-300. Total whole-generator strong error is at most

 5r^-118+r^-150+r^-300<r^-80

for all sufficiently large valid parameters. Every prefix seed, both sampler seeds and helper seed in every round, and all ignored padding bits are retained. The first c coordinates still force the unchanged shifted-agreement probability below r^-200 for every distinct input pair and fixed output shift. No rank claim about the nonlinear complete G is needed.

The number of rounds is bounded by an ABSOLUTE CONSTANT, making the external seed logarithmic rather than squared-logarithmic. If M1>0, successive floor bounds give

 ell>=n/(2*2176000), kappa>=ell/32000, M1>=kappa/16,
 M1>=n/2228224000000.

For any admitted2m<=r<=n, rounds<=1114112000000<2^41. The displayed universal constant suffices also for the actual m=floor(r/125). Because s2=O(w), the exact external seed

 2s+rounds*(2*1024*w+5s2)=O(w)                          (B5.6)

has a fixed constant independent of r and the header. There is no hidden search or unbounded seed-recycling depth.

## Actual data and uniform evaluation

The prefix modulus of degree1000w, sampler modulus of degree w and helper dimension4 field/matrix are constructed canonically by C1, consuming ZERO input bits. C1 proves exact good-data tests, existence and uniform polynomial(r) search cost; no unsupported efficient primitive-element algorithm is used. This source-specific bounded search claim is not extended to the four growing-dimension internal packages.

Scalar planning still checks all sampler/A1/error/seed guards and never allocates2^n. Binding computes the canonical data and the weak CE design. Per-round sampling, fixed-dimension helper streaming and weak hashing remain polynomial in r. The total canonical search cost is polynomial in r with large fixed exponents as proved in C1. The evaluator never computes conditional probabilities or comparison-distribution repairs.

# B6. A robust function for almost independent uniform rows

We supply a specialized bounded-independence argument below. No general
AC0 fooling theorem or non-oblivious bit-fixing extractor is assumed.
This proof uses the partition architecture of Chattopadhyay--Zuckerman,
with B5's much smaller shifted agreement and an explicit two-level
Bonferroni comparison. Only constant separation from zero and one is
needed for dispersal.

## Exact layout and bias parameter

Let there be M0=2^d0 available rows, each with at least u output bits. Choose the largest integer
r with

 m=floor(delta*r), b=floor(floor(sqrt(S_header))*ceil(log2 r)/(30*125^2)),
 u=2b+2ceil(log2 r)+10,
 2^(m+b)<=M0.                                         (B6.1)

S_header is the actual fixed input scale for this hardness argument.
The implementation finds r by integer binary search, using
m+b<=d0, an equivalent comparison. Put R=2^r, B=2^b,
M=2^m,S=BM. The robust input now consists of M0 rows of at least u bits each. Use the first u bits of EACH of the first S rows as one group. No group combines different rows.
For a fixed sufficiently large header put C_S=floor(sqrt(S_header))/(30*125^2) and require C_S>=2000. Then b=C_S*ceil(log2 r)+O(1), so B=Theta_S(r^C_S),r/B tends to zero,b=Theta(log r)=o(r), and
u=Theta(log r). Hence the maximal layout has r=Theta(log M0).
B7 chooses one header meeting this condition and the logarithmic B5 seed bound; no bias estimate is asserted for arbitrary small headers. At an adjacent integer r, m increases by at most one. The ceil-logarithm
w increases by at most one, and b increases by at most C_S*w, for a
fixed constant C_S. Maximality of the GROUP-count capacity therefore implies

 M0*2^(-C'_S*log r) < B*M <= M0.                   (B6.1a)

A constant-fraction coverage claim is NOT used. Instead the uncovered
factor is 2^poly(log r), subpolynomial in M0. We have

 b/m -> 0, r=Theta(log M0), u=Theta(log r).

For the sufficiently large fixed header of B7, B5's seed fits b.
All computation and the displayed comparisons use exact integers.

Choose an INTEGER J by bisection, as the largest 0<=J<2^u satisfying

 R*(2^(uB)-(2^u-J)^B)^M <= 2^(uBM).                      (B6.2)

Define g on a group to be one iff its positively weighted binary integer
is at least J. This is monotone in its bits. On uniform group bits let
p1=J/2^u, p2=(1-p1)^B, p3=(1-p2)^M, and alpha=Rp3.
Then

 alpha -> 1,    Mp2~ln R,    Bp1~ln M.                  (B6.3)

Here is a check that bisection hides no real constant or accuracy premise.
The continuous root with Rp3=1 has
p2=1-R^(-1/M)~ln R/M and p1~(ln M-ln ln R)/B.
The derivative of ln p3 with respect to p1 near that root is O(B ln R).
The grid step 2^-u is O(1/(B^2*r^2)), so its product with that derivative
tends to zero. The root lies between consecutive grid points chosen in
(B6.2), proving alpha->1 from below. These estimates also give 0<J<2^u
and the other two comparisons eventually. The implemented comparisons
are exact integer exponentiations on at most uBM+O(r) bits. Since
uBM<=u*M0. For each fixed header u=O_S(log M0), so this is polynomial time in M0. Uniformly over actual admitted headers, u<=ell<=N and M0<=2^(5L)<=(2N)^5, giving a fixed polynomial in N.

For each v in F_2^r the B seed outputs z_j=G(v,j) define a partition of
the S groups into M blocks

 P_v,w={(j,z_j+w): j in F_2^b},  w in F_2^m.

The group indexed (j,a) has array index j*M+a; addition is XOR. Each block
has size B and these blocks partition the groups. Blocks from DIFFERENT
partitions meet in at most eta*B groups by B5(ii), for ANY two shifts.
Define the monotone function

 H=AND_v OR_w AND_(group in P_v,w) g(group bits).         (B6.4)

Its loops have at most R*M*B iterations, with each group map computable in polynomial time in r. G includes
B5's one-time polynomial(r) weak-design construction and polynomial(r)
work per value. Since r=O(log M0), the total is polynomial time in M0. Its polynomial exponent is very
large; no modest-length efficiency assertion is made.

## Two elementary probability tools

Bonferroni follows by expanding product(1-I_i) and truncating the
alternating sum. For any finite event family, truncations of consecutive
parities bound the probability that none occurs; their gap is the
corresponding factorial moment sum S_t=sum_(|T|=t)Pr[all events in T].
If the events are independent with mean count lambda, S_t<=lambda^t/t!.
We also use t!>=(t/e)^t, obtained for example by integrating log x.

We record and prove the needed multiplicative form of Janson's inequality.
For conjunction events E_i in independent Bernoulli coordinates, let
p_i=Pr E_i<=tau<1, and let Delta=sum_(i<j, supports intersect) Pr[E_i E_j].
Then

 product_i(1-p_i) <= Pr[no E_i]
 <= product_i(1-p_i)*exp(Delta/(1-tau)).                 (B6.5)

Harris association is elementary here: induct on the last independent
coordinate; conditional covariances of increasing functions are nonnegative,
and their conditional expectations are increasing in that last coordinate.
It gives positive association of decreasing events too, proving the lower
bound. For the upper bound expose E_i in order. Write A for nonoccurrence
of earlier NONneighbors (their supports are disjoint from E_i), and B'
for nonoccurrence of earlier neighbors. A is independent of E_i. Since
E_i E_j is increasing and A decreasing, Harris gives
Pr[E_i E_j|A]<=Pr[E_i E_j]. Hence

 Pr[E_i|A B'] >= Pr[E_i B'|A]
              >= p_i-sum_(earlier neighbors j) Pr[E_i E_j].

Multiply the resulting upper bounds on conditional nonoccurrence
probabilities, and use 1+x<=exp x. This gives (B6.5). This proof also
applies to conjunctions of independent group variables with the same bias.

## Arbitrary bad coordinates

Assume at most q of the original M0 ROWS are bad. Each good row supplies u bits; every set of at most t good ROWS is jointly within gamma of independent uniform u-bit rows. Bad rows may depend arbitrarily on all good rows, including additional conditional randomness. Applying g to each row preserves this total-variation bound. It is enough that, eventually,

 q<=S^(1-delta), t>=A_t*r^2, gamma<=2^(-A_gamma*r^2*u).     (B6.6)

Here A_t,A_gamma are fixed finite positive constants, depending only on
delta and the fixed outer truncation order 100. Their existence is PROVED
by the finite degree and coefficient estimates below, not assumed as an
external extractor premise. B7 explains how the fixed total language
contains an actual input header large enough to meet all such constants.

A group is bad exactly when its supplying row is bad. There are at most q bad groups. Define H0,H1 by setting ALL bits of bad groups to zero or one,
respectively. Monotonicity gives H0<=H<=H1 pointwise. The good u-bit groups retain the asserted ROW-wise almost independence. On ideal independent good
bits, each good group has failure probability p1.

Call a partition exceptional if one of its blocks has more than
3epsilon*B bad groups, where epsilon is the strong error of B5. There
are fewer than 2^(2m)*M=2^(3m)<=R^(3delta) exceptional partitions.
Indeed, otherwise some shift w has at least 2^(2m) such v. Take their
uniform source of min-entropy 2m and apply the strong extractor G, keeping
the seed in the output. Its probability of landing in the bad-group set
exceeds 3epsilon, whereas the uniform probability is q/S<=epsilon
for large r; this contradicts strong error epsilon. The latter inequality
uses q/S<=S^-delta, which is exponentially small in r, whereas
epsilon=r^-80.

We will ignore exceptional partitions temporarily. Let F_v be the event
that partition v's OR-of-ANDs is zero. For any partition, H0 has at least
M-q wholly unaffected disjoint blocks. Under the independent model

 Pr[F_v for H0] <= (1-p2)^(M-q)=(1+o(1))/R.              (B6.7)

The last identity uses q/M->0 exponentially: from (B6.1),(B6.6),
log2(q/M)<= (1-delta)b-delta*m=-Theta(delta^2*r).
The same bound, with an additional negligible term, will be proved below
under the almost-independent good-bit distribution. Since H1's event
F_v is contained in H0's, discarding all exceptional partitions changes
either output probability by o(1), using their R^(3delta) count.

## Ideal independent probabilities on the remaining partitions

Fix any CONSTANT c<=100 distinct nonexceptional partitions and either
bad-group setting. For H0 discard blocks made impossible by a zero group;
at most q per partition are discarded. Other block events have probability
p2. For H1 each block loses at most 3epsilon*B random groups, and its
all-one probability is p2*(1+o(1)), uniformly. The product of individual
nonoccurrence probabilities over the c partitions is

 D=(p3)^c*(1+o(1)).

For H0 its logarithmic discrepancy is O(c*q*p2)=o(1). For H1 it is
O(c*epsilon*r^2)=o(1), using (B6.3). A block meets at most B blocks of
each other partition. The number of dependent pairs is at most c^2*M*B.
Their joint probabilities are at most
p2^(2-6epsilon-eta)=O(p2^2), since (epsilon+eta)r->0. Thus

 Delta=O(c^2*B*r^2/M)=o(1).

The ratio B/M tends to zero exponentially in r. Applying the supplied
inequality (B6.5) proves, uniformly in these partitions and the bad set,

 Pr[intersection of these c F_v]=(p3)^c*(1+o(1)).        (B6.8)

There are R-o(R) retained partitions. Their outer factorial moment sums
of each fixed order c<=100 tend to 1/c!, by alpha->1. The two outer
Bonferroni truncations at consecutive orders 99 and 100 consequently
place the retained function's independent-model expectation strictly
between .36 and .38 for sufficiently large r. These bounds follow also
by the elementary alternating series for exp(-1); the truncation gap
is at most 1/99!, far below .001.

## Comparison under limited independence

We give the details because pairwise independence of the good bits would
not suffice. Throughout, c is fixed, at most 100. Constants below may
depend on c and delta, but never on r, M0 or the chosen bad set.

Consider the no-block-event probability for these c partitions. For j
DISTINCT block events, j<=C*r, their union contains at least

 jB(1-3epsilon)-j^2*eta*B/2

random good groups; impossible H0 events simply have probability zero.
Using independent group bits, the joint block probability is at most
2*p2^j for large r. Therefore its factorial moment of order j is at most
2*(c*M*p2)^j/j!. Its mean parameter cMp2 is O(r). Choose an EVEN
J_inner=C*r, rounded up, with C sufficiently large that the moment at
that order is at most exp(-H*r), where H is any prescribed fixed constant.
The factorial estimate above permits an arbitrarily large H by increasing C.

To compare a joint block event under the two distributions, view it as
NO group failures on a union of at most jB groups. In the independent
model the failure count has mean at most jB*p1=O(r^2). Truncate its
Bonferroni expansion at an EVEN K_inner=C'*r^2. Choosing C' sufficiently
large makes the independent factorial remainder at most exp(-H'*r^2),
for any prescribed H'. The same remainder under the other distribution
can be compared termwise. Each term involves at most K_inner=O(r^2) good GROUPS, each supplied by one distinct row. Its expectation differs by at most gamma by the joint ROW assumption(B6.6). The u bits inside a row do not consume u different row-independence slots. The total number of terms is at most
exp(O(r^2*u)), since log(jB)=O(u) and K_inner=O(r^2). Consequently EACH joint block probability differs by at most

 2exp(-H'*r^2)+gamma*exp(O(r^2*u)).

There are at most (cM)^j=exp(O(r^2)) block tuples of order j<=J_inner.
Choose H' large enough before using the estimate. The differences of all
block factorial moments through J_inner are then exp(-H''r^2) plus
an error gamma*exp(O(r^2*u)). Apply the block-level Bonferroni bounds to
both distributions; their gap is controlled by the J_inner moment just
bounded. The probability of the c-fold F_v intersection differs by at
most exp(-H*r)+o(exp(-H*r)).

All required ROW degrees are at most C_deg*r^2, and the coefficient
counts are at most exp(C_count*r^2*u), for FIXED finite constants after
J_inner and K_inner above have been chosen. These constants are uniform
in the selected bad set and c<=100. Set A_t>=C_deg and choose
A_gamma*ln(2)>C_count+1. Then gamma*exp(C_count*r^2*u)<=exp(-r^2*u),
which is smaller than exp(-H*r^2) for every fixed H eventually. This
proves the existence and adequacy of the constants in (B6.6). They do
not depend on N, an input header, or the particular source. Choose H>1000. After weighting by
at most R^c=2^(cr) choices, every outer factorial-moment discrepancy
through order 100 tends to zero. Thus the same .36--.38 expectation
bounds hold under the actual good-bit distribution for the retained
H0 and H1.

The one-partition comparison for (B6.7) is identical with its impossible
H0 blocks discarded. All remaining blocks are disjoint and wholly good,
so their factorial moments have the same O(r) mean bound. It yields
Pr[F_v for H0]<2/R for large r, plus an error exp(-1000r). Union over
R^(3delta) exceptional partitions is o(1). Restoring those partitions
therefore puts BOTH E[H0] and E[H1] between .35 and .4 eventually.
Since H0<=H<=H1, both values of H have positive probability, regardless
of how the bad bits depend on the good ones. This proves the robust
non-oblivious bit-fixing disperser needed here.

# B7. Polylogarithmic design seed and a fixed scale header

## Fixed algorithm on every header

Let N be actual input length,L=ceil(log2 N),z=ceil(log2 L),C=2^30,D_int=max(2,z). The ordinary L-bit header encodes S-1,1<=S<=2^L<=2N. Exact scales are

 h=P2(S L^2), e=P2(S L^2 z), ell=P2(S^3 L^2 z^2),
 eprime=P2(S^5 L^4 z^2), d=P2(S^9 L^4 z^5),
 Y=P2(S^11 L^7 z^5), k=P2(S^13 L^9 z^5),
 a=2,b=d0,q0=(2h*a+C*e*z^2+10)ell.                      (B7.1)

Here h is ROW independence and e is error precision; they are deliberately different. Raw log2 k is checked before expansion. Rejected plans and N<128 output0. Initial A1 dimension11 is fixed; internal dimension is D_int. C5 decodes and C3 reconstructs the paid internal primitive packages and C1 constructs initial/B5 data canonically. Short E1 uses no data. Every header/data string is handled by the same total polynomial algorithm.

## Uniform estimates for each fixed scale

B2.2 supplies actual internal seed<25ze, paid by ell/e>=S^2*z/2 for fixed S>=8. Largest internal K remains O_S(L^6*z^8), exactly its field/advice premise. Short E1 seed O(eprime) fits d since d/eprime>=S^4*z^3/2. Entropy10d meets M+2eprime+10. Transcript precision eprime>=100h*a*ell+2e gains a power of S. Initial advice Y*k^(1/11)*L is o(k) by its strict L exponent7+9/11+1<9, so one call suffices; d0<4L and checked5L cap remain. Current prefix P=O_S(L^2z^2)=o(k).

Power-of-two rounding gives

 h*a*ell<=8S^4 L^4 z^2, eprime>=S^5 L^4 z^2,
 b*h*d<=4B0 S^10 L^7 z^5, Y>=S^11 L^7 z^5,
 hY<=4S^12 L^9 z^5, k>=S^13 L^9 z^5,                    (B7.2)

where b<=B0 L,B0=5. Fixed S>=3200 pays transcript precision,S>=8000B0 pays Y>=1000bhd,S>=16 pays hY<=k/4. Also q0<=8C*S^4 L^4 z^5 eventually because e>=h; hence q0/d<=8C/S^5 is at most1/100 for a sufficiently large fixed S. Thus d>=q0,d>=(2h*a+10)ell. Here2e=o(eprime).

All internal entropy thresholds are C*e*z^2 times output. q0-2h*(a-1)ell=(2h+C*e*z^2+10)ell>=C*e*z^2ell. Y>=1000bhd leaves at least.9Y and e*z^2*d/Y=O_S(z^3/L) tends to zero. Thus Y-100bhd>=C*e*z^2max(d,q0). The residual X entropy k-hY-100bh*a*ell exceeds10d+C*e*z^2ell: hY/k is made below1/4 by S and all other terms have smaller L powers. These are the EXACT B3--B4 thresholds and finite source guards. No h=e substitution is made. B2.2 proves all field/advice/seed guards.

## One fixed hard scale and row grouping

The new group-count layout has r<125(d0+1)<=c_r L,c_r=1000 under the checked5L cap, and B*M<=M0 rows. It needs ell>=u so that one row supplies a whole group; this holds eventually because ell is polynomial in L while u=O_S(log r). At fixed S, eventually ceil(log2 r)<=2z and u<=C0*sqrt(S)*z for an absolute constant C0.

B6 now requires h>=A_t*r^2 ROWS of almost independence, while gamma=100b*h^2*2^-e must be at most2^(-A_gamma*r^2*u). Choosing S>=2A_t*c_r^2 pays the first requirement. Choosing sqrt(S) above a fixed multiple of A_gamma*C0*c_r^2 pays e>=A_gamma*r^2u+log2(100bh^2), including its lower-order logarithm. B5's seed still fits once sqrt(S)/(30*125^2)>2C_boot; also require that ratio>=2000 for the bias proof. These are finitely many absolute scale requirements along with(B7.2). Choose ONE integer S0 meeting them. Executable guards h>=r^2,e>=100r^2u and ell>=u are resource checks, not a claim that every admitted header proves hardness. This S0 meets all required full constants eventually.

B4 gives every at most h GOOD ROWS joint distance at most gamma from independent uniform ell-bit rows. Selecting their first u bits and applying the biased group map are deterministic data processing, preserving that bound with the SAME number of rows. Thus the exact grouped-row hypothesis(B6.6) is met; no assumption on all u*h individual bits from different rows is needed. Initial bad row count is at most2M0/N^(1/28).

B6.1a gives S_groups>=M0*2^(-O_S0(log r)). Since r=O(L),
the initial bad-coordinate bound q<=2M0/N^(1/28) implies

 q/S_groups<=N^(-1/28+o(1)), S_groups^delta<=2^((4/125)L).

Here S_groups<=M0=2^d0 and d0<4L eventually, so the second bound follows with delta=1/125. Their product is N^(-13/3500+o(1)) and tends to zero. Thus q<=S_groups^(1-delta) eventually. The rebalanced precision E0 still supplies enough initial good rows: fixed-seed rank deficiency has fraction at most2*2^-E0<=2*N^(-1/28). No N^-2 initial-error claim is retained.
Every B6 hypothesis holds. This ONE S0 fits the header once2^L>=S0,
and every ready guard holds eventually, for ALL sufficiently large
actual integer lengths.

This is a lower-bound restriction, not a nonconstructive definition of
the language. The evaluator already computes all header values. It never
receives S0 externally. The same algorithm and language work at every
length and have no good-header promise.

## Actual prefix and uniform computation

After the L-bit header, the input cursor consumes only the four internal packages in order T,W,E2,E3. Initial A1 and short E1 consume zero bits; B5 consumes zero bits. The four packages each have t-a(t) bits,t=D_int*s; C5 decodes and C3 reconstructs their subfield and relative companion action, or returns a zero map at a failed rank guard. C1 constructs initial/B5 data canonically with proved good hypotheses and uniform bounded searches. All actual bits are charged and the cursor is checked. P=L+sum_internal(D_int*s-a(D_int*s))=O_S(L^2z^2)=o(k). The remaining N-P actual payload bits are padded with P zeros to length N. Good internal data exist and are fixed only as a paid hard restriction.

On arbitrary headers, scalar bit lengths are polynomial in N and large k rejects before expansion. Every admitted dimension is polynomial in N by explicit seed/source/prefix guards. Initial/helper A1 dimensions are fixed; internal dimension D_int=O(loglog N) streams only the requested monomials with polynomial delay as proved in B2.2. The number of calls is at most the source threshold; each call, all PV/ring/matrix arithmetic, CE designs, samplers and the robust loops are polynomial in their finite dimensions. The checked d0<=5L cap bounds seed enumeration by a fixed polynomial in N. C1 canonical searches cost a fixed polynomial in N at every admitted header, using the initial seed cap and uniform B5 parameter bounds. Rejected plans and bad internal admitted data define Boolean values in the same polynomial time. Hence the complete language is in P and has the stated empty-witness NP verifier.

The initial d0 estimate, bad-row precision and robust layout are retained from earlier attributed candidates. The fixed header chosen above gives r=O(L),u=O_S(z), correct full-seed extraction and shifted agreement for B5, and adequate independence/error in B6. Its affine-source robustness and gate application are unchanged. The exact actual-length loss is4P+8k=O_S(L^9 z^5). Both canonical claims assert that whole theorem and its specified eventual weakenings for this SAME family.

## Provenance of this revision

The scale-header architecture, A1/B1/B3--B7 base proofs and original code derive from9597bd82. The self-contained gate theorem is retained from28feddfa; charged-prefix restriction accounting traces to db5f9c6e (its legacy code is removed). Later published candidates510a59cf,19095349,b71edbca,b5cfccde,7fa479d9,b4e82f7b,681c4fbf and35570f1a supply the direct short hash, merged evaluator, initial resource accounting, half-output/overlap history, D11 balance, exact A1 reconstruction condition, adaptive B1 allocation and shared-multiplier weak hash. Unpublished9b6e4ab5 supplies half-output adaptation ancestry; its failed review is no mathematical premise. The logarithmic external-seed bootstrap, A2 and11/4 layout closure were supplied in25690f4d. The earlier variable-dimension revision reused them, changed the four internal maps to charged dimensionD_int B2.1 calls and closed the9,7 scales. The current whole-row revision separates h/e and closes9,5. Every required argument is supplied; predecessor verdicts and scores are not premises. Primary mathematical sources are credited in seeded.md and the retained supporting arguments.

The variable-dimension budget and9,7 closure are new relative to90a83107, which supplied the charged direct-internal-map method at fixed dimension1024 with901/100,4 loss.

C1 canonical initial/B5 searches and the exact prefix saving are new relative to61c0ea3e; the four internal packages remain charged.

C2 companion-package compression derives from the full-matrix precursor f72b2eba and removes sum_internal(D_int^2-D_int)s actual bits without changing seed, call, source or layout parameters. Full primitive-companion existence and exact loss accounting are supplied in C2.

C3 primitive binary packages and bounded subfield/relative-basis reconstruction now replace the C2 internal encoding. Compiler rank guards return zero maps on failure, all data bits remain charged and fixed-seed linearity survives. All extraction guarantees apply under one good primitive tuple per map, whose existence and representation-invariant A1 hypotheses are supplied in C3.

C4 further encodes a primitive reciprocal representative in D_int*s-2 bits, using its supplied direct polynomial-time decoder and a canonical t=2 exception. Failure returns zero. Every primitive reciprocal pair is represented, so the same universal good-data restriction exists. That preceding reciprocal-only representation saves sum_internal(s+2) versus C2. Current C5 uses its own explicit a(t) and prefix comparisons below.

The whole-row grouped robust input and separated h/e closure are new in this revision. Earlier encodings C1--C4 and their local exact savings are retained; only the robust row consumer/layout and resource scales change.

C5 now compresses the paid polynomial representation by excluding all small irreducible factors through a deterministic CRT unit code. The exact package length, complete good-primitive coverage, positive Euler-product counting bound and uniform polynomial decoder are supplied in seeded.md. The whole-row9,5 resources are unchanged.
