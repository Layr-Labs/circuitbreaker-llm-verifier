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

This is NOT claimed O(K+log n+e) for arbitrary parameters. In each actual application here e>=c*log n eventually for a fixed absolute c>0: initial precision has c=1/28 up to constant slack, and internal e=h grows faster than log n. Hence ceil(log2 n)/a is bounded by a fixed constant, making nprime=O(K+log n+e) in these uses. On all other admitted calls the displayed general finite-parameter bound remains polynomial. No ambient-length term is omitted from source or time accounting.

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

For the initial map, consider an affine extractor plan for source length n, entropy K, output
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

Initial extraction keeps D=11 FIXED. It has source k, output Y and precision E0=ceil(L/28)+10, with B7's k,Y and L=ceil(log2 N). At every fixed scale one call eventually suffices. Its advice order Y*k^(1/11)*L has L exponent9+1+1=11 but z exponent7+9/11=86/11<9. Therefore it is o(k); strict L-power slack is not required. The generic A1 argument applies unchanged at this fixed dimension.

For D=11, r=192 and r-D=181, E0=ceil(L/28)+10. Exact A1 rounding gives s0<=(192/1267+o(1))L. The adaptive B1 block exponent is a0=(1/28+o(1))L and its field-rounding target is(15/14+o(1))L. Therefore its seed is at most(15/7+o(1))L and its output is O(k). Combining12s0 with this seed yields

 d0<=(2304/1267+15/7+o(1))L<4L.

The explicit5L cap passes eventually. Initial122-element data, header and B5 modulus give P<=(24691/1267+o(1))L<20L. Each estimate is at fixed scale and no numerical threshold is claimed.

## B2.2 Direct half-output extraction with hash blocks

The four internal maps use HalfSeedPlan with integer overlap parameter t>=1,rho=2^t. For positive n,K,M,E require K<=n,M<=K/2 and K>=32B, where B=E+ceil(log2 M)+10. Apply B1 at precision E+10, giving exact length nprime. Put s=P2(8*(B+ceil(log2 nprime))),a=s,T=ceil(2a^2/t)+a. Impose the additional exact guard

 K*B>=4rho*s*M.                                       (B2.2a)

Use the proved deterministic degree-s tower field. The seed consists of the B1 seed, ONE independent shared field multiplier z, and T uniform design bits U. Input data length is zero; exact seed length is B1.s+s+T. Construct J=ceil(M/B) rows A_i of size a in [T], with

 sum_(j<i)2^|A_i intersection A_j|<=rho*i.              (B2.2b)

Existence and the exact algorithm follow the full CE calculation in B5. For a uniform a-subset A and fixed a-subset B0, expected2^|A intersection B0| is at most exp(a^2/(T-a+1))<=exp(t/2)<2^t. The exact series bound exp(1/2)<=5/3<2 proves the strict inequality at integer t. Summing over predecessors and choosing no-larger conditional expectations at every bit gives(B2.2b); the code checks the final integer sum. B5 has its separately proved5/4 design below.

Pack the condensed source into polynomial coefficients. Row i gets ONLY y_i=U|A_i in F_(2^s), and emits the first c_i=min(B,M-iB) coordinates of z*P_X(y_i). Thus the same z is used on every row. For any H-min-entropy source the full-seed collision/root proof in B2.3 gives strong block error

 (1/2)*sqrt(2^c_i*(2^-H+nprime/2^s))<2^-B

when H>=7B. The chosen s and c_i<=B imply the required field/error margins. This calculation applies to general conditional sources, not just affine sources.

Here is the changed conditioning proof, including the shared multiplier. Fix U outside the CURRENT y-row and later independent uniform hybrid blocks, but DO NOT FIX z. For each prior row, list the RAW field value P_X(y_j) for every assignment on its overlap with the current y-row. Each raw value has s bits, so this table has length

 s*sum_(j<i)2^overlap<=s*rho*i<s*rho*M/B<=K/4.

The strict middle inequality uses i<ceil(M/B), hence Bi<M, including the partial final block. The table is a function of X and the fixed outside U only: it is independent of current-row seed bits and of z. Given it, outside U, current-row bits and z, every prior output is determined by multiplication and projection. It is therefore sufficient for the data-processing step in the hybrid, even though prior outputs share z with the current block. Fixing a table value leaves (y_i,z) uniform and independent of the conditional source. No prior output is fixed in a way that prematurely fixes z.

Table values of probability below2^(-table_length-B) have total mass at most2^-B. Every other value leaves entropy at least3K/4-B>=23B>=7B. Apply the general block collision bound to this conditional source and the independent pair(y_i,z), retaining both seed components. Each step loses at most2*2^-B. All J<=M steps lose at most2M*2^-B<=2^(-E-9). Add the B1 bad-seed mass at most2^(-E-9). Total strong error is below2^-E, INCLUDING the B1 seed, shared z and every U bit. Larger affine inputs follow from parallel K-flat decomposition or the image-rank bound. Every fixed complete seed still defines a binary linear map.

All construction/evaluation operations have polynomial finite-parameter bit complexity by B5's exact binomial CE formula and tower-field Horner evaluation. The final block alone is truncated. For actual maps t=floor(log2 L),rho<=L and z/2<=t<=z with z=ceil(log2 L). Here E=h dominates logarithms of every admitted asymptotic source/output dimension. Hence s/B<=32 eventually and K=128rho*M meets(B2.2a). Exact guards are checked on every finite plan; no bound on s/B is silently applied to arbitrary headers. Seed cost is O(E^2/z+E), as before, with smaller row universe because a=s rather than2s. B3--B4 and B7 explicitly pay the larger source thresholds128rho times output. Short E1 still uses B2.3; only initial A1 uses input matrix data.

For an affine source and any fixed seed, output is uniform on an affine image. If it is not uniform on the full output space its distance is at least1/2. Thus strong error eta bounds bad-seed fraction by2eta; an independent weak seed of min-entropy a-t can hit that bad set with probability at most2^(t+1)*eta. These observations used in B3--B4 are unchanged.

## B2.3 Multi-bit polynomial hashing for short extraction

Let n,K,M,E be positive integers with K<=n, M<=K/2 and K>=M+2E+10. Put s=P2(8*(M+E+ceil(log2 n)+10)). Use the deterministic degree-s tower field of B1. Pack the n input bits injectively into coefficients of a polynomial P_X over that field, padding the last coefficient with zeros. For a uniform seed(y,z) in F_(2^s)^2 output the first M tower-basis coordinates of z*P_X(y). Since M<=s, this coordinate projection is surjective. The exact seed length is2s, input data length is zero, and every fixed-seed map is binary linear. Unused stipulated seed bits are ignored.

Here is the full strong-extraction proof for ANY source with min-entropy at least K. Independent source copies coincide with probability at most2^-K. For distinct source points the polynomial difference is nonzero and has degree below ceil(n/s), hence at most n roots. At every other y, multiplication by its value is a bijection on the field, so uniform z makes the M-bit difference exactly uniform. Averaging over the ENTIRE seed(y,z), the conditional output collision probability is at most2^-M+2^-K+n/2^s. Cauchy--Schwarz therefore bounds statistical distance between(seed,output) and(seed,uniform M bits) by

 (1/2)*sqrt(2^M*(2^-K+n/2^s)).

The entropy margin K>=M+2E+10 bounds the first summand after multiplication by2^M by2^(-2E-10). The chosen field degree is at least M+2E+ceil(log2 n)+10, so the second summand has the same bound. The total strong error is at most(1/2)*sqrt(2^(-2E-9))<2^-E. This argument includes all seed bits, does not condition on future outputs, and requires no rank-condenser or external extractor theorem. Affine sources are included, so B3--B4's fixed-seed rank and weak-seed estimates apply.

The evaluator uses tower multiplication and Horner's rule on ceil(n/s) coefficients and then projects the first M bits. It has polynomial bit complexity in n,s,M on every admitted call. For short E1, M=ell, E=Eprime, and ambient n=N; Eprime dominates ell and log N at every fixed scale. Thus exact seed2s=O(Eprime), fitting the stipulated d once the fixed scale is large enough. K=10d meets the entropy margin by the explicit scale comparisons in B7. No field, matrix, seed table or good construction is read as input data for this map.

# B3. Parallel affine look-ahead extraction

Let at most h streams Q_i be affine functions of Z, each padded to a common
ambient bound. Let a<=h be the look-ahead time limit, separate from the
number of streams. One Q_v is exactly uniform on q0=(2h*a+10)ell bits.
Let rho be the overlap weight in B2.2 and require q0-2h*(a-1)ell>=128rho*ell. Let H(X|Z) be at least h*a*ell+128rho*ell. The extractors W on X and T
on Q_i have output and seed length ell, source entropy threshold128rho*ell,
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
by at most h*ell. Before the last needed transition, q0-2h*(a-1)ell>=128rho*ell and
the remaining B dimension is at least128rho*ell. The initial uniform seed and the
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

Let k,h,ell,d be the parameters of B7 and let Y_i be the initial B2
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
* E2: ambient Y, entropy128rho*q0, output q0, seed ell,
  error eta=2^-e.
* E3: ambient Y, entropy128rho*d, output d, seed ell, error eta.
* W: ambient N, entropy128rho*ell, output ell, seed ell, error eta.
* T: ambient d, entropy128rho*ell, output ell, seed ell, error eta.

Short Q_i are padded at the END before T is evaluated. Construction B2
provides all five maps with these seed lengths, as verified in B7.
Assume the following conservative bounds, also checked eventually there:

 Y>=1000*b*h*d,
 k-h*Y-100*b*h*a*ell>=10d,
 d>=(2h*a+10)ell,
 e'>=100h*a*ell+2e,
 Y-100bhd>=128rho*max(d,q0),
 q0-2h*(a-1)ell>=128rho*ell,
 k-hY-100bh*a*ell>=10d+128rho*ell.                         (B4.1)

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
is at least128rho*max(d,q0) by (B4.1). This meets both larger E2 and E3 source thresholds. The residual B entropy meets the larger W threshold128rho*ell; the designated Q_v rank meets the T threshold by the displayed q0 guard. These replace the former constant output-ratio thresholds.

For j<b apply E3 to Y_v using the uniform target chosen R as an independent
seed. Fix this seed, and then fix the other current B-parts. Outside another
bad mass 2eta the resulting D_v,j+1 is linear in A and exactly uniform.
The other B-fixings leave this uniformity intact. Do NOT fix the next
seeds of still-same-prefix predecessors: those may be correlated with it.
They are fixed only on the next round, AFTER the target seed is fixed.
This gives precisely the induction assertion conditional on E_j next seeds.

All original Y_v rank losses in one round are bounded by 100h*d: the
seed observations and next-seed observations cost O(hd), while Q lists
and look-ahead A-parts cost O(h^2*a*ell)<=O(hd). Only O(h*a*ell) output bits
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

# B5. A shared-multiplier weak-design generator

Write w=ceil(log2 r),delta=1/125,m=floor(delta*r),s=1000w,c=250w and B=750w. A good actual-input modulus defines F_(2^s). Pack v in F_2^r injectively into coefficients of P_v over this field, degree below ceil(r/s). The complete unpadded seed consists of independent(y0,z0) in F_(2^s)^2, one further independent shared multiplier z in F_(2^s), and T uniform bits U. Put a=s,T=5a^2+a,J=ceil((m-c)/B). The output G is the c-bit prefix

 C(v;y0,z0)=first c coordinates of z0*P_v(y0),

followed by J blocks, truncated only at the end, whose i-th block is the first min(B,m-c-iB) coordinates of z*P_v(U|A_i). The y-rows A_i below have size s, not2s. Exact seed length is3s+T; larger stipulated seeds are padded with ignored bits. Every fixed complete seed gives a binary linear map of v. Only the modulus has input data; the design is computed.

For a nonzero v, P_v has at most r roots. At every other field point y, multiplication by P_v(y) makes uniform z exactly uniform. Thus every c-bit offset has probability at most

 2^-c+r/2^s.                                           (B5.1)

For any H-min-entropy source, the supplied collision/root calculation gives full-seed strong block error at most

 (1/2)*sqrt(2^b*(2^-H+r/2^s))                           (B5.2)

for b projected bits. This proof is exactly B2.3 with n=r and arbitrary conditional sources. For H>=1000w and b<=750w it is below r^-120 eventually: the two terms under the square root are at most2^-250w and r*2^-250w, respectively. For H=2m the c-bit prefix has error below r^-300 eventually, and(B5.1) is below r^-200. Here m grows linearly in r while w is logarithmic; c=o(m). No claim is made that toy or small admitted plans meet these eventual margins.

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

## Strong extraction with the shared multiplier

First take any source X of min-entropy at least7m/4 and independent NW seed(z,U). In hybrid step i fix U outside A_i and all later independent uniform hybrid outputs, but leave z uniform. For each preceding row list RAW P_X(y_j) for all assignments on its overlap with A_i. This table has at most

 s*sum2^overlap<=(5/4)*s*i<(5s/(4B))*(m-c)
                 <=5m/3

bits, because Bi<m-c for i<J, s/B=4/3, and the partial final block is included. It depends on X and outside U, not on z or the current y seed. Given the table, outside/current U and z, the previous outputs are determined. Thus conditional data processing reduces the step to the strong collision bound for the current independent pair(y,z), even though z is shared across rows.

Table values of probability below2^(-table_length-100w) have total mass at most2^-100w. All others leave entropy at least

 7m/4-5m/3-100w=m/12-100w>=1000w

eventually. Apply(B5.2) to the general conditional source. The step error is at most r^-120+2^-100w. There are J<=r steps, so the total NW error is at most2r^-99 eventually, retaining the ENTIRE shared seed and outside bits.

Now let X have entropy at least2m and include the independent prefix seed(y0,z0). For each fixed prefix seed, prefix values of probability below2^(-c-100w) have total mass at most2^-100w. Every other value leaves entropy at least2m-c-100w>=7m/4 eventually. The NW seed(z,U) remains independent of this conditional source. Apply the previous paragraph, then replace the prefix with an independent uniform c-bit word using its strong error below r^-300. Total full-seed error is at most

 2r^-99+r^-100+r^-300<r^-80

eventually. This includes every prefix seed, the shared multiplier and every design coordinate. It makes no future-output or unproved independence assumption. For affine sources, rank deficiency of a fixed-seed linear output has distance at least1/2 from uniform, so the bad-seed fraction is at most2r^-80.

For any distinct v,v' and any m-bit shift w0, the equation G(v,j)=G(v',j)+w0 forces its first c coordinates to agree with the fixed offset. By(B5.1) its probability over the complete seed is at most2^-c+r/2^s<r^-200 eventually. Extra ignored padding bits preserve this. This supplies the shifted-agreement bound used in B6 as well as the strong extraction/rank statement.

Construction uses polynomial(r) exact CE arithmetic. Evaluation uses field Horner rules at the prefix point and the J row points, with the same shared multiplier. All field degrees are O(log r), input data has exactly s modulus bits, and every operation is polynomial in r on good or bad admitted data. There is no field search and no design table in the input prefix. Exact seed3s+T<=C_design*w^2 for a fixed constant; this is the sole B5 seed requirement paid by B7.

# B6. A robust function for almost independent bits

We supply a specialized bounded-independence argument below. No general
AC0 fooling theorem or non-oblivious bit-fixing extractor is assumed.
This proof uses the partition architecture of Chattopadhyay--Zuckerman,
with B5's much smaller shifted agreement and an explicit two-level
Bonferroni comparison. Only constant separation from zero and one is
needed for dispersal.

## Exact layout and bias parameter

Let the available bit string have length M0=2^d0. Choose the largest integer
r with

 m=floor(delta*r), b=floor(floor(sqrt(S_header))*ceil(log2 r)^2/(30*125^2)),
 u=2b+2ceil(log2 r)+10,
 u*2^(m+b)<=M0.                                         (B6.1)

S_header is the actual fixed input scale for this hardness argument.
The implementation finds r by integer binary search, using
ceil(log2 u)+m+b<=d0, an equivalent comparison. Put R=2^r, B=2^b,
M=2^m, S=BM. Use the first u*S input bits as S consecutive u-bit groups.
For EVERY fixed header S_header, b=Theta(log(r)^2)=o(r), and
u=Theta(log(r)^2). Hence the maximal layout has r=Theta(log M0).
At an adjacent integer r, m increases by at most one. The ceil-logarithm
w increases by at most one, and b increases by at most C_S*w, for a
fixed constant C_S. The adjacent ratio of u is at most two eventually.
Maximality therefore implies

 M0*2^(-C'_S*log r) < u*B*M <= M0.                   (B6.1a)

A constant-fraction coverage claim is NOT used. Instead the uncovered
factor is 2^poly(log r), subpolynomial in M0. We have

 b/m -> 0, r=Theta(log M0), u=Theta(log(r)^2).

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
uBM<=M0, they require polynomial time in M0.

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

Assume at most q of the original M0 bits are bad, and the remaining bits
have the following almost independence: every set of at most t of them
is within gamma of uniform. The bad bits may have arbitrary joint dependence on the good bits,
including additional conditional randomness. It is enough that, eventually,

 q<=S^(1-delta), t>=A_t*r^2*u, gamma<=2^(-A_gamma*r^2*u).     (B6.6)

Here A_t,A_gamma are fixed finite positive constants, depending only on
delta and the fixed outer truncation order 100. Their existence is PROVED
by the finite degree and coefficient estimates below, not assumed as an
external extractor premise. B7 explains how the fixed total language
contains an actual input header large enough to meet all such constants.

Call a group bad if it contains any bad bit. There are at most q bad
groups. Define H0,H1 by setting ALL bits of bad groups to zero or one,
respectively. Monotonicity gives H0<=H<=H1 pointwise. The good groups
retain their asserted almost independence. On ideal independent good
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
can be compared termwise. Each term involves at most K_inner good groups,
hence at most K_inner*u=O(r^2*u) good underlying bits. Its expectation
differs by at most gamma by (B6.6). The total number of terms is at most
exp(O(r^2*u)), since log(jB)=O(u) and K_inner=O(r^2). Consequently EACH joint block probability differs by at most

 2exp(-H'*r^2)+gamma*exp(O(r^2*u)).

There are at most (cM)^j=exp(O(r^2)) block tuples of order j<=J_inner.
Choose H' large enough before using the estimate. The differences of all
block factorial moments through J_inner are then exp(-H''r^2) plus
an error gamma*exp(O(r^2*u)). Apply the block-level Bonferroni bounds to
both distributions; their gap is controlled by the J_inner moment just
bounded. The probability of the c-fold F_v intersection differs by at
most exp(-H*r)+o(exp(-H*r)).

All required bit degrees are at most C_deg*r^2*u, and the coefficient
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

Let N be the ACTUAL input length, L=ceil(log2 N), z=ceil(log2 L).
The first L input bits encode S-1, so 1<=S<=2^L<=2N. The algorithm
evaluates every S without knowing or searching for a good scale. Let
P2(t) be the least power of two at least positive integer t. The exact power-of-two parameters are

 k=P2(S^13*L^11*z^9), h=e=P2(S*L^2*z^2),
 ell=P2(S^3*L^4*z^3), e'=P2(S^5*L^6*z^5),
 d=P2(S^9*L^6*z^5), Y=P2(S^11*L^9*z^7),
 a=2, b=d0.                                              (B7.1)

All operations are exact integer arithmetic. Initial A1 dimension is 11 FIXED; short E1 uses B2.3 and the other four internal maps use B2.2, and B5 uses its supplied weak-design recycling
map rather than a design A1 package. The raw log2 k is checked before expansion. Rejected large-length
plans return zero; small lengths N<128 output0.
Bad field data on an admitted plan still give bounded ring/matrix evaluation.

Eligibility is proved for the ONE fixed algorithm on every header. The
lower-bound argument separately proves that one fixed integer S0 makes a
hard actual-input restriction at every sufficiently large length. Header
and data bits are charged; S0 is not an external input, advice or witness.

## Uniform estimates for each fixed scale

The four half-output maps use B2.2 with seed cost at most C_seed*(e^2/z+e) plus lower terms. Short E1 uses B2.3 with total seed at most C_short*eprime plus lower terms. Their deterministic tower field degrees are O(e) and O(eprime), respectively, at every fixed scale. The implied constants are independent of fixed S once N is large enough for it, since log S/z->0. Both entropy guards hold eventually. All five internal maps have zero data length.

Initial output Y is much smaller than k. At D=11 the advice ratio is O_S(z^(-13/11)), because the L powers cancel and7+9/11-9=-13/11. Thus one call eventually suffices. The preceding exact formulas give initial field degree<=(192/1267+o(1))L, combined seed d0<4L and full prefix P<20L eventually. The checked5L seed cap therefore passes. Since the index alphabet is BINARY, b=d0<=B0*L with B0=5.

Power-of-two rounding costs at most two. The key comparisons are

 e^2/t <=8S^2L^4z^3,     ell>=S^3L^4z^3,
 e' <=2S^5L^6z^5,         d>=S^9L^6z^5,
 h*a*ell <=8S^4L^6z^5,     e'>=S^5L^6z^5,
 bhd <=4B0S^10L^9z^7,     Y>=S^11L^9z^7,
 hY <=4S^12L^11z^9,        k>=S^13L^11z^9.                (B7.2)

The new overlap entropy budgets also hold eventually. Since rho<=L and h>=S*L^2*z^2, q0-2h*(a-1)ell=(2h+10)ell>=128rho*ell. The source bound Y>=1000bhd leaves Y-100bhd>=.9Y>=128rho*max(d,q0): q0<=d, b>=1 and h/rho tends to infinity. The remaining X budget is at least10d+128rho*ell, since hY/k can be made below1/4 by fixed S and all terms bh*a*ell, d and rho*ell have smaller L powers than k. Thus all four HalfSeedPlan calls satisfy K=128rho*M, the B3/B4 source thresholds and K>=32B eventually. No rho*M loss is charged as though it were O(M).

Each required constant ratio gains a power of S. A sufficiently large
fixed S meets all seed, advice, weak-seed and row-rank requirements. For
example S>=3200 leaves weak-seed room, S>=8000B0 leaves Y>=1000bhd
room, and S>=16 gives hY<=k/4. The other source loss bh*a*ell has
smaller L power than k, as does10d. Also q0=(4h+10)ell<=56S^4 L^6z^5, so q0/d<=56/S^5<=1/100 once S is a sufficiently large fixed constant. This is a constant-fraction bound, not an o(d) assertion. It implies d>=q0. Also2e=o(e').
The code retains the actual checks with h*a, not a mistaken h^2 alphabet.

## One fixed hard scale

The layout uses delta=1/125 and r<(d0+1)*125<=c_r L, c_r=1000, under the checked5L cap.
For every fixed S, eventually ceil(log2 r)<=2z and
u<=C0*sqrt(S)*z^2, for an absolute finite C0. B6 requires

 t>=A_t*r^2*u, gamma<=2^(-A_gamma*r^2*u).

Thus h=e>=S*L^2*z^2 meets both requirements, including the polynomial
prefactor in100bh^2*2^-e, once sqrt(S) exceeds the corresponding fixed
constants times C0*c_r^2. B5's seed also fits once sqrt(S) exceeds its
fixed C_design bound. Choose ONE integer S0 exceeding those finitely
many bounds and the seed/advice/coefficient bounds in (B7.2).
Executable checks h>=r^2*u,e>=100r^2*u are resource checks, not assertions
that every admitted header certifies dispersion; this S0 meets the full
proved constants eventually.

B6.1a gives S_groups>=M0*2^(-O_S0(log r))/u. Since r=O(L),
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

The first L bits are the header, followed by only the initial A1 package, then B5's hash modulus. All five internal maps use deterministic tower fields and have zero data bits. An A1 package has exactly (D^2+1)s bits at its stated dimension. After reading the header the engine computes
all lengths, starts its field cursor AFTER that header, and removes the
ENTIRE header-plus-data prefix before padding the payload to N coordinates.

Only the initial map has A1 data, with fixed dimension11 and field degree O(L). B5 adds O(log r)=O(loglog N) modulus bits and the header has L bits. Thus the complete prefix is P(N,S0)=O(L)=o(k(N,S0)), including every matrix and modulus.

The large condenser modulus is still computed deterministically. One good
tuple works simultaneously for all affine sources at the fixed scale and
length. Zero padding injects every k-flat on N-P actual payload bits into
the N-coordinate source. The payload restriction is an affine disperser.
All its actual input bits are charged and P+k<=N eventually.

Uniform eligibility also covers all other headers. S<=2N, and k>N plans
reject early. Every remaining scalar is polynomially bounded in N. The
component/source checks and total-prefix guard bound admitted fields and
matrices. The d0 cap bounds all seed enumeration; r=O(log N) bounds
robust-function loops by fixed polynomials. Only initial A1 uses fixed dimension11 and streams its first nprime coefficient monomials with polynomial delay. Internal B2.2 weak-design rows and the B2.3 multi-bit polynomial hash have polynomial finite-parameter cost. B5's weak design is explicitly computed with exact
conditional expectations in polynomial(r) time, not supplied as advice. Field degrees and all bit operations are polynomially
bounded. Rejected large plans give zero. Small finite cutoff lengths use
zero. Thus eligibility does not depend on
finding S0, or on any per-length advice or unproved efficient construction.

## Provenance

A1 supplies the Shaltiel--Umans reconstruction, B1 the GUV graph/binary
linearity and deterministic tower. B4 follows Li's index-based affine
architecture with a supplied binary alphabet, reduced prephase transcript,
correct A/B conditioning and every-round refresh. B6 follows the
Chattopadhyay--Zuckerman partition architecture with the supplied shifted-
agreement generator and moment comparison. No modern extractor theorem
or prior acceptance is a premise.

Primary provenance: [GUV](https://people.eecs.berkeley.edu/~venkatg/pubs/papers/PV-condenser.pdf),
[Li2015](https://eccc.weizmann.ac.il/report/2015/121/download), and
[Chattopadhyay--Zuckerman](https://annals.math.princeton.edu/wp-content/uploads/annals-v189-n3-p01-s.pdf).
The finite full-B2 gate proof is separately supplied in elimination.md.

The joint resource choice is D=11,E0=ceil(L/28)+10,delta=1/125 and checked d0<=5L. The derived eventual bounds d0<4L and P<20L include all header/data bits. Initial bad-row fraction<=2N^(-1/28), rather than N^-2, suffices because delta*d0<=(4/125)L leaves the strict (13/3500)L exponent margin above. The layout denominator and bias formula both use125; no old10^12 denominator remains in the defined family. This reduces the seed enumeration, prefix and robust-loop constants while keeping the11/9 remainder and leading(4,1) projection.
