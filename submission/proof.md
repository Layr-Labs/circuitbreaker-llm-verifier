# One fixed uniform family with lower bound 4N-O((log N)^25/loglog N)

The supplied evaluator defines one fixed total language in P, hence NP,
with the following unconditional theorem in the full pinned B2 model:
there exist fixed K>0 and N0 such that for EVERY integer N>=N0,

    C_B2(F_N)>=4N-K*(log_2 N)^25/(log_2 log_2 N).

Both structured claims assert this same complete theorem, its exact
version-2 Big-O declaration and its joint projection (4,1). The seven
evidence files supply every material construction and gate argument.
There is one family, not a separately chosen language for each length,
coefficient, exponent or source subspace. Every N is actual input length.

## Provenance and mathematical organization

The construction, primitive evaluator and full component arguments derive
from fbadbf8d-fdb1-495b-a11f-2aba6ad4a9d2, source commit
41f5c012ad9745b26125b387924c77861d873995. Their authorship is not claimed
here. Prior verdicts, score records and attribution decisions are not
premises. The inherited kernel traces to 11e8b36d-cbd4-44f6-810f-7a7f875b4ad0
and its fair-state framework to 237016f9-ece3-4a05-9933-dfcf6056e787.
The supplied manuscripts preserve the underlying literature attribution.
All arguments needed here are supplied anew with the changed parameters.

The new parameter allocation increases the short-extractor output,
uses a degree-64 orbit extension and degree floor(k^(1/32)) polynomial
encoding, and pays formerly strict logarithmic exponent margins with
explicit fixed constants. It also retains the saving in the number
of advice stages, which is order log N/loglog N. The proof below audits
the resulting borderline transcript budget rather than treating equal
powers as asymptotically separated.

The short extractor is a folded affine condenser followed by the basic
Shaltiel--Umans orbit construction (seeded.md). Fixed-seed maps are
linear, giving exactly uniform good rows on every affine source.
Indexed alternating extraction supplies approximately independent bounded
marginals (correlation.md), and a monotone polynomial-graph resilient
disperser handles arbitrary bad rows (resilient.md). The full-B2 kernel
4N-8(d+1) is reproduced in elimination.md. These are proof obligations,
not consequences of the finite primitive audit.

## 1. Exact total family

At N<2 return zero. Otherwise let

    L=floor(log_2 N), J=ceil(log_2(L+2)),
    k=floor(2^10000*L^25/J).

These are integer computations: L is N.bit_length()-1 and J uses integer
bit lengths. If k>N, return zero. In all remaining shapes put

    m=floor(k/(2^300*L^2)),
    h=the least power of two >=2^200*L^2,
    ell=2^500*L^4, C_L=2^4110,
    Qlen=(16h^2+16C_L)*ell,
    d=2^64*(h^3*ell+L)^2,
    E=2^200*L^2, E_1=2^20*h^3*ell+E.

For the short extractor, q1 is the least power of two at least N^16*k;
q2 is the least power of two at least N^128. Write w1=log_2 q1 and
w2=log_2 q2. Its orbit extension has dimension 64 over F_(q2), polynomial
encoding degree below floor(k^(1/32)), output length m, and total seed
length s=w1+65w2. Put M=2^s, l=log_2 h and b=ceil(s/l).

FamilyShape.build lists exact coefficient-capacity, reconstruction-advice,
rank-budget, code/design entropy and seed guards, the guard M<=N^16384,
and the resilient-order/error guards. Any failure returns zero. Thus
these conditions are part of a total algorithm, without any claim that
an unavailable shape is silently used. Section 3 proves they all pass
through every sufficiently large integer length.

When they pass, pack ALL N input bits as the payload. Enumerate the M
short seeds in increasing binary order, obtaining m-bit rows Y_i.
Their indices have s bits; split them into l-bit blocks and pad the
last block on the right. Starting with the first d bits of Y_i,
perform Xrefresh, Yrefresh and the h-slot alternating Xalt/Qalt chain
at each index block. Select the response at its digit's slot; Ynext
refreshes the next stage's d-bit seed from the ORIGINAL Y_i. The first
bit of the final selected response is coordinate i of an M-bit word.
The code fixes every truncation, ordering, field and fallback.

Apply resilient.md's monotone polynomial-graph function to that word.
Its fixed bad-coordinate parameter is 1/4096 and graph degree is
16384. Its field size q is the largest power of two below
floor(M/ceil(log_2(M+2))); its coordinate grid, selected graph offsets,
rational approximation to ln 2 and row ordering are fixed by the code.
The ignored suffix is explicitly accounted for in that component proof.
No good-set information or source-dependent advice is provided.

## 2. Uniform P and NP eligibility

All metadata operations terminate. The largest numerical constants above
are fixed integers. Computing L^25 and its constant shift uses a number
of bits polynomial in log N and a fixed constant; comparison against N
is effective. Once k<=N, m,h,ell,d,Qlen and code parameters are bounded
by fixed powers of N. A failed guard returns immediately.

Small-field searches have degree O(log N), and hence search spaces
polynomial in N. The degree-64 extension search has at most
q2^64<=2^64*N^8192 candidates. Generator factorization and enumeration
also take a fixed polynomial in N. This runtime is charged to original
N, not to the smaller condensed ambient length. seeded.md proves search
existence and the exact tests. Large code fields use deterministic
quadratic towers with polynomial bit arithmetic, not exhaustive searches
at their larger degrees. Explicit stacks avoid hidden recursion limits.

The short row enumeration has M<=N^16384 entries. Every row has b<=s
stages, h alternating slots per stage, and polynomially computable
linear code bits. The final graph map enumerates q^16384 rows with
q<=M, polynomially many graph blocks and O(log M) coordinates per block.
Every exponent and constant is fixed; the entire evaluator has time
at most c0*(N+1)^c1 for fixed finite c0,c1. Its polynomial exponent
and nontrivial onset are enormous, and no practical implementation at
nontrivial asymptotic lengths is claimed.

Define V(x,w)=1 exactly when w is empty and this total evaluator returns
one. Members have that empty witness, and acceptance of an allowed
witness implies membership by definition. This proves both directions
with witness bound zero. Witnesses, fields and internal enumerated seeds
are not extra inputs to the target circuits. No per-length advice,
truth-table oracle, hard table or sparse-length choice occurs.

## 3. Every sufficiently large integer passes every guard

All limits here are through ALL integers N. Then L~log_2 N,
J~log_2 L, k~2^10000*L^25/J and eventually k<=N. Floors change only
factors tending to one. We have

    m=Theta(L^23/J), h=Theta(L^2), ell=Theta(L^4),
    Qlen=Theta(L^8), d=Theta(L^20), E=2^200*L^2,
    s=Theta(L), b=Theta(L/log L).

The hidden constants are fixed. The following explicit bounds pay the
equal-power comparison B_0=O(L^23/log L) versus m=Theta(L^23/J).

Since k<=N, w1<=17log_2 N+1, w2<=128log_2 N+1. Therefore
s<=8337log_2 N+66<=2^15*L eventually, and M<=N^16384 eventually.
Also log_2 h>=2log_2 L, so b<=2^15*L/log_2 L eventually.
The rounding of h gives

    2^200*L^2<=h<=2^201*L^2.

Eventually h^2>=C_L. Then
Qlen<=32h^2*ell<=2^907*L^8. Also
h^3*ell<=2^1103*L^10, so d<=2^2272*L^20.
The sum d+Qlen+h*ell is at most 2d eventually. Consequently

    B_0=2^20*C_L*b*h*(d+Qlen+h*ell)
        <=2^6619*L^23/log_2 L.

On the other hand k>=2^9999*L^25/J and
m>=k/(2^301*L^2)>=2^9698*L^23/J eventually. Since
J<=2log_2 L eventually,

    m>=2^9697*L^23/log_2 L >4B_0.

The input joint budget has fixed slack as well:
h*m<=2^-99*k by the upper bound on h and the exact floor defining m.
Furthermore B_0/k=O(L^-2), so k>=4h*m+4B_0 eventually.
The remaining rank guards follow from m/d~constant*L^3/J tending to
infinity, and d having order L^20 while
4h*(Qlen+(h+1)ell) has order at most L^10. All constants are fixed.

For the short polynomial encoder, its coefficient capacity has order
k^2 log N, exceeding n'=k*w1. Its reconstruction advice is
O(m*k^(1/32)*log N). Dividing by k gives order
L^(-7/32)*J^(-1/32), up to fixed factors, which tends to zero.
Thus the exact advice guard passes. Its q-ary reconstruction errors
remain smaller than 1/(16q2^64): m<=N, the encoded degree is at most
N^(1/32), and rho=N^-8/(2m)>=N^-9/2 give per-curve bounds
O(N^(-119*128)) and O(N^(27-(4095/32)*128)). Both are much smaller
than N^-8192. seeded.md supplies the full counting and hybrid argument.

For every ordinary long call, requested output bit length is O(log L),
so its exact code w_target<=3E eventually. Code-field rounding gives
w<=6E; its index length is at most 12E, design-field size at most 24E,
and seed D<=288E^2<2^409*L^4. The supplied ell=2^500*L^4 pays D by
a fixed factor exceeding 2^91. The design capacity is at least a fixed
positive multiple of L^8194, exceeding every requested output, whose
largest order is L^20. Its fixed degree 4096 is eventually below the
code index length. The entropy threshold's dominant term is 2^4100
times output, D is far below the smallest output ell, and all other
terms are lower order. Therefore H_call<=C_L times output eventually.

For Xrefresh, E_1=(2^20+o(1))*h^3*ell. Its w_target<=4E_1 eventually,
so D<=512E_1^2<=2^51*h^6*ell^2. The supplied d>=2^64*h^6*ell^2
pays it by a fixed factor of at least 8192. Its actual larger entropy
threshold is charged to C_L*d. Its graph capacity also has ample
power slack. All exact internal contracts therefore pass.

The resilient shape has q=Theta(M/log M), z=Theta(log M),
w=Theta(log M) and q_0=Theta(q), so its basic shape guards pass.
For its actual integer moment order R=2^20*(16384z+1)+4 we have
z,w<=s<=2^15*L and R<=2^50*L eventually. Thus

    Rw<=2^65*L^2 < h.

Its exact coded error charge is

    20+2ceil_log2(b+1)+4ceil_log2(h+1)+ceil_log2(R+1)
    +R*(1+ceil_log2(q_0+1))+10*16384*z+100.

Since q_0<=2^z, its R term is at most 2^66*L^2. The remaining terms
are eventually below 2^40*L^2. Hence the whole charge is at most
2^67*L^2<E. These explicit fixed-factor bounds replace any unsupported
assertion that equal logarithmic exponents have a strict power margin.

There are finitely many guards, each covered above by a fixed positive
factor or a strictly lower growth order. Therefore one fixed finite
onset makes every guard pass at EVERY integer N beyond it. No field-
compatible length, prime-degree subsequence or sampled input is used.

## 4. Universal affine dispersal

Let X be uniform on ANY affine subspace of F_2^N of dimension at least
k. At sufficiently large N all source-independent guards pass.
seeded.md proves strong error at most 2N^-8 for Short, which is linear
for every fixed seed. On an affine source, a rank-deficient output
image has distance at least 1/2 from uniform. Thus a FIXED fraction at
least 1-4N^-8 of its M rows are exactly uniform m-bit images. Their
good set depends on the source direction, not its sample or offset.

correlation.md proves for every at-most-h subset of good indices that
their selected ell-bit responses are close to independent uniform
rows, with joint error at most

    gamma=2^20*b^2*h^4*2^-E.

It supplies the affine residual decomposition, seed-before-linear-output
ordering, deficient-seed bound, equal-prefix treatment, retained target
transcripts, and lexicographic first-difference argument. Xrefresh pays
its larger entropy threshold; no near-full extraction is assumed. The
final bounded-marginal assertion does not assert that all M rows are
independent or that their full nonlinear matrix is an affine source.

The remaining at most 4M*N^-8 coordinates may be arbitrary functions
of X. Since M<=N^16384, M^(-1/4096)>=N^-4; hence eventually
4N^-8<=M^(-1/4096). Their count is therefore at most M^(1-1/4096).
The fixed-factor rank and error guards just proved give precisely the
bounded-marginal hypotheses of resilient.md. That proof uses monotone
zero/one extremes, graph packing, explicit two-row moments and bounded
inclusion-exclusion; it allows arbitrary dependence of the bad bits.
Both output values have positive probability. Therefore this SAME F_N
is nonconstant on EVERY affine subspace of dimension at least k.
The source and affine offset were arbitrary, and no finite experiment
is used to infer this universal statement.

## 5. Gate count, exact length and full remainder

Apply elimination.md with n=N and d=k-1, a nonnegative integer at all
sufficiently large lengths. Its complete full-B2 argument gives

    C_B2(F_N)>=4N-8k.

It permits every binary truth table, unrestricted sharing and fan-out,
constant and unary effects, repeated pins and single terminal outputs.
It converts W=s+4I+2q back to original gate count; it does not confuse
a potential with gates or presume direct-sum additivity.

Every N bit is actual payload. No witness input, certificate prefix,
field encoding prefix or padded replacement length is charged. Since
L<=log_2 N and L+2>log_2 N, we have J>=log_2 log_2 N whenever that
quantity is positive. Thus the EXACT definition implies

    8k<=2^10003*(log_2 N)^25/(log_2 log_2 N).

Take fixed K=2^10003 and a common finite onset from Section 3. This
proves the complete declared version-2 Big-O guarantee through all
integer lengths. No bounded rational coefficient for this very large
K is claimed. The loss divided by N tends to zero, yielding the same
joint projection (4,1). The particular loss is proved directly rather
than inferred from the leading liminf.

## 6. Finite evidence and limitations

family.py's declared audit takes exactly {}. It checks bounded integer,
field, orbit, rank, graph, covariance and Bonferroni predicates, and
actual small-length zero guards. No nontrivial asymptotic instance is
claimed to have run. The universal source, transcript, resilient-output,
all-length and gate-elimination assertions are mathematical arguments
in the supplied evidence. A finite check or prior review cannot replace
them. The adaptive transcript induction is especially intricate and
requires independent scrutiny; neither a local audit nor a simulation
is an authority to accept it.
