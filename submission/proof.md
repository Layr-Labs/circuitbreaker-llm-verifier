# One fixed family with lower bound 4N-O((log N)^13/loglog N)

For the single total language defined by family.py, there exist fixed
K>0 and N0 such that for EVERY integer actual input length N>=N0,

    C_B2(F_N)>=4N-K*(log_2 N)^13/(log_2 log_2 N).

The evaluator is deterministic polynomial time, with an enormous fixed
degree and onset, hence gives an empty-witness NP verifier. Both
structured score claims assert this SAME full version-2 Big-O guarantee
and its joint finite projection (4,1). All material component arguments
are supplied. A finite experiment or previous verdict is not a premise.

## Provenance and the new component

The fair gate kernel, indexed correlation proof and resilient graph
argument derive from fbadbf8d-fdb1-495b-a11f-2aba6ad4a9d2, source commit
41f5c012ad9745b26125b387924c77861d873995, with their earlier sources
preserved in the manuscripts. Their authorship is not claimed here.
No prior scientific score or attribution decision authorizes this proof.
The supplied state kernel retains 4n-8(d+1) in the full B2 model.

The new seeded component replaces the previous full-field primitive
orbit and code-bit long extractors. It uses an element of only bounded
large order, found by a root-count search polynomial in ambient length,
and a diagonal orbit whose values form a Vandermonde information set.
Thus a full primitive-element search in a huge field is unnecessary.
A deterministic quadratic tower represents that field efficiently.
Its seed length is linear in error exponent plus logarithmic ambient
length; the entropy requirement is paid explicitly, not mislabeled
constant-fraction extraction. seeded.md supplies the complete condenser,
reconstruction, hybrid, field and runtime arguments.

These shorter seeds allow both h and ell to have order (log N)^2,
and d order (log N)^8. The transcript budget has order
(log N)^11/loglog N; the initial joint-row entropy budget costs another
factor h. Explicit constants pay the resulting borderline comparisons.
The changed extractor thresholds, all metadata guards and resilience
parameter are accounted for below and in correlation.md.

## 1. Exact total family

For N<2 return zero. Otherwise put

    L=floor(log_2 N), J=ceil(log_2(L+2)),
    k=floor(2^10000*L^13/J).

All are exact integer computations. Return zero if k>N. Otherwise set

    m=floor(k/(2^300*L^2)),
    h=the least power of two >=2^200*L^2,
    ell=2^300*L^2, C_L=2^4110,
    Qlen=64h^2*ell,
    d=2^64*(h^3*ell+L),
    E=2^200*L^2, E_1=2^20*h^3*ell+E.

Short is seeded.md's folded/diagonal extractor at ambient N, threshold
k, output m and error exponent E_short=16ceil(log_2(N+2)). Its field
width is the least power of two >=128*(E_short+bitlength(N)+
bitlength(m+2)+10). Its seed has s=130w bits. Put M=2^s,
l=log_2 h and b=ceil(s/l). Each internal call uses the same extractor
with its requested output/error and the exact threshold H_call from
seeded.md. Guards check that threshold, its advice and coefficient
capacity, available source entropy, seed lengths and the full rank
budgets. Failure of any listed FamilyShape.build condition returns zero.
The guard M<=N^(2^20) is part of the total definition.

When all guards pass, pack ALL N input bits as payload. Enumerate the
M Short seeds in increasing binary order to get m-bit rows Y_i.
Split each s-bit index into l-bit advice blocks, padding the last on
the right. Starting with the first d bits of Y_i, at each block run
Xrefresh, Yrefresh and the h-slot Xalt/Qalt alternating chain. Select
the response at that block's digit; Ynext refreshes the d-bit next seed
from the ORIGINAL Y_i. Use the first bit of the final selected response
as word coordinate i. All ordering and masks are fixed by family.py.

Apply resilient.md's monotone polynomial-graph map to that M-bit word,
with fixed delta=1/2^20 and graph degree a=4/delta=2^22. Its exact
field, rational ln2 approximation, grid, selected offsets and row order
are determined uniformly. Its ignored suffix is explicitly accounted for.
Its integer moment/order/error guards are checked before enumeration.
No good set, source direction, witness, advice or certificate prefix is
an input to this evaluator or to the target circuit.

## 2. Uniform polynomial runtime and NP membership

Every branch terminates. Numerical exponents and constants are fixed.
Metadata uses polynomially many bits and operations in actual N; guarded
parameters are bounded by fixed powers of N. A failed guard returns zero.
For each extractor call the exact entropy guard gives H<=ambient A.
Its bounded-order search uses U=max(A,h_SU^128)<=A^2, at most U^2+1
candidate field elements and U multiplications each: polynomial in A.
The field may have superpolynomial cardinality, but its tower degree
w is polynomially bounded in N and arithmetic has polynomial bit cost.
There is no enumeration of that field, no factorization of its group
order and no search for a transitive primitive orbit.

Every actual extractor call is polynomial in its ambient length, output
and w, hence polynomial in N here. Reconstruction is only a counting
proof, not an evaluator routine. Code stacks avoid unbounded host recursion.
The M<=N^(2^20) row enumeration has b<=s stages and h polynomially many
slots each. The final graph map enumerates q^(2^22) rows with q<=M,
polynomially many blocks per row and O(log M) positions per block.
Its small graph field has degree O(log M)=O(log N), with a fixed large
constant, so its exhaustive irreducible search is also polynomial in N.
Thus fixed finite c0,c1 bound time by c0*(N+1)^c1 on ALL inputs.

The degree and onset are enormous. No nontrivial asymptotic instance
is claimed practical or executed. Define V(x,w)=1 exactly when w is
empty and this total evaluator outputs one. A member has that empty
witness; an accepted permitted witness implies membership. This proves
both directions with witness bound zero. The same fixed language and
algorithms work at every length. Witnesses and internal seeds are not
circuit inputs, and no per-length hard table is chosen.

## 3. Every sufficiently large integer passes every guard

Let N range through ALL integers. Then L~log_2 N, J~log_2 L,
k~2^10000*L^13/J and eventually k<=N. In orders of L,

    m=Theta(L^11/J), h=Theta(L^2), ell=Theta(L^2),
    Qlen=Theta(L^6), d=Theta(L^8), E=Theta(L^2), E_1=Theta(L^8).

Floors change only factors tending to one; all hidden constants are
fixed. The following explicit bounds pay the borderline comparisons.

Short has m<=N and E_short<=16(log_2 N+2). Its width satisfies
w<=256*(E_short+bitlength(N)+bitlength(m+2)+10), hence
s=130w<=2^20 log_2 N eventually. Thus M<=N^(2^20), and
s<=2^21 L eventually. Since log_2 h>=2log_2 L,
b<=2^21 L/log_2 L eventually. Also

    2^200 L^2<=h<=2^201 L^2,
    Qlen<=2^708 L^6,
    h^3*ell<=2^903 L^8, d<=2^968 L^8.

The sum d+Qlen+h*ell is at most 2d eventually. Therefore

    B_0=2^20*C_L*b*h*(d+Qlen+h*ell)
        <=2^5321*L^11/log_2 L.

On the other hand k>=2^9999 L^13/J and eventually
m>=k/(2^301 L^2)>=2^9698 L^11/J. With J<=2log_2 L,

    m>=2^9697 L^11/log_2 L.

Thus m exceeds 4B_0 with a huge fixed factor of slack. Moreover
h*m<=2^-99 k exactly by the upper bound on h and the floor defining m.
The direct paired-seed defect bound also has fixed slack:
Qlen+(h+1)ell<=66h^2 ell for h>=1, hence
4h*(Qlen+(h+1)ell)<=264h^3 ell<d/2^55. Eventually h^2>=C_L,
and m/d has order L^3/J and tends to infinity.

For Short's encoding, its dimension 128 and degree floor(k^(1/64))
give coefficient capacity Theta(k^2), exceeding its k field coefficients.
Its advice is O(m*k^(1/64)*w). Relative to k this has order
L^(-51/64)*J^(-1/64): w=Theta(L) and m/k=Theta(L^-2).
Thus advice/k tends to zero. Its additional error/log term is O(L),
also o(k). All exact Short capacity and advice guards pass eventually.
The strong Short error is <=N^-16, hence <=2N^-8 as used below.

For ordinary internal calls w_call<=512E eventually, because requested
output bit length is O(log L) and ambient bit length is at most O(L).
Their seed length 130w_call is below 2^17 E; the supplied
ell=2^300 L^2 pays this by a fixed factor exceeding 2^83.
For Xrefresh, w_call<=512E_1 eventually, so its seed is below
2^17 E_1<=2^38 h^3 ell. The supplied d>=2^64 h^3 ell pays it
by a fixed factor at least 2^26.

The exact entropy thresholds H_call=ceil((2^20*output*w_call)^(33/32))
have growth powers

    H_Xrefresh,H_Ynext=Theta(L^(165/16)),
    H_Yrefresh=Theta(L^(33/4)),
    H_Xalt,H_Qalt=Theta(L^(33/8)).

Each lies below its ambient length: 165/16<11 for Ynext,
33/4<11 for Yrefresh, 33/8<6 for Qalt, and both X thresholds
are below 13. Their coefficient capacities grow as Theta(H_call^2),
and the exact advice guards pass by seeded.md's
(33/32)*(63/64)=2079/2048>1 estimate.
In particular H_Y=max(H_Yrefresh,H_Ynext)=o(m),
H_X=max(H_Xrefresh,H_Xalt)=o(k), and
H_alt=max(H_Xalt,H_Qalt)=o(h^2 ell). Combining these with the
fixed-factor B_0/m slack and h*m/k bound proves

    m>=4B_0+4H_Y, k>=4h*m+4B_0+4H_X,
    Qlen>=8h^2 ell+8H_alt,
    k-h*m-B_0>=4h^2 ell+4H_alt.

These are the actual guarded source budgets, not an unproved claim
that the new long extractor needs only O(output) entropy. All supplied
seeds, exact advice, capacity and entropy conditions therefore pass.

The resilient shape has q=Theta(M/log M), z,w_graph=Theta(log M),
and q_0=Theta(q), so its basic shape guards eventually pass.
Use a=2^22. Since z,w_graph<=s<=2^21 L, its integer moment order
R=2^20*(a z+1)+4 is at most 2^64 L eventually. Consequently

    R*w_graph<=2^85 L^2<h.

Its exact coded error charge is

    20+2ceil_log2(b+1)+4ceil_log2(h+1)+ceil_log2(R+1)
    +R*(1+ceil_log2(q_0+1))+10a z+100.

The R term is at most 2^86 L^2 since q_0<=2^z. The other terms
are eventually below 2^50 L^2. Thus the full charge is at most
2^87 L^2<E. Equal logarithmic powers are paid by these explicit
fixed factors, rather than by a nonexistent strict exponent margin.

There are finitely many guards. Each is covered by a fixed positive
factor or strictly smaller growth order. Thus one common fixed finite
onset makes all of them pass at EVERY integer N beyond it, including
all rounding effects. No special lengths, sampled sources or compatible
field-degree subsequence are used.

## 4. Universal affine dispersal

Take ANY affine source X on F_2^N of dimension at least k. At the
onset from Section 3 all guards pass independently of this source.
seeded.md proves Short strong error <=2N^-8 and fixed-seed linearity.
A rank-deficient affine image has distance >=1/2 from uniform, so a
FIXED fraction at least 1-4N^-8 of the M rows are exactly uniform m-bit
images. The good set depends only on the source direction, not its
sample or affine shift.

correlation.md gives for every at-most-h subset of those good indices
approximately independent selected ell-bit responses, with error

    gamma=2^20*b^2*h^4*2^-E.

Its affine decomposition, source-direction rank test, deficient-seed
bound and lexicographic induction are reproduced. The revised exact
H_call thresholds and guarded budgets pay each use. Equal-prefix rows
may coincide; their simultaneous responses are excluded from the
independence assertion. All seeds are fixed before linear data
observations; the target's own previous transcript is charged. No
claim of independence of the full nonlinear M-row matrix is made.

The other at most 4M*N^-8 coordinates may depend arbitrarily on X.
Since M<=N^(2^20), M^(-1/2^20)>=N^-1. Thus eventually
4N^-8<=M^(-1/2^20), so their number is at most M^(1-delta),
delta=1/2^20. The rank and error guards just proved meet precisely
the hypotheses in resilient.md. That proof permits non-oblivious bad
bits, using monotone zero/one extremes, graph packing, explicit
negative correlation, two-row moments and truncated inclusion-exclusion.
Both final outputs have positive probability. Therefore this SAME F_N
is nonconstant on EVERY affine subspace of dimension at least k.
The source, direction and affine offset were arbitrary.

## 5. Original full-B2 gates and exact remainder

Apply elimination.md at n=N,d=k-1, eventually a nonnegative integer.
Its complete kernel gives C_B2(F_N)>=4N-8k. It handles every binary
truth table, constant and unary effects, repeated pins, unrestricted
sharing, safe paired domains and terminal output cases. Its charged
potential is telescoped back to original gates; no direct-sum or
restricted-basis lower bound is assumed.

All N bits are actual payload. No witness, encoding prefix, field
certificate or padded replacement length is counted. As L<=log_2 N
and L+2>log_2 N, J>=log_2 log_2 N whenever positive. Hence

    8k<=2^10003*(log_2 N)^13/(log_2 log_2 N).

Take fixed K=2^10003 and the common onset. This proves the entire
version-2 Big-O loss through all integer lengths; no bounded rational
coefficient for this enormous K is asserted. The loss divided by N
tends to zero, giving the SAME joint projection (4,1). The specific
remainder is proved directly, never inferred from the leading liminf.

## 6. Computational evidence and limits

The declared family.py audit takes {}. It checks bounded integer,
field, orbit, folded-rank, graph, covariance and Bonferroni predicates,
plus bounded-order searches and small diagonal Vandermonde matrices.
Its actual full-family calls exercise only the explicit zero branch.
The older primitive-orbit toy has degree16, while the actual new
extractor has 128 diagonal coordinates; neither finite check is a
certification of the universal source theorem. No nontrivial asymptotic
instance is claimed run or feasible. The adaptive transcript induction
remains especially intricate and requires independent scrutiny. A
simulation, prior verdict or finite audit is not mathematical authority.
