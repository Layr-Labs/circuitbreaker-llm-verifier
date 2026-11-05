# One fixed uniform family with lower bound 4N-O((log N)^32)

## Claim, model and evidence

Put C=32, a fixed positive integer. The evaluator family.py
defines one total Boolean language L in P and therefore NP. For every
sufficiently large integer ACTUAL input length N its characteristic
function satisfies

    C_B2(F_N)>=4N-O((log_2 N)^C).

One fixed Big-O constant and one fixed onset work for this SAME family.
The full B2 model includes every binary truth table, constants, unary
functions, repeated pins, acyclic unrestricted sharing and one output.
Both structured source claims assert the same complete bound.json
expression and its joint finite projection (4,1), with identical evidence.
There is no separate function proving each coordinate.

The seven evidence files reproduce the complete affine disperser and
full gate-elimination proofs. The evaluator and its finite-audit entry
point share family.py to stay within the eight-file profile. The audit
checks only bounded predicates and zero-guard examples; no nontrivial
asymptotic instance is claimed to have been executed. The polynomial
runtime exponent and the nontrivial onset are enormous.

The source builds on the submitted 36f8d801-ded3-4b75-93e7-d5ee4c7e0ae5
for deterministic fields, linear code-bit extraction and affine residual
facts, and on 11e8b36d-cbd4-44f6-810f-7a7f875b4ad0 for the complete
4n-8(d+1) gate kernel. Its state framework traces to production
237016f9-ece3-4a05-9933-dfcf6056e787. All material is supplied anew;
prior scientific verdicts and attribution decisions are not premises.

The new construction replaces the polynomial-entropy tree by indexed
alternating extraction, organized after Xin Li, ECCC TR15-121, Sections
3--4. Its two otherwise advanced components are supported here: a folded
affine condenser followed by the basic Shaltiel--Umans orbit extractor
on the SMALL condensed ambient space (seeded.md), and an explicit
monotone polynomial-graph resilient disperser (resilient.md). The latter
uses only a proved covariance identity and two-row moments, rather than
an external non-oblivious-source extractor theorem. The entire adapted
conditioning and indexed separation argument is in correlation.md.

## Parameter redesign and provenance

The construction and full component proofs derive from
fbadbf8d-fdb1-495b-a11f-2aba6ad4a9d2, source commit
41f5c012ad9745b26125b387924c77861d873995. Its authorship and prior
judgments are not claimed as premises. This candidate changes five
fixed exponents together, retaining the exact metadata guards and
all primitive algorithms. The parameter choice is

    C=32, m=floor(k^(4/5)), h=pow2ceil(ceil(k^(1/15))),
    ell=ceil(k^(2/15)), E=ceil(ell^(12/25)).

Every rational power is implemented by an exact integer root of an
integer power. The short SU polynomial degree stays floor(k^(1/8)).
Increasing its output length m changes its charged advice bound to
O(k^(37/40)*log N), with total k-exponent 153/160<1. Its existing
reconstruction remains valid: m<=N^(4/5), rho>=N^(-44/5)/2. With
q>=N^128 the per-curve error bounds are O(N^(-119.2*128)) and
O(N^(26.4-127.875*128)), both much smaller than N^-2048.
The same seed length and strong short-extractor error therefore hold.

For the long calls the low-error code/design seed has order E^2,
which is ell^(24/25)=o(ell). Its graph capacity is far larger than
any requested output. The exceptionally low-error Xrefresh still
uses E_1=2^20*h^3*ell+E and d=2^64*(h^3*ell+logceil)^2;
the identical fixed-factor seed and entropy accounting applies.

The new full rank budget has k-exponent
7/15+4/15+1/32=367/480<4/5. The residual input also suffices since
h*m has exponent 13/15<1. At the final resilience step,
h=Theta((log N)^(32/15)) and E=Theta((log N)^(256/125));
both logarithmic exponents strictly exceed two. Thus h pays the
bounded-marginal order and E pays the complete logarithmic error
charge. The full all-guard analysis is below. The language changes
with these parameters; a smaller remainder is not asserted for the
old family by merely changing a score declaration.

## 1. The exact family

The code fixes every field, ordering, bit position, rounding, index and
fallback. Bit zero is stored in integer position zero. At length N<2
return zero. Compute

    L_N=ceil(log_2(N+2)), k=L_N^C

by capped integer exponentiation, returning zero if k>N. Capping computes
min(L_N^C,N+1) exactly, so this branch is a mathematical guard rather
than a hardware cutoff. In the other branch set

    m=floor(k^(4/5)),
    h=the least power of two >=ceil(k^(1/15)),
    ell=ceil(k^(2/15)), C_L=2^4110,
    Qlen=(16h^2+16C_L)*ell,
    d=2^64*(h^3*ell+L_N)^2,
    E=ceil(ell^(12/25)), E_1=2^20*h^3*ell+E.

All powers here are exact integer powers or integer roots. C and C_L
are fixed construction constants, not functions of N.

The initial short extractor uses q_1=2^w1>=N^16*k, with the least such
binary exponent w1, and q_2=2^w2>=N^128, likewise minimal. Its seed has

    s=w1+17w2 bits; M=2^s; l=log_2 h; b=ceil(s/l).

The numerical guards require the exact coefficient capacity for the
sixteen-variable encoding, its advice/entropy inequality, M<=N^4096,
the rank budgets in correlation.md, and every actual long-extractor
capacity, seed and entropy inequality. All metadata computations and
all conditions are listed in FamilyShape.build; their failures return
zero. The resulting function is defined on every binary input.

If all guards pass, pack ALL N input bits into the payload integer.
For each of the M short seeds in increasing order, compute its m-bit
linear image Y_i. Split the s-bit binary row index into l-bit advice
blocks, padding the last on the right. At each advice block perform
Xrefresh on the current d-bit seed, Yrefresh on the original Y_i,
then the h-step alternating Xalt/Qalt sequence. Choose the response at
the position specified by that advice block. Ynext refreshes the d-bit
seed for the next block. The final selected ell-bit response is W_i.
Take its first bit as coordinate i of an M-bit word. No rows are
selected using the source's unknown entropy or good set.

Apply the monotone polynomial-graph function in resilient.md to this
word. Its field size q is the largest power of two at most
floor(M/ceil(log_2(M+2))). Its row polynomials have degree at most the
fixed a=16384 and constant term zero. Its graph size w, rational ln2
approximation and selected offset count q_0 are computed exactly by
ResilientShape.build. Each row ORs its graph ANDs; the result ANDs all
rows. This order is exactly its code. The ignored suffix beyond w*q
coordinates is part of the specified map and the proof accounts for it.

The additional guards h>=R*w and E>=the computed error charge use
R=2^20*(a*z+1)+4, where z=log_2 q. They conservatively imply the bounded
marginal and error requirements in resilient.md. They are computed
before evaluating the enormous M rows. No oracle, advice, certificate
prefix, probabilistic choice or input-independent hard table appears.

## 2. Uniform polynomial-time evaluation and NP membership

Each guarded branch terminates. Capped powering uses O(log N)-bit
intermediates and a fixed number of exponent bits. Once k<=N, all h,
ell, d, Qlen, code degrees and output lengths are bounded by fixed
polynomials in N. Moreover w1<=17log_2N+1 and w2<=128log_2N+1, so
s<=2193log_2N+18. The explicit guard M<=N^4096 is eventually true and
bounds every row enumeration even at an exceptional nonzero length.

All binary field searches have degree O(log N), hence polynomial search
space in N. The short extractor's degree-16 extension is over F_(q2):
its search has at most q2^16=O(N^2048) candidates, a fixed polynomial.
Finding its multiplicative generator by factoring the finite order and
exhaustive element tests also costs a fixed polynomial in N. This step
would not be polynomial in the smaller condensed ambient length; the
runtime is honestly charged to the ORIGINAL N. seeded.md proves search
existence and the exact irreducibility and generator tests.

Larger code fields have power-of-two degrees and use explicit quadratic
towers, without exponential irreducible search. Their multiplication
uses an evaluation stack with polynomial bit cost. Polynomial-graph
design fields have degree O(log N), so the small exhaustive searches
are appropriate there. The coefficient enumeration has fixed dimension
sixteen and an explicit stack. No mathematical recursion relies on a
host recursion limit.

Each code bit, folded evaluation, multivariate evaluation and alternating
response is computable in polynomial time. There are M<=N^4096 short
rows, b<=s stages per row, and h polynomially many alternating slots.
The final resilient graph has q^a rows with fixed a=16384, q<=M,
polynomially many selected graph blocks per row, and O(log M) points
per block. Both its evaluation and memory are bounded by fixed powers
of N, albeit enormous ones. No simulation of the unbounded proof
reconstruction in seeded.md occurs during evaluation.

Thus fixed finite c_0,c_1 bound the full running time by c_0(N+1)^c_1
on ALL inputs, including every exceptional guard branch. Define the
NP verifier to accept exactly when the witness is empty and the total
evaluator returns one. A member has that permitted empty witness;
an accepted permitted witness implies membership by definition. Both
directions use this same uniform algorithm. Witnesses and internal
extractor seeds are not extra inputs to the target B2 circuit.

## 3. Every sufficiently large integer passes every guard

Let N tend through ALL positive integers. Then k=(log_2N)^C(1+o(1)),
and eventually k<=N. In powers of k the relevant lengths are

    m=k^(4/5+o(1)), h=Theta(k^(1/15)),
    ell=k^(2/15+o(1)), d=Theta(k^(2/3)),
    Qlen=Theta(k^(4/15)),
    b=Theta(log N/log log N)=k^(1/32+o(1)),
    B_0=k^(367/480+o(1)).

The power-of-two rounding of h changes only fixed factors. The logceil
term in d is negligible because 1/32<1/3. Since 367/480<4/5,
B_0=o(m). Also h*m=o(k), h^2 eventually exceeds C_L,
d exceeds 4h*(Qlen+(h+1)ell), and m/d tends to infinity.
The short polynomial coefficient capacity has order k^2 log N and
exceeds its k*w1 input bits. The advice is O(k^(37/40)log N)=o(k),
since 37/40+1/32=153/160<1. All short guards eventually pass.

For a low-error long call, E=Theta(k^(8/125)) dominates log N and
logarithms of all requested output lengths. Its actual code index
length is Theta(E), within fixed rounding factors, and its design
seed length D is Theta(E^2)=Theta(k^(16/125))=o(ell), since
16/125<2/15. With design order 4096, graph capacity is at least a
fixed positive multiple of k^(32776/125), far exceeding the largest
requested output scale k^(2/3). Eventually the graph degree is below
the evaluation length. The entropy threshold is dominated by
2^(4096+4) times the output scale, with fixed spare factors in C_L.
Thus all exact low-error contracts eventually pass.

For Xrefresh, E_1=(2^20+o(1))*h^3*ell. Eventually its w_target is at
most 4E_1; rounding the field degree and graph universe gives an actual
seed D<=512E_1^2<=2^51h^6ell^2. But d>=2^64h^6ell^2. This pays the
longer seed with a fixed factor of at least 8192. Its entropy threshold
is bounded by C_L*d, and the graph capacity again has a strict power
margin. It is not asserted to have the low-error row's smaller entropy
budget. All these checks concern the ACTUAL coded seeds and outputs.

The resilient metadata has q=Theta(M/log M) within fixed factors,
z=Theta(log M), w=Theta(log M), and q_0=Theta(q). Its fixed shape
conditions therefore eventually hold. Its conservative R*w is
O(log(M)^2)=O(log(N)^2), whereas h has order (log N)^(32/15). The computed
error charge is also O(log(N)^2), while E has order (log N)^(256/125).
Both final guards therefore pass. The series approximation of ln2 and
the integer q_0 rounding have errors tending to zero at every length,
not just a special subsequence.

There are finitely many guards. All comparison powers and constants
are fixed, and every asymptotic inequality above has a strict exponent
or fixed-factor margin. Hence one finite onset works for all of them
and for all sufficiently large integer N. Field-size and code-degree
rounding stay within their displayed fixed factors. No prime degree,
input-certificate shape or selected-length assumption occurs.

## 4. Universal affine dispersal

Take any affine subspace of F_2^N of dimension at least k and its uniform
distribution X. All guards pass at sufficiently large N independently
of this source. seeded.md proves that Short has strong error at most
2N^-8 and is linear for each fixed seed. Therefore at least a fraction
1-4N^-8 of its M row seeds give EXACT uniform m-bit Y_i on this source.
Its good set is fixed by the source direction space.

correlation.md supplies the complete indexed alternating argument on
any at-most-h subset of good rows. It fixes all seeds before exposing
linear outputs, splits affine residuals into independent components,
and pays for the larger Xrefresh entropy requirement. Equal-prefix
rows need not be independent and can have equal responses. The first
differing advice block puts every smaller separated row's selected
response before the target's slot. Fresh extraction from the retained
original row supplies the next seed. The resulting first bits on good
rows are gamma-close to h-wise independent uniform bits, where

    gamma=2^20*b^2*h^4*2^-E.

This is a statement about all bounded marginals, not an assertion that
the entire row matrix is affine or independent. The remaining at most
4M*N^-8 bits can have arbitrary dependence on the good bits. Since
M<=N^4096, their count is eventually at most M^(1-1/4096).
The final rank/error guards imply precisely the sufficient bounded
marginal conditions in resilient.md. That manuscript handles these
non-oblivious coordinates by monotone zero/one extremes. Its polynomial
graph packing bound leaves almost all rows light; a proved conditional
covariance calculation and two-row inclusion-exclusion retain positive
mass for both output values. Thus both values occur on X.

The source, direction space and affine offset were arbitrary. Consequently
this one F_N is nonconstant on EVERY affine subspace of dimension at
least k. All entropy losses, code errors, retained components and source
observations are internal to its supplied proof; none is an assumption
about the unknown input source.

## 5. Exact loss assembly and the joint projection

Apply elimination.md with n=N and d=k-1, a nonnegative integer eventually.
It yields

    C_B2(F_N)>=4N-8k.

All N coordinates are actual payload. No field prefix has been fixed,
no unused witness bit has been charged to the target, and no padded
length has replaced N. For N>=4,
ceil(log_2(N+2))<=2log_2N. Therefore

    8k<=2^(C+3)*(log_2N)^C.

Take K=2^(C+3) and a finite common onset from the component arguments.
Both K and that onset are fixed independently of N. This establishes
exactly the declared version-2 Big-O logarithmic remainder, including
its all-length quantifier. An explicit numerical loss coefficient is
not claimed: the schema uses Big-O, as appropriate. Since
(log N)^C/N tends to zero, the same full guarantee implies the joint
projection (4,1). No inference from a leading liminf to this particular
remainder is made in the opposite direction.

## 6. Finite evidence and practical limits

The declared family.py audit takes exactly {}. It checks bounded
integer roots and capped powering, finite field extension identities
and primitive order, a small orbit polynomial encoder, folded affine
rank cases, polynomial graph intersections, truncated-cube covariance
and Bonferroni identities. Its full-family evaluations at small lengths
all return through the k>N guard. It does not evaluate a nontrivial
asymptotic instance or certify the universal lemmas. Those obligations
are carried by the supplied proofs.

The algorithm has one enormous fixed polynomial runtime exponent and
an enormous onset. It is a uniform asymptotic construction, not a
practical cryptographic extractor or a claim about feasible input sizes.
A successful finite audit, prior submission judgment, graph edge or
arithmetic checker is not a substitute for the complete mathematical
arguments and their configured scientific review.
