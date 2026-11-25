# A format-ceiling iterated-log bound for full-B2 element distinctness

## Main theorem and attribution

Fix R=10^256-1, the integer whose decimal expansion consists of 256 nines.
One fixed total language L in P satisfies, for every sufficiently large integer N,

    C_B2(L_N) >= R N (log_2(log_2 N))^R.

The coefficient and iterated-log exponent are both R; the remainder is empty.
R is a fixed finite integer, independent of N. Its exponent-one projection is

    lim_(N -> infinity, N integer) C_B2(L_N)/N = +infinity.

Equivalently, for every real A>0 there is N_A such that every integer
N>=N_A satisfies C_B2(L_N)>=A*N. Circuits are deterministic, acyclic,
single-output full B2, allowing all sixteen binary truth tables, arbitrary
depth and fan-out, repeated pins and shared subcircuits. Gates are counted
once; inputs, wires, output selection and substituted constant nodes are
free. N is the entire input bit length.

This is a Boolean-circuit corollary of OpenAI, *The matrix multiplication
exponent is nine fourths*, November 24, 2026 edition, Theorem 2.2 (label
`thm:many`). The immutable source snapshot is
`openai/math@92b516ef8e2abf66c3db7aa31438d3010bc7956f`:

https://github.com/openai/math/tree/92b516ef8e2abf66c3db7aa31438d3010bc7956f/preprints/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026

The tensor-network lower bound and its geometric/counting machinery belong
to the cited OpenAI authors. The fixed language, Boolean-network conversion,
padding proof and selection of supplied excerpts are adapted from Gal Arnon
(Yukon solver galarnon42), CircuitBreaker submission
e70d8290-d404-403a-8e48-c77118440d3a, immutable candidate commit
6d40928003a7a83b62a7bcfe34c8255a49cdf1e0:

https://github.com/Layr-Labs/circuitbreaker-llm-verifier/pull/942

The immediate starting package is our submission
c8205ffb-2702-48b8-8547-e8d8043fa084 (PR947), candidate commit
81d61c568f84c9de47b48358b21de7526f3e6d0f:

https://github.com/Layr-Labs/circuitbreaker-llm-verifier/pull/947

It extends our submission 449e2336-6ab9-4519-ad1e-22019c89d545 (PR946):

https://github.com/Layr-Labs/circuitbreaker-llm-verifier/pull/946

We retain their exact consumer accounting, output fusion, quadratic packing
and source-context restorations. This revision changes only the declared
fixed iterated-log consequence and its elementary limiting proof.
It supplies no new tensor lower-bound technique or numerical positive
logarithmic exponent. The finite R uses the full 256-digit rational format:
within p=1,q=0 it is the largest representable positive coefficient and
iterated-log exponent, since every positive denominator is at least one.
This encoding observation is not an optimality theorem for circuit complexity.
All mathematical premises are supplied for fresh review; no previous
judgment or platform status is a premise.

The supporting source material is supplied in source-*.md with original
file names, theorem labels, source provenance and attribution. Provenance
hashes identify the original versions, not an approval. Source-1 and
source-3 are unchanged from the cited CircuitBreaker source. Source-2
clearly marks our quadratic-packing adaptation of the original compressor
and preserves the complete remaining argument and its hypotheses. Source-4
restores the opening context verbatim from the same immutable upstream
model.tex, and source-5 explicitly restores obstruction.tex's standing
finite-field assumption. The source theorem supplies a fixed positive exponent
on log N. The declared iterated-log bound works for every such exponent.

## 1. One explicit total language

For each integer w>=1 let n_w=2*w*2^w. On any input x of length N:

1. If N<4, output 0.
2. Otherwise let w be the largest integer with n_w<=N.
3. Split the first n_w bits into 2^w consecutive records of 2w bits each.
4. Output 1 iff every pair of records is different. Ignore the suffix.

There is no advice or promise. In particular, the suffix can be arbitrary.
A fixed deterministic evaluator finds w by successive integer increments
and compares every pair of records bit by bit. There are at most N^2
pairs; even locating each pair by scanning the input uses polynomial time.
A conservative O((N+2)^4) multitape time bound suffices, including counter
and parameter arithmetic. Thus L is in P.

For NP take witness bound W(N)=0 and a fixed verifier V(x,z) that rejects
nonempty z and otherwise runs the evaluator. If x is in L its empty
witness is accepted. If V accepts, z is empty and the evaluator accepted x.
Both directions hold at all lengths, including the small rejecting cases.

## 2. The exact external theorem and its model

A bit tensor of arity u over F2 is a table {0,1}^u -> F2. Each network
vertex carries a table. Internal wires pair argument positions; each
unpaired position is a distinct named exterior bit. For fixed exterior
bits, the represented value is the sum over internal wire assignments of
the product of the table values. Every vertex of a trivalent network has
arity at most three. Volume counts vertices and all wires, including
exterior wires.

Let D_w:{0,1}^{n_w}->F2 be the distinctness indicator on 2^w records of
2w bits. The cited Theorem 2.2 states that for every fixed finite field k
there are sigma>0 and w_* such that, for all integers w>=w_* and
1<=t<=2^w,

    vol_k(D_w^(tensor t)) > t*n_w*w^sigma.

Sharing between the t lists is allowed. Specializing to k=F2 and t=1 gives

    vol_F2(D_w) > n_w*w^sigma.                         (1)

By decreasing this fixed exponent if necessary, assume 0<sigma<=1.
This preserves (1), since w>=1. No numerical lower bound on sigma is used.

There is no hypothesis omega_k<9/4 in this theorem. That assumption
belongs to the separate compression theorem (Theorem 2.1), which this
entry does not use. Section 5 of the source, `list.tex`, proves (1)
by an exhaustive ports/no-ports distinction, using the geometric and
counting results identified in the supplied evidence. Neither specific
case for characteristic two is assumed.

The Exact batches theorem in source-4 has the standing hypothesis that
the characteristic admits ports. This is explicit in the restored opening
of P/model.tex. Its constants may depend on the fixed characteristic and
port constants. M/list.tex invokes it only inside the ports branch, after
fixing the field. The no-ports branch instead uses the geometric
uniform-destinations argument. Thus the exhaustive two-case proof neither
drops the ports hypothesis nor requires it to hold in characteristic two.

## 2A. Quadratic packing and its precise quantitative consequence

Source-2 proves the packing lemma with K=4d^2 and Q=3d^2 for every
power-of-two base degree d>=1. Distinct nonzero linear forms on
F2^(log_2 d+2) replace one independent bit per base color. Any two
distinct such forms remain independent after arbitrary lane translations,
including when either endpoint's outgoing-slot label is fixed. The complete
retained proof counts exactly Q demands per bundle and dQ occurrences of
each surviving role, moves every whole old factor once in either common
layout, separates its bases by tag offsets, and preserves remaining-lane
translation and slot indexing. The binary tests concern lane indices and
work over every tensor ground field. Slot-list bounds depend on the fixed
output degree or subdivision counts as well as d.

Here is the full additional numerical regularization argument. The
six-party masses are 8/31, 4/31, 40/93 and 16/31, so their minimum is
4/31. With K=4d^2 the new slot count is

    m_e=floor(q_e*d'/(4d^3)).

For any fixed 0<delta_0<1 choose the next power of two d' at or above
31d^3/delta_0. Every unrounded count z is at least 1/delta_0>1; hence
floor(z)>=z-1>=(1-delta_0)z and m_e>0. Thus all source mass-loss
conditions hold. The common upper-degree calculation still holds because
Q=3K/4. Rounding the degree to a power of two gives

    d' < (62/delta_0)*d^3.

For example H'=125/64 and delta_0=1/256 satisfy the source's strict
31/16<H'<2 and give rho=H'/(2*(1-delta_0))=50/51. In freely chosen
intermediate steps the degree can therefore satisfy d'<15872*d^3.
For A=62/delta_0, iteration gives

    log_2 d_i
      < 3^i*log_2 d_0 + ((3^i-1)/2)*log_2 A       (i>=1).

This follows by induction from log_2 d_(i+1)<3log_2 d_i+log_2 A.
The last output degree must still be a usable order in the port family;
no upper bound on that next usable order is supplied here. The initial
port degree, uniform port constant, number of shrinking steps, motif maps
and some local-record bounds also remain nonnumerical. The new packing
preserves the full qualitative compressor and exact-descent proofs but
does not supply a numerical positive logarithmic exponent. Inequality
(1) and the remainder of this proof continue to use only a fixed sigma>0.

## 3. Exact consumer accounting in the Boolean conversion

**Lemma.** If a Boolean function depends essentially on each of its q>=3
inputs and has a full-B2 circuit with s gates, its predicate has an exact
trivalent F2 network of volume at most 5s-2q. For essential q=2 inputs,
volume at most 5s-2q+2 suffices. For general circuits, including constant
and unary outputs, volume at most 5s+2q+1 suffices.

Here an input is essential when changing only that bit changes the output
for at least one assignment. The zero-offset bound requires q>=3; it is
not asserted for a one-gate two-input function or a zero-gate primary
input function.

**Proof.** First suppose q>=1 and the computed function is nonconstant.
Substitute every free constant input pin into its binary truth table and
reconnect that now-ignored pin to x_1. This is a legal full-B2 gate at
unchanged cost. It preserves the function without adding a computational
basis. Trim to the output's ancestor subgraph after this normalization.
Write r<=s for its retained gate count and u for the number of primary
input signals having at least one consumer. Then u>=1. Count pin
occurrences separately even when a gate uses the same predecessor twice.

Use the ternary table T_g(a,b,c)=1 iff c=g(a,b) at each retained gate.
Temporarily attach a unary output test P(c)=c. A logical signal with d>=2
consumers is distributed by a binary tree of d-1 ternary equality tables
E(a,b,c)=1 iff a=b=c. If d=1, use a direct wire.

There are exactly 2r+1 consumer occurrences, including the output test.
Every retained gate-output signal has a consumer, because it is in the
output's ancestor subgraph. Exactly r+u logical signals therefore have
positive consumer count. Thus the number of equality vertices is exactly

    c = sum_(signals with d>=1) (d-1)
      = (2r+1)-(r+u) = r-u+1.

Connect each of the q-u unused primary inputs to a unary all-one tensor,
so every original input remains a distinct named exterior bit. With the
temporary output test, the number of vertices is

    v = r+c+(q-u)+1 = 2r+q-2u+2.

Gate and equality vertices are ternary, and the remaining vertices unary.
There are exactly q exterior wires. If i counts internal wires, the
incidence identity is

    2i+q = 3r+3c+(q-u)+1,
    i = 3r-2u+2.

Consequently this network has volume

    v+i+q = 5r+2q-4u+4.                              (2)

For every primary input assignment, topological evaluation determines all
gate signals and every equality copy uniquely. The constraints therefore
have exactly one internal assignment before the output test, and it
contributes one iff the circuit accepts. The contraction equals its
zero-one predicate over F2, without cancellation. Undirected cycles or
parallel wires arising from shared descendants and repeated pins do not
change this uniqueness.

If the output is a retained gate, its only consumer after trimming is the
output test. Indeed any other consumer cannot be an ancestor of that
terminal gate in an acyclic circuit. Contract its output wire with P:

    sum_(c in {0,1}) T_g(a,b,c)*P(c) = g(a,b).

Replace the terminal ternary gate table and the unary test by this binary
table. This removes one vertex and one internal wire, reducing (2) by two.
It is a tensor contraction, not an extra Boolean operation or an uncharged
Boolean basis. Every original retained Boolean gate is still counted.

If all q inputs are essential, each has to occur among the normalized
ancestors, so u=q. For q>=2 the output cannot be a primary input or a
constant. The terminal contraction therefore applies, and its volume is

    5r-2q+2 <= 5s-2q+2.

For q>=3 there is one more contraction. This network is connected:
every retained gate lies on a path to the output and every primary input
is used. It has no unused-input unary components. Every tensor is ternary
except the terminal binary tensor just formed. The latter has an internal
neighbor; otherwise both its legs would be exterior, and connectedness
would give at most two exterior bits. Let its neighbor be ternary and
let k>=1 be the number of joining internal wires. Contract all these
wires simultaneously. There are at most two because one tensor is binary.
The resulting tensor has arity 2+3-2k<=3. No named exterior index is
identified or summed. The contraction deletes one vertex and k internal
wires, reducing volume by 1+k>=2. Finite tensor contraction is associative,
so it preserves the entire predicate exactly, even when parallel wires
occur. It does not need a new uniqueness claim after fusion. Thus

    vol_F2(f) <= 5r-2q <= 5s-2q                    (q>=3).

The separate q=2 bound is necessary: one binary table has volume three,
whereas 5s-2q would equal one for s=1,q=2.

For a general nonconstant function, (2) and u>=1 give volume at most
5s+2q even without the last contraction. A primary-input output is covered
directly by its unary output table and q-1 all-one input tables, of volume
2q. A retained negated-projection gate is one ordinary counted B2 gate.
For a constant output use its scalar zero or one and q unary all-one
input tables, of volume 2q+1. This includes q=0. These cases establish the
general 5s+2q+1 bound as well. QED.

The construction permits all sixteen truth tables, arbitrary nonuniform
wiring, unbounded fan-out, shared subgraphs, repeated pins and unrestricted
depth. Only the direction Boolean circuit -> tensor network is used.

## 4. Essential inputs and the distinguished lengths

Every one of the n_w bits of D_w is essential. To see this, choose any
record and any bit within it. Set this record and one other record to
the all-zero word. Let e be the word having only the chosen bit equal to
one. Give the remaining 2^w-2 records distinct words avoiding both 0 and e.
This is possible because 2^(2w)-2>=2^w-2 for w>=1. Before the chosen bit
is flipped there is a duplicate zero record. After it is flipped, all
records are distinct. The output changes. This also covers w=1, when
there are no remaining records.

Thus q=n_w>=4 meets the sharper lemma's hypotheses. Let s_w=C_B2(D_w).
Combining it with (1), for all sufficiently large integer w,

    5s_w-2n_w > n_w*w^sigma,
    s_w > n_w*(w^sigma+2)/5.                         (3)

In particular s_w/n_w tends to infinity for one fixed sigma>0.

## 5. Every integer length and the declared fixed iterated logarithm

For N>=4 take the w selected by the evaluator. Maximality and the exact
ratio of consecutive distinguished lengths give

    n_w <= N < n_(w+1) = 2*(w+1)*n_w/w,
    n_w > N*w/(2*(w+1)).

Fix the ignored suffix of any circuit for L_N to zero. Its truth tables
can absorb these constants, with ignored pins reconnected to an active
input if necessary, without adding gates. It now computes D_w, so
C_B2(L_N)>=s_w. From (3),

    C_B2(L_N)
      > N*w*(w^sigma+2)/(10*(w+1))
      = N*w^sigma/10
        + N*(2w-w^sigma)/(10*(w+1)).                 (4)

We fixed sigma<=1, so w^sigma<=w for w>=1. The last term in (4)
is at least N*w/(10*(w+1))>=N/20. Therefore, for every sufficiently
large integer N, with its selected w,

    C_B2(L_N) > N*w^sigma/10 + N/20.                  (5)

This uses the positive linear term in (3) to cover the exact padding
ratio, with no terminal constant loss. It does not infer an empty
remainder from a bare liminf.

For large w, N<4n_w=8w*2^w gives
log_2 N < w+log_2 w+3 <= 3w. Let B_w=log_2(3w)>0. Then
0<log_2(log_2 N)<=B_w eventually. Equation (5) implies

    C_B2(L_N)/(R*N*(log_2(log_2 N))^R)
      > w^sigma/(10*R*B_w^R).                        (6)

Here is a complete elementary divergence proof for this enormous but
fixed integer R. For w>=3, ln(3w)<=2 ln w, so B_w<=2 ln w/ln 2.
Put t=ln w>0 and use the (R+1)-st term of the exponential series:

    w^sigma/(10*R*B_w^R)
      >= (ln 2)^R*exp(sigma*t)/(10*R*2^R*t^R)
      >= [(ln 2)^R*sigma^(R+1)/(10*R*2^R*(R+1)!)]*t
      -> infinity.

The bracketed constant is fixed and strictly positive: sigma>0 and R
is finite. Its size and the resulting threshold need not be practical.
The selected w tends to infinity with EVERY integer input length N.
Consequently (6) eventually exceeds one, proving the declared inequality
with coefficient R and empty remainder at every sufficiently large N.
The evaluator, verifier and sigma were fixed before R or any threshold;
no advice or family choice depends on N.

The same argument with any fixed real k>0 and integer j>k proves
C_B2(L_N)/(N*(log_2(log_2 N))^k)->infinity. This observation explains
the encoding choice; no numerical value of sigma is inferred from it.

For the declared lane endpoint, fix any real A>0. Choose an integer
w_0>=max(2,w_*) large enough for (3) and with w_0^sigma>=10A.
Set N_A=n_(w_0). For every integer N>=N_A the selected w is at least
w_0 and N>=16. Equation (5) gives C_B2(L_N)>A*N. The same evaluator
and empty-witness verifier work for every A; only N_A depends on A.

There is no direct-sum assumption, sparse-length substitution, change of
family with A, or extrapolation from finite experiments. The first score
coordinate `positive_infinity` is the typed endpoint transport value,
not a numerical gate coefficient. The second is the exact rational 1.
Both structured claims assert the same complete theorem.

## 6. Evidence and verification scope

The mathematical input is the source theorem, supported by its supplied
companion arguments. The Boolean conversion and padding proof above
make the eligible instantiation explicit. Retain the authorship and
source provenance when reusing this entry.

This is a proof-only candidate with no experiments or checker claims.
Structural validation and an empty execution plan are not mathematical
verification. The supplied evidence and every imported material premise
must be reviewed on their merits. No platform acceptance, novelty of
the underlying lower bound, or Lean verification of this particular
Boolean corollary is asserted.
