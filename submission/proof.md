# Superlinear full-B2 circuits for padded element distinctness

## Main theorem and attribution

One fixed total language L in P satisfies, for every sufficiently large integer N,

    C_B2(L_N) >= N log_2(log_2 N).

The coefficient is one and the remainder is empty. Its exponent-one projection is

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
to the cited authors. This entry specifies the total family and proves its
uniform P/NP eligibility, Boolean model transfer and all-length endpoint
quantifiers. No prior CircuitBreaker result or accepted judgment is a premise.
The supporting source material is supplied in source-*.md with original
file names, theorem labels, source provenance and attribution. Provenance
hashes identify versions. The source theorem supplies a fixed positive
exponent on log N. The declared iterated-log bound works for every such exponent.

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

There is no hypothesis omega_k<9/4 in this theorem. That assumption
belongs to the separate compression theorem (Theorem 2.1), which this
entry does not use. Section 5 of the source, `list.tex`, proves (1)
by an exhaustive ports/no-ports distinction, using the geometric and
counting results identified in the supplied evidence. Neither specific
case for characteristic two is assumed.

## 3. Constant-overhead conversion from arbitrary Boolean circuits

**Lemma.** A full-B2 circuit with s gates and q input bits has an exact
trivalent F2 network for its acceptance predicate of volume at most
8s+3q+5.

**Proof.** For q>=1 and nonconstant output, discard gates outside the
output's ancestor subgraph. If a gate has a constant input pin, substitute
the constant into its binary truth table and reconnect that now-ignored
pin to x_1. This keeps a legal full-B2 gate at the same cost and preserves
its Boolean function. It handles constant nodes without charging an
uncharged extra computational basis.

For each retained binary gate g use the ternary table

    T_g(a,b,c) = 1 if c=g(a,b), and 0 otherwise.

For each logical signal with d>=2 consumers, use a binary distribution
tree with d-1 equality vertices E(a,b,c)=1 iff a=b=c. For d=1 use a
direct wire. Count two input-pin occurrences separately even when they
have the same predecessor. Add one unary output pin P(c)=c. Across all
signals there are at most 2s+1 consumers, so at most 2s+1 copy vertices.
Each primary input occurs as one exterior bit. Connect each unused input
to a unary all-one tensor, requiring at most q additional vertices.

There are at most v=3s+q+2 vertices, every one of arity at most three,
and exactly q exterior wires. If i is the number of internal wires,
2i+q<=3v. Therefore

    volume = v+i+q <= (5v+q)/2
           <= (15s+6q+10)/2 <= 8s+3q+5.

For every fixed input, acyclicity uniquely determines all original
signals, and each equality tree uniquely determines their copies. The
gate and equality constraints therefore have exactly one consistent
internal assignment before imposing the output pin. The pin keeps it
iff the circuit accepts. The contraction is exactly its zero-one
predicate over F2, with no cancellation or probabilistic argument.
Undirected cycles arising from sharing do not change this uniqueness.

If the output is constant, use its scalar tensor and one unary all-one
table for each input, of total volume 2q+1. This includes q=0. If the
output is an unnegated primary input, use its unary output pin and the
other all-one tensors. A retained negated-projection gate is treated as
one ordinary B2 gate. Thus constants, unary and terminal outputs are
all included. QED.

This construction covers arbitrary nonuniform wiring, unbounded fan-out,
repeated pins, shared descendants, all sixteen gates and unrestricted
depth. It is not a formula, monotone or formal-arithmetic simulation.
Only the direction Boolean circuit -> tensor network is used.

## 4. The distinguished input lengths

Let s_w=C_B2(D_w). By the lemma and (1), eventually

    8s_w+3n_w+5 > n_w*w^sigma,
    s_w > (n_w*(w^sigma-3)-5)/8.                      (2)

Since sigma>0, s_w/n_w tends to infinity. In fact
s_w=Omega(n_w*(log_2 n_w)^sigma), since log_2 n_w=w+log_2(2w)=Theta(w).
The next step transfers this rate to every sufficiently large input length.

## 5. Every integer input length and every linear coefficient

For N>=4 take the w from the evaluator. Consecutive lengths satisfy

    n_(w+1)/n_w = 2*(w+1)/w <= 4.

Maximality gives n_w<=N<n_(w+1)<=4n_w. Fix the ignored suffix of any
circuit for L_N to zero. The circuit now computes D_w without increasing
its gate count; substituting constants and retaining or simplifying
truth tables is allowed in the pinned B2 model. Hence C_B2(L_N)>=s_w.

For w^sigma>3, (2) and n_w>N/4 give

    C_B2(L_N) > N*(w^sigma-3)/32 - 5/8.              (3)

Since N<4n_w=8w*2^w, we have log_2 N < w+log_2 w+3 <= 3w
for all sufficiently large w. Also w tends to infinity with N. Thus
log_2(log_2 N) <= log_2(3w). Divide (3) by N log_2(log_2 N).
Its right-hand side after division tends to infinity, since
w^sigma/log_2(3w) tends to infinity for each fixed sigma>0. In
particular the ratio exceeds one at every sufficiently large integer N.
This proves the full declared inequality, including its coefficient and
empty remainder. The same calculation gives divergence after dividing
by N*(log_2(log_2 N))^k for each fixed real k>0.

Fix any real A>0. Choose an integer w_0>=max(1,w_*) with
w_0^sigma>=32A+8 and put N_A=n_(w_0). For every integer N>=N_A,
its w is at least w_0, so (3) implies

    C_B2(L_N) > A*N + 5N/32 - 5/8 >= A*N,

using N>=4. The same language, evaluator and empty-witness verifier work
for every A. Only the threshold varies. This proves exactly the declared
lane endpoint at exponent 1, with all integer lengths covered.

There is no direct-sum assumption, sparse-length substitution, limsup
argument, change of family with A, or numerical extrapolation. The first
score coordinate `positive_infinity` is the required typed endpoint
transport value, not a numerical gate coefficient. The second is the
exact rational 1. Both claims assert this entire theorem.

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

