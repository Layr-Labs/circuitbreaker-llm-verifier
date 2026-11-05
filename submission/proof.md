# One fixed uniform NP family with a polylogarithmic loss

For the single total language defined below, at every sufficiently large
integer ACTUAL length N,

 C_B2(F_N) >= 4N - 32*(log2 N)^28*(log2(log2 N))^10.

In particular the same family satisfies C_B2(F_N)>=4N-N^epsilon eventually
for every fixed real epsilon>0. Its native leading projection in the selected
v1 scoring policy is (4,1). The displayed remainder is a supplied upper bound; no
claim of an optimal affine extractor or a larger leading score is made.

The exact eventual loss is R(N)=4P(N)+8k(N), where k and the ENTIRE data
prefix length P are specified by the executable integer plan. This proof
uses only the supplied full-B2 finite theorem in elimination.md Section12
and the supplied extraction arguments in disperser.md A1 and B1--B7.
The finite experiments are supporting bounded checks, not their proof.

## 1. The total language

On every bit string x of length N, evaluate family.py. At lengths with
N.bit_length()<1024 it returns the fixed legacy evaluator. Otherwise it
constructs engine.PolyPlan(N). If ready is false it again returns that same
legacy value. If ready is true it reads exactly PolyPlan.prefix bits as
construction data, pads the remaining payload with zeros to length N, and
returns engine.PolyEngine(plan,x).evaluate(). This is ONE algorithm for all
strings and lengths. It performs no search for good input data.

For all sufficiently large lengths put L=ceil(log2 N), z=ceil(log2 L),
and let P2(t) be the least power of two at least the positive integer t.
The exact executable parameters are

 k=P2(L^28*z^10), h=e=P2(L^3*ceil(z^(1/10))),
 ell=P2(L^6*ceil(z^(7/2))), d=P2(L^21*z^10),
 Y=P2(L^25*ceil(z^(19/2))), e'=P2(L^15*z^4).

Y is the independently chosen initial row width. The fractional ceilings
use exact integer roots, including sqrt(z^7) and sqrt(z^19). The program
compares ceil(log2(L^28*z^10)) with floor(log2 N) before allocating k.
Small lengths or failed guards return the fixed legacy value. All guards
pass at every sufficiently large integer length, not a sparse subsequence.

## 2. Construction data are actual input bits

The six strong linear extractor plans are listed in disperser.md B4:
initial extraction, the Q-side and X-side look-ahead maps, and E1,E2,E3.
Each consists of a deterministic polynomial condenser followed by A1's
seeded map, repeated with fresh independent seeds when needed. Its modulus
and D-by-D traversing matrix are read from the actual prefix: dimension
10000 initially and dimension z for the five internal maps. The robust-
function partition generator additionally reads its
hash-field modulus and its own A1 package, with D=100. The exact prefix
order is initial,Q-side,X-side,E1,E2,E3,hash modulus,partition A1.

Every modulus has its monic leading coefficient fixed to one and all its
other coefficients supplied as input bits. A1's matrix entries are also
input bits. Their lengths, order, and bit order depend ONLY on N. Repeated
calls share these data. Seed strings are computed by the stated algorithm;
they are not additional existential data. The payload has exactly N-P bits.

Good data exist simultaneously for all required components: finite fields
of every degree exist, and multiplication by a generator of a finite-field
extension gives A1's traversing matrix. These facts, including cyclicity,
are justified in A1. For a hardness restriction we fix ONE such tuple of
prefix bits. It works for all affine sources of the stated dimensions,
rather than being selected separately for each source.

The large condenser modulus would require order-N data if supplied as an
arbitrary polynomial. It is instead COMPUTED by B1's explicit deterministic
Artin--Schreier tower and relative-minimal-polynomial algorithm. Its binary
field degree is a power of two; the extension degree is a power of two at
least two. The supplied elementary proof establishes every field and relative
modulus. No field-finding oracle, advice, random search or witness is used.
Its construction takes polynomial time in its polynomially bounded degrees.

Bad prefix data remain valid inputs to the LANGUAGE. They still define
monic quotient-ring multiplication and fixed matrix evaluation. The evaluator
does not test their irreducibility or cyclicity and has no promise that they
are good. It returns a bit on every input. Prefix restrictions are used only
in the lower-bound proof and every restricted bit is included in P.

## 3. What the advanced evaluator computes

Initial extraction enumerates all d0-bit seeds, where its computed seed
length satisfies the explicit guard d0<=10^9 L. It produces Y-bit rows Y_i.
For good construction data and any affine input of entropy k, all but at
most 2M0/N^2 of these rows are exactly uniform, M0=2^d0.

Write each seed index as a d0-bit word, append RIGHT zero padding to make
b=ceil(d0/log2 h) base-h digits, and add one to each digit. This representation
is injective and ordered. On EVERY round the row generator refreshes a
look-ahead source using E1 and E2, runs the specified alternating extraction
up to its selected digit, and uses E3 for the next round's seed. On the last
round it still refreshes and outputs the first bit of the selected row.
This last-round refresh is part of our explicit variant. The full adaptive
conditioning proof, including additional target A-parts and the distinction
between diverged and same-prefix predecessors, is supplied in B3--B4.

For ANY h good indices their selected ell-bit rows are jointly within
100bh^2*2^-e of independent uniform rows. Taking their first bits preserves
this bound. The resulting vector therefore consists of h-wise almost
independent good bits and at most 2M0/N^2 arbitrary bad bits. The bad set
may depend on the affine SOURCE distribution, but is a fixed coordinate
set for that source, as required by the robust-function proof.

B5 supplies a strong partition generator with an explicit code/hash prefix
that makes shifted agreement between every two distinct messages at most
2^-r^.09. B6 supplies the monotone robust function applied to the first
uS row bits. Its block partitions and bias threshold are determined by N.
The threshold is computed by exact INTEGER bisection; its comparison is

 2^r*(2^(uB)-(2^u-J)^B)^M <= 2^(uBM).

There is no real-number oracle or hidden balancing advice. Bit lengths are
at most uBM+O(r)<=M0+O(r). Its two-level Bonferroni comparison and proved
Janson variant establish positive probability for both output bits under
the asserted almost independence, even for arbitrary dependence of bad
bits. No general AC0-fooling theorem, modern affine extractor, or resilient
function is used as an unexplained black box.

## 4. Guards hold at all sufficiently large integer lengths

All proofs of eventual inequalities are in B7, including their rounding
coverage. The seed budget retains the extractor dimension explicitly. For the
condensed input length n', put H=D*ceil((n')^(1/D)), R=16D+16. The
supplied bound is

 C(log n+e+log K)+C(D+1)*R*H*(e+log K+log H+log R)^2.

Only the five internal maps use D=z. They have n'=O(K), log K=O(z),
H=O(z), R=O(z) and field degree O(precision). Their general half-output
seed cost is O(z^3*e^2+log N)=O(L^6*z^(16/5))=o(ell), since16/5<7/2.
E1 needs one SHORT call: its advice is O(ell*z^2*e')=O(L^21*z^(19/2))
=o(d), and its seed cost O(z*e'+log N)=o(d). Initial extraction keeps
D=10000 fixed; its one-call advice Y*k^(1/10000)*L has scale
O(L^26.0028*z^9.501)=o(k), so d0<=10^9 L remains valid. Internal longer
field seeds are not enumerated.

The retained exact entropy checks include

 Y>=1000*b*h*d,
 k-h*Y-100*b*h^2*ell>=10d,
 e'>=100*h^3*ell+2e.

Now b=O(L/z). Row-rank loss is O(L^25*z^(91/10))=o(Y), while hY is
O(L^28*z^(48/5))=o(k). The residual loss bh^2ell=O(L^13*z^(27/10))
is smaller. Weak-seed loss h^3ell=O(L^15*z^(19/5)) is o(e').
All compared log powers have strict gaps: 91/10<19/2, 48/5<10,19/5<4.

For the robust layout r=O(L),ceil(log2 r)=Theta(z). Put

 a_r=ceil((ceil(log2 r))^(1/20)).

The guards are h>=r^3*a_r and e>=100r^3*a_r. They pass because h=e
has order L^3*z^(1/10), while the right sides have order at most
L^3*z^(1/20). The row joint error is at most2^(-r^3*a_r) eventually.
All needed observations have O(r^3) bits and coefficient counts
exp(O(r^3)); a_r tending to infinity controls their fixed constants.
The groups satisfy S>=M0/poly(log M0). From d0<=10^9 L and delta=10^-12,
q/S<=10u/N^2 and S^delta<=2^(.001L), giving q<=S^(1-delta) eventually.

Every actual prefix bit is included. The five internal data packages have
D=z; their largest field-degree scale is e'=Theta(L^15*z^4). Thus their
total size is O(L^15*z^6). Initial data are O(L), and design data smaller.
Therefore P=O(L^15*z^6)=o(k), and P+k<=N eventually.
Every other dimension and validity guard follows from the same supplied
component comparisons. Because z differs from log2 L by less than one,
all statements hold as N ranges over EVERY sufficiently large integer,
including at each rounding discontinuity.

## 5. Totality, uniform polynomial time, and NP

The scalar plan is computed by integer arithmetic and bounded halving and
binary-search loops. Its initial k exponent is checked before allocation;
its d0 cap is checked before any seed enumeration. The other integer values
and field dimensions have polynomially many bits in N. Prefix+k<=N and the
explicit source/output/seed checks bound all admitted allocations and loops.

A1's monomials are streamed with polynomial delay by an explicit iterative
weak-composition enumeration. Only the first n' coefficient monomials
are used. The initial dimension is fixed; internal dimensions are
z=O(loglog N), with polynomial delay in z and the degree bound. Tower-field
multiplication is evaluated with an explicit stack, not interpreter recursion.
All other polynomial, ring and matrix arithmetic is by bounded loops.
The deterministic condenser fields have degree polynomial in N, even when
the field's cardinality is much larger; the algorithm never enumerates that
field for an internal call.

There are at most N^(2*10^9) initial seeds from the stated cap (the extra
factor accounts harmlessly for L rounding). The index rounds and look-ahead
steps are polynomially bounded. In the robust layout r<=10^12(d0+1),
so its R=2^r loops and all partition/block loops are also bounded by a fixed
power of N. This power is enormous, but is independent of the input and
length. The integer threshold search has polynomial bit cost in M0. No
practical speed claim is made for the eventual branch at its huge threshold.

At failed lengths the fixed legacy evaluator also halts: its prime search
terminates by existence of primes, and all subsequent ring, powering and
symmetric-polynomial loops are finite. Since only finitely many lengths
fail the new guards, these exceptional computations contribute a finite
constant to a polynomial-time bound. This proves a single uniform total
polynomial-time evaluator on ALL bit strings, including invalid prefixes.
Older engine routines remain solely as named finite-audit fixtures and are
not called by the language's advanced branch.

For NP take the relation V(x,w) which rejects every nonempty witness and
otherwise accepts exactly when family.evaluate(x)=1. Its witness length is
zero and its running time is polynomial. Both directions of membership
follow immediately. The prefix is part of x, not part of w.

## 6. Full-B2 lower bound and exact accounting

Fix good prefix data at an actual sufficiently large length N. For every
k-dimensional affine flat in the N-P payload coordinates, zero padding
embeds it injectively into the N-coordinate ambient source. B1--B7 therefore
make the restricted payload function an affine disperser of dimension k.
The finite theorem in elimination.md Section12 applies to this function
on m=N-P variables:

 C_B2(payload function) >=4m-8k.

It permits every one of the 16 binary Boolean tables, arbitrary DAG sharing,
constants, aliases, unused inputs, and all terminal output cases in the
pinned model. The separate promised-pair proof and its four-object refinement
are supplied there; no gate-basis restriction is introduced by this family.
Fixing input bits in any circuit for F_N yields a circuit for that payload
function without increasing gate count. Hence

 C_B2(F_N)>=4(N-P)-8k=4N-R(N),   R(N)=4P+8k.

Every construction bit is charged in that expression. Since P=o(k),
R<=9k eventually, and k<2L^28*z^10. Since L/log2 N->1 and
z/log2(log2 N)->1, eventually R<=32*(log2 N)^28*(log2(log2 N))^10.
This proves the explicit bound over all sufficiently large integer N. For
every fixed epsilon>0 the SAME family has R<=(log2 N)^(28+epsilon)
eventually, and hence also R<=N^epsilon eventually.

The canonical score-source claims are generated from bound.json with the
same sorted evidence. They both assert this complete joint theorem. The
projection (4,1) and the full family theorem are asserted together, rather
than assembling independently incompatible coefficient and exponent claims.

## 7. Verification scope and provenance

The declared finite audit supports only its named bounded predicates.
Neither small legacy evaluations nor mocked traversal tests execute the
entire eventual construction, and no such execution is claimed. The universal
mathematical statements rest on the supplied arguments. Independent critique
is advisory; configured-provider and human judgment are separate.

Gate provenance is Mishra--Santhanam revision1, November3,2026, for the
coefficient-four architecture, followed by our supplied Section12 finite
refinement. Extraction provenance and all needed components are specified
in disperser.md. The retained older gate lemmas and finite fixtures are
credited to the declared baseline and previous submissions in the public
note; their earlier acceptance is not a mathematical premise.
