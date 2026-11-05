# A fixed uniform NP family with subpolynomial gate remainder

For the one total language defined here, at every sufficiently large integer
actual length N, we prove C_B2(F_N)>=(7/2)N-2^(200*(log2 N)^(3/4)*log2(log2 N)). More strongly, for
every fixed epsilon>0 the remainder can be N^epsilon for the SAME family.
The pinned leading-term projection remains (7/2,1).

The finite gate theorem is supplied in elimination.md Section 10. The
large-length disperser adapts Shaltiel's 2011 challenge-response construction,
with supplied proofs of its specialized extraction ingredients in
disperser.md A1--A6. A1 adapts Shaltiel--Umans's low-degree construction;
A2 adapts Ben-Sasson--Kopparty's cubic coefficient argument; A3 supplies
Haramaty--Shpilka's bounded-rank derivative argument; A4 adapts Rao's affine
somewhere-random condenser; A5--A6 adapt Shaltiel's block extractor and
challenge-response analysis. Our substitutions use smaller outputs and
larger thresholds, and specify every construction field and matrix as
charged input data. No cited specialized component is assumed without
its supplied argument.

Primary provenance: [Shaltiel 2011](https://www.cs.haifa.ac.il/~ronen/online_papers/AffineDisp.pdf),
[Shaltiel--Umans](https://www.cs.haifa.ac.il/~ronen/online_papers/SUlong-rev2.ps),
[Haramaty--Shpilka revision 1](https://eccc.weizmann.ac.il/report/2009/080/revision/1/download),
[Rao](https://homes.cs.washington.edu/~anuprao/pubs/lowweight.pdf), and
Ben-Sasson--Kopparty, Affine Dispersers from Subspace Polynomials, 2010 author
version, Sections 3 and 4. These credits are mathematical provenance;
Yukon submission dependencies are separately declared in the public note.

## Exact fixed evaluator

Keep the existing family.py evaluator as legacy_evaluate. On a binary
string x of actual length N, return legacy_evaluate(x) if N.bit_length()
is less than 1024. Otherwise compute Plan(N) in engine.py. If its ready
flag is false, return legacy_evaluate(x). If ready, split off exactly
Plan.prefix bits as construction data and run Engine(plan,x).evaluate().
All integer operations, ordering choices, seed interpretation and padding
are specified by engine.py. This is one deterministic evaluator for all
strings and lengths; it chooses no data by an unspecified search.

Plan uses L=ceil(log2 N), z=ceil(log2 L), k=2^(100*ceil(L^(3/4))*z), t=2^ceil(L^(1/4)),
ell=2^ceil(sqrt L), pmax=2^ceil(L^(3/4)), q=ceil(sqrt L)*z. All these roots
are computed by integer comparisons: e.g. ceil(L^(3/4)) is the least integer
u with u^4>=L^3, not a floating-point operation. The padded ambient
length n=t^ceil(L/log2 t) satisfies N<=n<Nt. The tree is truncated
when its node length is less than floor(k^.4). The exact active index
limit is hA=floor(h/q). The remaining quantities are the exact integer
parameters in the component proofs, specialized as follows.

Let K=2^floor(log2 k/2), r=2^floor(log2 K/2),
D=ceil(100 ceil(log2 n)/log2 K), h0=D*ceil(n^(1/D)). E0 is the seeded
extractor on n bits, output r, error 1/8, field degree selected by the
explicit q-bound. E1 uses the same output but error 2^-floor(K^.001).
Let T=2^floor(log2 r/10), b=1+ceil(d0/log2 T), where d0 is E0's seed
length. Every use of BE has b blocks, each padded to n bits. Its initial
2^d0 rows are grouped into T-row chunks at each reduction. Final
chunks are padded by repeating their last row. The output is the first
pmax bits of the final r-bit row.

For the somewhere-random chunk extractor, the exact row/width sequence
is t_j=T/2^j and r_j=r/2^j for j<log2 T. Let w_j=floor(r_j^.1),
p_j=2w_j+1 (odd). Its base extractor uses p_j field coordinates
and computes the trace coordinates of the seventh power, on two slices
padded with zeros. The near-half extractor uses ambient t_j*r_j,
D_j=100000, m_j=2^floor(999 log2 r_j/1000), error
2^-floor(r_j^.00001), and the accompanying explicit seeded parameters.
It concatenates ceil((r_j/2)/m_j) fresh seeds and truncates the output
to exactly r_j/2 bits. The base extractor supplies exactly that many
seed bits. (At large lengths m_j divides r_j/2; the finite algorithm
also defines the other case by truncation.)

## Every construction datum is an actual input bit

Prefix order is E0's modulus and matrix, E1's modulus and matrix, then
for each j in increasing order the base field modulus and the near-half
extractor's modulus and matrix. A degree-s modulus is specified by s
bits for its lower coefficients, with leading coefficient 1 implicit.
A D-by-D matrix uses D^2 consecutive s-bit field elements in row order.
Least-significant bits come first. There are no additional hardwired
length-dependent field or matrix choices. Repeated uses at the same
scale share these input data. Every matrix entry, every modulus bit,
and the entire prefix count towards N.

Good choices exist: finite fields of every size 2^s and their extensions
exist by the Frobenius root argument, and multiplicative groups of finite
fields are cyclic. Multiplication by a generator of F_(2^(sD)) in a basis
over F_(2^s) is a D-by-D matrix cycling through all nonzero vectors.
These existence statements supply prefixes on which the component proofs
apply. The evaluator never searches for these moduli or matrices. It
does not validate them, either. For every other prefix it runs the same
bounded quotient-ring and positive matrix-power algorithms. A monic
modulus always defines reduction and multiplication; even a reducible
modulus therefore gives a total Boolean output. The trace loop's
constant coefficient is used, agreeing with the trace bit on good data.

## Uniformity and polynomial time, including bad prefixes

All prefix dimensions and loop bounds depend only on N. Plan has no
field-data search. It returns ready=false immediately when k>N, before
constructing any seeded-extractor parameters. Computing this comparison
is polynomial: 100*ceil(L^.75)*ceil(log2 L)<=10000*L for L>=2,
so even this intermediate k has only O(log N) bits. On the remaining
branch k<=N, all later scalar dimensions are polynomial in N. Each p_j=2w_j+1 is computed directly; no prime-distribution theorem
or search is used in the large-length branch. Integer roots are found by
binary search using fixed or O(log N) integer exponents. D<=O(log N),
while each near-half D_j is the fixed constant 100000. All scalar integer
bit lengths in parameter planning are polynomial in N, including the
powers of two used to specify seed spaces. The field degrees and data
lengths are polynomial in N before any prefix allocation.

Once ready is true, guards give d0<=1000000 ceil(log2 n), so enumeration
of E0 seeds uses at most n^2000000 possibilities. The guard
h*t^b<=ell bounds tuple enumeration per tree node, and the complete
tree has at most 2n nodes. The stored matrix sizes, all row lengths,
all field degrees, all output counts, all polynomial exponents and
all concatenation counts are polynomial in N. The low-degree extension
stores or streams only its first n monomials; enumerating these in total
degree and lexicographic order has polynomial delay. For each monomial,
compute powers by binary powering and multiply at most D field values.
There is no enumeration of the field or a reconstruction search in the
evaluator: reconstruction is solely a counting proof of extraction.

Field/matrix arithmetic is polynomial even on invalid prefixes. The
number of tree-node visits is the fixed finite tree size. Both tree
traversal and weak-composition enumeration are iterative; they do not
require an unbounded interpreter recursion limit. Thus one fixed
polynomial in N bounds the whole evaluator, with an enormous constant
exponent allowed. Failed guards invoke the already polynomial legacy
evaluator. The witness verifier accepts exactly the empty witness and
runs this evaluator. This gives a total uniform P language, hence a
uniform NP language with witness bound zero, without advice.

## All large actual integer lengths take the new branch

Here L=ceil(log2 N) tends to infinity through integers and log n=L+O(L^.25).
Thus the asymptotic comparisons in the component proofs hold along ALL
large integer N, including at rounding discontinuities. They imply:

* E0 and E1 satisfy the explicit seeded entropy inequalities. Their
  monomial degree is K^(.01+o(1)); E0's field degree is O(log K+log D),
  and E1's is O(K^.001+log K). Thus their advice terms are respectively
  K^(.51+o(1)) and K^(.511+o(1)), both o(K).
* d0/log n is bounded by an absolute constant less than 1000000:
  from s<=8*(21+log2(16D+16)+log2 h0+4*(6+log2 r)),
  D=ceil(100 ceil(log2 n)/log2 K), log2 h0<=.01 log2 K+log2(2D),
  log2 r<=.5 log2 K, one gets d0<10000 log2 n eventually.
* T>=2, b<=t, hA>=1, ell*pmax<=r,
  h*t^b<=ell, and ell^(6hA+6)<=pmax. The exponents are respectively
  O(L^.5/log2 L) for log2(t^b), L^.5 for log2 ell,
  O(L^.75/log2 L) for the largest log2 confidence, and L^.75+O(1)
  for log2 pmax. The only nonvanishing entropy-loss ratio in the
  structural descent is h*log2(4b)/log2 k, which tends to 1/400
  and leaves the stated k^.77 and k^.75 margins; the other relevant
  ratio log2(t^q)/log2 k also tends to 1/100.
* Each near-half extractor satisfies its explicit advice/entropy bound
  with threshold r_j/4. Its concatenated seed length is r_j^(.00101+o(1)),
  whereas the base extractor supplies p_j^(1/8)=r_j^(.0125+o(1)) bits.
  All scale lengths r_j are at least r/T>=r^.9 and tend to infinity,
  so the inequalities hold uniformly over every one of O(log r) rounds.
* The chunk extractor's final width r/T>=r^.9 exceeds E1's seed
  length K^(.001+o(1)).

The number P of prefix bits is

 (D^2+1)(s0+s1)+sum_j [p_j+(100000^2+1)*s_j].

Here p_j=O(r_j^.1), s_j=O(r_j^.00001+log r_j), and there are O(log r)
terms. Therefore P=K^(.05+o(1))=k^(.025+o(1))=N^o(1). In particular
P+k<=N eventually. Every ready condition is consequently met at all
sufficiently large lengths. The fixed cutoff N.bit_length()<1024
excludes only finitely many lengths.

## Circuit lower bound with actual prefix cost

Fix any sufficiently large N and good values of all P prefix bits.
Restricting a full B2 circuit to these constants cannot increase its
number of gates in the pinned normalization. The remaining m=N-P bits
are embedded injectively into n coordinates by zero padding. The supplied
challenge-response argument makes the new disperser nonconstant on
every affine subspace of dimension k in those n coordinates. Therefore
its restriction to the payload is also nonconstant on every dimension-k
affine subspace of F2^m, since injective linear embedding preserves
dimensions. The universal finite theorem of elimination.md Section 10
then gives, for every circuit computing the WHOLE new family,

 C_B2(F_N)>= (7/2)(N-P)-(13/2)k
          = (7/2)N-R(N),  R(N)=(7/2)P+(13/2)k=N^o(1).

This conclusion charges the entire prefix and the exact integer N.
No claim that the entire input function is itself an affine disperser
is needed. For every real epsilon>0, R(N)<=N^epsilon eventually, hence
C_B2(F_N)>=(7/2)N-N^epsilon at all sufficiently large integer N. In
particular the stronger explicit fixed remainder is 2^(200*(log2 N)^(3/4)*log2(log2 N)), and
the leading pair is still (7/2,1). The stronger remainder changes neither
the exponent 1 nor the leading coefficient 7/2 in the pinned v1 policy.

The last explicit bound follows from
log2 k=100*ceil(L^.75)*ceil(log2 L) and P=k^(.025+o(1)):
log2 R(N)=(100+o(1))*L^.75*log2 L. Replacing L=ceil(log2 N)
by log2 N changes this asymptotic ratio by 1+o(1). Therefore
R(N)<=2^(200*(log2 N)^(3/4)*log2(log2 N)) eventually.
All limits and margins hold along every sufficiently large integer N.
