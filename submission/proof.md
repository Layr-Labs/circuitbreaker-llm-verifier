# Coefficient 5: the fixed BK/Shoup family and ternary-pair gate elimination

This submission proves the finite lower bound (5,1) for ONE fixed total
uniform family in the pinned full B2 model. The ingredients are the
Ben-Sasson--Kopparty affine disperser, reproduced in `disperser.md`; the
ternary-pair gate-elimination theorem with potential W=s+5I+3q, proved in
`elimination.md`; and the Li--Yang B2 conventions cited below. No
experiment, advice table, nonconstructive existence statement or unproved
efficiency claim enters the argument.

## 1. Fixed total family

Set f_0(empty)=0. For n>=1 let k be the least prime strictly above
2*n^(3/5), r=ceil(n/k), t=ceil(sqrt(r)). Fix once a deterministic
procedure from Section 3 realizing Shoup-preprint Theorem 3.2, with
deterministic choices fixed by that procedure, to produce a monic
irreducible q_k over F_2 of degree k. Let K=F_2[Z]/(q_k), with polynomial
basis 1,Z,...,Z^(k-1). This fixes one language, not per-length advice.

Append r*k-n zeros to x and split it into r blocks. Interpret each block
in K with its first bit the constant coefficient; this is an injective
F_2-linear map phi. Let pi take the constant coefficient, so pi is linear
and pi(1)=1. For 1<=i<=r put y_i=phi_i(x)^(2^(i+1)-1). Define

    f_n(x)=pi(Sym_r^t(y_1,...,y_r)),

where Sym_r^t sums the products indexed by t-element subsets. This uses
the precise ascending powers in BK-preprint Theorem 2.2, equation (4).
Because the powers depend only on the position i and not on n or on any
choice made later, the construction is uniform in the actual input
length.

## 2. Uniform P and NP membership

One deterministic algorithm computes this function at all lengths. Find
k by scanning and primality trial division, using q^5>32*n^3 for the
strict threshold. Bertrand's bound makes k=O(n^(3/5)), so this scan is
polynomial in the actual input length n. Shoup's unconditional
fixed-prime procedure is polynomial in k. Binary powering uses only O(i+1)
field multiplications for each exponent, whose bit length is i+1. Field
operations in the polynomial basis have polynomial bit cost.

Compute the symmetric polynomial without characteristic-two division:
start E_0=1 and E_j=0 for 1<=j<=t; after reading y_i update
E_j=E_j+y_i*E_(j-1) in descending j from min(i,t) to 1. Splitting subsets
by inclusion of i proves the recurrence; it takes O(r*t) field
operations. Output pi(E_t). Thus one polynomial-time evaluator A computes
the total family. Use V(x,w)=1 iff w is empty and A(x)=1, with witness
bound zero. Both directions are immediate, establishing L in P, hence NP.
The witness is not an input to the target circuits.

## 3. Disperser property and the dimension d(n)

`disperser.md` reproduces the binary case of Ben-Sasson--Kopparty
Theorem 7.3 and proves it from the field and polynomial preliminaries
rather than invoking it as a black box. Its Lemma states: for prime k,
1<=t<r with k*t>=r(r+1), and any nonzero linear pi:K->F_2, the function
of (a_1,...,a_r) defined with the weights w_i=r+2-i is nonconstant on
every affine subspace of K^r of dimension strictly greater than
kr/t + kt + r(r+1). The construction of Section 1 uses exactly these
weights, because 2^(i+1)-1 = 2^(r+2-i)-1 under the substitution i -> r+1-i
and the symmetric sum is over all t-element subsets.

With the parameters of Section 1 one has, for every sufficiently large n,
k<=3*n^(3/5), r<=n^(2/5)+1 and t<=n^(1/5)+1, while k>2*n^(3/5) gives
r<n^(2/5)+1. Substituting these native bounds into the three terms gives
kr/t+kt+r(r+1)+1 <= 6*n^(4/5) for every sufficiently large integer n.
Set d(n)=floor(6*n^(4/5))+1. Eventually 1<=d(n)<=n and
d(n)/n<=6*n^(-1/5)+1/n tends to zero. BK-preprint Theorem 2.2 applies to
arbitrary integer input length through the injective map, not to a
divisibility assumption or a sparse subsequence. Zero padding is inside
the theorem's injection; the circuit bound remains measured against the
n actual input bits.

The same argument at n=0 is the declared base case f_0(empty)=0. Thus the
hypothesis needed in Section 4 holds for every sufficiently large integer
length.

## 4. Ternary-pair gate elimination

`elimination.md` proves the following theorem, self-contained in the
pinned model. Let d be a nonnegative integer and suppose f:F_2^n->F_2 is
nonconstant on every affine subspace of dimension strictly greater than
d. Then every circuit for f in the full binary basis has at least

    5n - 10(d+1)

gates, including arbitrary sharing, repeated pins, constants, unary
tables, degenerate gates and a single terminal output.

The proof runs a fair affine system A v=N u+c, marks disjoint three-point
pairs with one excluded corner each, and charges the potential
W=s+5I+3q, where s counts equations plus acyclic gates, I counts used
unpaired parameters and q counts marked pairs. Its legal operations are
restrictions by injective affine maps, parameter/auxiliary exchanges and
literal rerouting; none creates a gate. The deletion rows D1-D4 and the
affine move M are nonincreasing in W. The paired rows P0, P2, PC, PE and
the AND-elimination rows AP1, AP2, GP satisfy Delta W<=10 Delta a, with
Delta a in {0,-1}, and the tabulated bounds are

    AP1: Delta W<=-11,  AP2: Delta W<=-10,  GP: Delta W<=-9-k<=-10.

The two strengthened rows shorten the argument as follows. In AP2 the
constant gate G has deletable children; a fourth distinct deletion is
forced either by a child of G outside the two partner gates, or by the
constant descent through the partner gate that reads G and its fixed
partner. In GP the unique child of G in the k=1 case is not the output,
so it has a child; the constant-or-literal split of that child either
deletes the child directly or transfers the identification to an
ancestor, and the finite acyclic chain ends with one further distinct
deletion. The rows therefore delete 1, 2, 3 and 2+k distinct objects
respectively, and the batch counting rule makes repeated pins, common
parents, shared descendants and cyclic affine dependencies harmless.

Summing with W_0<=|C|+5n and W_T>=0 gives

    |C|+5n >= W_T+10(n-a_T) >= 10(n-d-1),

and a_T<=d+1 because coverage holds at every loop state and the integer
m+s strictly decreases. Rearranging gives the theorem. Each affine row
was introduced one-for-one from a charged gate, exchanges changed roles
without extra computation, and restrictions never added a gate, so the
inequality bounds the original gate count and not a potential mistaken
for it.

## 5. The bound at actual input length

Apply the theorem to the family of Section 1 with the d(n) of Section 3.
For every sufficiently large integer n,

    C_B2(f_n) >= 5n - 10(d(n)+1)
              = 5n - 10*floor(6*n^(4/5)) - 20
              = 5n - o(n).

Hence liminf over all integer n of C_B2(f_n)/n >= 5, and for every
epsilon>0 there is N_epsilon with C_B2(f_n)>=(5-epsilon)*n for every
integer n>=N_epsilon. The declared guarantee is the polynomial/log
expression c*n^p*(log2 n)^q*(log2 log2 n)^r with c=5, p=1, q=0, r=0,
and the declared projection is the same finite pair (5,1); the full
expression and the projection are both established by the displayed
inequality.

## 6. Model conventions and accounting

The full B2 model includes all 16 binary truth tables, unbounded sharing
and depth, one output and all input assignments; inputs and wires are not
gates. The source allows constant nodes during substitutions and handles
repeated pins and unary effects by normalization. Output cases are
covered: the elimination proof excludes constant, unary and terminal
literals on the restricted domains before charging them, and a retained
negated projection costs one, as required. Negations on edges are
absorbed into the child's truth table in `elimination.md`, and the affine
system is uniquely solvable everywhere. The length n=|x| is the actual
input bit length; blocks, field encodings, padding and every linear loss
are accounted for in d(n). Witness bits are not circuit inputs and are
never charged.

The same family serves the coefficient and the exponent jointly: there is
no exponent parameter in the construction, and the bound is linear at
this one declared exponent. The P evaluator of Section 2 is fixed for all
lengths and does not depend on the circuits or on any review outcome.

## 7. Sources

- LY-full: ECCC TR21-023, front-page date 2021-02-20; Sections 1.1, 2.1,
  2.4 (footnotes 3-4) and 3.3 for the full B2 model, Definition 3.6 and
  Theorem 4.1 for the model's counting conventions;
  https://eccc.weizmann.ac.il/report/2021/023/download/;
  SHA256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062
- BK-preprint: Author preprint 2010-03-12; Theorem 2.2, equation (4),
  printed pp.6-7 (PDF pp.7-8);
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf;
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf
  (binary case of Theorem 7.3, reproduced and proved in `disperser.md`)
- Shoup-preprint: Author preprint 1989-01-31; Section 3, Theorem 3.2,
  printed p.3 and the fixed-prime consequence;
  https://www.shoup.net/papers/detirred.pdf;
  SHA256 00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644
- Mishra--Santhanam: ECCC TR26-260, 2026-10-19, Sections 2-4 for the fair
  affine system, role exchanges and deletion machinery adapted in
  `elimination.md`. The three-point-pair domain, the weights W=s+5I+3q
  and the affine-disperser variant are proved in `elimination.md` itself
  rather than invoked from that paper.

The two supporting documents are part of the supplied evidence. Their
statements, constants and row bounds are the ones used above; no external
premise is left to a link or to a retrieval step.
