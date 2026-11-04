# Coefficient 3 for the fixed BK/Shoup affine disperser family

Claim: for the single total family f_n defined below (the same family as the
published 31/10 seed), with d(n)=floor(6*n^(4/5)) the Ben-Sasson–Kopparty
dimension threshold, C_B2(f_n) >= 3n - 4*d(n) - 4 for every sufficiently large
integer n; in particular C_B2(f_n) >= 3n - o(n), liminf C_B2(f_n)/n >= 3, and
the declared pair is (3,1).

Structure. Part A defines the family and proves NP membership (unchanged
from the seed). Part B states the one external premise, the Ben-Sasson–Kopparty
affine-disperser theorem. Part C is a COMPLETE, self-contained gate-elimination
proof that every B2 circuit for a function that is nonconstant on all affine
subspaces of dimension > d has at least 3n - 4d - 4 gates. Part C is our own
write-up of the idea of Demenkov–Kulikov (MFCS 2011, LNCS 6907, "An elementary
proof of a 3n-o(n) lower bound on the circuit complexity of affine dispersers");
it is reproduced in full and does not rely on that paper. No new circuit
lower bound beyond the literature is claimed; the 31/10 bound of Li–Yang is
stronger, but this proof replaces their long case analysis by a short one.
Part D combines A–C and handles the quantifiers.

## Part A. Fixed total family and NP membership

Set f_0(empty)=0. For n>=1, let k be the least prime strictly above 2*n^(3/5),
r=ceil(n/k), t=ceil(sqrt(r)). Fix once a deterministic procedure from Section 3
realizing Shoup-preprint Theorem 3.2, with deterministic choices fixed by that
procedure, to produce a monic irreducible q_k over F_2 of degree k. Let
K=F_2[Z]/(q_k), with polynomial basis 1,Z,...,Z^(k-1). This fixes a single
language, not per-length advice.

Append r*k-n zeros to x and split it into r blocks. Interpret each block in K
with its first bit the constant coefficient; this is an injective F_2-linear
map phi. Let pi take the constant coefficient, so pi is linear and pi(1)=1.
For 1<=i<=r put y_i=phi_i(x)^(2^(i+1)-1). Define
f_n(x)=pi(Sym_r^t(y_1,...,y_r)), where Sym_r^t sums the products indexed by
t-element subsets. This uses the precise ascending powers in BK-preprint
Theorem 2.2, equation (4).

Uniform polynomial time. One deterministic algorithm computes this function at
all lengths. Find k by scanning and primality trial division, using q^5>32*n^3
for the strict threshold q>2*n^(3/5). Bertrand's postulate gives k=O(n^(3/5)),
so the scan is polynomial in the actual input length n. Shoup's unconditional
fixed-prime procedure is polynomial in k. Binary powering uses O(i+1) field
multiplications for each exponent 2^(i+1)-1. Field operations in the polynomial
basis have polynomial bit cost. Compute the symmetric polynomial without
division: start E_0=1 and E_j=0 for 1<=j<=t; after reading y_i update
E_j=E_j+y_i*E_(j-1) in descending j from min(i,t) to 1. Splitting t-subsets by
whether they contain i proves the recurrence; it takes O(r*t) field operations.
Output pi(E_t). Thus one polynomial-time algorithm A computes the total family.

NP membership. Let L={x : f_|x|(x)=1}. Use the verifier V(x,w)=1 iff w is the
empty string and A(x)=1, with witness bound W(n)=0. If x in L then V(x,empty)=1;
if V(x,w)=1 for some w with |w|<=0 then w is empty and A(x)=1, so x in L. Both
directions hold, L is in P and hence in NP. Witness bits are never circuit
inputs; the circuits below have exactly the n input bits of x.

## Part B. The affine-disperser premise (the only external premise)

Definition. g:{0,1}^n->{0,1} is nonconstant above dimension d if g is
nonconstant on every affine subspace of F_2^n of dimension at least d+1.

Premise (BK-preprint Theorem 2.2, equation (4), with p=2). With the above
k, r, t, phi, pi and powers, f_n is nonconstant on every affine subspace of
F_2^n of dimension strictly greater than 6*n^(4/5).

Hypotheses, checked one by one. (i) The base field is F_2 (p=2) and K is the
degree-k extension given by the irreducible q_k, with k prime and k>2*n^(3/5).
(ii) The number of blocks is r=ceil(n/k) and the threshold is t=ceil(sqrt(r)).
(iii) Block i is raised to the power 2^(i+1)-1 for i=1..r (ascending, as in
equation (4)), the symmetric function Sym_r^t is applied and the result is
projected by the nonzero F_2-linear functional pi. (iv) The dimension bound is
stated in terms of the actual input length n. The theorem is applied to arbitrary input length n via the
injective F_2-linear zero-padding map phi: an affine subspace of F_2^n of
dimension D is mapped by phi injectively onto an affine subspace of K^r of the
same dimension D, so no divisibility of n by k and no sparse subsequence is
needed. If one reads the BK threshold instead as a function of the padded
length r*k <= n+k = n+O(n^(3/5)), it is still at most 6*(n+k)^(4/5) = o(n),
and Part D below only uses d(n)=o(n); the coefficient 3 is unaffected.

Consequence. Put d(n)=floor(6*n^(4/5)). An integer dimension D satisfies
D>6*n^(4/5) iff D>=d(n)+1. So f_n is nonconstant above dimension d(n).

## Part C. Self-contained gate elimination: at least 3n-4d-4 gates

Theorem C. Let n>=d+2 and let g:F_2^n->F_2 be nonconstant above dimension d.
Every B2 circuit computing g has at least 3n-4d-4 gates.

Model. B2 as in Li–Yang (ECCC TR21-023, Sections 1.1, 2.1, 2.4, 3.3): an
acyclic graph with n input nodes x_1..x_n, free constant nodes 0 and 1, and
gates; each gate has two ordered pins (repeated pins allowed) and any of the 16
binary truth tables; unbounded fan-out and depth; one output node; exact
computation on all inputs. Every gate is counted, including degenerate,
unary, negated-projection or unused gates.

### C.1 Generalized circuits over an affine subspace

Fix an affine subspace S of F_2^n of dimension m. An affine function on S is
the restriction to S of a map x -> <a,x>+b. A nonconstant affine function l on
S takes value a on an affine subspace {s in S : l(s)=a} of dimension exactly
m-1, for each a in F_2.

An S-circuit is an acyclic graph with: input nodes, each labelled by an affine
function on S; constant nodes 0,1; gates with two ordered pins and a truth
table; one output node; and an output bit beta in F_2. On s in S each input
node takes the value of its label at s, nodes are evaluated in topological
order, and the circuit computes (value of the output node) XOR beta. Let G be
the number of gates and I the number of input nodes; the measure is mu=G+I.
Constant nodes cost zero. Clearly mu>=0. A successor of a node w is a gate
having w as a pin; out(w) is the number of distinct successors.

Start. Given a B2 circuit with G_0 gates computing g, take S=F_2^n, label x_i
by the coordinate function, beta=0. This S-circuit computes g and has
mu = G_0 + n.

Restriction. If an S-circuit computes h on S and S' is an affine subspace of
S, restricting every label to S' gives an S'-circuit computing h restricted to
S', with the same graph and the same mu. (Every node takes the same value at
each s in S'.)

Why this adds no uncharged computational basis. The generalized model is used
only to prove a LOWER bound on mu, and only via mu_start = G_0 + n exactly:
the start circuit has beta=0 and literally the gates of the given B2 circuit.
Every input node, including any affine label created later, is charged 1 in mu;
N4 exchanges one counted gate for one counted input, so a free XOR never lowers
mu. The bit beta only records negations already present in counted gates and
never computes anything. No potential is identified with a gate count: the
conclusion is G_0 >= mu_start - n and mu_start >= 4(n-d-1).

### C.2 Normalization rules

Rules may be applied in any order chosen by the prover; each application is
justified on its own (repeated pins of B2 are handled by N3).

Each rule preserves the function computed on S and does not increase mu.

- N0 (constant label). An input node u whose label is constant c on S:
  redirect every pin at u to constant node c, if u is the output make constant
  node c the output, delete u. I decreases by 1.
- N1 (unused input). An input node with no successors that is not the output:
  delete it. I decreases by 1.
- N2 (merge). Input nodes u != v with label(v)=label(u)+e, e in F_2: redirect
  every pin at v to u, replacing the table of each such gate by its composition
  with x -> x+e on that pin; if v is the output, make u the output and add e to
  beta. Delete v. I decreases by 1.
- N3 (degenerate gate). A gate A is degenerate if a pin is a constant node, or
  its two pins are the same node, or its truth table does not depend on both
  arguments. Then A's value equals a constant c, or equals h(w) for one
  non-constant node w that is a pin of A and h in {identity, negation}. In the
  first case redirect every pin at A to constant node c (and if A is the
  output, make c the output). In the second case redirect every pin at A to w,
  composing each successor's table with h on that pin; if A is the output make
  w the output and add [h=negation] to beta. Delete A. G decreases by 1.
  (This is the Li–Yang normalization: negations are absorbed into successors'
  tables or the free output bit; no uncharged computational basis is added.)
- N4 (top XOR). A gate A whose pins are two distinct input nodes u,v, whose
  table is (a,b) -> a+b+c, and with label(u)+label(v) nonconstant on S: A
  computes the affine function l=label(u)+label(v)+c. Create a new input node
  labelled l, redirect every pin at A to it (and the output, if A is the
  output), delete A. G decreases by 1 and I increases by 1: mu unchanged.
- N5 (redundant intermediate). A gate B with pins {u, A}, where u is an input
  node, A is a gate whose pins are {u, v} with v an input node, B is the only
  successor of A, and A is not the output. B's value is a function
  H(value(u), value(v)) obtained by substituting A's table into B's. Make B's
  pins (u,v) with table H and delete A. G decreases by 1.
- N6 (dead gate). A gate with no successors that is not the output: delete it.
  G decreases by 1.

Termination. N0, N1, N2 keep G and decrease I; N3, N4, N5, N6 decrease G.
So every rule strictly decreases the pair (G, I) lexicographically, and any
sequence of rule applications is finite. A circuit to which no rule applies is
normalized. Since no rule increases mu, normalizing never increases mu.

Facts about a normalized S-circuit (immediate from the rules):
(F1) every label is nonconstant (N0), and no two labels are equal or
complementary (N2), so label(u)+label(v) is nonconstant for distinct inputs;
(F2) every gate has two distinct non-constant pins and a table depending on
both arguments (N3); such a table is either XOR-type (a,b)->a+b+c or AND-type
(a,b)->((a+alpha)(b+gamma))+c for some alpha,gamma,c in F_2 (the 16 binary
functions are 2 constants, 4 that depend on one argument, 2 XOR-type and 8
AND-type);
(F3) every gate that is not the output has at least one successor (N6).

### C.3 The elimination step

Lemma. Let g be nonconstant above dimension d. Let S be an affine subspace of
dimension m>=d+2 and let C be a normalized S-circuit computing g on S. Then
there is an affine subspace S' of S of dimension m-1 and an S'-circuit
computing g on S' with mu at most mu(C)-4.

Proof. Since m>=d+1, g is nonconstant on S, so the output is not a constant
node. If the output were an input node u, g would equal label(u)+beta on S,
hence be constant on {s in S: label(u)(s)=0}, which has dimension m-1>=d+1, a
contradiction. So the output is a gate and C has gates. Let A be a gate that
is first in a topological order. Its pins are not gates and (F2) not
constants, so they are two distinct input nodes u,v. By N4 and (F1), A is not
XOR-type, so by (F2) A is AND-type: A=((u+alpha)(v+gamma))+c. If A were the
output, g would be constant on {s in S: label(u)(s)=alpha} (dimension m-1>=d+1),
a contradiction. So by (F3) A has at least one successor.

Case 1: out(u)>=2 or out(v)>=2; by symmetry out(u)>=2. Let B be a successor
of u with B != A. A has a successor C with C != B: otherwise B is the only
successor of A, B has the two distinct pins u and A, A has pins {u,v} and A is
not the output, so N5 applies, contradicting normalization. Let
S'={s in S: label(u)(s)=alpha}, of dimension m-1 by (F1). Restrict to S'.
Apply N0 to u (I drops by 1); now A has the constant pin alpha, so A has the
constant value c and N3 deletes A, redirecting its pins to constant c; B has
a constant pin (from u) and N3 deletes B; C has a constant pin (from A) and N3
deletes C. A,B,C are distinct (B != A by choice; C != A by acyclicity; C != B
by choice), and deleting one of them does not remove the constant pin of
another, so G drops by at least 3. Total drop of mu: at least 4.

Case 2: out(u)=out(v)=1, so A is the only successor of both. Let C be a
successor of A. Let S'={s in S: label(u)(s)=alpha} and restrict to S'. Apply
N0 to u (I drops by 1). A now has constant value c; N3 deletes A (G drops by
1). C has a constant pin; N3 deletes C (G drops by 1). C's pins are A and some
node other than u and v (each of u, v has A as its only successor), so deleting C does not redirect
anything to v; deleting A (a constant) does not either. Thus v now has no
successor and is not the output, and N1 deletes v (I drops by 1). Total drop:
at least 4.

In both cases further normalization does not increase mu. QED.

### C.4 Proof of Theorem C

Let a B2 circuit with G_0 gates compute g. The start S-circuit (S=F_2^n,
m=n) has mu=G_0+n; normalize it (mu does not increase). While m>=d+2, apply
the Lemma and then normalize: m decreases by exactly 1 and mu by at least 4.
This runs for exactly n-d-1 steps (from m=n down to m=d+1; n>=d+2 makes the
count positive). Since mu>=0 at the end,
0 <= G_0 + n - 4(n-d-1), i.e. G_0 >= 3n - 4d - 4. QED.

Remarks on coverage. Repeated pins, constant nodes created by substitution,
unary/negated gates, dead gates, an output that is an input, a constant or a
negated node are all handled by N0–N6 and the free output bit beta, which only
records the fixed function g and is not a gate. Shared descendants and
coincidences among A, B, C are treated explicitly in Case 1 and Case 2. The
original circuit's gates are all counted in G_0; the generalized model is only
an accounting device and is at least as permissive as B2, so a lower bound on
its mu yields the stated lower bound on B2 size.

## Part D. Combining, with quantifiers

For n>=1 let d(n)=floor(6*n^(4/5)). By Part B, f_n is nonconstant above
dimension d(n) for every sufficiently large n (if the BK statement carries its
own large-n condition, take n beyond it as well; this only enlarges the
threshold below). For n>=248832=12^5 we have n^(1/5)>=12, so
6*n^(4/5) <= n/2 and d(n)+2 <= n/2+2 <= n. Theorem C applies at every such
sufficiently large integer n and gives

    C_B2(f_n) >= 3n - 4*floor(6*n^(4/5)) - 4 >= 3n - 24*n^(4/5) - 4.

Given epsilon>0, (24*n^(4/5)+4)/n <= epsilon for all n >= N_epsilon for an
explicit N_epsilon (e.g. N_epsilon = max(N_0, ceil((28/epsilon)^5)), where
N_0>=248832 is the threshold beyond which Part B applies), so
C_B2(f_n) >= (3-epsilon)n at EVERY integer n>=N_epsilon. Thus
liminf_{n->inf} C_B2(f_n)/n >= 3 over all integer lengths, which is the
declared polynomial_log guarantee c=3, p=1, q=r=0 and the finite projection
(3,1). The same fixed family, verifier and witness bound serve both
coordinates. Lengths below N_0 are finitely many and do not affect the
liminf. (Under the padded-length reading of Part B the same argument with
d(n)<=6*(n+k)^(4/5)<=12*n^(4/5) for large n gives loss <=52*n^(4/5)=o(n).)
No linear loss occurs: the only loss term 24*n^(4/5)+4 is o(n).

Check of the epsilon bound: for n >= (28/epsilon)^5, n^(1/5) >= 28/epsilon,
so 24*n^(4/5)+4 <= 28*n^(4/5) = 28*n/n^(1/5) <= epsilon*n.

## Sources

- BK-preprint: E. Ben-Sasson, S. Kopparty, "Affine dispersers from subspace
  polynomials", author preprint 2010-03-12; Theorem 2.2, equation (4), printed
  pp.6–7 (PDF pp.7–8); journal version SIAM J. Comput. 41(4), 2012;
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf;
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf
- Shoup-preprint: V. Shoup, "New algorithms for finding irreducible polynomials
  over finite fields", author preprint 1989-01-31; Section 3, Theorem 3.2,
  printed p.3 and fixed-prime consequence; https://www.shoup.net/papers/detirred.pdf;
  SHA256 00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644
- Model conventions only: Li–Yang, ECCC TR21-023 (2021-02-20), Sections 1.1,
  2.1, 2.4 (footnotes 3–4), 3.3. No theorem of this paper is used.
- Idea provenance only (proof reproduced above, not cited as a premise):
  Demenkov–Kulikov, MFCS 2011, LNCS 6907, pp. 256–265.
