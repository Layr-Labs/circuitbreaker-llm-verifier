# A full-B2 family bound with error constants 14 and 13

## Main theorem

There is one fixed total language L in P such that, writing F_N for its
N-bit Boolean function, every circuit over all sixteen binary Boolean
tables, with arbitrary sharing and fan-out and one output, satisfies

    C_B2(F_N) >= 4N - 14 N^(4/5)       for every integer N >= 2^20,
    C_B2(F_N) >= 4N - 13 N^(4/5)       for every integer N >= 2^25.

More precisely, for every N >= 2^20 the proof gives

    C_B2(F_N) >= 4N - E(N),

    E(N) = (8 sqrt(2)+1) N^(4/5) + 16 N^(3/5)
           + 12 N^(2/5) + 4 sqrt(2) N^(1/5) + 24.

The coefficient of N is four, and the polynomial exponent is one.
Both logarithmic exponents are zero. The two structured claims in
submission.json state this same complete family theorem.

The proof has four parts: a full-model gate bound, a fully proved algebraic
affine-disperser lemma, elementary prime/field support, and exact actual-
length accounting. No DD hypothesis, old transport surplus, direct-sum
circuit hypothesis, or unproved structural charge is used.

## One fixed total evaluator

On an input z of length N:

1. If N < 2^20, output zero.
2. Let b be the least positive integer with b^5 >= N^3. Thus b=ceil(N^(3/5)).
3. Let ell be the first prime strictly larger than b, found by trial division.
4. Put r=floor(N/ell)-1 and t=ceil(sqrt(r)), using integer arithmetic.
5. Interpret the first ell bits as coefficients p_0,...,p_(ell-1) of the monic
   polynomial P(X)=X^ell+sum_(j=0)^(ell-1) p_j X^j over F2.
6. Test P for irreducibility by checking
      X^(2^ell) == X mod P,  and  gcd(P,X^2-X)=1.
   If either test fails, output zero.
7. In K=F2[X]/(P), interpret the next m=r ell bits as r field elements
   x_1,...,x_r using the polynomial basis in the same coefficient order.
8. Put w_i=r+2-i and y_i=x_i^(2^w_i-1). Compute the elementary symmetric
   polynomial Sym_t(y_1,...,y_r) by the descending dynamic program
      E_0=1, E_1=...=E_t=0;
      for each y_i, for j=min(i,t),...,1:
          E_j <- E_j+y_i E_(j-1).
9. Output Tr_K/F2(E_t)=sum_(j=0)^(ell-1) E_t^(2^j).
10. The remaining rho=N-(r+1)ell bits are ignored. They still count in N.

The symmetric-polynomial update is correct because each subset either
contains the new index or does not; descending j preserves the preceding
iteration's E_(j-1).

An exponent 2^w-1 is evaluated without writing an exponential-size object:
start p_1=x and use p_(a+1)=p_a^2*x, for a=1,...,w-1.

This is one deterministic algorithm, not a family chosen separately for
different error constants. Part IV proves its polynomial running time,
that all required parameters exist, and that its hard restriction exists
at every sufficiently large integer length. Define the NP witness bound
W(N)=0. The fixed verifier V(z,w) accepts exactly when w is empty and
this evaluator returns one. If z belongs to L, the empty witness is
accepted. Conversely, an accepted witness implies that the evaluator
returns one, so z belongs to L. This proves both directions of the NP
definition without advice or witness bits as circuit inputs.

## Provenance of the algebraic template

The trace of an elementary symmetric polynomial of consecutive all-ones
field powers is the Ben-Sasson--Kopparty construction. Reference:
Eli Ben-Sasson and Swastik Kopparty, Affine Dispersers from Subspace
Polynomials, March 12, 2010, Sections 3, 6, and 7:
https://www.math.toronto.edu/swastik/affine-disperser.pdf

The supporting arguments needed here are proved below. The sharper
entropy selection uses delta_i=max(d_i-w_i,0), and a weighted allocation
argument excludes cross terms. The source citation is provenance, not an
unproved premise of the family theorem.

## Source credit and evidence scope

The certificate-prefix method and the complete algebraic proof developed
in submissions f4f2c464-3045-4c32-a6e4-32478c05f343 and
2633389c-f681-4e93-b3c8-a86793dcc031 supplied the starting construction.
The three-point product domain and gate case framework came from the
public note of submission 5e9fcf7f-8338-4729-9faa-a18994661af6. The affine
system representation and role exchanges follow the parameter method of
Mishra and Santhanam, ECCC TR26-260, 2026-10-19. All needed reductions
are proved in Part I; no published theorem is imported as a premise.

This proof changes the potential to s+4I+2q, proves the fourth AP2
deletion, sharpens the algebraic dimension estimate, and uses a new fixed
choice of input-length parameters. The source submissions and their
judgments are not evidence in place of the arguments below. Contribution
dependencies are declared separately; judgment reuse is not used.

This is a proof-only submission. No finite experiment or checker result
is a premise of the theorem. All proof evidence is in this file.

# Part I. Full-B2 gate bound from ordinary affine dispersal

Theorem. If f:F2^n -> F2 is nonconstant on every affine flat of dimension
greater than an integer d>=0, then C_B2(f)>=4n-8(d+1).

The potential is V=s+4I+2q. Its surplus is different from the surplus of
the earlier W=s+5I+3q histories; no earlier token or AND-mass inequality is
transferred to it.

## 1. State and domain

Keep M independent binary parameters. There are q disjoint parameter pairs; each pair {p,h} ranges over F_2^2 minus one excluded corner (e_p,e_h). All other parameters range freely. Let D be this product domain.

Keep an affine injection Phi:F_2^M -> F_2^n, defined and injective on the whole parameter cube. The represented computation equals f(Phi(u)) for u in D. Equality outside D is not required.

The computation consists of:

* a square invertible affine system A v = N u + c, where each row has at most three distinct nonconstant terms;
* acyclic B2 gates, with no gate output entering the affine system;
* a selected output that is a constant or a literal of a parameter, auxiliary, or gate value.

Let s count all retained rows and gates, one charge each. A parameter's children are the distinct rows and gates containing it; the output reference is not a child. Let I count used free (unpaired) parameters. Define

    a = M-q,     V = s+4I+2q.

Initially M=a=n, q=0, Phi is the identity, s=|C|, and V<=|C|+4n.

Fixing one endpoint of each pair to the complement of its excluded coordinate leaves its partner free. Such a safe fixing picks an affine a-flat inside D. Its image under Phi is an affine a-flat of the same dimension. Therefore a constant represented output is impossible whenever a>d.

While a>=d+2, an affine represented output is also impossible: on a safe a-flat it has a constant subflat of dimension at least a-1>d. The same is true of an AND-type output whose parents are affine values. On the safe flat, fix a nonconstant affine factor to its blocking value; if there is no such factor, the output is affine or constant. Thus a minimal AND cannot be the selected output during the active loop.

## 2. Legal coordinate operations

The following preserve the domain form and affine injectivity, with the corresponding affine substitution also made in Phi.

1. Replace a free coordinate x by w+ell(other coordinates), where ell is affine. This is a bijection on each free two-point fiber, even if ell uses paired coordinates. It adds no gate.
2. Fix a free coordinate to a constant: M and a each fall by one; q is unchanged.
3. Safely fix one paired endpoint: M and q each fall by one, the partner becomes free, and a is unchanged.
4. Restrict a pair to any affine line contained in its three-point domain: parametrize the line by one free coordinate; M and q fall by one, with a unchanged.
5. Restrict two free coordinates to three corners of their square by declaring them a new pair: q rises by one and a falls by one; M and Phi do not change.

No nonlinear change of parameter is used.

## 3. Affine movement and elementary deletions

After absorbing negations into literal edges/tables, a nondegenerate binary gate is XOR-type or AND-type. All other tables are constant or depend on at most one parent.

Move an affine gate having no AND ancestor into the affine system, in topological order. Add its output as an auxiliary and its defining equation as a row. Its old gate charge becomes its row charge; nothing is added to s. The enlarged auxiliary matrix is block triangular with new diagonal entry one. Its row has at most three terms. The move preserves M,q,a,I,V.

Use the following elementary deletions.

**D1.** Delete a constant or one-parent-literal gate and reroute its uses. This includes fixed pins and repeated pins. If an AND has affine parents with dependent linear parts, those parents are constant or equal up to a constant, so its output is likewise constant or a literal and can be deleted.

**D2.** Delete a gate with no children unless it is the selected output.

**D3.** A row with at most two nonconstant terms contains an auxiliary, since A is invertible. Solve it as v=U+c, where U is absent or one other value, substitute, and remove its row and auxiliary column. Substitution does not increase row width. Adding the solved row to each row containing v and expanding on its column proves that the remaining auxiliary matrix is invertible.

For D1 and D3, any parameter receiving a rerouted use was already used in the deleted object; thus I cannot increase. D2 also cannot increase I. Each deletion removes a distinct retained charge and decreases V by at least one.

**D4.** Suppose a free parameter x has sole child row E. The column N_x is nonzero, so some auxiliary Q has x-coefficient one in A^{-1}N_x. Write Q=x+ell(other parameters). Exchange x and Q in role: Q becomes the new free parameter, and x becomes auxiliary. Substitute x=Q+ell in Phi. The new auxiliary matrix replaces A's Q-column by N_x and is invertible, since

    det(A with Q-column replaced by N_x)
      = det(A) * (A^{-1}N_x)_Q = 1.

No equation's terms change, so row width is unchanged. The old x now occurs as an auxiliary only in E, and it is not the selected output by the active terminal test. Remove x and E; expanding along its single-occurrence column again proves invertibility. The old used free x is replaced by at most one used free Q, so I does not increase. One charge is removed, and V drops by at least one.

The same role exchange can prepare an auxiliary AND parent Q without deleting a row. In that case the old free x is used in the system, and the new free Q is used by the AND; I is unchanged. The other statistics are unchanged. This is an invertible free-coordinate substitution, not an encoder.

## 4. Pair cleanup

At the start of each iteration, first move all available affine gates into the system. Then apply elementary deletions whenever possible, followed by the first applicable pair rule below. After any completed action, restart the iteration. In particular, finish the required pair cleanup after affine movement and before preparing an AND source. Propagate all forced constants and literals using D1/D3. When counting, count distinct old objects: an object already deleted supplies its one charge, never a second charge. Fixing a pin or shortening a row eventually removes the affected object, even if other propagation removes it first.

**P0: an unused paired endpoint.** If p is unused, safely fix its partner h and remove their pair. The surviving free p is unused. Then

    Delta a=0, Delta q=-1, Delta I<=0,
    V_before-V_after >= 2.

Any further deletions improve the bound.

**P2: an endpoint with at least two children.** Safely fix that endpoint and remove its pair. Its two distinct affected children are removed. Only its partner can become a new used free parameter. Therefore

    Delta a=0, Delta s<=-2, Delta I<=1, Delta q=-1,
    Delta V <= -2+4-2=0.

After P0/P2 fail, every paired endpoint has exactly one child.

**PC: partners with the same child B.** If B is a row, its third term is an auxiliary and its output is an affine function of the two endpoints. If B is a gate, its output is a Boolean function of them. Among three allowed pair points, two give the same output value. Those two points form an affine line inside the allowed domain. Restrict to it and remove the constant B. Both old endpoints occurred only in B, so the surviving line coordinate is unused. Hence

    Delta a=0, Delta s<=-1, Delta I<=0, Delta q=-1,
    Delta V<=-3.

**PE: a paired endpoint p has sole child row E.** Its partner h has a different sole child B, since PC is absent. Safely fix h and remove B and the pair. The free p still occurs only in E. If E has already disappeared in propagation, it supplies the second distinct deletion. Otherwise apply D4 to remove E. At most one newly free parameter is used, including a replacement parameter introduced by D4. Thus

    Delta a=0, Delta s<=-2, Delta I<=1, Delta q=-1,
    Delta V<=0.

When these rules all fail, every paired endpoint has one gate child, partners have distinct children, and paired parameters have zero columns in N.

## 5. Prepare a minimal AND

After affine movement and the required elementary and pair cleanup, if no acyclic gate remains, the selected output is affine, contrary to the active terminal test. Hence some AND G has no AND ancestor. Its parents are affine values. Elementary cleanup ensures their linear parts are nonzero and independent.

Paired columns of N are zero, since PE is absent, so all auxiliaries depend only on free parameters. For an auxiliary parent Q of G, exchange it with a free parameter on which it depends. If the other parent is already a free parameter, independence supplies a nonzero coefficient outside that parameter; if it is paired, Q has no dependence on it. If both parents are auxiliary, after the first exchange independence supplies a free pivot distinct from the first new parent. Thus at most two exchanges make both parents parameters.

Recheck elementary and pair rules after each exchange and complete any resulting action. If none arises, continue preparing the fixed G. This needs no endless source reselection. G has a child: it is not the selected output by the terminal test and it is not dead by D2.

## 6. Source restrictions and their charges

The parents of G are distinct parameters. Every AND-type binary table has one exceptional corner and one majority value at its other three points.

### GP: both source inputs free

Let G read x,y and have k>=1 distinct gate children. Restrict (x,y) to its three majority points and record them as a new pair. G is constant. Remove G and all its k children by fixed-pin propagation. These are 1+k distinct old gates, even if a child also reads x or y or one child is an ancestor of another. Each of x,y was used and now leaves I; old free parameters cannot acquire a first use through rerouting.

    Delta M=0, Delta q=1, Delta a=-1,
    Delta I<=-2, Delta s<=-(1+k).

Consequently

    Delta V <= -(1+k)-8+2 = -7-k <= -8.

### AP1: one source input paired

Let G read paired p and free x; h is p's partner. Safely fix h and remove its pair. Then fix x to block G. The now-free p has no remaining use, since G was its sole child. Let B be h's old sole child. PC guarantees B!=G. Remove G,B.

If G has a child J!=B, delete J as the third distinct object. Otherwise B is G's sole child and reads both G,h. Both inputs become constant, so B becomes constant and cannot be the selected output: the resulting safe dimension is a-1>d. A child of B supplies the third distinct old gate. Acyclicity makes it different from G,B.

This count is independent of propagation order. If B first becomes a literal of G, its later rerouted uses are still fixed when G becomes constant.

    Delta M=-2, Delta q=-1, Delta a=-1,
    Delta I<=-1, Delta s<=-3.

Hence

    Delta V<=-3-4-2=-9<=-8.

### AP2: both source inputs paired; four deletions

Let G read p,r, whose respective partners are h,k. The pairs differ; otherwise PC would apply to G. Safely fix h,k and remove both pairs. Fix the newly free p to block G. Then r becomes unused, since its sole child was G. Let B_h,B_k be the old sole children of h,k. Neither equals G, since PC is absent.

**Case 1: B_h!=B_k.** Removing G,B_h,B_k gives three distinct charges. If G has a child outside {B_h,B_k}, it is a fourth. Otherwise choose a child inside that set, say B_h. Its parents are G,h, which both become constant, so B_h becomes constant. If it has a child outside {G,B_h,B_k}, that child is a fourth deletion. It cannot feed G, whose parents are p,r, or itself, by acyclicity. If its only remaining child is B_k, then B_k reads B_h and k, both constant, and becomes constant too. Nonconstant output and acyclicity force a child outside the three-object set, supplying the fourth deletion. A selected output inside this constant chain is impossible after restriction, since a-1>d.

**Case 2: B_h=B_k=B.** The gate B reads h,k, so after fixing them B becomes constant and does not read G. G reads p,r, so neither G nor B reads the other. Both constant nodes have a child, since the restricted output is nonconstant. If they have different children, those children together with G,B give four distinct deletions. Otherwise their common child reads both constants, becomes constant, and has a child distinct from these three by acyclicity and nonconstant output. This is again a fourth deletion.

Every counted object is an old charged gate that is actually removed by the selected restrictions and propagation. If a previous substitution removed it already, its existing deletion counts once. There is no addition of two alternative restrictions or of hypothetical deletion sets.

The state changes are

    Delta M=-3, Delta q=-2, Delta a=-1,
    Delta I<=0, Delta s<=-4.

Therefore

    Delta V<=-4-4=-8.

All source restrictions end at safe dimension a-1>=d+1>d. Thus each use of nonconstant restricted output is valid. Further constant, literal, or dead propagation removes charges and never introduces a first use of an old free parameter, so it only improves all bounds.

## 7. Coverage, termination, and summation

Run while a>=d+2. Every active iteration performs finite affine movement and elementary/pair cleanup, or prepares a fixed minimal AND using at most two legal exchanges and applies GP/AP1/AP2. These exhaust its possible parameter-parent types.

Every completed action preserves the representation, never increases M or s, strictly decreases the nonnegative integer M+s, and obeys

    V_before-V_after >= 8(a_before-a_after).

Preparatory moves have zero changes in these quantities and cannot create an infinite preparation loop by the fixed-source protocol above. The process therefore terminates at a_T<=d+1. Since V_T>=0,

    |C|+4n >= V_0
             >= V_T+8(n-a_T)
             >= 8(n-d-1),

which yields

    |C| >= 4n-8(d+1).

If the loop never starts, n<=d+1 and the claimed bound is nonpositive. This completes the proof for every n.

Equivalently, let U_0=n-I_0 and let S_8 be the actual sum of V-surpluses. The exact identity is

    |C| = 4n+4U_0+S_8+V_T-8a_T.

No retained gate or row is counted twice, and no charged object is created. All sixteen binary tables are covered by constant/literal cleanup and the XOR/AND classification. Sharing and fan-out are unrestricted throughout.


# Part II. Algebraic disperser with a sharper threshold

## 1. Precise theorem

Let ell be an odd prime, let K be any field of size 2^ell, and let r,t be integers with 1<=t<=r and r+1<ell. Put

    w_i = r+2-i,  alpha_i = 2^(w_i)-1,  i=1,...,r,
    W = sum_i w_i = r(r+3)/2.

Define the Boolean function on K^r

    f(x_1,...,x_r) = Tr( sum_{I subset [r], |I|=t} product_{i in I} x_i^(alpha_i) ),

where Tr(z)=sum_{j=0}^{ell-1} z^(2^j). Under any F2-linear identification K^r=F2^(r ell), f is nonconstant on every affine flat of dimension

    D > W + ell * max(r/t, t-1).                       (T)

In particular, for t=ceil(sqrt r),

    D > r(r+3)/2 + ell sqrt r                         (T')

suffices, since both r/t<=sqrt r and t-1<sqrt r.

The construction has a deterministic polynomial-time evaluator given an irreducible degree-ell modulus: compute x_i^(2^w_i-1) by repeated squaring/multiplication and compute the elementary symmetric polynomial by the standard descending dynamic program. No sum with exponentially many terms is evaluated explicitly.

## 2. Elementary field/polynomial facts

Let q=2^ell. Reduction modulo X_j^q-X_j gives a unique polynomial of degree<q in each variable representing the same function on K^r. Such a reduced polynomial is the zero function only if it is the zero polynomial: prove this first in one variable by the root bound, then induct on the number of variables. Thus a nonzero coefficient of a nonconstant reduced monomial certifies a nonconstant function.

A polynomial sum_{j=0}^{ell-1} a_j X^(2^j) is F2-linear on K. Conversely every F2-linear map K->K has exactly such a representation. There are q^ell=2^(ell^2) such polynomials and the same number of linear maps. The map from polynomials to functions is injective by the root bound (degree at most 2^(ell-1)<q), proving the converse.

The trace takes values in F2 because Tr(z)^2=Tr(z), using z^q=z. It is nonzero since Tr(1)=1 for odd ell.

For a positive exponent that is a sum of h powers of two, reduction modulo X^q-X does not increase binary weight. Reduce the power positions modulo ell. Whenever two units occupy one position, merge them into one unit in the next cyclic position. Each merge decreases the number of units by one and preserves the represented monomial function. The process terminates with distinct occupied positions. A nonzero exponent never becomes exponent zero. Therefore the reduced binary weight is at most h, and equality h occurs only when the chosen positions were distinct (after cyclic position reduction). A constant term supplies zero units.

## 3. Image-subspace polynomials and their coefficient gaps

For the zero subspace use the linearized map Q_0=0. In the following image-polynomial construction assume 1<=d<=ell.

For an F2-linear subspace V<=K of dimension d define

    P_V(X)=product_{v in V}(X+v).

This is a monic linearized polynomial of degree 2^d. Induct on d: if V=U directsum F2 b, then

    P_V(X)=P_U(X) P_U(X+b)
          =P_U(X)^2 + P_U(b) P_U(X).

Its roots are distinct, and its linear coefficient is nonzero.

Write W0=P_V(K), a subspace of dimension ell-d, and put Q_V=P_{W0}. Then

    Q_V(P_V(X))=X^q+X,

because both sides are monic degree-q polynomials with all elements of K as roots. Also P_V(Q_V(X))=X^q+X: composing their difference with the nonconstant P_V gives the zero polynomial, since the coefficients of P_V lie in K. A nonzero polynomial cannot become identically zero after composition with a nonconstant polynomial.

Consequently Q_V(K) is contained in V. Its kernel is W0 of dimension ell-d, so its image has dimension d and equals V. Write

    Q_V(X)=sum_{j=0}^{ell-d} a_j X^(2^j),

and extend a_j cyclically modulo ell by setting a_j=0 for ell-d<j<ell. Here a_0 and a_(ell-d) are nonzero, and Q_V has exactly 2^(ell-d) roots in K.

### Gap lemma

If d>=2 and ell is prime, the cyclic coefficient sequence of Q_V has:

1. no d consecutive zero entries;
2. no d-1 consecutive zero entries except the designated tail
   ell-d+1,...,ell-1.

For (1), a cyclic Frobenius shift of d consecutive zero coefficients to the highest positions would give a nonzero polynomial of degree at most 2^(ell-d-1), with the same 2^(ell-d) roots as Q_V, a contradiction. If d=ell, the claim simply says the coefficient sequence is not identically zero.

For (2), a second gap of length d-1 can be cyclically shifted into the designated tail by a nonzero shift k modulo ell. Thus

    Q'(X)=Q_V(X)^(2^k) mod (X^q-X)

has degree at most 2^(ell-d), and it has exactly the same roots as Q_V. Comparing leading coefficients and using the root bound gives Q'=c Q_V for a nonzero c. Hence every y in V satisfies y^(2^k)=c y. Choose beta in V nonzero. Then every z=y/beta satisfies z^(2^k)=z. Since gcd(k,ell)=1 and the Frobenius automorphism has its ell-th power equal to the identity, a Bezout combination gives z^2=z. Thus z is in F2, so V is contained in beta F2 and has dimension at most one, a contradiction.

In particular every zero gap lying strictly between two nonzero coefficients in indices 0,...,ell-d has length at most d-2.

## 4. The independent-product coefficient lemma

Let V_1,...,V_s<=K have dimensions d_i, and let integers w_i satisfy

    2<=w_i<d_i.

Let Q_i be their image-subspace polynomials. Put

    e_i=w_i-1,  h_i=d_i-w_i>=1,  b_i=ell-d_i.

If sum_i h_i>=ell, then for every mu in K nonzero the reduced polynomial

    R(X_1,...,X_s)=Tr(mu product_i Q_i(X_i)^(2^w_i-1))

has a nonzero coefficient on a monomial product_i X_i^(beta_i) in which the binary weight of every beta_i is exactly w_i.

### Coefficient formula

Write Q_i(X)=sum_j a_(i,j) X^(2^j), with indices modulo ell. For a set B_i of w_i positions, beta_i=sum_{b in B_i}2^b, define M_i by

    M_i(u,b)=a_(i,b-u)^(2^u),  0<=u<=e_i.

Expanding the trace, the target coefficient is

    sum_{tau=0}^{ell-1} mu^(2^tau)
        product_i Perm M_i([0,e_i], B_i-tau)^(2^tau).   (C)

Indeed Q_i^(2^w_i-1) is a product of w_i Frobenius powers of Q_i. To obtain an exponent of weight w_i, one must choose distinct power positions; the assignments of the w_i factors to B_i are exactly the permutations in the permanent. Terms with repeated positions have smaller weight by the cyclic-carry fact. The trace rotates every chosen position by tau and raises its coefficient to 2^tau, giving (C).

We choose B_i so that the tau=0 product is nonzero, while every other product has a factor with a zero column.

### Triangular minors anchored at a desired position

Fix one i and an anchor z with 0<=z<=b_i.

If a_(i,z) is nonzero, take the ordered columns

    B_i=(z, b_i+1,...,b_i+e_i).

The minor is lower triangular. Its first diagonal entry is a_(i,z); its other diagonal entries are Frobenius powers of a_(i,b_i). All entries above the diagonal lie in the zero tail.

If a_(i,z)=0, let u<z<v be the immediately surrounding nonzero positions. Both exist because positions 0 and b_i are nonzero. Let g=v-u-1 and a=min(g,e_i). Take

    B_i=(v, u+1,...,u+a, b_i+a+1,...,b_i+e_i).

These are w_i distinct positions. For the first a+1 columns, the entries above the diagonal lie in the zero gap (u,v), and the diagonal entries are the nonzero a_(i,v), followed by Frobenius powers of a_(i,u). The remaining columns have their above-diagonal entries in the designated zero tail and their diagonal entries equal to Frobenius powers of a_(i,b_i). Thus this minor is also lower triangular with nonzero diagonal.

### Which shifts a selected column kills

A selected column k gives an entirely zero column in the tau-shifted minor whenever

    tau in [k+1,k+h_i] modulo ell.                    (Z)

For tau=k+a, 1<=a<=h_i, the coefficient indices are -a,-a-1,...,-a-e_i modulo ell. They all lie in the designated tail b_i+1,...,ell-1.

For the B_i anchored at z above, its intervals (Z) cover the whole integer interval

    [z+1,z+h_i].                                     (A)

If a_(i,z) is nonzero, the column z proves this immediately. If a_(i,z)=0 and g<e_i, the selected columns contain z, again proving it. If g>=e_i, the intervals from u+1,...,u+e_i cover [u+2,u+d_i-1]. The gap lemma gives v<=u+d_i-1, so adding the interval [v+1,v+h_i] gives an interval containing [z+1,z+h_i].

### Cover all nonzero shifts

Reorder the pairs (V_i,w_i) so d_1=max_i d_i, and let b_1=ell-d_1. Let s0 be the first index for which sum_(i<=s0)h_i>=b_1+1. Such an index exists since sum h_i>=ell.

For i<=s0, choose the anchor

    z_i=sum_(j<i)h_j.

Minimality gives z_i<=b_1<=b_i, so these anchors are valid. Their intervals (A) cover [1,b_1+1]. The first anchor is zero, and its selected columns are 0,b_1+1,...,b_1+e_1. The latter columns' intervals additionally cover [b_1+2,ell-1]. Together all nonzero shifts 1,...,ell-1 are covered. For i>s0, simply choose anchor zero; their minors remain triangular.

Every tau!=0 therefore has a zero column in some factor of (C), while tau=0 has all factors nonzero. The coefficient in (C) is nonzero. This proves the independent-product lemma.

## 5. Affine triangularization

Let A=gamma+V<=K^r be any F2-affine space of dimension D. Elementary coordinate-block elimination gives linear subspaces Y_i<=K and linear maps sigma_ij:Y_j->K for j<i such that

    x_i=gamma_i+y_i+sum_(j<i)sigma_ij(y_j),  y_i in Y_i,
    sum_i dimY_i=D.

For completeness, project V to its first block, choose a linear section of that projection, and apply the same construction to the kernel (whose first block is zero). Repeating gives the displayed triangular decomposition, with unique independent y_i coordinates and additive dimensions.

Let Q_i map K onto Y_i as in Section3, explicitly taking Q_i=0 whenever Y_i={0}. Every sigma_ij composed with Q_j is a linearized polynomial Q_ij in X_j. Therefore A is the image of the polynomial map whose coordinates are

    R_i(X_1,...,X_r)=Q_i(X_i)+sum_(j<i)Q_ij(X_j)+gamma_i.

No claimed circuit saving or free nonlinear encoder is involved: this is only a proof parameterization of a flat.

## 6. A weighted packet lemma isolates one symmetric summand

Weights w_i=r+2-i are positive and strictly decreasing with i.

Let I,J subset[r] each have t elements. Imagine that source j in J supplies w_j units and each unit can be sent only to an index i<=j (or discarded). Suppose every i in I receives at least w_i units. Give a unit received at i in I weight 1/w_i; ignore discarded units and all others. Required target mass is at least |I|=t. Every source-j unit has weight at most1/w_j, since i<=j implies w_i>=w_j. Total source mass is at most sum_(j in J)w_j/w_j=t.

All inequalities must be equalities. Strict decrease of the w_i forces every unit to go to i=j, with none discarded or sent outside I. Hence J=I, and every source uses only its own index.

This applies directly to reduced monomials in a trace of product_(j in J)R_j^(2^w_j-1). Expand it as w_j Frobenius-linear factors per sourcej. A chosen nonconstant term of such a factor is a scalar times a power-of-two monomial in some X_i with i<=j, so it supplies one unit there. A constant term supplies none. Polynomial reduction can only reduce a variable's binary weight. Thus a target monomial with weight w_i on i in I and no other variables requires at least w_i received units there. The packet lemma implies that only J=I can contribute, and then all factors must come from the own-variable term Q_i(X_i), without triangular cross terms or constants.

Consequently any nonzero independent-product coefficient with weight pattern (w_i)_(i in I) remains the same coefficient after triangularization and after adding all the other t-subset summands.

## 7. Sharper entropy selection and completion of the theorem

Let d_i=dimY_i, so sumd_i=D, and define

    delta_i=max(d_i-w_i,0).

Then

    sum_i delta_i >= D-W.

Under (T), this sum exceeds both ell*r/t and ell*(t-1). Each delta_i is strictly less than ell. Therefore at least t of the delta_i are positive. Choose I to contain t largest delta_i. Averaging gives

    sum_(i in I) delta_i >= (t/r)sum_i delta_i > ell.

For every selected i, d_i>w_i>=2. Thus the independent-product lemma applies to the selected Q_i,w_i. It supplies a nonzero reduced coefficient of a monomial M with weight exactly w_i in every selected variable. Section6 shows that the same coefficient occurs in the full composed trace-symmetric polynomial f(R_1,...,R_r), and cannot cancel with another summand.

That reduced polynomial is nonconstant, so f is nonconstant on A. This proves (T), including all affine translations and all cross-block dependencies. No direct-sum assumption about circuits is used.


# Part III. Elementary prime and field support

## 1. Bertrand's interval for every integer x >= 512

Let P(y) be the product of primes p <= y. Then P(y) <= 4^y for every real y >= 1.

It suffices to prove P(N) <= 4^N for integers N. The cases N=1,2 are immediate. For even N=2n >=4, 2n is composite, so P(2n)=P(2n-1)<=4^(2n-1). For odd N=2n+1, every prime in (n+1,2n+1] divides binomial(2n+1,n), since it appears in the numerator (n+2)...(2n+1) and not in n!. Therefore

  P(2n+1) <= P(n+1) binomial(2n+1,n) < 4^(n+1) 4^n.

The binomial inequality follows because the two middle coefficients of (1+1)^(2n+1) are equal and each is strictly less than half the sum 2^(2n+1). This completes the strong induction. For real y, apply the integer result to floor(y).

Now fix an integer x>=512 and suppose there is no prime in (x,2x]. Put B=binomial(2x,x). Since B is the largest of 2x+1 binomial coefficients summing to 4^x,

  B >= 4^x/(2x+1).

For any prime p, Legendre's formula gives

  v_p(B)=sum_{j>=1}(floor(2x/p^j)-2 floor(x/p^j)).

Every summand is 0 or 1, so p^v_p(B)<=2x. All primes p<=sqrt(2x) therefore contribute at most (2x)^sqrt(2x).

For p>sqrt(2x), the exponent in B is at most one. Primes with 2x/3<p<=x contribute zero, because the floors in the displayed formula are 2 and 1. There are no primes above x by supposition. Thus the remaining prime contribution is at most P(2x/3)<=4^(2x/3).

Combining bounds gives

  4^(x/3) <= (2x)^sqrt(2x) (2x+1),

or

  2/3 <= sqrt(2/x) log_2(2x) + log_2(2x+1)/x.

The right side is decreasing for x>=512. For the first term, differentiate (log(2x))/sqrt(x), whose derivative has sign 1-(log(2x))/2. For the second, differentiate log(2x+1)/x; its derivative is negative because log(1+t)>t/(1+t) for t>0, which follows by differentiating the difference from t=0.

At x=512 the right side is less than

  (320+11)/512 =331/512 <2/3.

This is a contradiction. Therefore every integer x>=512 has a prime ell with x<ell<=2x.

This proof uses no finite computational prime table. To select ell deterministically, test the successive integers x+1,...,2x by trial division up to their square roots and take the first prime. When x is polynomially bounded in the actual input length, this is a polynomial-time algorithm.

## 2. Existence of a monic irreducible polynomial of every prime degree ell

Take a finite splitting field for X^(2^ell)-X over F2. Existence is elementary: repeatedly adjoin a root of an irreducible factor by taking a polynomial quotient field, until the polynomial splits. Its derivative is 1, so it has 2^ell distinct roots.

The roots form a field: Frobenius shows closure under addition and multiplication, and if a nonzero root a satisfies a^(2^ell)=a, its inverse also satisfies that equation. Thus the root field has exactly 2^ell elements.

The Frobenius map a->a^2 has order dividing ell on this root field. For prime ell, each orbit has size 1 or ell. Its fixed points satisfy a^2=a and are exactly 0,1. Since 2^ell>2, an element with orbit size ell exists.

For such an element alpha, the product of (X-alpha^(2^j)) over j=0,...,ell-1 is fixed coefficientwise by Frobenius, so its coefficients lie in F2. Any nonzero F2 polynomial vanishing at alpha vanishes at all ell distinct orbit elements. Hence its degree is at least ell. The orbit polynomial is therefore the minimal polynomial, and is irreducible and monic of degree ell.

## 3. A polynomial-time prime-degree irreducibility test

For a monic P(X) in F2[X] of prime degree ell, the following conditions are necessary and sufficient for irreducibility:

  X^(2^ell) == X mod P,
  gcd(P, X^2-X)=1.

Necessity: F2[X]/(P) is a field with 2^ell elements. Every nonzero element has exponent dividing 2^ell-1 by the elementary group-order theorem, so the first identity holds. P has no linear factor.

Sufficiency: the first condition says P divides X^(2^ell)-X, which is squarefree and whose roots are the field of 2^ell elements constructed above. Every irreducible factor therefore has degree 1 or ell, by the Frobenius-orbit argument. The gcd condition excludes degree 1. Since deg(P)=ell, P must itself be irreducible.

Compute the first condition by ell successive polynomial squarings modulo P. With naive binary polynomial arithmetic this takes O(ell^3) bit operations. Euclid computes the gcd in polynomial time. Thus parsing an ell-bit coefficient prefix as

  P(X)=X^ell+sum_{i=0}^{ell-1} p_i X^i

and rejecting invalid prefixes is polynomial time. For a valid prefix, quotient-field addition, multiplication, squaring, and absolute trace are all polynomial time as well. An exponent with O(r) binary digits can be evaluated by O(r) squarings and multiplications.

For the lower bound, existence lets one FIX the prefix to any valid P. This is a restriction of one uniform total evaluator, not advice to the evaluator. All ell prefix bits must still be counted in actual input length.


# Part IV. Complete actual-length analysis

## 1. Parameter hypotheses and uniformity

For N>=2^20, b=ceil(N^(3/5))>=4096, so Part III gives

    b < ell <= 2b.

The chosen ell is odd. The integer b can be computed by binary search
using the exact comparison b^5>=N^3. Trial division finds ell using
polynomially many operations in N. No real-number oracle is used.

Let A=ell/N^(3/5). Since N^(3/5)>=4096,

    1 < A <= 2+2/N^(3/5) <= 33/16.

Hence N/ell >= (16/33) N^(2/5) >= 4096/33 > 124,
so r=floor(N/ell)-1>=123. Also

    r+1=floor(N/ell) <= N/ell < N^(2/5) < ell.

Thus 2<=t=ceil(sqrt(r))<r, and all hypotheses of Part II hold.

For a valid prefix, field arithmetic is polynomial time. With naive binary
polynomial arithmetic, each multiplication or modular squaring costs
O(ell^2) bit operations. The irreducibility test uses ell modular squarings
and Euclid. The field-power computations use O(r^2) field operations, the
symmetric dynamic program uses O(rt), and trace uses O(ell). Since r ell<=N,
all these costs are polynomial in N. Invalid prefixes are handled by the
same polynomial-time algorithm and return zero. The language is total.

The prime-degree irreducible-polynomial existence proof in Part III supplies
at least one valid prefix at every such length. For the lower bound, fix
the prefix to any one valid P and fix the rho unused bits arbitrarily.
Every remaining m=r ell-bit payload is a legal tuple in K^r, and the
restricted function is exactly the Part II disperser on its entire cube.

This fixed-prefix choice is a restriction in the lower-bound proof. It is
not advice to the evaluator. The disperser theorem holds for every valid
modulus, so no special modulus or hidden selector is needed.

## 2. The hard restriction and one-count circuit accounting

Put

    T = r(r+3)/2 + ell max(r/t,t-1),
    d = floor(T).

Part II proves that the payload function is nonconstant on every affine
flat of dimension greater than d: an integer dimension greater than
floor(T) is greater than T. Part I therefore gives

    C_B2(payload) >= 4m-8(floor(T)+1).

Restricting a full-B2 circuit for F_N to the fixed prefix and unused bits
never increases its number of gates. Thus

    C_B2(F_N) >= 4m-8(floor(T)+1).                    (Exact bound)

This is one gate count for one restricted circuit. There is no addition
of overlapping deletion histories or separately counted encoder gates.
The nonlinear Q_i maps in Part II parameterize affine flats only for the
algebraic proof; they are not uncharged circuit reparameterizations.

Since t=ceil(sqrt(r)), both r/t<=sqrt(r) and t-1<sqrt(r). Hence

    T <= r(r+3)/2 + ell sqrt(r).

Write N=(r+1)ell+rho, where 0<=rho<ell. All prefix and padding bits are
accounted for by

    N-m=ell+rho<2ell.

The exact bound consequently implies

    C_B2(F_N) >= 4N - [8ell+4r^2+12r+8ell sqrt(r)+8]. (1)

## 3. A sharper explicit error expression

Let u=N^(1/5), so N=u^5, ell=A u^3, and u>=16. We have

    1 < A <= A_max=2+2/u^3,
    r <= N/ell = u^2/A <=u^2.

The sum of the two leading terms in (1) obeys

    4r^2+8ell sqrt(r) <= u^4 g(A),
    g(A)=4/A^2+8sqrt(A).

Its derivative is

    g'(A)=4 A^(-1/2)-8 A^(-3),

which changes from negative to positive exactly once. The maximum on
[1,A_max] is therefore at an endpoint. At A=1, g(1)=12. At A=A_max>=2,

    g(A_max) <= 1+8sqrt(2) sqrt(1+u^(-3))
               <= 1+8sqrt(2)+4sqrt(2)u^(-3),

using sqrt(1+x)<=1+x/2. The latter bound also exceeds 12, so it bounds
both endpoints and hence g(A) throughout the interval.

Also 8ell<=16u^3+16, and 12r<=12u^2. Substitution in (1) proves

    C_B2(F_N) >= 4N-E(N),

    E(N)=(8sqrt(2)+1)N^(4/5)+16N^(3/5)+12N^(2/5)
          +4sqrt(2)N^(1/5)+24.                       (2)

This holds for every integer N>=2^20.

## 4. Error constants 14 and 13 with explicit thresholds

Use sqrt(2)<17/12, whose square is the inequality 288<289.

For u>=16, the leading coefficient in E(N)/u^4 is less than 37/3,
16/u<=1, and

    12/u^2 + 4sqrt(2)/u^3 + 24/u^4
      < 3/64 + 17/12288 + 3/8192
      = 1195/24576 < 1/16.

Therefore E(N)/N^(4/5)<37/3+1+1/16=643/48<14, proving

    C_B2(F_N)>=4N-14N^(4/5)  for every N>=2^20.

For u>=32 (equivalently N>=2^25), 16/u<=1/2, and

    12/u^2 + 4sqrt(2)/u^3 + 24/u^4
      < 3/256 + 17/98304 + 3/131072
      = 4685/393216 < 1/64.

Thus E(N)/N^(4/5)<37/3+1/2+1/64=2467/192<13, proving

    C_B2(F_N)>=4N-13N^(4/5)  for every N>=2^25.

The same single evaluator supplies both statements. In particular the
proof also gives the asymptotic form

    C_B2(F_N)>=4N-(8sqrt(2)+1)N^(4/5)-O(N^(3/5)).

## 5. Scope of the completed result

The family is fixed, total, and deterministically polynomial time; P is
contained in NP via the empty-certificate verifier. Every sufficiently
large integer input length is covered. All field-modulus prefix bits and
unused padding bits are included in N and its explicit error term.

The gate theorem covers arbitrary sharing, fan-out, and all sixteen binary
tables, with one output. Gates moved to affine rows retain their one old
charge. No encoder is made free and no alternate restriction histories
are added together.

For every real epsilon>0, let

    N_epsilon=max(2^20, ceil((14/epsilon)^5)).

Then every integer N>=N_epsilon satisfies

    C_B2(F_N)>=(4-epsilon)N.

Thus liminf over all integer lengths of C_B2(F_N)/N is at least four.
For the canonical guarantee g(N)=4N and relative tolerance eta>0,
use epsilon=4 eta. The evaluator and witness bound remain fixed as
these thresholds vary. The finite projection is exactly (4,1).
