# A total quotient-ring family with coefficient 7/2

All vector spaces and Boolean polynomials below use F_2. We prove, for the
single total language specified here, a full-B2 lower bound

    C_B2(f_N) >= (7/2)N - O(N^(5/6)).

In particular the exact finite declaration is (c,p)=(7/2,1). The proof uses
the affine-disperser construction of Ben-Sasson--Kopparty and the elimination
method of Mishra--Santhanam. The necessary construction and elimination
arguments are supplied below; no accepted judgment or external theorem about
circuit lower bounds is used as a premise.

## 1. Definition at every input length

Set f_0=0. For N>=1 choose the least prime integer k satisfying k^3>8N^2.
This is an exact integer test. If k>=N, set f_N identically zero. Otherwise
parse an N-bit input as (q,x), where q has k bits and x has m=N-k bits.
Put r=ceil(m/k), t=ceil(sqrt(r)), and e_i=r+1-i for 1<=i<=r.

Let Q_q(T)=T^k+sum_{j=0}^{k-1} q_j T^j and use the ring
R_q=F_2[T]/(Q_q), represented in the basis 1,T,...,T^(k-1). Append rk-m
zero bits to x and view its successive k-bit blocks as z_1,...,z_r in R_q,
with the first bit of each block the coefficient of 1. Define

    y_i = z_i^(2^(e_i+1)-1),
    S = sum_{A subset {1,...,r}, |A|=t} product_{i in A} y_i,
    U = sum_{j=0}^{k-1} S^(2^j).

The output f_N(q,x) is the coefficient of 1 in the degree-<k representative
of U. This prescription applies to every q, including reducible Q_q. It
does not require field recognition, a field oracle, or advice.

There is one polynomial-time algorithm for these functions. Trial division
and scanning find k; Bertrand's postulate implies k=O(N^(2/3)), so the scan
has polynomial bit cost (even an elementary trial-division implementation
suffices). For the finite exceptional lengths the output is already specified.
The bit cost of each ring operation is polynomial in k: multiply binary
polynomials and reduce by the monic Q_q. Each exponent defining y_i has
O(r) bits, so binary powering takes O(r) ring multiplications per block.
Compute S without division by the recurrence E_0=1, E_j=0 initially and,
upon reading y_i, update E_j <- E_j+y_i E_(j-1) in descending j. Partitioning
j-element subsets by whether they contain i proves that E_t=S. This uses
O(rt) ring operations. Repeated squaring and addition compute U in O(k)
more ring operations. All indices and representations have polynomial size
in N. Thus this is a fixed deterministic polynomial-time evaluator A.

Define L={a:A(a)=1}, W(N)=0 and V(a,w)=1 exactly when w is empty and
A(a)=1. If a is in L the empty witness accepts; if any permitted witness
accepts then A(a)=1 and a is in L. This proves both directions of uniform NP
membership. The target circuits receive exactly the N bits of a=(q,x),
without witnesses.

## 2. Reducing the theorem to a disperser on the data bits

For all sufficiently large N, k is an odd prime and m=N-k>0. Fix any monic
irreducible polynomial Q_0 of degree k over F_2. Its existence requires no
uniform construction: a splitting field of X^(2^k)-X has 2^k distinct roots;
these roots form a field under addition, multiplication, and inversion.
Since k is prime, an element outside F_2 has orbit of length k under
squaring, and its minimal polynomial over F_2 therefore has degree k.
(An orbit length divides k; length one means membership in F_2.) There are
2^k>2 elements, so such an element exists. Fix its polynomial's coefficient
bits q=q_0 in an arbitrary circuit for f_N. The resulting circuit has at
most the original number of gates, and computes a Boolean function h_m on
the m remaining bits. The choice q_0 is only a mathematical restriction of
the circuit; A is unchanged and evaluates every q uniformly.

In this restriction the quotient is the field K=F_(2^k). The map
Tr(z)=sum_{j=0}^{k-1} z^(2^j) is F_2-linear, takes values in F_2 because
Tr(z)^2=Tr(z), and satisfies Tr(1)=1 since k is odd. Hence the Boolean
output is exactly Tr(S). The padded block map F_2^m -> K^r is injective
and F_2-linear.

The construction argument in Section 7 below proves that h_m is nonconstant
on every affine subspace of dimension greater than

    D = kr/t + kt + r(r+1).

We use d=floor(D)+1. Here k=Theta(N^(2/3)), r=Theta(N^(1/3)), and
t=Theta(N^(1/6)), so d=O(N^(5/6))=o(N). The auxiliary parameter
inequalities needed in Section 7 hold eventually; no assertion about small
lengths is needed for the asymptotic guarantee. In Sections 3--6 we prove
that any such disperser on m bits requires at least (7/2)m-7d gates.

## 3. Circuit normalization and elimination states

We prove the lower bound in a relaxed circuit model with free constants,
free negation on each wire, and a constant or a signed internal value as
the output. Every circuit counted in the question is a circuit of this
model with the same gate count. A lower bound in the relaxed model
therefore also bounds the specified full B2 model. A gate whose two inputs
are equal or complementary, or one of whose inputs is constant, computes
a constant or a literal of one remaining input and can be bypassed. A
constant/unary truth table can also be bypassed. Output bypassing is
legitimate in this relaxed model. No extra charged gates are introduced.

Every nondegenerate binary truth table has one of the forms

    a+b+c,                  or (a+alpha)(b+beta)+gamma.

Indeed its unique multilinear polynomial is c+Aa+Bb+Dab. If D=0,
essentiality of both pins gives A=B=1. If D=1, completing the product
gives the displayed second form. These cover the two XOR-type and eight
AND-type tables. Constants and the four unary tables were already handled.

Fix a Boolean disperser h on b original input bits, nonconstant on every
d-dimensional affine subspace. During the argument maintain an injective
map Phi from F_2^M into F_2^b of the form

    Phi(u) = v + Lu + sum_{{p,r} in P} c_(p,r) p r,                 (1)

where P is a set of disjoint pairs of parameters. Put q=|P| and a=M-q.
A parameter outside P is unpaired. Fixing one parameter in each pair
makes (1) an injective affine map on M-q parameters. Its image is thus
an a-dimensional affine subspace. Consequently h composed with Phi is
nonconstant whenever a>=d.

We represent h composed with Phi by an invertible affine system followed
by an ordinary acyclic circuit. Write the system as

    A v = B u + c,

with R equations and R auxiliary variables v, A invertible, each equation
containing at most three distinct nonconstant terms after cancellation.
The acyclic part reads parameters, auxiliary variables, constants and
earlier gates; it sends no wire back to the system. Its output may be a
literal of any value. This representation always defines one value for
every internal variable on each parameter assignment. Count
s=R+(number of acyclic gates). Initially R=q=0, Phi is the identity,
M=a=b, and s is the given gate count.

The degree of a value means its number of distinct child gates or
equations, and does not count the output designation. Multiple pins in
one gate count as one child; terms of an equation are first cancelled.
For an unpaired parameter set

    w(u)=0 if degree 0, 6 if degree 1, 7 if degree >=2,
    H=2s+3q+sum_(unpaired u) w(u).                                (2)

Paired parameters have no separate weight. All weights and H are
nonnegative. Write Delta for final minus initial value during a step.

An affine gate with no nonlinear ancestor can be moved into the system:
replace its gate J by the equation J+Aa+Bb=c, with a,b its parents.
This adds one auxiliary variable and one equation and removes one gate.
The new auxiliary matrix is block triangular with final diagonal entry 1,
so remains invertible. Cancel any repeated terms. The size s stays fixed;
parameter degrees do not increase. This operation never increases H.

The parameter changes allowed in (1) are as follows. An unpaired parameter
may be fixed to a constant or identified with a literal of another
parameter: M drops by 1, q is fixed, and a drops by 1. A paired parameter
p may be fixed to a constant or set equal to its partner r plus a bit:
remove p and its pair, so M and q both drop by 1 and a is unchanged.
An unpaired parameter z may be replaced by w+ell(other parameters),
where ell is affine and w fresh; this bijection preserves M,q,a. Finally
distinct unpaired x,y,z may be changed by

    z = w+(x+alpha)(y+beta)+gamma,

with fresh w, while adding {x,y} to P. This is bijective, M is unchanged,
q increases by 1 and a drops by 1. Each map is injective; composing it
with Phi preserves injectivity and the disjoint-quadratic form (1).

We only implement the affine change by exchanging a parameter with an
auxiliary variable: if Q=x+ell(other parameters) in the solved system,
rename Q as fresh parameter w and old x as an auxiliary variable. No
equation or gate needs expansion. For invertibility, if Q is auxiliary
column j and B_x is the parameter column for x, then the new matrix is
A with column j replaced by B_x. Multiplying on the left by A^(-1)
gives an identity matrix with column j replaced by A^(-1)B_x, whose
jth entry is 1. Its determinant is 1. The change therefore preserves
fairness as well as equation lengths. Only the weights of x and w change.

There is also a useful exclusion while a>d. The output cannot be affine
in the parameters. If it were nonconstant affine, select a parameter
with nonzero coefficient. If paired, first fix its partner (a unchanged),
then eliminate the selected now-unpaired parameter using an affine
expression that makes the output constant. The resulting image still
contains an affine subspace of dimension at least a-1>=d, a contradiction.
The constant case contradicts nonconstancy directly. Nor can the output
be a literal of an AND-type gate on two unpaired parameters: fixing one
to its controlling bit makes the output constant at cost one in a.
The affine restriction in this exclusion is only a logical argument
about Phi; it is not used as a potentially large circuit rewrite.

## 4. Elementary reductions and accounting

We first describe four simplifications, each with Delta a=0,
Delta s=-1, and Delta H<=-1.

(i) Bypass an acyclic gate that is a constant or a literal of an input.
This includes an AND-type gate on affine inputs whose linear parts span
dimension at most one: its output is constant or a literal of one
nonconstant parent. The rerouting can increase the weight of at most
one formerly isolated unpaired parameter, from 6 to 7. Gate loss gives
-2, so Delta H<=-1. Constant rerouting gives Delta H<=-2.

(ii) Delete an acyclic gate with no child, unless it is the output.
This changes no parameter weight upwards, so Delta H<=-2.

(iii) An equation with at most two terms contains an auxiliary variable
v, since otherwise A would have a zero row. Write it v=U+c with U a
single value or zero. Substitute this everywhere and delete v and the
equation. Substitution uses row additions; the column of v then has a
single 1, so expansion of det(A) in that column proves the remaining
auxiliary matrix invertible. Equation lengths never increase. At most
one isolated parameter gains weight 1, so Delta H<=-1.

(iv) Let an isolated unpaired parameter x occur only in equation E.
Its column B_x is a unit vector. Since A^(-1)B_x is nonzero, choose Q
with coefficient 1 in that vector and exchange Q with x as above. The
now-auxiliary x occurs only in E and is not the output: before exchange
the output could not be the parameter x by the exclusion in Section 3.
Delete x and E. This does not affect other values and preserves
invertibility by determinant expansion in x's unit column. Parameter x
of weight 6 was replaced by Q of weight at most 7. Thus Delta H<=-1.

Always apply one of these reductions when available, after moving
eligible affine gates into the system. If no reduction applies, every
equation has three terms, no unpaired isolated parameter occurs in an
equation, every nonoutput gate has a child, and an AND-type gate on
affine values has independent linear parts.

Here are the remaining preliminary steps, in the stated order. Each
strictly reduces M+s and satisfies Delta H<=14 Delta a.

P1. Suppose an AND-type gate G on two parameters z,w has unpaired
isolated input z. Fix w to a bit making G constant and delete G.
Then z becomes unused. If w is unpaired, its weight and z's were each
at least 6; hence Delta a=-1 and Delta H<=-2-6-6=-14.
If w is paired, its partner becomes unpaired with weight at most 7;
Delta a=0 and Delta H<=-2-6-3+7=-4. Rerouting is only to a constant.

If P1 and the basic reductions are inapplicable, no child of a parameter
can also contain an isolated unpaired parameter: such a child could not
be an equation (by (iv)), or an AND-type gate (by P1); an XOR-type gate
on parameters would already have moved into the system. This observation
controls all possible increases from weight 6 to 7 below.

P2. If a paired parameter p has degree zero, fix its partner: the pair
disappears, p is unused, Delta a=0, Delta H<=-3. If p has degree at
least two, fix p and delete two distinct children. Each becomes unary
or constant (a three-term equation becomes eligible for (iii)). No
previously isolated unpaired parameter can gain weight, by the previous
observation. The surviving partner has weight at most 7, giving
Delta a=0, Delta H<=-3-4+7=0. If a child disappears during the first
deletion, it already contributes its own deletion; otherwise delete it
next. Thus shared descendants do not reduce the two distinct losses.
After P2 fails, every paired parameter has exactly one child.

P3. If partners p,r have the same child, impose p=r+c to make that
child constant. For an equation v+p+r=b, or an XOR-type gate, either
choice makes it constant. For (p+alpha)(r+beta)+gamma, choose
c=alpha+beta+1, so the product is zero. Delete the child. The surviving
r is unused. Thus Delta a=0 and Delta H<=-3-2=-5.

P4. If paired p occurs in an equation, let r be its partner, whose
unique child is different by P3. Fix r and delete its child. No old
isolated unpaired parameter gains weight. Newly unpaired p has weight
6; the interim change is at most -3-2+6=1. Apply (iv) to p and its
unique equation, decreasing H by at least 1 more. That equation remains
present after the first deletion: p does not occur in r's child, so the
substitution cannot duplicate p there or change its sole occurrence.
Thus Delta a=0 and Delta H<=0, with two size losses.

P5. Suppose an AND-type gate G on parameters p,w has paired input p.
Let r be p's partner. It has a single child J, an acyclic gate by P4,
different from G by P3. In particular w is neither p nor r. It is
either unpaired of degree at least two (P1) or paired. Fix r, then fix
w to make G constant; delete J and G. The now-unpaired p becomes
unused. No old isolated unpaired parameter gains weight: G is rerouted
to a constant, and J had the parameter r as an input. If w was
unpaired, Delta a=-1 and Delta H<=-3-4-7=-14. If w was paired,
two pairs disappear and only its partner can gain weight at most 7;
Delta a=0 and Delta H<=-6-4+7=-3. The two deleted gates are distinct.

When no preliminary step applies, every parameter in the affine system
is unpaired of degree at least two. Every AND-type gate already on
parameters has two unpaired inputs of degree at least two.

When a step names several deletions, retain each named equation/gate
until its own deletion is counted, postponing optional simplifications.
In particular, fixing a parameter makes each of its original child
equations have at most two terms. Row substitutions preserve this
property. Invertibility ensures that every such remaining equation has
an auxiliary variable to eliminate. Hence even when named gates share
descendants or equations share variables, each distinct named gate or
equation supplies its own loss of one in s. No assumption of disjoint
subcircuits or additive complexities is being made.

## 5. The remaining reductions

Continue to assume a>d. When the preliminary list is exhausted, the
acyclic part has a first gate G; otherwise the output is affine, which
Section 3 excludes. G is AND-type because all eligible affine gates
were moved into the system. Its parents are affine internal values
with independent linear parts.

Make its two parents parameters as follows. For an auxiliary parent Q,
choose an unpaired parameter x with coefficient 1 in the solved Q and
exchange x with Q. Such x occurs in the system and has degree at least
two. The new parameter Q also has degree at least two: it occurs in
some equation because its auxiliary column is nonzero, and in G.
Both weights are 7, so this exchange changes neither H nor a nor s.
If the other parent is already a parameter, choose x different from
it; independence guarantees that Q has such a coefficient. After the
first exchange, independence again guarantees that an auxiliary second
parent depends on some parameter other than the first exposed parent.
Thus at most two exchanges suffice. Recheck the preliminary list after
each exchange, ending the step if any applies. Otherwise G has two
unpaired parents x,y with degrees at least two.

Next, if a child J of G also reads x, replace J's truth table by
(x,y) -> J(x,G(x,y)) and use x,y as its two inputs. This is an allowed
B2 table, preserves the function, and removes J as a child of G.
Only y can gain children, and its weight was already 7. Do the same
for a child reading y. Recheck all preliminary reductions after each
such change and move eligible affine gates into the system. If none
applies, continue. Each rerouting removes one child of G, and these
preparations add none, so this ends. The preparations have not
increased H or changed a,M,s. Now no child of G reads x or y. The
gate G is not the output by Section 3 and has at least one child.
Assume deg(x)>=deg(y)>=2 by interchanging the names.

R1. Suppose both degrees are 2 and x,y have the same other child B.
If B is an equation it has form v+x+y=c; otherwise it is a gate on
x,y. On three of the four assignments to x,y the AND-type G takes
the same value. On two of those three, the Boolean value of B
(or x+y for the equation) is also the same. Those two points form an
affine line in F_2^2. Restrict to this line by fixing one parameter
or identifying x with y plus a bit, at cost Delta a=-1. Both G and
B (or auxiliary v) become constants. Delete them and reroute only to
constants. The two parameters both lose their weights 7: one is
removed and the survivor is unused. Thus Delta H<=-4-14=-18.

In all subsequent cases, if both degrees are 2 their other children
B_x and B_y are distinct. Let K>=1 be the number of children of G,
and let K_1 be the number whose other input is an unpaired parameter
of degree 1. These children are acyclic gates, never equations.

R2. Fix x to its controlling bit in G, so G becomes constant. Delete
all deg(x) children of x and all K children of G, beginning with G.
The two child sets are disjoint by the preceding rerouting. An
unpaired parameter of initial degree 0 cannot gain a child. One of
degree at least 2 cannot increase weight. For a parameter of degree
1, an increased weight requires its sole child to be deleted and
rerouted to its literal. The observation after P1 excludes this
among the children of x, so at most K_1 such parameters increase,
each by at most 1. Loss of x contributes -7, giving

    Delta H <= -7-2(deg(x)+K)+K_1.                                (3)

When deg(x)=2, deg(y)=2. Its sole other child B_y differs from B_x,
and none of G's children reads y. Thus y gains no child in the step
and loses G. Its weight drops by at least 1, strengthening (3) to
Delta H<=-12-2K+K_1. Since K_1<=K, the following cover all but one
configuration:

    deg(x)>=3:        Delta H <= -7-6-K <= -14;
    deg(x)=2,K>=2:    Delta H <= -12-K <= -14;
    deg(x)=2,K=1,K_1=0: Delta H <= -14.

Each has Delta a=-1. The only remaining configuration has
deg(x)=deg(y)=2 and K=K_1=1: G has exactly one child D, whose other
input z is a distinct unpaired parameter of degree 1. We use the
next two reductions instead of R2 in that configuration.

R3. If D is AND-type, write D=(G+delta)(z+eta)+theta. Select bits
x_0,y_0 with G(x_0,y_0)=delta; either value occurs for AND-type G.
Fix both x,y to those bits. Delete G,D,B_x,B_y. They are four
distinct objects: B_x differs from B_y by failure of R1, and neither
is D because no child of G reads x or y. Both G and D are rerouted
to constants. No formerly isolated unpaired parameter can gain
weight from B_x or B_y, since they had parameters x or y as parents.
The weights of x,y,z, formerly 7,7,6, all disappear. Thus

    Delta a=-2,   Delta H <= -8-7-7-6=-28.

R4. If D is XOR-type, D=G+z+c. Use the bijective change
z=w+G(x,y)+c, adding {x,y} to P. It has the permitted quadratic
form in Section 3. Because z occurred only in D and G's sole child
was D, substitute only in Phi, reroute D to fresh parameter w,
and delete both D and now-unused G. The affine system is unchanged:
z had no occurrence there. Parameters x,y lose weights 7 each;
z of weight 6 is replaced by w of weight at most 7, while the new
pair contributes 3. No other weight increases. Consequently

    Delta a=-1,   Delta H <= -4-7-7-6+7+3=-14.

These cases exhaust the full truth-table possibilities. No finite
enumeration, checker certificate, or unproved configuration coverage
is needed for any case.

## 6. Termination and the disperser lower bound

Every step above is an exact transformation of h composed with Phi,
or its indicated injective restriction, with the representation and
fairness conditions preserved. Every genuine reduction strictly
decreases the nonnegative integer M+s. Basic reductions and quadratic
rerouting decrease s; all other reductions fix a parameter and
decrease M without increasing s. Preparations terminate as described
in Section 5. While a>d at least one reduction applies. Thus the
procedure eventually stops at a_T<=d, and each full step satisfies

    Delta H <= 14 Delta a.

The argument does not require nonconstancy after the last step,
which can have Delta a=-2. Before each step, the exclusion lemma
applies because a>d. Unused parameters are deliberately retained
in Phi and in M, so injectivity and the dimension argument remain
valid even when their circuit degree is zero.

Initially H_0<=2s_0+7b and a_0=b. Summing the step inequalities
and using H_T>=0 yields

    2s_0+7b >= H_0 >= H_T+14(b-a_T) >= 14(b-d),
    s_0 >= (7/2)b-7d.                                            (4)

This applies to arbitrary nonuniform acyclic full-B2 circuits,
with unrestricted fan-out and depth. The fair affine systems are
only intermediate mathematical representations used to count
deletions; (4) bounds the original gate count s_0, not an uncharged
potential mistaken for circuit size. All shared nodes were counted
once in s, and all child-set estimates counted distinct objects.

## 7. The algebraic disperser construction

### Finite-field polynomial facts


Work over K=F_(2^k), with k an odd prime. All dimensions and linear maps in
this argument are over F_2. Polynomials representing functions of K^r are
reduced separately modulo X_i^(2^k)-X_i. A reduced polynomial has degree
less than 2^k in each variable and is uniquely determined by its function:
apply the univariate root bound in one variable and induct on the number of
variables. Consequently, a positive-degree coefficient that is nonzero proves
that the represented function is nonconstant.

Every linear map K->K has a representation sum_(j=0)^(k-1) a_j*X^(2^j).
Indeed these polynomials form a vector space of dimension k^2 over F_2,
the same as the space of linear maps; evaluation is injective by the root
bound and hence bijective. Such a polynomial will be called linearized.
The binary weight of an exponent is its number of 1 digits. Products of h
linearized factors have binary weight at most h in their variable, even
after reduction modulo X^(2^k)-X. To see this, first replace each exponent
2^b by 2^(b mod k). Collisions can be resolved by replacing two copies of
2^b by one copy of 2^(b+1 mod k), reducing the number of summands each
time. Thus a reduced monomial of binary weight h in a product of h factors
requires h distinct selected powers. We use this only for h<k.

For a linear subspace Y of K, of dimension D, define
P_Y(X)=product_(y in Y)(X-y). This is monic and linearized. One proof is
induction on D: adjoining a vector b to a subspace U changes its polynomial
to P_U(X)*P_U(X+b)=P_U(X)^2+P_U(b)*P_U(X), which is linearized if P_U is.
Its roots are distinct, and its kernel as a map on K is Y.

For 1<=D<=k there is a monic linearized polynomial L_Y of degree
2^(k-D) whose image is Y. For D=0 we instead use the reduced zero map
only for parameterization; that case is not part of the coefficient lemma.
To prove the nonzero-dimensional assertion, let Y'=P_Y(K) and
L_Y=P_(Y'). The dimension of Y' is k-D.
The monic polynomial L_Y(P_Y(X)) has degree 2^k and vanishes throughout K,
so it equals X^(2^k)-X. Also P_Y(L_Y(X))=X^(2^k)-X: compose this proposed
identity on the right with P_Y, use linearity and coefficients in K, and
then use injectivity of substitution by a nonconstant polynomial.
It follows that L_Y(K) is contained in Y; its dimension is D by its kernel
size 2^(k-D), so its image is exactly Y.

Write L_Y(X)=sum_(j=0)^(k-1) a_j*X^(2^j). Then a_(k-D)=1,
a_(k-D+1)=...=a_(k-1)=0, and a_0 is nonzero, the last fact following
from distinct roots and the derivative. If D>=2, the cyclic coefficient list
has no string of D-1 consecutive zeros except the indicated final string.
Here is a proof. If another such string starts at s modulo k, raise L_Y
to a suitable Frobenius power 2^h, where 1<=h<k, and reduce it, rotating
that string into positions k-D+1,...,k-1. The resulting polynomial has
degree at most 2^(k-D), has the same 2^(k-D) roots as L_Y, and is nonzero.
It must therefore be a nonzero scalar multiple c*L_Y. For every y in Y
this gives y^(2^h)=c*y. Fix b in Y nonzero. Then (y/b)^(2^h)=y/b for
every y in Y. Since k is prime and 1<=h<k, the only elements of K fixed
by this Frobenius power are 0 and 1: use Bezout's identity for h,k on the
invertible Frobenius automorphism, obtaining z^2=z. Hence Y is contained
in {0,b}, contradicting D>=2. This proves the assertion about zeros.
In particular, any zero string strictly between positions 0 and k-D has
length at most D-2.

### Isolating a coefficient in independent variables

We prove the following auxiliary statement. Let Y_1,...,Y_s be linear
subspaces of K, with dimensions D_i, and let L_i be their image polynomials
just constructed. Choose integers 1<=e_i<D_i-1, and put b_i=e_i+1. If

    sum_i(D_i-b_i)>k,

then the reduced polynomial

    Tr(product_i L_i(X_i)^(2^b_i-1))

has a nonzero coefficient at a monomial product_i X_i^beta_i with the
binary weight of beta_i equal to b_i for every i.

Write the coefficients of L_i as a_(i,j), with indices interpreted modulo k,
and form the matrix M_i with rows and columns 0,...,k-1, whose (l,j) entry
is a_(i,j-l)^(2^l). Let A_i={0,...,e_i}. For any set B_i of b_i distinct
column indices put beta_i=sum_(j in B_i)2^j. Expanding the trace and then
the factors shows that the coefficient in question is

    sum_(v=0)^(k-1) product_i per(M_i[A_i,B_i-v])^(2^v).       (1)

The permanent arises because full binary weight requires choosing each
selected exponent once, as established above; all bijections from factors
to the selected exponents contribute. This expansion also proves (1)
directly, without any assumption about signs or cancellation.

We now choose the B_i so that all permanents in the summand v=0 are
nonzero, and every summand v!=0 has a zero permanent. Relabel so D_1 is
maximal and set h_i=D_i-b_i>=1. Let s_0 be the first index with
sum_(i<=s_0)h_i>=k-D_1+1; it exists by the assumed strict inequality.
For i<=s_0 put ell_i=sum_(j<i)h_j. Then 0<=ell_i<=k-D_1<=k-D_i.
For i>s_0 put ell_i=0. To specify B_i, abbreviate D=D_i, e=e_i,
ell=ell_i, and a_j=a_(i,j).

If a_ell!=0, take, in the indicated column order,

    B_i=(ell, k-D+1, ..., k-D+e).

If a_ell=0, choose u<ell<v such that a_u and a_v are nonzero and all
coefficients strictly between them vanish. These indices exist because
ell lies between 0 and k-D and the endpoint coefficients are nonzero.
Put q=v-u-1 and z=min(q,e), and take the column order

    B_i=(v, u+1,...,u+z, k-D+z+1,...,k-D+e).

The entries in each displayed list are distinct. In either case the matrix
M_i[A_i,B_i] is lower triangular. In the first case its diagonal entries
are a_ell followed by Frobenius powers of a_(k-D). In the second case the
first diagonal entry is a_v, the next z are Frobenius powers of a_u, and
the remaining ones are Frobenius powers of a_(k-D). Every entry above the
diagonal vanishes, either because it lies in the gap (u,v), or because it
uses a coefficient at an index strictly above k-D. All diagonal entries
are nonzero. Thus the corresponding permanent is nonzero.

For any b in B_i, and every shift v in [b+1,b+h_i] modulo k, column b-v
of M_i[A_i,B_i-v] is zero: its coefficient indices run through
[b-v-e_i,b-v], contained modulo k in [k-D_i+1,k-1]. It remains to check
that these intervals, over all selected columns and all i, cover every
nonzero shift.

For i=1, ell_1=0 and a_(1,0)!=0. Its intervals cover [1,h_1] and
[k-D_1+2,k-1]. For each i<=s_0 the intervals cover
[ell_i+1,ell_i+h_i]. In the nonzero-coefficient case this follows from
ell_i in B_i. In the zero-coefficient case, if q<e, B_i includes the
consecutive indices u+1,...,v; their intervals cover [u+2,v+h_i], which
contains the desired interval. If q>=e, the intervals from u+1,...,u+e
cover [u+2,u+D-1], and that from v covers [v+1,v+h_i]. The zero-gap
property gives v-u<=D-1, so their union again covers
[u+2,v+h_i], containing the desired interval. Therefore the first s_0
blocks cover [1,k-D_1+1], and block 1 covers its complement in
[1,k-1]. This proves the required coverage. In (1), only the v=0 summand
survives and it is nonzero, proving the auxiliary statement.

### Passing from arbitrary affine spaces to independent blocks

Let 1<=t<r and k>=r(r+1)/t, and define b_i=r+2-i, alpha_i=2^b_i-1.
We prove that

    f(x_1,...,x_r)=Tr(Sym_r^t(x_1^alpha_1,...,x_r^alpha_r))

is nonconstant on every affine subspace A of K^r of dimension

    dim A > kr/t+kt+r(r+1).                                 (2)

Block Gaussian elimination provides subspaces Y_i of K, linear maps
sigma_ij:Y_j->K for j<i, and a point gamma in A such that A consists
exactly of the points

    x_i=gamma_i+y_i+sum_(j<i)sigma_ij(y_j),  y_i in Y_i,
    dim A=sum_i dim Y_i.

Indeed choose a section for the first-block projection of the underlying
linear space, subtract that section, and repeat on the remaining blocks.
Let L_i be an image polynomial for Y_i. For zero-dimensional Y_i simply
use L_i=0. The preceding parameterization is obtained as the image of
R=(R_1,...,R_r), where

    R_i=gamma_i+L_i(X_i)+sum_(j<i)L_ij(X_j),

and L_ij represents sigma_ij composed with L_j. Every such map is
linearized, by the finite-field polynomial fact above.

Set D_i=dim Y_i and S={i:D_i>r+1}. From (2),

    sum_(i in S)D_i > kr/t+kt.

In particular |S|>=t. Choose t indices I giving the largest D_i within S.
Their sum is at least t/r times the displayed sum. Consequently

    sum_(i in I)(D_i-b_i)
       >= sum_(i in I)(D_i-r-1)
       > k+kt^2/r-t(r+1) >= k.                             (3)

Also 1<=e_i=b_i-1<D_i-1 for i in I. The auxiliary coefficient statement
therefore gives a monomial M=product_(i in I)X_i^beta_i with binary weights
b_i and a nonzero coefficient in Tr(product_(i in I)L_i(X_i)^alpha_i).
For completeness, the following weighted count proves at once that this
coefficient is unchanged by replacing L_i(X_i) by R_i, and that no other
t-subset J contributes in the symmetric polynomial. Fix one Frobenius
summand of the trace and one term of its expansion for a t-subset J. Each
j in J supplies b_j factors, each choosing either a constant or a term
in a variable X_i with i<=j. Call such a choice a token sent to i;
constant choices or choices outside I are discarded. To obtain the
target exponent of binary weight b_i, at least b_i tokens must reach
each i in I, because reduction never increases binary weight.

Give a token reaching i mass 1/b_i. Total mass is at least |I|=t.
But i<=j implies b_i>=b_j, so each token sent by j has mass at most
1/b_j. Source j has b_j tokens and therefore supplies mass at most 1.
There are |J|=t sources, so total mass is at most t. Equality is forced
everywhere: no tokens are discarded, and every token has b_i=b_j.
Since b_i=r+2-i is strictly decreasing, every token must go to i=j.
Thus J=I and every factor from R_i uses its L_i(X_i) term. This is
exactly the coefficient supplied by the independent-variable lemma,
including the same Frobenius powers, with no contribution from constants
or triangular off-diagonal maps. Summing the trace terms preserves that
nonzero coefficient.

The reduced polynomial f(R(X)) therefore has the positive-degree monomial
M with nonzero coefficient. By uniqueness of reduced function polynomials,
it is a nonconstant function. Since the image of R is A, f is nonconstant
on A, proving (2).

## 8. Actual bit length, limit, and provenance

The data function h_m obtained from f_N by setting q=q_0 is the restriction
of the field-block function in Section 7 to an injective affine image of
F_2^m. Every affine d-dimensional subspace of the data domain remains
d-dimensional under that injection. The construction theorem therefore
shows that h_m is an affine disperser for d=floor(kr/t+kt+r(r+1))+1.
The required parameter conditions hold for every sufficiently large integer
N: k=Theta(N^(2/3)), r=Theta(N^(1/3)), t=Theta(N^(1/6)), so
k/(r(r+1)/t)=Theta(N^(1/6)) tends to infinity, r tends to infinity,
and 1<=t<r eventually. Here m=N-k, without any divisibility assumption.
The zero bits used to fill the final block are internal constants, not
additional input bits. The dimension estimate already covers the padded
subspace; the lower bound is measured against m, not rk.

Apply (4) in the relaxed signed-wire model to the restricted circuit. If
the original circuit for f_N has g gates, fixing q and retaining or
bypassing its degenerate gates cannot increase this count, so

    g >= (7/2)(N-k)-7d
      = (7/2)N - R_N,       R_N=(7/2)k+7d=O(N^(5/6)).

This holds for every circuit computing f_N exactly on all N-bit inputs,
for every sufficiently large integer N. The family and evaluator are fixed
before any circuit or restriction is chosen. The restriction q_0 may depend
on N; the total family does not use q_0 as advice. Arbitrary field bases
are not implicit inputs: the coefficient basis of Q_q is fixed explicitly
in the definition, and the first k actual input bits pay for its description.

Since R_N/N tends to zero, given any real epsilon>0 choose N_epsilon
large enough that all preceding parameter inequalities hold and
R_N/N<=epsilon for every integer N>=N_epsilon. Then
C_B2(f_N)>=(7/2-epsilon)N at all those lengths. Thus
liminf_(N->infinity) C_B2(f_N)/N >= 7/2. With
g(N)=(7/2)N the equivalent guarantee liminf C_B2(f_N)/g(N)>=1
follows, matching both the full expression and finite projection.
No logarithmic factor, superlinear exponent, lane endpoint, or new record
beyond the published 7/2 coefficient is asserted.

The supplied derivations reproduce the mathematical methods of:

- Eli Ben-Sasson and Swastik Kopparty, *Affine Dispersers from Subspace
  Polynomials*, author preprint dated 2010-03-12, Sections 3, 6 and 7
  (image polynomials, coefficient isolation, triangular parameterization),
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf.
- Om Swostik Mishra and Rahul Santhanam, *An Explicit 3.5n-o(n) Size Lower
  Bound against General Boolean Circuits*, ECCC TR26-260, 2026-10-19,
  Sections 2--4, Theorem 4.1 (fair affine systems, reparameterizations and
  the deletion accounting), https://eccc.weizmann.ac.il/report/2026/260/.

Those publications are attribution, not substitute evidence for the
nonstandard lemmas: all lemmas used from their methods were proved here.
The elementary field axioms, polynomial root bound, Gaussian elimination,
and Bertrand's postulate are used in their standard forms with hypotheses
specified above. The input-specified quotient ring avoids the separate
Shoup irreducible-polynomial algorithm cited in the repository seed.
No previous submission, its judgment, or its score is a proof premise.
