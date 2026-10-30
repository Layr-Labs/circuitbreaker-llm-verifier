# A complete coefficient four lower bound for one uniform Boolean family

**Joint theorem.** The total family F_n defined in Section 1 is computed by
one fixed deterministic polynomial-time algorithm, and therefore defines one
language in NP with the empty-witness verifier specified there. In the full
Li-Yang B2 model of rules.md, at actual input length n, it satisfies

\[
\liminf_{n\to\infty}\frac{C_{B_2}(F_n)}{n}\ge4.
\]

More explicitly, for every integer n>=1024, the parameters k,m,r,t and d(n)
defined below satisfy

\[
 C_{B_2}(F_n)\ge4(n-k)-9(d(n)+1)
             \ge4n-\frac{627}{4}n^{4/5}.
\]

The full proof is in this file. The construction and affine-disperser
argument occupy Sections 1--5, the circuit lower bound is proved in Section
6, and Section 7 combines them and verifies the all-length quantifiers.
The whole n-bit function is not claimed to be an affine disperser: invalid
field headers give zero. The disperser property is used only for a fixed
valid-header restriction, with the header length explicitly charged.

The coefficient-four class argument is adapted from the public submissions
credited in Section 8. This package does not claim independent discovery
of that coefficient. The complete mathematical arguments are supplied here;
no prior score, judgment, audit, or acceptance is used as a premise.

This section is a self-contained version, over the binary field, of the
sublinear-dimensional affine-disperser argument of Ben-Sasson and Kopparty,
author preprint dated 2010-03-12, Theorem 7.3 and its supporting lemmas in
Sections 3, 6 and 7. The construction below uses descending powers and field
trace. Its field-modulus header is an input to one fixed polynomial-time
algorithm. No algorithm for uniformly constructing irreducible polynomials is
assumed.

## 1. Total definition and efficient evaluation

For every input length n<1024 define F_n=0. For n>=1024 let k be the least
prime satisfying k^5>32 n^3, put m=n-k, r=ceil(m/k), and t=ceil(sqrt(r)). All
these choices are deterministic integer computations. The first k bits of the
input are q_0,...,q_{k-1}, representing

    q(Z)=Z^k+sum_{j=0}^{k-1} q_j Z^j.

The last m bits are the data. If q is reducible over F_2, output zero.
Otherwise let K=F_2[Z]/(q), append kr-m zeros to the data, and interpret its r
blocks in the polynomial basis 1,Z,...,Z^{k-1}, with the first bit of a block
its constant coefficient. Denote the resulting elements by z_1,...,z_r.
For i=1,...,r set

    e_i=r+1-i,  w_i=e_i+1=r+2-i,  a_i=2^{e_i+1}-1,
    y_i=z_i^{a_i}.

Output

    F_n(q,data)=Tr_K/F_2(sum_{I subset [r], |I|=t} product_{i in I} y_i).

Trace means Tr(u)=sum_{h=0}^{k-1} u^{2^h}. Its value belongs to F_2, because
squaring permutes the summands and u^{2^k}=u. Thus this is a total Boolean
function.

Here are elementary details of the algorithms and their running times.
Bertrand's postulate, that an integer a>1 has a prime between a and 2a,
gives k<4 n^{3/5}+2. Thus searching consecutive integers and testing primality
by trial division is polynomial time in n. Also k>2 n^{3/5}. Since n>=1024
implies n^{2/5}>=16,

    k<4 n^{3/5}+2 <= n/4+2 < n/3.

Consequently m>0, r>=3, and t<r. In particular the definition needs no
nonexistent block or negative padding. The integer square root and fifth-power
comparison avoid any reliance on floating-point computation.

Because k is prime and greater than two, irreducibility can be tested by the
following particularly simple deterministic procedure. Require q(0)=q(1)=1;
compute Z^{2^k} modulo q by k modular squarings; accept q precisely when this
remainder equals Z. To justify the test, the polynomial Z^{2^k}-Z has no
repeated roots, since its formal derivative is one, and its roots in a
splitting field form the field of 2^k elements. Indeed the roots are closed
under addition, multiplication and inversion, and there are exactly 2^k of
them. An irreducible factor of this polynomial has degree dividing k: a root
generates a subfield, and the tower law for dimensions gives this divisibility.
As k is prime, its factors therefore have degrees one or k. Excluding the
two linear factors and requiring q to divide Z^{2^k}-Z forces q itself to be
irreducible. Conversely an irreducible degree-k q defines a field of size
2^k, whose elements satisfy u^{2^k}=u, and it has no linear factors. This
proves both directions of the test.

It also proves existence of an irreducible degree-k polynomial. Construct a
splitting field of Z^{2^k}-Z by adjoining roots; its 2^k roots form the
field just described. An element outside F_2 exists. Its minimal polynomial
has degree dividing the prime k and greater than one, so has degree k.
This is only an existence argument, not a step of the evaluator.

With schoolbook arithmetic, a multiplication or reduction of degree-k binary
polynomials has polynomial bit cost, for example O(k^2). The exponents a_i
have at most r+1 binary digits, so all y_i take O(r^2) field multiplications
by repeated squaring. To compute the elementary symmetric polynomial, start
E_0=1 and E_j=0 for 1<=j<=t, and for each y_i update

    E_j <- E_j+y_i E_{j-1}

in descending order j=min(i,t),...,1. Partitioning subsets according to
whether they contain i proves this recurrence. It uses O(rt) field
operations and involves no division in characteristic two. Trace takes k
successive squarings and additions. Thus one fixed deterministic algorithm A
evaluates F_n in polynomial time at every n. Let

    V(x,w)=1 iff w is empty and A(x)=1,    W(n)=0.

Both implications in x in L iff there exists |w|<=W(n) with V(x,w)=1 now
hold by the definition L={x:A(x)=1}. Thus L is one fixed language in P and
therefore NP. The modulus header belongs to x; it is neither advice nor a
witness.

## 2. Elementary polynomial facts

We give all nonstandard premises needed for the affine-disperser property.
In the rest of this section k is an odd prime, K is any field of size 2^k,
and Q=2^k. Every polynomial function on K^s has a unique representative of
degree less than Q in each variable, obtained by reduction modulo
X_j^Q-X_j. Uniqueness follows from the univariate root bound, applied one
variable at a time. Therefore a reduced polynomial is a constant function
exactly when it is a constant polynomial.

For an exponent b in [0,Q-1], write wt(b) for the number of ones in its
binary expansion. A monomial's binary degree in X_j is wt of that
variable's exponent. A product of h monomials of the form c X_j^{2^v},
after reduction modulo X_j^Q-X_j, has binary degree at most h. To see this,
reduce v modulo k, and whenever two equal powers occur replace them by the
next power, cyclically modulo k. Every such carry reduces the number of
terms by one. The process ends in distinct powers, representing the unique
positive reduced exponent; it cannot convert a positive exponent into the
zero exponent. Constants contribute zero binary degree. Consequently
Frobenius powers preserve binary degree, and multiplication or summation
cannot create a monomial whose binary degree exceeds the number of
linearized factors contributing to it.

A linearized polynomial is a finite sum of terms b_j X^{2^j}. Its
evaluation is F_2-linear. Conversely every F_2-linear map K->K has a unique
such polynomial using just j=0,...,k-1:
evaluation on the k^2-dimensional coefficient space is injective by the
root bound (degree at most 2^{k-1}<Q), and both coefficient space and the
space of F_2-linear maps have dimension k^2 over F_2.

For a d-dimensional linear subspace Y of K with 1<=d<=k, define

    P_Y(X)=product_{y in Y}(X-y).

It is monic linearized of degree 2^d. Induct on d: if Y=U+F_2 v with v
outside U, then

    P_Y(X)=P_U(X)P_U(X+v)=P_U(X)^2+P_U(v)P_U(X),

which preserves the linearized form. The roots are exactly Y and simple.
Let J=P_Y(K). It has dimension k-d. Set Q_Y=P_J. Then

    Q_Y(P_Y(X))=X^Q-X.

Both sides are monic degree-Q polynomials with every element of K as a
root. Moreover P_Y(Q_Y(X))=X^Q-X as well: composing the preceding identity
on the left by P_Y and using
P_Y(X^Q-X)=P_Y(X)^Q-P_Y(X) shows
(P_Y o Q_Y-(X^Q-X)) o P_Y=0; composition with a nonconstant polynomial is
injective. It follows that Q_Y(K) is contained in Y. Its kernel is J, so
rank-nullity gives equality. Thus Q_Y is a monic linearized polynomial of
degree 2^{k-d} whose image is Y and which has exactly 2^{k-d} roots in K.
We call it the image polynomial of Y.

Write its coefficients as Q_Y(X)=sum_{j=0}^{k-1} b_j X^{2^j}, with indices
understood modulo k. The following two facts about zero coefficients are
central:

1. No d consecutive cyclic coefficients can all vanish. For d=k this is
   immediate from Q_Y=X. For d<k, if they did, a
   Frobenius power of Q_Y, reduced modulo X^Q-X, would have degree at most
   2^{k-d-1} but retain its 2^{k-d} roots, a contradiction. In particular
   b_0 and b_{k-d} are nonzero, and b_{k-d+1},...,b_{k-1} are zero.
2. If d>1, the last displayed run is the only cyclic run of d-1 zero
   coefficients. An additional such run, rotated to the last d-1
   positions by a Frobenius power 2^v with 0<v<k, would give another
   degree-at-most-2^{k-d} polynomial with the same roots. It must be c Q_Y
   for a nonzero c. Thus y^{2^v}=c y for every y in Y. Fix a nonzero beta
   in Y. Every y/beta is fixed by the 2^v power map. Since this map and
   the 2^k power map fix y/beta and gcd(v,k)=1, the 2 power map does too
   (use Bezout's identity in the cyclic group generated by Frobenius).
   Therefore Y/beta is contained in F_2, contradicting d>1.

The first fact also covers d=k; in that case Q_Y=X. The second argument
uses only primeness of k and not a classification of arbitrary subfields.

## 3. A coefficient lemma for independent blocks

Let Y_1,...,Y_s be linear subspaces of K, with dimensions d_i>1. For each
i choose an integer e_i with 1<=e_i<d_i-1; write w_i=e_i+1 and
a_i=sum_{j=0}^{e_i}2^j. Let Q_i be the image polynomial of Y_i. Suppose

    sum_i(d_i-w_i)>k-max_i d_i+1.                 (1)

We will prove that the reduced polynomial

    R(X_1,...,X_s)=Tr(product_i Q_i(X_i)^{a_i})

has a nonzero coefficient at a monomial product_i X_i^{b_i} with
wt(b_i)=w_i for every i.

For Q_i(X)=sum_j c_{ij}X^{2^j}, define a k-by-k matrix

    M_i[u,v]=c_{i,v-u}^{2^u},

where indices are modulo k. This entry is the coefficient of X^{2^v}
in Q_i(X)^{2^u} reduced modulo X^Q-X. Put A_i={0,...,e_i}. For a set
B_i of w_i distinct indices put b_i=sum_{v in B_i}2^v. The coefficient
of product_i X_i^{b_i} in R is

    sum_{h=0}^{k-1} product_i Perm(M_i[A_i,B_i-h])^{2^h}.          (2)

Indeed expand trace into its k Frobenius powers and expand each a_i into
its w_i distinct binary powers. A term with maximal binary degree w_i
must choose w_i distinct exponents from those w_i factors: a constant or
a repeated exponent would reduce binary degree. Bijections between the
factors and the chosen exponents give exactly the permanent. This also
justifies formula (2) in characteristic two without assuming that reduced
exponents behave like ordinary degrees.

It remains to select the B_i so that the h=0 summand of (2) is nonzero
and every other summand is zero. Permuting the independent factors, assume
d_1=max_i d_i. Write h_i=k-d_i and s_i=d_i-e_i-1>0. Let T be the first
index with sum_{i<=T}s_i>=k-d_1+1; it exists by (1). For i<=T put
l_i=sum_{j<i}s_j. Thus

    l_1=0,   0<=l_i<=k-d_1<=h_i,
    sum_{i<=T}s_i <= k-d_1+s_T < k.

For each i<=T choose B_i as follows.

If c_{i,l_i} is nonzero, take

    B_i={l_i,h_i+1,...,h_i+e_i},

in that column order. If c_{i,l_i}=0, let b<l_i<j be the nearest nonzero
coefficient indices on either side; they exist between 0 and h_i. Put
q=j-b-1 and v=min(q,e_i), and take, in the displayed order,

    B_i={j} union {b+1,...,b+v} union {h_i+v+1,...,h_i+e_i}.

These are distinct indices. For i>T take the simpler choice
B_i={0,h_i+1,...,h_i+e_i}. With the rows in order 0,...,e_i, every
M_i[A_i,B_i] is lower triangular with nonzero diagonal. In the first
choice the diagonal coefficients are c_{i,l_i},c_{i,h_i},...,c_{i,h_i};
all entries above the diagonal use the terminal zero run. In the second
choice the first diagonal coefficient is c_{i,j}, the next v diagonal
coefficients use c_{i,b}, and the remaining ones use c_{i,h_i}; entries
above the diagonal use either the gap b+1,...,j-1 or the terminal zero
run. Frobenius powers preserve nonzero values. Hence the h=0 term of
(2) is nonzero.

For a column index z in B_i, the shifted column in M_i[A_i,B_i-h] is
all zero whenever

    h in [z+1,z+s_i] modulo k,                                (3)

because every coefficient index z-h-u, for 0<=u<=e_i, then lies in the
terminal zero run h_i+1,...,k-1. We check that these intervals cover all
h=1,...,k-1. From B_1 they cover

    [1,s_1] union [k-d_1+2,k-1].

For 1<i<=T, they cover [l_i+1,l_i+s_i]. If c_{i,l_i} is nonzero,
this follows from l_i in B_i. If it is zero and q<e_i, the gap indices
together with j are in B_i, including l_i. Finally suppose q>=e_i.
The nonterminal gap has q<=d_i-2 by the second zero-run fact, so
j=b+q+1<=b+d_i-1. The intervals (3) from b+1,...,b+e_i cover
[b+2,b+d_i-1], reaching j; the interval from j continues through
j+s_i. As b<l_i<j, this union contains [l_i+1,l_i+s_i]. Thus the
first T factors cover [1,sum_{i<=T}s_i], which joins the interval
[k-d_1+2,k-1] from the first factor. Every nonzero shift h has an
all-zero column in at least one factor. Its term in (2) vanishes, proving
the coefficient lemma.

## 4. The affine-disperser property for dependent blocks

Let r>=3, t=ceil(sqrt r)<r, and suppose k>=r(r+1)/t. Define

    H(z_1,...,z_r)=Tr(sum_{I subset [r], |I|=t}
                        product_{i in I} z_i^{2^{r+2-i}-1}).

We prove H is nonconstant on every F_2-affine subspace A of K^r with

    dim A > D := kr/t+kt+r(r+1).                              (4)

Choose gamma in A. Gaussian elimination, grouping basis vectors by their
first nonzero block, gives linear subspaces Y_i of K and linear maps
sigma_ij:Y_j->K for j<i such that the points of A have the unique form

    z_i=gamma_i+v_i+sum_{j<i}sigma_ij(v_j),   v_i in Y_i,
    dim A=sum_i d_i,   d_i=dim Y_i.

For completeness, each group of basis vectors has an independent set of
leading values in its first block; let Y_j be their span. The later
blocks of the same basis vectors define the maps sigma_ij. Distinct
groups together are a basis, proving both the parametrization and its
dimension statement.

For d_i>0 let Q_i be the image polynomial of Y_i, and for Y_i={0}
put Q_i=0. Extend each sigma_ij linearly from Y_j to K. The polynomial
representing sigma_ij o Q_j is linearized, by the elementary
representation proved in Section 2. We consequently have a polynomial
parametrization A={ (R_1(x),...,R_r(x)):x in K^r }, where

    R_i(X)=Q_i(X_i)+sum_{j<i}Q_ij(X_j)+gamma_i.

Each R_i involves only X_1,...,X_i, and each of its nonconstant monomials
has binary degree one.

Let S={i:d_i>r+1}. Equation (4) implies

    sum_{i in S} d_i > kr/t+kt,

since the omitted dimensions total at most r(r+1). In particular S has
at least t members. Choose the t largest d_i among S, with index set I.
Their sum is at least (t/r)sum_{i in S}d_i, so

    sum_{i in I}(d_i-(r+1))
        > k+kt^2/r-t(r+1) >= k.                             (5)

The last step uses kt>=r(r+1). For i in I take e_i=r+1-i and
w_i=r+2-i. Then 1<=e_i<d_i-1 and, by (5),
sum_{i in I}(d_i-w_i)>k. The coefficient lemma supplies a monomial

    M=product_{i in I} X_i^{b_i},   wt(b_i)=w_i,

with nonzero coefficient in Tr(product_{i in I}Q_i(X_i)^{a_i}).

That coefficient is unchanged on replacing Q_i(X_i) by R_i(X). To see
this, expand trace and each a_i into its w_i binary powers, and begin
with the largest i in I. Earlier-index factors do not involve X_i.
Obtaining binary degree w_i in X_i therefore requires taking a term in
X_i from every one of the w_i factors belonging to R_i; these terms
come exactly from Q_i(X_i). Having selected them, descend to the next
largest index. This proves the assertion for each trace summand, hence
for their sum.

No subset J of size t other than I can contribute the same monomial M
to H(R_1,...,R_r). Here is a direct argument that also rules out all
cancellations between different subsets. An expanded term from subset
J selects one constant or linearized monomial from each of its w_j
factors for each j in J. Let a_ji count the selections whose variable
is X_i. Only i<=j can occur; sum_i a_ji<=w_j. In order to yield M,
binary-degree subadditivity requires sum_{j in J}a_ji>=w_i for each
i in I. Variables outside I cannot occur because positive exponents
never reduce to zero. As w_i=r+2-i, the condition i<=j implies
w_i>=w_j. It follows that

    t <= sum_{i in I} sum_{j in J} a_ji/w_i
      <= sum_{j in J} sum_{i in I} a_ji/w_j <= t.

All inequalities must be equalities. Any positive a_ji with i<j would
make the middle inequality strict. Thus all selections have i=j.
Every i in I needs positive supply, so I is contained in J; both have
size t, hence J=I. The coefficient of M in the complete reduced
polynomial H(R_1,...,R_r) is consequently nonzero. It has positive
degree, so Section 2's uniqueness fact proves that H is nonconstant
on A. This establishes (4), including arbitrary dependencies between
blocks.

## 5. Parameters, actual input length and restriction

For the parameters of Section 1,

    r=ceil(n/k)-1 < n/k < (1/2)n^{2/5}.

Since t>=sqrt r,

    r(r+1)/t <= r^{3/2}+sqrt r
      < (1/(2 sqrt 2))n^{3/5}+(1/sqrt 2)n^{1/5}
      < 2 n^{3/5}<k.

We already proved r>=3 and t<r. Thus Section 4 applies for every
n>=1024. Define the integer

    d(n)=ceil(kr/t)+kt+r(r+1)+1.

It is strictly greater than D. For any irreducible header q, padding
the m data bits defines an injective F_2-linear map from F_2^m to K^r.
It preserves affine dimensions. The restricted data function therefore
is nonconstant on every affine subspace of dimension at least d(n).
This statement is vacuous when d(n)>m, which causes no asymptotic issue.

The estimates k=Theta(n^{3/5}), m=n-k=Theta(n),
r=Theta(n^{2/5}), and t=Theta(n^{1/5}) imply

    d(n)=O(n^{4/5}),   k/n -> 0,   d(n)/n -> 0.

These statements hold along all integer lengths, not a subsequence.
An explicit prime exists and an irreducible header exists at every
n>=1024, as proved above. For each such length fix any one irreducible
header only for the purpose of restricting a hypothetical circuit.
The evaluator never makes this choice and never receives advice.

Substituting the fixed header bits cannot increase the number of gates
in a full B2 circuit. One may first retain the same gates with constants
on some pins, then process them in topological order. A constant output
is propagated. A projection is bypassed. A negation is absorbed into
the truth tables of all immediate descendants; fan-out causes no gate
cost. If its output is the final output, it can instead be retained as
the existing unary gate, or absorbed into the preceding gate. A unary
gate on a lone remaining input can use that input on both pins, again
retaining at most the existing gate. Repeated pins and coincident
subgraphs are allowed throughout. Thus the restricted data function
has circuit complexity no larger than C_B2(F_n).


## 6. A complete gate-elimination proof on three-point domains

### C.0. Statement and credit

**Theorem.** Let `d` be a nonnegative integer, and let
`f : F_2^n -> F_2` be nonconstant on every affine subspace of dimension at
least `d`. Every deterministic acyclic, one-output circuit for `f` over
all sixteen binary Boolean truth tables, with arbitrary sharing, has at
least

    4n - 9(d+1)

gates. Input variables, wires, and output selection have zero cost; each
retained gate, including a trivial gate, has cost one. This theorem is a
class theorem; eligibility and an all-length family instantiating it must
be supplied separately.

This section adapts the complete reconstruction in public submission
`33f76827-c1af-454c-bd5b-304495c3db1a`, archive member
`submission/elimination.md`, candidate commit
`8d2c1d9d8bd13915d64b84a4853b0f951d799908`.

The three-point-pair-domain idea and the potential `s+5I+3q` are taken
from the public note titled "Coefficient 4: three-point pair
domains give 4n - 9(d+1) for every affine disperser", associated with
submission `5e9fcf7f-8338-4729-9faa-a18994661af6`. That note credits GPT 6
Astra, GPT 6.1 Sol, and Claude Opus 5.5 for the idea, accounting, and
writing. The proof below reconstructs and proves the claim independently;
neither its submission status nor its reported audits are premises.
It uses the fair affine-system/parameter framework of Mishra--Santhanam,
ECCC TR26-260 (19 October 2026), Sections 2--4. Everything needed from
that framework is proved below. Three-point fibers were also used in
weighted gate elimination by Golovnev--Kulikov, ECCC TR15-170 (ITCS 2016);
no theorem from that work is a premise here.

All arithmetic and affine maps below are over `F_2`.

### C.1. Representation and the domain invariant

We work in a relaxation of the circuit model, starting from the given
circuit without increasing its size. There are `M` independent Boolean
parameters, a collection `P` of `q` disjoint pairs of these parameters,
and a domain `D` in parameter space. For every pair `{p,h}` exactly one
corner `(alpha,beta)` is excluded; the allowed values of that pair are
the other three points of `F_2^2`. Unpaired parameters range freely.
The domain is the product of these pair domains and the free bits.

There is an affine map

    Phi : F_2^M -> F_2^n

which is injective on the **entire** parameter cube, not just on `D`.
There are auxiliary variables `v` specified by an affine system

    A v = N u + c,

where `A` is square and invertible. Each equation has at most three
distinct nonconstant terms after cancellation. An acyclic collection of
binary Boolean gates follows the system. Gate inputs can be parameters,
auxiliaries, or earlier gates. No gate is an input to an equation.
The output is a constant or a literal of a parameter, auxiliary, or
gate. Permitting the final output literal for free is a relaxation of
the stated model and therefore strengthens the lower-bound theorem.

The representation is required to compute `f(Phi(u))` for `u in D`.
Its behavior outside `D` is unrestricted. Auxiliaries nevertheless have
unique affine values on the entire cube because `A` is invertible.
The representation's size `s` is its number of equations plus retained
gates. Initially `M=n`, `Phi` is the identity, `q=0`, there are no
equations, and `s` is the original circuit size.

A value's degree is the number of distinct equations in which it occurs
as a nonzero term, plus the number of distinct gates having it as a
parent; repeated pins count once. Equations have no distinguished
output term for this definition. The output reference is not counted.
Let `I` be the number of **used unpaired
parameters**, namely unpaired parameters of positive degree. Put

    a = M-q,                 W = s + 5I + 3q.

Always `W>=0`, and initially `W<=s_initial+5n`. We retain unused
parameters in `M` unless a specified restriction removes them; syntactic
nonuse does not by itself justify a dimension decrease or increase.

For a pair excluding `(alpha,beta)`, fixing `p=1+alpha` leaves `h`
free, and fixing `h=1+beta` leaves `p` free. We call these the safe
values. Fixing one member of every pair safely gives an affine subspace
of the parameter cube, contained in `D`, of dimension `M-q=a`.
Its image under `Phi` is an affine subspace of exactly that dimension.
In particular, a constant output is impossible whenever `a>=d`.

We perform reductions only while

    a >= d+2.                                      (guard)

Two useful terminal exclusions follow directly, without requiring that
the following affine slices themselves be implemented as reductions.
On a safe affine subspace of dimension `a`, any affine output is
constant on a subspace of dimension at least `a-1`. Thus a constant or
affine output, including a literal of a parameter or auxiliary, is
impossible under the guard. An AND-type operation on two affine values
also has a constant affine slice of codimension at most one: write it
`(L+b)(R+c)+e`. If `L` is nonconstant, set `L=b`; if `L` is constant,
the operation is either constant or affine in `R`, to which the previous
argument applies. Consequently a gate of this form with affine parents
cannot be the output under the guard. Each resulting slice has dimension
at least `a-1>=d`, and `Phi` preserves its dimension.

### C.2. Boolean tables and legal simplifications

The sixteen two-input tables consist of two constants, four literals
of one input, two tables `U+V+e`, and eight tables

    (U+b)(V+c)+e.

These last eight are called AND-type. Their majority value is `e`,
and their sole exceptional input is `(1+b,1+c)`. Constant nodes introduced
by substitution are free. Negations are absorbed into descendant truth
tables or into affine equation constants; a terminal negation is covered
by the relaxed output-literal convention. No retained negation or
projection gate is silently omitted from the size count.

Here are four simplifications. Each deletes at least one gate or
equation, keeps `M,q,a` unchanged, and does not increase `I`.

**D1: a degenerate gate.** A gate which is a constant or a literal of
one actual parent is rerouted to that value and deleted. This includes
constant or unary tables, fixed pins, and repeated pins. It also applies
when both parents are affine values with dependent linear parts: they
are constant or equal up to a constant, so composing the two-input table
gives a constant or a literal of a nonconstant actual parent. A rerouted
parameter was already used in the deleted gate, so it cannot newly
enter `I`.

**D2: a dead gate.** Delete a gate of degree zero unless it supplies the
output. This cannot increase parameter usage.

**D3: a short equation.** An equation with at most two distinct terms
contains an auxiliary, since otherwise its row in `A` would be zero.
Solve it as `v=U+c` with `U` either one other value or zero. Substitute
in all other equations, gates and the output, and delete this equation
and `v`. Substitution replaces a single term by at most one term, so
equation lengths do not increase. Ordinary row elimination followed by
expansion along the `v` column proves that the remaining auxiliary
matrix is invertible. If `U` is a parameter, it was already used in the
deleted equation, so `I` does not increase.

**D4: an isolated unpaired parameter in an equation.** Suppose an
unpaired parameter `x` has just one child, an equation `E`. Some
auxiliary `Q` has coefficient one on `x` in its affine solution, because
the column `N_x` is nonzero and `A^{-1}N_x` is nonzero. Exchange the
roles of `x` and `Q`, leaving the equations' terms unchanged. Taking
`w=Q` as the new unpaired parameter, write

    Q=x+ell(other parameters),   x=w+ell(other parameters),

and make this affine substitution in `Phi`. The coordinate change is
a bijection on the full cube. It preserves the product domain: all
paired coordinates stay fixed by the change, and for every choice of
them the old free coordinate `x` is exchanged bijectively with the new
free coordinate `w`. Hence `Phi` stays affine and injective.

The new auxiliary matrix replaces the column `A_Q` by `N_x`. Its
determinant is nonzero because the `Q` coordinate of `A^{-1}N_x` is
one. The old parameter `x`, now auxiliary, occurs only in `E`; it is
not the output, since that would be affine under the guard. Delete it
and `E`. Expanding the determinant on its single-occurrence column
proves invertibility of the remaining matrix. No gate or equation has
been added. Before deletion a used unpaired parameter was exchanged
for another used unpaired parameter: `Q` occurs in some equation
because the original column `A_Q` is nonzero. Thus `I` cannot increase.

All four simplifications give `Delta W<=-1` and `Delta a=0`.
The deletion arguments count original distinct objects. In particular,
if propagation already deletes a second designated object, it counts
as deleted; if it survives, a short equation stays short under D3 and
a fixed-pin gate stays degenerate until it is deleted. These observations
will justify simultaneous deletion counts below.

An affine gate whose parents are parameters or auxiliaries can instead
be moved into the system: add its defining equation and make its output
a new auxiliary. The new matrix is block triangular with a final
diagonal entry one; the equation has at most three terms. This removes
one gate and adds one equation, so `s` does not change, and `I` cannot
increase. Repeat in topological order whenever possible. This is a
transfer of one charged object to another, not an additional free
computation.

The role exchange used in D4 is also legal without the final deletion,
whenever an unpaired `x` has coefficient one in the solution for an
auxiliary `Q`. It preserves `s,M,q,a,I`, for the same reasons. We call
this an affine role exchange. The possibly dense affine expression is
substituted only in `Phi`; the computational equations retain their
original terms. No dense affine gate is inserted for free.

### C.3. Clean up the paired parameters

At the start of each iteration, move eligible affine gates and perform
any available D1--D4 simplification. A simplification is a successful
iteration, after which we restart. Otherwise try the following pair
cases in order, restarting after any successful case.

**P0: an unused member.** If paired `p` is unused, fix its partner `h`
safely, and remove `h` and this pair. The surviving free `p` is unused.
Thus `Delta M=Delta q=-1`, `Delta a=0`, `Delta I<=0`, and
`Delta W<=-3`.

**P2: a member of degree at least two.** Fix such a member `p` safely,
remove it and the pair, and delete at least two distinct children of
`p`. A gate child has a fixed pin and is deleted by D1; an equation
child becomes short and is deleted by D3. The other member becomes an
unpaired parameter, so it is the only possible new member of `I`.
Consequently

    Delta s<=-2,  Delta I<=1,  Delta q=-1,  Delta a=0,
    Delta W<=-2+5-3=0.

When neither P0 nor P2 applies, every paired member has exactly one
child.

**PC: partners have the same child.** Suppose partners `p,h` have the
same unique child `B`. If `B` is a gate, its Boolean table has the same
value on two of the three allowed points, by the pigeonhole principle.
Those two points form an affine line contained in the allowed domain.
Restrict the pair to that line. The line is expressible by fixing one
coordinate or identifying the coordinates up to a constant; remove one
parameter and the pair. Now `B` is constant and can be deleted. If `B`
is an equation, its other term is an auxiliary and it has form
`v+p+h=c`. Among the three allowed points two have the same `p+h`;
on their line `v` is constant, and D3 deletes the equation. The
surviving parameter is unused. Therefore

    Delta s<=-1,  Delta I<=0,  Delta q=-1,  Delta a=0,
    Delta W<=-4.

The restricted `Phi` is still affine and injective. This includes the
case that the deleted gate supplies the output: if the resulting
constant output conflicts with the domain invariant, the supposed
configuration simply cannot occur under the guard.

**PE: a paired member occurs in an equation.** Let `p` have sole child
`E`, an equation, and let its partner `h` have sole child `B`. By PC,
`E` and `B` are distinct. Fix `h` safely and remove it and the pair.
Delete `B` using D1 or D3. It contained no `p`, since `p` had only
child `E`; none of these substitutions creates new occurrences of `p`.
If `E` survives, `p` is an unpaired parameter isolated in `E`, so D4
deletes `E`. If `E` has already disappeared, it is already the second
distinct deletion; `p` may remain used after rerouting. In either case
at most one parameter, `p` or its D4 replacement, can newly enter `I`. Thus

    Delta s<=-2,  Delta I<=1,  Delta q=-1,  Delta a=0,
    Delta W<=0.

If none of these pair cases applies, every paired parameter has exactly
one child, that child is a gate, and partners have distinct children.
In particular no paired parameter occurs in any equation. All auxiliary
values therefore depend only on unpaired parameters.

### C.4. Expose a first AND-type gate

Assume none of the preceding reductions applies under the guard. The
acyclic gate part is nonempty, since otherwise the output would be
affine. Its first gate `G` has only parameter or auxiliary parents.
It is neither degenerate nor affine, by D1 and the topological moves,
and hence is AND-type. Its parents' affine linear parts are independent,
or D1 would apply.

If one parent `Q` is auxiliary, its affine form depends on an unpaired
parameter `x`: it is nonconstant, and paired parameters occur in no
equation. Use an affine role exchange to make `Q` an unpaired parameter.
If the other parent is still auxiliary, linear independence guarantees
that its form involves some unpaired parameter other than the first
new parent; exchange using that parameter. Thus at most two exchanges
make both parents of `G` parameters. Paired coordinates and all their
occurrences are unchanged. Newly exposed degeneracies or eligible affine
gates can be simplified before continuing; otherwise these exchanges
have zero cost and the same `G` has two distinct parameter parents.

The terminal exclusion in Section C.1 shows that `G` is not the output.
D2 consequently gives at least one gate child of `G`. Gates cannot
have equation children. No preparation which reroutes children of `G`
away from its input parameters is needed.

Exactly three cases remain, according to the pair status of its parents.

**GP: both parents are unpaired.** Write

    G=(x+b)(y+c)+e.

Restrict `(x,y)` to the three-point majority fiber, excluding
`(1+b,1+c)`, and record this as a new pair. No parameter or coordinate
map is changed. On the new domain `G=e`, so delete `G` and each of
its original children by constant substitution and D1. There is at
least one child, and it is distinct from `G`. The previously used
unpaired parameters `x,y` cease to contribute to `I`. Ordinary
rerouting cannot create another used parameter, since any parameter
rerouted from a deleted gate was already used there. Hence

    Delta s<=-2,  Delta I<=-2,  Delta q=1,  Delta a=-1,
    Delta W<=-2-10+3=-9.

The new computation need only be correct on the new domain. In
particular `G` is not asserted constant on the full parameter cube.
The coordinate map remains affine and injective on that full cube.

**AP1: exactly one parent is paired.** Let the parents be paired `p`
and unpaired `w`, and let `h` be the partner of `p`. Let `H` be the
unique child of `h`; by the pair cleanup `H` is a gate distinct from
`G`. Fix `w` to the value forcing the AND-type `G` to its majority
constant. Also fix `h` safely in its pair. Remove `w,h` and that pair.
The surviving `p` is free and becomes unused after `G` is removed.
Delete `G` and `H`: the former is constant, and the latter has its
`h` pin fixed. These are two distinct deletions, even if `H` is a
child of `G` or they share descendants. No third deletion is needed.
Removing the used unpaired `w` loses one member of `I`; the only
newly unpaired coordinate is unused `p`. Thus

    Delta s<=-2,  Delta I<=-1,  Delta q=-1,  Delta a=-1,
    Delta W<=-2-5-3=-10<=-9.

This slightly weaker estimate than the public note's `-11` avoids any
need to count a successor of `H` and suffices for the theorem.

**AP2: both parents are paired.** They are in different pairs, since
PC excludes partners sharing `G`. Write the pairs as `{p,h}` and
`{r,k}`, with `G` reading `p,r`. Let `H,K` be the unique children
of `h,k`, respectively. They are gates, each distinct from `G`.
Choose a valid point of the first pair with `p` equal to a forcing
value for `G`, and fix that whole pair to this point. Such a point
always exists: either value of `p` occurs in the three-point domain.
The second member `r` becomes unused when the now constant `G` is
deleted; fix its partner `k` safely, so `r` remains a free unused
parameter. These restrictions remove three parameters and two pairs.

Delete `G,H,K`. If `H!=K`, these are three distinct gates. If `H=K`,
their common gate reads `h` and `k`: it has two distinct pins and both
parameters have it as their unique child. Therefore this common gate
cannot also read `G`. Since `G` has a child, choose one such child `J`;
then `G,H,J` are three distinct gates. The pins `h,k` are fixed, so
`H` is deletable, and the `G` pin of `J` is constant, so `J` is
deletable. This proves at least three distinct deletions in every
overlap pattern. The only surviving newly unpaired parameter is `r`,
which is unused; no old unpaired parameter becomes used. Hence

    Delta s<=-3,  Delta I<=0,  Delta q=-2,  Delta a=-1,
    Delta W<=-3-6=-9.

All these restrictions preserve the product form of the domain and
affine injectivity of `Phi`: a safe single-coordinate fixation removes
one parameter and one pair, a whole-pair fixation to a valid point
removes two parameters and one pair, and an unpaired fixation removes
one parameter. In particular their displayed changes in `a` are exact.

### C.5. Coverage, termination, and the inequality

Every successful reduction satisfies

    Delta a in {0,-1},        Delta W <= 9 Delta a.

The zero-dimensional-loss reductions are D1--D4 and P0, P2, PC, PE;
the unit-loss reductions are GP, AP1, AP2. Each preserves all
representation invariants, never increases `M` or `s`, and strictly
decreases the nonnegative integer `M+s`. Each preparation is finite:
affine moves transfer existing gates into equations without introducing
new gates, and at most two role exchanges expose the selected first
AND-type gate. If a simplification becomes applicable during a
preparation, perform it as the successful reduction and restart.

Coverage is exhaustive. After simplifications and pair cleanup, a
nonempty acyclic part has a first nondegenerate nonaffine gate, whose
parents can be exposed as parameters. It cannot be the output, and
the status of its two parents is exactly one of GP, AP1, AP2.
If the acyclic part were empty, the output would be affine, which is
excluded by the guard. Output constants or literals produced by a
reduction are handled by the same domain invariant; if such an output
would conflict with that invariant, the purported pre-reduction state
was impossible. No estimate presumes an extra successor of an output.

There are at most `n+s_initial` successful iterations. Consequently the
process must leave the guard, with a final `a_T<=d+1`. If `n<=d+1`
initially, the asserted lower bound is already nonpositive and is
immediate. Otherwise sum the displayed inequalities over the process:

    s_initial + 5n >= W_0
                      >= W_T + 9(n-a_T)
                      >= 9(n-d-1),

since `W_T>=0`. Therefore

    s_initial >= 4n - 9(d+1),

as claimed.

The size in this inequality is the original number of charged gates.
Every equation was introduced one-for-one from an existing charged
gate. Role exchanges change the representation and coordinate map
without adding computational objects. Every other step is a restriction
or a deletion. The argument thus covers all sixteen tables, constants,
unary and repeated-pin gates, terminal outputs, arbitrary fan-out,
shared descendants, and unrestricted nonuniform circuits in the stated
full-B2 model.

## 7. Combination, quantitative remainder, and the two score coordinates

Fix an integer n>=1024. An irreducible degree-k header exists by Section 1.
Fix any such header in an arbitrary B2 circuit for F_n. By Section 5 the
resulting m-input data function h has circuit complexity no greater than
that of F_n, and is nonconstant on every affine subspace of dimension d(n)
whenever d(n)<=m. The class theorem of Section 6, with its input count
replaced by m, gives

\[
 C_{B_2}(F_n)\ge C_{B_2}(h)\ge4m-9(d(n)+1)
                  =4n-\bigl(4k+9d(n)+9\bigr).                 \tag{J1}
\]

If d(n)>m, the asserted lower bound 4m-9(d(n)+1) is negative and follows
from nonnegativity of circuit size. Thus (J1)'s bound on F_n holds at
every integer n>=1024, including lengths with a vacuous disperser threshold.

Here is an explicit estimate of the remainder. Bertrand's bound from
Section 1 gives k<4n^{3/5}+2<=6n^{3/5}. From Section 5,
r<n/k<(1/2)n^{2/5}. Since sqrt(r)<=t<=sqrt(r)+1 and a ceiling adds at
most one, we have

\[
\begin{aligned}
d(n)&\le kr/t+kt+r(r+1)+2\\
    &\le 2k\sqrt r+k+r^2+r+2\\
    &\le 2\sqrt{kn}+k+r^2+r+2\\
    &\le \frac{21}{4}n^{4/5}+6n^{3/5}
                 +\frac12n^{2/5}+2 .                         \tag{J2}
\end{aligned}
\]

The last line uses 2sqrt(6)<5 and r^2<(1/4)n^{4/5}.
Consequently, using n^b<=n^{4/5} for 0<=b<=4/5 and n>=1,

\[
\begin{aligned}
R(n):=4k+9d(n)+9
 &\le\frac{189}{4}n^{4/5}+78n^{3/5}
                    +\frac92n^{2/5}+27\\
 &\le\frac{627}{4}n^{4/5}.                                  \tag{J3}
\end{aligned}
\]

Hence R(n)/n tends to zero over all integer lengths. For every real
epsilon>0 one may choose

\[
 N_\epsilon=\left\lceil\max\left\{1024,
                    \left(\frac{627}{4\epsilon}\right)^5\right\}\right\rceil.
\]

Every integer n>=N_epsilon then satisfies C_B2(F_n)>=(4-epsilon)n.
The family, evaluator, verifier and witness bound were fixed before
epsilon was chosen. This proves the required all-length liminf theorem.

The polynomial_log declaration has c=4, p=1, and both logarithmic
exponents zero. Thus g(n)=4n and liminf C_B2(F_n)/g(n)>=1. Its finite
projection is exactly (4,1). The first two structured manifest claims
assert the same entire theorem; their scores are its coefficient and
exponent. No weaker local lemma is used to supply a score coordinate.

## 8. Sources and provenance

- Om Swostik Mishra and Rahul Santhanam, *An Explicit 3.5n-o(n) Size Lower
  Bound against General Boolean Circuits*, ECCC TR26-260, initial version
  dated 19 October 2026. Sections 2--4, especially Theorem 4.1, supply the
  affine-system parameter framework adapted in Section 6 above. The
  coefficient-four class argument is credited separately below. PDF:
  https://eccc.weizmann.ac.il/report/2026/260/download/ . SHA-256:
  `38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc`.
- Eli Ben-Sasson and Swastik Kopparty, *Affine Dispersers from Subspace
  Polynomials*, author preprint dated 12 March 2010. Section 3, Theorem
  6.1 with Claims 6.2--6.3, and Theorem 7.3 with Lemmas 7.1--7.2 and
  Claim 7.4 supply the construction method. Sections 2--4 above give
  the needed arguments in full for characteristic two; the normalized
  counting argument in Section 4 is a direct proof of the needed
  uniqueness property. PDF:
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf . SHA-256:
  `922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf`.
- The repository's starting proof invokes Li and Yang, ECCC TR21-023,
  version dated 20 February 2021, for coefficient 31/10. That baseline
  identifies the comparison and the circuit conventions; its LP and
  unproved citations are not premises of the new proof. The source
  convention is implemented directly in Section 6. No optimality or
  stronger coefficient is inferred from the old LP.

- Our prior submission `904ae6cd-1d63-4223-a94f-5db5b6e799da`:
  Sections 1--5 of its `submission/proof.md` supply the fixed family,
  evaluator, empty-witness NP verification and full disperser proof retained
  here. Its parameter estimates are retained in Section 7. Its earlier
  coefficient 7/2 and judgment are not premises of the new circuit theorem.
- Public submission `33f76827-c1af-454c-bd5b-304495c3db1a`, candidate
  commit `8d2c1d9d8bd13915d64b84a4853b0f951d799908`, archive member
  `submission/elimination.md`: the complete coefficient-four class proof
  adapted in Section 6, including sparse affine equations, domain-preserving
  exchanges, pair cleanup, AP2 overlap handling and the potential telescope.
  The supplied proof is evaluated anew; no judgment is imported.
- Public submission `5e9fcf7f-8338-4729-9faa-a18994661af6`, candidate
  commit `4a61f024a4993592d97b49ca7d00979029b7fece`: its public note
  and the state/potential section of `submission/elimination.md` supply the
  original three-point-pair idea, safe affine core and weights s+5I+3q.
  We consulted these directly and retain this attribution in addition to
  the full reconstruction above.
- Alexander Golovnev and Alexander Kulikov, *Weighted gate elimination*,
  ECCC TR15-170 (ITCS 2016): historical antecedent for three-point fibers
  and weighted elimination. No theorem from it is a premise here.
- Public submission `08913fb8-2e14-449b-9644-36f8ffc9ac80`: its public
  note was consulted and records the PE clarification that deleting an
  equation need not make its parameter unused. Section C.3 states only
  the needed bound of at most one newly used unpaired parameter, and our
  independent final review checked this accounting directly. Its family,
  finite constant 52 and reported judgment are not used here.

The only basic general results used without a separate textbook development
are Bertrand's postulate, elementary polynomial algebra (including splitting
fields and the polynomial root bound), and elementary finite-dimensional
linear algebra. The specialized disperser and circuit arguments are proved
above. No checker, computation, copied dossier, or previous acceptance
replaces any mathematical premise.
