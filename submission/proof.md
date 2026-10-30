# A complete coefficient 7/2 lower bound for one uniform Boolean family

**Joint theorem.** The total family F_n defined in Section 1 is computed by
one fixed deterministic polynomial-time algorithm, and therefore defines one
language in NP with the empty-witness verifier specified there. In the full
Li-Yang B2 model of rules.md, at actual input length n, it satisfies

\[
\liminf_{n\to\infty}\frac{C_{B_2}(F_n)}{n}\ge\frac72.
\]

More explicitly, for every integer n>=1024, the parameters k,m,r,t and d(n)
defined below satisfy

\[
 C_{B_2}(F_n)\ge\frac72(n-k)-7d(n)
             \ge\frac72n-\frac{469}{4}n^{4/5}.
\]

The full proof is in this file. The construction and affine-disperser
argument occupy Sections 1--5, the circuit lower bound is proved in Section
6, and Section 7 combines them and verifies the all-length quantifiers.
The whole n-bit function is not claimed to be an affine disperser: invalid
field headers give zero. The disperser property is used only for a fixed
valid-header restriction, with the header length explicitly charged.

This applies published mathematical ideas; it claims no new circuit
lower-bound record. The principal sources are identified in the proof and
Section 8. No prior submission or prior judgment is used as a premise.

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


## 6. A complete gate-elimination proof

**Lemma.** If a Boolean function \(f:\mathbb F_2^n\to\mathbb F_2\)
is nonconstant on every affine subspace of dimension \(d\), where
\(1\le d\le n\), then every circuit in the challenge's full \(B_2\)
model computing \(f\) has at least

\[
                 \frac72 n-7d                                      \tag{C1}
\]

gates. All dimensions and all circuit sizes in this lemma are integers.

The argument below independently presents the parameter method of Om Swostik
Mishra and Rahul Santhanam, *An Explicit 3.5n-o(n) Size Lower Bound against
General Boolean Circuits*, ECCC TR26-260, initial version dated 19 October
2026, especially Sections 2--4, Theorem 4.1. Source PDF:
https://eccc.weizmann.ac.il/report/2026/260/download/ . Its SHA-256 is
`38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc`.
The definitions, transformations, exhaustive cases, and accounting needed for
(C1) are all supplied here; the citation is attribution, not a missing premise.

### C.1. What is counted, including boundary cases

Start with any allowed circuit, of size \(S\), including all gates whose
operations are constant or unary. Repeated pins are allowed. An operation of
two Boolean arguments either ignores both, depends on just one, or belongs to
one of these two classes:

\[
 a+b+c,\qquad (a+\alpha)(b+\beta)+\gamma,
 \quad c,\alpha,\beta,\gamma\in\mathbb F_2.                       \tag{C2}
\]

There are respectively two and eight operations in these classes. This follows
by writing the unique multilinear polynomial for a Boolean truth table as
\(c+b_1a+b_2b+b_{12}ab\): if both arguments matter and \(b_{12}=0\),
both linear coefficients are one; if \(b_{12}=1\), completing the product
gives the second form. We call these classes affine and product gates.

During the lower-bound argument only, we permit constants and complemented
wires and permit a constant or a complemented existing value as output.
This is a relaxation of the challenge model. In particular, it can only make
the minimum size smaller, so proving (C1) for this representation proves it
for the original model. No retained gate in the original circuit is declared
free. If a gate becomes a constant or a literal of one of its inputs, substitute
that value in its descendants, absorb complements into their binary truth
tables, and remove that gate. An output gate that simplifies can similarly
be removed in this relaxed representation. All such removals are charged
explicitly below. Initial output selection and inputs require no extra gates.
Thus unary terminal cases do not require an uncharged computational basis in
the original circuit. The prohibition on a simple output while the iteration
is active is proved in C.3.

Equivalently, one can start with the stronger analysis representation from
the outset: every original \(S\)-gate circuit embeds in it with exactly
\(S\) counted gates. The argument proves its numerical inequality for that
embedded circuit by valid transformations. A terminal complement need not be
translated back into a new original-model circuit, so no terminal gate can
reduce the lower bound and no additive transfer loss is required.

### C.2. States: sparse affine equations and independent parameters

A state has \(m\) independent Boolean parameters \(u=(u_1,\ldots,u_m)\)
and a map into the original \(n\)-bit input space:

\[
 \Phi(u)=b+Lu+\sum_{\{p,r\}\in P}c_{p,r}pr .                    \tag{C3}
\]

Here all coefficients are vectors or matrices over \(\mathbb F_2\), \(P\)
is a collection of disjoint pairs of parameters, and \(\Phi\) is injective.
Pairs with a zero vector coefficient are permitted. Write \(q=|P|\) and
\(a=m-q\). A parameter in a pair is paired; the others are unpaired.

The computing representation has a square affine system

\[
                         Av=Nu+c                              \tag{C4}
\]

with \(r\) auxiliary values \(v\), invertible \(A\in\mathbb F_2^{r\times r}\),
and at most three nonzero variable occurrences in each equation. It is
followed by an acyclic list of ordinary binary gates. A gate may read
parameters, auxiliary values, constants, or earlier gates. The acyclic part
never feeds back into (C4). The output is a constant or a literal of a value.
Set \(s=r+\text{number of remaining gates}\). The affine equations are
accounting objects, not newly free computations: each costs one in \(s\).
They may have cyclic syntactic dependencies; invertibility gives their unique
values for every assignment of the parameters. The central invariant is

\[
              \text{represented output at }u=f(\Phi(u)) .       \tag{C5}
\]

Initially \(\Phi\) is the identity, all inputs are unpaired parameters,
the affine system is empty, and the ordinary circuit is unchanged. Hence
\(m=a=n, q=0, s=S\).

An equation is a child of each parameter or auxiliary value occurring in it
with coefficient one. A gate is a child of each distinct value on its input
pins. Degree means the number of these distinct children, not the number of
pins; output selection is not a child. Thus coincidences of two pins never
inflate a degree. A parameter is isolated if its degree is one and unused if
its degree is zero. An unpaired parameter has weight

\[
 w(t)=\begin{cases}0&\deg(t)=0,\\6&\deg(t)=1,\\7&\deg(t)\ge2.
 \end{cases}
 \qquad W=2s+3q+\sum_{t\text{ unpaired}}w(t).                  \tag{C6}
\]

The potential is nonnegative. Paired parameters have no separate weights.
We retain unused parameters in the domain of \(\Phi\); ignoring their
occurrences does not erase dimensions from the invariant.

An affine gate whose two inputs are parameters or auxiliary values can be
replaced by the equation recording its value. For a new auxiliary \(J\),
this equation is \(J+b_1t_1+b_2t_2=c_0\), after collecting constants.
The new auxiliary matrix is block triangular with diagonal blocks \(A\) and
the scalar one, hence remains invertible. One gate becomes one equation, so
\(s\) is unchanged. Parameter degrees and weights cannot increase. Repeat in
topological order until no affine gate without a product-gate ancestor remains
outside the system. No inputs to these equations come from the remaining
acyclic part.

### C.3. Parameter changes and the disperser invariant

The following changes preserve the form (C3) and its injectivity.

* An unpaired parameter \(x\) may be replaced by \(w+\ell(u_{-x})\),
  where \(w\) is fresh and \(\ell\) is affine. This is a bijection of
  parameter spaces. It preserves \(m,q,a\).
* For three distinct unpaired parameters \(x,y,z\), replace \(z\) by
  \(w+(x+\alpha)(y+\beta)+\gamma\) and add the pair \(\{x,y\}\).
  This is also bijective. It preserves \(m\), increases \(q\) by one, and
  decreases \(a\) by one.
* An unpaired parameter \(x\) may be fixed to an affine function of the
  others. This restricts to an injectively parameterized subset and decreases
  \(m\) and \(a\) by one. The actual circuit operations below use only a
  constant or a literal of one other parameter in this restriction.
* For a pair \(\{p,r\}\), fix \(p\) to a constant or to \(r+c_0\),
  and remove \(p\) and that pair. The product \(pr\) then becomes affine
  in \(r\), using \(r^2=r\). Both \(m\) and \(q\) decrease by one, so
  \(a\) is unchanged. Either partner can play the role of \(p\).

Each new map is the old map composed with an injective substitution, so it
remains injective. Replacing an unpaired parameter changes only the affine
part of (C3), except for the one explicitly introduced disjoint pair in the
second operation. Fixing one partner destroys its pair's only quadratic
term. These observations establish the claimed preservation of (C3).

Fixing one parameter from each pair to zero makes (C3) affine in the remaining
\(m-q=a\) parameters, still injectively. Its image is therefore an affine
subspace of dimension \(a\). Consequently (C5) is nonconstant whenever
\(a\ge d\).

In fact, while \(a>d\), the output cannot be any affine function of the
parameters. A constant would already contradict the preceding paragraph.
For a nonconstant affine output, choose a parameter with nonzero coefficient.
If it is paired, first fix its partner to a constant, without changing \(a\).
Then solve for the now unpaired parameter so that the output is constant.
This last affine restriction decreases \(a\) by one, leaving \(a\ge d\),
a contradiction. Similarly the output cannot be a literal of a product gate
on two unpaired parameters: setting an appropriate one of those parameters
to a constant makes that output constant and costs just one unit of \(a\).
These are mathematical restrictions used to prove the output assertions;
they do not assert a gate-saving implementation for an arbitrary affine
restriction.

### C.4. Sound rewrites and their exact accounting

When a gate becomes a constant or a literal of an existing input, replace all
uses and delete it. This removes one gate and changes no parameter or pair.
Only an isolated unpaired input that is forwarded to several children can
gain weight, and the gain is at most one, from six to seven. Thus this
operation has \(\Delta s=-1,\Delta a=0,\Delta W\le-1\). If the
forwarded value is constant, \(\Delta W\le-2\). This applies to unary
operations, a constant pin, and repeated pins. It also applies to a product
gate with affine-valued inputs whose linear parts are dependent: over
\(\mathbb F_2\), one input is constant or the two are equal up to a
constant, so the product is constant or a literal of one input. Testing this
uses the solved expressions \(v=A^{-1}Nu+A^{-1}c\).

An unused gate other than the designated output can be deleted, with
\(\Delta s=-1,\Delta a=0,\Delta W\le-2\).

An equation with at most two nonconstant terms must contain an auxiliary
variable, since otherwise it gives a zero row of \(A\). Write it as
\(v=U+c_0\), where \(U\) is zero, another auxiliary value, or a parameter.
Replace \(v\) everywhere and delete this equation and \(v\). The number of
terms in any other equation does not increase. This preserves invertibility:
adding the eliminated equation to rows containing \(v\) is an invertible row
operation; subsequently the column of \(v\) has only its pivot, and deleting
that row and column leaves an invertible minor. There is again at most one
weight gain of one, so \(\Delta s=-1,\Delta a=0,\Delta W\le-1\).

There is also a deletion for an isolated unpaired parameter \(x\) occurring
in one equation \(E\). The column \(N_x\) is a unit column; thus
\(A^{-1}N_x\ne0\). Select an auxiliary \(Q\) whose solved expression has
coefficient one on \(x\): \(Q=x+\ell(u_{-x})\). Use the first parameter
change in C.3, taking \(Q\)'s value as the new parameter \(w\), and regard
the old \(x\) as an auxiliary instead. No equation or gate changes except
the roles of these two symbols; the three-term limit is unchanged.

To check this exchange, let \(j\) index \(Q\) and let \(t=A^{-1}N_x\),
so \(t_j=1\). Replacing column \(j\) of \(A\) by \(N_x\) gives \(AB\),
where \(B\) is the identity with column \(j\) replaced by \(t\).
Its determinant is \(t_j=1\). The new affine system is therefore invertible.
In the new parameterization the old \(x\) occurs only in \(E\), and nowhere
else, including the output by C.3. Delete this auxiliary and its single
equation; determinant expansion along its unit column again proves
invertibility of the remaining system. The old weight six is replaced by
at most seven, and one equation is deleted. Hence
\(\Delta m=\Delta q=\Delta a=0,\Delta s=-1,\Delta W\le-1\).
The map (C3) is updated by the affine exchange, so no lost computation is
concealed in this deletion.

More generally the same exchange between \(x\) and \(Q\) is possible
whenever \(x\) is unpaired and \(Q=x+\ell(u_{-x})\). It preserves
\(m,q,a,s\); its potential change is \(w'(Q)-w(x)\). We will use it
only in the preceding deletion, or when both degrees are at least two,
in which case its potential change is zero.

Constant and literal restrictions of parameters in C.3 leave \(A\)
unchanged before deletions and cannot increase an equation's term count.
They require no additional gates. Rerouting a deleted acyclic gate cannot
affect (C4). These facts together with the determinant arguments justify all
matrix manipulations below, even when auxiliary equations have cyclic
dependencies.

For simultaneous deletions of children of fixed values below, retain the
identities of the specified gates and equations until they are deleted.
Restrict parameters first. Every equation targeted for a two-term deletion
then has at most two terms; eliminating another such equation cannot increase
this number, since it substitutes a constant or one literal. Every targeted
acyclic gate has a fixed input and thus is removable;
forwarding a predecessor preserves its represented constant or literal.
Delete targets in a compatible gate order, using the stated equation
substitutions. A target removed sooner is counted once. If its surviving form
becomes even simpler, it is still removed once. Additional deletions can be
deferred. Thus shared descendants, overlapping fan-out, and coincidences do
not create extra distinct deletions or invalidate the counts.

For clarity, an equation's identity is its row, not any particular auxiliary
name on that row. Substituting a two-term equation performs row additions,
then removes only its chosen pivot row and its auxiliary column. Keep every
other row with its original identity. A targeted remaining row still has at
most two terms; the invertible remaining auxiliary matrix ensures that it has
an auxiliary variable and is still deletable. No second targeted row disappears
silently through cancellation: a zero auxiliary row would contradict that
invertibility. This is why the counts remain valid when several targets share
auxiliary variables.
The isolated-parameter deletion has its own pivot argument above and does
not require its row to have only two terms.

### C.5. Preliminary cases

We work while \(a>d\). Move eligible affine gates into (C4), then take the
first applicable case in the following list. One complete iteration will
strictly decrease \(m+s\), with

\[
          \Delta a\in\{0,-1,-2\},\qquad \Delta W\le14\Delta a .\tag{C7}
\]

**A. A deletion from C.4 is available.** Perform it. This decreases \(s\),
does not change \(a\), and lowers \(W\). If no such deletion is available,
all gates have distinct nonconstant pins, every remaining equation has exactly
three terms, nonoutput gates have a child, no unpaired isolated parameter
occurs in (C4), and product gates with affine inputs have independent linear
parts.

**B. A product gate \(G\) on parameters \(z,w\) has an isolated unpaired
input \(z\).** Fix \(w\) to the value making \(G\) constant, then delete
\(G\) and forward that constant. Parameter \(z\) becomes unused. If \(w\)
was unpaired, its weight and \(z\)'s weight were each at least six, giving
\(\Delta a=-1,\Delta W\le-6-6-2=-14\). If \(w\) was paired, remove
its pair; its surviving partner acquires weight at most seven. Then
\(\Delta a=0,\Delta W\le-3-6-2+7=-4\). In both alternatives \(m\)
decreases and \(s\) does not increase.

When A and B are both unavailable, we use this consequence repeatedly:

> No child has both an isolated unpaired parameter and a different parameter
> among its parents.                                                   (R)

For an equation this follows from A. For a gate the affine case would already
have moved into (C4), and the product case would be B. Also an isolated
parameter cannot itself be the output by C.3. Thus its unique occurrence is
indeed in a gate with a nonparameter other input.

**C. A paired parameter has degree zero or at least two.** If \(p\) has
degree zero, fix its partner \(r\), removing the pair. The surviving \(p\)
is unused, so \(\Delta a=0,\Delta W\le-3\). If \(\deg(p)\ge2\),
fix \(p\), remove its pair, and delete two distinct children of \(p\),
whether equations or gates. They are removable because \(p\) is fixed.
By (R) none had an initially isolated unpaired parameter as another parent,
so these two deletions cannot increase any previously unpaired weight.
Only \(r\)'s new weight, at most seven, must be paid for. Hence
\(\Delta a=0,\Delta W\le-3-4+7=0\). These actions decrease \(m\).
After excluding C, every paired parameter has exactly one child.

**D. The two partners \(p,r\) have the same child \(H\).** If \(H\) is
an equation, it has form \(v+p+r=c_0\). Set \(p=r+c_1\), so \(v\)
becomes constant. If \(H\) is affine, the same substitution makes its
value constant. If \(H=(p+\alpha)(r+\beta)+\gamma\), choose
\(c_1=\alpha+\beta+1\), making the two factors complementary and
their product zero. Remove \(p\) and its pair and delete \(H\), forwarding
the constant. The survivor \(r\) is unused. Thus
\(\Delta a=0,\Delta W\le-3-2=-5\), and \(m\) decreases.

**E. A paired parameter \(p\) occurs in (C4).** Its sole child is an
equation \(H_p\); its partner \(r\) has a distinct sole child \(H_r\)
by D. Fix \(r\), remove the pair, and delete \(H_r\). Property (R)
excludes an increase of an old unpaired weight. The survivor \(p\) has
weight six and still has the unique equation child \(H_p\), so the change
so far is at most \(-3-2+6=1\). Apply the isolated-parameter equation
deletion from C.4 to \(p,H_p\); it saves at least one further unit of
potential. Hence \(\Delta a=0,\Delta W\le0\), and \(m+s\)
strictly decreases. Deleting \(H_r\) cannot introduce an occurrence of
\(p\), since \(H_r\ne H_p\) and \(H_p\) was its unique child.
If \(H_r\) is an equation, its substitution contains no \(p\), so the
coefficient of \(p\) in \(H_p\) stays one and it stays zero in every other
row. The row \(H_p\) may shorten, but the isolated-parameter deletion is
valid for a row of any length; its proof used only that the parameter column
is a unit column and that the auxiliary matrix is invertible. This covers
dependencies between the two equations.

**F. A product gate \(G\) on parameters \(p,w\) has a paired input \(p\).**
Let \(r\) be its partner. The sole child of \(r\) is a gate \(H_r\ne G\)
by C--E. The other input \(w\) differs from \(p,r\); it is either an
unpaired parameter of degree at least two (B is unavailable) or belongs to a
different pair. Fix \(r\), then fix \(w\) to make \(G\) constant; delete
the two distinct gates \(H_r,G\). The survivor \(p\) is unpaired and
unused. No old unpaired weight rises, by (R) for \(H_r\) and the constant
forwarding of \(G\). If \(w\) was unpaired, remove its weight seven and
obtain \(\Delta a=-1,\Delta W\le-3-7-4=-14\). If \(w\) had another
partner, two pairs are removed, and only that other partner can acquire
positive weight, at most seven. Then
\(\Delta a=0,\Delta W\le-6-4+7=-3\). Here \(m\) decreases by two.

If no preliminary case applies, all parameters occurring in the affine system
are unpaired and have degree at least two. Furthermore, every product gate
whose inputs are parameters has two unpaired inputs of degree at least two.

### C.6. Exposing and preparing a product gate

There must be a gate outside the affine system: otherwise the output is affine
in the parameters, contrary to C.3. Take the earliest remaining gate \(G\).
It is a product gate since all earlier affine computations were moved into
(C4). Its parents are parameters or auxiliaries, with independent linear
parts by A.

Suppose an input \(Q\) is auxiliary. It depends on a parameter \(x\):
otherwise its linear part is zero. Such \(x\) occurs in an equation, since
a zero column of \(N\) is also a zero column of \(A^{-1}N\). Thus \(x\)
is unpaired of degree at least two. Exchange \(x\) and \(Q\) as in C.4.
The new parameter \(Q\) has degree at least two: an auxiliary column of
invertible \(A\) is nonzero, and \(Q\) also feeds \(G\). Both the old
and new weights are seven, so this exchange preserves \(W,m,s,q,a\).
If one input is already a parameter, choose a nonzero coefficient outside that
parameter for the exchange of the other input; independence of their linear
parts ensures that one exists. Thus at most two exchanges make both inputs
parameters, without exchanging away an exposed input.

After an exchange, again move newly eligible affine gates and inspect A--F.
If a case becomes available, perform it and end the iteration, with (C7).
Otherwise continue. Write the exposed inputs as \(x,y\). They are unpaired
and satisfy \(\deg(x),\deg(y)\ge2\).

If a child \(J\) of \(G\) also reads \(x\), its value is a Boolean
function \(h(x,G(x,y))\) of just \(x,y\). Replace \(J\)'s operation by
this binary truth table and give it inputs \(x,y\). This removes its edge
from \(G\). It does not add a gate. The only possibly increased parameter
degree is that of \(y\), whose weight is already seven. The same rewrite
applies to a child also reading \(y\). After each rewrite, move eligible
affine gates and check A--F; an applicable case ends the iteration. If none
applies, each such rewrite strictly decreases the number of children of \(G\),
so the preparatory process ends. These rewrites preserve the function for
every assignment, including coincident descendants, since each rewritten
gate computes exactly its former value.

At the end \(G\) has no child also reading \(x\) or \(y\). It is not the
output by C.3, so it has \(k\ge1\) distinct children. Reorder its inputs
so that \(h=\deg(x)\ge\deg(y)\ge2\). All preparation preserved
\(m,s,a\) and did not increase \(W\).

### C.7. The remaining four cases

**G. Both degrees are two and the second children coincide.** Besides \(G\),
the parameters \(x,y\) have one common child \(B\). If it is an equation,
write it as \(v+x+y=c_0\) and consider the function \(x+y\). If it is
a gate, consider its binary truth table. The product operation \(G(x,y)\)
has one value on three of the four assignments of \(x,y\). On at least
two of those three assignments the Boolean function associated with \(B\)
also agrees. Those two points form an affine line. Over \(\mathbb F_2^2\)
this line is \(x=c_1\), \(y=c_1\), or \(x=y+c_1\); each is a permitted
single unpaired restriction. On it \(G\) and \(B\)'s value (or \(v\))
are constant. Delete the two distinct objects and forward the constants.
One parameter is removed; the other is unused. Thus
\(\Delta a=-1,\Delta W\le-7-7-4=-18\). This accounts for the
coinciding second-child possibility before any separate-child count is used.

For the rest of the proof, if both degrees are two their other children
\(B_x,B_y\) are distinct.

**H. The ordinary constant restriction.** Fix \(x\) so that \(G\) is
constant. Delete all \(h\) children of \(x\), starting with \(G\), and
all \(k\) children of \(G\). These are \(h+k\) distinct objects because
no child of \(G\) reads \(x\). All are removable after the restriction:
equations containing \(x\) have at most two terms, its gate children have
a constant pin, and every child of \(G\) has its now constant value as pin.
The accounting in C.4 permits further simplifications without double-counting
shared descendants; it suffices to delete these specified objects.

Let \(k_1\le k\) count children of \(G\) whose other parent is an isolated
unpaired parameter. The deleted weight of \(x\) saves seven. A parameter
already of weight seven cannot gain, and an unused parameter cannot acquire
an occurrence through these rewrites. An initially isolated unpaired
parameter can gain only if its sole child is deleted and forwards it.
Property (R) excludes such a parameter from any child of \(x\); therefore
only the \(k_1\) indicated parameters can gain, each by at most one. Thus

\[
                    \Delta W\le-7-2(h+k)+k_1.                 \tag{C8}
\]

If \(h=2\), both input degrees are two. The only children of \(y\) were
\(G,B_y\). Since \(B_x\ne B_y\), and no child of \(G\) reads \(y\),
none of the other deleted objects contains \(y\). None can forward \(y\).
It loses \(G\) and gains no new child, saving at least one more potential
unit. Therefore the right side of (C8) improves to \(-8-2(h+k)+k_1\).
Using \(k_1\le k\), the cases are:

\[
\begin{array}{c|c}
\text{condition}&\text{upper bound on }\Delta W\\\hline
h\ge3&-7-2h-k\le-14\\
h=2,\ k\ge2&-12-2k+k_1\le-14\\
h=2,\ k=1,\ k_1=0&-14.
\end{array}                                                    \tag{C9}
\]

Each of these has \(\Delta a=-1\) and decreases \(m+s\).
Use H only in these circumstances. The single excluded situation is

\[
 \deg(x)=\deg(y)=2,\quad k=k_1=1.                              \tag{C10}
\]

Here \(G\)'s sole child is a gate \(D\) whose other input is an isolated
unpaired parameter \(z\); \(z,x,y\) are distinct. The second children
\(B_x,B_y\) are distinct from each other, \(G\), and \(D\). These
facts follow respectively from G's exclusion and from the rewiring in C.6.

**I. Configuration (C10), with a product gate \(D\).** Write
\(D=(G+\delta)(z+\eta)+\theta\). Choose \(x_0,y_0\) with
\(G(x_0,y_0)=\delta\), possible since a product operation takes both
values. Fix both \(x=x_0\) and \(y=y_0\). Then \(G,D\) are constant.
Delete \(G,D,B_x,B_y\). They are four distinct objects; \(B_x,B_y\)
are removable because their respective parameter input is fixed, whether
they are gates or equations. The isolated \(z\) loses its sole child \(D\)
and becomes unused. None of these rewrites increases an old unpaired weight:
\(G,D\) forward constants, and (R) excludes an isolated unpaired parent
from \(B_x,B_y\). Consequently

\[
 \Delta a=-2,\qquad \Delta W\le-7-7-6-8=-28=14\Delta a .     \tag{C11}
\]

**J. Configuration (C10), with an affine gate \(D=G+z+c_0\).** Make the
bijective parameter change \(z=w+G(x,y)+c_0\), introducing the disjoint
pair \(\{x,y\}\) in (C3). This is the quadratic parameter operation
from C.3. Since \(z\) occurs only at \(D\), no equation or other gate
has to implement this expression. Instead \(D\)'s value becomes exactly
the fresh parameter \(w\). Forward \(w\) in place of \(D\) and delete
\(D\). Now \(G\) is unused and is not the output, so delete it too.
The representation still satisfies (C5) and the affine system is unchanged.
Here \(m\) is unchanged, \(q\) rises by one, \(a\) falls by one, and
\(s\) falls by two. The old weights of \(x,y,z\) are \(7,7,6\);
the new pair costs three and the new parameter \(w\) weighs at most seven.
No other parameter gets an additional occurrence. Therefore

\[
 \Delta W\le-4-7-7-6+3+7=-14=14\Delta a.                    \tag{C12}
\]

In particular the substituted product is not an uncharged gate. Its value is
incorporated into the injective description of the original input space;
the actual computing representation loses exactly the two stated gates.
The lost unit of affine dimension is recorded by the new pair and by \(a\).

### C.8. Exhaustiveness, termination, and the lower bound

For every live state \(a>d\), either a preliminary case A--F applies,
or C.6 prepares a product gate after finitely many zero-cost changes.
Preparation uses at most two exchanges, at most the current number of gates
as affine moves, and at most one rewrite per child of that chosen gate.
If it reveals a preliminary case, that case supplies the reducing step.
Otherwise G covers the common-second-child configuration; H covers all
other possibilities except (C10); and I,J cover its two possible gate
classes. Thus the cases are exhaustive for every valid state, without a
bounded-neighborhood or tree assumption.

Each completed iteration has (C7). A preliminary deletion or J decreases
\(s\); all other cases decrease \(m\). Neither quantity ever increases.
Hence the nonnegative integer \(m+s\), initially \(n+S\), decreases
strictly at each iteration. The process must reach \(a_T\le d\) in at
most \(n+S\) iterations; otherwise the exhaustive next step would give an
infinite descending sequence of nonnegative integers. Terminal constant,
unary, or empty-circuit states create no gap: while \(a>d\) C.3 forbids
their output functions, and after \(a\le d\) no further reduction is
required.

Initially \(W_0\le2S+7n\), because \(q=0\) and each of the \(n\)
parameters has weight at most seven. Summing (C7) yields

\[
 W_T-W_0\le14(a_T-n),\qquad
 2S+7n\ge W_0\ge W_T+14(n-a_T)\ge14(n-d).
\]

The last inequality uses only \(W_T\ge0\) and \(a_T\le d\).
Rearrangement proves \(S\ge7n/2-7d\). Since the original circuit was
arbitrary, this is the required minimum-size bound. The proof applies to
nonuniform circuits of arbitrary depth and unrestricted fan-out and makes
no direct-sum assumption. Its potential is used through the displayed
telescoping inequality, never identified with a gate count.

## 7. Combination, quantitative remainder, and the two score coordinates

Fix an integer n>=1024. An irreducible degree-k header exists by Section 1.
Fix any such header in an arbitrary B2 circuit for F_n. By Section 5 the
resulting m-input data function h has circuit complexity no greater than
that of F_n, and is nonconstant on every affine subspace of dimension d(n)
whenever d(n)<=m. The circuit lemma (C1), with its n replaced by m, gives

\[
 C_{B_2}(F_n)\ge C_{B_2}(h)\ge\frac72m-7d(n)
                  =\frac72n-\left(\frac72k+7d(n)\right).       \tag{J1}
\]

If d(n)>m, the right side of the first inequality's lower bound is negative,
so the asserted lower bound follows directly from nonnegativity of circuit
size. Thus (J1)'s bound for F_n holds at every integer n>=1024, including
any length where the disperser threshold is vacuous.

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

For the last line use 2sqrt(6)<5 and r^2<(1/4)n^{4/5}.
Combining this with k<=6n^{3/5}, and using n^b<=n^{4/5} for
0<=b<=4/5 and n>=1, yields

\[
\begin{aligned}
R(n):=\frac72k+7d(n)
 &\le\frac{147}{4}n^{4/5}+63n^{3/5}
                    +\frac72n^{2/5}+14\\
 &\le\frac{469}{4}n^{4/5}.                                  \tag{J3}
\end{aligned}
\]

Hence R(n)/n tends to zero over all integer lengths. For any real
epsilon>0 one may choose

\[
 N_\epsilon=\left\lceil\max\left\{1024,
                    \left(\frac{469}{4\epsilon}\right)^5\right\}\right\rceil.
\]

Then every integer n>=N_epsilon satisfies
C_B2(F_n)>=(7/2-epsilon)n. The family, evaluator, verifier and witness
bound were fixed before epsilon was chosen. This proves the required
liminf theorem, not a sparse-length or limsup assertion.

In the polynomial_log declaration, c=7/2, p=1, and both logarithmic
exponents are zero. Thus g(n)=(7/2)n and the proved conclusion is
liminf C_B2(F_n)/g(n)>=1. Its finite projection is exactly (7/2,1).
The first two structured manifest claims assert this same entire theorem;
their score fields are its coefficient and exponent. Neither coordinate
comes from a separate lemma. There is no superlinear or endpoint claim.

## 8. Sources and provenance

- Om Swostik Mishra and Rahul Santhanam, *An Explicit 3.5n-o(n) Size Lower
  Bound against General Boolean Circuits*, ECCC TR26-260, initial version
  dated 19 October 2026. Sections 2--4, especially Theorem 4.1, supply the
  parameter representation and gate-elimination method independently
  detailed in Section 6 above. PDF:
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

The only basic general results used without a separate textbook development
are Bertrand's postulate, elementary polynomial algebra (including splitting
fields and the polynomial root bound), and elementary finite-dimensional
linear algebra. The specialized disperser and circuit arguments are proved
above. No checker, computation, copied dossier, or previous acceptance
replaces any mathematical premise.
