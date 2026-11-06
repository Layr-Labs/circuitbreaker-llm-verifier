# Full-B2 lower bound of 4n-10d for an affine disperser

This file reconstructs the size argument in Mishra--Santhanam, ECCC
TR26-260, **Revision 1**, November 3, 2026, Sections 2--4. A function on n
bits is an affine disperser for dimension d if it is nonconstant on every
d-dimensional affine subspace, and hence on every larger affine subspace.

**Theorem.** Every full-B2 circuit computing such a function has at least
4n-10d gates.

No quadratic term is ever inserted into the parameter map in this Revision-1
argument. Pairs below are bookkeeping constraints, not quadratic terms.

## 1. Circuit normalization and output convention

A B2 gate has one of the 16 binary truth tables. Of these, two are constant,
four are unary, two have the form u+v+c over F_2, and eight have the form
(u+alpha)(v+beta)+gamma. Call the last two classes XOR-type and AND-type;
both types together also include the constant and unary degeneracies in our
normalization. Negations on internal wires are absorbed into descendant truth
tables. Inputs and wires cost zero, and every original gate is counted once,
including trivial or repeated-pin gates.

Start with an arbitrary circuit, not necessarily minimum or a formula. We may
allow its output to be a constant or a literal of any internal value without
adding a gate; this only makes the lower bound stronger. Every later deletion
is charged against the initial gate count. Unrestricted sharing and depth are
allowed throughout.

## 2. States, pair constraints, and fair semi-circuits

A state consists of an injective affine map, a set of disjoint pairs, and a
fair semi-circuit. There are m independent parameters u and an injective map

    Phi(u)=b+L*u,

where b is a constant vector and L has linearly independent columns. If P is
a set of q disjoint parameter pairs, define

    Z(P)={u : p*r=0 for every {p,r} in P},    a=m-q.

Thus each pair promises that at least one member is zero. Initially Phi is the
identity, P is empty, and m=a=n.

The computational representation has an affine system and an acyclic part.
The system is

    A*v=N*u+c,

with r auxiliary variables v and square invertible A. Each row has at most
three nonconstant terms. The acyclic part is topologically ordered and may
read constants, parameters, auxiliary values, and earlier gates, but no wire
returns to the affine system. The output is a constant or literal of any
internal value. Fairness means that every parameter assignment has a unique
solution v=A^(-1)N*u+A^(-1)c. Let s be the number of equations plus the
number of acyclic gates; initially s=|C0|.

The maintained semantic invariant is

    C(u)=C0(Phi(u)) for every u in Z(P),

and every state operation below preserves it. The number of direct children
of a value is its degree: count each equation or gate once even if a pin is
repeated. An unused parameter has degree zero, and an isolated one has degree
one.

### Four parameter operations

(AR) Affine reparameterization replaces an unpaired x by
w+phi(other parameters), where w is fresh. This is bijective and preserves
m,q,a. It is invoked only when an auxiliary Q has the solved expression
Q=x+phi(other parameters); Section 4 then exchanges the roles of x and Q
without changing any equation row. No arbitrary dense substitution is made
into the computational representation.

(UF) Fix an unpaired x to a constant. This decreases m and a by one and
leaves q unchanged. Only constant UF instances are used below.

(PF) For a pair {p,r}, restrict to one of p=0, r=0, or p=r+1, remove the
restricted parameter, and remove that pair. These are exactly the three affine
lines inside pr=0. Both m and q decrease by one, so a is unchanged.

(PC) Add a new pair {x,y} of distinct unpaired parameters. This leaves Phi
and m unchanged, increases q by one, and decreases a by one. On the enlarged
pair set, xy=0.

Each operation composes Phi with an affine map tau. In AR, tau is bijective;
in UF and PF it is injective on the remaining parameters; in PC it is the
identity. Each new Z(P') maps into the old Z(P), and Phi composed with tau is
injective and affine. Constant UF and the three listed PF lines cannot
increase a row's term count, while admitted AR exchanges leave rows unchanged.
Thus all invariants are preserved by the operations defined above.

### Why simple outputs cannot occur while a>d

The following covers constants, unary outputs, and terminal AND-type outputs.
Choose one member from every pair and fix it to zero. This gives an affine
space U0 contained in Z(P) of dimension a.

If the output is affine in the parameters, either it is constant on U0, or one
nonconstant affine equation makes it zero on an affine subspace U of U0 of
dimension a-1. If the output is a literal of

    (ell1+alpha)(ell2+beta)+gamma

with affine ell1,ell2, either both factors are constant on U0, or forcing one
nonconstant factor to zero makes the output constant on a dimension-(a-1)
subspace. In every case we obtain U contained in Z(P), of dimension at least
a-1, on which the output is constant.

If a>d, this dimension is at least d. Since Phi is injective and affine,
Phi(U) is an affine subspace of the original input space of dimension at
least d on which C0 is constant. That contradicts the affine-disperser
property. Therefore, for an exact circuit, no such simple output can occur
while a>d. This argument does not overlook a constant circuit or a terminal
AND gate.

## 3. Weight and potential

An unpaired parameter has weight 0 if unused and weight 6 otherwise. Paired
parameters have no individual weight. If I is the number of used unpaired
parameters, set

    W=s+4q+6I.

Always W>=0. Initially q=0, s=|C0|, and at most n parameters are used, so

    W0<=|C0|+6n.

Rerouting means replacing every occurrence of a value that became constant or
a literal. In each deletion below, literal rerouting targets a value already
present in the deleted object. Hence a previously unused parameter cannot
become used, so no unpaired weight can rise.

## 4. Preserving the affine system

### Moving affine gates

If an affine gate has no AND-type ancestor, move it in topological order into
the affine system. A gate J=b1*t1+b2*t2+c becomes the three-term equation
J+b1*t1+b2*t2=c and a new auxiliary variable. The enlarged auxiliary matrix
is block triangular with old A and a new diagonal 1, hence invertible. One
charged gate becomes one charged equation, so m,q,a,s do not change and W
does not increase. Repeat until every affine gate without an AND ancestor has
moved.

### AR as an exchange

Suppose an auxiliary Q depends with coefficient one on an unpaired parameter
x in its solved affine expression:

    Q=x+phi(other parameters).

Use AR on Phi, rename Q as the new parameter w, and regard x as auxiliary.
In the affine system this exchanges Q's column in A with x's column in N,
leaving every equation term unchanged. Let N_x be x's column and j Q's
coordinate. The j entry of A^(-1)N_x is the coefficient of x in Q, namely 1.
The new auxiliary matrix is A times the identity with column j replaced by
A^(-1)N_x; its determinant is therefore nonzero. No gate, equation, m,q,a,s,
or pair changes. Both old x and new w occur in the system, so both have weight
6 and W is unchanged.

A general affine substitution could make a row dense; this exchange avoids
that by changing only Phi and the roles of two variables.

### Rerouting an equation

If an equation has at most two nonconstant terms, invertibility of A implies
that it contains an auxiliary v and can be written v=U+c, where U is zero or
one other value. Substitute this expression everywhere and delete v and its
row. Substitution is an elementary row operation on A; afterwards v's column
has one nonzero entry, so determinant expansion along that column proves the
smaller matrix invertible. Term counts cannot increase. The same determinant
argument applies when deleting an auxiliary v and the only equation in which
it occurs.

## 5. Four deletion rules

(D1) Delete a gate computing a constant or a literal of one parent and reroute
its uses. This covers constant and unary tables, fixed or repeated pins, and
an AND on affine parents with dependent linear parts: those values are equal
up to a constant, making the table constant or unary. Also delete a gate that
became constant on Z(P). Then Delta a=0, Delta s=-1, Delta W<=-1.

(D2) Delete a non-output gate with no children. Again Delta a=0 and
Delta W<=-1.

(D3) Delete a row with at most two nonconstant terms as in Section 4. Then
Delta a=0, Delta s=-1, Delta W<=-1.

(D4) Let unpaired isolated x have a unique child equation E. Since x's column
N_x is nonzero, A^(-1)N_x is nonzero. Choose an auxiliary Q with nonzero
coordinate in A^(-1)N_x; its solved expression depends on x with coefficient
one. Exchange x and Q by AR. The former x is now auxiliary and still occurs
only in E. Delete E and auxiliary x by the determinant argument above. The
old parameter of weight 6 is replaced by Q of weight at most 6, while m,q,a
are unchanged and s falls by one. Thus Delta W<=-1.

Apply D1--D4 and affine-gate moves until none applies. From then on, every
affine-system row has three terms, every non-output gate has a child, and
every gate with two distinct nonconstant parents is nondegenerate.

## 6. Paired-parameter reductions

Use these cases in order. Each preserves all invariants and has Delta a=0.

**Case 1: a paired p is unused or has at least two children.** Let r be its
partner. If p is unused, PF-fix r=0 and erase the pair. The surviving p is
unused, m,q each fall by one, s and a stay fixed, and Delta W<=-4. If
p has at least two children, PF-fix p=0 and delete two distinct children by
D1 or D3. If a child disappears during rerouting, it is already one of the
counted deletions. Only the newly unpaired r can gain weight, by at most 6.
Thus Delta m=-1, Delta q=-1, Delta a=0, Delta s<=-2, and

    Delta W<=-4-2+6=0.

After Case 1 no longer applies, each paired parameter is isolated.

**Case 2: partners p,r have the same unique child H.** If H is a gate, its
value on the three points (0,0),(0,1),(1,0) agrees at two points. The line
through them is one of p=0, r=0, or p+r=1. PF-restrict to that line, making
H constant, then D1-delete it. If H is an equation, D3's absence gives three
terms and invertibility gives an auxiliary v, so H is v+p+r=c. PF-impose
p=r+1, reducing it to v=c+1, and D3-delete it. In either case reroute only to
constants; the surviving partner becomes unused. Therefore Delta a=0,
Delta s=-1, Delta q=-1, and Delta W<=-4-1=-5.

After this case, partners have distinct children.

**Case 3: paired p occurs in an equation Hp.** Its partner r has a different
unique child Hr. PF-fix r=0 and D1/D3-delete Hr. Since p does not occur in
Hr, p now has only Hp, so it is an unpaired isolated parameter in an equation.
D3-delete Hp if it became short, and otherwise D4-delete it with p. The new
parameter introduced by D4 has weight at most 6. Altogether Delta m=-1,
Delta q=-1, Delta a=0, Delta s<=-2, and

    Delta W<=-4-2+6=0.

When Cases 1--3 are exhausted, every paired parameter is isolated in the
acyclic part, and its partner's unique child is a different gate.

## 7. Preparing the first gate

Assume a>d and none of the above reductions applies. If there were no acyclic
gate, the output would be affine and Section 2 would give a contradiction.
Take the first acyclic gate G in topological order. It is AND-type, because
otherwise it would have moved into the affine system. D1's absence implies
that its parents are parameters or auxiliary values with independent linear
parts.

If a parent is auxiliary Q, choose any parameter x on which its solved affine
expression depends. Preliminary Case 3 shows x is unpaired. AR-exchange Q and
x at zero potential cost. The new parameter Q has degree at least two: it
occurs in the system and is a parent of G. Repeat for the second parent if
needed. Independence supplies a parameter distinct from the first new one.
After each exchange, redo affine moves and Cases 1--3; if one applies, that
finishes the iteration. Otherwise at most two exchanges make both parents
parameters x,y. Write

    G=(x+alpha)(y+beta)+gamma.

If G or its literal is the output, Section 2 again gives a contradiction.
Otherwise G has at least one child by D2. All of its children are in the
acyclic part.

## 8. Eliminating G: the last two cases

**Case 4: both x,y are unpaired.** Complement x or y in Phi and at every use
if necessary, absorbing the negation into constants and truth tables, so that
G=xy+gamma. Complementation changes no count or degree. Add pair {x,y} by PC.
Then xy=0 throughout Z(P), so G is constantly gamma. D1-delete G, and choose
one of its children and D1-delete that too. If rerouting the child yields a
literal, its other parent was already used, so no previously unused parameter
gains weight. Both x,y lose their individual weight 6 and the new pair adds
4. Thus

    Delta m=0, Delta q=1, Delta a=-1, Delta s<=-2,
    Delta W<=-2+4-6-6=-10.

**Case 5: at least one input is paired.** By Cases 1--3, each paired input
has G as its only child, and its partner has a different unique gate child.

If exactly y is paired, let r be its partner and Hr r's unique child. PF-fix
r=0, D1-delete Hr, and UF-fix x=alpha. This makes G=gamma; D1-delete G.
Now y is unused. Only already-used literal targets occur, and x loses weight
6. Hence

    Delta m=-2, Delta q=-1, Delta a=-1, Delta s<=-2,
    Delta W<=-2-4-6=-12.

If both x,y are paired, let their partners be rx,ry and choose Hrx, the unique
child of rx. The pairs are distinct. PF-fix rx=ry=0, D1-delete Hrx, then
UF-fix x=alpha and D1-delete G=gamma. As before y becomes unused and no old
unpaired weight rises. Therefore

    Delta m=-3, Delta q=-2, Delta a=-1, Delta s<=-2,
    Delta W<=-2-8=-10.

These cases are exhaustive: after AR, each of the two inputs is either paired
or unpaired. Shared descendants are handled by rerouting every occurrence,
not by assuming tree fan-out. Any gate deleted early is counted once.

## 9. Termination and the gate count

Every preparatory move preserves m,s,a and does not increase W. Every
deletion has Delta a=0 and Delta W<=-1. Cases 1--3 have Delta a=0 and
Delta W<=0; Cases 4--5 have Delta a=-1 and Delta W<=-10. Every complete
iteration therefore satisfies

    Delta a in {0,-1},       Delta W<=10*Delta a.

It also strictly decreases m+s: deletions and Case 4 decrease s, while the
other cases decrease m; AR and affine moves preserve both. Preparations are
finite, with at most two AR exchanges per restructuring, so at most n+|C0|
complete iterations occur.

Now assume contrary to the theorem that |C0|<4n-10d. After any sequence of
reductions ending in state j, summation gives

    0<=W_j<=W0+10(a_j-n)
           <(4n-10d)+6n+10(a_j-n)=10(a_j-d).

Consequently a_j>d. The next iteration cannot terminate merely because
a<=d. But while a>d, Section 2 says a simple output would force C0 to be
constant on an affine subspace of dimension at least d, contradicting the
disperser hypothesis. All other states admit one of the exhaustive
reductions. Thus the finite decrease of m+s is impossible, and no such C0
exists. This proves |C0|>=4n-10d. If the displayed right side is negative,
the bound is of course immediate.

The proof needs no direct-sum additivity, no identification of potential with
a gate count, and no quadratic disperser assumption. The potential is only a
weighted accounting device; the final inequality converts its initial and
terminal values back to the original charged gates.

## 10. Source scope

Source: Mishra--Santhanam, ECCC TR26-260, Revision 1, Theorem 4.1 and its
proof, including Lemma 2.2 and Lemma 4.2.

https://eccc.weizmann.ac.il/report/2026/260/revision/1/download/

PDF SHA256:
a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921.

The algorithmic counterexample conclusion of the published theorem is not
needed for `proof.md`; only its size consequence is claimed here. The
original report's 3.5n-7d analysis and the old quadratic substitutions are
not premises of this proof.
