# Coefficient 4 for a certificate-prefix BK family: parameter-space gate elimination

## Main theorem

For the one fixed total family $(F_N)$ defined in Section 1, which is
computable in deterministic polynomial time and hence lies in NP with
empty witnesses, every circuit in the full Li-Yang $B_2$ model
computing $F_N$ at actual input length $N=|x|$ satisfies, for every
integer $N\ge 1$,

$$
C_{B_2}(F_N)\;\ge\;4N-R(N),\qquad
R(N)=140\,N^{4/5}+56\,N^{3/5}+50\,N^{1/5}+38.
$$

This is exactly the declared `eventual_lower_bound`: the leading term
is $L(N)=4N$ and the explicit remainder is $R(N)$, with threshold
$N=1$ (the inequality is proved for every $N\ge1$, not eventually).
Since $R(N)/N\to0$, it also gives the finite projection $(c,p)=(4,1)$:
$\liminf_N C_{B_2}(F_N)/N\ge 4$.

The proof has three parts, and no part is a cited black box.
Section 1 defines the family (a Ben-Sasson–Kopparty payload whose
field-defining polynomial is supplied and *verified* in the input,
so no irreducible-polynomial construction theorem is needed), proves
uniform membership in P, and derives the payload's disperser property
from Lemma D of `disperser.md`, where the Ben-Sasson–Kopparty
nonconstancy argument is reproduced in full. Sections 3–9 prove,
completely and self-containedly, the general gate elimination theorem
(Theorem G): every affine disperser for dimension $d$ requires $B_2$
circuit size at least $4n-10d$; this is also the unscored lemma
declared as a string claim. The parameter-pair gate elimination
argument is due to Mishra–Santhanam (Section 10 gives provenance);
it is reproduced here in full. The Li-Yang $31/10$ theorem and its LP
are *not* used, and neither are the BK Theorem 2.2 statement or any
result of Shoup: every external mathematical input has its supporting
argument reproduced in the submitted evidence.

## 1. The fixed family, uniform membership in P and NP, disperser property

**Parameters.** For an input length $N\ge1$ let $k$ be the least prime
with $k^5>32N^3$ (equivalently the least prime $>2N^{3/5}$; the integer
form is what the evaluator uses). If $k\ge N$ the evaluator outputs
$0$. Otherwise put $m=N-k$, $r=\lceil m/k\rceil$,
$t=\lceil\sqrt r\,\rceil$, and $w_i=r+2-i$ for $i=1,\dots,r$.

**The function.** Read the first $k$ input bits $b_0,\dots,b_{k-1}$ as
$Q(Z)=Z^k+\sum_{j=0}^{k-1}b_jZ^j\in\mathbf F_2[Z]$. Test whether $Q$ is
irreducible (exact test below). If the test fails, output $0$. If it
passes, let $K=\mathbf F_2[Z]/(Q)$ with polynomial basis
$1,Z,\dots,Z^{k-1}$, zero-pad the remaining $m$ input bits to $rk$
bits, split them into blocks $a_1,\dots,a_r\in K$ (first bit of each
block the constant coefficient), put $y_i=a_i^{2^{w_i}-1}$, and output
$\pi(\mathrm{Sym}_r^t(y_1,\dots,y_r))$, where $\pi$ extracts the
constant coefficient and $\mathrm{Sym}_r^t$ sums the products indexed
by $t$-element subsets. For $N=0$, $F_0=0$. This is one fixed total
function per length; $Q$ is part of the input, not advice.

**The irreducibility test is exact.** *Claim: for prime $k$, a monic
degree-$k$ polynomial $Q\in\mathbf F_2[Z]$ is irreducible iff
$Z^{2^k}\equiv Z\pmod Q$ and $\gcd(Q,Z^2+Z)=1$.*

*Proof.* In an algebraic closure, $Z^{2^k}-Z$ is separable (its
derivative is $1$) and its root set is the field $\mathbf F_{2^k}$:
the roots are closed under addition, multiplication and inversion and
number $2^k$. The first condition says $Q\mid Z^{2^k}-Z$ in
$\mathbf F_2[Z]$. Every irreducible factor of $Z^{2^k}-Z$ has degree
dividing $k$: a root $\alpha$ has $\mathbf F_2(\alpha)\subseteq
\mathbf F_{2^k}$, so $[\mathbf F_2(\alpha):\mathbf F_2]$ divides $k$.
Since $k$ is prime, factor degrees are $1$ or $k$; the degree-$1$
factors are exactly $Z$ and $Z+1$ (roots $0,1$), whose product is
$Z^2+Z$. If $Q$ is irreducible of degree $k\ge3$, its root generates
$\mathbf F_{2^k}$, so $Q\mid Z^{2^k}-Z$, and $Q\notin\{Z,Z+1\}$ gives
the gcd condition. Conversely, if both conditions hold, every
irreducible factor of $Q$ divides $Z^{2^k}-Z$ (so has degree $1$ or
$k$) and none has degree $1$ (a degree-$1$ factor would divide
$\gcd(Q,Z^2+Z)$); so every factor has degree $k$, and since
$\deg Q=k$, $Q$ is itself irreducible. $\square$

Both test steps are polynomial time: $Z^{2^k}\bmod Q$ by $k$ repeated
squarings with reduction (schoolbook $O(k^2)$ bit operations each),
the gcd by the Euclidean algorithm. Note the evaluator never
*constructs* an irreducible polynomial; it only verifies one supplied
in the input. This is what removes any need for a deterministic
construction theorem.

**A valid prefix exists.** The roots of $Z^{2^k}-Z$ in
$\mathbf F_{2^k}$ not lying in $\mathbf F_2$ number $2^k-2$; each
monic irreducible of degree $k$ has exactly $k$ roots there (one
Frobenius orbit, of size $k$ since the root is not in $\mathbf F_2$
and orbit sizes divide the prime $k$); orbits partition the set, so
the number of monic degree-$k$ irreducibles is $(2^k-2)/k\ge1$ for
$k\ge2$. Fix $Q_0$ to be the lexicographically least valid prefix;
it is used only in the lower-bound proof, never by the evaluator.

**Uniform P and NP.** One deterministic algorithm computes $F_N$ at
all lengths: find $k$ by scanning integers upward and testing
$p^5>32N^3$ and primality by trial division; by Bertrand's postulate
applied to $M=\lfloor2N^{3/5}\rfloor+1\ge3$ there is a prime in
$(M,2M)$-range giving

$$k\;\le\;4N^{3/5}+2\qquad(N\ge1), \tag{1}$$

so the scan is polynomial in $N$. Parsing, the irreducibility test,
and field arithmetic in the polynomial basis are polynomial in $k$.
Binary powering for the exponent $2^{w_i}-1$, whose bit length is
$w_i\le r+1$, takes $O(w_i)$ field multiplications per $y_i$. The
symmetric polynomial is computed without characteristic-two division:
start $E_0=1$, $E_j=0$ for $1\le j\le t$; after reading $y_i$ update
$E_j=E_j+y_iE_{j-1}$ in descending $j$ from $\min(i,t)$ to $1$
(splitting $t$-subsets by membership of $i$ proves the recurrence);
$O(rt)$ field operations. Output $\pi(E_t)$. Thus one polynomial-time
evaluator $A$ computes the total family. With $V(x,w)=1$ iff $w$ is
empty and $A(x)=1$ and witness bound $W\equiv0$, both NP directions
are immediate, so $L\in\mathrm P\subseteq\mathrm{NP}$. The witness is
not an input to the target circuits, whose input length is exactly
$N=|x|$; the zero padding is inside the payload embedding and adds no
circuit inputs; the $k$ prefix bits are genuine inputs and are
charged to $N$ in the bound.

**Restriction to the payload.** Fix $N$ and a $B_2$ circuit $C$
computing $F_N$. Substitute the constant bits of $Q_0$ for the first
$k$ inputs. No gate is added and no gate count grows: the same gate
list, with some inputs now constant, computes the payload function
$f:\mathbf F_2^m\to\mathbf F_2$,
$f(a)=\pi(\mathrm{Sym}_r^t(a_1^{2^{w_1}-1},\dots,
a_r^{2^{w_r}-1}))$ over $K=\mathbf F_2[Z]/(Q_0)$, where $a$ is read
through the zero-padding embedding $\iota:\mathbf F_2^m\to K^r$.
Hence $C_{B_2}(f)\le C_{B_2}(F_N)$: input fixing is free in the $B_2$
cost model (inputs cost $0$; retained gates are counted once).

**Disperser property of the payload.** Call $N$ *good* if $k<N$,
$t<r$, and $kt\ge r(r+1)$. For good $N$, Lemma D of `disperser.md`
applies verbatim to $f$ (with $\pi$ the constant-coefficient
functional, nonzero since $\pi(1)=1$): $f$ is nonconstant on every
affine subspace of $K^r$ of dimension $>D_0(N)$, where

$$D_0(N)\;=\;\frac{kr}{t}+kt+r(r+1).$$

Since $\iota$ is $\mathbf F_2$-linear and injective, affine subspaces
of $\mathbf F_2^m$ map to affine subspaces of $K^r$ of the same
dimension, so $f$ is an affine disperser for dimension
$d(N)=\lfloor D_0(N)\rfloor+1$ on $\mathbf F_2^m$.

**Remainder accounting.** For all $N\ge1$:

$$r\;\le\;\frac{m}{k}+1\;\le\;\frac12 N^{2/5}+1, \tag{2}$$

using $k>2N^{3/5}$. Hence, with $t\ge\sqrt r$ and
$\sqrt r\le\sqrt{\tfrac32}\,N^{1/5}\le1.25\,N^{1/5}$,

$$\frac{kr}{t}\le k\sqrt r\le 5N^{4/5}+2.5N^{1/5},\qquad
kt\le 5N^{4/5}+4N^{3/5}+2.5N^{1/5}+2,$$
$$r(r+1)\le \tfrac{15}{4}N^{4/5}\le4N^{4/5},$$

(the second using $t\le\sqrt r+1$), so
$D_0(N)\le14N^{4/5}+4N^{3/5}+5N^{1/5}+2$ and, with (1),

$$4k+10\,d(N)\;\le\;4k+10D_0(N)+10\;\le\;R(N)
=140N^{4/5}+56N^{3/5}+50N^{1/5}+38. \tag{3}$$

**Bad lengths are few and trivial.** If $N$ is not good, then
$N<10^6$: (a) if $k\ge N$, then by (1)
$N^{2/5}\le4+2/N^{3/5}\le6$, so $N\le6^{5/2}<89$; (b) if $t\ge r$ then
$r\le2$ (for $r\ge3$, $\lceil\sqrt r\,\rceil\le\sqrt r+1<r$), so
$N\le3k\le12N^{3/5}+6$ and $N^{2/5}\le18$, $N<1375$; (c) if
$kt<r(r+1)$ with $r\ge3$, then $t\ge\sqrt r$ gives
$k<r^{3/2}+r^{1/2}\le2r^{3/2}$, so
$N=m+k\le kr+k\le2r^{5/2}+2r^{3/2}\le4r^{5/2}$; if moreover
$N\ge10^6$ then $N^{2/5}\ge251$ and (2) gives
$r\le N^{2/5}(\tfrac12+\tfrac1{251})\le0.504\,N^{2/5}$, whence
$N\le4(0.504)^{5/2}N<0.73\,N$, a contradiction. For every
$N\le10^6<35^5$ we have $140N^{4/5}\ge4N$ (equivalently
$N^{1/5}\le35$), so $R(N)\ge4N$ and the Main Theorem's inequality is
trivial at bad lengths since $C_{B_2}(F_N)\ge0\ge4N-R(N)$.

## 2. The $B_2$ model and gate types

Circuits are deterministic, acyclic, single-output, with two ordered input
pins per gate and all 16 binary Boolean truth tables allowed; depth and
fan-out are unrestricted; the circuit must compute $f_n$ exactly on every
input. Inputs, wires and output selection cost zero; each retained gate
costs one, including constant and unary gates. Negations are absorbed
into the truth tables of descendant gates; a constant node used during a
substitution costs zero. These are the conventions of rules.md, following
Li-Yang Sections 1.1, 2.1, 2.4 (footnotes 3-4) and 3.3. Consequently,
"rerouting" a gate $J$ that computes a constant $c$ (replacing every use
of $J$ by the constant node $c$) or a literal of one of its inputs
(replacing every use of $J$ by that input, with the negation, if any,
absorbed into each reading gate's truth table) removes $J$ and adds no
gate. This is the only cost convention the elimination steps use.

Identify $\{0,1\}$ with $\mathbf F_2$. Every binary operation has a unique
multilinear form $h(t_1,t_2)=c+b_1t_1+b_2t_2+b_{12}t_1t_2$. Of the 16
operations, 2 are constant and 4 are unary (depending on at most one
input); of the remaining 10, the 2 *XOR-type* operations are
$t_1+t_2+c$ and the 8 *AND-type* operations are
$(t_1+\alpha)(t_2+\beta)+\gamma$. XOR-type and more generally affine
operations ($b_{12}=0$) are called *affine*. A gate is constant, unary,
XOR-type or AND-type according to its operation.

## 3. The general theorem and its application

**Theorem G (gate elimination for affine dispersers).** *Let
$f:\mathbf F_2^n\to\mathbf F_2$ be an affine disperser for dimension
$d\ge 1$, i.e. $f$ is nonconstant on every $d$-dimensional affine
subspace. Then every $B_2$ circuit computing $f$ has size at least
$4n-10d$.*

Sections 4-9 prove Theorem G. Application: fix $N\ge1$. If $N$ is
not good (Section 1), the bound is trivial by the bad-length analysis
there. If $N$ is good, Section 1 gives a payload function
$f:\mathbf F_2^m\to\mathbf F_2$, $m=N-k$, which is an affine disperser
for dimension $d(N)$, and $C_{B_2}(F_N)\ge C_{B_2}(f)$. If
$d(N)\ge m$, then $4m-10d(N)<0$ and there is nothing to prove;
otherwise Theorem G gives $C_{B_2}(f)\ge4m-10d(N)$. In both cases,
by (3) of Section 1,

$$
C_{B_2}(F_N)\;\ge\;4m-10d(N)\;=\;4N-\bigl(4k+10d(N)\bigr)
\;\ge\;4N-R(N).
$$

No integer rounding is used: the inequality is proved for the
real-valued right-hand side, and the declared remainder is exactly
$R(N)=140N^{4/5}+56N^{3/5}+50N^{1/5}+38$.


## 4. States: parameters, the map $\Phi$, pairs, and the restricted domain

Fix a $B_2$ circuit $C_0$ computing the disperser $f$ of Theorem G.
Gate elimination maintains a *state* $(D,\Phi,P)$:

- $D$ is the current fair semi-circuit (Section 5) on parameters
  $\vec u=(u_1,\dots,u_m)$;
- $\Phi:\mathbf F_2^m\to\mathbf F_2^n$ is an injective affine map,
  $\Phi(\vec u)=\vec b+L\vec u$;
- $P$ is a set of pairwise disjoint *pairs* of parameters; a parameter in
  a pair is *paired*, the other member is its *partner*; all others are
  *unpaired*.

The *restricted domain* is

$$
Z(P)=\{\vec u\in\mathbf F_2^m : pr=0\ \text{for every }\{p,r\}\in P\},
$$

i.e. for each pair, at least one member is $0$. The invariant is

$$
D(\vec u)=C_0(\Phi(\vec u))\quad\text{for all }\vec u\in Z(P). \tag{$\ast$}
$$

Initially $m=n$, $\Phi$ is the identity, $P=\varnothing$ (so
$Z(P)=\mathbf F_2^n$) and $D=C_0$. Write $q=|P|$ and $a=m-q$. Note $a$ is
the dimension of the largest affine subspace of the form used in
Section 6 that $Z(P)$ is guaranteed to contain: fixing one member of
each pair to $0$ leaves $m-q$ free coordinates.

Four kinds of change are used. Each replaces $\Phi$ by $\Phi\circ\tau$,
where $\tau$ expresses old parameters in terms of new ones; $\tau$ is
injective and affine, so $\Phi\circ\tau$ is again injective affine of
the form $\vec b+L\vec u$, and each operation is arranged so that
$\tau$ maps the new restricted domain into the old one, preserving
($\ast$) for the transformed semi-circuit $D'$ defined by
$D'(\vec u)=D(\tau(\vec u))$ on the new domain.

- **(AR) Affine reparameterisation.** For an unpaired parameter $x$ and
  an affine map $\varphi$, substitute $x\leftarrow w+\varphi(\vec u_{-x})$
  with $w$ fresh. $\tau$ is a bijection preserving paired coordinates,
  so it maps $Z(P)$ onto $Z(P)$; $\operatorname{im}\Phi$ is unchanged.
  $\Delta m=\Delta q=\Delta a=0$.
- **(UF) Fixing an unpaired parameter.** For unpaired $x$ and affine
  $\psi$, impose $x=\psi(\vec u_{-x})$ (including a constant). $\tau$ is
  injective, paired coordinates are unchanged, and the image lies in
  the old domain. $\Delta m=-1$, $\Delta q=0$, $\Delta a=-1$.
- **(PF) Fixing a paired parameter.** For $\{p,r\}\in P$, impose one of
  the three lines contained in $pr=0$: $p=0$, $r=0$, or $p=r+1$; remove
  the substituted parameter and the pair. The remaining member becomes
  unpaired. By the choice of line, $\tau$ maps every new assignment
  into $Z(P)$. $\Delta m=-1$, $\Delta q=-1$, $\Delta a=0$.
- **(PC) Pair creation.** For distinct unpaired $x,y$, add $\{x,y\}$ to
  $P$. $\Phi$ and $\tau=\mathrm{id}$ are unchanged, but
  $Z(P)$ shrinks to $Z(P)\cap\{xy=0\}$. $\Delta m=0$, $\Delta q=+1$,
  $\Delta a=-1$.

A restriction only ever shrinks (or, for (AR), preserves) the set of
parameter assignments on which ($\ast$) is claimed. This is legitimate
because the conclusion drawn at the end (Section 9) is the existence of
*one* affine subspace of the final restricted domain on which $C_0$ is
constant; shrinking the domain during the argument cannot create such a
subspace, it can only discard assignments.

## 5. Fair semi-circuits, moving affine gates, and the invariants

A *fair semi-circuit* on $\vec u$ consists of an affine system followed by
an acyclic part. The affine system is $A\vec v=N\vec u+\vec c$ in
auxiliary variables $\vec v=(v_1,\dots,v_r)$, where $A$ is $r\times r$
*invertible* (fairness: every parameter assignment extends to a *unique*
solution $\vec v$), and every row of $(A\mid N)$ has at most three
nonconstant terms. The acyclic part is an ordinary circuit whose gates,
in topological order, read constants, parameters, auxiliary variables or
earlier gates; there are no wires back to the affine system. The output
is a constant or a literal of a parameter, an auxiliary variable, or a
gate. Parameters, auxiliary variables and gates are *internal values*;
the *parents* of a gate or equation are the internal values it reads or
mentions, and it is a *child* of each of them; the *degree*
$\deg(t)$ of an internal value is its number of distinct children. A
parameter is *unused* if its degree is $0$ and *isolated* if its degree
is $1$. The *size* $s$ of a semi-circuit is the number of equations plus
the number of acyclic gates. Initially the system is empty and
$s=|C_0|$: every original gate is charged exactly once, either as an
acyclic gate or, after the move below, as one equation.

Every auxiliary variable is an affine function of the parameters,
$\vec v=A^{-1}N\vec u+A^{-1}\vec c$.

*Moving an affine gate into the system.* Let $J$ be an acyclic gate with
affine operation $h(t_1,t_2)=b_1t_1+b_2t_2+c$ whose inputs are constants,
parameters or auxiliary variables. Replace $J$ by the equation
$J+b_1t_1+b_2t_2=c$ and regard $J$ as a new auxiliary variable. The new
row has at most three nonconstant terms; the new auxiliary matrix is
$\bigl(\begin{smallmatrix}A&0\\ *&1\end{smallmatrix}\bigr)$, hence
invertible; the unique solution assigns $J$ its old gate value, so the
computed function is unchanged, and $s$ is unchanged (one gate becomes
one equation). Applying this in topological order moves every affine
gate with *no AND-type ancestor* into the system; an affine gate with an
AND-type ancestor stays in the acyclic part. In particular the *first*
gate of the acyclic part in topological order, whose inputs are
necessarily constants, parameters or auxiliary variables, is never
affine after this process: if it were affine it would have been moved.

*Preservation of the two invariants* (at most three nonconstant terms
per row; $A$ invertible). Fixing a parameter to a constant deletes its
column in $N$ and updates $\vec c$, leaving $A$ unchanged. Setting
$x=y+c_0$ for distinct parameters replaces one parameter term by
another in every row, leaving $A$ and all term counts unchanged. A
general affine substitution is *never* expanded inside the system;
it is performed only in $\Phi$ by the exchange (AR) of Section 7, which
permutes the roles of total variables and therefore leaves the set of
equations, as polynomials in all variables, literally unchanged.
Rerouting an auxiliary variable $v$ given by an equation
$v=U+c_0$ (with $U$ a constant, a parameter or an auxiliary variable)
substitutes $U+c_0$ for $v$ everywhere: in other equations this is an
elementary row operation on $(A\mid N\mid\vec c)$, preserving
$\det A\ne 0$, and it cannot increase any row's term count since $v$
(one term) is replaced by at most one term. Afterwards the column of
$v$ in $A$ has a single nonzero entry, in the row of its own equation;
expanding $\det A$ along that column shows the matrix obtained by
deleting that row and column is still invertible. The same determinant
argument justifies deleting an auxiliary variable whose only occurrence
is a single equation, together with that equation: every other equation
is unchanged, the remaining auxiliary variables keep their old values
(the restricted old solution satisfies the remaining system, and
invertibility gives uniqueness), and no other internal value used that
variable.

## 6. The simple-output lemma

**Lemma S.** *Suppose $a>d$ and the output of $D$ is either (i) an
affine function of the parameters, or (ii) a literal of an AND-type gate
whose two parents are affine functions of the parameters. Then there is
an affine subspace $U\subseteq Z(P)$ of dimension at least $a-1\ge d$
and a bit $b$ with $D(\vec u)=b$ for all $\vec u\in U$.*

*Proof.* Fix one member of each pair in $P$ to $0$; the remaining
$m-q=a$ coordinates are free, giving an affine subspace
$U_0\subseteq Z(P)$ of dimension $a$, represented by the point $\vec 0$
and the standard basis vectors of the unfixed coordinates. An affine
function $\ell(\vec u)=\sum_j\lambda_j u_j+\kappa$ is constant on $U_0$
iff $\lambda_j=0$ for every free coordinate $j$; if it is not constant,
some free coordinate $j_0$ has $\lambda_{j_0}=1$, and
$\{\vec u\in U_0:\ell(\vec u)=0\}$ is an affine subspace of dimension
$a-1$ (a point and explicit basis are obtained from the coefficients,
as $\ell(\vec u)=0$ expresses $u_{j_0}$ affinely in the other free
coordinates). In case (i) apply this to $\ell$ itself: $U=U_0$ if
$\ell$ is constant there, else the level set, with $b$ the constant
value. In case (ii) the output is
$(\ell_1+\alpha)(\ell_2+\beta)+\gamma$ with $\ell_1,\ell_2$ affine in
the parameters. If both $\ell_1+\alpha$ and $\ell_2+\beta$ are constant
on $U_0$, so is the output and $U=U_0$. Otherwise one of them, say
$\ell_1+\alpha$, is nonconstant on $U_0$; on its level set of dimension
$a-1$ the output equals the constant $\gamma$. Since $a>d$ and $a,d$
are integers, $\dim U\ge a-1\ge d$. $\square$

By ($\ast$), $S=\Phi(U)$ is then an affine subspace of $\mathbf F_2^n$
of dimension at least $d$ (injectivity of $\Phi$ preserves dimension)
with $C_0(\vec x)=b$ for all $\vec x\in S$.


## 7. Potential, deletions, and affine reparameterisation

The weight of an unpaired parameter is
$\mathrm{wt}(x)=6$ if $\deg(x)\ge 1$ and $\mathrm{wt}(x)=0$ if
$\deg(x)=0$. Paired parameters have no individual weight. The potential is

$$
W \;=\; s \;+\; 4q \;+\; \sum_{x\ \text{unpaired}} \mathrm{wt}(x)\;\ge\;0.
$$

For a quantity $X$, $\Delta X$ is its change over a step (before to
after). Two different facts must not be conflated. *Weights* are
binary (a parameter is used or it is not): whenever a deletion below
reroutes a value to a literal of a parameter $x$, that parameter
already occurs in the deleted gate or equation, so it already has that
gate or equation as a child; hence a parameter that is unused before a
deletion remains unused after it, and no unpaired parameter's weight
can increase except where a case explicitly says so. Rerouting to a
constant never increases any weight. *Degrees*, by contrast, are not
preserved by rerouting (a value rerouted to can acquire all of the
deleted item's children), and no step below relies on a degree bound
surviving a deletion: each case is entered only after the preliminary
phase has re-established that no deletion applies, and within a case,
every degree-sensitive fact used after a deletion is re-derived at
that point in the text. Finally, the output selection counts as a use
of the value it references: rerouting a value also reroutes the output
if the output referenced it, and ``child'' and ``degree'' never count
the output reference itself.

**Deletions.** Each deletion has $\Delta m=\Delta q=\Delta a=0$ and
$\Delta W\le -1$.

- **(D1)** A gate that computes, on all assignments in $Z(P)$, a constant
  or a literal of one of its inputs: reroute it to that constant or
  literal and delete it; $\Delta s=-1$. This covers every constant or
  unary gate, every gate with a constant input, every gate whose two
  inputs are identical, and every gate that has become constant on
  $Z(P)$. It also covers an AND-type gate whose two inputs have
  *dependent linear parts* as affine functions of the parameters:
  writing that gate as $(\ell+\alpha)(\ell'+\beta)+\gamma$, if
  $\ell'=\ell+c_0$ then, using $t^2=t$ on $\mathbf F_2$, the product is
  an affine function of the single value $\ell$, namely a constant or
  the literal $\ell$ or $\ell+1$ of the first input; if $\ell'$ has
  linear part $0$ it is a constant function and the gate is likewise a
  constant or a literal of the first input; and symmetrically for
  $\ell$.
- **(D2)** A gate with no children that is not the output: delete it;
  $\Delta s=-1$.
- **(D3)** An equation with at most two nonconstant terms: it must
  contain an auxiliary variable, for otherwise its row of $A$ would be
  zero, contradicting invertibility. Write it as $v=U+c$ with $v$
  auxiliary and $U$ a parameter, an auxiliary variable, or $0$;
  substitute $v=U+c$ everywhere and delete $v$ and the equation (the
  determinant argument of Section 5 preserves invertibility);
  $\Delta s=-1$.
- **(D4)** An unpaired isolated parameter $x$ whose unique child is an
  equation $E$: the column $N_x$ of $N$ has a single nonzero entry, in
  the row of $E$, so $A^{-1}N_x\ne\vec 0$. Choose an auxiliary variable
  $Q$ whose coordinate in $A^{-1}N_x$ is $1$; then
  $Q=x+\varphi(\vec u_{-x})$ for an affine $\varphi$. Apply (AR) below
  to exchange $x$ and $Q$. The equations as polynomials are unchanged,
  so the new auxiliary variable $x$ still occurs only in $E$; delete
  $x$ together with $E$ (Section 5's determinant argument). The old
  parameter $x$, of weight $6$, is replaced by the new parameter $Q$,
  of weight at most $6$; $\Delta s=-1$ and $\Delta W\le -1$.

**Affine reparameterisation (AR).** Let $Q$ be an auxiliary variable
depending on an unpaired parameter $x$, i.e.
$Q=x+\varphi(\vec u_{-x})$. Substitute $x\leftarrow w+\varphi(\vec u_{-x})$
in $\Phi$, regard $Q$ as the new parameter $w$, and regard $x$ as an
auxiliary variable. Concretely this exchanges the column of $x$ in $N$
with the column of $Q$ in $A$: the new auxiliary matrix is $A'=AB$,
where $B$ is the identity with its $j$-th column (the index of $Q$)
replaced by $\vec t=A^{-1}N_x$; the coefficient of $x$ in $Q$ is
$t_j=1$, so $\det B=t_j=1$ and $A'$ is invertible. No equation or gate
changes; $\Delta m=\Delta q=\Delta a=\Delta s=0$. Moreover $x$ occurs
in some equation (otherwise $N_x=0$, hence $A^{-1}N_x=0$, contradicting
that $Q$ depends on $x$), so $\mathrm{wt}(x)=6$; and $w$ occurs in some
equation (the column of $Q$ in the invertible matrix $A$ is nonzero),
so $\mathrm{wt}'(w)=6$. Hence $\Delta W=0$.

Moving an affine gate into the system (Section 5) has
$\Delta m=\Delta q=\Delta a=\Delta s=0$ and does not increase any
parameter's degree (a gate child becomes an equation child of the same
parents), so $\Delta W\le 0$.

## 8. The iteration lemma

**Lemma I.** *Suppose $a>d$. Then a finite sequence of the operations
above either reaches a state whose output has one of the two forms of
Lemma S (and Section 9 draws the conclusion), or transforms the state
into a new state with $m'+s'<m+s$, $\Delta a\in\{0,-1\}$ and
$\Delta W\le 10\,\Delta a$.*

*Proof.* **Preliminary phase.** First check the output directly: if it
is a constant or a literal of a parameter or an auxiliary variable
(auxiliary variables are affine functions of the parameters), or a
literal of an AND-type gate whose two parents are constants, parameters
or auxiliary variables, then Lemma S applies and we stop. Otherwise
apply (D1)-(D4) and move affine gates into the system until none
applies. This phase is finite: use the measure
$s+\#\{\text{acyclic gates}\}$. A gate deletion lowers $s$ by $1$ and
does not increase the number of acyclic gates; an equation deletion
lowers $s$ by $1$; moving an affine gate into the system keeps $s$
fixed and lowers the number of acyclic gates by $1$. No step of the
phase increases the measure, each step strictly decreases it, and it
is a nonnegative integer, even though a deletion can make a previously
unmovable affine gate movable (that only gives the phase another
measure-decreasing move to perform). Then perform the first
applicable of Cases 1-3; for each case, all earlier cases are assumed
inapplicable.

**Case 1: a paired parameter has degree $0$ or at least $2$.** Let
$\{p,r\}\in P$. If $\deg(p)=0$, use (PF) to fix $r=0$; $p$ becomes an
unpaired parameter of degree $0$, hence weight $0$. Then $\Delta m=-1$,
$\Delta s=0$, $\Delta q=-1$, $\Delta a=0$, and
$\Delta W\le -4$ (only the $4q$ term changes; no weight increases, by
the rerouting remark of Section 7). If $\deg(p)\ge 2$, use (PF) to fix
$p=0$ and delete two distinct children of $p$: a child gate now has the
constant input $0$, so (D1) deletes it; a child equation loses the term
$p$ and has at most two nonconstant terms, so (D3) deletes it. The two
deletions are sequential and do not interfere: the substitution $p=0$
is global, the first deletion never reintroduces $p$ (its rerouting
targets are constants or literals of other values), and the second
child therefore still has its constant $p$-input, respectively still
lacks the term $p$, when its turn comes. Then
$\Delta m=-1$, $\Delta s=-2$, $\Delta q=-1$, $\Delta a=0$. The newly
unpaired $r$ can contribute weight at most $6$; no other unpaired
parameter increases its weight (deleted gates and equations are
rerouted to constants or to literals of values that already had the
deleted item as a child). Hence
$\Delta W\le -4-2+6=0=10\,\Delta a$.

If Case 1 does not apply, every paired parameter is isolated.

**Case 2: two partners have the same child.** Let $\{p,r\}\in P$ with a
common child $H$, the unique child of each. If $H$ is a gate, its two
inputs are exactly $p$ and $r$; its operation takes the same value at
two of the three points $(0,0),(0,1),(1,0)$ of $pr=0$, and the line
through those two points is one of $p=0$, $r=0$, $p+r=1$. Restrict by
(PF) to that line; $H$ becomes constant on the new restricted domain
and (D1) deletes it. If $H$ is an equation, it has exactly three
nonconstant terms ((D3) is inapplicable) and one of them is an
auxiliary variable (otherwise its row of $A$ would be zero), so
$H$ is $v+p+r=c$; impose $p=r+1$ by (PF), reducing $H$ to $v=c+1$,
and delete it by (D3). In both subcases the surviving member of the
pair loses its only child and becomes unused; rerouting is to
constants only. Hence $\Delta m=-1$, $\Delta s=-1$, $\Delta q=-1$,
$\Delta a=0$ and $\Delta W\le -4-1=-5$.

If Case 2 also does not apply, every paired parameter is isolated and
has a different child from its partner.

**Case 3: a paired parameter occurs in the affine system.** Let
$p\in\{p,r\}\in P$ occur in the system; its unique child $H_p$ is an
equation, and its partner $r$ has a unique child $H_r\ne H_p$ (Case 2
inapplicable), so $p$ does not occur in $H_r$. Use (PF) to fix $r=0$
and delete $H_r$: by (D1) if $H_r$ is a gate (it now has a constant
input) and by (D3) if it is an equation (it lost a term). The deletion
of $H_r$ creates no new child of $p$, since $p$ does not occur in
$H_r$ and any value rerouted from $H_r$ already had $H_r$ as a child,
whereas $p$'s only child is $H_p$. Now $p$ is an unpaired isolated
parameter with unique child the equation $H_p$. If $H_p$ has at most
two nonconstant terms, delete it by (D3); then $p$ becomes unused. If
$H_p$ has exactly three nonconstant terms, apply (D4) to $p$ and
$H_p$; (D4) does not increase the weight of any initially unpaired
parameter, and its replacement parameter has weight at most $6$.
Overall $\Delta m=-1$, $\Delta s=-2$, $\Delta q=-1$, $\Delta a=0$. The
weight ledger, with the emergent pair made explicit, is:
$-4$ from $\Delta q=-1$; $-1$ for deleting $H_r$ and $-1$ for deleting
$H_p$ in $s$; $+6$ because $p$, formerly paired with no individual
weight, becomes an unpaired parameter whose unique child is $H_p$;
then $-6$ because $p$ ends unused (via (D3), giving the branch total
$-6$) or is replaced by the
(D4) parameter of weight at most $6$ in place of its $6$ (via (D4),
net $\le 0$). Hence $\Delta W\le -4-2+6=0$, the worse of the two
branch totals, as before.

If none of Cases 1-3 applies, every paired parameter is isolated,
occurs only in the acyclic part, and has a different child from its
partner. The minimal-AND analysis of Section 9 then applies.

In Cases 1-3 the required inequality $\Delta W\le 10\,\Delta a=0$
holds, $m$ strictly decreases, and $s$ does not increase.


## 9. The minimal AND gate, Cases 4-5, and the endgame

Assume none of Cases 1-3 applies. The acyclic part has at least one
gate: otherwise the output would be a constant or a literal of a
parameter or auxiliary variable, and the preliminary phase would have
stopped via Lemma S. Let $G$ be its first gate in topological order.
By Section 5, $G$ is an AND-type gate: an affine first gate would have
been moved into the system, since its inputs are constants, parameters
or auxiliary variables. Since (D1) does not apply, both inputs of $G$
are parameters or auxiliary variables and, by the dependent-linear-parts
clause of (D1), their linear parts as affine functions of the parameters
are independent.

*Making both inputs parameters.* Suppose one input of $G$ is an
auxiliary variable $Q$. It is a nonconstant affine function of the
parameters (independence), so it depends on some parameter $x$. That
$x$ occurs in some equation: if it did not, its column in $N$, and
hence in $A^{-1}N$, would be zero, contradicting the dependence of
$Q$ on $x$. Since Case 3 does not apply, no paired parameter occurs in
the affine system, so $x$, occurring in an equation, is unpaired.
Apply (AR) to exchange $Q$ and $x$. The new parameter $w$ (formerly
$Q$) has degree at least $2$: it occurs in an equation (its column in
the invertible matrix $A$ is nonzero) and it is a parent of $G$. By Section 7, $\Delta W=0$ for this exchange (the degree $\ge2$ of
$w$ is recorded here because the weight computation for (AR) uses
that $w$ is used, and because it shows $w$ is not isolated, so the
restructuring cannot strand the new input of $G$). After the exchange,
re-run the preliminary phase and re-inspect Cases 1-3; if a deletion
or a case applies, perform it and the iteration ends, with the
accounting of Section 8 (deletions only add $\Delta a=0$,
$\Delta W\le-1$ steps, so they cannot hurt the ledger). If none
applies, repeat for the second input if it is also an auxiliary
variable: the first input is now literally the parameter $w$, a
single coordinate, so if the second input depended only on $w$, its
linear part would be a scalar multiple of $w$'s coordinate vector,
i.e. of the first input's linear part, contradicting the independence
of the two inputs, which the bijective reparameterisation preserves;
hence the second input depends on some parameter other than $w$, and
the same exchange and re-inspection follow. The persistence needed
for Cases 4-5 is now justified rather than assumed: an (AR) exchange
involves only the unpaired parameter $x$ and the auxiliary variable
$Q$; it does not touch any paired parameter, any gate, or the
polynomial of any equation — it only re-labels which of $x,Q$ is
called a parameter, with $w$ occurring in equations exactly where the
auxiliary $Q$ occurred and the new auxiliary $x$ occurring exactly
where the parameter $x$ occurred (the columns of $A$ and $N$ are
exchanged, Section 7). Consequently every paired parameter keeps the
same children it had before, remains outside the affine system (its
occurrences in equations are unchanged, since neither $x$ nor $Q$ is
paired), and keeps a child distinct from its partner's; and term
counts of equations are unchanged, so the (D3)-inapplicability facts
used in Section 8 are preserved. The standing consequences of Cases
1-3 failing - paired parameters isolated, only in the acyclic part,
with children distinct from their partners' - therefore persist
throughout the restructuring.

After at most two exchanges, both inputs of $G$ are parameters, call
them $x$ and $y$, and
$G=(x+\alpha)(y+\beta)+\gamma$ for some
$\alpha,\beta,\gamma\in\mathbf F_2$. If the output is a literal of
$G$, Lemma S applies (its parents $x,y$ are affine functions of the
parameters) and we stop. Otherwise $G$ has at least one child, for
(D2) does not apply; every child of $G$ is an acyclic gate, since there
are no wires back to the affine system.

**Case 4: both $x$ and $y$ are unpaired.** Complementing a parameter -
substituting $x\leftarrow w+1$ in $\Phi$ and in $D$, with the negation
absorbed into the constant term of each equation mentioning $x$ and
into the operation of each gate reading $x$ - changes no quantity
$m,q,a,s$ and preserves degrees, hence $W$; apply it to $x$ if
$\alpha=1$ and to $y$ if $\beta=1$, so that $G=xy+\gamma$. Now apply
(PC) to add the pair $\{x,y\}$. On the new restricted domain,
$xy=0$, so $G=\gamma$ throughout $Z(P)$; delete $G$ by (D1),
rerouting it to the constant $\gamma$, and delete one child of $G$ by
(D1), since that child now has a constant input (its rerouting is to a
constant or to a literal of its other input, which already had that
child as a child, so no other unpaired parameter's weight increases).
Both $x$ and $y$ had weight $6$ (each had $G$ as a child and was
unpaired); as paired parameters they lose those weights, while the new
pair contributes $+4$. Altogether $\Delta m=0$, $\Delta s\le -2$,
$\Delta q=+1$, $\Delta a=-1$, and
$\Delta W\le -2+4-6-6=-10=10\,\Delta a$.

**Case 5: at least one of $x,y$ is paired.** Each paired input of $G$
has $G$ as its unique child, and its partner has a unique child that
is an acyclic gate different from $G$ (Cases 1-3 inapplicable, and
Case 2 gives distinctness from $G$ for the partner's child).

*Subcase $k=1$:* exactly one input is paired; say $x$ is unpaired, $y$
is paired with partner $r$, whose unique child is the acyclic gate
$H_r\ne G$. Use (PF) to fix $r=0$, removing the pair $\{y,r\}$. Then
use (UF) to fix the still-unpaired $x$ to the constant $\alpha$, so
$(x+\alpha)=0$ and $G$ computes the constant $\gamma$; delete $G$ by
(D1), rerouting it to $\gamma$, after which $y$, whose only child was
$G$, is unused. Delete $H_r$ by (D1), since it now has the constant
input $r=0$; its rerouting is to a constant or to a literal of its
other input, which already had $H_r$ as a child, so no initially
unpaired parameter increases its weight. Now $\Delta m=-2$ ($r$ and
$x$ are fixed), $\Delta s\le -2$ ($G$ and $H_r$), $\Delta q=-1$,
$\Delta a=-1$. The parameter $x$ loses its weight $6$; for $y$ the
emergent pair is explicit: after (PF) removes its pair, $y$ is
unpaired with unique child $G$ and so transiently has weight $6$
($+6$), and after $G$ is deleted $y$ is unused, weight $0$ ($-6$), a
net contribution of $0$ — $y$ had no individual weight before (paired)
and has weight $0$ after. Hence $\Delta W\le -2-4-6=-12\le 10\,\Delta a$.

*Subcase $k=2$:* both inputs are paired, with partners $r_x$ and
$r_y$; the two pairs are distinct and disjoint, for if $x$ and $y$
were partners of each other then $G$ would be a common child of the
pair $\{x,y\}$ and Case 2 would apply, contrary to assumption. Let $H_{r_x}$ be the
unique child of $r_x$, an acyclic gate different from $G$. Use (PF)
twice to fix $r_x=0$ and $r_y=0$, removing both pairs. Then use (UF)
to fix the now-unpaired $x$ to $\alpha$, so $G$ computes $\gamma$;
delete $G$ by (D1), after which $y$ is unused ($y$'s unique child is
$G$, by the Case 5 preamble, so no other gate, in particular not
$H_{r_x}$, reads $y$). Delete $H_{r_x}$ by
(D1), since it now has the constant input $r_x=0$. (The other
partner's child is left in place; it now has a constant input and
will be removed by (D1) in a later preliminary phase.) Rerouting
increases no initially unpaired parameter's weight, by the same
argument as above. Now $\Delta m=-3$ ($r_x,r_y,x$ are fixed),
$\Delta s\le -2$ ($G$ and $H_{r_x}$), $\Delta q=-2$, $\Delta a=-1$.
The parameters $x,y$ had no individual weight before (paired) and
contribute none after ($x$ removed, $y$ unused), so
$\Delta W\le -2-8=-10=10\,\Delta a$. (Ledger detail: after the two
(PF) steps, $x$ and $y$ are transiently unpaired with child $G$,
$+6$ each; $x$ is then removed by (UF) and $y$ becomes unused when
$G$ is deleted, $-6$ each; the two emergent pairs cancel, and neither
$x$ nor $y$ contributes net weight, exactly as counted above.)

This completes the proof of Lemma I: deletions in the preliminary
phase contribute $\Delta a=0$ and $\Delta W\le -1$ each; the
preparatory moves and (AR) exchanges contribute $\Delta W\le 0$ with
$m,s,a$ unchanged in aggregate accounting; and exactly one of
Cases 1-5 then applies (Cases 4-5 are exhaustive for the minimal AND
gate once Cases 1-3 fail), each satisfying
$\Delta a\in\{0,-1\}$ and $\Delta W\le 10\,\Delta a$. Cases 1-3 and 5
strictly decrease $m$; deletions and Case 4 strictly decrease $s$; no
step increases $m+s$. Hence each iteration strictly decreases the
nonnegative integer $m+s$.

*Proof of Theorem G.* Let $C_0$ compute the disperser $f$ and suppose
$|C_0|<4n-10d$. If $d\ge n$ then $4n-10d<0$, contradicting
$|C_0|\ge 0$ under the supposition; so assume $d<n$. Initially
$m=a=n$, $q=0$, $s=|C_0|$, and
$W_0\le |C_0|+6n<4n-10d+6n=10(n-d)$. Apply Lemma I repeatedly. Suppose
the first $t$ iterations each returned a new state, with values
$a_t,W_t$. Summing $\Delta W\le 10\,\Delta a$ over those iterations
gives $0\le W_t\le W_0+10(a_t-n)<10(a_t-d)$, hence $a_t>d$, so Lemma I
applies to the state after $t$ iterations as well. Since $m+s$
strictly decreases at every iteration and is a nonnegative integer
starting at $n+|C_0|$, the iterations cannot return new states
forever; some iteration must stop in the preliminary phase or at the
minimal AND gate via Lemma S. That gives $U\subseteq Z(P)$ of
dimension at least $d$ and a bit $b$ with $D=b$ on $U$, hence by
($\ast$) the affine subspace $S=\Phi(U)\subseteq\mathbf F_2^n$ of
dimension at least $d$ has $C_0(\vec x)=b$ for all $\vec x\in S$.
Inside $S$, take a point of $S$ together with $d$ linearly independent
direction vectors of $S$: they span a $d$-dimensional affine subspace
$S'\subseteq S$. The disperser $f$ is nonconstant on $S'$, so some
$\vec x\in S'$ has $f(\vec x)\ne b=C_0(\vec x)$, contradicting that
$C_0$ computes $f$. Therefore $|C_0|\ge 4n-10d$. $\square$

Note what the argument never uses: it never assumes additivity under
sharing (all counting is inside the one circuit $D$, with gates
counted once in $s$), never passes to a restricted basis (every step
is valid for all 16 $B_2$ operations, classified exhaustively as
constant, unary, XOR-type or AND-type), and never selects sparse
lengths (Theorem G is applied at the given $n$).

## 10. Sources, provenance, and scope

*The payload disperser.* The nonconstancy argument reproduced in
`disperser.md` is due to Ben-Sasson–Kopparty, *Affine Dispersers from
Subspace Polynomials* (author preprint 2010-03-12),
https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf,
SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf:
the preliminaries, kernel/image polynomials and structural lemma are
their Sections 3.2–3.4 (Lemmas 3.6, 3.7, 3.10, 3.13), the block
decomposition and selection are their Section 7 (Lemma 7.2 and
Theorem 7.3). No statement from that paper is used as a premise: the
exact statement consumed here (Lemma D) is proved in `disperser.md`
from the definitions, specialised to $\mathbf F_2$ and prime field
degree, with the exponents $2^{w_i}-1$, $w_i=r+2-i$, in the decreasing
order the argument uses. The two places where the specialisation
simplifies the original are flagged in `disperser.md`: prime degree
makes the structural lemma's subfield hypothesis automatic, and the
final isolation uses a direct weighted-capacity count. The
certificate-prefix wrapper (field polynomial supplied and verified in
the input) is the standard way to avoid a deterministic
irreducible-polynomial construction; the irreducibility test and the
existence of a valid prefix are proved elementarily in Section 1, and
no result of Shoup is used anywhere.

*The model.* Li-Yang, ECCC TR21-023 (front-page date 2021-02-20),
Sections 1.1, 2.1, 2.4 and 3.3, for the $B_2$ conventions only,
https://eccc.weizmann.ac.il/report/2021/023/download/,
SHA256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062.
The Li-Yang lower-bound theorem and linear program are not used, and
no `lp-certificate` check is declared.

*The gate elimination argument.* Sections 4-9 reproduce, in the form
used here, the parameter-space gate elimination proof of
Mishra-Santhanam, "An Explicit $4n-\mathcal O(\log n)$ Size Lower
Bound against General Boolean Circuits", STOC 2027 submission 1241,
https://stoc2027.sigact-ec.workers.dev/submissions/1241.pdf
(Theorem 4.1 and Lemma 4.2 there; the state operations, fair
semi-circuits, potential $W=s+4q+\sum\mathrm{wt}$, deletions
(D1)-(D4), reparameterisation (AR) and Cases 1-5 correspond to their
Sections 2-4). That work's earlier version, "An Explicit
$3.5n-o(n)$ Size Lower Bound against General Boolean Circuits",
ECCC TR26-260 (19 October 2026),
https://eccc.weizmann.ac.il/report/2026/260/download/,
SHA256 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc,
proves a weaker bound by the same parameter-substitution method. No
authority is borrowed from either source: Sections 4-9 are a
complete, self-contained proof of Theorem G as used here, and the
sources are cited for attribution and provenance only. Two
adaptations are made, and only these: (1) the disperser used is the
payload family of Section 1, with the explicit dimension threshold
$D_0(N)$ and dimension $d(N)=\lfloor D_0(N)\rfloor+1$ of Lemma D,
rather than an $O(\log n)$-dimension disperser — the general theorem
is stated and proved for an arbitrary dimension $d$, and the declared
explicit remainder $R(N)$ is exactly the loss $4k+10\,d(N)$ majorised
as in Section 1, so no $O(\log n)$ disperser construction is needed
or claimed; (2) no constructive refuter or polynomial-time search
claim is made — Lemma I and Theorem G are used purely existentially
(the final step exhibits one miscomputed input, whose existence, not
its efficient computation, contradicts $C_0$ computing $f$).

*Scope.* The bound proved is the declared one and no more:
coefficient $4$ at exponent $1$ with the stated explicit remainder,
for the single family of Section 1, at every integer length $N\ge1$.
Nothing is claimed about any other family, about exponent $>1$, about
an endpoint, or about upper bounds for $F_N$.
