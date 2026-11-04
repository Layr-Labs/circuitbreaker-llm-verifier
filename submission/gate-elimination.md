# The Mishra–Santhanam gate-elimination theorem, reproduced

**Theorem A (Theorem 4.1 of ECCC TR26-260).**
For every affine disperser $f:\{0,1\}^n\to\{0,1\}$ for dimension $d$, every
$B_2$ circuit $C$ computing $f$ has size at least
$$
|C|\;\ge\;\frac{7}{2}\,n-7d .
$$

This file reproduces the complete proof of Theorem A from
O. S. Mishra and R. Santhanam, *An Explicit 3.5n−o(n) Size Lower Bound against
General Boolean Circuits*, ECCC TR26-260 (19 October 2026), Sections 2–4. The
exposition is rewritten but follows the source's definitions, invariants,
potential function, operations and ten cases exactly. Nothing outside
elementary linear algebra over $\mathbb F_2$ is used. We write
$\mathbb F=\mathbb F_2=\{0,1\}$, so addition is XOR, multiplication is AND, and
negation is addition of $1$.

---

## 1. Circuits

A **fan-in-two Boolean circuit** is an acyclic directed graph in which every
vertex has indegree zero or two, and the edges entering a vertex are ordered.
Vertices of indegree zero are **input vertices** and are labelled by variables
$X_1,\dots,X_n$ or by constants $0,1$. Vertices of indegree two are **gates**.
The **fan-out** of a vertex is its outdegree. The two wires entering a gate may
leave the same vertex. Fan-out is unrestricted. Each gate is labelled by a
Boolean operation $g:\mathbb F^2\to\mathbb F$; the set of all sixteen
operations is $B_2$. One vertex of outdegree zero is the **output**. The
**size** $|C|$ of $C$ is its number of gates.

Identify $\mathbb F$ with the field of two elements. Every operation
$h:\mathbb F^2\to\mathbb F$ has a unique multilinear form
$h(t_1,t_2)=c+b_1t_1+b_2t_2+b_{12}t_1t_2$. The **XOR-type** operations are
$t_1+t_2+c$ (two of them); the **AND-type** operations are
$(t_1+\alpha)(t_2+\beta)+\gamma$ (eight of them); the remaining six are
constant or unary operations. An **affine operation** is
$b_1t_1+b_2t_2+c$ with either coefficient allowed to vanish. A gate is
constant, unary, XOR-type or AND-type according to its operation.

We may assume without loss of generality that the circuit has no redundant
gate or input: no gate without a path to the output, no input of outdegree
zero. Removing them does not change the function computed and does not
increase the size. Hence a lower bound on reduced circuits is a lower bound on
all circuits.

---

## 2. Parameters and the parameter map

During the proof we maintain an injective map
$$
\Phi:\mathbb F^m\longrightarrow\mathbb F^n,\qquad
\Phi(\vec u)=\vec b+L\vec u+\sum_{\{p,r\}\in P}\vec c_{p,r}\,p\,r,
\tag{2.1}
$$
where $\vec b\in\mathbb F^n$, $L\in\mathbb F^{n\times m}$, $P$ is a set of
**disjoint** pairs of coordinates of $\vec u$, and $\vec c_{p,r}\in\mathbb F^n$
for each pair. The coordinates $u_1,\dots,u_m$ are **parameters**. A parameter
in a pair of $P$ is **paired**, the other element of its pair is its
**partner**, and all remaining parameters are **unpaired**.

Since $pr$ becomes affine once $p$ or $r$ is fixed to a constant, or $p$ is
identified with a literal of $r$, a pair contributes only one degree of
freedom. Paired parameters may still occur in the affine part $L\vec u$; we
allow $\vec c_{p,r}=0$.

The **main invariant** of the whole argument is
$$
C(\vec u)=f(\Phi(\vec u))\qquad\text{for every }\vec u\in\mathbb F^m,
\tag{2.2}
$$
where $C$ is the current circuit (defined below). Initially $\Phi$ is the
identity: $m=n$, $P=\varnothing$, the parameters are the original inputs, and
$C$ is the original circuit.

**Proposition 2.1.** The image of any injective map of the form (2.1) contains
an affine subspace of dimension $a:=m-|P|$. Consequently $f\circ\Phi$ is
nonconstant whenever $a\ge d$.

*Proof.* Choose one parameter from each pair in $P$ and fix it to a constant.
The pairs are disjoint, so exactly $|P|$ parameters are fixed; every quadratic
term of (2.1) then becomes affine in the remaining $m-|P|=a$ parameters. The
restricted map is the restriction of an injective map, hence injective and
affine, so its image is an affine subspace of dimension $a$. If $a\ge d$, the
affine disperser $f$ is nonconstant on it. $\square$

---

## 3. Fair semi-circuits and the affine system

Let $\vec u=(u_1,\dots,u_m)$ be the parameters. A **fair semi-circuit** on
$\vec u$ consists of three parts.

**(i) Affine system.** $r_0\ge0$ equations in the parameters and in $r_0$
further variables $\vec v=(v_1,\dots,v_{r_0})$, the **auxiliary variables**:
$$
A\vec v=N\vec u+\vec c,\qquad A\in\mathbb F^{r_0\times r_0},\
N\in\mathbb F^{r_0\times m},\ \vec c\in\mathbb F^{r_0}.
\tag{3.1}
$$
$A$ is the **auxiliary matrix**, $N$ the **parameter matrix**. Every row of
$(A\mid N)$ has at most three nonzero entries. **Fairness** means that for
every $\vec u$ there is a unique $\vec v$; equivalently $A$ is invertible.

**(ii) Acyclic part.** An ordinary circuit whose gates are listed in a
topological order. Each gate is labelled by an operation from $B_2$, and each
of its two inputs is a constant, a parameter, an auxiliary variable, or an
earlier gate. There are no wires from the acyclic part back to the affine
system.

**(iii) Output.** A constant, or a literal (a variable or its negation) of a
parameter, an auxiliary variable, or a gate. No variable of an equation is
designated as an output; this is what allows us to exchange an auxiliary
variable with a parameter later.

Parameters, auxiliary variables and gates are **internal values**. By
fairness, every assignment $\vec u$ extends uniquely to all internal values;
in particular every auxiliary variable is an affine function of the
parameters. We write $s$ for the **size** of the semi-circuit: the number
$r_0$ of equations plus the number of gates in the acyclic part. Initially
$r_0=0$ and the acyclic part is the original circuit, so $s=|C|$.

**Moving affine gates into the system.** Suppose a gate $J$ of the acyclic part
is labelled by an affine operation $h(t_1,t_2)=b_1t_1+b_2t_2+c$ and each of
$t_1,t_2$ is a constant, a parameter, or an auxiliary variable. Replace $J$ by
the equation
$$
J+b_1t_1+b_2t_2=c
\tag{3.2}
$$
and regard $J$ as a new auxiliary variable. The new equation has at most three
nonconstant terms, and the new auxiliary matrix is
$\begin{pmatrix}A&0\\ *&1\end{pmatrix}$, which is invertible. For every
parameter assignment, the solution of the old system extends uniquely with $J$
taking the value the gate computed, so the function is unchanged; one gate is
replaced by one equation, so $s$ is preserved. Doing this in topological order
moves every affine gate with no AND-type ancestor into the affine system;
affine gates with an AND-type ancestor stay in the acyclic part.

**Rerouting and representation.** When a gate or auxiliary variable becomes a
constant or a literal of another internal value, we replace all its
occurrences by that constant or literal (**rerouting**). Rerouting a gate of
the acyclic part does not touch the affine system. If an auxiliary variable
$v$ occurs in an equation with at most two nonconstant terms, written
$v=U+c_0$ with $U$ a parameter, an auxiliary variable, or $0$, rerouting $v$
replaces it by $U+c_0$ and deletes the equation; this is an elementary row
operation on $(A\mid N\mid\vec c)$ and hence preserves $\det A$ up to sign,
and it cannot increase the number of nonconstant terms in an equation. So both
invariants — at most three nonconstant terms per row and invertible $A$ — are
preserved. (The detailed determinant argument is in §6, where it is needed.)

---

## 4. States and the four kinds of change of parameters

At every step the **state** is $(C,\Phi,P)$. We track
$$
s=|C|,\qquad m=\text{number of parameters},\qquad q=|P|,\qquad a=m-q .
$$
We often keep parameters in $\Phi$ even when $C$ no longer reads them;
$m$ counts them.

Let $\vec u-x$ denote all parameters other than $x$, and $w$ a fresh parameter
not occurring in $\Phi$. There are four changes.

**(AR) Affine reparameterisation.** For an unpaired parameter $x$ and an
affine map $\varphi$, substitute
$$
x\;\leftarrow\;w+\varphi(\vec u-x).
$$
This does not change $m,q,a$, because $w$ replaces $x$ and no pair changes.

**(QR) Quadratic reparameterisation.** For distinct unpaired parameters
$x,y,z$ and constants $\alpha,\beta,\gamma\in\mathbb F$, substitute
$$
z\;\leftarrow\;w+(x+\alpha)(y+\beta)+\gamma,
\tag{4.1}
$$
and add the pair $\{x,y\}$ to $P$. This keeps $m$, increases $q$ by one, so
decreases $a$ by one.

**(UF) Fixing an unpaired parameter.** For unpaired $x$ and affine $\psi$
(possibly constant), impose $x=\psi(\vec u-x)$. This decreases $m$ by one,
keeps $q$, so decreases $a$ by one.

**(PF) Fixing a paired parameter.** For $\{p,r\}\in P$, either fix $p$ to a
constant or impose $p=r+c$ for a constant $c$, then remove $p$ and the pair.
The same can be done for $r$, the other parameter of the pair. This decreases
both $m$ and $q$ by one, so keeps $a$.

Each change replaces $\Phi$ by $\Phi\circ\tau$, where $\tau$ expresses the old
parameters in terms of the new ones. In (AR) and (QR), $\tau$ is a bijection,
so $\operatorname{im}\Phi$ is unchanged. In (UF) and (PF), $\tau$ fixes one
parameter as a function of the others, hence is injective and
$\operatorname{im}\Phi$ can only shrink. In all cases $\Phi\circ\tau$ is
injective and again has the form (2.1). We make the same change in $C$, so the
new circuit computes $C\circ\tau=f\circ(\Phi\circ\tau)$ and (2.2) is preserved.

The algorithm uses (UF) only to fix a parameter to a constant or to impose
$x=y+c$. Thus each occurrence of a restricted parameter becomes a constant or
a literal of another parameter; this adds no gate, and the gates and equations
that it simplifies are counted in $s$ until they are deleted. Substituting the
expressions of (AR) or (QR) directly into the circuit could add gates, so those
operations are implemented by exchanging roles in the affine system, as
follows.

**Implementing (AR).** Suppose an auxiliary variable $Q$ has the solved
expression $Q=x+\varphi(\vec u-x)$ with $x$ unpaired. Substitute
$x\leftarrow w+\varphi(\vec u-x)$ in $\Phi$, rename $Q$ as the new parameter
$w$, and regard $x$ as an auxiliary variable. The equations are unchanged; we
are only exchanging the column of $x$ in $N$ with the column of $Q$ in $A$.
Let $j$ be the column index of $Q$ in $A$ and $t=A^{-1}N_x$ (the coefficient
vector of $x$ in the solved auxiliary variables). The coefficient of $x$ in
$Q$ is $t_j=1$. The new auxiliary matrix is $A'=AB$, where $B$ is the identity
with its $j$-th column replaced by $t$; $\det B=t_j=1$, so $A'$ is invertible.

**Implementing (QR).** We apply (QR) only in the configuration where $z$ is
isolated (degree one), its unique child is an XOR-type gate $D$, and the other
input of $D$ is an AND-type gate $G=(x+\alpha)(y+\beta)+\gamma$ whose unique
child is $D$. After the substitution (4.1),
$D=G+(w+G+c)+c=w$, so we reroute every use of $D$ to the literal $w$ and
delete $D$. Since $z$ occurred in no equation and had no other child, the
affine system is unchanged; $G$ becomes unused and is deleted. Two gates are
deleted, one pair is added, and the new parameter replaces $z$.

---

## 5. Affine dispersers and non-triviality

An **affine subspace** of $\mathbb F^n$ is a translate $b+V$ of a linear
subspace $V\subseteq\mathbb F^n$. Recall:

**Definition 5.1.** $f:\mathbb F^n\to\mathbb F$ is an **affine disperser for
dimension $d$** if it is nonconstant on every affine subspace of dimension
$d$ (equivalently, of dimension at least $d$).

We use only Proposition 2.1 for non-triviality. The following consequence is
needed to exclude degenerate outputs.

**Lemma 5.2 (Lemma 2.3 of the source).** If $a>d$, the output of the circuit
cannot be an affine function of the parameters, and cannot be a literal of an
AND-type gate whose two parents are unpaired parameters.

*Proof.* A constant output would contradict Proposition 2.1 (the image
contains an affine subspace of dimension $a\ge 1$ on which $f\circ\Phi$ would
be constant).

Suppose the output is a nonconstant affine function of the parameters, and
choose a parameter $x$ with nonzero coefficient. If $x$ is paired, let $y$ be
its partner; fix $y$ to a constant. This is a (PF), so $a$ is unchanged, and
now $x$ is unpaired. Substitute for $x$ an affine expression in the other
parameters chosen so that the output becomes constant. This (UF) decreases
$a$ by one, so the restricted image still contains an affine subspace of
dimension at least $a-1\ge d$, on which $f$ is constant, contradicting
Proposition 2.1.

If the output is a literal of an AND-type gate on two unpaired parameters
$x,y$, fix $y$ to a constant making the gate constant. Again the output is
constant while $a$ decreases by only one, so $a-1\ge d$, a contradiction.
$\square$

---

## 6. Degrees, weights and the potential

**Parents and children.** The **parents** of a gate are its inputs (constants
or internal values); the gate is a **child** of each internal value it reads.
An equation is a child of each internal value appearing in it with nonzero
coefficient, and the internal values in an equation are its parents. A gate
reading the same value twice counts as one child. The **degree** $\deg(t)$ of
an internal value $t$ is its number of distinct children. A parameter is
**unused** if $\deg=0$ and **isolated** if $\deg=1$. Degree is syntactic: it
counts only direct occurrences.

**Weight and potential.** For an unpaired parameter $x$,
$$
\operatorname{wt}(x)=
\begin{cases}
0,&\deg(x)=0,\\
6,&\deg(x)=1,\\
7,&\deg(x)\ge2,
\end{cases}
\qquad
W=2s+3q+\sum_{x\ \text{unpaired}}\operatorname{wt}(x).
\tag{6.1}
$$
Paired parameters contribute only through $3q$. The potential is always
nonnegative. For any tracked quantity $X$ we write $X'$ for its value after an
operation and $\Delta X=X'-X$.

**Rerouting.** When a gate or auxiliary variable becomes a constant or a
literal of another internal value, replace all its occurrences by that
constant or literal (**rerouting**). Rerouting can increase the degree of a
parameter that was a parent of the replaced value; in particular it can turn
an isolated parameter into one of degree at least two.

**Preservation of the representation.** Fixing $x=c_0$ deletes the column of
$x$ in $N$, updates $\vec c$, and leaves $A$ unchanged. Setting $x=y+c_0$
likewise leaves $A$ unchanged and does not increase the number of nonconstant
terms in an equation. For (AR), let $j$ be the column of $Q$ in $A$,
$N_x$ the column of $x$ in $N$, and $t=A^{-1}N_x$. The coefficient of $x$ in
$Q$ is $t_j=1$ (because $Q=x+\varphi(\vec u-x)$), and $A'=AB$ with $B$ the
identity modified in column $j$ by $t$; $\det B=t_j=1$, so $A'$ is invertible
and every equation keeps exactly the same terms. For (QR), $z$ occurs in no
equation (it is isolated and its only child is a gate $D$), so the affine
system is unchanged. For rerouting an auxiliary variable $v$ that occurs in an
equation $v=U+c_0$: replacing $v$ by $U+c_0$ in another equation is an
elementary row operation on $A$ (it adds the row of $v=U+c_0$), hence preserves
$\det A$; after the substitution the column of $v$ in $A$ has a single nonzero
entry, and expanding the determinant along it shows the reduced auxiliary
matrix is invertible. The same argument applies when deleting an auxiliary
variable whose only use is in a single equation. Thus every operation below
keeps the two structural invariants.

---

## 7. Deletions

In each deletion $\Delta m=\Delta q=\Delta a=0$ and $\Delta s=-1$. Hence the
$2s$ term of $W$ contributes $-2$; the rest is the change of weights.

**D1: a gate depending on at most one input.** This covers every constant or
unary gate, every gate with a constant input, every gate whose two inputs are
identical, and every AND-type gate whose two inputs are parameters or
auxiliary variables with dependent linear parts. Reroute every use of the gate
$J$ to the constant or the literal of its input, and delete $J$. An unpaired
parameter that was a parent of $J$ increases its weight from $6$ to $7$ only
if it was isolated and has at least two children after rerouting. Hence
$\Delta W\le-2-6+7=-1$ in that case; if $J$ is rerouted to a constant no
unpaired parameter can increase its weight, so $\Delta W\le-2$.

**D2: an unused gate.** A gate with no children that is not the output can be
deleted without changing the function. $\Delta W\le-2$.

**D3: an equation with at most two nonconstant terms.** Such an equation must
contain an auxiliary variable, otherwise its row in $A$ would be zero. Write
it $v=U+c$ with $U$ a parameter, an auxiliary variable, or $0$; substitute
$v=U+c$ everywhere, and delete $v$ and the equation. As in D1, at most one
isolated unpaired parameter can increase its weight from $6$ to $7$, so
$\Delta W\le-1$ (and $\le-2$ if no weight increases).

**D4: unit parameter deletion.** Let $x$ be an unpaired isolated parameter
whose unique child is an equation $E$. Write the system as $A\vec v=N\vec u+
\vec c$. Since $x$ occurs only in $E$, the column $N_x$ has a single nonzero
entry, so $A^{-1}N_x\neq0$. Choose an auxiliary variable $Q$ whose coordinate
in $A^{-1}N_x$ is $1$; then $Q=x+\varphi(\vec u-x)$ for some affine
$\varphi$. Apply (AR) exchanging $x$ and $Q$. No equation changes, so the new
auxiliary variable $x$ still occurs only in $E$. Delete $x$ and $E$. The old
parameter $x$ of weight $6$ is replaced by the new parameter $Q$ of weight at
most $7$, no other unpaired parameter increases its weight, and $s$ drops by
one, so $\Delta W\le-2-6+7=-1$.

---

## 8. Reparameterisations

Both (AR) and (QR) are bijective changes of parameters that preserve both the
image and the form of $\Phi$. We record their effect.

**AR.** Let $Q$ be an auxiliary variable depending on an unpaired parameter
$x$, i.e. $Q=x+\varphi(\vec u-x)$. ($Q$ and $x$ need not share an equation,
but $x$ occurs in some equation, otherwise the column of $x$ in $N$, hence in
$A^{-1}N$, would be zero and $Q$ could not depend on $x$.) Substitute
$x\leftarrow w+\varphi(\vec u-x)$ in $\Phi$, regard $Q$ as the new parameter
$w$ and $x$ as an auxiliary variable. No equation or gate changes; every
parameter other than $x,w$ has the same children. Hence
$\Delta m=\Delta q=\Delta a=\Delta s=0$ and $\Delta W=\operatorname{wt}'(w)-
\operatorname{wt}(x)$. Since $w$ has exactly the children of $Q$,
$\deg'(w)=\deg(Q)\ge1$ (the column of $Q$ in $A$ is nonzero). If $Q$ has at
least one other child and $\deg(x)\ge2$, then both weights are $7$ and
$\Delta W=0$. This is the only situation in which we use (AR), apart from D4.

**QR.** We apply (QR) only as in §4: $x,y,z$ distinct unpaired, $G=(x+\alpha)
(y+\beta)+\gamma$ AND-type with only child $D$, $D=G+z+c$ XOR-type, and $z$
isolated with only child $D$. Then $z\leftarrow w+G(x,y)+c$, add the pair
$\{x,y\}$, reroute $D$ to $w$, delete $D$, and delete the now unused $G$. We
have $\Delta m=0$, $\Delta q=1$, $\Delta a=-1$, $\Delta s=-2$. For the
potential: the two deletions give $-4$; the new pair gives $+3$; $x,y$ lose
their weights; the new parameter $w$ has exactly the children of $D$ and
weight at most $7$, replacing the isolated $z$ of weight $6$; no other
unpaired parameter gains a child. Thus
$$
\Delta W\le -4+3+7-6-\operatorname{wt}(x)-\operatorname{wt}(y)
=-\operatorname{wt}(x)-\operatorname{wt}(y).
$$
We apply QR only when $\deg(x),\deg(y)\ge2$, so both weights are $7$ and
$\Delta W\le-14=14\Delta a$.

---

## 9. The algorithm and Lemma 4.2

The algorithm runs until $a\le d$. Each iteration first does preparatory
operations that preserve $m,s,a$ and do not increase $W$:

* move every affine gate with no AND-type ancestor into the affine system
  (each move keeps $s$ and does not increase any degree, so $\Delta W\le0$);
* at most two (AR) exchanges to turn the inputs of the minimal AND into
  parameters;
* reroutings to remove children of the minimal AND that also read one of its
  inputs.

Then it performs the first applicable reduction among the cases below.

**Lemma 9.1 (Lemma 4.2 of the source).** Suppose $a>d$. In polynomial time the
algorithm transforms the state so that $m+s$ strictly decreases,
$\Delta a\in\{0,-1,-2\}$, and $\Delta W\le14\Delta a$.

*Proof.* The cases below each verify the three claims. The preparatory
operations preserve $m,s,a$ and do not increase $W$, so the claims for an
iteration follow from the claims for the chosen reduction. Every reduction
either deletes gates/equations or restricts/reparameterises parameters, never
increases $s$, and strictly decreases $m+s$; termination follows.

### 9.1 Preliminary reductions

**Case 1: a deletion is available.** Apply D1–D4 repeatedly until none
applies. Each has $\Delta a=0$, $\Delta s=-1$, $\Delta W\le-1$.

If no deletion applies, the state has these properties:

* (from D1) every gate is XOR-type or AND-type, has two distinct inputs, and
  neither input is a constant; if both inputs of an AND-type gate are
  parameters or auxiliary variables, their linear parts are independent;
* (from D2, D3) every gate other than the output has a child, and every
  equation has exactly three nonconstant terms;
* (from D4) no unpaired isolated parameter occurs in the affine system.

**Case 2: an AND-type gate on parameters has an unpaired isolated input.**
Let $G$ be an AND-type gate with inputs the parameters $z$ and $w$, where $z$
is unpaired and isolated. Choose a constant $c$ such that $G$ becomes constant
when $w=c$, and restrict $w$ to $c$: this is (UF) if $w$ is unpaired and (PF)
if $w$ is paired. Delete $G$ by D1 and reroute all occurrences of $G$ to that
constant. Since $G$ was the only child of $z$, $z$ is now unused; no other
degree changes because $G$ is rerouted to a constant.

If $w$ was unpaired, both $z$ and $w$ had weight at least $6$. Thus
$$
\Delta m=-1,\quad\Delta s=-1,\quad\Delta q=0,\quad\Delta a=-1,
\quad\Delta W\le-6-6-2=-14.
$$
If $w$ was paired, let $r$ be its partner. The pair is removed and $r$ becomes
unpaired with weight at most $7$. Thus $\Delta m=-1$, $\Delta s=-1$,
$\Delta q=-1$, $\Delta a=0$, and
$$
\Delta W\le-3-6-2+7=-4.
$$

**Consequence of Cases 1–2.** For every isolated unpaired parameter $z$: it is
not the output (Lemma 5.2), it does not occur in the affine system (Case 1
inapplicable), so its unique child is a gate $G$ of the acyclic part. The other
input of $G$ cannot be a parameter: if $G$ were XOR-type with two parameter
inputs it would have been moved into the affine system, and if $G$ were
AND-type then Case 2 would apply. Hence no gate or equation has both an
isolated unpaired parameter and a different parameter as parents.

**Case 3: a paired parameter has degree zero or at least two.**
Let $p$ be paired with partner $r$.

If $\deg(p)=0$, apply (PF) to fix $r$ to a constant. The pair is removed and
$p$ becomes unpaired but unused, hence contributes $0$. Thus
$\Delta m=-1$, $\Delta s=0$, $\Delta q=-1$, $\Delta a=0$, and
$\Delta W\le-3$.

If $\deg(p)\ge2$, apply (PF) to fix $p$ to a constant, and delete two children
of $p$ using D1–D4. These deletions do not increase the weight of any
initially unpaired parameter; only $r$, newly unpaired, can reach weight at
most $7$. Thus $\Delta m=-1$, $\Delta s=-2$, $\Delta q=-1$, $\Delta a=0$, and
$$
\Delta W\le-3-4+7=0.
$$
If Case 3 is inapplicable, every paired parameter is isolated and has exactly
one child.

**Case 4: two partners have the same child.** Let $p$ be paired with $r$, and
let $H$ be a common child. If $H$ is an equation, it reads $v+p+r=c$ for an
auxiliary variable $v$. If $H$ is a gate, its two inputs are $p,r$. Use (PF) to
impose $p=r+c'$, choosing $c'$ so that $H$ becomes constant: any $c'$ works
for an equation or an XOR-type gate; for an AND-type gate
$(p+\alpha)(r+\beta)+\gamma$, take $c'=\alpha+\beta+1$. Delete $H$ and reroute
its occurrences to the constant. The parameter $r$ survives but is now unused,
since $H$ was its only child. Thus $\Delta m=-1$, $\Delta s=-1$,
$\Delta q=-1$, $\Delta a=0$, and $\Delta W\le-3-2=-5$.

**Case 5: a paired parameter occurs in the affine system.** Let $p$ be paired
with $r$ and suppose $p$ occurs in the affine system. Since Case 3 does not
apply, $p$ has a unique child $H_p$; since $p$ occurs in the affine system,
$H_p$ is an equation. Similarly $r$ has a unique child $H_r$, and $H_r\neq
H_p$ because Case 4 does not apply.

Apply (PF) to fix $r$ to a constant, and delete $H_r$ (by D1 if it is a gate,
by D3 if it is an equation). Since $r$ is a parent of $H_r$ and Cases 1–2 do
not apply, $H_r$ has no isolated unpaired parameter as a parent; hence deleting
$H_r$ does not increase the weight of any initially unpaired parameter. After
$r$ is fixed, $p$ becomes unpaired with unique child $H_p$ and gets weight $6$.
Thus removing the pair contributes $-3$, deleting $H_r$ contributes $-2$, and
the new weight of $p$ contributes $+6$, so
$\Delta W\le-3-2+6=1$. Now $p$ is unpaired and isolated with unique child the
equation $H_p$, so D4 applies and gives $\Delta W\le-1$. Altogether
$\Delta W\le0$, while $\Delta m=-1$, $\Delta s=-2$, $\Delta q=-1$,
$\Delta a=0$.

**Case 6: an AND-type gate on parameters has a paired input.** Let $G$ be an
AND-type gate with inputs the parameters $p$ and $w$, where $p$ is paired with
$r$. Since Case 3 does not apply, $G$ is the unique child of $p$. The partner
$r$ has a unique child $H_r$. Since Case 5 does not apply, $r$ does not occur
in the affine system, so $H_r$ is a gate. Since Case 4 does not apply,
$H_r\neq G$ and $w\neq r$. Also $w\neq p$ because Case 1 does not apply.
Since Case 2 does not apply, $w$ is not unpaired isolated; hence $w$ is either
unpaired of degree at least two, or paired.

Apply (PF) to fix $r$ to a constant, then fix $w$ to a constant so that $G$
becomes constant (UF if $w$ is unpaired, PF if paired); delete $G$ and $H_r$
by D1. After $r$ is restricted, $p$ is unpaired and unused (its only child was
$G$). No initially unpaired parameter increases its weight: $G$ is rerouted to
a constant, and $H_r$ has no isolated unpaired parameter as a parent (Cases 1
and 2 inapplicable).

If $w$ was unpaired, its degree was at least two, so it contributed $7$ and
now contributes $0$. Then $\Delta m=-2$, $\Delta s=-2$, $\Delta q=-1$,
$\Delta a=-1$, and
$$
\Delta W\le-7-3-4=-14=14\Delta a.
$$
If $w$ was paired with partner $r'$, two pairs are removed; only $r'$ can
become newly unpaired with positive weight, at most $7$. Then $\Delta m=-2$,
$\Delta s=-2$, $\Delta q=-2$, $\Delta a=0$, and
$$
\Delta W\le-6-4+7=-3.
$$

**Consequence of Cases 1–6.** If none applies, then every parameter occurring
in the affine system is unpaired of degree at least two, and for every
AND-type gate $G$ with two parameter inputs, both parameters are unpaired and
have degree at least two.

### 9.2 Restructuring the minimal AND

Assume $a>d$ and none of Cases 1–6 applies. The acyclic part has at least one
gate: otherwise the output is affine in the parameters, contradicting
Lemma 5.2. Choose the first gate $G$ in the acyclic part in topological order.
It must be AND-type: all affine gates whose two inputs are already in the
affine system were moved at the start of the reduction. Call $G$ the
**minimal AND**. By Case 1 inapplicability, both inputs of $G$ are parameters
or auxiliary variables with independent linear parts.

**Turning the inputs into parameters.** Suppose an input of $G$ is an
auxiliary variable $Q$. Since the affine system is fair, $Q$ is an affine
function of the parameters. Its linear part is nonconstant because the inputs
of $G$ are independent, so $Q$ depends on some parameter $x$. That $x$ must
occur in some equation: otherwise the column of $x$ in $N$, and hence in
$A^{-1}N$, would be zero, contradicting that $Q$ depends on $x$. Since no
preliminary case applies, every parameter occurring in the affine system is
unpaired of degree at least two; in particular $x$ is unpaired with
$\deg(x)\ge2$. Apply (AR) to exchange $Q$ and $x$. The new parameter $w$ has
degree at least two (the column of $Q$ in $A$ is nonzero and $Q$ is a parent
of $G$), and both $x$ and $w$ have weight $7$, so
$\Delta m=\Delta s=\Delta q=\Delta a=\Delta W=0$.

Move all affine gates into the system again, and inspect Cases 1–6 in order;
if one applies, perform it and end the iteration. Otherwise repeat for the
other input if it is an auxiliary variable. There is a parameter it depends on
other than $w$, because the two inputs of $G$ were linearly independent and an
invertible affine change of parameters preserves linear independence. After at
most two applications of (AR), both inputs of $G$ are parameters, denoted
$x,y$. Since Cases 2 and 6 do not apply, $x$ and $y$ are unpaired of degree at
least two.

**Rerouting children that read $x$ or $y$.** Suppose a child $J$ of $G$ also
has parameter $x$ as a parent. Since there are no wires from the acyclic part
to the affine system, $J$ is a gate of the acyclic part. The value at $J$ is
$h(x,G(x,y))$ for some operation $h\in B_2$ depending only on $x$ and $y$.
Relabel $J$ by the operation $(x,y)\mapsto h(x,G(x,y))$ and reroute its inputs
to $x$ and $y$. The function is unchanged; $J$ is no longer a child of $G$ but
is still a child of $x$; $y$ may gain occurrences but had degree at least two,
so its weight stays $7$. Hence $\Delta m=\Delta s=\Delta q=\Delta a=0$ and
$\Delta W\le0$. Repeat with $x,y$ interchanged. Each rerouting removes a child
of $G$ and creates no new child of $G$, so the process terminates. If it ends
with no preliminary case applicable, then no child of $G$ reads $x$ or $y$.
By Lemma 5.2, $G$ is not the output and, since Case 1 does not apply, it has
at least one child.

### 9.3 Eliminating the minimal AND

Assume no child of $G$ reads $x$ or $y$. Swap the names so that
$\deg(x)\ge\deg(y)\ge2$, and let $k\ge1$ be the number of children of $G$.

**Case 7: $\deg(x)=\deg(y)=2$ and $x,y$ have a common second child.**
Let $B\neq G$ be the common second child.

If $B$ is an equation, it has exactly three nonconstant terms (Case 1
inapplicable) and at least one auxiliary variable (invertibility of $A$), so
$B$ reads $v+x+y=c$ for an auxiliary variable $v$ and a constant $c$. If $B$
is not an equation, it is a gate with inputs $x,y$. Let $h(x,y)=x+y$ in the
first case and the operation of $B$ in the second.

The AND-type gate $G$, as a function of $(x,y)\in\mathbb F^2$, takes the same
value on three of the four points. Since $h$ takes at most two values, two of
those three points have the same $h$-value and the same $G$-value. Two points
of $\mathbb F^2$ lie on a line $\alpha x+\beta y=\gamma$ with $(\alpha,\beta)
\neq(0,0)$. Restrict to this line by one (UF): fix $x=\gamma$ if $\beta=0$,
fix $y=\gamma$ if $\alpha=0$, or impose $x=y+\gamma$ if $\alpha=\beta=1$. On
this line $G$ is constant and $B$ is constant if it is a gate; if $B$ is an
equation it becomes $v=c_0$. Delete $G$ by D1 and $B$ by D1 or D3, rerouting
each to its constant. One of $x,y$ is removed and the other becomes unused,
since all occurrences of $x,y$ were in $G$ and $B$; both lose weight $7$, and
no other parameter increases its weight because we reroute only to constants.
Thus $\Delta m=-1$, $\Delta s=-2$, $\Delta q=0$, $\Delta a=-1$, and
$$
\Delta W\le-7-7-4=-18=18\Delta a\quad(\le14\Delta a).
$$
From now on assume this case is not applicable; in particular, if
$\deg(x)=\deg(y)=2$, the second children $B_x,B_y$ of $x,y$ are distinct.

**Case 8: $\deg(x)\ge3$, or $G$ has at least two children, or its unique
child has no unpaired isolated input.** Let $k_1$ be the number of children of
$G$ whose other parent is an unpaired isolated parameter, so $0\le k_1\le k$.
Fix $x$ to a constant making $G$ constant (UF), then delete all $\deg(x)$
children of $x$, starting with $G$, and then all $k$ children of $G$. These
$\deg(x)+k$ internal values are distinct because no child of $G$ reads $x$
(established in §9.2). Once $x$ is fixed, every child of $x$ is deletable by
D1 if a gate and D3 if an equation; once $G$ is rerouted to a constant, every
child of $G$ is a gate of the acyclic part, deletable by D1.

The parameter $x$ had degree at least two, so it loses weight $7$. Each of the
$\deg(x)+k$ deletions contributes $-2$. An unpaired parameter of degree at
least two already has the maximum weight $7$ and cannot increase. An
initially unused parameter remains unused. An initially isolated unpaired
parameter increases its weight only if its unique child is deleted and
rerouted to a literal of it; by the consequence of Cases 1–2 its unique child
is not a child of $x$, so it must be one of the $k_1$ children of $G$ whose
other parent is an unpaired isolated parameter. Hence at most $k_1$ parameters
increase, each by at most one. Therefore
$$
\Delta m=-1,\quad\Delta s=-(\deg(x)+k),\quad\Delta q=0,\quad\Delta a=-1,
$$
and
$$
\Delta W\le-7-2(\deg(x)+k)+k_1. \tag{9.1}
$$

If $\deg(x)=2$, then $\deg(y)=2$ as well. The children of $y$ are $G$ and
$B_y$; $y$ is not a parent of $B_x$ because $B_x\neq B_y$, and no child of $G$
has $y$ as a parent. As $G$ is rerouted to a constant, $y$ gains no new child
and loses the child $G$, so it becomes isolated: its weight drops from $7$ to
at most $6$, an extra decrease of $1$. Combining,
$$
\Delta W\le
\begin{cases}
-7-2\deg(x)-k\le-14, & \deg(x)\ge3,\\[2pt]
-12-2k+k_1\le-14, & \deg(x)=2,\ k\ge2,\\[2pt]
-14+k_1, & \deg(x)=2,\ k=1.
\end{cases}
$$
(The second line uses $k_1\le k$ and $k\ge2$; the first uses $k\ge1$.)

Since $\Delta a=-1$, Case 8 satisfies $\Delta W\le14\Delta a$ — except when
$\deg(x)=\deg(y)=2$ and $k=k_1=1$. In that remaining situation $G$ has a
unique child $D$, whose other parent is an unpaired isolated parameter $z$;
this is the **bottleneck configuration**. Since Case 1 does not apply, $D$ is
either AND-type or XOR-type. We test for this configuration before applying
Case 8 and, if it occurs, use Case 9 or Case 10 instead.

**Case 9: bottleneck configuration with $D$ AND-type.** Write
$D=(G+\delta)(z+\eta)+\theta$ with constants $\delta,\eta,\theta$. Since $G$
is nonconstant, choose $(x_0,y_0)$ with $G(x_0,y_0)=\delta$. Fix $x=x_0$ and
$y=y_0$ by two (UF) restrictions. Then $G=\delta$ and $D=\theta$. Delete $G$
and $D$ by D1, rerouting them to their constants, and delete $B_x$ and $B_y$
(by D1 if a gate, by D3 if an equation), which is possible since $x,y$ are
fixed. These four internal values are distinct: $B_x\neq B_y$, and neither
equals $D$ because no child of $G$ reads $x$ or $y$. The only child of $z$ was
$D$, which is rerouted to a constant, so $z$ becomes unused. An unpaired
parameter can increase its weight only if it is initially isolated and its
unique child is deleted and rerouted to a literal of it; this does not happen:
$G,D$ are rerouted to constants, and $B_x,B_y$ have $x,y$ as parents, so by the
consequence of Cases 1–2 they have no isolated unpaired parameter as parent.
Thus $x,y,z$ lose weights $7,7,6$ and the four deletions contribute $-8$:
$$
\Delta m=-2,\quad\Delta s=-4,\quad\Delta q=0,\quad\Delta a=-2,
\quad\Delta W\le-7-7-6-8=-28=14\Delta a.
$$

**Case 10: bottleneck configuration with $D$ XOR-type.** Write $D=G+z+c$
with $c\in\mathbb F$. Apply (QR) as in §4:
$$
z\leftarrow w+G(x,y)+c,
$$
add the pair $\{x,y\}$, reroute $D$ to $w$ and delete $D$, then delete the now
unused $G$. Since the only child of $z$ was $D$, no equation and no other gate
is affected. Two gates are deleted ($-4$); the new pair contributes $+3$;
$x,y$ lose their weights $7,7$; the new parameter $w$ has the former children
of $D$ (weight at most $7$) and replaces the isolated $z$ of weight $6$. No
other unpaired parameter gains a child. Hence
$$
\Delta m=0,\quad\Delta s=-2,\quad\Delta q=1,\quad\Delta a=-1,
\quad\Delta W\le-4-7-7-6+3+7=-14=14\Delta a.
$$

**End of the proof of Lemma 9.1.** Every possible state with $a>d$ falls into
one of Cases 1–10 or the restructuring of §9.2. The preparatory operations
preserve $m,s,a$ and do not increase $W$; each reduction has
$\Delta a\in\{0,-1,-2\}$ and $\Delta W\le14\Delta a$, and strictly decreases
$m+s$ (deletions and QR decrease $s$ and keep $m$; the other reductions
decrease $m$; none increases $s$). An iteration performs at most $s$ gate
moves, at most two affine exchanges, and at most one rerouting per child of
the minimal AND. We use Gaussian elimination over $\mathbb F$ to express the
auxiliary variables as affine functions of the parameters, which lets us test
linear dependence and perform (AR); all other tests and updates involve lists
of gates and equations and constant-size truth tables. Hence each iteration
runs in polynomial time. $\square$

---

## 10. Proof of Theorem A

Initialise the state with $\Phi(\vec u)=\vec u$, no pairs, an empty affine
system, and the input circuit as the acyclic part. Then $m=a=n$, $q=0$,
$s=|C|$, and since every input is an unpaired parameter of weight at most $7$,
$$
W_0\le 2|C|+7n. \tag{10.1}
$$

While $a>d$, apply Lemma 9.1: try the preliminary reductions, and if none
applies follow §9.2 and §9.3. By Lemma 9.1 each iteration strictly decreases
$m+s$; since $m+s=n+|C|$ initially and $m+s\ge0$ throughout, the algorithm
stops with $a\le d$ after at most $n+|C|$ iterations. The map $\Phi$ is stored
through $\vec b,L,P,\vec c_{p,r}$ in (2.1), all updated in polynomial time;
since $m\le n$ and $s\le|C|$ throughout, each iteration is polynomial in
$n$ and $|C|$, so the whole algorithm is polynomial time.

Let $T$ be the number of iterations, so $a_T\le d$. Summing
$\Delta W\le14\Delta a$ over all iterations,
$$
W_T-W_0\le14(a_T-n).
$$
Since $W_T\ge0$ and $a_T\le d$,
$$
2|C|+7n\ \ge\ W_0\ \ge\ W_T+14(n-a_T)\ \ge\ 14(n-d).
$$
Rearranging gives
$$
|C|\ \ge\ \frac{14(n-d)-7n}{2}\ =\ \frac72\,n-7d.
$$
This proves Theorem A. $\square$

---

## 11. Remark on the leading coefficient

For potentials of the form $V=s+\beta q+\alpha I-\sigma\ell$, where $I$ is the
number of used unpaired parameters and $\ell$ the number of isolated unpaired
parameters, with $\alpha\ge\sigma\ge0$ and $\beta\ge0$, checking the estimates
above shows that $\alpha=7/2$, $\beta=3/2$, $\sigma=1/2$ are optimal: no
choice of weights improves the coefficient $7/2$. This is a statement about
this family of potential functions only. It is not a barrier, and we do not
claim that $7/2$ is optimal for unrestricted $B_2$ circuits.
