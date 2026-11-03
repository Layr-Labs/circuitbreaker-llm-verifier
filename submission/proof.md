# Coefficient 2, exponent 1 for the threshold-2 family: a self-contained gate-elimination proof

## Theorem claimed

Let $L$ be the set of binary strings containing at least two $1$s.
For $n \ge 0$ let $f_n(x_1,\dots,x_n)=1$ if and only if
$\sum_{i=1}^n x_i \ge 2$, and $f_n=0$ otherwise. Thus $f_0=0$,
$f_1=0$, $f_2=x_1 \land x_2$, and $f_3$ is majority of three bits.
In the full Li-Yang B2 model, with $C_{B2}(f_n)$ the minimum number
of gates in a circuit computing $f_n$ exactly on every input in
$\{0,1\}^n$,

$$
C_{B2}(f_n) \ge 2n-4 \qquad\text{for every integer } n \ge 0.
$$

In particular, with $g(n)=2\cdot n^1$,

$$
\liminf_{n\to\infty} \frac{C_{B2}(f_n)}{g(n)} \ge 1,
$$

so the declared finite pair $(c,p)=(2,1)$ holds. The proof below
uses no external circuit theorem, no literature construction, and
no unproved premise: every step is proved here. It is a
gate-elimination proof for one fixed threshold, not a claim of
coefficient above $2$, of exponent above $1$, of a lane or global
endpoint, or of the $31/10$ baseline. Tightness of $2n-4$ is not
claimed.

## 1. One fixed total family, in P and hence in NP

The language is fixed once, for all lengths, by the single sentence
"$x$ contains at least two $1$s". There is no per-length advice,
no truth-table oracle, no block, field, or padding encoding, and
no choice that depends on $n$ other than the length itself. The
function is total: every binary string has a defined value,
including the empty string and one-bit strings, for which the
value is $0$.

Define a verifier $V$ and a witness bound $W$ as follows. For every
$n$, $W(n)=0$, a constant, hence polynomial, bound. On input
$(x,w)$, $V$ accepts if and only if $w$ is the empty string and a
left-to-right scan of $x$, stopping as soon as two $1$s have been
seen, finds at least two $1$s. Scanning $x$ takes time linear in
$|x|$, and checking that $w$ is empty takes constant time, so $V$
runs in polynomial time in $|x|$.

Both directions of NP membership hold. If $x \in L$, then the empty
witness $w$ is accepted, because the scan finds two $1$s in $x$.
Conversely, if $V(x,w)=1$, then $w$ is empty and the scan found
two $1$s in $x$, so $x \in L$. The same scan, with no witness, is
a deterministic polynomial-time evaluator for the family, so $L$
is in P, hence in NP. The witness is part of the eligibility
argument only. Circuits in Sections 3-5 receive $x$ alone, never
a witness, and the witness bound is $0$ at every length.

For $n \ge 3$, $f_n$ and its negation $\lnot f_n$ both depend on
every variable. Indeed, fix all variables other than $x_i$ so that
exactly one of them is $1$ and the rest are $0$. Then $f_n=0$ when
$x_i=0$ (total $1$) and $f_n=1$ when $x_i=1$ (total $2$); negation
reverses the two values but still flips. This dependence is used
repeatedly below. It also shows, for $n \ge 3$, that no circuit
computing $f_n$ or $\lnot f_n$ can have an input or a constant as
its output source: such a source depends on at most one variable.
Hence, for $n \ge 3$, a correct circuit for $f_n$ or $\lnot f_n$
must have a gate as its output source.

## 2. The circuit model and what is counted

Fix $n$ and a circuit in the model of rules.md: deterministic and
acyclic, one output, two input pins per gate, every one of the
$16$ binary Boolean truth tables allowed, unrestricted depth and
fan-out, exact computation on every input in $\{0,1\}^n$, and
nonuniform circuits allowed. The argument below is made for an
arbitrary such circuit at each $n$, so nonuniformity is no
obstacle: uniformity is required of the family and verifier in
Section 1, not of circuits.

Counting follows rules.md. Each retained gate is counted once.
Inputs, wires, and output selection cost zero. A constant source,
if the model represents constants as nodes, costs zero and is not
a gate. Repeated pins are allowed. Shared subgraphs are allowed:
one gate may feed many descendants, and is still counted once.
No monotone, bounded-depth, XOR-free, formula, or other
restricted-basis assumption is used, and no transfer from a
restricted model is needed: the argument is made directly in full
B2, for all $16$ truth tables. In particular, no step assumes a
gate table is monotone, an AND/OR/XOR table, or anything other
than a binary table.

$C_{B2}(f_n)$ is the minimum gate count over all such circuits
computing $f_n$. The minimum is over a nonempty set, hence is
attained: every Boolean function on $n$ variables has a B2
circuit, for example via its disjunctive normal form using AND,
OR, and negated-projection gates, each of which is one of the
$16$ binary tables. No explicit upper bound is needed for the
lower bound claimed here, and none is claimed.

A node is an input node, a constant node, or a gate. A gate has
two pins; each pin's source is a node, and acyclicity means no
gate is a pin source, directly or indirectly, of itself. The
output source is one node, selected at zero cost.

## 3. Eliminable gates

A gate $g$ is called eliminable if at least one of the following
syntactic conditions holds:

- the gate's table depends on at most one of its two pin values
  (a constant table, or a projection or negated projection);
- at least one pin source is a constant node;
- both pins have the same source node.

In each case the output of $g$ is a unary function of a pin
source $a$ of $g$: if both pins have the same source $s$, take
$a=s$ and $u(t)$ to be the table evaluated at $(t,t)$; else if
a pin source is a constant node $c$, take $a$ to be the other
pin source and $u(t)$ to be the table evaluated with $c$ in the
constant position; else the table depends on at most one pin,
so take $a$ to be the pin source of that pin (either pin source
if the table depends on neither, with $u$ constant) and $u$ to
be that unary function. In every case $a$ is a pin source of
$g$, so $a \ne g$ by acyclicity and $a$ precedes $g$
topologically. This syntactic requirement is essential: a gate
is not called eliminable merely because its output happens to
equal a unary function of some later node, which would not
precede $g$ and could not be substituted for $g$ without a cycle.

A unary function here is not a gate and is never inserted as one.
It is only a description of what a binary gate computes in one
of the cases above. In particular, a gate whose table is, say,
AND, but whose one pin is a constant node, is eliminable in this
sense, even though the AND table as a function of two variables
depends on both variables. Eliminability is a property of the
gate in the circuit, with its actual pin sources, not of its
table in isolation.

Lemma 1 (single elimination). Let $D$ be a circuit computing a
nonconstant function $F$, and let $g$ be an eliminable gate of
$D$. Then there is a circuit $D'$ with exactly one fewer gate
than $D$, computing $F$ or $\lnot F$.

Proof. Write the output of $g$ as $u(a)$ for a pin source $a$
of $g$ and a unary $u$, as in the definition above. Then
$a \ne g$ by acyclicity and $a$ precedes $g$ topologically.

If $g$ is not the output gate, construct $D'$ as follows. Delete
$g$. For every successor gate $h$ that had $g$ as a pin source,
replace that pin by $a$ and replace the table of $h$ by the
composed table: the new table, as a function of $a$ (or of the
constant value, if $a$ is a constant node) and the other pin of
$h$, equals the old table as a function of $u(a)$ and the other
pin. If both pins of $h$ were $g$, both are replaced in the same
way. The composed table is still a binary Boolean table, so $h$
remains one B2 gate. If $g$ had no successors, $D'$ is just $D$
with $g$ deleted. Acyclicity is preserved, because $a$ precedes
$g$ and hence precedes every successor of $g$. The values of all
remaining gates, and therefore the output, are unchanged, so $D'$
computes $F$. The gate count drops by exactly one: $g$ is gone,
and every successor is still counted once, with a changed table
only.

If $g$ is the output gate, first apply the same successor
replacement and table composition as above to every successor
gate of $g$ (such successors, if any, are dead with respect to
the output, since a successor of the output gate cannot also be
an ancestor of it without a cycle, but their pins must still
not be left pointing at the deleted gate), and then delete $g$
and take $a$ as the new output source. If
$u=\mathrm{id}$, the new output equals the old output, so $D'$
computes $F$. If $u=\lnot$, the new output is the negation of
the old output, so $D'$ computes $\lnot F$. If $u$ is constant,
or if $a$ is a constant node with $u$ equal to $\mathrm{id}$ or
$\lnot$, then the old output would be constant, contradicting
the assumption that $F$ is nonconstant. Hence this case cannot
occur for a nonconstant $F$, and $D'$ computes $F$ or $\lnot F$
with exactly one fewer gate. ∎

Two consequences are used. First, in a circuit computing a
nonconstant $F$, repeatedly applying Lemma 1 strictly decreases
the gate count at each step, so the process terminates. Second,
and more directly, an optimal circuit cannot contain an
eliminable gate, as shown in Corollary 1.

Corollary 1. Let $n \ge 3$ and let $F$ be $f_n$ or $\lnot f_n$.
Let $C$ be an optimal circuit for $F$, of minimum size
$s=C_{B2}(F)$. Then $C$ contains no dead gate and no eliminable
gate.

Proof. A dead gate is a gate that is not an ancestor of the
output gate (for $n \ge 3$ the output source is a gate, by
Section 1). Deleting a dead gate changes no gate value that
affects the output, so it gives a smaller circuit for the same
$F$, contradicting optimality. Hence $C$ has no dead gate.

Suppose $C$ contains an eliminable gate. Lemma 1 gives a circuit
with $s-1$ gates computing $F$ or $\lnot F$. If it computes $F$,
this contradicts optimality for $F$ directly. If it computes
$\lnot F$, convert it to a circuit for $F$ with the same number
$s-1$ of gates as follows: if its output source is a gate, flip
that gate's table (replace $\varphi$ by $1-\varphi$, which is
again one of the $16$ tables); if its output source is an input
or a constant, then $\lnot F$ would depend on at most one
variable, contradicting the dependence of $F$ (and hence of
$\lnot F$) on all $n \ge 3$ variables proved in Section 1. Either
way, optimality for $F$ is contradicted. Hence $C$ has no
eliminable gate. ∎

The same argument, with $F$ and $\lnot F$ interchanged, shows
that an optimal circuit for $\lnot F$ also has no dead or
eliminable gate, and that flipping the output gate's table gives
$C_{B2}(f_n)=C_{B2}(\lnot f_n)$ for $n \ge 3$: flipping the output
gate of an optimal circuit for one gives a same-size circuit for
the other, and vice versa. This symmetry is not needed as a
separate lemma in the induction, because the induction hypothesis
is stated for both $F$ and $\lnot F$ and the elimination steps
are allowed to flip the sign.

## 4. First gate and fanout

Fix $n \ge 3$ and $F \in \{f_n, \lnot f_n\}$, and let $C$ be an
optimal circuit for $F$. By Corollary 1, $C$ has no dead or
eliminable gate.

A first gate is a gate none of whose pin sources is a gate. Such
a gate exists: start from any gate and repeatedly follow a pin
source that is a gate; acyclicity and finiteness force this walk
to stop at a gate with no gate pin source. Let $g$ be a first
gate. Its pin sources are therefore inputs or constants. Since
$g$ is not eliminable, it cannot have a constant pin source,
cannot have both pins from the same source, and cannot have a
table depending on at most one pin value. Hence the two pin
sources of $g$ are two distinct input variables, call them
$x_i$ and $x_j$. The table of $g$ is otherwise unrestricted.

Lemma 2 (fanout). At least one of $x_i, x_j$ is a pin source of
a gate $h \ne g$.

Proof. Suppose neither $x_i$ nor $x_j$ feeds any gate other than
$g$. Let $R$ be the set of the remaining $n-2$ variables and let
$m=n-2 \ge 1$. Form a circuit $C_y$ by deleting $g$, introducing
a fresh input node $y$, and replacing every pin that was $g$ by
$y$, keeping all other gates and tables unchanged. We claim that
for every assignment to $(x_i,x_j,R)$, the output of $C$ equals
the output of $C_y$ on $(g(x_i,x_j),R)$. This is proved by
induction on a topological order of the gates other than $g$:
such a gate has no pin source equal to $x_i$ or $x_j$ by the
supposition, and its gate pin sources are earlier gates whose
values, by the induction hypothesis, depend on $(x_i,x_j)$ only
through $g$, or are $g$ itself. Hence every gate value other
than $g$'s depends on $(x_i,x_j)$ only through $g(x_i,x_j)$,
and the output does as well.

Now consider the three assignments $(x_i,x_j)=(0,0)$, $(0,1)$,
$(1,1)$. Their $g$-values lie in $\{0,1\}$, so two of the three
assignments, say $p \ne q$, share the same $g$-value $v$. Then
for every $R$, the output of $C$ on $(p,R)$ equals
$C_y(v,R)$ and also equals the output of $C$ on $(q,R)$. Hence
the two restriction functions of $F$ obtained by fixing
$(x_i,x_j)$ to $p$ and to $q$ are identical as functions of $R$.

But for $F=f_n$, the three restrictions are pairwise distinct.
They are: for $(0,0)$, the function $\mathrm{Th}_m^2(R)$, which
is $1$ exactly when $R$ contains at least two $1$s; for
$(0,1)$, the function $\mathrm{Th}_m^1(R)$, which is $1$ exactly
when $R$ contains at least one $1$; and for $(1,1)$, the
constant function $1$. These are pairwise distinct for
$m \ge 1$: $\mathrm{Th}_m^2$ and $\mathrm{Th}_m^1$ differ on an
$R$ with exactly one $1$ ($0$ versus $1$); $\mathrm{Th}_m^1$ and
$1$ differ on the all-$0$ $R$ ($0$ versus $1$); and
$\mathrm{Th}_m^2$ and $1$ also differ on the all-$0$ $R$.
For $F=\lnot f_n$, the three restrictions are the negations of
these three functions, which are still pairwise distinct,
because negation is injective: if $\lnot A = \lnot B$ as
functions, then $A=B$. This contradicts the identity of the
restrictions for $p$ and $q$. Therefore the supposition is
false, and at least one of $x_i, x_j$ feeds a gate $h \ne g$. ∎

## 5. Induction: two gates are eliminated at each step

Theorem. For every integer $n \ge 3$,
$C_{B2}(f_n) \ge 2n-4$ and $C_{B2}(\lnot f_n) \ge 2n-4$.

Proof by induction on $n$, proving both statements together at
each $n$.

Base case $n=3$. Let $F \in \{f_3, \lnot f_3\}$. A circuit of
size $0$ has an input or constant output source, depending on at
most one variable, contradicting dependence of $F$ on all three
variables. A circuit of size $1$ has a single gate whose pins are
inputs or constants, so its output depends on at most two
distinct input variables, again contradicting dependence on all
three variables. Hence every circuit for $F$ has size at least
$2 = 2\cdot 3 - 4$.

Inductive step. Let $n \ge 4$ and assume the theorem for $n-1$,
for both $f_{n-1}$ and $\lnot f_{n-1}$. Let
$F \in \{f_n, \lnot f_n\}$ and let $C$ be an optimal circuit for
$F$, of size $s$. We show $s \ge 2n-4$.

By Corollary 1, $C$ has no dead or eliminable gate. Let $g$ be a
first gate of $C$, with distinct input pin sources $x_i$ and
$x_j$, as in Section 4. By Lemma 2, at least one of $x_i, x_j$
feeds a gate $h \ne g$. Interchange the names of $x_i$ and $x_j$
if necessary so that $x_i$ feeds $h$; $x_i$ is the variable that
will be set to $0$.

Substitute $x_i=0$. Let $D$ be the circuit obtained from $C$ by
replacing the input node $x_i$ by a constant $0$ node, keeping
every gate and every table unchanged. $D$ still has $s$ gates.
Its function $F'$ on the remaining $n-1$ variables is the
restriction of $F$ to $x_i=0$: since $x_i$ contributes $0$ to the
count, $F'$ is $f_{n-1}$ if $F=f_n$, and $\lnot f_{n-1}$ if
$F=\lnot f_n$. For $n-1 \ge 3$, $F'$ is nonconstant (all-$0$
gives $0$, all-$1$ gives $1$, for $f_{n-1}$; reversed for its
negation).

In $D$, the gate $g$ is eliminable: its pins are the constant
$0$ and $x_j$, so its output is a unary function of $x_j$.
Also, the gate $h$ is eliminable in $D$: one of its pins is the
constant $0$, so its output is a unary function of its other
pin. Note that $h \ne g$, and that $D$ is a circuit computing
the nonconstant function $F'$, so Lemma 1 applies to $D$.

Apply Lemma 1 to $g$ in $D$. This gives a circuit $D_1$ with
$s-1$ gates, computing $F'$ or $\lnot F'$. The construction in
Lemma 1 deletes $g$ and leaves every other gate in place, with
successor tables composed as described; in particular, $h$ is
still a gate of $D_1$, even if it has become dead after $g$ was
deleted, because $D_1$ is defined to keep all gates other than
$g$.

We claim $h$ is still eliminable in $D_1$. Its one pin is still
the constant $0$ node: elimination of $g$ never replaces a
constant pin, because a constant node is not a gate. Its other
pin is either unchanged, or, if it was $g$, it is now the node
$a$ from the elimination of $g$, which is $x_j$ or a constant
node (since $g$'s pins in $D$ are $0$ and $x_j$). In every case,
$h$'s output in $D_1$ is a unary function of a single node (or
is constant), so $h$ is eliminable in $D_1$.

Apply Lemma 1 to $h$ in $D_1$. This gives a circuit $D_2$ with
$s-2$ gates, computing $F'$ or $\lnot F'$: each elimination step
preserves the function or flips it to its negation, and two such
steps still leave $F'$ or $\lnot F'$. The constant-output case
of Lemma 1 cannot occur here, because $F'$ and $\lnot F'$ are
nonconstant for $n-1 \ge 3$.

By the induction hypothesis for $n-1$, every circuit for
$f_{n-1}$ or $\lnot f_{n-1}$ has size at least $2(n-1)-4$.
$D_2$ is such a circuit, so $s-2 \ge 2(n-1)-4$, and therefore
$s \ge 2n-4$. Since $F$ was arbitrary in $\{f_n, \lnot f_n\}$,
both statements hold at $n$. This completes the induction. ∎

For $n \le 2$, the claimed inequality $C_{B2}(f_n) \ge 2n-4$
holds trivially because $C_{B2}(f_n) \ge 0$ and $2n-4 \le 0$.
Combined with the theorem for $n \ge 3$,

$$
C_{B2}(f_n) \ge 2n-4 \qquad\text{for every integer } n \ge 0,
$$

which is the unscored string claim in submission.json.

## 6. From the linear bound to the declared coefficient and exponent

For $n \ge 3$, $C_{B2}(f_n)/n \ge (2n-4)/n = 2-4/n$, which tends
to $2$ as $n\to\infty$ over all integers $n$. Hence

$$
\liminf_{n\to\infty} \frac{C_{B2}(f_n)}{n} \ge 2.
$$

In the epsilon form required by rules.md, let $\varepsilon>0$
and set $N_\varepsilon=\max(3,\lceil 4/\varepsilon\rceil)$.
For every integer $n \ge N_\varepsilon$, $C_{B2}(f_n) \ge 2n-4$
and $\varepsilon n \ge 4$ (since $n \ge 4/\varepsilon$), so
$2n-4 \ge (2-\varepsilon)n$. Thus
$C_{B2}(f_n) \ge (2-\varepsilon)n$ for every integer
$n \ge N_\varepsilon$. With $g(n)=2\cdot n^1\cdot
(\log_2 n)^0\cdot(\log_2\log_2 n)^0$, the two logarithmic
factors are $1$ for every $n \ge 2$ where they are defined, so
$g(n)=2n$ here, and $\liminf C_{B2}(f_n)/g(n) \ge 1$ follows:
indeed, for every $\delta>0$, $C_{B2}(f_n) \ge (1-\delta)2n$
whenever $n \ge \max(3,\lceil 2/\delta\rceil)$, because
$2n-4 \ge 2n-2\delta n$ is equivalent to $\delta n \ge 2$.
The finite projection $(c,p)=(2,1)$ therefore follows. The
quantifiers cover every sufficiently large integer length, not
a sparse subsequence: the bound $2n-4$ holds at every $n \ge 3$,
and the small lengths $n=0,1,2$ are defined and handled in
Sections 1 and 5.

No remainder subtlety arises. The loss relative to $2n$ is
exactly $4$, and $4/n\to 0$; it is not a positive fraction of
$n$. No direct-sum additivity, potential-to-gate identification,
sparse-length, limsup, or finite-experiment argument is used.
No monotone, formula, or other restricted-basis bound is used,
and no transfer from such a bound is needed, because every step
is in full B2.

## 7. Checklist against the challenge rules

- One fixed family: threshold $2$, defined once in Section 1,
  for all lengths, with $f_0=f_1=0$ and $f_2$ the AND of two bits.
- NP eligibility: fixed polynomial-time verifier $V$ and witness
  bound $W(n)=0$ in Section 1, with both directions proved. The
  family is in fact in P, so NP membership does not depend on
  any unproved construction.
- Model: full B2, direct argument, all $16$ truth tables,
  unrestricted sharing and depth, exact computation, nonuniform
  circuits allowed, inputs and wires uncharged, each retained
  gate counted once.
- Constant, unary, and terminal-output cases: constant sources
  cost zero and are vertices, not gates; a gate with a constant
  pin, repeated pins, or a table depending on at most one pin is
  eliminable and is handled by Lemma 1, never inserted as a free
  unary gate; a correct circuit for $n \ge 3$ cannot have an
  input or constant as its output source, so the output gate
  exists, and flipping an output gate's table is used only to
  convert between $F$ and $\lnot F$ at the same gate count.
- Bit accounting: $n$ is the actual input bit length. There are
  no blocks, fields, encodings, or padding, so there is no linear
  loss to account for.
- Joint theorem: the coefficient and the exponent are proved
  for the same family, the same verifier, and the same model.
  No separate lemma supplies one coordinate and a different
  lemma the other.
- What is not claimed: no coefficient above $2$, no exponent
  above $1$, no lane or global endpoint, and no $31/10$ result.
  This submission does not reproduce or rely on the Li-Yang/BK
  baseline construction.
