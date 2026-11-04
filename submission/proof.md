# Threshold-two: a self-contained $2n-4$ lower bound in the full B2 model

This proof is self-contained. No external theorem is used as a premise.
Every step needed for the scored claim is proved below, in the model
and with the quantifiers stated in `rules.md`.

## 1. The family, and why it is in NP

Let $L$ be the set of binary strings containing at least two $1$s.
For $n\ge 0$ let $f_n:\{0,1\}^n\to\{0,1\}$ be the characteristic
function of $L$ at length $n$; thus

- $f_0$ on the empty string is $0$;
- for $n\ge 1$, $f_n(x)=1$ if and only if at least two coordinates
  of $x$ are $1$.

In particular $f_1$ is identically $0$ and $f_2$ is AND.
The family is one fixed total family, defined at every length,
including the small exceptions, by the same rule.

One deterministic algorithm evaluates every $f_n$: scan the
$n$ input bits, count the $1$s, and output $1$ exactly when the
count is at least $2$. This runs in time linear in $n=|x|$,
with no padding, no advice, and no oracle. It is a single
polynomial-time evaluator for all lengths.

Define the verifier $V$ by $V(x,w)=1$ if and only if $w$ is the
empty string and the evaluator above accepts $x$; otherwise
$V(x,w)=0$. Take the witness bound $W(n)=0$ for every $n$.
Then $x\in L$ if and only if there is $w$ with $|w|\le W(|x|)$
and $V(x,w)=1$: if $x\in L$, the empty witness works; if the
empty witness works, $x\in L$ by the definition of $V$; and
no non-empty witness works, because $V$ rejects every
non-empty $w$. Hence $L$ is in P, and therefore in NP.
The witness is not an input to the circuits below, and the
circuit input length is the actual length $n=|x|$.

For $n\ge 3$, $f_n$ depends essentially on every variable.
Indeed, fix $i$, set exactly one other coordinate to $1$ and
all remaining coordinates to $0$. Then $f_n$ is $0$ when
$x_i=0$ (the total number of $1$s is $1$) and $1$ when
$x_i=1$ (the total is $2$). So flipping $x_i$ flips the
output for this assignment to the other variables. In
particular $f_n$ is not constant and not a function of a
single variable (or its negation) for $n\ge 3$.

For $m\ge 0$ and $k\in\{0,1,2\}$ write $\mathrm{Thr}_k^m$ for
the function on $m$ bits that is $1$ exactly when at least
$k$ bits are $1$. Thus $\mathrm{Thr}_2^m=f_m$,
$\mathrm{Thr}_1^m$ is OR, and $\mathrm{Thr}_0^m$ is the
constant $1$ function. Fixing one coordinate of $f_n$ to $0$
gives $\mathrm{Thr}_2^{\,n-1}$ on the remaining coordinates,
because the total number of $1$s is then the number of $1$s
among the remaining coordinates.

## 2. The circuit model

A circuit in the model of `rules.md` is a finite directed
acyclic graph. Its nodes are input nodes $x_1,\dots,x_n$,
constant nodes with fixed values $0$ or $1$ (cost zero),
and gate nodes. Each gate node has exactly two pins, each
pin pointing to a node earlier in a topological order, and
a truth table, one of the $16$ binary Boolean functions.
The value of a gate is its truth table applied to the
values at its two pins. One node is designated the output;
choosing the output costs zero. The size of a circuit is
the number of its gate nodes; input nodes, constant nodes,
wires, and output selection cost zero. Fan-out and depth
are unrestricted, subgraphs may be shared, and circuits
may be nonuniform: the argument below works for an
arbitrary single circuit at each length $n$, with no
uniformity assumption on circuits.

A gate is *simplifiable* if its output, as a function of
the values at its two pins, is constant, or is equal to
one pin value, or is equal to the negation of one pin
value, for every assignment to the pin values. This
covers every way a gate can fail to use both pins
essentially: a degenerate truth table (one of
$0,1,u,\bar u,v,\bar v$ as a function of pin values
$u,v$); two pins pointing to the same node, in which
case the output is a unary function of that node; and a
pin pointing to a constant node, in which case the output
is a unary function of the other pin. Conversely, a gate
that is not simplifiable has two distinct pin nodes,
neither constant, and a truth table depending
essentially on both pins.

Before simplifying a circuit, trim it: delete every gate
that is not an ancestor of the output node. Deleted gates
do not affect the output value, so the computed function
is unchanged and the size does not increase. After
trimming, every remaining gate is an ancestor of the
output, and the output node, if it is a gate, has no
successor gate. Normalisation repeats this trimming
before each removal; see Lemma 3.2.

## 3. Removing one simplifiable gate

Lemma 3.1. Let $C$ be a trimmed circuit computing a
function that depends on at least two variables. Let $G$
be a simplifiable gate of $C$. Then there is a circuit
$C'$ computing the same function, with exactly one fewer
gate than $C$, obtained from $C$ by removing $G$ only
(and rewiring its successors as described below), with
no gate added and no uncharged computational basis
introduced.

Proof. By simplifiability, the value of $G$ is, for every
assignment, a constant $c$, or a pin value $w$, or
$\bar w$, where $w$ is the value at one of its pins.
The pins of $G$ point to earlier nodes, so $w$ is not
$G$ itself.

If $G$ is not the output node, remove $G$ and rewire
every successor gate $H$ that reads $G$ as follows.
If $G$ is constant $c$, replace each pin of $H$ that
pointed to $G$ by a constant node $c$. If $G$ equals
$w$, replace each such pin by the node that supplied
$w$. If $G$ equals $\bar w$, replace each such pin by
the node that supplied $w$ and negate the corresponding
argument of $H$'s truth table: replace the table
$t$ by the table $t'$ with $t'(\dots,\bar w,\dots)=
t(\dots,w,\dots)$ in that argument (and similarly in
both arguments if $H$ read $G$ twice). The table $t'$
is still one of the $16$ binary truth tables, so $H$
remains a B2 gate and is still counted once; no new
gate is created, and negation is absorbed into an
already-counted gate, not supplied as a free gate.
Pointwise, $H$ then computes the same value as before,
because the value supplied to its (possibly negated)
pin is exactly the value $G$ had. All other gates are
unchanged, and the output node is unchanged. By
induction on topological order, every node computes
the same value as before, so $C'$ computes the same
function, with one fewer gate.

If $G$ is the output node, then, because $C$ is
trimmed, $G$ has no successor gate. Its value is the circuit output.
If its value is constant, or equals an input node or
its negation, or equals a constant node, then the
computed function is constant or a function of one
variable, contradicting the hypothesis that it depends
on at least two variables. Hence this case cannot
occur unless the pin node $w$ is a gate node. In that
remaining case, remove $G$ and designate $w$ as the
new output node (output selection costs zero). If
$G$ equals $w$, designate $w$ unchanged. If $G$
equals $\bar w$, replace $w$'s truth table by its
negation, so $w$ now computes $\bar w$, designate
$w$ as output, and compensate every other successor
gate of $w$ (other than the removed $G$) by negating
each argument of its truth table corresponding to a
pin that reads $w$ (both arguments, if it reads $w$
twice), so that each such successor still receives
the original value of $w$; negated tables remain among
the $16$ B2 tables and no gate is added. Thus $C'$
computes the same function with one fewer gate.
This covers the constant, unary, and terminal-output
cases. ∎

Lemma 3.2 (normalisation). Let $C$ compute a function
depending on at least two variables. The following
trim-then-remove loop terminates in a circuit $C^*$
computing the same function, which is trimmed and
has no simplifiable gate, and with at most as many
gates as $C$: at each step, first trim the current
circuit to the ancestors of its output node, then,
if the trimmed circuit has a simplifiable gate,
apply Lemma 3.1 to it; repeat. Trimming deletes
only non-ancestor gates, adds no gate, rewires no
pin, and leaves the computed function unchanged.

Proof. Trimming never increases the gate count, and
each application of Lemma 3.1 to a trimmed circuit
removes exactly one gate and adds none, so each
loop iteration that removes a gate strictly decreases
the gate count; the process must therefore stop
after finitely many steps. It stops only at a
trimmed circuit with no simplifiable gate, since
otherwise another iteration would apply. The
function is preserved by every trim and every
Lemma 3.1 step. The gates of $C^*$ are a subset of
the gates of $C$, because no step adds a gate. ∎

In a normalised circuit, therefore, every gate has two
distinct pin nodes, neither a constant node, and a
truth table depending essentially on both pins. In
particular no gate reads the same node twice, no gate
reads a constant node, and the output node is a gate.
For $f_n$ with $n\ge 3$, every input variable feeds
at least one gate in a normalised circuit: otherwise
the computed function would be independent of that
variable, contradicting essential dependence in
Section 1.

## 4. Fixing a variable that feeds two gates

Lemma 4.1. Let $n\ge 4$, let $C$ be a normalised
circuit computing $f_n$, and let $x_i$ feed at least
two distinct gates of $C$. Fix $x_i=0$ and normalise
the resulting circuit. The result computes
$f_{\,n-1}$ on the remaining $n-1$ variables and has
at most $|C|-2$ gates, where $|C|$ is the number of
gates of $C$.

Proof. Replace the input node $x_i$ by a constant
node $0$, leaving every gate unchanged. Call this
circuit $C_0$. For the remaining variables $y$,
$C_0$ outputs $f_n$ with $x_i=0$, which is
$\mathrm{Thr}_2^{\,n-1}(y)=f_{\,n-1}(y)$, by
Section 1. The function $f_{\,n-1}$ depends on at
least two variables, since $n-1\ge 3$, so Lemmas
3.1 and 3.2 apply at every step of normalising $C_0$.

Let $S$ be the set of gates of $C$ that read $x_i$;
$|S|\ge 2$ by hypothesis. In $C_0$ each gate in $S$
reads the constant node $0$ on one pin, so each is
simplifiable. Normalisation by Lemma 3.2 proceeds
by trimming and Lemma 3.1 removals, adds no gate,
and its final circuit's gates are a subset of the
gates of $C_0$. We claim no gate in $S$ is present
in that final circuit. Suppose a gate $G\in S$ were
present there. Its pin that read $x_i$ still points
to a constant node: trimming rewires no pin, and a
pin is replaced in a Lemma 3.1 step only when it
points to a removed gate, and a constant node is not
a gate, so that pin is never replaced. Hence $G$
still reads a constant on one pin in the final
circuit, so $G$ is still simplifiable, contradicting
the fact that the final circuit, by Lemma 3.2, is
trimmed and has no simplifiable gate. The other pin
of $G$ may have been rewired by earlier removals,
but that does not remove the constant pin, so the
argument is unaffected; a gate in $S$ may also have
been deleted by trimming, which likewise counts as
removed. Thus all of $S$, at least two distinct
gates, are absent from the final circuit. Hence
the final circuit has at most $|C|-2$ gates
($|C_0|=|C|$) and, by Lemma 3.1 applied at each
step and trimming preserving the function,
computes the same function as $C_0$, namely
$f_{\,n-1}$. ∎

## 5. A variable in the first gate feeds another gate

Lemma 5.1. Let $n\ge 4$ and let $C$ be a normalised
circuit computing $f_n$. Let $G$ be a gate of $C$
with no gate ancestor (a topologically first gate).
Then the pins of $G$ are two distinct input variables
$x_i$ and $x_j$, and at least one of $x_i,x_j$ feeds
a gate distinct from $G$.

Proof. Such a gate $G$ exists: $C$ has at least one
gate, since $f_n$ is non-constant, and a gate of
minimal topological rank has no gate ancestor. Its
pins cannot point to gates, by minimality, and by
normalisation they are distinct non-constant nodes,
so they are distinct input variables $x_i,x_j$.

Suppose, for contradiction, that neither $x_i$ nor
$x_j$ feeds any gate other than $G$. Let
$R$ be the set of the remaining $n-2$ variables and
let $g$ be the value computed by $G$, as a function
of $(x_i,x_j)$.

Claim: every gate $H\ne G$ of $C$, and the output,
computes a function depending only on $R$ and on $g$.
Indeed, prove this by induction on topological order.
For a gate with no gate predecessors, its pins are
input variables, none of which is $x_i$ or $x_j$ by
the supposition (except for $G$ itself), so its value
depends only on $R$. If the output gate is $G$
itself, its value is $g$. In general, a gate
$H\ne G$ cannot read $x_i$ or $x_j$ directly,
by the supposition. Its pins point to $R$ variables
or to earlier gates. Earlier gates are $G$ itself,
whose value is $g$, or gates already covered by the
induction hypothesis, whose values depend only on
$R$ and $g$. Hence the value of $H$, obtained by
applying its truth table to such values, depends
only on $R$ and $g$. The output node is a gate, so
the same holds for the circuit output.

Consequently, if $(a,b)$ and $(a',b')$ are assignments
to $(x_i,x_j)$ with $G(a,b)=G(a',b')$, then the
restrictions of $f_n$ to $R$ obtained by fixing
$(x_i,x_j)$ to those assignments are equal as
functions on $R$.

But the three assignments $(0,0)$, $(1,0)$, $(1,1)$
give three pairwise distinct restrictions, since
$n-2\ge 2$:

- fixing $(0,0)$ gives $\mathrm{Thr}_2^{\,n-2}$;
- fixing $(1,0)$ gives $\mathrm{Thr}_1^{\,n-2}$ (OR);
- fixing $(1,1)$ gives $\mathrm{Thr}_0^{\,n-2}$
  (constant $1$).

They are distinct: on $R$ with exactly one $1$,
$\mathrm{Thr}_2^{\,n-2}=0$ while
$\mathrm{Thr}_1^{\,n-2}=1$; on $R$ all $0$,
$\mathrm{Thr}_1^{\,n-2}=0$ while
$\mathrm{Thr}_0^{\,n-2}=1$; and on $R$ all $0$,
$\mathrm{Thr}_2^{\,n-2}=0$ while
$\mathrm{Thr}_0^{\,n-2}=1$.

Since $G$ takes only the two values $0$ and $1$,
two of the three assignments share the same $G$
value, so by the claim their restrictions must be
equal, contradicting distinctness. Hence the
supposition is false: $x_i$ or $x_j$ feeds another
gate. ∎

## 6. Induction: $C_{B2}(f_n)\ge 2n-4$ for $n\ge 3$

We prove, by induction on $n$, that every B2 circuit
computing $f_n$ has at least $2n-4$ gates, for every
integer $n\ge 3$.

Base $n=3$: the bound is $2$. Suppose a circuit
computing $f_3$ had at most one gate. Normalising
it (Lemma 3.2) gives a circuit with at most one gate
computing the same function. With $0$ gates the
output is an input node or a constant node, so the
function depends on at most one variable. With $1$
gate, by normalisation neither pin is a constant
node and the pins are distinct, so that gate's pins
are distinct input nodes (there is no other gate),
at most two distinct variables, so the function
depends on at most two variables. Both
contradict essential dependence of $f_3$ on all
three variables (Section 1). Hence every circuit
for $f_3$ has at least $2$ gates.

Inductive step: let $n\ge 4$ and let $C$ be any
circuit computing $f_n$. Normalise $C$ by Lemma 3.2
to a normalised circuit $C^*$; $C^*$ computes $f_n$
and $|C^*|\le |C|$. Apply Lemma 5.1 to $C^*$ to get
a variable $x_i$ feeding at least two distinct gates
of $C^*$ (the variable supplied by Lemma 5.1,
together with the first gate $G$). Apply Lemma 4.1
to $C^*$, fixing $x_i=0$: we obtain a circuit for
$f_{\,n-1}$ with at most $|C^*|-2$ gates. By the
induction hypothesis applied to $n-1\ge 3$, that
circuit has at least $2(n-1)-4$ gates. Hence
$|C^*|-2\ge 2(n-1)-4$, i.e. $|C^*|\ge 2n-4$, and
therefore $|C|\ge |C^*|\ge 2n-4$. This completes
the induction.

Thus, for every integer $n\ge 3$,
$C_{B2}(f_n)\ge 2n-4$, where $C_{B2}$ counts gates
in the full B2 model of Section 2.

## 7. The declared guarantee and its quantifiers

Let $g(n)=2\cdot n^1$, with log and loglog exponents
$0$. For $n\ge 3$, $(2n-4)/g(n)=1-2/n$, which tends
to $1$ as $n\to\infty$. Hence
$\liminf_{n\to\infty} C_{B2}(f_n)/g(n)\ge 1$.
Equivalently, for every $\epsilon>0$, taking
$N_\epsilon=\max(3,\lceil 4/\epsilon\rceil)$ gives,
for every integer $n\ge N_\epsilon$,
$C_{B2}(f_n)\ge 2n-4\ge (2-\epsilon)n$,
since $\epsilon n\ge 4$. The loss $4$ is $o(n)$,
not a positive fraction of $n$; the coefficient is
exactly $2$, not $2$ reduced by a linear remainder.
The bound is proved for every sufficiently large
integer length, not for a sparse subsequence, a
limsup, or finitely many lengths: the induction in
Section 6 covers every integer $n\ge 3$, for the
one fixed family of Section 1, with the same
evaluator, verifier, and witness bound throughout.

The finite projection is therefore $(c,p)=(2,1)$:
coefficient $2/1$, exponent $1/1$, with log
exponent $0/1$ and loglog exponent $0/1$ in the
full guarantee. No stronger rate is claimed.

## 8. Model transfer and accounting

The circuits counted above are exactly the circuits
of `rules.md`: deterministic, acyclic, one output,
two input pins per gate, all $16$ binary truth
tables, unrestricted depth and fan-out, exact
computation on every input, nonuniform circuits
allowed. Normalisation and substitution never leave
this model: successor truth tables stay among the
$16$ tables, negation is absorbed into an
already-counted successor gate (or into the output
gate in the terminal case of Lemma 3.1), no gate
in the final normalised circuit reads a constant
node (unused constant nodes, if any, cost zero
and do not affect the gate count), and no gate is
ever added without being counted. Repeated pins and shared
subgraphs are handled by normalisation: a gate
with repeated pins is simplifiable and is removed,
and sharing is irrelevant to the counting, which
counts distinct gates removed, never a potential
identified with a gate count. The first-gate
argument in Lemma 5.1 uses only the DAG structure
and the distinctness of the three restrictions; it
does not assume anything about depth, fan-out, or
uniformity. Witness bits are never circuit inputs.
All parameter dependencies are explicit: the only
constants are the base-case bound $2$ and the
loss $4$, and $N_\epsilon$ above depends only on
$\epsilon$.

No premise from outside this proof is used. The
argument is the classical gate-elimination pattern,
but every lemma it needs — essential dependence,
one-gate removal, normalisation, substitution,
the first-gate factoring argument, and induction —
is proved here from the definitions in Sections 1
and 2.
