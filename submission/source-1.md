# Source evidence 1. Spectral spaces and curvature

Excerpts from OpenAI's mathematics preprints at commit 92b516ef8e2abf66c3db7aa31438d3010bc7956f. The dependency map, bibliography and Apache 2.0 license are in source-5.md. The authors' result labels are preserved.

## G / prices.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/prices.tex, Git blob 408773173318cb2dc09b5d288d021c969b7fed76. The fixed-placement Zariski-closure lemma and proof are omitted.

~~~latex
\section{The calculus of prices}\label{sec:spectra}
Every spectral point on a fixed set of parties gives a pseudometric on
those parties. The distance between two parties is the logarithm of the
spectral point's value on their binary wire. We first show how the same invariant evaluates tensors on subsets of the
parties and controls the cost of moving their indices. All logarithms are to
base two unless indicated otherwise.

\subsection{Spectral points and prices}
Fix a field $k$ and a finite set $J$ of at least two parties. A tensor has
one finite-dimensional vector space at each party. We use bases throughout,
and contraction is bilinear. Write $X\le Y$ when local linear maps, one at
each party, take $Y$ to $X$. Such a transformation is a \emph{restriction}.
Two operations combine tensors on $J$. In a tensor product we group the
corresponding factors at each party. In a direct sum $X\oplus Y$ we take the
direct sum of the local spaces at every party, so the two summands occupy
separate blocks. Up to mutual restriction, these operations make the tensors
a commutative semiring. Its integer $r$ is the diagonal tensor
\[
                 \Delta_J(r)=\sum_{s=1}^r\bigotimes_{j\in J}e_s.
\]
A nonzero tensor $T$ of rank $R(T)$ satisfies $1\le T\le R(T)$.

A \emph{spectral point} $F$ assigns a nonnegative real number to each tensor.
It preserves $0$ and $1$, is monotone under restriction, and satisfies
\[
       F(X\otimes Y)=F(X)F(Y),\qquad F(X\oplus Y)=F(X)+F(Y).
\]
In particular $F(\Delta_J(r))=r$. Matrix rank after flattening across a
fixed nontrivial cut of $J$ is one example. We call $H(T)=\log_2 F(T)$
the \emph{price} of a nonzero tensor. Tensor products then add prices,
and local maps can only decrease them. The spectral points form a compact
space in the topology of pointwise convergence. Indeed, in standard
coordinate spaces their axioms define a closed subset of the product of
the intervals $[0,R(T)]$.

\subsection{Placement}
For a nonempty $S\subseteq J$, let $\Delta_S(q)$ be the diagonal tensor
sharing a common $q$-valued label among the parties in $S$, with
one-dimensional spaces elsewhere. A singleton diagonal is a nonzero local
vector. For $a\ne b$, the two-party diagonal $E_{ab}(q)=\Delta_{\{a,b\}}(q)$
is a wire of dimension $q$. A wire whose two ends lie at the same party
is also a nonzero local vector. Set
\[
     \delta(S)=H(\Delta_S(2)),\qquad d(a,b)=H(E_{ab}(2)).
\]
Tensoring diagonals multiplies their dimensions, and restriction can
reduce a dimension to any smaller value. Comparing powers of $q$ with
powers of two therefore gives $H(\Delta_S(q))=\delta(S)\log_2 q$.
For a flattening rank, $d(a,b)$ is one when the cut separates $a$ and $b$
and zero otherwise.

For every spectral point, $d$ is a pseudometric. The triangle inequality
comes from joining two wires by contraction at their common endpoint.
A diagonal restricts to a wire between any two of its parties. Conversely,
suppose a tree joins the sites of $S$, possibly through auxiliary vertices.
At each vertex a local map equates its incident wire labels, retaining the
common label at vertices in $S$. Applying these maps to the product of
binary wires produces $\Delta_S(2)$. Hence
\begin{equation}\label{eq:span-price}
 \operatorname{diam}(S)\ \le\ \delta(S)\ \le\
 \sum_{ab\text{ in the tree}} d(a,b).
\end{equation}

A \emph{placement} $x=(x_1,\ldots,x_r)$ assigns the indices of a nonzero
tensor $T$ to parties. We group the indices at each party and write $H(T;x)$ for the
resulting price. Changing a placement keeps the coefficient array of $T$
fixed. To move its $i$th index from $x_i$ to $y_i$, supply a wire of the
index dimension $q_i$ and contract one end with that index at $x_i$.
The same construction moves it back. Thus
\begin{equation}\label{eq:move}
 |H(T;x)-H(T;y)|\ \le\ \sum_{i=1}^r d(x_i,y_i)\log_2 q_i .
\end{equation}
Sites at distance zero may therefore be identified. A nonzero tensor
whose indices are all at one site has price zero.

A \emph{spectral space} is a metric space equipped with prices $H(T;x)$
for all finite placements over $k$. Inserting trivial slots leaves prices
unchanged, and a binary wire has price equal to the distance between its
ends. On each finite set $S$ of positive diameter, all prices must have
the form $\delta(S)\log_2 F_S(T)$ for one spectral point $F_S$ on $S$.
The scale $\delta(S)$ is the price of $\Delta_S(2)$, since
$F_S(\Delta_S(2))=2$. Multiplying all nonzero tensor prices and all
distances by the same positive constant gives another spectral space.
We set $H(0)=-\infty$. Nonzero tensors at a single site and nonzero
scalars on no legs have price zero.

\paragraph{Normalization on subsets.}
Fix a spectral point $F$ on $J$ and a subset $S$ of at least two parties.
Extend tensors on $S$ to $J$ by adjoining trivial slots, and denote the
resulting evaluation by $G$. We have $G(\Delta_S(q))=q^e$ for
$e=\delta(S)$. If $e=0$, the rank bound on $S$ gives $G(T)=1$ for
every nonzero $T$. If $e>0$, the evaluation $G^{1/e}$ is a spectral
point on $S$. Only additivity needs proof, since direct sum on $S$ and
direct sum on $J$ separate blocks at different sets of parties.

Group the terms of $(U\oplus V)^{\otimes n}$ by the number $j$ of
$U$-factors. After local changes of coordinates the $j$th block is
$U^{\otimes j}V^{\otimes(n-j)}\Delta_S(\binom nj)$. Each block is a
restriction of the full tensor. Conversely, the full tensor is the
ordinary sum of the blocks in the ambient spaces. Since an ordinary sum
is a restriction of a direct sum on $J$, $G(U\oplus V)^n$ lies between
the largest block value and $n+1$ times that value. The $1/e$ powers
of the block values are precisely the terms of the binomial expansion
of $(G(U)^{1/e}+G(V)^{1/e})^n$. Their sum lies between the largest
term and $n+1$ times that term. Taking $1/(en)$ powers in the bounds
for $G(U\oplus V)^n$ and letting $n\to\infty$ therefore gives
$G(U\oplus V)^{1/e}=G(U)^{1/e}+G(V)^{1/e}$.

The same normalization works for a list of parties in which several
indices occupy one site. Grouping the indices of a direct sum at that
site introduces additional zero blocks. Deleting those blocks by local
projection gives the direct sum of the grouped tensors, and local
inclusions recover the padded tensor. Equation~\eqref{eq:move} allows
indices at distance zero to move to a common representative without
changing their price. The resulting spectral space on the metric
quotient therefore retains all the original log prices.

\paragraph{Ordinary sums.}
If a product of tensors supplies $T$ by local maps and contractions,
then $H(T)$ is at most the sum of their prices. Suppose instead that
several constructions produce tensors $T_i$ with the same output
spaces on a finite set $S$ of positive diameter. Their auxiliary sites
may differ. Addition in the output spaces is a restriction of a direct
sum on $S$, so the normalized spectral point on $S$ gives
\begin{equation}\label{eq:price-ordinary-sum}
 H\!\left(\sum_i T_i\right)
 \le \delta(S)\log_2\sum_i 2^{H(T_i)/\delta(S)}
 \le \max_i H(T_i)+\delta(S)\log_2\#\{i\}.
\end{equation}
Here $2^{-\infty}=0$. At a single output site every nonzero sum has price
zero. These bounds also apply when terms cancel.


% Editorial omission of the unused fixed-placement Zariski-closure lemma and proof.
\subsection{Limits of spectral spaces}\label{sec:limits}
Two limiting operations preserve these data. For metric completion take
the limit of every fixed tensor price along convergent placements,
using \eqref{eq:move}. For a pointed metric ultralimit, possibly after
rescaling, take the ultralimit of prices of the representing placements.
For each fixed nonzero tensor these numbers are bounded by \eqref{eq:move}
because its indices are a bounded distance from the basepoint.
The same bound gives independence of representatives and agreement when
two parties merge. On a finite limit set of positive diameter, the
multipliers $\delta(S)$ are bounded above and away from zero by
\eqref{eq:span-price}. Dividing by them and passing to the limit preserves
all spectral point axioms. Thus the limits carry simultaneous values
for every tensor, including tensors over uncountable fields.


\subsection{Turning prices into restrictions}
Monotonicity gives a necessary condition for a tensor restriction.
A strict inequality at every spectral point has a converse for powers
of the tensors. We use the following form of Strassen's asymptotic
spectral theory \cite{Strassen1988} (see also \cite{Fritz,WZ}).
\begin{lemma}[Strict spectral comparison]\label{lem:spectral-comparison}
Suppose that $X,Y$ are nonzero and $F(X)<F(Y)$ for every spectral
point on $J$. There is an integer $K\ge 1$ such that
$K Y^m\ge X^m$ for every $m\ge 1$.
\end{lemma}
\begin{proof}
The proof separates an element from the cone of restrictions by a
spectral point, then clears the denominators in the resulting cone
inequality. Pass to the additive group completion of the tensor semiring
$S_J$ and localize at all
nonzero tensors, obtaining a commutative $\mathbb Q$-algebra
$A$. Nonzero tensors form a multiplicative set, and flattening rank
extends to this localization. Let $P$ be the cone of finite sums of
$(a-b)/s$ with $a\ge b$ and $s$ a nonzero tensor. This rational
convex cone is closed under multiplication by tensor fractions $t/s$.
Write $v\le_P w$ if $w-v\in P$. Every element of $A$ is between
integer bounds in this order, since $0\le_P t/s\le_P t$.
Evaluating at a flattening rank shows that $-1\notin P$.

We claim that $z\notin P$ implies the existence of a spectral point
whose extension to $A$ takes a nonpositive value at $z$. Indeed,
$P'=P-\mathbb Q_{\ge0}z$ also excludes $-1$, since $-1=p-qz$
with $p\in P$ would give $q>0$ and $z=(1+p)/q\in P$.
Extend the identity on $\mathbb Q1$ to a rational-linear functional
on $A$ positive on $P'$ (where positive means nonnegative on the
cone). For the extension by a new direction $h$, prescribe a value
between the supremum on old vectors bounded above by $h$ and the
infimum on old vectors bounded below by $h$, using the order of $P'$.
Integer bounds make the two sets nonempty, positivity orders their
values, and rational rescaling verifies the extension. Well-ordering $A$
gives the functional, with value at $z$ at most zero.

The rational-linear functionals positive on $P$ with value one at the
unit form a compact convex set, by the integer bounds. Minimize
evaluation at $z$, then (taking intersections at limit stages) minimize successively on a well-ordered list of
all coordinates in the preceding minimizing faces. The nested compact
faces have a unique common member $\ell$, which is extreme and has
$\ell(z)\le0$. For any nonzero tensor $s$, put $u=s/(1+s)$ and
$v=1/(1+s)$. Multiplication by either preserves $P$, and $u+v=1$.
If $s\le r$, both $u,v$ are at least $1/(r+1)$ in the cone order,
using $s\ge1$. Thus the functionals $\ell(u\,\cdot)$ and
$\ell(v\,\cdot)$ normalize to states and decompose $\ell$ with
positive coefficients. Extremality gives
$\ell(uh)=\ell(u)\ell(h)$ and $\ell(vh)=\ell(v)\ell(h)$.
As $vs=u$, it follows that $\ell(sh)=\ell(s)\ell(h)$.
The restriction of $\ell$ to tensors is the desired spectral point,
and the same identity determines its extension to fractions.

By compactness choose an integer $b$ with
$b F(Y)>(b+1)F(X)$ at every spectral point. The claim gives
$bY-(b+1)X\in P$. Collecting denominators and clearing the
localization gives $cbY-c(b+1)X=a-a'$ in the additive group
completion, for some nonzero tensor $c$ and tensors $a\ge a'$.
Equality in the completion allows addition of a common tensor in
$S_J$, so $cbY+w\ge c(b+1)X+w$ for some $w$ (absorb $a'$ into
$w$). Choose $q$ with $qcY\ge w$, possible since $cY\ge1$.
Trading $q$ times while retaining $w$ gives
\[
 q(b+1)cY\ \ge\ qbcY+w\ \ge\ q(b+1)cX+w\ \ge\ q(b+1)cX .
\]
Consequently $c'Y\ge c'X$ with $c'=q(b+1)c$. Reusing $c'$ gives
$c'Y^m\ge c'X^m$ for every $m$. Any integer $K\ge c'$ works.
\end{proof}


Strict comparison also lets us change the ground field without losing
prescribed spectral values (see \cite{Strassen1988,CVZ}).
\begin{lemma}\label{lem:comparison}
For $k\subseteq K$, the possible joint values of spectral points on
any family of tensors defined over $k$ are the same over the two fields.
\end{lemma}

\begin{proof}
Restricting spectral points over $K$ to tensors over $k$ gives a compact
subset of the spectrum over $k$. If this subset misses $F$, there are
finitely many tensor evaluations whose values at $F$ do not occur
jointly in the subset. Choose a rational polynomial in those evaluations
which is negative at $F$ and positive on the compact image. For example,
approximate the squared distance to the tuple at $F$, minus a small
positive constant. After clearing denominators, the polynomial is the
difference in evaluation of two nonzero tensors $Y,X$. Monomials
correspond to tensor products, and adding $1$ to both parts ensures
that $X$ and $Y$ are nonzero.

By Lemma~\ref{lem:spectral-comparison} over $K$ there is a fixed integer
$C$ such that $X^m\le C Y^m$ for all $m$. Fix $m$. The finitely many
coefficients of the local maps generate a $k$-algebra. Taking its
quotient by a maximal ideal specializes the restriction to a finite
extension $L/k$ by Zariski's lemma. The $n$th power of this restriction
is still given by one linear map per party over $L$. Expand each map
in a $k$-basis of $L$ and apply a $k$-linear functional $L\to k$ that
takes $1$ to $1$. This expresses $X^{mn}$ over $k$ as an ordinary sum
of at most $[L:k]^{|J|}$ scaled restrictions of $(CY^m)^n$. Hence
\[
              F(X)^{mn}\le [L:k]^{|J|} C^nF(Y)^{mn}.
\]
Taking $n$th roots and letting $n\to\infty$ gives
$F(X)^m\le C F(Y)^m$. Since $C$ is independent of $m$, this contradicts
$F(X)>F(Y)$.
\end{proof}
\section{Geodesic extension}
\label{sec:geodesic}
Can one insert a point between two existing sites while preserving the
prices of all tensors on the old sites? The next theorem gives such an
extension for every division of the distance between the sites. Its
second assertion adds a site at the center of a star and identifies the
least possible weighted sum of the new wire lengths. We then use an
exterior-algebra identity to control tensor prices as one or more
indices move between sites.
\subsection{Extending a spectral point by one party}

For an additional party $p$, let $\iota$ embed old tensors by the
trivial leg at $p$. We say $F'$ \emph{extends} $F$ if
$F'(\iota T)=F(T)$ for every old tensor. This condition fixes the
old spectrum values, including those of tensors with many
nontrivial legs.

\begin{theorem}[Exact extension]\label{thm:extension}
Let $F$ be a spectral point on $J$ over $k$, and write
$d_F(a,b)=\log_2 F(E_{ab}(2))$.
\begin{enumerate}
\item For $a\ne b$ in $J$ and any nonnegative $\lambda_a,\lambda_b$ with
$\lambda_a+\lambda_b=d_F(a,b)$ there is an extension by $p$ such
that $d_{F'}(a,p)=\lambda_a$ and $d_{F'}(p,b)=\lambda_b$.
\item For integers $w_i\ge 0$ on $J$, put
\[
 g_F(w)=\sup_{h\ge1,\ T\ne0}
                  \frac{\log_2 F(T)}h ,
 \qquad T\in\bigotimes_{i\in J} k^{2^{hw_i}},
 \quad h\in\mathbb Z.
\]
Then
\begin{equation}\label{eq:spectral-star}
                 g_F(w)=\min_{F'\text{ extends }F}
                             \sum_{i\in J} w_i d_{F'}(p,i).
\end{equation}
For each width profile $w$ the minimum is attained by an extension
whose center $p$ may depend on $w$. Taking $w$ to be the indicator
of a subset gives the equal-dimension star.
\end{enumerate}
\end{theorem}
\begin{proof}
Use the embedded nonzero old tensors as generators with prescribed
targets $f(\iota T)=F(T)$. For the first assertion add generators
$E_{pa}(2),E_{pb}(2)$ with targets $2^{\lambda_a},2^{\lambda_b}$.
For the second add the single generator
$S_w=\prod_i E_{pi}(2^{w_i})$ with target $2^{g_F(w)}$.
The rank bound gives $0\le g_F(w)\le\sum_i w_i$.
Targets multiply on formal monomials in the selected generators.

Contract the new leg of such a monomial $M$ by a covector. The old
spectral value of the result is at most $f(M)$. The old factors pass
unchanged through the contraction. In the first case the remaining
tensor is a matrix of rank at most $2^{\min(u,v)}$, where $u,v$ are
the respective wire multiplicities. Its value is at most
$2^{d_F(a,b)\min(u,v)}\le 2^{u\lambda_a+v\lambda_b}$.
In the second case a covector on the center of $S_w^u$ gives exactly
an arbitrary tensor of dimensions $2^{uw_i}$, whence the bound.
This also shows $g_F(w)\le\sum_i w_i d_{F'}(p,i)$ for every
extending $F'$.

Let $C_m$ be the supremum of the old spectral values of contractions of
$M^m$ at $p$. In the first case with $u=v$, choosing the identity
matrix gives $C_m=f(M)^m$. In the star case write
$M=(\iota T_0)S_w^u$, combining the embedded old factors into $T_0$.
Fix $h\ge1$ and a nonzero tensor $T$ of local dimensions $2^{hw_i}$.
For every $m$, take $\lfloor um/h\rfloor$ copies of $T$ and pad the
remaining dimensions by a nonzero simple tensor. The result is a
contraction of $S_w^{um}$, so
\[
  F(T_0)^m F(T)^{\lfloor um/h\rfloor}\le C_m\le f(M)^m .
\]
Taking the lower limit of $m$th roots gives
$\liminf_m C_m^{1/m}\ge F(T_0)F(T)^{u/h}$. Taking the supremum over
$h,T$ gives $C_m^{1/m}\to f(M)$, including when $u=0$.

We have bounded the spectral value of every contraction by its monomial
target. The bound is sharp in exponential rate for a monomial with balanced
wire multiplicities in the first case and for every star monomial in the
second. We now use strict spectral comparison to show that the targets
are attained simultaneously.
Suppose no spectral point attains them all. Compactness then gives a finite list $G_1,\ldots,G_s$, including
the added generators, and $\delta>0$ such that every spectral point
$F'$ on $J\cup\{p\}$ has some ratio $F'(G_i)/f(G_i)$ or its inverse
at least $1+\delta$. Let $M=\prod_i G_i^e$ and let
$M_1,\ldots,M_{2s}$ be the variants obtained by changing exactly one
exponent either from $e$ to $0$ or from $e$ to $2e$. Take $e$ large
enough that $(1+\delta)^e>8s$. Choose positive integers $L,c_j$ by
rational approximation so that
\[
 \frac1{8s}<\frac{c_j f(M_j)}{L f(M)}<\frac1{4s}
 \qquad(1\le j\le 2s).
\]
For any $F'$, the violating index and the choice of sign give an
index $j$ with $c_jF'(M_j)>L F'(M)$. By
Lemma~\ref{lem:spectral-comparison} there is a fixed $K$ with
\[
 (L M)^m\ \le\ K\biggl(\bigoplus_{j=1}^{2s} c_jM_j\biggr)^m
 \qquad(m\ge1).
\]
On the left apply at $p$ the all-ones covector on the $L^m$-copy
factor, times a covector on $M^m$ whose contraction has old spectral
value greater than $C_m/2$. Such a covector exists since $0<C_m<\infty$. The
copy factor becomes the full $L^m$-fold diagonal on the old parties.
The equal multiplicities required in the first case hold for $M$.
The resulting old values therefore approach $L f(M)$ after taking
$m$th roots. On the right, pull the covector back through the restriction at
$p$ and restrict it to each direct-sum branch. On a branch the input
is a monomial $(\iota T_*)V$, with $T_*$ the product of its old
tensors. Its leg at $p$ is just that of $V$, so an arbitrary covector
$\xi$ there gives the old tensor $T_*\,\xi(V)$ before applying the
old-party maps. Those maps commute with contraction at $p$ and can
only decrease the value. The contraction bound above applies jointly
to $V$, using its total wire multiplicities. A star power gives a
single tensor of the prescribed local dimensions and the split wires
give a single matrix of bounded rank. Thus each branch is bounded by
its monomial target without a factorization assumption on $\xi$.
Summing over branches by ordinary subadditivity bounds the old value by
$K(\sum_j c_j f(M_j))^m$. Since $\sum_j c_j f(M_j)<L f(M)/2$,
this is impossible for large $m$. The spectral point matching the
targets proves both assertions.
\end{proof}

Theorem~\ref{thm:extension} extends to arbitrary spectral spaces.
Apply the finite theorem, with the appropriate scale, on each finite
subset containing the specified sites. For the star give weight zero
to any additional sites and use prices in place of log spectral
values. Take the pointed ultralimit from Section~\ref{sec:limits}
along an ultrafilter containing all tails of this directed family.
Take a specified site as basepoint, choosing one of positive weight
for a star. The prescribed distance or the bound on the star sum keeps
the new point at bounded distance from the basepoint. The limit therefore
contains the new point and preserves every old price. Repeated sites
in a star are handled by adding their weights. For zero weights
throughout, or a singleton space, any old site serves as center.
\subsection{Relocating one tensor index}

To interpolate the price of a tensor, we first move a single index.
On an auxiliary exterior algebra, the sum of two operator products
in opposite orders recovers the original tensor index. Related fermionic constructions occur in
the crossing identity \cite{Anthropic,Khoruzhii} used in
Section~\ref{sec:crossing}.
Fix a spectral point $F$ and write $d=d_F$. Consider a tensor with a
specified local factor $k^q$. Write $T_A$ for its coefficient array
with that factor placed at $A$ and every other factor in its original
position. Moving the factor from $B$ to $A$ is a restriction after
supplying a $q$-dimensional wire between $A$ and $B$. Consequently
\[
                         F(T_A)\le q^{d(A,B)}F(T_B).
\]
Two indices of a wire at the same party give a nonzero local vector,
equivalent by restriction to the unit, even if the characteristic
divides $q$. Thus the estimate holds also for repeated parties, and
moving a factor a distance zero preserves its value.

On the exterior algebra $\bigwedge k^Q$, for integers $Q\ge2$, write
$L_i$ for left exterior multiplication by the $i$th basis vector
and $D_j$ for contraction by the $j$th dual basis vector. Explicitly,
\[
 D_j(v_1\wedge\cdots\wedge v_b)
   =\sum_{a=1}^b(-1)^{a-1}(v_a)_j\,
      v_1\wedge\cdots\wedge\widehat{v_a}\wedge\cdots\wedge v_b .
\]
In particular $D_j(e_i\wedge v)=\delta_{ij}v-e_i\wedge D_jv$, or
\begin{equation}\label{eq:wire-exterior}
                       D_j L_i+L_i D_j=\delta_{ij}I .
\end{equation}
This identity holds over the integers in the ordered subset basis, so
also in characteristic two. On $\bigwedge^b k^Q$ write
$R=\binom Qb$ and $R_\pm=\binom Q{b\pm1}$. For each fixed
$s\in[-1,1]$ we may take $1\le b\le Q-1$ so that, as $Q\to\infty$,
\begin{equation}\label{eq:wire-binomial}
                   R_+/R=Q^{s+o(1)},\qquad R_-/R=Q^{-s+o(1)} .
\end{equation}
For $s>0$ use $b=\lfloor Q^{1-s}\rfloor$ for sufficiently large
$Q$, for $s<0$ take the complementary level for $-s$, and for $s=0$
use a middle level. The ratios are $(Q-b)/(b+1)$ and $b/(Q-b+1)$.

\begin{lemma}[One index]\label{lem:wire-one-index}
Let $0\le u\le1$ and suppose $d(A,C)=u\,d(A,B)$ and
$d(C,B)=(1-u)d(A,B)$. Every nonzero tensor with a marked index
satisfies $F(T_C)\le F(T_A)^{1-u}F(T_B)^u$.
\end{lemma}
\begin{proof}
Put $t=d(A,B)$. The index-moving estimate handles $t=0$ and
$u\in\{0,1\}$. A one-dimensional index is trivial, so assume $q\ge2$
and $0<u<1$. In the $n$th power the marked index has dimension $Q=q^n$.
Using the matrices in Equation~\eqref{eq:wire-exterior}, replace the
marked index $j$ by a fresh index $i$ of the same dimension, with row
and column indices $\beta,\alpha$ on $\bigwedge^b k^Q$ at $B,A$.
The two arrays have entries
\[
   \sum_j(D_jL_i)_{\beta\alpha}(T^n)_{\ldots j\ldots}
   \quad\text{and}\quad
   \sum_j(L_iD_j)_{\beta\alpha}(T^n)_{\ldots j\ldots}.
\]
Their sum, with $i$ at $C$, is $T_C^n E_{AB}(R)$.
The two arrays have local implementations with the following wire
dimensions. In both rows $\alpha$ starts at $A$, $\beta$ ends at $B$,
and the new index $i$ is at $C$.
\[
\begin{array}{ccccc}
 A && C && B\\[2pt]
 \alpha & \xrightarrow{\ R\ } & L_i & \xrightarrow{\ R_+\ } & D_j\\[2pt]
 D_j & \xrightarrow{\ R_-\ } & L_i & \xrightarrow{\ R\ } & \beta
\end{array}
\]
Contraction by $D_j$ uses the old index of $T_B^n$ in the first row
and that of $T_A^n$ in the second. Both implementations are local
restrictions even when sites coincide. Subadditivity, $d(A,C)=ut$ and $d(C,B)=(1-u)t$
give
\[
 F(T_C)^n\le F(T_B)^n(R_+/R)^{(1-u)t}
                    +F(T_A)^n(R_-/R)^{ut}.
\]
Set $s=t^{-1}\log_q(F(T_A)/F(T_B))$, which lies in $[-1,1]$ by the
index-moving estimate. Use Equation~\eqref{eq:wire-binomial} and take
$n$th roots. The two candidate limits for the maximum are
$F(T_B)q^{(1-u)ts}$ and $F(T_A)q^{-uts}$, both equal to
$F(T_A)^{1-u}F(T_B)^u$.
\end{proof}

\subsection{Moving several indices together}
Consider a binary wire whose ends both move from $a$ to $b$ and occupy
the same site at fraction $t$ of the journey. Its price there is zero.
Interpolating one end at a time compares its price with the average
for two independent endpoint choices, namely $2t(1-t)d(a,b)$. Choosing
the same endpoint for both ends gives the sharp bound zero. Such
synchronized choices are available for every tensor. More generally,
the endpoint choices can have any joint law with the prescribed marginals.
\begin{proposition}[Joint convexity of placement]\label{prop:joint-convexity}
For each marked index $i$ let $z_i$ have distances
$t_i d_i,(1-t_i)d_i$ to $a_i,b_i$, respectively, where
$d_i=d(a_i,b_i)$ and $0\le t_i\le1$. For any probability law $\mu$
on endpoint choices choosing $b_i$ with probability $t_i$
and $a_i$ otherwise, with $x_i$ the chosen position,
\[
                 H(T;z)\le \mathbb E_\mu H(T;x).
\]
Unmarked indices keep their positions.
\end{proposition}

\begin{proof}
The assertion is immediate for $T=0$ or if every move is constant.
Omit constant moves and one-dimensional marks, and normalize prices by
the scale on the finite set of all positions. Successive applications
of Lemma~\ref{lem:wire-one-index} give the bound when the endpoint
choices are independent. Under any coupling with the same marginals,
a function of the form $c_0+\sum_i f_i(x_i)$ has the same expectation.
We use the exterior identity to decompose $T^{\otimes n}$ into terms
whose bounds at the endpoints approach functions of this form.

Divide log prices by $n$. Each term has the original marked positions
and some auxiliary indices at fixed sites. We record an upper bound
$N(x)$ for its normalized price at each endpoint configuration $x$,
after subtracting the normalized prices of the auxiliary wires adjoined
along its branch. These bounds may have an error $o(1)$ as $n\to\infty$.
The initial table for $T^{\otimes n}$ is $N(x)=H(T;x)$.

Choose one mark of original dimension $q$, with endpoints $a,b$ at
distance $d>0$. Its current dimension is $Q=q^n$. On exterior
degree $l$ of $k^Q$, multiplication by a basis vector $i$ and
contraction by the dual vector $j$ satisfy
\[
                 D_jL_i+L_iD_j=\mathbf 1_{i=j} I .
\]
Adjoin the $ab$-wire of dimension $R=\binom Ql$. Use this identity on
it, summing against the old index $j$ and leaving the new index $i$
marked. The two terms form two child nodes. At either endpoint a child
can be made using both operators there and the $R$-wire. In the first
term, with the new index at $a$, one may instead put creation at $a$
and the old index with contraction at $b$, using a wire of dimension
$R_+=\binom Q{l+1}$. In the second term the reversed placement uses
$R_-=\binom Q{l-1}$ between the operators.

Write $N_0,N_1$ for the two slices at $a,b$. For any
$|c|\le d\log_2 q$ choose $l$ so that
$(d/n)\log_2(R_\pm/R)=\pm c+o(1)$. For $s=c/(d\log_2 q)>0$ take
$l=\lfloor Q^{1-s}\rfloor$, for $s<0$ take the complementary level
for $-s$, and for $s=0$ take a middle level. After subtracting the
price of the adjoined wire the child tables may therefore be taken as
\[
 \bigl(\min\{N_0,N_1+c\},N_1\bigr),\qquad
 \bigl(N_0,\min\{N_1,N_0-c\}\bigr),
\]
with entrywise minima. The passage from $N_1$ to $N_1+c$, for instance, is precisely the alternative $R_+$ implementation at $a$.
A zero child may be discarded.
Choose $c$ midway between the extreme entries of $N_0-N_1$. Writing $D=N_0-N_1$, the new differences are $\min(D,c)$ and
$\max(D,c)$. Thus the range of differences starts in
$[-d\log_2 q,d\log_2 q]$ by \eqref{eq:move}, and each split reduces
its width by at least half. The bounds on differences along every other
coordinate are preserved, since entrywise minimum is monotone and
commutes with addition of a constant.

After sufficiently many rounds with every mark, the variation of each
coordinate difference in any leaf table is at most $\epsilon$.
Consequently each such table with $r$ marks lies within $r\epsilon$
of a function of the form $c_0+\sum_i f_i(x_i)$, by telescoping
coordinate changes from one corner. The means of this latter function
under the product law and $\mu$ agree. Apply independent relocation to
the leaf tensor, with auxiliary indices fixed, and subtract the accumulated
wire prices. Its price is bounded by the $\mu$-mean of its table plus
$2r\epsilon+o(1)$, hence by the same expression using the original table,
as entries only decreased. At each split the new wire's price cancels
and ordinary subadditivity costs at most $1/n$ in the normalized log.
For fixed $\epsilon$ the depth is fixed. Let $n\to\infty$ and then
$\epsilon\downarrow0$.
\end{proof}


Write $[a,b]_t$ for a point at distances $t\,d(a,b)$ from $a$ and
$(1-t)d(a,b)$ from $b$.
\begin{corollary}[Geodesic extension]\label{cor:spectral-geodesic}
Every spectral space embeds with prices unchanged into a complete uniquely geodesic spectral space.
Prices are convex along simultaneous geodesics. In particular in any spectral space,
if the two interval points exist,
\begin{equation}\label{eq:wire-busemann}
 d([a,b]_t,[c,e]_t)\le (1-t)d(a,c)+t d(b,e)
                       \qquad (0\le t\le1).
\end{equation}
\end{corollary}
\begin{proof}
Apply Proposition~\ref{prop:joint-convexity} to the ends of a binary
wire with the same endpoint choice for both ends. At a common fraction
$t$ the distance between interval points is bounded by
$(1-t)d(a,a')+t d(b,b')$ for their endpoint pairs $(a,b),(a',b')$.
This implies uniqueness of each interval point. Use
Theorem~\ref{thm:extension} and its infinite-space version to insert midpoints for all pairs, taking
successive extensions and their union for an infinite set of pairs.
Repeat in countably many rounds, then take the completion. The midpoint
estimate extends midpoints to the completion, and iterated bisection
gives constant-speed segments. Propagating prices through these steps as
above gives a complete geodesic spectral space $X$. On such a space
the synchronized endpoint choice in Proposition~\ref{prop:joint-convexity}
gives convexity of every price along simultaneous geodesics.
\end{proof}
~~~

## G / curvature.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/curvature.tex, Git blob 9269654e26f361ff7b130ae3adcc8007db719a04. Through the CAT(0) and four-point inequalities. Later flats and the illustrative figure are omitted.

~~~latex
\section{From crossing decompositions to nonpositive curvature}
\label{sec:crossing}

The geodesic extension of Section~\ref{sec:geodesic} gives convexity of
wire distances. To obtain nonpositive curvature, we need a quadratic
bound on those distances. Its algebraic source is a decomposition of two
wires into two tensors that factor across different cuts of the four
endpoints. The cost of the internal wire on one cut can be traded
against its cost on the other. We use this tradeoff first to prove the
parallelogram law for normed planes realized by wire distances and then
to prove the squared midpoint inequality in every geodesic spectral
space.

\subsection{The integral crossing identity}\label{app:crossing}

The crossing construction originates in work of Anthropic~\cite{Anthropic}.
Khoruzhii~\cite{Khoruzhii} gives an integral construction in which two
products of fermionic operators factor through states below the same
energy bound. Extracting a shifted coefficient changes the two bounds
by equal and opposite amounts.
\begin{theorem}[Shifted crossing]\label{thm:crossing}
Let $t,u>0$, put $S=t^2+u^2$, and take $-S<\delta<S$. There are
integers $d_X,d_Y$, depending on $n$, with
$\log_2 d_X=nt+o(n)$ and $\log_2 d_Y=nu+o(n)$, for which the
two-wire tensor has the ordinary sum decomposition
\[
 E_{X_+X_-}(d_X)\,E_{Y_+Y_-}(d_Y)=T_++T_-.
\]
On the cuts $X_+Y_+\mid X_-Y_-$ and $X_+Y_-\mid X_-Y_+$ respectively,
$T_+$ and $T_-$ have flattening ranks at most
\[
  2^{n\sqrt{S+\delta}+o(n)}
  \quad\text{and}\quad
  2^{n\sqrt{S-\delta}+o(n)} .
\]
One can choose the tensors and their two matrix factorizations over
$\mathbb Z$, so the same choices work over every field.
\end{theorem}

The two terminal pairs $(X_+,X_-)$ and $(Y_+,Y_-)$ carry the original
wires. Each summand can be represented by joining two local factors
with a single internal wire, as in Figure~\ref{fig:crossing-cuts}.
In the integral construction, the internal wire records the intermediate
state of a product of operators together with auxiliary summation indices.
There are polynomially many choices of the auxiliary indices at the given
energy cutoff. Reversing the order of the operators exchanges the two cuts. The shift $\delta$ lets us choose how much dimension to
allow on each internal wire.


% Editorial omission of an illustrative figure.


\input{crossing-appendix}

\subsection{Rigidity of normed planes}

Place the four terminals in a normed plane. The two auxiliary sites in
each factorization may be placed anywhere. If the norm violates the
parallelogram law, we can place these sites so that both summands cost
strictly less at exponential scale than the original pair of wires.
Subadditivity of a spectral point rules this out.
\begin{lemma}[Norm rigidity]\label{lem:norm}
Let $V$ be a real normed plane. If each finite subset of some
nonempty open set in $V$ is realized at some scale by a spectral wire
metric over the fixed field, then the norm comes from an inner product.
\end{lemma}

\begin{proof}
Any finite subset of the plane can be put inside the open subset by a
translation and a positive dilation. Let $x,y$ be unit vectors at
which the norm is differentiable, and let $f,g$ be its derivatives
there, elements of the dual unit sphere. We prove that for $t,u>0$,
\begin{equation}\label{eq:crossing-dual}
   \|tf+ug\|_*^2+\|tf-ug\|_*^2\le 2(t^2+u^2).
\end{equation}
Suppose the inequality fails. Put $S=t^2+u^2$ and
$A_\pm=\|tf\pm ug\|_*$. Both are positive and
$A_\pm^2\le(t+u)^2\le 2S$. Thus one can choose $-S<\delta<S$
such that $S-A_-^2<\delta<A_+^2-S$. Take unit vectors $z_\pm$
on which $tf\pm ug$ attain their respective positive maxima.
Place the terminals at $x,-x,y,-y$. On the first cut one center serves
$x,y$ and the opposite center serves $-x,-y$. On the second cut the
partners change to $x,-y$ and $-x,y$. Starting both centers at the
origin gives each placement the same terminal cost. Moving them apart
shortens the terminal wires at a rate controlled by $tf+ug$ or $tf-ug$,
respectively, and requires a wire between the two centers. Choose the
centers $\epsilon z_+,-\epsilon z_+$ on the first cut and
$\epsilon z_-,-\epsilon z_-$ on the second, with $0<\epsilon<1$.

For any matrix factorization of a four-terminal tensor across a cut,
give each terminal an identity wire of its dimension to the center on
its side, and join the centers by a wire of the factorization size.
Contracting the indices at the centers against the two factors
produces the tensor by local restriction, with trivial output at
both centers. Coincident auxiliary roles may be grouped at one party.
Apply this observation to $T_+$ and $T_-$ of
Theorem~\ref{thm:crossing}. Realize the terminals and both pairs of
centers simultaneously, with wire distances $\sigma$ times their
norm distances for some $\sigma>0$. By multiplicativity, restriction
monotonicity, and the formula $F(E_{ab}(q))=q^{d(a,b)}$ of
Section~\ref{sec:spectra}, the log base two of $F(T_+)$, divided by
$2n\sigma$, is at most
\[
 t\|x-\epsilon z_+\|+u\|y-\epsilon z_+\|
                 +\epsilon\sqrt{S+\delta}+o(1).
\]
Here and below a tensor of value zero satisfies the corresponding
exponential bound directly. For $T_-$ the same expression applies
with $y$ replaced by $-y$, $z_+$ by $z_-$ and $\delta$ by $-\delta$.
As $\epsilon\downarrow0$ the leading expressions are
\[
 t+u-\epsilon\big(A_\pm-\sqrt{S\pm\delta}\big)+o(\epsilon),
\]
because the norm derivative at $-y$ is $-g$. Both linear savings
are strictly positive. Fix $\epsilon$ small enough for both, and
then fix the joint realization, independently of $n$. Subadditivity
under ordinary tensor sums now contradicts the wire formula for the
target, whose log value divided by $2n\sigma$ tends to $t+u$.
This proves \eqref{eq:crossing-dual}.

The closed convex hull of all derivatives at smooth unit vectors
is the dual unit ball. Indeed the norm is differentiable on a dense
set by finite-dimensional convexity. For any unit vector $v$, take
smooth unit vectors $v_j\to v$ with derivatives $f_j$. Then
$f_j(v)\ge f_j(v_j)-\|v-v_j\|\to1$, so these derivatives have the
same support function as the dual unit ball. Joint convexity and
continuity in $(f,g)$ now extend \eqref{eq:crossing-dual} to arbitrary
$f,g$ in that ball. Taking $t,u$ to be arbitrary dual vector norms
gives
\[
   \|a+b\|_*^2+\|a-b\|_*^2\le 2(\|a\|_*^2+\|b\|_*^2).
\]
Applying this inequality to $a+b,a-b$ proves its reverse. The
parallelogram law gives an inner product on the dual, and duality in
finite dimension gives the assertion for the original norm.
\end{proof}

\subsection{Nonpositive curvature}\label{sec:curvature}
On a spectral space, distances between corresponding points of two
geodesics are convex by \eqref{eq:wire-busemann}. Norm rigidity upgrades
this convexity to the squared midpoint bound characteristic of
$\operatorname{CAT}(0)$ spaces. The resulting extension preserves the prices of all
tensors on the original sites.
\begin{theorem}[Curvature]\label{thm:geometry}
Every complete geodesic spectral space $X$ is $\operatorname{CAT}(0)$.
In particular every spectral space embeds with prices in such a space, called below a \emph{Hadamard spectral space}.
\end{theorem}
We first examine two geodesics from a common origin whose separation
grows exactly linearly. They span a Euclidean sector, possibly
degenerate. We then rescale an arbitrary pair of geodesics at their
origin. Convexity makes the separation linear in this limit, so the
sector calculation gives an inequality for the original distances.

Fix a complete geodesic spectral space $X$.
Consider two constant-speed segments $a,b$ starting at $o$, defined on $[0,R]$, or two rays on $[0,\infty)$. Suppose $d(a(r),b(r))=cr$ throughout. On the coordinate triangle $u,v\ge0$, $u+v\le R$, or on the whole quadrant for rays, define
\[
 H(u,v)=[a(u+v),b(u+v)]_{v/(u+v)},\qquad H(0,0)=o.
\]
This map takes affine segments to geodesics with their affine parametrizations. Indeed, let $(u,v)$ be the affine interpolation of two parameter pairs, and $w$ the corresponding interpolation of their images. Comparing with the same interpolation on each boundary ray, \eqref{eq:wire-busemann} gives $d(w,a(u+v))\le cv$ and $d(w,b(u+v))\le cu$. Their sum is bounded by the distance $c(u+v)$ between these boundary points. Equality and uniqueness of geodesics force $w=H(u,v)$. The formula and \eqref{eq:wire-busemann} also give continuity, including at the origin.

The pullback distances of $H$ are given by a seminorm on the coordinate plane. To see this, let $N_p(v)$ be the speed of $H(p+hv)$ as a function of $h$ near $0$, for an interior point $p$. For any two interior points $p,q$ choose $s<1$ and an interior point $z$ with $q=(1-s)p+sz$. Comparing $q+hv$ to $q$ by interpolating from $p+hv/(1-s)$ and $p$ to the same $z$ shows $N_q(v)\le N_p(v)$. Interchanging $p,q$ gives equality. The common speed is absolutely homogeneous, and the triangle inequality on small segments gives subadditivity. The segment property and continuity extend the distance formula to the boundary. In rank two the seminorm is a norm, and $H$ embeds the interior of its parameter domain isometrically into $X$ with this norm. Lemma~\ref{lem:norm} therefore makes the norm Euclidean. In rank at most one it is already induced by a positive semidefinite quadratic form. Thus the sector, with all its interpolations, has the distances of a possibly degenerate Euclidean sector.

\begin{lemma}[Squared midpoint inequality]\label{lem:catzero}
For $x,y,z\in X$ and $m=[x,y]_{1/2}$,
\begin{equation}
 d(z,x)^2+d(z,y)^2\ge 2d(z,m)^2+\tfrac12 d(x,y)^2.
 \label{eq:cur-mid}
\end{equation}
At equality, the two sectors spanned at $m$ by the segments to $x,z$ and to $y,z$ have Euclidean distances, allowing degeneracy.
\end{lemma}

\begin{proof}
Take the three segments from $m$ to $x,y,z$, each on $[0,1]$. Take a pointed ultralimit at $m$ after rescaling distances by integers $j\to\infty$, still geodesic (use representative segments), carrying spectral structure by Section~\ref{sec:limits}. Take the limit of the three segments with domain rescaled by the same $j$, obtaining three rays (possibly constant). For two original segments $a,b$, the ratio $d(a(s),b(s))/s$ is nondecreasing by \eqref{eq:wire-busemann} and bounded by the sum of their speeds. Its limit at zero therefore gives a constant $c$ for which the equal-parameter distances of the limiting rays are $cr$.

Put $A=d(x,m)=d(y,m)$ and $B=d(z,m)$, initially both positive. The rays for $x$ and $y$ point in opposite directions on one geodesic. The preceding sector argument gives Euclidean angle cosines $u,v$ for their respective pairs with the third ray. The triangle inequality between the first two rays at parameter $T$, through the third at parameter $1$, reads
\[
 2TA\le\sqrt{(TA)^2+B^2-2TABu}+\sqrt{(TA)^2+B^2-2TABv}.
\]
Subtract $2TA$ and let $T\to\infty$ to find $u+v\le0$. Each original distance $d(z,x),d(z,y)$ is at least the corresponding limiting distance at parameter $1$, by monotonicity of the ratios. The sum of their squares is therefore at least $2(A^2+B^2)-2AB(u+v)\ge2(A^2+B^2)$, proving \eqref{eq:cur-mid}. Equality forces each original ratio at parameter $1$ to equal its limit at zero. The ratio is then constant on $(0,1]$, and the sector argument applies to each pair of original segments. If $A$ or $B$ vanishes, the assertions follow directly on a segment.
\end{proof}

Applying \eqref{eq:cur-mid} at dyadic fractions and passing to the limit gives the CAT(0) inequality in its squared-distance form
\begin{equation}
 d(z,[x,y]_s)^2\le (1-s)d(z,x)^2+s d(z,y)^2-s(1-s)d(x,y)^2 .
 \label{eq:cur-cn}
\end{equation}
The midpoint characterization \cite[Ch.~II.1]{BH} now gives the curvature
assertion of Theorem~\ref{thm:geometry}. Corollary~\ref{cor:spectral-geodesic}
gives the embedding assertion. To understand when a finite configuration
spans a Euclidean rectangle, we use the four-point consequence of
\eqref{eq:cur-mid}. For $A,B,C,D\in X$ it reads
\begin{equation}
 d(A,C)^2+d(B,D)^2
 \le d(A,B)^2+d(B,C)^2+d(C,D)^2+d(D,A)^2 .
 \label{eq:cur-quad}
\end{equation}
To derive this bound, apply \eqref{eq:cur-mid} on $AC$ with $z=B,D$, and use $d(B,D)^2\le 2d(B,m)^2+2d(m,D)^2$. At equality, $m$ is also the midpoint of $BD$, and \eqref{eq:cur-mid} is sharp on either diagonal with either remaining vertex. Consequently the two side squares incident at any vertex sum to $(d(A,C)^2+d(B,D)^2)/2$, and opposite sides have equal lengths.
~~~

## G / crossing-appendix.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/crossing-appendix.tex, Git blob 188af506e2a3abde8732967a0a6f69ea22e77970. Complete shifted-crossing proof.

~~~latex
\begin{proof}[Proof of Theorem~\ref{thm:crossing}]
The construction has three steps. We identify an integral pairing
of operators, factor it through states with two different energy
bounds, and count the states allowed by each bound.

\smallskip\noindent\emph{The pairing.}
Use the exterior algebra over $\mathbb Z$ on generators indexed by
$i\ge 0$. Its basis states are finite subsets of these indices,
written in increasing order. Let $c_i$ wedge the $i$th generator on
the left and let $a_i$ be its signed deletion. Creation operators
anticommute with one another, as do deletion operators, and
$a_i c_j+c_j a_i=\mathbf 1_{i=j}\mathrm{Id}$. Give a basis state $I$
energy $h(I)=\sum_{i\in I}(i+1/2)$. For background on the free fermion and vertex superalgebras see \cite{KacWang}. Introduce the fields
\[
 \psi(z)=\sum_{i\ge0}(c_i z^i+a_i z^{-i-1}),
 \qquad \psi_j(z)=\partial^{(j)}_z\psi(z),
 \qquad Y_I(z)={:}\prod_{i\in I}\psi_i(z){:}.
\]
Here $\partial^{(j)}_z z^b=\binom bj z^{b-j}$ for every integer $b$,
using the coefficient of $v^j$ in $(z+v)^b$. The product is ordered
by increasing $i$. The colons instruct us to expand and move creations
to the left of deletions using permutation signs only. Such a normal
product is skew under interchange of two fields. A coefficient $[z^b]Y_I$
shifts energy by $b+h(I)$. Indeed each summand of a differentiated
field shifts energy by its exponent plus the index of differentiation
plus $1/2$.

For equally large nonempty sets $I,J$, put $m=|I|=|J|$, $h=h(I)$,
$h'=h(J)$ and $L=h+h'-1$, a nonnegative integer. We claim
\begin{multline}\label{eq:crossing-field}
 (z-w)^L\big(Y_I(z)Y_J(w)-(-1)^{m^2}Y_J(w)Y_I(z)\big)
   \\= C_{IJ}\sum_{s\in\mathbb Z}z^{-s-1}w^s\,\mathrm{Id},
\end{multline}
where
\[
 C_{IJ}=(-1)^{\binom m2+\sum_{i\in I}i}
       \det\left[\binom{i+j}{i}\right]_{i\in I,\,j\in J}.
\]
To see this, normal order the product of the two blocks. When a
deletion from the first block meets a creation from the second, the
anticommutation relation permits either a signed interchange or a
scalar contraction. Summing the scalar for fields $i,j$ gives the
expansion in $w/z$ of
\[
  \partial^{(i)}_z\partial^{(j)}_w\frac{1}{z-w}
   =\frac{(-1)^i\binom{i+j}{i}}{(z-w)^{i+j+1}},
\]
since the undifferentiated contraction is
$\sum_{r\ge0}z^{-r-1}w^r$. Thus the expansion is indexed by partial
matchings between the blocks. Each summand is the product of the
matched scalars times the normal product of the unmatched fields.
Its sign is that of the permutation listing the matched pairs first,
in first-field/second-field order within each pair, followed by the
unmatched fields in their original order. In fact the sign for a
contraction is the sign of bringing the two operators together.
Moving the resulting pair as a block incurs no further sign.

For the reversed product, use the same ordering of pairs and the
same final order of the unmatched fields. Interchanging the original
blocks contributes $(-1)^{m^2}$. Reversing the orientation of any
pair contributes a further minus sign in the permutation, precisely
canceled by reversing its rational kernel, since
$(-1)^j/(w-z)^{i+j+1}=(-1)^{i+1}/(z-w)^{i+j+1}$. Consequently after
multiplication by $(-1)^{m^2}$ the reversed summands expand the same
rational kernels in $z/w$. For an incomplete matching an unmatched
field in each block accounts for at least one unit of missing pole
order. Its pole order is therefore at most $L$, and multiplication
by $(z-w)^L$ makes the two expansions equal. For a full matching the
pole order is $L+1$. Interleaving the two increasing lists contributes
$(-1)^{\binom m2}$, permuting partners contributes the determinant
sign, and the kernels contribute the remaining row signs in $C_{IJ}$.
The two expansions of $1/(z-w)$ differ by
$\sum_{s\in\mathbb Z}z^{-s-1}w^s$. This proves
\eqref{eq:crossing-field}.

We justify the formal products coefficientwise. Truncate the modes at
$N$ and first
perform the two normal-ordering expansions with their truncated scalar
kernels. Fix input and output states and a monomial $z^b w^{b'}$.
In the first product the intermediate energy is the input energy plus
$b'+h'$. In the reversed product the increment is $b+h$. There are
only finitely many states at either energy. Each individual normal
block between specified states uses deletion modes from its input
and creation modes from its output. Hence its matrix coefficient
stabilizes as a Laurent polynomial. The same argument applies to
the unmatched normal product, even in two variables. For each of its
finitely many monomials, prescribing the total $w$-power bounds all
contracted modes in the first expansion: a creation in mode $r$
differentiated $j$ times has $w$-power $r-j\ge0$ whenever its
coefficient is nonzero. The reversed expansion uses the $z$-power in
the same way. Both Wick expansions therefore pass coefficientwise to
the infinite-mode limit, where the rational-kernel cancellation above
applies. Multiplication by the fixed polynomial $(z-w)^L$ preserves
these limits. All coefficients are integers, including the divided
derivatives, so reduction modulo $2$ also preserves the identity.

\smallskip\noindent\emph{The factorizations.}
Let $D(H)$ count the states of energy at most $H$, with $D(H)=0$ for
$H<0$. Fix $H,E\ge1$ and an integer $s$. Use nonempty $I,J$ of
energy at most $H$, and input and output states $\beta,\alpha$ of
energy at most $E$. For a block of given $m,h,h'$ take
$[z^{-s-1}w^s]$ in \eqref{eq:crossing-field}. For the multiplier
monomial $z^{L-a}w^a$, $0\le a\le L$, the first product acts on
$\beta$ first by $[w^{s-a}]Y_J$, reaching energy at most
$E+h'+s-a\le E+H+s$. The reversed product acts first by
$[z^{-s-1-L+a}]Y_I$, reaching energy at most
$E-s-h'+a\le E+H-s$. Sum over intermediate states in each matrix
product. This writes $C_{IJ}\mathbf 1_{\alpha=\beta}$ as a sum of
two matrices, each factored over $\mathbb Z$. On the cut
$(I,\alpha)\mid(J,\beta)$ the first factors through at most
$(L+1)D(E+H+s)$ coordinates, indexed by $a$ and the intermediate
state. On the cut $(I,\beta)\mid(J,\alpha)$ the second factors
through at most $(L+1)D(E+H-s)$ coordinates. The integer coefficients
of the multiplier and the sign of the reversed product can be placed
on either factor.

Extend each block by zero and sum over $m,h,h'$. The number of blocks
and the largest $L+1$ are bounded by a polynomial in $H+1$, since
all energies are half-integers and $m\le 2H$. Define $C_{IJ}=0$ for
$|I|\ne |J|$. The resulting target is the full matrix
$C\otimes\mathrm{Id}_{D(E)}$ on the four indices. The two factorizations
have respective sizes at most
\[
  \operatorname{poly}(H+1)D(E+H+s),\qquad
  \operatorname{poly}(H+1)D(E+H-s).
\]
To turn $C$ into an identity, let $P_{i\ell}=\binom{i}{\ell}$
for nonnegative indices up to the largest mode under consideration.
Vandermonde's identity gives
$(PP^{\mathsf T})_{ij}=\sum_\ell\binom{i}{\ell}\binom{j}{\ell}
=\binom{i+j}{i}$. For fixed $m$, let $P^{[m]}$ be the matrix of
$m$-rowed minors restricted to sets of energy at most $H$ on both
indices. A minor of $P$ vanishes unless its increasing column list
is componentwise bounded by its increasing row list, since $P$ is
lower triangular. Thus restricting to this energy downset discards
no term of Cauchy--Binet, and the corresponding block of $C$ is
$P^{[m]}(P^{[m]})^{\mathsf T}$ up to row signs. In lexicographic
order $P^{[m]}$ is unit lower triangular. Every block is therefore
unimodular, and applying $C^{-1}$ at the $I$ index gives the
identity of dimension $D(H)-1$ without increasing either factorization
size.

\smallskip\noindent\emph{Counting the intermediate states.}
Write $\ln$ for the
natural logarithm. With
$A=\int_0^\infty\ln(1+e^{-x})\,dx>0$, the partition function is
\[
  \mathcal Z(\lambda)=\prod_{i\ge0}(1+e^{-\lambda(i+1/2)}),
  \qquad \ln\mathcal Z(\lambda)\sim A/\lambda
  \quad(\lambda\downarrow0).
\]
Under probabilities proportional to $e^{-\lambda h(I)}$ the occupation
variables are independent, with means $1/(e^{\lambda(i+1/2)}+1)$.
The mean energy is asymptotic to
$\lambda^{-2}\int_0^\infty x/(e^x+1)\,dx=A/\lambda^2$,
by integration by parts. Its variance is $O(\lambda^{-3})$ by
summing $(i+1/2)^2e^{-\lambda(i+1/2)}$. These Riemann-sum
estimates have exponentially decreasing tails. From
$D(H)\le e^{\lambda H}\mathcal Z(\lambda)$ at
$\lambda=\sqrt{A/H}$ we obtain the upper asymptotic bound. For the
lower bound take $0<\epsilon<1/2$ and
$\lambda=\sqrt{A/((1-\epsilon)H)}$. Chebyshev's inequality puts
probability $1-o(1)$ on energies between $(1-2\epsilon)H$ and $H$.
The number of states in this interval is at least
$(1-o(1))\mathcal Z(\lambda)e^{\lambda(1-2\epsilon)H}$.
Letting $\epsilon\downarrow0$ gives
\[
     \log_2 D(H)\sim c\sqrt H,\qquad c=2\sqrt A/\ln 2 .
\]
Finally take $H=(nt/c)^2$, $E=(nu/c)^2$ and
$s=\lfloor\delta n^2/c^2\rfloor$. Both shifted energy bounds are
positive for large $n$. Set $d_X=D(H)-1$, $d_Y=D(E)$ and identify
$(I,J,\alpha,\beta)$ with $(X_+,X_-,Y_+,Y_-)$. The polynomial
factors cost $o(n)$ in the logarithms, proving the theorem.
\end{proof}
~~~
