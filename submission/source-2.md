# Source evidence 2. Ports and exact compression

Excerpts from OpenAI's mathematics preprints at commit 92b516ef8e2abf66c3db7aa31438d3010bc7956f. The dependency map, bibliography and Apache 2.0 license are in source-5.md. The authors' result labels are preserved.

## G / port-graphs.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/port-graphs.tex, Git blob 418a2e0442817bb69cbbd22abbf8b3009624c289. Port definition, equivalence and normalization. The subsequent cube-exclusion implication is unused.

~~~latex
Let $G$ be a symmetric stochastic matrix on a finite vertex set $V$
and assign a point $x_i$ of a spectral space to each vertex $i$. The
mean squared length of a step of $G$ and the mean squared distance
between two independent uniform vertices are, respectively,
\[
 \mathcal E_G(x)=\frac1{|V|}\sum_{i,j\in V}G_{ij}d(x_i,x_j)^2,
 \qquad
 \mathcal V(x)=\frac1{|V|^2}\sum_{i,j\in V}d(x_i,x_j)^2.
\]

\begin{definition}\label{def:ports}
A field \emph{has ports} if it admits a family of regular, symmetric
multigraphs of unbounded order and bounded degree, and a constant
$\gamma$ such that
\[
 \mathcal V(x)\le\gamma\mathcal E_G(x)
\]
for every graph in the family and every assignment into a spectral space
over the field. Here $G$ denotes normalized adjacency. A diagonal arc
keeps the walk at its current vertex.
\end{definition}
This property depends only on the characteristic, since it concerns
only joint prices of binary wires (Lemma~\ref{lem:comparison}).
By Theorem~\ref{thm:geometry} one may test either finite spectral
points or complete $\operatorname{CAT}(0)$ spectral spaces. Testing
an assignment that is constant on each of two components and takes two
distinct values shows that every port graph is connected.

For use in Section~\ref{subsec:refresh}, we may take the graphs to be
loopless, of arbitrarily large power-of-two order and of one fixed
power-of-two degree. Given a graph of order $M$, duplicate enough
vertices to reach the next power of two. Keep the old edges between
the original representatives, join each duplicate to its representative,
and pad with holding arcs to a fixed power-of-two degree. To compare any
pair, follow the path through its two representatives. A uniform endpoint
chooses any particular representative with probability at most $2/M$.
The mean squared distance from a vertex to its representative is controlled
by a constant times the new edge energy, since each such pair is an edge
and the degree is bounded. The triangle inequality therefore gives a
uniform Poincar\'e bound after duplication. To remove loops, take two
copies of this graph. Replace each holding arc by a crossing to the
corresponding vertex in the other copy, and join every such pair of vertices
by as many additional parallel edges as the previous degree. The degree
doubles. The same comparison through representatives again gives a
uniform bound.
~~~

## G / compressor.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/compressor.tex, Git blob a5368884147dca2dac0442f715ac8727ed60bd29. Complete mathematical text. The illustration is omitted.

~~~latex
\part{Descent and contradiction}
\input{port-graphs}

\section{Fixed wire identities from metric inequalities}
\label{sec:motifs}

A product of wires can be replaced by a sum of products on a different
graph. For compression, we need a bound on the total wire width in
every summand after its endpoints are regrouped. Theorem~\ref{thm:conversion} turns a comparison
of squared wire distances into a bound on the widths in each summand.
Lemma~\ref{lem:motif} supplies this comparison for a six-party graph
that admits two ways of regrouping its endpoints. Both results hold
over an arbitrary field, denoted here and in Section~\ref{sec:compressor}
by $\mathbb F$.

Fix a finite set of at least two parties, a multiset $D$ of demand edges,
and a connected spanning multigraph $S$ of supply edges. Each edge has distinct endpoints.
For nonnegative integer widths $b=(b_e)_{e\in D}$,
let $W_D(b)$ be the product of identity tensors of dimensions $2^{b_e}$
on their respective edges. Define $W_S(c)$ similarly. Give the edges
positive rational masses $p_e$ on $D$ and $q_e$ on $S$, and set
\[
  \|b\|_p^2=\sum_{e\in D}\frac{b_e^2}{p_e},\qquad
  \|c\|_q^2=\sum_{e\in S}\frac{c_e^2}{q_e}.
\]
For a spectral point $F$ on these parties over $\mathbb F$, write $d_e$ for its
wire distance on the endpoints of $e$. In particular,
$\log_2 F(W_D(b))=\sum_{e\in D}b_e d_e$.
We use this comparison on six parties to shrink the total width and later on a block of ports to restore bounded degree.
The squared-distance energies in Equation~\eqref{eq:motif-rayleigh} are dual to the width energies above. For
example, the maximum of $\sum_{e\in D}b_e d_e$ over real $b\ge0$ with
$\|b\|_p\le1$ is $(\sum_{e\in D}p_e d_e^2)^{1/2}$.
The next theorem turns this comparison into exact identities for
integer widths.

\begin{theorem}[Conversion of a fixed motif]
\label{thm:conversion}
Let $H\ge0$ and suppose that for every such spectral point,
\begin{equation}
 \sum_{e\in D}p_e d_e^2\ \le\ H\sum_{e\in S}q_e d_e^2 .
 \label{eq:motif-rayleigh}
\end{equation}
For every $H'>H$ and all $b$, there is an ordinary sum identity
\begin{equation}
 W_D(b)=\sum_\alpha
       \left(\bigotimes_v L_{v,\alpha}\right) W_S(c_\alpha),
 \qquad
 \|c_\alpha\|_q^2\le H'\|b\|_p^2+O(1).
 \label{eq:motif-conversion}
\end{equation}
The constant in the energy estimate is independent of $b$. The identities
use finitely many scalar parameters $\kappa$ subject to fixed integer
polynomial equations with a solution in $\mathbb F$. They hold over every
commutative unital ring with any assignment satisfying these equations.

With $M=\max(\{0\}\cup\{b_e:e\in D\})$, each index $\alpha$ consists of
a fixed number of integers of magnitude $O(M+1)$. Admission and supply
widths are computable in $\operatorname{poly}(\log(2+M))$ bit operations.
Given the index and its incident bit strings, an entry of $L_{v,\alpha}$
is specified as zero or as a monomial in $\kappa$ of degree $O(M+1)$.
Its description is computable in time polynomial in $M+2$ with an
absolute exponent.
\end{theorem}

\begin{proof}
Choose a rational $R$ strictly between $\sqrt H$ and $\sqrt{H'}$.
The case $b=0$ uses one unit term. Otherwise round $b/\|b\|_p$
coordinatewise upwards to a sufficiently fine fixed rational grid. This
gives one of finitely many vectors $z$ with $\|z\|_p\le 1+\epsilon$,
where $(1+\epsilon)\sqrt H<R$. Weighted Cauchy--Schwarz and
Equation~\eqref{eq:motif-rayleigh} give
\[
  \sum_{e\in D}z_e d_e
   <R\left(\sum_{e\in S}q_e d_e^2\right)^{1/2}
    =\max_{w\ge0,\,\|w\|_q\le R}\sum_{e\in S}w_e d_e .
\]
Indeed the supply norm on the right is positive. If it vanished, a
product of binary wires along a spanning tree would have value one,
whereas it restricts to the two-copy diagonal, of value two. Rational
approximation inside the supply ball and compactness of the spectrum
give a finite list $\mathcal W$ of nonnegative rational $w$ that suffices
for the strict inequality for every $z$ and $F$.

Choose a large common denominator multiple $B$ for both grids. For each
$z$, every spectral point strictly orders $W_D(Bz)$ below
$Y=\bigoplus_{w\in\mathcal W}W_S(Bw)$. By strict spectral comparison
(Lemma~\ref{lem:spectral-comparison}), some fixed $K_z$ copies of $Y^{\otimes h}$
restrict to $W_D(hBz)$ for every $h\ge1$. Expand these direct sums into
an ordinary sum of images of supply tensors under local maps. Each
base term has supply-width vector of norm at most $hBR$ by the triangle
inequality. Its index fits in
\[
   s=\left\lceil\log_2(\max_z K_z)+h\log_2|\mathcal W|\right\rceil
\]
bits. We may make $s/(hB)$ as small as desired by choosing $B$ large
first and then $h$ large. Fix $B,h$ at these values. Introduce one
parameter for each entry of each base local map, with the base
identities themselves imposed entrywise as polynomial equations.
These finitely many equations have integer coefficients and a simultaneous
solution in $\mathbb F$.

For the given nonzero $b$, take $u=\lceil\|b\|_p/(hB)\rceil$ copies of
the base identity at $z$. Its demand widths $uhBz$ dominate $b$. Fix
all surplus demand bits to zero at both ends. To avoid listing a separate
outer term for every sequence of $u$ base selectors, group sequences
by their total supply payload $v=(v_e)_{e\in S}$. Allow every integral
$v\ge0$ with $\|v\|_q\le uhBR$, including values that come from no
sequence. For this one outer index, use width $c_e=v_e+us$ on supply
edge $e$.

At each party, read the $us$ selector bits on every incident supply edge.
Require the lists to agree, require every selector to be valid, and
require the sequence's total payload vector to be the stated $v$.
Use the selectors in order to split each incident payload into the
corresponding base fields. Return the product of the selected base-map
entries, with the demand padding prescribed above. Failure of any check
gives a zero entry. Since the supply graph is connected, a surviving
assignment of selector copies is determined by a single global sequence.
Conversely every sequence gives exactly one such assignment. Summing
the payload bits therefore reproduces the product identity, with each
sequence counted once in its $v$-term. This proves the identity over every
ring satisfying the base equations.

Writing $\mathbf 1$ for the all-ones vector on $S$, the resulting widths
satisfy
\[
  \|c\|_q\le uhB\left(R+\frac{s}{hB}\|\mathbf 1\|_q\right).
\]
Take the selector fraction small enough that the factor in parentheses
is strictly less than $\sqrt{H'}$. Since $uhB\le\|b\|_p+hB$, squaring
gives the bound in Equation~\eqref{eq:motif-conversion} because the strict gap in the
quadratic coefficient absorbs the linear term up to a constant.

The outer index $v$ has a fixed number of coordinates of size $O(M+1)$.
Deciding the grid cell of $b/\|b\|_p$, computing $u$, and testing its
norm bound use integer comparisons after clearing fixed denominators and
squaring nonnegative quantities. All operands have $O(\log(2+M))$ bits.
To evaluate one map entry, scan the $u=O(M+1)$ fixed-size selectors and
their fields, incrementing one parameter exponent per copy. All local
dimensions of the base identities are constants. Even rescanning all
incident strings to locate each successive field takes only
$O((M+2)^2\operatorname{polylog}(M+2))$ bit operations with a fixed
multiplicative constant. This gives the stated uniform exponent.
\end{proof}

\subsection{Two layouts on six parties}
The endpoints $A,B$ represent sites where several vertices have been
packed together. We pair their demand wire $AB$ with a second demand
wire $CD$. The placement of $D$ may already be fixed by another
occurrence of the motif. We retain the freedom to place $C$ at either
$A$ or $B$, and use two auxiliary parties $M,N$ to make a different
group of supply wires local in each case. For a related packing
construction, see \cite[Appendix~D]{falseeth}.

Consider the six-party graph with demands $AB,CD$ and supply
groups
\[
 G_1=\{AC,AM,CM,DN\},\qquad G_2=\{BC,BN,CN,DM\}.
\]
Under the first layout, identify the parties within $ACM$ and within
$DN$. Under the second layout, identify the parties within $BCN$ and
within $DM$. The two demand wires remain part of the local maps. In each
layout every supply wire in the corresponding group can be summed at a
single position, as indicated in Figure~\ref{fig:motif-layout}.


% Editorial omission of an illustrative figure.


\begin{lemma}[Six-party energy saving]
\label{lem:motif}
Take $t=1/2$ and give $AB$ and $CD$ masses
$r=t^2/(1-t^2)=1/3$ and $1$, respectively. There are rational positive
supply masses, of total $1+r$ in each group, for which
\begin{equation}
 r d_{AB}^2+d_{CD}^2\ \le\ \frac{31}{16}\sum_{e\in G_1\cup G_2}q_e d_e^2
 \label{eq:motif-six}
\end{equation}
holds for every six-point configuration in a $\mathrm{CAT}(0)$ space,
and hence for every assignment in a spectral wire pseudometric over $\mathbb F$.
\end{lemma}

\begin{proof}
Embed the wire distances into a $\operatorname{CAT}(0)$ space by Corollary~\ref{cor:spectral-geodesic} and Lemma~\ref{lem:catzero}. Write $XY$ for the distance between $X$ and $Y$. Give the pairs $(AC,BC)$, $(AM,BN)$, $(CM,CN)$ and $(DM,DN)$ effective coefficients $1/2,1/4,5/6,1$, respectively. Let $E$ be their weighted sum of squared lengths. Put $U=[C,M]_{1/3}$ and $V=[C,N]_{1/3}$. Equation~\eqref{eq:cur-cn} gives
\[
  \tfrac12 AC^2+\tfrac14 AM^2\ge\tfrac34 AU^2+\tfrac16 CM^2
\]
and the corresponding bound with $B,N,V$. Adding the remaining terms leaves $CM^2+CN^2+DM^2+DN^2$ in addition to $\tfrac34(AU^2+BV^2)$. At the midpoint $O$ of $MN$, the midpoint inequality bounds those four squares below by $2CO^2+2DO^2+MN^2\ge CD^2+MN^2$. Hence
\[
 E\ge\tfrac34(AU^2+BV^2)+CD^2+MN^2 .
\]
Convexity along geodesics gives $UV\le MN/3$, so $AB\le AU+BV+MN/3$. Weighted Cauchy--Schwarz gives
\[
 \tfrac13 AB^2\le\tfrac{25}{27}\bigl(\tfrac34(AU^2+BV^2)+MN^2\bigr).
\]
Thus $E$ bounds $\tfrac13 AB^2+CD^2$ from above. The effective coefficients sum to $31/12$ in each layout group. Dividing them by $31/16$ gives the required rational masses, of total $4/3$ per group.
\end{proof}

In the packing used in Section~\ref{sec:compressor}, a bit at each
endpoint marks it thick or thin. A fraction
$t^2=1/4$ of wires have two thick endpoints, and the packing distributes
their total mass over the remaining fraction $3/4$. This fixes $r=1/3$.
Fix any $31/16<H'<2$ in Theorem~\ref{thm:conversion} for this motif.
After choosing one term at each occurrence in a collection of motifs,
one common layout suffices to discard at least half of their total
supply energy. Indeed the surviving costs under the two layouts sum to
the entire supply energy. Before rounding widths to regularize the
degrees, the leading coefficient on the demand energy is therefore
$H'/2<1$. The construction in Section~\ref{sec:compressor} realizes both identifications consistently
when different occurrences share positions and preserves the strict
saving through rounding.
\section{Compression of a wire network}\label{sec:compressor}

Theorem~\ref{thm:compressor} gives an exact expansion of a wire network.
After $T$ blocks, the sum of the surviving widths in each term is
bounded by $2^{-T}$ times a quantity determined by the initial network.
A block has two stages. Several applications of the six-party identity
reduce the leading coefficient of a quadratic energy bound while
increasing the degree. The port graph then restores the original degree
with a multiplicative loss independent of the intermediate degree.
This independence permits iteration.

\subsection{Wire networks and local identities}

Let the base be a finite indexed directed multigraph with distinct endpoints
on every arc and with $d$ incoming and $d$ outgoing arcs at every vertex,
where $d$ is a power of two. Number the slots at each vertex separately in
each direction. Give each base arc $e:x\to y$ a displacement $\tau_e$ in a
binary vector space of lanes. The positions are the pairs $(x,g)$ of base
vertex and lane. The actual arcs over $e$ are
\[
       (x,g)\longrightarrow (y,g+\tau_e).
\]
Directions and slots index the wires. An actual arc of width $b_e$ equates
two strings of $b_e$ bits. The integer width, possibly zero, can vary from one
lane to another. With $P$ positions there are $Pd$ actual arcs. Give each
mass $p=1/(Pd)$ and put
\[
   \mathcal L=\sum_{e\,\text{actual}} b_e^2/p,
   \qquad \sum_{e\,\text{actual}} b_e\le\sqrt{\mathcal L}.
\]
To bound the energy of a concatenated wire, let $p_i$ and $b_i$ be
the masses and widths of its constituents. Cauchy--Schwarz gives
\[
        \frac{(\sum_i b_i)^2}{\sum_i p_i}\le\sum_i b_i^2/p_i.
\]
Thus a concatenated wire can inherit the total mass of its constituents
without increasing energy. Splitting it across integral slots incurs only
the rounding error estimated below.

Write $a_u$ for the string of incident half-wire indices at position $u$,
and $W((a_u)_u)$ for the product of the wire equalities. A substep gives an
identity of the form
\begin{equation}\label{eq:kernel}
 W((a_u)_u)=\sum_{\iota}\ \sum_{z}
   \prod_{v} K_{v,\iota}
       \big((a_u)_{\pi_\iota(u)=v},\,z|_v\big).
\end{equation}
Here $\iota$ is a list of admitted motif terms. The new topology, widths,
placement $\pi_\iota$ of old positions, and kernels are determined by this
list. The variable $z$ chooses one string per new actual arc, used at both
ends. A kernel consists of local ties, sums over local indices, and entries
of the maps in Theorem~\ref{thm:conversion}. After multiplying by arbitrary old
vertex factors and summing the $a_u$, every whole old factor belongs to the new
vertex $\pi_\iota(u)$. In this sum the two old half-wire indices can be
chosen separately, as their former equality has been expanded. Any uncontracted
indices of an old factor move with it. All identities below hold in every
commutative unital ring satisfying the fixed equations on the map entries.

\begin{theorem}[Iterated compression]\label{thm:compressor}
Suppose there are ports over $\mathbb F$ in the sense of
Definition~\ref{def:ports}. Use the family of connected loop-free symmetric regular
graphs of fixed degree and arbitrarily large power-of-two orders from
Section~\ref{sec:ports}. Its squared-distance Poincar\'e bound holds
with one constant for every spectral wire pseudometric over $\mathbb F$.
There exist a
power-of-two degree $d_0$, a positive integer $b_0$, a constant $C_2\ge2$
and finitely many motif templates, all fixed, with the following property.
Let a network of the stated form have $P_0$ positions of degree $d_0$,
maximum width $W_0$, and at least $Tb_0$ lane coordinates, for $T\ge1$.
It has an exact expansion by $T$ blocks of local identities
\eqref{eq:kernel}, ending at degree $d_0$. The same finite system of
parameter equations, satisfiable in $\mathbb F$, suffices for every network and $T$.
For every sequence of admitted terms,
\begin{equation}\label{eq:compress-bound}
  \sum_{e\,\mathrm{final}} b_e\le\sqrt{\mathcal L_T}
  \le 2^{-T}\big(\sqrt{\mathcal L_0}+O(P_0)\big).
\end{equation}
All intermediate wire and kernel fields have length at most
$C_2^T(W_0+1)$. At each final position, the recursive expansion of its
vertex tensor uses at most $C_2^T$ occurrences of vertices, map entries
and fields. This counts each initial factor and each map entry as an
unevaluated operation in the local contraction, with one occurrence
shared by all values of its hidden indices. For fixed arguments, a nonzero map entry contributes
a monomial of parameter degree at most $C_2^T(W_0+1)$.
\end{theorem}

\input{compressor-appendix}

\subsection{The energy estimate for one shrinking step}\label{subsec:shrink}

Use the masses of Lemma~\ref{lem:motif}. Write $t=1/2$, $r=1/3$, give the
two demands $AB,CD$ relative masses $r,1$, and give each supply type $e$ its
positive rational relative mass $q_e$. Each of the two supply groups has
total relative mass $1+r$. Fix $31/16<H'<2$ in Theorem~\ref{thm:conversion},
using the metric inequality of Lemma~\ref{lem:motif}, and choose $0<\delta_0<1$ such that
\[
                  \rho=\frac{H'}{2(1-\delta_0)}<1.
\]
With current degree $d$, the step works for every sufficiently large
power of two $d'$ as its output degree. For every term its energy satisfies
\begin{equation}\label{eq:compress-shrink}
        \mathcal L_{\mathrm{new}}\le \rho\mathcal L+O_{d,d'}(P^2).
\end{equation}
The fixed field and motif data may enter all implicit constants.

\paragraph{Packing and assignment.}
Apply Lemma~\ref{lem:compressor-slots}. There are $P'=P/K$ new
positions, with $K=2^{2d+1}d$. Each bundle concatenates $K/4$ thick--thick
wires between the packed endpoints $A,B$. Divide this payload into
$Q=3K/4$ integral pieces as evenly as possible and pair each piece
with one of the $Q$ demands $CD$ assigned to the bundle.

Give these two demands actual masses $rp,p$. Since $rQ=K/4$, the pieces
collectively inherit the mass of the bundled wires. A bundle of total width $b$
would cost $b^2/(Qrp)=b^2/((K/4)p)$ under perfectly equal division. By
Cauchy--Schwarz this is at most the original energy of its thick--thick
wires. Dividing an integer $b$ evenly into $q$ integers raises their sum of
squares above $b^2/q$ by at most $q/4$: their deviations are centered and
the variance of the two consecutive values is at most $1/4$.
Thus rounding costs at most $Q/(4rp)$ per bundle. There are $P'd$ bundles,
and each $CD$ wire keeps its energy. The total demand energy is therefore
at most $\mathcal L+O_d(P^2)$. Apply the wire identity of
Theorem~\ref{thm:conversion} with the masses of Lemma~\ref{lem:motif} to
each motif. After scaling masses by $p$, its additive cost is $O(1/p)$ per motif,
and there are $P'dQ=O_d(P)$ motifs. In every list of terms the total supply
energy is at most $H'\mathcal L+O_d(P^2)$.

\paragraph{Choosing the placement.}
The two layouts send each thin position to opposite ends of its bundle.
Because Equation~\eqref{eq:motif-conversion} is an identity before any
sites are identified, we can select a common layout after choosing the
motif terms.

The first supply group $G_1$ consists of the triangle on $A,C,M$ and the
edge $DN$. The second, $G_2$, consists of the triangle on $B,C,N$ and the
edge $DM$. Choose a single layout for the entire list of terms, retaining the
group of smaller total energy, with deterministic tie breaking. Each layout
first defines the placement $\pi$ of whole old positions. A thick position
goes to its packed site. A thin position goes to the $A$ of its assigned
bundle in layout~1 and to its $B$ in layout~2. This places every $C$ as
required. In every motif the $D$ map goes to $D^*=\pi(D)$, even if that old
position plays $C$ or $D$ elsewhere. Its old vertex factor thus moves only
once.

In layout~1 put the fresh $M$ at $A$ and the fresh $N$ at $D^*$. All edges
of $G_1$ become local. Of the surviving edges, $BC$ joins $B,A$, $BN$ joins
$B,D^*$, and $CN,DM$ join $A,D^*$. Layout~2 interchanges $A,B$ and $M,N$,
making $G_2$ local. The choice retains at most half of the total supply energy.

\paragraph{Regularizing the degree.}
Lemma~\ref{lem:compressor-slots} arranges that $A,B,D^*$ have distinct
bases in each motif and that each of these three roles occurs exactly $dQ$
times at every new position, in either layout.
The new arc mass is $p'=1/(P'd')=pKd/d'$. For every surviving type $e$ use
$m_e=\lfloor p q_e/p'\rfloor$ parallel arcs per occurrence, oriented by its
two roles in a fixed order. Take $d'$ large enough that $m_e>0$ and
$m_e p'\ge(1-\delta_0)p q_e$ for all types. Since each role occurs $dQ$
times, both the in-count and the out-count at every position are
\[
       dQ\sum_{e\,\text{surviving type}}m_e
       \le dQp(1+r)/p'=d'.
\]
Fill the common deficit with zero-width arcs as in
Lemma~\ref{lem:compressor-slots}. The slots commute with translations of the
remaining lane, so the new graph again has fixed displacements on its base
arcs. This symmetry concerns only the unweighted slots. Widths and map
entries may differ on every lane.

Split the surviving payload of width $b$ as evenly as possible among its
$m_e$ slots. The new cost is at most $b^2/(m_ep')+m_e/(4p')$.
The multiplicative loss relative to its supply energy is at most
$1/(1-\delta_0)$. For these fixed degrees there are $O_d(P)$ occurrences,
$m_e$ is bounded, and $1/p'=O_{d,d'}(P)$. Summing this estimate after the
choice of layout proves \eqref{eq:compress-shrink}.

\paragraph{Exact contraction and size of a kernel.}
Concatenating the thick--thick fields and slicing them identically at both
ends rewrites the old wire tensor as the product of the joint $AB,CD$
demands. Expand each with its selected motif identity, summing over lists
of selections. Put the maps for $C,D$ at the images of these old endpoints
and place all other maps by the selected layout. Every absorbed wire is
summed once inside its new position. The remaining wire bits are exactly the
bits allocated to parallel slots. This proves \eqref{eq:kernel}, including
when the layout depends on the term list. A new vertex receives boundedly many
old positions, since the number of slots and consumed coordinates is fixed.
It contains boundedly many maps and wire fields, by the role counts. If the
largest old width is $W$, each new or local field has length $O_{d,d'}(W+1)$,
by the bundle bound and Theorem~\ref{thm:conversion}. Concatenations and
slices are described by bounded lists of field identifiers, offsets and
lengths, even where two slice boundaries differ. Each map needs only the
compact term description and widths from the conversion theorem.
\subsection{Restoring the degree}\label{subsec:refresh}
Let $d_H$ and $\gamma_H$ be the degree and Poincar\'e constant of the
loop-free port graphs constructed in Section~\ref{sec:ports}.

Fix a power of two $d_0\ge \max\{4,d_H+1\}$.
Choose the last shrinking step so that its output degree $D$ is an order in the port family. Use the graph of that order on the $D$ ports. We give a step that
keeps all $P$ positions and returns their in- and out-degrees to $d_0$.
The arc mass changes from $p=1/(PD)$ to $p'=1/(Pd_0)=Dp/d_0$.

Consume $\log_2 D$ lane coordinates as a port index $j$, and include $j$ in
the new base label. For each old base vertex $x$ and remaining lane $g$, its
$D$ ports form one block. An old base arc $e:x\to y$ sends the $D$ ports in a source
block bijectively to those in a single destination block, by translation of both
the port and the remaining lane. Concatenate these $D$ wires, in source-port order,
into one external wire. Its endpoints are the port numbered by the outgoing slot of
$e$ at $x$ and the port numbered by its incoming slot at $y$. Record the same
ordering at the destination. Every port has exactly one outgoing and one incoming
external arc. If the old widths in this bundle were $b_1,\ldots,b_D$, its new
energy is bounded by
\[
 \frac{(\sum_j b_j)^2}{p'}\le \frac{D}{p'}\sum_j b_j^2
     =d_0\sum_j\frac{b_j^2}{p}.
\]

Inside each block, tie every old half-wire to its segment of the external
wire at the designated port. A tie within one port is local. The other ties
form two copies of all ordered pairs of distinct ports, one copy for the
outgoing slots and one for the incoming slots. Each tie has the width of its
old wire and gets mass $p$, so the total internal demand energy over all blocks
is at most $2\mathcal L$. Supply these ties with the directed expander arcs of
the block, each of mass $p'$. The demands test all port pairs, while the
supply tests only expander neighbors. The Poincar\'e inequality compares these
two averages with a constant independent of the number of ports. This is
exactly the metric comparison needed to convert the dense demand into
bounded-degree supply. In detail, for any spectral wire metric $d$ on the ports,
the common Poincar\'e bound gives
\begin{align*}
 \sum_{e\,\text{internal demand}}p\, d_e^2
   &= 2p\sum_{u,v}d(u,v)^2
    \le \frac{2\gamma_H Dp}{d_H}\sum_{e\,\text{expander arc}}d_e^2 \\
   &= \frac{2\gamma_H d_0}{d_H}\sum_{e\,\text{expander arc}}p'\,d_e^2 .
\end{align*}
Apply Theorem~\ref{thm:conversion} to the fixed $D$-port motif with all
masses divided by $p$, and use $H'=2\gamma_Hd_0/d_H+1$.
In each term the internal supply energy increases by a factor at most $H'$
up to an additive $O_D(1/p)$ per block. Summing the external and internal
estimates yields
\begin{equation}\label{eq:compress-refresh}
 \mathcal L_{\mathrm{new}}\le C_*\mathcal L+O_D(P^2),
 \qquad C_*=d_0+2H'.
\end{equation}
In particular $C_*$ is independent of $D$.

The internal identities and the external wires produce the
kernel in Equation~\eqref{eq:kernel}. Before the internal replacements,
summing an external wire with its segment ties restores the old wire equalities.
Each port now has $1+d_H$ arcs in each direction. Pad the deficit with
zero-width arcs from each port to the next in cyclic order. The external
arcs have distinct old base endpoints. Internal and padding arcs have distinct
port endpoints. All slots commute with translations of the remaining lane.
One old position stays at each new position, the number of local maps and
fields there is bounded in terms of $D$, and their lengths are
$O_D(W+1)$ when the incoming maximum width is $W$.

\subsection{Iteration}

\begin{proof}[Proof of Theorem~\ref{thm:compressor}]
Choose a positive integer $s_0$ so that $C_*\rho^{s_0}\le1/4$. Starting
from $d_0$, choose in succession $s_0$ sufficiently large power-of-two
output degrees for the shrinking steps. Choose the last to be an expander
order $D$ from this family, then restore the degree to
$d_0$. This order of choices is possible because $\rho$ and $C_*$ are
independent of the chosen degrees. The resulting block uses a fixed
number $b_0$ of lane coordinates. Each of its finitely many templates has
equations solvable in $\mathbb F$. Use disjoint names for their parameters and take
the union of the equations.

Every shrinking step divides the position count by its own $K$, while the
refresh preserves the count. Write $\kappa_0>2$ for the product of these
fixed $K$'s in a block. At the start of block $i$, numbered from zero,
there are $P_i=P_0/\kappa_0^i$ positions. Composing
\eqref{eq:compress-shrink} and \eqref{eq:compress-refresh} within a block
gives a fixed $c$ with
\[
          \mathcal L_{i+1}\le\tfrac14\mathcal L_i+cP_i^2 .
\]
Consequently
\[
   \mathcal L_T\le 4^{-T}
     \left(\mathcal L_0+4cP_0^2
                \sum_{i=0}^{T-1}(4/\kappa_0^2)^i\right).
\]
The geometric sum is uniformly bounded, giving \eqref{eq:compress-bound}.

At each of the fixed substeps a bound on the incoming field lengths changes
from $W$ to at most a constant times $W+1$, including the internal fields.
At each new position the numbers of immediate preimages, maps and field
descriptions are bounded by a constant. Iterating these two estimates down
the tree of local contractions gives the two bounds with $C_2$ sufficiently
large. This constant may also absorb a fixed number of description items per
occurrence and the constant in the linear bound on the degree of a map entry
from Theorem~\ref{thm:conversion}. Once final arc bits and term choices
are fixed, each final vertex sums only its own hidden indices, because each
substep sends each entire old vertex to one new position.
\end{proof}

\subsection{Computing the local contraction records}\label{app:descriptors}

The width estimate controls the summands in the expansion. To implement
their sum, we must also recover the contractions specified by a sequence
of terms. A \emph{descriptor} records the offset and length of a bit
field, together with the map or old vertex that uses it. These records
can be constructed before assigning any wire bits.

\begin{lemma}[Compressor descriptors]\label{lem:descriptors}
Under the hypotheses of Theorem~\ref{thm:compressor}, suppose that
$n\ge P=P_0$, $W_0\le n^{O(1)}$, and $T=O(\log n)$. There is a
nondeterministic RAM procedure with $O(\log n)$-bit words that takes
the initial slot records and widths as input. Its successful paths are
in bijection with the sequences of admissible term indices in the
theorem. Each such path constructs the following records at every level.
\begin{enumerate}
\item Each actual incident slot gives access to its arc width and both
endpoint slots. The base slots specify lane displacements.
\item Each new position has its full list of immediate old positions,
its local maps and their term data, and the slices, ties and summed fields
needed to evaluate its kernel in Equation~\eqref{eq:kernel}. Each position
has boundedly many records, accessible by their addresses.
\end{enumerate}
The procedure takes $O(P\operatorname{polylog} n)$ steps. With
comparisons, bit operations, shifts and ordinary word arithmetic, the
polylogarithmic exponent can be taken to be an absolute constant
(for instance, ten).
\end{lemma}

\begin{proof}
Store the lemma's records in arrays indexed by level, slot and position.
To reconstruct a map argument, identify each field by its source, offset
and length. For instance, when
concatenated fields are repartitioned, intersect the two lists of
intervals in the concatenation to obtain the slices to copy. Each local
map description contains its motif, role, fixed template, demand widths,
and the term index, including the numerical data for checking the total
payload widths on \emph{all} edges of that motif. At this stage a map
entry is represented by this description, without writing out the
selectors or the bit positions of its copies.

All degrees and templates in one compressor block are fixed. Theorem~\ref{thm:compressor}
bounds the number of immediate preimages, fields and maps per position
by a constant for each substep. Thus each position needs only a bounded
number of records per level. A record uses a bounded number of words.
Levels, slots, base vertices, tags and lanes admit direct product indexing
with $O(\log n)$ bits per index, since the numbers of sites and levels
are polynomial in $n$. Arrays for successive levels can occupy disjoint
address ranges. Blank space needs no initialization. 

To implement the slots and placements of Lemma~\ref{lem:compressor-slots},
a shrinking substep uses three scans.
\begin{enumerate}
\item In base order, compute the greedy colors and tag offsets. Looking
up a base vertex's neighbors, their incident displacements and earlier
assignments takes boundedly many accesses. Identify a bundle by its base
arc, source tag and remaining lane, and a motif by its $CD$ arc. Scan the
actual arcs, appending each thick--thick arc and each $CD$ arc to the
appropriate lists. Record the packing of every thick position and both
possible placements of each thin position. Order the bounded lists using
identifiers in which lane coordinates are written relative to the bundle's
source lane (by xor). This fixes the order consistently under translation
of the remaining lane. Divide each bundle into the prescribed pieces and
record their slices at both ends.
\item In a fixed order guess one term index per motif, using fixed-length
binary encodings and rejecting out-of-range or inadmissible indices. Each
admitted index occurs once, including indices whose terms may eventually
contribute zero. Accumulate the two possible surviving energies while
scanning, then choose the common layout with the prescribed tie rule. The
common mass factor cancels in this comparison. The other coefficients are
fixed rationals, and the scores are sums of polynomially many squared
integer widths of polynomial magnitude. The numerators and denominators
therefore use $O(\log n)$ bits after clearing a fixed denominator.
The grid, payload and rounding tests for an index likewise use fixed
rational tables and integer comparisons.
\item Store the chosen layout once for the level. Scan old positions,
appending each position exactly once to the preimage list at its selected
destination. For a thin position this takes the recorded endpoint of
its one assigned bundle, even when the position carries roles in several
motifs. Record each role incidence separately at this same destination.
Generate the surviving slots for each motif together with
endpoints, role edge and subdivision number. Include zero-width slots and
generate the prescribed dummy slots. At each position order the outgoing
and incoming slots separately by origin, using the motif's identity
relative to that position's lane, its role edge and subdivision number,
or the dummy slot's identity. These keys use only the combinatorial
placement in the chosen layout. In particular they are independent of
payload widths. Give each generated actual arc its own record indexed
by its absolute motif ID, role edge and subdivision number (or its dummy
ID). Both endpoint lists carry that address. After each bounded list is
sorted, write its local slot number back to the corresponding endpoint
field in the arc record. A scan of the outgoing slots at lane zero can
now read the destination slot and xor displacement directly. Since
translation preserves the relative keys, these are the base records for
every lane. The individual widths and term indices stay in the actual
records. Split surviving widths as prescribed.
Record the argument slices for the maps, the internal identifications of
absorbed supplies and the slices occupying the new incident arcs.
\end{enumerate}
For a degree refresh the port blocks have fixed size $D$. For each old
base arc and source block, list its $D$ actual arcs in source-port order,
using the known translated port at the other end. This describes both
ends of the concatenated external arc. In each block guess one internal
motif term as above and record its demand, supply and map arguments.
The new incoming and outgoing lists each have one external slot. Reserve
the other slots for the indicated internal and dummy arcs. These direct
correspondences specify the slices and ties at every port.

Each scan visits $O(P)$ records, since the number of positions never
increases. Each visit needs at most polylogarithmically many elementary
operations and probes. There are $O(\log n)$ substeps, so the entire
descriptor phase costs $O(P\operatorname{polylog} n)$. In particular this phase pays for numerical widths
through their binary lengths.
\end{proof}
~~~

## G / compressor-appendix.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/compressor-appendix.tex, Git blob b620f6fd0aa309d4bb79ee3386833265c99819ed. Complete placement and packing proof.

~~~latex
\subsection{Packing the wire slots}\label{sec:compressor-slots}

The lane coordinates let us replicate each base arc in four endpoint
configurations. We mark each endpoint thick or thin. Thick positions
will merge at a fixed site, whereas thin positions can move to either
end of an assigned bundle. The assignment must move every occurrence
of a given old position together, even when several motifs meet there.

\begin{lemma}[Packing and placement of the wire slots]\label{lem:compressor-slots}
For every power-of-two base degree $d$, set $K=2^{2d+1}d$ and $Q=3K/4$.
There is a packing of the wire slots, using a number of lane coordinates
depending only on $d$, with the following properties. The two layouts
refer to the six-party motif of Figure~\ref{fig:motif-layout}.
\begin{enumerate}
\item Packing the thick positions produces $P/K$ new positions. Each
thick--thick bundle contains $K/4$ wires. Every other wire is assigned as
a $CD$ demand to the bundle of its thin $C$ endpoint, with $Q$ per bundle.
\item In either layout the packed endpoints $A,B$ and the image $D^*=\pi(D)$
of the other endpoint have distinct bases in each motif. Each of these three
roles occurs exactly $dQ$ times per new position.
\item Translation of the remaining lane preserves the assignments, placements
and unweighted slots. Slot numbering, including the padding
to any admissible output degree, can be determined by bounded-size lists
per new position after a greedy scan of the base.
\end{enumerate}
\end{lemma}

\begin{proof}
Color the base greedily with $2d+1$ colors. Split the lane as $(U,h,j,g)$,
where $U$ supplies one bit per color, $h$ indexes the $d$ outgoing slots,
$j$ is a tag in a binary vector space of order $J$, and $g$ is the remaining
lane. We determine $J$ below. Call a position thick if its base-color bit
is $1$, and thin otherwise. Put $K=2^{2d+1}d$, the number of pairs $(U,h)$.
Along a base arc the endpoint colors are distinct. With $j,g$ fixed at one
end, the two tests use distinct free bits of $U$ up to translations. They
are independent fair tests, also independently of $h$ at either end.
Consequently each of the four endpoint types occurs $K/4$ times.

For given $x,j,g$ pack all thick positions there into a single new one, with
base $(x,j)$ and lane $g$. There are $P'=P/K$ new positions. Bundle the
$K/4$ thick--thick wires of each base arc for fixed $j,g$ at its source by
concatenating them in a common order at both endpoints. The bundle runs from
a packed position $A$ to a packed position $B$. For every wire outside these
bundles take a thin end for $C$ and the other end for $D$. For thin--thin
wires take their old direction $C\to D$.

Give each base vertex $x$ an offset $\phi_x$. Assign each thin position
$(x,U,h,j,g)$ to the bundle with source $(x,j+\phi_x,g)$ on the outgoing
base arc of slot $h$. Every wire on which this position plays $C$ is assigned
to that bundle. To count these wires, fix the base, tag, remaining lane and
slot of its possible thin sources. The $d$ originally outgoing incidences
contribute $d(K/(2d))=K/2$ wires. The $d$ originally incoming incidences
from thick positions contribute $d(K/(4d))=K/4$. Thus each bundle receives
$Q=3K/4$ wires of type $CD$.

We choose the tag offsets so that $A,B,D^*$ always have distinct new bases.
The old base components of $A,B$ differ. For $C$ at old base $x$ with tag
$j$, and $D$ at $y$ on their old base arc $e$ in either direction, write
\[
          L_x=\{0\}\cup\{\tau_f^{\mathrm{tag}}: f\text{ leaves }x\}.
\]
The new tags of $A,B$ lie in $j+\phi_x+L_x$. A thick $D$ has tag
$j+\tau_e^{\mathrm{tag}}$, and the possible tags of a thin $D^*$ belong to
$j+\tau_e^{\mathrm{tag}}+\phi_y+L_y$. It suffices to impose, at each
incidence,
\[
      \phi_x\notin\tau_e^{\mathrm{tag}}+L_x,
      \qquad \phi_x+\phi_y\notin\tau_e^{\mathrm{tag}}+L_x+L_y.
\]
Greedy choice in base order works for any power of two
$J>2d((d+1)+(d+1)^2)$, since at most this many values are forbidden by
the first conditions and the already assigned neighbors in the second.

\smallskip
At every new position $d$ bundles start and $d$ bundles end, so each of
the roles $A,B$ occurs $dQ$ times. The thick $D$ endpoints at this position
contribute $2dK/4=dK/2$ occurrences of $D^*$, counting mixed arcs in both
old directions. The others are destinations of thin--thin arcs. Before
placement, for fixed base, tag and remaining lane, these contribute $dK/4$
occurrences, equally distributed over destination slot labels $h$. In
layout~1 the tag offset permutes the tags. In layout~2 each position receives
through $d$ incoming base arcs, each with a uniquely determined source tag
and remaining lane, and each carrying $K/4$ such occurrences from its
assigned slot. In both cases the total number of $D^*$ roles is
$dK/2+dK/4=dQ$.

Packing, bundle assignment, and the two choices of image depend on the
remaining lane only through the given arc translations. The common choice
of layout therefore preserves this symmetry. For each surviving occurrence,
label its parallel slots by the original $CD$ arc identifying the motif, the
surviving type, and the subdivision index. When a common deficit remains
after the degree calculation in Subsection~\ref{subsec:shrink}, add zero-width
arcs to the next tag cyclically at the same old base and remaining lane.
Number incoming and outgoing slots separately by sorting their identifiers,
expressing all remaining lanes relative to the vertex being numbered. The
orders transport from lane zero to every lane. Widths and term-dependent map
entries need never enter these identifiers.

Greedy coloring and tag selection each use a single scan of the indexed
base and bounded neighbor lists. After these scans the assignments can be
formed by listing the actual arcs by bundle, and the regularized slots by
their endpoints. The list sizes depend only on the fixed degrees and lane
coordinates consumed in this step. The sorted identifiers give the two
corresponding slot numbers and their displacement, as required for the
next level of the construction.
\end{proof}
~~~

## P / compress.tex

OpenAI, preprints/Circuit-lower-bounds-for-Hamiltonian-cycle-and-the-permanent-November-24-2026/build/compress.tex, Git blob 4b0887da937eaf4c512232571fae596093a6ee12. Complete expansion, enumerable contraction and routing argument.

~~~latex
\section{Compression from ports}\label{sec:compr}

Assume throughout this section that the characteristic admits ports.
To enumerate a contraction directly, one assigns a bit to every wire
and multiplies the entries of all tensors. The compressor replaces
this enumeration by a sum in which fewer bits are shared between
positions. Once those bits are fixed, each position can sum over its
own indices independently. We describe the expansion from
\cite[Theorems~\ext{G}{thm:compressor} and~\ext{G}{thm:conversion},
Lemma~\ext{G}{lem:descriptors}]{oaG} and use it to evaluate a closed
network in subexponential time. Section~\ref{sec:gap} will use the
separate bounds on shared and private work.

\subsection{An exact expansion with fewer shared bits}
The compressor applies to a network indexed by base vertices and lanes.
The base is a directed loop-free multigraph of fixed in- and out-degree
$d$, a power of two with $d\ge4$. Each pair $(u,x)$ of a base vertex
and a lane $x\in\{0,1\}^h$ is a \emph{position}. A base arc
$e:u\to v$ with mask $\tau_e$ induces an actual arc from $(u,x)$ to
$(v,x+\tau_e)$ in each lane, with addition bitwise modulo two. Each
actual arc carries a string of bits shared by its two endpoints. The
width may vary from lane to lane and may be zero. Write $P$ for the
number of positions, $D$ for the maximum width and $A$ for the
contraction of the factors at the positions. A factor may also
depend on parameters, including advice bits.

One compression block applies a fixed sequence of local tensor
identities. Each identity is a sum, and a \emph{history} records one
admissible summand at every application. Given $t\ge1$ blocks with
$b_0t\le h$, where $b_0$ is fixed, the final network in a history
$\lambda$ has bounded degree and at most $P$ positions. Let $z$ assign
bits to its arcs. At a final position $v$, let $y_v$ assign the bits
private to that position and let $K_{\lambda,v}$ be the product of its
local entries and equality tests. The compressor gives the identity
\[
 A=\sum_{\lambda}\sum_z
          \prod_v\left(\sum_{y_v}
             K_{\lambda,v}(z|_v,y_v)\right),
\]
where $z|_v$ denotes the strings on arcs incident to $v$. Each original
factor moves intact to one final position, with all of its arguments.
Three estimates determine the cost of implementing this identity.

\begin{enumerate}
\item The shared string $z$ has length $O(2^{-t}P(D+1))$ in every
history. The construction controls a quadratic energy of the wire
widths. Local tensor identities derived from a comparison of
Dirichlet energies reduce the leading coefficient in the energy bound
while raising the degree. The ports restore bounded degree at a fixed
cost. Grouping enough shrinking steps into each block gives the factor
$2^{-t}$ in the bound on total width.

\item The local contraction at any one $v$ uses at most $C^t$
occurrences of bit strings, map entries and original factors across
all its recursive levels, for a fixed $C>1$. Each string has length
at most $C^t(D+1)$. Here an original factor counts as one occurrence.
The maps use a fixed finite list of scalars satisfying fixed integer
polynomial equations. Applying the companion theorems over the prime
field gives a solution there, and the identities hold over any
commutative unital ring containing a solution. For given bit strings
a map entry is zero if a forbidden pattern occurs. Otherwise it is a
monomial in these scalars, of degree at most $C^t(D+1)$. The tests
and exponents take binary work polynomial in the string widths, with
an absolute exponent.

\item A history specifies the layout before any private bits are
read. Suppose $P,D\le V$ and $t=O(\log V)$. Given the initial arc
widths, indexed endpoint slots and masks, a nondeterministic
random-access procedure uses $O(P\operatorname{polylog}(2V))$ bit
operations and guessed bits per path. Its successful paths are in
bijection with admissible histories. At each level the procedure
records the fields, slices and ties needed for the local sums, and a
description of each map test of length $\operatorname{polylog}(2V)$.
These records can be accessed by address.
\end{enumerate}

\subsection{Evaluating closed networks}

The application to integer multiplication in \cite{oaT} uses both the cost
of evaluating the contraction and the cost of constructing its formula.

\begin{theorem}[Enumerable contraction]\label{thm:closed}
Fix a characteristic admitting ports. There is an $\eta>0$ such that
every closed bit network of volume at most $s$ and maximum degree at
most $s^{1/20}$, over any field of this characteristic, has a formula
for its contraction using $2^{O(s^{1-\eta})}$ arithmetic operations
and local entry queries. The topology determines the formula, which
can be generated deterministically in $2^{O(s^{1-\eta})}$ bit
operations and uses fixed prime-field constants. In characteristic
zero, for integer entries of polynomial bit length, evaluation takes
the same number of bit operations plus the time for the queries.
\end{theorem}

\begin{proof}
\emph{Routing the input network.}
For large $s$ choose the least $h$ with $P=h2^h\ge s^{9/10}$.
Then $P=\Theta(s^{9/10})$ and $h=\Theta(\log s)$. Use the positions
$(j,x)$ with $j\in\mathbb Z/h\mathbb Z$ and $x\in\{0,1\}^h$. From
$(j,x)$ put arcs to $(j+1,x)$ and $(j+1,x+e_j)$, where $e_j$ is the $j$th unit
vector. Add parallel slots of width zero to reach base degree $d$.

Place the input tensors in the $2^h$ lanes at $j=0$. Give each
tensor weight one plus its degree and always place the next tensor
in a lane of least current load. The total weight is at most $2s$.
The largest final load $q$ is at most the average load plus the
largest weight, so
\[
       q\le 2s/2^h+1+s^{1/20}=O(s^{1/10}\log s).
\]
Orient the original wires arbitrarily, and form a bipartite graph on
two copies of the lanes by joining the source lane of each wire to
its destination lane. Both sides have maximum degree at most $q$.
Pair their unused slots to complete a $q$-regular bipartite
multigraph. Successive perfect matchings split its edges into $q$
permutations.

Each permutation can be routed in one forward and one reverse trip
around the cursor $j$, by the Bene\v s construction~\cite{Benes}.
To see this, pair lanes differing in the first bit at both source
and destination. The permutation induces a degree-two bipartite
multigraph on the pairs. Color each cycle alternately. The switches
at the first and last stages can now send one color to each half, and
within either half the colored edges specify a permutation on one
fewer bit. Induction routes each permutation with at most one item
per lane on each pass. The reverse trip traverses the same arcs in
the opposite direction, which is allowed for bit equalities. Give
the two passes distinct bit slots and discard the padding wires.
Each original bit assignment has a unique extension satisfying the
equalities along these routes.
The maximum actual arc width $D$ and the local load, including the
tensors placed at $j=0$ and the routing ties, are both $O(q)$.

\emph{Enumerating the compressed expansion.}
Take $T=\lfloor a\log_2 s\rfloor$ blocks for a sufficiently small
fixed $a>0$ so that $b_0T\le h$. For each history the shared arc
bits have total length
\[
        O(2^{-T}P(D+1))=O(s^{1-a}\log s).
\]
The choices of local summands themselves use $O(PT\log s)$ bits.
There are $O(P)$ applications per substep and boundedly many substeps
per block. Each choice is a fixed number of integers of polynomial
magnitude. To enumerate histories deterministically, enumerate
the guess tapes for the layout procedure, checking choices in order
and requiring unused bits to be zero. This counts each successful
path once and costs $O(P\operatorname{polylog}s)$ bits per tape.

At each final position there are at most $C^T$ private fields of
length at most $C^T(D+1)$. Expanding the initial position factors
adds at most $O(C^T q)$ original tensors and routing ties. By taking
$a$ small enough, the number of private bits and the local load are
each at most $s^{1/2}$ for large $s$. For each pair $(\lambda,z)$
form the sum over $y_v$ at each position, then multiply the local
sums. For fixed $(\lambda,z)$ the work is the sum of the local enumeration
costs, because each private assignment occurs in a single local sum.
Thus the logarithm of the total work is at most
\[
 O\bigl(P\operatorname{polylog}s+s^{1-a}\log s+s^{1/2}\bigr)
       =O(s^{1-\eta})
\]
for some $\eta>0$. The layout records generate all indices, tests
and scalar monomials in polynomial work per assignment. This gives
both the arithmetic formula and the claimed bound for listing it.

In characteristic zero choose the fixed map scalars in $\mathbb Q$
and let $m$ clear their denominators. Each term of a local sum is
a product of polynomially many input entries and scalars with
polynomial total degree. Its denominator divides $m^{\mathrm{poly}(s)}$.
A local sum has at most $2^{\mathrm{poly}(s)}$ terms. Multiplying
polynomially many such sums and then adding over the outer enumeration
still gives polynomial bit length, since a common denominator at each
stage is a polynomial power of $m$. This argument also works for
inputs in a fixed ring of cyclotomic integers by using a fixed
integral basis. The polynomial cost per operation is absorbed into
the displayed bound.
\end{proof}
~~~
