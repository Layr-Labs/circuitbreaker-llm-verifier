# Source evidence 3. Metric entropy and product closure

Excerpts from OpenAI's mathematics preprints at commit 92b516ef8e2abf66c3db7aa31438d3010bc7956f. The dependency map, bibliography and Apache 2.0 license are in source-5.md. The authors' result labels are preserved.

## P / entropy.tex

OpenAI, preprints/Circuit-lower-bounds-for-Hamiltonian-cycle-and-the-permanent-November-24-2026/build/entropy.tex, Git blob 267da065e125912693ffa5dce25b6b50fd25a91d. Complete mathematical text. The illustration is omitted.

~~~latex
\section{A logarithmic lower bound for the Poincaré constant}\label{sec:poincare}

The price comparison in Section~\ref{sec:comparison} requires a
single placement of the row sites that works for every balanced cut.
To construct it, we shall look for placements that make two random
changes of state have very different mean squared lengths. The
simplest comparison is on a regular graph. One step follows an edge,
and the other goes to an independent uniform vertex. If a class of
Hadamard spaces has no ports, the constant that bounds the second
energy by the first is at least of order \(\log n\) on \emph{every}
\(n\)-vertex graph of bounded degree. This bound will allow us to
distinguish two walks by the information they generate.

For a finite Markov kernel \(P\) with stationary law \(\pi\) and a
placement \(x=(x_i)\) of its states, put
\[
 \mathcal E_P(x)=\sum_{i,j}\pi_iP_{ij}d(x_i,x_j)^2,\qquad
 \mathcal V_\pi(x)=\sum_{i,j}\pi_i\pi_j d(x_i,x_j)^2 .
\]
For a pair law \(D\), the weights in \(\mathcal E_D(x)\) are its joint
probabilities. Write \(\mathcal V\) when \(\pi\) is uniform. States
of zero stationary mass may be discarded. A \(d\)-regular multigraph
has a symmetric matrix of nonnegative integral arc counts with every
row sum \(d\). A holding arc counts once. Dividing the counts by
\(d\) gives the walk kernel with uniform stationary law.

Suppose one large graph had too small a Poincar\'e constant. Placing
a copy of it at each vertex of a preceding graph would multiply the
order at every iteration. To keep the Poincar\'e constants bounded,
we average the walk on the preceding graph before inserting the
copies. The estimate for this averaging is a form of nonlinear
spectral calculus \cite{MendelNaor,MendelNaorBarycentric}. We give
the barycentric argument for Hadamard spaces
\cite{BH,Sturm}.

\begin{lemma}[Smoothing a walk]\label{lem:power}
Let \(P\) be a finite Markov kernel with stationary law \(\pi\), and let
\(t\ge2\) be even. Put
\[
 A=\tfrac12(I+P),\qquad S=\tfrac1t\sum_{h=1}^t A^h .
\]
For each placement \(x\) in a Hadamard space \(X\) there is a placement
\(y\) in \(X\) such that
\begin{equation}\label{eq:hadamard-smoothing}
 \sum_i\pi_i d(x_i,y_i)^2+t\mathcal E_P(y)
 \le 16\mathcal E_S(x).
\end{equation}
If \(D\) has both marginals \(\pi\) and
\(\mathcal E_D\le\beta\mathcal E_P\) for all placements in \(X\), then
\(\mathcal E_D\le C(1+\beta/t)\mathcal E_S\) for all placements in \(X\),
with an absolute constant \(C\).
\end{lemma}

\begin{theorem}[Logarithmic Poincaré obstruction]\label{thm:loggraph}
Let \(P\) be the walk on a \(d\)-regular multigraph with \(n\ge2\)
vertices, \(d\ge1\), and uniform stationary law \(\pi\). Suppose
that \(\mathcal V_\pi(x)\le\beta\mathcal E_P(x)\) for every
placement \(x\) in every member of a class of Hadamard spaces
without ports. Then
\[
                 \beta\ \ge\ c\,\frac{\log n}{\log(d+1)}
\]
for an absolute constant \(c>0\).
\end{theorem}

\begin{proof}[Proof of Lemma~\ref{lem:power}]
The barycenter \(b\) of a finite probability measure \(\mu\) in \(X\)
satisfies the variance inequality
\(\int d(z,w)^2\,d\mu(w)\ge d(z,b)^2+\int d(b,w)^2\,d\mu(w)\).
Indeed the midpoint inequality makes any minimizing sequence Cauchy,
and squared-distance convexity along the segment from the minimizer
to \(z\) gives the displayed inequality.
Starting with \(x^0=x\), take \(x_i^{g+1}\) to be the barycenter of the
\(x_j^g\) with weights \(A_{ij}\). Put
\[
 q_g=\sum_{i,j}\pi_iA_{ij}d(x_i^{g+1},x_j^g)^2,\qquad
 u_g=\sum_i\pi_i d(x_i,x_i^g)^2 .
\]
Run the stationary \(A\)-walk for \(h\) steps. Conditional on the first
\(h-1\) steps, apply the variance inequality to the last point with
the starting point as \(z\). Replace the last point by its conditional
barycenter and repeat. At the \(g\)-th replacement, starting at
\(g=0\), stationarity makes the mean variance term equal to \(q_g\).
Consequently
\[
 \mathcal E_{A^h}(x)\ge u_h+\sum_{g=0}^{h-1}q_g .
\]
On averaging over \(1\le h\le t\), each \(q_g\) with \(g<t/2\)
occurs in at least \(t/2\) summands. Since \(u_0=0\),
\[
 \sum_{0\le g<t/2}\bigl(u_g+(t/2)q_g\bigr)\le t\mathcal E_S(x).
\]
Choose a summand of size at most \(2\mathcal E_S(x)\).
Comparison along \(x_i^g,x_i^{g+1},x_j^g\) gives
\(\mathcal E_P(x^g)\le4q_g\), because \(A\) averages a holding step
and a \(P\)-step. Taking \(y=x^g\) proves \eqref{eq:hadamard-smoothing}.
Finally
\(\mathcal E_D(x)\le3(2u_g+\mathcal E_D(y))\) by the two marginal
identities. Apply the assumed comparison at \(y\).
\end{proof}

\begin{proof}[Proof of Theorem~\ref{thm:loggraph}]
A class without ports contains a space with two distinct
points. Choose one of the points independently and uniformly at
each vertex, and average the asserted inequality. The expected pair
energy is \((1-1/n)\) times half their squared distance, and the edge
energy is at most half their squared distance. Thus \(\beta\ge1/2\).
Write \(b=\max(1,\beta)\) and call the given multigraph \(H\).

Fix an even \(t\ge2\) with \((2d+1)^t\le n/2\). Given a previous
graph of degree \(2d\) with walk \(K\) and Poincaré constant at most
\(\gamma\), form the averaged kernel \(S\) of Lemma~\ref{lem:power}
from \(K\). The new pair comparison for \(S\) has constant at most
\(C(1+\gamma/t)\) throughout the class.
We need a symmetric approximation to \(S\) whose probabilities have
denominator \(n\). In each row set \(m_{uv}=\lceil n S_{uv}/2\rceil\).
A row of \(S\) has support at most \((2d+1)^t\), since holds pad shorter
walks to length \(t\). Hence its integer row sum is at most
\(n/2+(2d+1)^t\le n\). Pad each diagonal entry to reach \(n\).
The counts remain symmetric and \(R_{uv}=m_{uv}/n\ge S_{uv}/2\).
It follows that the Poincaré constant of \(R\) is at most
\(\rho=2C(1+\gamma/t)\).

Replace each old vertex \(u\) by a copy of \(H\). Partition its \(n\)
vertices into sets of sizes \(m_{uv}\), and match the set indexed by
\(v\ne u\) bijectively to the set indexed by \(u\) in the copy at \(v\).
Fix the vertices in the set indexed by \(u\) itself. This gives an
involution whose projected step from a uniform vertex of the copy
has distribution \(R_{u\bullet}\).
The new graph has the \(d\) internal arcs of \(H\) at each vertex and
\(d\) copies of its involution arc. Its mean edge energy is
\((E_{\rm in}+E_{\rm match})/2\), with the two terms the corresponding
mean energies.

For a placement \(x\) in any one of the test spaces, let \(z_u\)
be the barycenter of the points in the copy at \(u\). Let \(v_0\) be
the mean squared distance to one's own barycenter. Minimization against
a uniformly chosen point in each copy and the comparison on \(H\)
give \(v_0\le b E_{\rm in}\). Apply the triangle inequality first from
one barycenter to the next via a matched pair, then between independent
vertices via their two barycenters. Uniformity of both endpoints of an
involution step gives
\[
 \mathcal E_R(z)\le3(2v_0+E_{\rm match}),\qquad
 \mathcal V(x)\le3(2v_0+\mathcal V(z)).
\]
Here each pair energy uses the uniform law on the indicated graph.
Since \(\mathcal V(z)\le\rho\mathcal E_R(z)\), the new graph has a
Poincaré bound \(B b(1+\gamma/t)\) for an absolute \(B\ge1\).

If \(Bb/t\le1/2\), iteration starting from one vertex with \(2d\)
loops gives graphs on \(n^k\) vertices of fixed degree \(2d\) with
bounds at most \(2Bb\). These graphs would be ports for the class.
Hence every admissible even \(t\) satisfies \(t<2Bb\). If
\(\log(n/2)/\log(2d+1)\ge4\), the largest admissible even \(t\)
is at least half this ratio. In the remaining range
\(\log n/\log(d+1)\) is bounded absolutely and \(\beta\ge1/2\)
suffices. Since \(b\le2\beta\), this proves the assertion.
\end{proof}

We shall apply the bound to walks with unequal transition
probabilities. Only the number of possible destinations matters.

\begin{corollary}[Symmetric kernels with sparse rows]\label{cor:support}
Suppose \(P\) is a symmetric stochastic matrix on \(n\ge2\) states
with at most \(d\) positive entries per row, where \(d\ge1\) is
an integer, and let \(\pi\) be uniform. If
\(\mathcal V_\pi(x)\le\beta\mathcal E_P(x)\) for every placement \(x\)
in every member of a class of Hadamard spaces without ports, then
\(\log n\le C\beta\log(d+1)\) for an absolute constant \(C\).
\end{corollary}
\begin{proof}
Take \(\lceil dP_{ij}\rceil\) arcs in each position and add holding
arcs to make every row sum \(2d\). The walk on the resulting symmetric
regular multigraph dominates \(P/2\) entrywise, so its Poincaré
constant is at most \(2\beta\). Apply Theorem~\ref{thm:loggraph}.
\end{proof}
\section{Observable entropy}\label{sec:entropy}

A graph of small degree gives a walker few choices at each step.
A walk with many possible destinations can also transmit little
information if a few choices carry most of the probability, or if
many destinations share the same observable label. We compare
two walks on the same state space. One, the \emph{demand} walk,
specifies the pairs of states we want far apart in a placement.
The other, the \emph{physical} walk, specifies the moves we want
short. A uniform comparison of their squared lengths bounds the
amount of information one demand step can generate about the
destination label, given the full starting state. The bound measures
the physical walk by the number of bits per observed label along a
long stationary trajectory.

Given a labeling $L$ of a finite set $V$ and a stationary $P$-walk
$(X_t)_{t\ge0}$ on $V$, put
\[
 H_j=H\big(L(X_0),\ldots,L(X_{j-1})\big),\qquad
 a=\lim_{j\to\infty}\frac{H_j}{j}=\inf_{j\ge1}\frac{H_j}{j}.
\]
Here \(H\) denotes Shannon entropy in bits, and all logarithms in this
section are in base two \cite{ct}. The equality follows from
stationarity and subadditivity of block entropy.

\begin{theorem}[Observable entropy]\label{thm:entropy}
Let $P,D$ be Markov kernels on a finite set with common stationary
law $\pi$. Suppose that $\beta\ge1$ and that
$\mathcal E_D(f)\le\beta\mathcal E_P(f)$ for every placement $f$ in
a class of Hadamard spaces without ports. For any labeling $L$, let $a$ be its entropy rate
for the stationary $P$-walk and draw $(X,Y)$ from the stationary
$D$-step. Then
\[
 H(L(Y)\mid X)\le C\beta(1+a)
\]
for an absolute constant $C$.
\end{theorem}

In particular, if the demand information exceeds \(C\beta(1+a)\)
for some \(\beta\ge1\), one of the target spaces admits a placement with
\(\mathcal E_D(f)>\beta\mathcal E_P(f)\). Section~\ref{sec:product}
uses the identity labeling on a walk of words modulo cyclic rotation.
The full state is then observable and the entropy rate is its
one-step conditional entropy. In Section~\ref{sec:lamps} the label
records only an array of lamps. The physical walk also keeps auxiliary
records, so the cost in Theorem~\ref{thm:entropy} is the entropy rate
of the observed arrays over time.

The additive constant accounts even for deterministic motion. If \(P\)
alternates between two states and \(D\) draws a uniform state
independently, the identity labeling gives \(H_j=1\) for every \(j\).
Here \(a=0\), yet \(H(Y\mid X)=1\) and
\(\mathcal E_D=\mathcal E_P/2\).

\begin{proof}
We may assume that the stationary $D$-step is exchangeable. Indeed,
replace its pair law by the average of that law and its transpose.
The squared energies are unchanged, and concavity of conditional
entropy in the joint law shows that the new $H(L(Y)\mid X)$ is at
least half the original value. The stationary $P$-walk is unchanged.

Write $h=H(L(Y)\mid X)$ for this symmetrized step and assume $h>0$.
We use \(r\) copies of the original state, called the reservoir, to
write entries in an array of \(m\) lamps. At each writing step the
value added to one lamp is a function of the reservoir labels.
From a set of starts of mass at least one half, a demand move in all
\(r\) cells yields an increment whose law dominates a fixed fraction
of the uniform law on \(\lfloor rh/4\rfloor\)-bit words. We shall use the assumed
comparison to charge a refresh of all the lamps to \(O(mr\beta)\)
physical steps. Corollary~\ref{cor:support} will then require the net
change over those steps to carry order \(mrh\) bits. A description
of the same change using only the observed reservoir labels and the
writing schedule gives the upper bound in terms of \(a\).

\smallskip\noindent
\emph{Extracting lamp values.}
Fix \(m\ge1\). Take \(r\) sufficiently large that \(r\ge2\) and
\(k=\lfloor rh/4\rfloor\ge1\).
Set $A=\mathbb F_2^k$ and $q=|A|$. Write $L^r$ for the coordinatewise
labeling. There is a function $b$ from $r$-tuples of labels to $A$
and a set $B\subseteq V^r$ of $\pi^r$-mass at least $1/2$ such that
for $U'\sim D^{\otimes r}(u,\cdot)$ one has
\begin{equation}\label{eq:lamp-hash}
 \Pr\{b(L^r(U'))=z\mid U=u\}\ge \frac1{4q}
 \qquad(u\in B,\ z\in A).
\end{equation}
To see this, draw $r$ independent stationary $D$-steps. The negative
logarithm of the conditional probability of their destination label
tuple, given the starting reservoir $U$, is the sum of $r$ independent
random variables of mean $h$. The law of large numbers implies that,
for a set of starts of mass $1-o(1)$, atoms of conditional probability
at most $2^{-3rh/4}$ carry at least half the conditional mass. Assign
each label tuple a uniformly random color in $A$, independently. At
any such start, the mass of these small atoms assigned to a given
color has mean at least $1/(2q)$ and variance at most $2^{-3rh/4}/q$.
Chebyshev's inequality and a union bound give failure probability at
most $16q^2 2^{-3rh/4}=o(1)$ for the lower bound $1/(4q)$ in all
colors. Averaging over the coloring gives~\eqref{eq:lamp-hash}.

\smallskip\noindent
\emph{Comparing two walks on lamps and reservoirs.}
The physical kernel $G$ has stationary law $\pi^r$ times the uniform
law on $A^m$. With probability $1-1/r$ it advances a uniformly chosen
reservoir cell by $P$. With the remaining probability it chooses a
uniform lamp $J$ and adds $b(L^r(U))$ to that lamp. The first operation
preserves $\pi^r$ by stationarity of \(P\), and the second permutes the
uniform array for every fixed reservoir.
The demand kernel $D_*$ keeps the reservoir fixed. It replaces all
lamps by independent uniform values if $U\in B$ and holds otherwise.
For every placement \(F\) of this enlarged state space in a member of
the class, we claim that
\begin{equation}\label{eq:lamp-compare}
 \mathcal E_{D_*}(F)\le Cmr\beta\,\mathcal E_G(F).
\end{equation}

We first compare a refresh of just one lamp with the three moves in
Figure~\ref{fig:lamp}. Write \(s\in A^m\) for the lamp array and
\(e_Jz\) for the array equal to \(z\) in position \(J\) and zero
elsewhere. At a good starting reservoir \(u\in B\), draw
\(u'\sim D^{\otimes r}(u,\cdot)\). Add \(z=b(L^r(u'))\) to lamp \(J\) and return the reservoir to \(u\).
The endpoint is \((u,s+e_Jz)\). By \eqref{eq:lamp-hash}, every
increment has probability at least \(1/(4q)\). Thus the endpoint
law dominates one quarter of the uniform replacement law.


% Editorial omission of an illustrative figure.


To pay for the outward reservoir move, we need the squared cost of
changing independent cells to be bounded by the sum of their
individual costs. If $(S_i,S'_i)$ are independent exchangeable
pairs and $S^{(i)}$ denotes $S$ with $S_i$ replaced by $S'_i$, then
\begin{equation}\label{eq:product-exchange}
 \mathbb E\,d(F(S),F(S'))^2
 \le\sum_i\mathbb E\,d(F(S),F(S^{(i)}))^2.
\end{equation}
For two blocks of coordinates form the quadrilateral obtained by
choosing the original or changed value in each block. The CAT$(0)$
quadrilateral inequality bounds the sum of squared diagonals by the
sum of squared sides. Exchanging one block preserves the joint law,
so the diagonals have equal expected squares, and each opposite pair
of sides does too. Dividing by two and iterating over the blocks
proves~\eqref{eq:product-exchange}.

Make independent $D$-moves in all $r$ reservoir cells, keeping the
lamps fixed. Under stationarity the squared cost of this change is at
most $2r\beta\,\mathcal E_G(F)$ by~\eqref{eq:product-exchange} and
the assumed comparison on each cell, since $1-1/r\ge1/2$. The squared
cost of a single random lamp addition is at most $r\mathcal E_G(F)$.
It remains to charge the return move. Integrate the path bound over
\(u\in B\), a uniform lamp array \(s\), and a uniform address \(J\).
The squared endpoint distance is at most three times the sum of the
squared lengths of the three legs. We bound this sum by integrating
over all starts. The first and last legs then have equal average cost,
since $s+e_Jz$ is uniform even conditional on $u,u',J$. The middle leg
has the stationary lamp-addition cost because \(u'\) has law \(\pi^r\).
Thus replacing one random lamp at good starts, with holding elsewhere,
costs at most $Cr\beta\,\mathcal E_G(F)$. Applying
\eqref{eq:product-exchange} to the $m$ independent pairs of old and
new lamp values at each good start gives~\eqref{eq:lamp-compare}.

\smallskip\noindent
\emph{Entropy forced by spreading the lamps.}
Take $T=2\lceil mr\beta\rceil$ and apply Lemma~\ref{lem:power}
to~\eqref{eq:lamp-compare}, using the Ces\`aro average of the lazy
$G$-walk up to time $T$. For placements depending only on the lamps,
the demand energy is $\pi^r(B)\mathcal V(f)$, where $\mathcal V(f)$
is the mean squared distance between images of two independent
uniform elements of $A^m$. The physical displacement of the lamps is
a random translation by some $Z\in A^m$, independent of their uniform
starting value $s$. All reservoir moves and choices of lamps are
independent of $s$. We obtain, with an absolute $K\ge1$,
\begin{equation}\label{eq:lamp-translation}
 \mathcal V(f)\le K\,\mathbb E\, d(f(s),f(s+Z))^2
 \qquad\text{for all placements } f \text{ in the class} .
\end{equation}

Each fixed translation has squared energy at most $4\mathcal V(f)$,
by the triangle inequality through an independent uniform point.
Put $\eta=1/(8K)$. Replace by zero every value of $Z$ whose probability
is less than $2^{-H(Z)/\eta}$. The replaced mass is at most $\eta$
by Markov's inequality for the information of the sampled atom (if
$H(Z)=0$ the support is a singleton). The discarded translations
contributed at most
\(4\eta\mathcal V(f)\) to the expectation in
\eqref{eq:lamp-translation}. Since \(4K\eta=1/2\), the resulting
translation kernel satisfies the same comparison with constant \(2K\).
It is symmetric because $A^m$ has exponent two, and each of its rows
has at most $1+2^{H(Z)/\eta}$ atoms. Corollary~\ref{cor:support}
therefore gives
\begin{equation}\label{eq:lamp-bits}
 mk=\log |A^m|\le C(1+H(Z)).
\end{equation}

\smallskip\noindent
\emph{Coding the physical displacement.}
In a trajectory of the lazy Ces\`aro walk, ignore lazy holds. There
are at most $T$ steps left. Let $S$ be the number of lamp additions,
so $\mathbb E S\le T/r$. Record $S$ in unary and the $S$ lamp
addresses in $\lceil\log m\rceil$ bits each. Within each of the $S+1$
intervals before, between and after additions, record in unary the
number of $P$-steps taken by each reservoir cell. The lengths of
these $r(S+1)$ unary words sum to at most $T+r(S+1)$ bits.
This prefix code for the schedule data has expected length at most
\[
 C\big(T+r+(T/r)\log(2m)\big).
\]
The order of moves in different cells inside an interval is
immaterial, since no lamp is read or written there.

Given the recorded data, the successive labels in a cell have the law
of a stationary $P$-block of the indicated length $j$, including its
initial label. This follows by constructing independent stationary
$P$-trajectories for the cells before sampling the schedule. Since
$H_j/j\to a$, there is $M_0<\infty$, depending only on $P,\pi,L$,
such that $H_j\le (1+a)j+M_0$ for every $j$. The $r$ lengths sum to
at most $T+r$. The schedule and these label blocks determine the
reservoir labels at every lamp addition and hence determine $Z$.
The chain rule for entropy now gives
\begin{equation}\label{eq:lamp-code}
 H(Z)\le C\big((T+r)(1+a)+rM_0+(T/r)\log(2m)\big).
\end{equation}
Combine~\eqref{eq:lamp-bits} and~\eqref{eq:lamp-code} and divide by
$mr$. Let $r\to\infty$ with $m$ fixed. The lamp-address contribution
tends to zero, and $k/r\to h/4$. Next let $m\to\infty$ to remove
$M_0/m$. As $T/(mr)\to 2\beta$ in the first limit, this proves the
theorem.
\end{proof}

\subsection{Conditional entropy of weighted edge laws}\label{sec:cursor-entropy}

For a graph with weighted edges, the demanded pairs and the
physical edges can have different vertex marginals. Symmetrize each
pair law and add suitable holding steps. The resulting walks have a
common stationary law and both squared energies are halved. Applying Theorem~\ref{thm:entropy} with the
identity labeling then bounds the one-step demand entropy.

\begin{corollary}[Pair laws on a graph]\label{cor:graph-entropy}
Let $\mu$ and $\alpha$ be laws on ordered pairs from a finite set, and
suppose that $\mathcal E_\mu(f)\le\beta\mathcal E_\alpha(f)$
for every placement \(f\) in every member of a class of Hadamard
spaces without ports, with \(\beta\ge0\). Put
$\alpha_s=(\alpha+\alpha^\top)/2$. Then
\[
 H_\mu(Y\mid X)\le C(1+\beta)\big(1+H_{\alpha_s}(Y\mid X)\big).
\]
In particular, if the undirected support of $\alpha$ has at most
$d$ neighbors at each vertex, $d\ge2$, then
$H_\mu(Y\mid X)\le C'(1+\beta)\log d$.
\end{corollary}
\begin{proof}
Put $\mu_s=(\mu+\mu^\top)/2$, and write $\rho,\sigma$ for the
marginals of $\mu_s,\alpha_s$, respectively. Write $\Delta_\rho$ for the
law of $(X,X)$ with $X\sim\rho$. The stationary pair laws
$(\mu_s+\Delta_\sigma)/2$ and $(\alpha_s+\Delta_\rho)/2$ define
kernels $D,P$ with common stationary law $(\rho+\sigma)/2$ on its
support. Their energies are half the original energies. Concavity of
conditional entropy gives $H_D(Y\mid X)\ge H_\mu(Y\mid X)/4$.
For the physical step, record a bit indicating which component of the
mixture was used. Conditional on the diagonal component the destination
is known, so $H_P(Y\mid X)\le 1+H_{\alpha_s}(Y\mid X)/2$.
This conditional entropy is the entropy rate for the identity labeling
of the stationary $P$-chain. The first assertion follows by applying
Theorem~\ref{thm:entropy} with $\max(1,\beta)$. In the bounded-degree
case at most $d$ neighbors and a possible holding move give
$H_{\alpha_s}(Y\mid X)\le\log(d+1)\le2\log d$.
\end{proof}
~~~

## P / product.tex

OpenAI, preprints/Circuit-lower-bounds-for-Hamiltonian-cycle-and-the-permanent-November-24-2026/build/product.tex, Git blob a23099e8f89e0e9e2754c9776abe6dc96578fcc0. Complete theorem and proofs.

~~~latex
\section{Products inside a class of Hadamard spaces}\label{sec:product}

Two placements of the same labels can be combined in an $\ell_2$ product.
For $0<\lambda<1$, scale the factors by $\sqrt{\lambda}$ and $\sqrt{1-\lambda}$. Then
the squared distances of the combined placement are the weighted averages
of the two squared distance tables. Section~\ref{sec:lamps} uses this
operation twice, first to place the lamps under individual edge constraints,
then to obtain one placement that gives a saving for every balanced cut.
Both applications require the product to belong to the class of test
spaces for the characteristic.

Recall that for a regular multigraph $G$ and a placement of its vertices,
$V$ is the mean squared distance of two independent uniform vertices and
$E_G$ is the mean squared distance along a uniform directed arc. A class
has \emph{ports} if some graphs of fixed degree and unbounded order
satisfy $V\le C E_G$ for every placement in the class, with a common $C$.

\begin{theorem}[Product theorem]\label{prod:main}
Let $\mathcal C$ be a nonempty isometry-invariant class of Hadamard spaces,
closed under positive rescaling, closed convex subspaces and pointed
ultralimits. If $\mathcal C$ has no ports, it contains a nonempty subclass
with those same properties and no ports that is also closed under finite
$\ell_2$ products.
\end{theorem}

The proof works with finite configurations. A \emph{table} on $A$ is
the array of distances of an $A$-indexed tuple, whose points may coincide.
A table is \emph{allowed} if some space of $\mathcal C$ realizes it.
At each arity the allowed tables are closed under coordinatewise limits,
by closure under pointed ultralimits. Conversely, if all finite tables
of a Hadamard space $X$ are allowed, then $X\in\mathcal C$. Fix a
basepoint in $X$, realize each finite subset containing it, and take an
ultralimit over those subsets along an ultrafilter containing every tail.
This gives an isometric copy of $X$. Its image is closed by completeness
and convex by uniqueness of geodesics. Throughout we allow pointed
ultralimits over arbitrary index sets.

Call a table $d$ on $A$ \emph{power-realizable} if for every $s\ge1$
the table on $A^s$ with squared distances $\sum_{i=1}^s d(a_i,b_i)^2$
is allowed. We first show that no fixed-degree family of unbounded order can expand
on all power-realizable tables. The entropy obstruction supplies grids of words
with large pair variance per coordinate. Averaging pairs of slices makes
these grids flat in the limit. In a flat grid several coordinates can
still move in a common Euclidean direction. Lemma~\ref{prod:extract}
separates that motion from a power-realizable table. If the Euclidean part satisfies a proposed bound on $V/E_G$, a
profile violating that bound must violate it on the power-realizable
part. We extend these tables to entire Hadamard spaces
and use minimality to pass from powers of one table to products of two
different tables.

\subsection{Flat grids and extraction of powers}

A change in one coordinate and a change in a disjoint block trace the
four sides of a quadrilateral. For four points in cyclic order $(a,b,e,c)$,
write $\Delta$ for the sum of the four squared sides minus the sum of
the two squared diagonals. The squared midpoint inequality gives
$\Delta\ge0$ and
\begin{equation}\label{prod:drop}
 d\bigl(\operatorname{mid}(a,c),\operatorname{mid}(b,e)\bigr)^2
 \le \tfrac12\bigl(d(a,b)^2+d(c,e)^2\bigr)-\tfrac14\Delta .
\end{equation}
Thus taking pointwise midpoints of two parallel slices lowers squared
lengths on average, strictly so wherever $\Delta>0$.

Call a quadrilateral \emph{flat} if $\Delta=0$. Its sides and the points
on them have the distances of a Euclidean parallelogram, possibly
degenerate. Indeed, four times the squared distance between the diagonal
midpoints is at most $\Delta$. When $\Delta=0$ the diagonals have a common
midpoint, and equality holds in each comparison of a vertex with the
midpoint of the opposite diagonal. The Euclidean median identities give
equal opposite sides and hence a parallelogram distance table.
Equality at the midpoint implies equality everywhere along that diagonal
by comparison on subintervals. First variation at an endpoint then
gives equality of the angle with its comparison angle, so the corner
triangle is flat. The four corner angles sum to $2\pi$ and the flat
quadrilateral theorem applies. The collinear case follows from equality
in the triangle inequality
\cite[Chapter II.2, Proposition 2.9, Exercise 2.10 and Theorem 2.11]{BH}.

We also use two averaging inequalities. The Cartan barycenter minimizes
mean squared distance and the distance between two barycenters is at most
the mean distance under any coupling of the two laws, by Jensen's
inequality on the product space \cite{Sturm}. For independent exchangeable
pairs $(X_i,X'_i)$, write $X^{(i)}$ for $X$ with coordinate $i$ replaced
by $X'_i$. The coordinate-exchange inequality, Equation~\eqref{eq:product-exchange},
bounds the cost of changing all coordinates by
\begin{equation}\label{prod:exchange}
 \mathbb E d(f(X),f(X'))^2
 \le\sum_i\mathbb E d(f(X),f(X^{(i)}))^2
\end{equation}
for every Hadamard-valued placement $f$.

A grid of dimension $s$ is a placement of $A^s$ in a space of
$\mathcal C$. To average slices while keeping track of limits, include
in its distance table all formal iterated midpoints of the grid vertices.
There are countably many such labels. Their distances are bounded by
the diameter of the grid, by convexity. At fixed dimension the enlarged
tables of diameter at most $R$ form a compact set in the coordinatewise
topology, since any limit is realized in a pointed ultralimit.

A \emph{system} assigns a law of these tables to each dimension. Its
laws are invariant under permutations of coordinates. They are also
consistent in the following sense. Fixing a coordinate to an independent
uniform letter in an $(s+1)$-grid gives the law of the $s$-grid.
Uniform tightness of diameters at each fixed dimension therefore permits
subsequential limits of systems by diagonal extraction.

Let $V_s$ be the expected uniform pair variance of the $s$-grid. Given
a coupling of the uniform law on $A^j$ with itself, draw $(u,v)$ from
that coupling independently of a $(j+s)$-grid. Take pointwise midpoints
of its two $s$-dimensional slices at $u,v$. Doing this in each dimension
defines a new system. By consistency both slices have the old law,
so Equation~\eqref{prod:drop} shows that every $V_s$ weakly decreases.
If a finite $V_s$ stays unchanged, its nonnegative mean deficit vanishes.
Almost surely the square for every supported $(u,v)$ and every transition
in the other $s$ coordinates is then flat. A system is \emph{flat} when
almost surely this holds for every pair of disjoint blocks and every
choice of transitions in them.

\begin{lemma}\label{prod:extract}
If the one-coordinate distance profile of a flat system belongs with positive
probability to a compact set $\mathcal K$, some $d\in\mathcal K$ admits
a decomposition
\[
 d(a,b)^2=\rho(a,b)^2+\|v_b-v_a\|^2 ,
\]
where $\rho$ is power-realizable and the $v_a$ belong to a Euclidean space.
\end{lemma}

Homogeneous coordinates in a large flat grid have two kinds of motion.
Aligned increments in the common direction contribute a term quadratic
in their number to the squared distance. The remaining contributions
add coordinatewise. Auxiliary coordinates will cancel the common
direction and realize powers of the remaining table.

\begin{proof}
On a flat grid a distance depends only on the coordinates that change and
their transitions, since opposite sides of a flat square are equal. For a
transition $u=(a,b)$ in coordinate $i$ write $d_i(u)$ for its length.
For two transitions in distinct coordinates $i<j$ define $B_{ij}(u,v)$
by subtracting $d_i(u)^2+d_j(v)^2$ from the squared distance of the joint
change and dividing by two. Reversing either transition changes its sign.
For any list of transitions the squared distance is
\begin{equation}\label{prod:gram}
 \sum_i d_i(u_i)^2+2\sum_{i<j} B_{ij}(u_i,u_j).
\end{equation}
To see this, partition the active coordinates into two nonempty blocks.
The parallelogram identity says that the squared length of the joint
change plus its squared length with one block reversed equals twice
the sum of the two block lengths squared. The right-hand side of
Equation~\eqref{prod:gram} satisfies the same identity. By induction on
the number of active coordinates, the error consequently changes sign
on reversing any proper nonempty block. For three or more coordinates,
reversing two of them together changes its sign, while reversing them
one at a time preserves it. The error is zero. For one or two coordinates
the formula holds by definition.

Let $p_0>0$ be the probability that a single-coordinate profile lies
in $\mathcal K$. An $m$-grid has $mp_0$ such coordinates in expectation,
so there exist flat grids with arbitrarily many of them. Their
one- and two-coordinate coefficients are bounded by compactness of
$\mathcal K$ and the triangle inequality. Color each coordinate by the
mesh cell of its profile and each pair $i<j$ by that of its ordered
coefficient array. Finite Ramsey theory gives, for every $m$, $m$
coordinates on which each of these colors is constant at mesh size
$1/m$. Take a convergent subsequence of the representative coefficient
arrays. The limits are a single $d\in\mathcal K$ and a single kernel $B$.
For each fixed dimension the restricted grids converge to the table
in Equation~\eqref{prod:gram} with $d_i=d$ and $B_{ij}=B$ for
all $i<j$. Closedness realizes that grid.

The kernel $B$ on ordered transitions is symmetric. In
Equation~\eqref{prod:gram} use $u,v,-u,-v$ in four successive blocks of equal
length $l$, where minus denotes reversal. The coefficient of $l^2$ is
$2(B(u,v)-B(v,u))$ and must be nonnegative. Exchanging $u,v$ gives the
reverse inequality. For any integer combination of transitions, use $l$
times the absolute multiplicities and reverse the transitions with negative
multiplicity. The $l^2$ coefficient is the corresponding quadratic form
of $B$. Homogeneity and density now imply positive semidefiniteness over
the reals. Represent $B$ by Gram vectors $v_{ab}$, with $v_{ba}=-v_{ab}$.
The change $(a^l,b^l,c^l)\to(b^l,c^l,a^l)$ is a product of two
block transpositions. Each transposition has zero $l^2$ coefficient
because its two transitions are reverses. Its length is $O(\sqrt l)$.
The triangle inequality therefore forces
$v_{ab}+v_{bc}+v_{ca}=0$ in the direct change. Hence $v_{ab}=v_b-v_a$.
Two opposite transitions on distinct coordinates show that
$\rho(a,b)^2=d(a,b)^2-\|v_b-v_a\|^2$ is nonnegative. Equation~\eqref{prod:gram} is now the sum of the residual squares in each coordinate
plus the squared norm of the sum of its vector increments.

It remains to cancel the shared vector increment. In one auxiliary
coordinate replace the transition $a\to b$ by its geodesic interval
parametrized by $[0,1]$. Suppose the background coordinates change with
residual squared cost $R$ and vector increment $W$. The two copies of
the interval bound a flat parallelogram. Its diagonal formula gives
the scalar product of the background side with the interval side as
$\langle W,v_b-v_a\rangle$. Interpolating on the two sides at parameters
with signed difference $t$ therefore gives squared distance
\[
 R+t^2\rho(a,b)^2+\|W+t(v_b-v_a)\|^2 .
\]
This formula satisfies the parallelogram identity on disjoint changes,
so we may interpolate further auxiliary coordinates in turn.

Set $v_o=0$ at one reference letter. For each original coordinate and
each letter $a\ne o$, adjoin $M$ intervals directed from $a$ to $o$.
When that coordinate reads $a$, use parameter $1/M$ on its $M$
intervals for $a$ and zero on the others. Relative to the zero-parameter
settings, the auxiliary increment is the negative of the original
word's vector. Thus the shared term vanishes between any two words. At fixed $s$ the total
residual cost of the auxiliary changes is $O(1/M)$, since $O(M)$
intervals change by $1/M$. Closedness as $M\to\infty$ realizes the
table with squared distances $\sum_{i=1}^s\rho(a_i,b_i)^2$ on $A^s$.
In particular $\rho$ is a pseudometric.

\end{proof}

\subsection{The graph obstruction on power-realizable tables}

\begin{lemma}\label{prod:power-obstruction}
If $\mathcal C$ has no ports, no fixed-degree graph family of unbounded
order has a common Poincar\'e bound on its power-realizable tables.
\end{lemma}
\begin{proof}
Suppose a family had such a bound. At the same orders superpose
bounded-degree scalar expanders, available as regular multigraphs at all
sufficiently large orders \cite{HLW}. Their scalar inequality extends
to Euclidean placements by summing over coordinates. The new graph
energies dominate those of either constituent up to the ratio of the
degrees. Thus there is a common bound $C_0$ on both Euclidean and
power-realizable tables. Fix one resulting graph $G$ of degree $D$ on
$A$, with $n=|A|$ large, and put
\[
 h=c\frac{\log n}{\log(2D)}
\]
for a sufficiently small absolute $c>0$. We will construct a flat
system with $V_1\ge h$ and expected one-coordinate $G$-energy at
most one. To retain this lower bound after averaging away defects, we
first need grids with pair variance proportional to their dimension.
In the entropy calculation logarithms are to base two.

On words in $A^p$, update a uniformly chosen coordinate by one $G$-step.
Identify words that differ by a cyclic rotation. Lifting a placement
of these orbits back to words gives every coordinate the same energy,
while passing to orbits loses at most $\log p$ bits of word entropy.
Let $\pi$ be the orbit law induced by the uniform word, and let $\beta$
be the supremum, over placements in $\mathcal C$, of the ratio of
independent $\pi$-pair energy to stationary one-step energy, enlarged
to one if necessary. The quotient is connected, so $\beta$ is finite.

For $t=2\lceil\beta\rceil$, Lemma~\ref{lem:power} bounds the same
comparison by an absolute constant against the lazy Ces\`aro walk
of at most $t$ steps. Apply Theorem~\ref{thm:entropy} to this smoothed
walk, independent demand and the identity labeling of the orbits.
The demand has conditional entropy $H(\pi)\ge p\log n-\log p$.
For the walk, the entropy rate is the one-step conditional entropy.
Given a starting orbit, fix one representative word. To describe any
endpoint, record the number of updates in each coordinate and the
local arcs in their order within that coordinate. Updates in different
coordinates commute, and holds may be omitted. There are at most
$\binom{p+t}{p}$ count vectors and at most $D^t$ choices of arcs
for each vector. Theorem~\ref{thm:entropy} therefore gives
\[
 p\log n-\log p\ \le\
 C\left(1+\log\binom{p+t}{p}+t\log D\right).
\]
Using $\binom{p+t}{p}\le2^{p+t}$ and $t=2\lceil\beta\rceil$ now gives
$\beta\gtrsim p\log n/\log(2D)$, first taking $n$ sufficiently large
and then $p$ large. Lift a placement attaining at least half this ratio
to a rotation-invariant placement of $A^p$. Rescale it so that the
mean $G$-arc energy in each coordinate, averaged over the background
word, is one. For the choice of $c$ above, its uniform pair variance
is at least $4ph$.

Passage to a system requires uniform integrability of squared distances
on fixed-dimensional slices. Smooth a placement $f$ by independently
resampling each letter with probability $\epsilon$ and assigning to
each word $x$ the barycenter of the resulting images. At a fixed
one-coordinate transition the old mean squared distance over background
words is $O_G(1)$, by following a fixed path between the two letters and
using the arc-energy normalization. Equation~\eqref{prod:exchange}
then bounds the expected squared distance between the old image of a
uniform word and its resampled image by $O_G(\epsilon p)$. Jensen's
inequality gives the same bound for the displacement to its barycenter.
By the triangle inequality in $L^2$, choosing $\epsilon>0$ sufficiently
small at fixed $G$ leaves pair variance at least $ph$.

To compare two smoothed words differing in one coordinate, resample
their other letters identically. In the changed coordinate either keep
both old letters or replace both by the same fresh letter. Barycenter
contraction bounds the new transition distance by the scalar noise
operator on the background coordinates applied to the old distance.
Consequently the expected $G$-arc energy in each coordinate is at most one.

On the fixed alphabet the noise operator
$T_\epsilon f=(1-\epsilon)f+\epsilon\mathbb E f$ satisfies
$\|T_\epsilon f\|_{2+\eta}\le\|f\|_2$ for some $\eta>0$
(see \cite{odon} for background). By positivity it suffices to verify
this for nonnegative functions. The $L^2$ contraction is strict away from constants,
uniformly on compact subsets of the unit sphere disjoint from constants.
Near the constant $1$, factor out the mean and expand the squared
output norm at $1+g$ with $\mathbb E g=0$. The quadratic coefficient
of $\mathbb E g^2$ is $(1+\eta)(1-\epsilon)^2<1$ for small $\eta$,
and the remainder is uniform on the finite alphabet. This proves
the contraction, which tensorizes by Minkowski's inequality. Every
smoothed one-coordinate transition thus has bounded $(2+\eta)$-moment
uniformly in the word length $p$.

Randomly permute the coordinates and fix the unused letters uniformly.
At every fixed dimension the diameter has bounded $(2+\eta)$-moment,
by connecting grid vertices along the finitely many one-coordinate edges.
Take $p\to\infty$ along factorials and extract a limit at every dimension.
Before passage to the limit, partition the coordinates into blocks of
size $k$, for any $k$ dividing $p$. Equation~\eqref{prod:exchange}
gives $V_p\le(p/k)V_k$. Uniform integrability on the restricted grids
gives $V_k\ge kh$ in the limit. The system also retains the upper bound
of one on the expected single-coordinate $G$-energy and uniform
$(2+\eta)$-moment bounds for the single-coordinate transition lengths.

Among systems obeying these three bounds, minimize
$\sum_{k\ge1}2^{-k}V_k$. The diameter moment bounds give tightness
of the grid laws in each dimension, and consistency passes to limits.
Uniform integrability preserves second moments, while lower
semicontinuity preserves the upper $(2+\eta)$-moment bounds. Since
$V_k\le kV_1=O_G(k)$ by Equation~\eqref{prod:exchange}, the objective
is continuous on this compact set.

In the minimizing system take pointwise midpoints of slices at two
independent uniform settings on $j$ extra coordinates. Convexity and
Jensen's inequality preserve the two upper bounds. For the lower bound,
let $V'_l$ be the new variance in dimension $l$. Couple the new grid
with its first slice. At a uniform word the point moves by half the
length of the segment between the slices, so its mean squared
displacement is $V_j/4$. Comparing two uniform words gives
$\sqrt{V'_l}\ge\sqrt{V_l}-\sqrt{V_j}$ by the triangle inequality
in $L^2$. For each $k$, partition the $l$ coordinates into blocks of
size $k$ to get $V'_k\ge(k/l)V'_l$. Let $l\to\infty$ through multiples
of $k$. Since $V_l\ge lh$, this gives $V'_k\ge kh$.

The smoothed system is therefore a competitor in the minimization.
Every variance weakly decreases, so minimality forces equality in
every dimension. The two slice settings have full support. By
Equation~\eqref{prod:drop}, every square with one change in those $j$
coordinates and one change in the others is flat almost surely.
Varying $j$ gives a flat system.

Take $n$ sufficiently large that $h>C_0$. Since $V_1\ge h$ and
the expected $G$-energy is at most one, the one-coordinate table $d$
satisfies $V(d)>C_0E_G(d)$ with positive probability. Some compact
subset of that event still has positive probability.
Lemma~\ref{prod:extract} decomposes one table in this subset into
squared distances of a power-realizable table and a Euclidean table.
Both satisfy the comparison with constant $C_0$. Adding their
inequalities gives the same comparison for $d$, a contradiction.
\end{proof}

\subsection{Midpoint extensions and minimality}

To find a subclass of spaces using only power-realizable tables, we
must enlarge each such table to a complete convex space. Restrictions,
repeated labels, positive rescalings and coordinatewise limits preserve
power-realizability. The issue is adding midpoints. For the midpoint of
$a,b$ we know the two half lengths in advance. Its distances to the
other labels must be chosen consistently across every power.

\begin{lemma}\label{prod:midextend}
Any power-realizable table admits a power-realizable extension by a midpoint
of any specified pair of its labels.
\end{lemma}
\begin{proof}
Adjoin a letter $z$ to the old alphabet. Prescribe all old distances
and $d(a,z)=d(z,b)=d(a,b)/2$. Call two letters \emph{linked} if they
are equal or their distance is prescribed. Two words are linked if
their letters are linked coordinatewise. Prescribe the squared distance of each linked pair of words as the sum
of its prescribed coordinate squares.

These requirements can be realized in every dimension. Start with a
power of the old table and insert $z$ in one coordinate at a time as
the midpoint of the $ab$ segment for each setting of the other
coordinates. If two such settings are linked, the six prescribed
distances at the ends of their $ab$ segments give a Euclidean rectangle,
possibly degenerate. Interpolating on these segments therefore gives
the required distances between two new midpoints and between a midpoint
and either linked endpoint.

Permute coordinates uniformly and pass to a system limit. The diameter
of each grid is bounded by a constant times its dimension, by connecting
words through linked changes. Among systems realizing the partial table,
minimize $\sum_{k\ge1}2^{-k}V_k$. The diameter bound gives compactness
and continuity as in the proof of Lemma~\ref{prod:power-obstruction}.

We apply midpoint averaging twice to this minimizer. First, on $j$
extra coordinates choose a coupling of uniform settings supported on
linked pairs, with positive probability on every such pair. To construct
it, put a small equal mass on each off-diagonal linked pair and fill
the row marginals on the diagonal. Symmetry gives the column marginals.
For a linked change in the remaining coordinates, the four slice
vertices form a rectangle by the prescribed distances. Averaging
along its opposite sides preserves their separation. Thus the smoothed
system still realizes the partial table. By minimality all variances are unchanged. Equation~\eqref{prod:drop}
therefore forces flatness in the original system for every square with
one linked change and \emph{any} disjoint change.

Second, average slices of the same original minimizer at independent
uniform settings. For each linked change in the remaining coordinates,
the preceding flatness gives a parallelogram, so the distance between
midpoints is again the required length. This smoothing also preserves
the partial table. Equality of all variances, now with full support
on settings, gives a flat system.

Extract homogeneous flat grids as in Lemma~\ref{prod:extract}. Linked
distances survive restriction and limits. Repeating a linked transition
$u=(x,y)$ in two coordinates gives squared length $2d(x,y)^2$.
Equation~\eqref{prod:gram} thus gives $B(u,u)=\|v_y-v_x\|^2=0$.
Linked pairs connect the alphabet, so all the shared vectors are equal.
The one-coordinate table itself is power-realizable. In any realization,
the two prescribed half lengths place $z$ at the midpoint of $a,b$.
\end{proof}

Iterate Lemma~\ref{prod:midextend}, eventually adjoining a midpoint
for every pair of labels so obtained. Realize all these nested tables
together in a pointed ultralimit. The closure of their union is
midpoint-closed by continuity and hence convex by dyadic approximation.
All its finite tables are power-realizable by approximation. The
subclass of $\mathcal C$ consisting of spaces whose every finite table
is power-realizable therefore realizes \emph{every} power-realizable
table. It is closed under isometry, positive rescaling, closed convex
subspaces and pointed ultralimits, by the corresponding properties of
the tables. Lemma~\ref{prod:power-obstruction} now shows that this
subclass has no ports.

\begin{proof}[Proof of Theorem~\ref{prod:main}]
Choose a minimal subclass $\mathcal M\subseteq\mathcal C$ with the
stated closure properties and without ports. We justify this choice
using Zorn's lemma on the systems of allowed finite tables, which
determine their classes. Consider a chain decreasing by inclusion.
For a table common to the chain, realize the tuple at every stage
and take a pointed ultralimit along an ultrafilter containing every
tail towards smaller classes, with basepoint in the tuple. The table
is realized in the limit, whose finite tables all belong to every
stage by closedness. Thus the intersection class realizes every table
in the intersection of the systems.

Every class in the chain realizes the unit two-point table, since it
contains a nontrivial space and is closed under rescaling. This remains
true in the intersection. For a fixed connected $D$-regular graph on
$n$ vertices, Theorem~\ref{thm:loggraph} bounds the Poincar\'e constant of each class below
by $c_0\log n/\log(1+D)$. The tables of edge energy one are bounded
by connectedness. Compactness gives a
witness for this lower bound in each class, and nested compactness
gives one in the intersection. A disconnected graph fails any finite
comparison on the unit two-point table, by placing its components
at different points. Thus the intersection has no ports.

Apply the preceding construction of the subclass with power-realizable
tables inside $\mathcal M$. It has no ports, so minimality makes it
all of $\mathcal M$. This deals with powers of a single table. To
combine two different tables, we first put their shapes approximately
in the same space. For each nonzero table choose an open neighborhood
of its diameter-one normalization among tables of its arity. The
spaces in $\mathcal M$ avoiding that neighborhood among their normalized
tables form a proper subclass with the same closure properties. For
ultralimits this follows because a table of positive diameter with
normalization in the open set is approximated by tables whose
normalizations also lie there. By minimality each avoiding subclass
has ports.

We need a common port family for the two avoiding subclasses. In
any class with ports, take a graph of degree $D$ and order $q$ in its
port family with the common comparison bound $b\ge1$. Choose $q$ so
large that an even $t\ge2Bb$ satisfies $(2D+1)^t\le q/2$, with $B$
as in the proof of Theorem~\ref{thm:loggraph}. The iteration in that
proof gives graphs of degree $2D$ with bounded comparison on every
order $q^k$. To obtain every intermediate order $m$, take the least
such power $N\ge m$ and distribute its vertices onto $m$ fibers of
sizes $\lfloor N/m\rfloor$ and $\lceil N/m\rceil$. The pushforward
of the uniform measure dominates half the uniform measure on $m$
vertices. Retain the arc multiplicities between the fibers and pad
with loops to degree $2qD$, which suffices since $N<qm$ for $m>1$.
For any placement on the quotient the new arc energy is $N/(qm)$
times the arc energy of its pullback, and the new pair variance is
at most four times that of the pullback. The comparison bound
increases by a factor of at most $4q$. We may therefore take port
graphs at the same orders for both subclasses. Superposing them gives
a common port family, since each constituent arc energy is bounded
by the new arc energy times the ratio of degrees. An empty subclass
imposes no constraints.

If the two avoiding subclasses covered $\mathcal M$, their common
family would give ports for $\mathcal M$. Hence some space in
$\mathcal M$ realizes tuples approximating both prescribed shapes,
possibly at different scales. Take their union as one finite table
and use its power-realizability. In a large power vary copies of the
first tuple in one block of coordinates and copies of the second in
a disjoint block. This realizes the product of the two tuples, with
their squared scales multiplied by the respective block lengths.
The ratio of those lengths can approximate any positive number.
A common rescaling and a limit therefore give any desired ratio of
diameters. Shrinking both neighborhoods gives the exact product table
by closedness. A zero-diameter factor only repeats labels.

Every finite table in a product of two spaces of $\mathcal M$ is the
restriction of a product of two allowed tables, obtained by projecting
its labels to the two factors. The product Hadamard space thus belongs
to $\mathcal M$, and induction gives all finite $\ell_2$ products.
\end{proof}
~~~
