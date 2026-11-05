# The payload disperser argument

This file proves, from elementary finite-field algebra, the one
algebraic fact that `proof.md` consumes: the payload function is an
affine disperser with an explicit dimension threshold. Nothing in this
file is cited as a black box. The argument is the subspace-polynomial
proof of Ben-Sasson–Kopparty (BK), *Affine Dispersers from Subspace
Polynomials*, specialised to $\mathbf F_2$ and prime field degree:
the preliminaries and structural lemma are BK Sections 3.2–3.4
(Lemmas 3.6, 3.7, 3.10, 3.13), the block decomposition and selection
are BK Section 7 (Lemma 7.2, Theorem 7.3), and the final isolation is
BK's coefficient argument with the weighted-capacity count of
Section 7 below. Prime degree is used in exactly one place: the
subfield exclusion in the structural lemma, where it makes the
hypothesis automatic.

**Lemma D.** *Let $k$ be prime, let $K=\mathbf F_{2^k}$, and let
$r,t$ be integers with $1\le t<r$ and $kt\ge r(r+1)$. Let
$\pi:K\to\mathbf F_2$ be a nonzero $\mathbf F_2$-linear map. For
$i=1,\dots,r$ put $w_i=r+2-i$ and define*

$$F(a_1,\dots,a_r)=\pi\Bigl(\mathrm{Sym}^t\bigl(a_1^{2^{w_1}-1},
\dots,a_r^{2^{w_r}-1}\bigr)\Bigr).$$

*Then $F$ is nonconstant on every affine subspace $A\subseteq K^r$
(over $\mathbf F_2$) with*

$$\dim A \;>\; \frac{kr}{t}+kt+r(r+1).$$

## 1. Preliminaries

Write $q=2^k$. A function $K^r\to K$ has a unique representative
polynomial with each variable of degree $<q$: reduce modulo
$X_i^q-X_i$; uniqueness follows by induction on $r$ from the fact that
a nonzero univariate polynomial of degree $<q$ cannot vanish on all
$q$ points of $K$. In particular, a nonzero coefficient of a
positive-degree monomial in the reduced representative proves the
function is nonconstant.

A *linearized* polynomial is $\sum_{j=0}^{k-1} c_jX^{2^j}$. These are
exactly the $\mathbf F_2$-linear maps $K\to K$ given by polynomials:
each is linear by the Frobenius identity
$(x+y)^{2^j}=x^{2^j}+y^{2^j}$; the coefficient space and the space of
linear maps both have $\mathbf F_2$-dimension $k^2$; and the
coefficient-to-function map is injective by the root bound (a nonzero
linearized polynomial has degree $\le 2^{k-1}<q$), hence surjective.
Consequently every $\mathbf F_2$-linear map $K\to K$ has a unique
reduced linearized representative, a fact used for the mixing maps in
Section 5.

The trace $\mathrm{Tr}(z)=\sum_{h=0}^{k-1} z^{2^h}$ maps $K$ to
$\mathbf F_2$ and is not the zero function (its polynomial has degree
$2^{k-1}<q$). For $\mu\ne0$, $z\mapsto\mathrm{Tr}(\mu z)$ is a nonzero
$\mathbf F_2$-linear functional, and $\mu\mapsto\mathrm{Tr}(\mu\,\cdot)$
is an injective linear map from the $k$-dimensional space $K$ to the
$k$-dimensional dual of $K$, hence a bijection: every nonzero linear
$\pi:K\to\mathbf F_2$ is $\pi(z)=\mathrm{Tr}(\mu z)$ for a unique
$\mu\ne0$.

For an integer $b$, $\mathrm{wt}(b)$ is the number of ones in its
binary expansion. An exponent arising from a product of $h$ linearized
factors is, before reduction, a sum of $h$ powers of two; reduction
modulo $X^q-X$ replaces $2^k$ by $1$ (wrap) and merges equal powers
$2^s+2^s=2^{s+1}$ (carry). Wraps preserve the number of power-of-two
summands and carries strictly decrease it. Hence a product of $h$
linearized factors has, after reduction, all monomials of binary
weight at most $h$, and a monomial of binary weight exactly $h$ arises
only from summand residues that are pairwise distinct modulo $k$.
This is the only degree bookkeeping used below.

## 2. Kernel and image polynomials

For a linear subspace $V\subseteq K$ of dimension $d$, let
$P_V(X)=\prod_{v\in V}(X-v)$, the *kernel polynomial* of $V$. Its root
set in $K$ is exactly $V$.

**$P_V$ is a monic linearized polynomial of degree $2^d$.** Induction
on $d$: $P_{\{0\}}=X$; if $V=U\oplus\langle b\rangle$ with
$\dim U=d-1$, then once $P_U$ is known linearized,
$P_U(X+b)=P_U(X)+P_U(b)$, so
$P_V(X)=\prod_{u\in U}(X-u)(X-u-b)
=P_U(X)\bigl(P_U(X)+P_U(b)\bigr)=P_U(X)^2+P_U(b)P_U(X)$,
which is monic linearized of degree $2^d$ with root set
$U\cup(U+b)=V$.

Let $V'=P_V(K)$, the image of $P_V$; it is a subspace of dimension
$k-d$. Define the *image polynomial* of $V$ by $Q_V=P_{V'}$.

**Properties.** (i) $Q_V(P_V(X))=X^q-X$. (ii) $P_V(Q_V(X))=X^q-X$.
(iii) $Q_V(K)=V$ and $\ker Q_V=V'$, so $Q_V$ has exactly $2^{k-d}$
roots in $K$ and is monic linearized of degree $2^{k-d}$.

*Proof.* (i): $Q_V\circ P_V$ is monic of degree $2^{k-d}2^d=q$ and
vanishes on $K$ (for $x\in K$, $P_V(x)\in V'$ and $Q_V$ vanishes on
$V'$), so it equals $X^q-X$, the unique monic degree-$q$ polynomial
vanishing on $K$. (ii): let $g(Y)=P_V(Q_V(Y))-(Y^q-Y)$. Substituting
$Y=P_V(X)$ and using (i) gives
$g(P_V(X))=P_V(X^q-X)-\bigl(P_V(X)^q-P_V(X)\bigr)$, and the two terms
are equal: $P_V(X^q-X)=P_V(X^q)-P_V(X)$ by linearity of $P_V$, and
$P_V(X^q)=\sum_jc_jX^{2^{j+k}}=\sum_jc_j^{q}X^{2^{j+k}}=P_V(X)^q$
because the coefficients $c_j$ lie in $K=\mathbf F_{2^k}$, so
$c_j^q=c_j$. Hence $g(P_V(X))=0$ as a polynomial; since $P_V$ is
nonconstant, $g=0$ (composing a nonzero polynomial with a nonconstant
polynomial gives a nonzero polynomial), proving (ii). (iii): by (ii), for
$P_V(Q_V(x))=x^q-x=0$, so $Q_V(x)\in\ker P_V=V$, i.e.
$Q_V(K)\subseteq V$; $\ker Q_V$ is the root set of $P_{V'}$, namely
$V'$, of dimension $k-d$, so $\dim Q_V(K)=k-(k-d)=d=\dim V$ and
$Q_V(K)=V$. $\square$

Write $Q_V(X)=\sum_{j=0}^{k-1}a_jX^{2^j}$ with indices read modulo $k$
and $a_j=0$ for $k-d<j<k$; monicity is $a_{k-d}=1$. The indices
$k-d+1,\dots,k-1$ are the *standard zero run*, of length $d-1$.

## 3. The structural lemma (zero runs)

**Rotation.** *For an integer $s$, $Q(X)^{2^s}\bmod (X^q-X)$ has
coefficient $a_{j-s}^{2^s}$ at $X^{2^j}$ (indices mod $k$), and its
value at $x\in K$ is $Q(x)^{2^s}$.*

*Proof.* $Q(X)^{2^s}=\sum_j a_j^{2^s}X^{2^{j+s}}$; reduce
$X^{2^{j+s}}$ modulo $X^q-X$ by wrapping exponents at $k$ and reindex.
Values: both sides are functions $K\to K$ agreeing by the same
computation, since $x^{2^k}=x$. $\square$

**Structural lemma.** *Let $\dim V=d$ and $Q_V$ as above.*
1. *No $d$ cyclically consecutive coefficients of $Q_V$ are all zero.
   In particular $a_0\ne0$ and $a_{k-d}\ne0$.*
2. *If $d\ge2$, no $d-1$ cyclically consecutive coefficients are all
   zero, except the standard run.*

*Proof.* (1) If $d=k$ then $Q_V=X$ and both parts are immediate from
$a_0=1$, so assume $d<k$. Suppose $a_j=0$ for $j$ in a cyclic interval
$J$ of length $d$. Apply the rotation with $s=k-(r+d)$ where $J=\{r,\dots,
r+d-1\}$: by the rotation formula the rotated polynomial $Q'$ has
coefficients $a_{j-s}^{2^s}$; for $j\in[k-d,k-1]$ the index $j-s$ lies
in $J$, so all top $d$ coefficients of $Q'$ vanish and
$\deg Q'\le 2^{k-d-1}$. But $Q'(x)=Q(x)^{2^s}$, so $Q'$ vanishes on
the $2^{k-d}$ roots of $Q$, and $Q'\ne0$ (rotation permutes
coefficients up to Frobenius, and $a_{k-d}=1$ survives). A nonzero
polynomial of degree $\le2^{k-d-1}$ cannot have $2^{k-d}$ roots:
contradiction. The "in particular" follows because indices
$k-d+1,\dots,k-1$ already form a zero run of length $d-1$: a further
zero at $k-d$ or at $0$ would extend it to a run of length $d$.

(2) Suppose $a_j=0$ on a cyclic interval $J=\{r,\dots,r+d-2\}$ of
length $d-1$ that is not the standard run. Set $s=k-(r+d)+1$; since
$J$ is not the standard run, $s\not\equiv0\pmod k$, and since
$|J|=d-1\le k-1$ we may take $1\le s\le k-1$. Rotate by $s$: the
resulting $Q'$ has, for $j\in[k-d+1,k-1]$, coefficients
$a_{j-s}^{2^s}$ with $j-s\in J$, hence zero; so
$Q'=cX^{2^{k-d}}+$ (lower terms) with $c=a_{r-1}^{2^s}$. If $c=0$ then
$\deg Q'\le2^{k-d-1}$ and part (1)'s root count gives $Q'=0$,
impossible. So $c\ne0$; then $Q'-cQ_V$ has degree $\le2^{k-d-1}$
(leading terms cancel) and vanishes on all roots of $Q_V$ (as
$Q'(x)=Q_V(x)^{2^s}$), hence $Q'=cQ_V$ identically. Therefore for
every $x\in K$, $Q_V(x)^{2^s}=c\,Q_V(x)$; since $Q_V(K)=V$, every
$z\in V$ satisfies $z^{2^s}=cz$. Fix $\beta\in V$, $\beta\ne0$
(possible as $d\ge2\ge1$). For any $z\in V$,
$(z/\beta)^{2^s}=cz/(c\beta)=z/\beta$, so $z/\beta$ is fixed by the
$2^s$-Frobenius, i.e. $z/\beta\in\mathbf F_{2^{\gcd(s,k)}}=\mathbf
F_2$ because $k$ is prime and $1\le s\le k-1$. Hence $z\in\{0,\beta\}$
and $d\le1$, contradicting $d\ge2$. $\square$

Note where primality entered: only through
$\mathbf F_{2^{\gcd(s,k)}}=\mathbf F_2$. For composite degree the same
proof gives the conclusion only for subspaces not contained in a
constant multiple of a proper subfield; for prime $k$ every subspace
of dimension $\ge2$ qualifies automatically, since a constant multiple
of a proper subfield has dimension $\le1$.

## 4. A nonzero coefficient for independent blocks

Let $V_1,\dots,V_t\subseteq K$ be subspaces, $d_i=\dim V_i>1$, and let
$Q_i=Q_{V_i}=\sum_j a_{(i,j)}X^{2^j}$. Fix integers
$1\le h_i<d_i$ and put $e_i=h_i-1$, $\delta_i=d_i-h_i>0$. Assume

$$\sum_i\delta_i \;>\; k-\max_i d_i+1. \tag{$\dagger$}$$

**Claim.** *The reduced polynomial
$\mathrm{Tr}\bigl(\mu\prod_i Q_i(X_i)^{2^{h_i}-1}\bigr)$ has a nonzero
coefficient at some monomial $\prod_i X_i^{b_i}$ with
$\mathrm{wt}(b_i)=h_i$ for all $i$.*

*Proof.* Since $2^{h_i}-1=\sum_{j=0}^{e_i}2^j$,
$Q_i(X_i)^{2^{h_i}-1}=\prod_{j=0}^{e_i}Q_i(X_i)^{2^j}$.
Define the $k\times k$ matrix $T_i$ by
$T_i[u,v]=a_{(i,v-u)}^{2^u}$ (indices mod $k$): this is the coefficient
of $X^{2^v}$ in $Q_i(X)^{2^u}$. Expanding the $h_i$ Frobenius factors,
a monomial $X_i^{b_i}$ with $\mathrm{wt}(b_i)=h_i$,
$b_i=\sum_{v\in B_i}2^v$, $|B_i|=h_i$, receives contributions only
from choices of distinct residues (Section 1), i.e. from bijections
between factors and columns; its coefficient in $Q_i^{2^{h_i}-1}$ is
$\mathrm{perm}\,T_i[A_i,B_i]$ with $A_i=\{0,\dots,e_i\}$. Expanding
the trace as $\sum_{\rho=0}^{k-1}\mu^{2^\rho}(\,\cdot\,)^{2^\rho}$ and
shifting indices gives the coefficient of $\prod_iX_i^{b_i}$ in the
traced polynomial as

$$\sum_{\rho=0}^{k-1}\mu^{2^\rho}\prod_i
\bigl(\mathrm{perm}\,T_i[A_i,B_i-\rho]\bigr)^{2^\rho}. \tag{1}$$

We choose the $B_i$ so that the $\rho=0$ summand is a product of
nonzero permanents and every $\rho\ge1$ summand has a zero factor.

*Triangular minors at $\rho=0$.* Order the sources so $d_1$ is
maximum. Let $J$ be the first index with
$\sum_{i\le J}\delta_i\ge k-d_1+1$ (it exists by $(\dagger)$), and for
$i\le J$ put $\ell_i=\sum_{j<i}\delta_j$ (so $\ell_1=0$ and
$\ell_i\le k-d_1\le k-d_i$, with no wraparound at the first crossing:
the crossing sum is at most $k-d_1+\delta_J\le k-d_1+(d_J-1)\le k-1$);
for $i>J$ put $\ell_i=0$. Fix one source, write
$d,e,\ell,a_j$ for its data, and let $B$ be chosen as follows.

- If $a_\ell\ne0$: $B=(\ell,\ k-d+1,\ \dots,\ k-d+e)$.
- If $a_\ell=0$: $\ell$ lies in a zero gap $a_{l+1}=\dots=a_{j-1}=0$
  with $a_l,a_j\ne0$, $l<\ell<j\le k-d$; the gap is not the standard
  run, so by the structural lemma its length $z=j-l-1$ is at most
  $d-2$. Put $v=\min(z,e)$ and
  $B=(j,\ l+1,\dots,l+v,\ k-d+v+1,\dots,k-d+e)$ (last interval empty
  if $v=e$).

In both cases $|B|=h$ and the entries are distinct: in the second
case, if $z<e$ then $j\le k-d<k-d+v+1$ separates the blocks; if
$z\ge e$ there is no final interval. Pairing rows $A=\{0,\dots,e\}$
with $B$ in the listed order, $T[A,B]$ is lower triangular with
nonzero diagonal: in the first case the diagonal is
$a_\ell,a_{k-d},\dots,a_{k-d}$ (note $a_{k-d}=1$) and every entry above
it has coefficient index in the standard run $[k-d+1,k-1]$. In the
second case the diagonal is $a_j$, then $a_l$ repeated ($v$ times,
from columns $l+1,\dots,l+v$ paired with rows $1,\dots,v$), then
$a_{k-d}$ repeated; entries above the diagonal have indices either in
the chosen gap $(l,j)$ or in the standard run. A lower-triangular
matrix's permanent is the product of its diagonal (any other
bijection uses an above-diagonal entry), so
$\mathrm{perm}\,T_i[A_i,B_i]\ne0$ for every $i$.

*Killing $\rho\ge1$.* A column $b\in B_i$ becomes an all-zero column
of $T_i[A_i,B_i-\rho]$ whenever
$\rho\in[b+1,\,b+\delta_i]\pmod k$: indeed the entries of that column
are $a_{(i,\,b-\rho-u)}$ for $u\in A_i=\{0,\dots,e_i\}$, and
$b-\rho$ ranges over $[-\delta_i,-1]$, so $b-\rho-u$ ranges inside
$[k-\delta_i-e_i,\,k-1]=[k-d_i+1,\,k-1]$, the standard run, where all
$a_{(i,\cdot)}$ vanish. The intervals $[b+1,b+\delta_i]$ over all
$b\in B_i$ and all $i$ cover every $\rho\in\{1,\dots,k-1\}$: source
$1$ has $\ell_1=0$ and $a_0\ne0$ by the structural lemma, so it is in
the first case; its column $0$ covers $[1,\delta_1]$ and its columns
$k-d_1+1,\dots,k-d_1+e_1$ have intervals
$[k-d_1+j+1,\,k-d_1+j+\delta_1]$ ($j=1,\dots,e_1$) whose union is the
tail $[k-d_1+2,k-1]$ (consecutive intervals overlap since
$\delta_1\ge1$, and the last ends at $k-d_1+e_1+\delta_1=k-1$); for $1<i\le J$,
if $a_{\ell_i}\ne0$ then $\ell_i\in B_i$ covers
$[\ell_i+1,\ell_i+\delta_i]$; if $a_{\ell_i}=0$ and $z<e_i$, then
$B_i\supseteq[l+1,j]\ni\ell_i$, same conclusion; if $z\ge e_i$, the
columns $l+1,\dots,l+e_i$ cover $[l+2,l+e_i+\delta_i]$ contiguously
and column $j$ extends coverage through $j+\delta_i\ge
\ell_i+\delta_i$ (the two intervals meet because
$j+1=l+z+2\le l+e_i+\delta_i$, using $z\le d_i-2$ from the
structural lemma), while $l+2\le\ell_i+1$, so
$[\ell_i+1,\ell_i+\delta_i]$ is covered. Together the first $J$
sources cover $[1,\sum_{i\le J}\delta_i]$, which meets source $1$'s
tail because $\sum_{i\le J}\delta_i\ge k-d_1+1$. Sources after $J$ are
irrelevant. Hence for each $\rho\ge1$ some factor in (1) has a zero
column, so its permanent is $0$.

Therefore only the $\rho=0$ summand of (1) survives; it is $\mu$ times
a product of nonzero permanents, hence nonzero. $\square$

## 5. Decomposing an affine subspace into blocks

Let $A\subseteq K^r$ be affine with direction space $\tilde A$. For
$i=1,\dots,r$ let $V_i$ be the projection to block $i$ of
$\{v\in\tilde A: v_1=\dots=v_{i-1}=0\}$, and choose, for each $i$, an
$\mathbf F_2$-linear section $s_i:V_i\to\tilde A$ of that projection
(choose a basis). Fixing $\gamma\in A$, every point of $A$ is uniquely
$\gamma+\sum_i s_i(z_i)$, $z_i\in V_i$: read off $z_1$ from block $1$,
subtract, and continue; uniqueness is triangularity, and counting
gives $\dim A=\sum_i\dim V_i$. Writing $d_i=\dim V_i$ and
$\sigma_{ij}:V_j\to K$ for block $i$ of $s_j$ ($j<i$),

$$A=\Bigl\{\bigl(\gamma_i+z_i+\textstyle\sum_{j<i}\sigma_{ij}(z_j)
\bigr)_{i=1}^r : z_i\in V_i\Bigr\}. \tag{3}$$

Let $Q_i$ be the image polynomial of $V_i$ (the zero polynomial as a
function if $V_i=\{0\}$; only nonzero $V_i$ are used below), and let
$Q_{ij}$ be a reduced linearized polynomial representing the linear
map $\sigma_{ij}\circ Q_j:K\to K$ (Section 1). Set
$R_i(X)=\gamma_i+Q_i(X_i)+\sum_{j<i}Q_{ij}(X_j)$. The map
$R:K^r\to K^r$ has image exactly $A$.

## 6. Selecting $t$ high-dimension blocks

Assume $\dim A=D>kr/t+kt+r(r+1)$. Let $S=\{i:d_i>r+1\}$. Since
$D=\sum_id_i\le|S|\,k+r(r+1)$ (each $d_i\le k$), we get
$|S|>(D-r(r+1))/k>t$, so $|S|\ge t$. Also
$\sum_{i\in S}d_i\ge D-r(r+1)>kr/t+kt$. Choose $I_0\subseteq S$,
$|I_0|=t$, with the largest $d_i$; since $|S|\le r$,
$\sum_{i\in I_0}d_i\ge\frac{t}{|S|}\sum_{i\in S}d_i\ge
\frac{t}{r}(kr/t+kt)=k+kt^2/r$, hence

$$\sum_{i\in I_0}\bigl(d_i-(r+1)\bigr)\;\ge\;
k+\frac{kt^2}{r}-t(r+1)\;\ge\;k, \tag{4}$$

using $kt\ge r(r+1)$, i.e. $kt^2/r\ge t(r+1)$.
For $i\in I_0$ set $h_i=w_i=r+2-i$: then $1\le e_i=r+1-i$ (as
$i\le r$), $h_i\le r+1<d_i$, so the claim of Section 4 applies to the
sources $V_i$, $i\in I_0$, once $(\dagger)$ holds:
$\sum_{i\in I_0}(d_i-h_i)\ge\sum_{i\in I_0}(d_i-(r+1))>k\ge
k-\max_{i\in I_0}d_i+1$. The claim therefore gives a monomial
$M=\prod_{i\in I_0}X_i^{b_i}$, $\mathrm{wt}(b_i)=w_i$, with nonzero
coefficient in
$\mathrm{Tr}\bigl(\mu\prod_{i\in I_0}Q_i(X_i)^{2^{w_i}-1}\bigr)$.

## 7. Mixing and isolation

**Mixing does not change the coefficient.** $M$ has the same
coefficient in
$\mathrm{Tr}\bigl(\mu\prod_{i\in I_0}R_i(X)^{2^{w_i}-1}\bigr)$.
Indeed expand each $R_i^{2^{w_i}-1}$ as a product of $w_i$ Frobenius
factors $R_i^{2^s}$, $s=0,\dots,w_i-1$; a factor of block $i$ mentions
only variables $X_j$ with $j\le i$. Process the selected indices in
decreasing order. For the largest selected $i$, no other selected
block's factors mention $X_i$ except block $i$'s own (blocks below
$i$ do not reach $X_i$; there is no selected block above $i$), and
reaching binary weight $w_i$ in $X_i$ forces every factor of block
$i$ to choose its $Q_i(X_i)$ term; those terms mention no other
variable. Inductively, when processing selected $i$, all selected
blocks above $i$ have already been forced onto their own diagonal
terms and contribute nothing to $X_i$, blocks below $i$ do not
mention $X_i$, so again all $w_i$ factors of block $i$ must choose
$Q_i(X_i)$ terms. Hence in every contributing choice the whole
product agrees with the unmixed product, for each trace summand
separately. (Constants $\gamma_i$ and mixing terms $Q_{ij}$ never
survive in a contributing choice.)

**Isolation.** No summand of $\mathrm{Sym}^t$ other than $I_0$
contributes to $M$. A summand is indexed by a $t$-subset $J$ and
equals $\mathrm{Tr}\bigl(\mu\prod_{j\in J}R_j(X)^{2^{w_j}-1}\bigr)$.
In an expansion choice let $n_{ji}$, for $j\in J$ and $i\in I_0$, be
the number of Frobenius factors of block $j$ that choose a monomial
summand mentioning $X_i$ (a mixing term $Q_{ji}(X_i)$ is itself a sum
of powers of $X_i$; each chosen summand contributes exactly one power
of two to the $X_i$-degree before carries). A contributing choice must satisfy

$$\sum_{i\in I_0}n_{ji}\le w_j,\qquad
\sum_{j\in J}n_{ji}\ge w_i\ \ (i\in I_0),\qquad
n_{ji}>0\Rightarrow i\le j, \tag{5}$$

the second by binary-degree subadditivity (Section 1: carries only
lower binary weight, also after reduction). Since $w_i$ is strictly
decreasing in $i$,

$$t\;\le\;\sum_{i\in I_0,\,j\in J}\frac{n_{ji}}{w_i}\;\le\;
\sum_{i\in I_0,\,j\in J}\frac{n_{ji}}{w_j}\;\le\;t, \tag{6}$$

the first inequality from the demands in (5), the last from the
supplies. If some $n_{ji}>0$ had $i<j$, then $w_i>w_j$ and the middle
inequality would be strict — impossible. So every chosen term has
$i=j$; every $i\in I_0$ has positive demand, hence $i\in J$; equal
cardinalities give $J=I_0$.

**Conclusion.** The coefficient of $M$ in the reduced polynomial of
$F\circ R$ equals the nonzero coefficient of Section 6. As
$\mathrm{wt}(b_i)=w_i\ge2$, $M$ has positive degree, so $F\circ R$ is
nonconstant; since $R(K^r)=A$, $F$ is nonconstant on $A$. This proves
Lemma D. $\square$

## 8. What is used from outside this file

Only the definition-level facts of Section 1 (unique reduced
representatives, linearized polynomials = linear maps, trace duality,
binary-weight bookkeeping), each proved above. The existence of a
field $K=\mathbf F_{2^k}$ for the specific $k$ and the evaluation
costs are supplied in `proof.md`; this file needs no irreducible
polynomial, no search procedure, and no choice: it speaks about an
arbitrary field of order $2^k$ and an arbitrary nonzero $\pi$.
