# The Ben-Sasson–Kopparty affine disperser, reproduced

This file states and proves the block affine-disperser theorem used in
`proof.md`, following E. Ben-Sasson and S. Kopparty, *Affine dispersers from
subspace polynomials* (preprint 12 March 2010; SIAM J. Comput. 41(4):880–914,
2012), Theorems 2.2 and 7.3 and Sections 3, 6, 7. We specialise to $p=2$, which
is the only case needed. Everything is proved from elementary finite-field
algebra; the only external computational statement is Shoup's theorem
(Theorem 12 below), used solely for uniform polynomial-time computability.

Throughout, $k$ is a prime, $F=\mathbb F_{2^k}$, and we identify
$\mathbb F_2^k$ with $F$ through a fixed $\mathbb F_2$-basis (ultimately the
polynomial basis of `proof.md`). We write $2^k=q$ often, but keep $k$ for the
field degree to avoid confusion with the block count $r$. The **trace** is
$\operatorname{Tr}(y)=\sum_{j=0}^{k-1}y^{2^j}:F\to\mathbb F_2$; it is
$\mathbb F_2$-linear and nonzero, and every nonzero $\mathbb F_2$-linear map
$F\to\mathbb F_2$ is $y\mapsto\operatorname{Tr}(\beta y)$ for a unique
$\beta\in F\setminus\{0\}$.

A polynomial is **linearized** if it has the form
$\sum_{j=0}^{k-1}a_jX^{2^j}$; it is $\mathbb F_2$-linear as a map $F\to F$.
It is **affine linearized** if it is linearized plus a constant. Reduction
"mod $(X^{2^k}-X)$" means replacing $X^{2^k}$ by $X$ repeatedly, which does
not change the values on $F$.

---

## 1. Subspace polynomials

For an affine subspace $A\subseteq F$ of dimension $d$ (as an
$\mathbb F_2$-vector space), the **kernel-subspace polynomial** is
$$
P_A(X)=\prod_{\alpha\in A}(X-\alpha).
$$
It is monic, affine linearized, of degree $2^d$, and its roots are exactly
$A$. The following is the key existence result.

**Lemma 1 (image-subspace polynomial; Lemma 3.7 of the source).**
For every affine subspace $A\subseteq F$ of dimension $d$ there is a monic
affine-linearized polynomial $Q_A$ of degree $2^{k-d}$ such that
$Q_A(F)=A$. Moreover
$$
P_A(Q_A(X))=Q_A(P_A(X))=X^{2^k}-X.
\tag{1.1}
$$
Consequently $Q_A$ has exactly $2^{k-d}$ distinct roots in $F$, and those
roots form an affine subspace of dimension $k-d$.

*Proof.* Let $A_0=P_A(F)$, the image of $P_A$. Since $P_A$ is affine
linearized, $A_0$ is an affine subspace; its dimension is $k-d$ because the
kernel of $P_A$ is $A$. Define $Q_A=P_{A_0}$, the kernel-subspace polynomial
of $A_0$; it is monic, affine linearized, of degree $2^{k-d}$.
Then $Q_A(P_A(X))$ is monic of degree $2^k$ and vanishes on $F$, hence equals
$X^{2^k}-X$. Thus $P_A(Q_A(P_A(X)))=P_A(X^{2^k}-X)=P_A(X)^{2^k}-P_A(X)$.
Let $g(Y)=P_A(Q_A(Y))-(Y^{2^k}-Y)$. We have $g(P_A(X))=0$. Since
$\deg(g(P_A(X)))=\deg g\cdot\deg P_A$ unless $g=0$, we get $g=0$, i.e.
$P_A(Q_A(Y))=Y^{2^k}-Y$ for all $Y\in F$. In particular $P_A(Q_A(y))=0$ for
every $y\in F$, so $Q_A(F)\subseteq A$. The kernel of $Q_A=P_{A_0}$ is $A_0$,
of dimension $k-d$, so its image has dimension $k-(k-d)=d$; since $A$ also has
dimension $d$, $Q_A(F)=A$. The root set of $Q_A$ is $A_0$, of dimension $k-d$.
$\square$

For a **linear** subspace $L$ we write $Q_L(X)=\sum_{j=0}^{k-1}a_jX^{2^j}$
(no constant term). The next lemma is the structural fact that drives
everything.

**Lemma 2 (Lemma 3.10 of the source).** Let $L\subseteq F$ be a linear
subspace of dimension $d$ and $Q_L(X)=\sum_{j=0}^{k-1}a_jX^{2^j}$ its
image-subspace polynomial.

1. For every integer $r$ and every set $J=\{(r+j)\bmod k:0\le j\le d-1\}$ of
   $d$ consecutive indices, some $a_j$ with $j\in J$ is nonzero. In
   particular $a_0\neq0$ and $a_{k-d}\neq0$.
2. Suppose $L$ is not contained in any constant multiple of a proper subfield
   of $F$, i.e. $L\not\subseteq c\,\mathbb F_{2^{k'}}$ for every $c\in F$ and
   every proper divisor $k'$ of $k$. Then for every integer $r\neq k-d+1$ and
   every set $J=\{(r+j)\bmod k:0\le j\le d-2\}$ of $d-1$ consecutive
   indices, some $a_j$ with $j\in J$ is nonzero.

*Proof.* (1) Suppose $a_j=0$ for all $j\in J$. By Lemma 1, $Q_L$ has exactly
$2^{k-d}$ distinct roots in $F$. Let $t\equiv k-(r+d)\pmod k$ and
$Q'(X)=Q_L(X)^{2^t}\bmod(X^{2^k}-X)$. Applying
$(Q_L(X))^{2^t}\equiv\sum_j a_{(j-t)\bmod k}^{2^t}X^{2^j}\pmod{X^{2^k}-X}$
(which follows by expanding $Q_L^{2^t}$), the coefficient of $X^{2^j}$ for
$j\in\{k-d,\dots,k-1\}$ is $a_{(j+r+d)\bmod k}$; writing $j=k-d+s$ with
$0\le s\le d-1$, this index is $(s+r)\bmod k$, which lies in
$J=\{r,\dots,r+d-1\}\bmod k$, hence is zero; so $Q'$ has degree at most
$2^{k-d}-1$. But for every root $\alpha$ of $Q_L$,
$Q'(\alpha)=Q_L(\alpha)^{2^t}=0$, so $Q'$ has at least $2^{k-d}$ roots, a
contradiction for the nonzero polynomial $Q'$ (nonzero because
$\mathbb F_2[X]/(X^{2^k}-X)$ is a product of fields, hence reduced, and
$Q_L\neq0$). Since an image-subspace polynomial has degree exactly
$2^{k-d}$, we already know $a_j=0$ for $j>k-d$ and $a_{k-d}\neq0$; applying the
just-proved statement to the $d$ consecutive indices
$\{k-d+1,\dots,k-1,0\}$ (take $r=k-d+1$) then gives $a_0\neq0$.

(2) Suppose $a_j=0$ for all $j\in J$ with $|J|=d-1$. If $r>k-d+1$, then $J$
contains $0$ and $a_0\neq0$ by part (1), so the claim is immediate. Hence
assume $r\le k-d$ (the case $d=1$ is vacuous). Let $h=k-(r+d)+1$; then
$1\le h<k$ because $r\le k-d$ and $r+d-1\ge1$. Put
$Q'(X)=Q_L(X)^{2^h}\bmod(X^{2^k}-X)$. Then $Q'$ is nonzero of degree at most
$2^{k-d}$, and it vanishes on all $2^{k-d}$ roots of $Q_L$; since a nonzero
polynomial of degree $<2^{k-d}$ cannot have $2^{k-d}$ roots, $Q'$ has degree
exactly $2^{k-d}$. Both $Q'$ and $Q_L$ have degree $2^{k-d}$ and vanish
exactly on the same $2^{k-d}$ distinct roots (the kernel of $Q_L$), so each is
a scalar multiple of $\prod_{\alpha\in\ker Q_L}(X-\alpha)$; since $Q_L$ is
monic, $Q'=c\,Q_L$ for some $c\in F\setminus\{0\}$. Hence
$Q_L(X)^{2^h}=c\,Q_L(X)\pmod{X^{2^k}-X}$. For $\alpha\in F$, this reads
$Q_L(\alpha)^{2^h}=c\,Q_L(\alpha)$. Since $Q_L(F)=L$, we get
$\lambda^{2^h}=c\lambda$ for all $\lambda\in L$. The multiplicative group
$F^*$ has order $2^k-1$, and the map $x\mapsto x^{2^h-1}$ is surjective
because $\gcd(2^h-1,2^k-1)=2^{\gcd(h,k)}-1=1$ (here $k$ is prime and
$0<h<k$). Hence $c=b^{2^h-1}$ for some $b\in F^*$, and
$\lambda^{2^h}=c\lambda$ forces $(\lambda/b)^{2^h-1}=1$, i.e.
$\lambda\in b\,\mathbb F_{2^h}$. We obtain $L\subseteq b\,\mathbb F_{2^h}$.
Intersecting with $F$ gives $L\subseteq b\,\mathbb F_{2^{\gcd(h,k)}}=b\,\mathbb
F_2$, a constant multiple of a proper subfield, contradicting the hypothesis.
$\square$

---

## 2. Coefficients of products of linearized polynomials

For a linearized $Q(X)=\sum_{i=0}^{k-1}a_iX^{2^i}$ define its **associated
matrix** $M_Q\in F^{k\times k}$ by $M_Q(i,j)=(a_{j-i})^{2^i}$, indices mod $k$.
Then the $(i,j)$-entry of $M_Q$ is the coefficient of $X^{2^j}$ in
$Q(X)^{2^i}\bmod(X^{2^k}-X)$. Let $M_Q^0$ be the $0/1$ matrix recording which
entries are nonzero. For sets $A,B\subseteq\{0,\dots,k-1\}$ let $M[A,B]$ be
the submatrix, and $B+r=\{(b+r)\bmod k:b\in B\}$.

**Lemma 3 (Lemma 3.13 of the source, with the roles of the two index sets
named explicitly).** Let $\beta\in F\setminus\{0\}$. For each $i$ let
$P_i,M_i\subseteq\{0,\dots,k-1\}$ be the **power set** and the **monomial
set**, with $|P_i|=|M_i|>0$; put
$$
\tau_i=\sum_{l\in P_i}2^l\quad(\text{the exponent applied to }Q_i),\qquad
\sigma_i=\sum_{j\in M_i}2^j\quad(\text{the exponent of the monomial}).
$$
Let $Q_1,\dots,Q_t$ be linearized and
$$
R(X_1,\dots,X_t)=\operatorname{Tr}\Bigl(\beta\prod_{i=1}^t Q_i(X_i)^{\tau_i}
\Bigr)\bmod\bigl(X_i^{2^k}-X_i\bigr)_{i\in[t]}.
$$
Then the coefficient $c_M$ of $M=\prod_{i=1}^tX_i^{\sigma_i}$ in $R$ is
$$
c_M=\sum_{r=0}^{k-1}\beta^{2^r}\prod_{i=1}^t
\operatorname{Perm}\bigl(M_{Q_i}[P_i+r,M_i]\bigr)
=\sum_{r=0}^{k-1}\beta^{2^r}\prod_{i=1}^t
\operatorname{Perm}\bigl(M_{Q_i}[P_i,M_i-r]\bigr)^{2^r}.
\tag{2.1}
$$

*Proof.* By definition of the trace,
$\operatorname{Tr}(\beta\prod_iQ_i^{\tau_i})=\sum_{r=0}^{k-1}
\beta^{2^r}\prod_iQ_i^{\tau_i2^r}$. For fixed $r$,
$Q_i^{\tau_i2^r}=\prod_{l\in P_i}Q_i^{2^{l+r}}$. Expanding each factor
$Q_i^{2^{l+r}}$, a contribution to $X_i^{\sigma_i}$ corresponds to a bijection
$h:P_i\to M_i$ choosing from the factor indexed by $l$ the monomial
$X_i^{2^{h(l)}}$; its coefficient is
$\prod_{l\in P_i}(a_{i,h(l)-(l+r)})^{2^{l+r}}$, which is the permanent of the
submatrix $M_{Q_i}[P_i+r,M_i]$ (rows indexed by $P_i+r$, columns by $M_i$).
Summing over $r$ gives the first form; using
$(a_{j-i})^{2^{i+r}}$ for $i=l+r$, the same product equals
$\operatorname{Perm}(M_{Q_i}[P_i,M_i-r])^{2^r}$, giving the second. $\square$

**Remark 4 (Remark 3.14 of the source).** In the notation of Lemma 3, suppose
that (i) each $M_{Q_i}^0[P_i,M_i]$ is, after reordering rows and columns,
triangular with nonzero diagonal, and (ii) for every $r\in\{1,\dots,
k-1\}$ there is some $i$ such that $M_{Q_i}^0[P_i,M_i-r]$ has an all-zero
column. Then $c_M\neq0$: the $r=0$ term of (2.1) is a product of nonzero
permanents of triangular matrices (a triangular matrix has permanent equal to
the product of its diagonal entries), and every $r\ge1$ term has a zero
permanent.

---

## 3. The independent-source theorem

**Theorem 5 (Theorem 6.1 of the source, $p=2$).** Let
$A_1,\dots,A_t\subseteq F$ be **linear** subspaces of dimensions
$d_1,\dots,d_t>1$, none contained in a constant multiple of a proper subfield
of $F$. Let $Q_i$ be the image-subspace polynomial of $A_i$. Let
$\beta\in F\setminus\{0\}$ and let $1\le e_i<d_i-1$. Put
$s_i=\sum_{j=0}^{e_i}2^j$ and
$$
R(X_1,\dots,X_t)=\operatorname{Tr}\Bigl(\beta\prod_{i=1}^tQ_i(X_i)^{s_i}
\Bigr)\bmod\bigl(X_i^{2^k}-X_i\bigr)_{i\in[t]}.
$$
If $\sum_{i=1}^t\bigl(d_i-(e_i+1)\bigr)>k-\max_i d_i+1$, then $R$ has a
monomial $M=\prod_{i=1}^tX_i^{\sigma_i}$ with
$\operatorname{wt}_2(\sigma_i)=e_i+1$ whose coefficient is nonzero. In
particular $R$ is nonconstant on $F^t$.

Because $k$ is prime, the only proper subfield of $F=\mathbb F_{2^k}$ is
$\mathbb F_2$; hence every linear subspace of dimension at least $2$ (in
particular every $A_i$ used below) automatically satisfies the
non-containment hypothesis.

(Here $\operatorname{wt}_2$ is the binary digit sum, so
$\operatorname{wt}_2(\sigma_i)=|M_i|$ whenever $\sigma_i=\sum_{j\in M_i}2^j$
is a sum of distinct powers.)

*Proof.* Let $D_i=\{0,1,\dots,e_i\}$ (the **power set** of exponents applied
to $Q_i$; this is distinct from the subspace $A_i$), so $|D_i|=e_i+1$. It
suffices to construct **monomial sets** $B_1,\dots,B_t$ with
$|B_i|=e_i+1$ and exponents $\sigma_i=\sum_{l\in B_i}2^l$ such that
the monomial $M=\prod_iX_i^{\sigma_i}$ has nonzero coefficient. We apply
Lemma 3 and Remark 4 with $P_i=D_i$ (power set) and $M_i=B_i$ (monomial set).

Write $Q_i(X)=\sum_{j=0}^{k-1}a_{i,j}X^{2^j}+\hat a_i$. Define
$\ell_1=0$ and $\ell_{i+1}=\ell_i+d_i-(e_i+1)\pmod k$ for $1\le i<t$; thus
$\ell_i=\sum_{i'<i}(d_{i'}-(e_{i'}+1))\pmod k$.

*Case $a_{i,\ell_i}\neq0$:*
$$
B_i=\{\ell_i\}\cup\{k-d_i+1,k-d_i+2,\dots,k-d_i+e_i\}.
$$
*Case $a_{i,\ell_i}=0$:* let
$j_i=\min\{j>\ell_i:a_{i,j}\neq0\}$ and
$\hat\jmath_i=\max\{j<\ell_i:a_{i,j}\neq0\}$; these exist by Lemma 2(1).
Let $q_i=j_i-\hat\jmath_i-1$ and $s_i'=\min\{q_i,e_i\}$. Put
$$
B_i=\{j_i\}\cup[\hat\jmath_i+1,\hat\jmath_i+s_i']
\cup\{k-d_i+s_i'+1,\dots,k-d_i+e_i\}.
$$

**Claim 5.1 (Claim 6.2).** For every $i$, $M_{Q_i}^0[D_i,B_i]$ is lower
triangular with nonzero diagonal.

*Proof.* An entry of $M^0[D_i,B_i]$ in row $u\in D_i$ and column $b$ records
whether $a_{i,b-u}$ is nonzero.

*First case.* Order the columns as $\ell_i,k-d_i+1,\dots,k-d_i+e_i$. An entry
strictly above the diagonal is a pair $(u,v)$ with $u<v$. For $v=0$ there is
no such $u\ge0$. For $v\ge1$ the column index is $b=k-d_i+v$, so
$b-u\ge k-d_i+v-u\ge k-d_i+1$ (because $u\le v-1$), and
$b-u\le k-d_i+e_i\le k-1$; hence $b-u$ lies in the zero block
$\{k-d_i+1,\dots,k-1\}$ of Lemma 2(1) and the entry is zero. The diagonal
entries are $(0,0)$, which is $a_{i,\ell_i}\neq0$, and $(v,v)$ ($v\ge1$),
which is $a_{i,k-d_i}\neq0$ by Lemma 2(1). Thus the minor is lower triangular
with nonzero diagonal.

*Second case.* The leading $(s_i'+1)\times(s_i'+1)$ block, with columns
$j_i,\hat\jmath_i+1,\dots,\hat\jmath_i+s_i'$ and rows $0,\dots,s_i'$, is lower
triangular: for $v\ge1$ and $u<v$, the entry index is
$\hat\jmath_i+v-u\in[\hat\jmath_i+1,\hat\jmath_i+s_i']\subseteq
[\hat\jmath_i+1,j_i-1]$ (as $s_i'\le q_i=j_i-\hat\jmath_i-1$), and all
coefficients with indices strictly between $\hat\jmath_i$ and $j_i$ vanish by
the definitions of $j_i,\hat\jmath_i$. Its diagonal entries are
$a_{i,j_i}\neq0$ and $a_{i,\hat\jmath_i}\neq0$. The remaining
$e_i-s_i'$ columns are the last columns of the first case, and the same
argument applies. $\square$

**Claim 5.2 (Claim 6.3).** For every $r\in\{1,\dots,k-1\}$ there is an $i$
such that $M_{Q_i}^0[D_i,B_i-r]$ has an all-zero column.

*Proof.* For $m\in B_i$, the entries of row $u\in D_i$ in the column $m-r$ of
$M^0[D_i,B_i-r]$ are $a_{i,(m-r)-u}$, $u=0,\dots,e_i$. If
$r\in[m+1,\,m+d_i-(e_i+1)]\pmod k$, then as $u$ runs over $0,\dots,e_i$ the
indices $(m-r)-u\bmod k$ run over
$[m-r-e_i,\,m-r]\pmod k\subseteq[k-d_i+1,\,k-1]$, the zero block of
Lemma 2(1); so that column is zero. Thus it suffices to show that the intervals
$$
\bigcup_{i=1}^{t}\bigcup_{m\in B_i}[m+1,\,m+d_i-(e_i+1)]\pmod k
\tag{5.4}
$$
cover $\{1,\dots,k-1\}$. We may assume $d_1=\max_i d_i$.

First, $B_1=\{0\}\cup\{k-d_1+1,\dots,k-d_1+e_1\}$. The element $m=0$ gives the
interval $[1,d_1-(e_1+1)]$. The elements $m=k-d_1+v$, $v=1,\dots,e_1$, give
the intervals $[k-d_1+v+1,\,k-d_1+v+d_1-(e_1+1)]$, whose union is
$[k-d_1+2,\,k-1]$.

Choose the least $t_0$ with
$\sum_{i=1}^{t_0}(d_i-(e_i+1))\ge k-d_1+1$; it exists because
$\sum_i(d_i-(e_i+1))>k-\max_i d_i+1=k-d_1+1$ by hypothesis. For
$2\le i\le t_0$, $\ell_{i+1}=\ell_i+d_i-(e_i+1)\pmod k$, and we claim that
$\bigcup_{m\in B_i}[m+1,m+d_i-(e_i+1)]$ contains $[\ell_i+1,\ell_{i+1}]$
(the source's equation (18)). The verification splits into two cases.

If $a_{i,\ell_i}\neq0$, then $\ell_i\in B_i$, and the interval starting at
$\ell_i$ is exactly $[\ell_i+1,\ell_i+d_i-(e_i+1)]=[\ell_i+1,\ell_{i+1}]$.

If $a_{i,\ell_i}=0$, then either $q_i<e_i$, in which case $B_i$ can be
rewritten as $[\hat\jmath_i+1,j_i]\cup\{k-d_i+q_i+1,\dots,k-d_i+e_i\}$ with
$j_i=\hat\jmath_i+q_i+1$, so $\ell_i\in B_i$ and we are back in the first
case; or $q_i\ge e_i$, in which case
$B_i=[\hat\jmath_i+1,\hat\jmath_i+e_i]\cup\{j_i\}$. Since $i\le t_0$ and
$d_1$ is maximal, $\hat\jmath_i<\ell_i\le k-d_1\le k-d_i$, so
$\hat\jmath_i+1\neq k-d_i+1$. The coefficients with indices
$\hat\jmath_i+1,\dots,j_i-1$ vanish, so the zero run has length
$q_i=j_i-\hat\jmath_i-1$. If $q_i\ge d_i-1$, then the $d_i-1$ consecutive
indices $\hat\jmath_i+1,\dots,\hat\jmath_i+d_i-1$ (whose start is not the
excluded alignment $k-d_i+1$) would all be zero, contradicting Lemma 2(2)
with the non-containment hypothesis. Hence $q_i\le d_i-2$, i.e.
$j_i-\hat\jmath_i\le d_i-1$; since also $j_i>\ell_i>\hat\jmath_i$ gives
$j_i\ge\hat\jmath_i+2$, we have
$j_i\in[\hat\jmath_i+2,\hat\jmath_i+d_i-1]$. Taking all but the last element
of $B_i$, the intervals
$\bigcup_{m\in[\hat\jmath_i+1,\hat\jmath_i+e_i]}[m+1,m+d_i-(e_i+1)]$ contain
$[\hat\jmath_i+2,\hat\jmath_i+d_i-1]\ni j_i$; reinserting $j_i$ as the last
element of $B_i$ covers $[\ell_i+1,\ell_{i+1}]$.

Thus $[1,\ell_{t_0+1}]\supseteq[1,k-d_1+1]$, and together with
$[k-d_1+2,k-1]$ this covers $\{1,\dots,k-1\}$. This proves (5.4). $\square$

With Claims 5.1 and 5.2, Remark 4 gives $c_M\neq0$ for the monomial whose
exponents are the sums of the powers indexed by $B_i$. Since each $B_i$ has
$e_i+1$ elements, $\operatorname{wt}_2(\sigma_i)=e_i+1$. $\square$

---

## 4. Bourgain's decomposition and invariance

**Lemma 6 (Lemma 7.1 of the source).** Let $A\subseteq F^r$ be an affine
subspace and $\xi\in A$. There exist linear subspaces $Y_1,\dots,Y_r
\subseteq F$ and linear maps $\psi_{ij}:Y_j\to F$ ($j<i$) such that
$$
A=\Bigl\{(x_1,\dots,x_r):\exists y_i\in Y_i,\
x_i=\xi_i+y_i+\sum_{j<i}\psi_{ij}(y_j)\Bigr\},
\qquad \dim A=\sum_{i=1}^r\dim Y_i.
\tag{4.1}
$$
*Proof.* Put the rows of a basis of the linear subspace underlying $A$ in
echelon form with respect to the block coordinates $1,\dots,r$. The pivot
elements of block $i$ form $Y_i$; the earlier free coordinates are expressed
through $\psi_{ij}$. The number of pivots is $\sum_i\dim Y_i=\dim A$. $\square$

**Lemma 7 (Lemma 7.2 of the source).** Let $P_i(X_i)$ be linearized
($i\in[r]$) and $P_{ij}(X_j)$ linearized for $j<i$. Let
$\xi\in F^r$, let $I_0=\{i_1<\dots<i_{t'}\}\subseteq[r]$, let
$e_{i_1},\dots,e_{i_{t'}}\ge1$, put $s_i=\sum_{l=0}^{e_i}2^l$, and let
$\beta\in F\setminus\{0\}$. Define
$$
g=\operatorname{Tr}\Bigl(\beta\prod_{i\in I_0}\bigl(P_i(X_i)
+\sum_{j<i}P_{ij}(X_j)+\xi_i\bigr)^{s_i}\Bigr)\bmod(X_i^{2^k}-X_i)_i,
$$
$$
g_0=\operatorname{Tr}\Bigl(\beta\prod_{i\in I_0}P_i(X_i)^{s_i}\Bigr)
\bmod(X_i^{2^k}-X_i)_i.
$$
Then for every choice of exponents $(\sigma_i)_{i\in I_0}$ with
$\operatorname{wt}_2(\sigma_i)=e_i+1$, the coefficient of
$\prod_{i\in I_0}X_i^{\sigma_i}$ in $g$ equals its coefficient in $g_0$.

*Proof.* Write $\hat P_i(X_{<i})=\sum_{j<i}P_{ij}(X_j)+\xi_i$, and expand
$$
g=\sum_{r=0}^{k-1}\beta^{2^r}\prod_{i\in I_0}
\bigl(P_i(X_i)+\hat P_i(X_{<i})\bigr)^{s_i2^r}.
$$
For a fixed $i$ and $r$, split
$\bigl(P_i+\hat P_i\bigr)^{s_i2^r}$ according to how many of the factors
$P_i^{2^{r+l}}$ versus $\hat P_i^{2^{r+l}}$ are chosen. The term using all
factors $P_i$ has $\mathbb F_2$-degree $e_i+1$ in $X_i$; every term using at
least one $\hat P_i$ factor has $\mathbb F_2$-degree at most $e_i$ in $X_i$
(because $P_i$ and $\hat P_i$ are linearized, hence of $\mathbb F_2$-degree
$1$, and $\hat P_i$ does not involve $X_i$). Since the target monomial has
$\mathbb F_2$-degree exactly $e_i+1$ in each $X_i$, only the all-$P_i$
terms can contribute. Multiplying over $i$ and summing over $r$ leaves the
coefficient in $g_0$. $\square$

---

## 5. The block theorem

**Theorem 8 (Theorem 7.3 of the source, $p=2$).** Let $t<r$ be integers and
$k$ be prime with $k\ge r(r+1)/t$. Put $e_i=r+1-i$,
$s_i=\sum_{j=0}^{e_i}2^j$ ($1\le i\le r$), and let
$\beta\in F\setminus\{0\}$. Define
$$
f(x_1,\dots,x_r)=\operatorname{Tr}\Bigl(\beta\sum_{I\subseteq[r],\,|I|=t}
\prod_{i\in I}x_i^{s_i}\Bigr).
$$
Then for every affine subspace $A\subseteq F^r$ with
$$
\dim A>\frac{kr}{t}+kt+r(r+1),
\tag{5.1}
$$
$f$ is nonconstant on $A$.

*Proof.* Let $d=\dim A$ and fix $\xi\in A$. Apply Lemma 6: there are linear
$Y_i$ and linear $\psi_{ij}$ with (4.1) and $\sum_i d_i=d$, where
$d_i=\dim Y_i$. Let $Q_i$ be the image-subspace polynomial of $Y_i$ (linear,
so linearized), and for $j<i$ let $Q_{ij}$ be the linearized polynomial with
$Q_{ij}(x)=\psi_{ij}(Q_j(x))$. Put
$$
R_i(X_1,\dots,X_r)=Q_i(X_i)+\sum_{j<i}Q_{ij}(X_j)+\xi_i .
$$
Then the map $H=(R_1,\dots,R_r):F^r\to F^r$ has image exactly $A$. Define
$$
h=f\circ H=\operatorname{Tr}\Bigl(\beta\sum_{|I|=t}\prod_{i\in I}R_i^{s_i}
\Bigr)\bmod(X_i^{2^k}-X_i)_i .
$$
By the reduction property, $h(F^r)=f(A)$; it suffices to show that $h$ is
nonconstant.

Let $S=\{i:d_i>r+1\}$. Since $d_i\le k$,
$$
|S|\ge t
\quad\text{and}\quad
\sum_{i\in S}d_i\ \ge\ \sum_{i\in[r]}\bigl(d_i-(r+1)\bigr)
=d-r(r+1)>\frac{kr}{t}+kt .
\tag{5.2}
$$
Hence there is $I_0\subseteq S$ with $|I_0|=t$ and
$$
\sum_{i\in I_0}\bigl(d_i-(r+1)\bigr)
\ \ge\ \frac{t}{r}\sum_{i\in S}d_i-(r+1)t
\ >\ \frac{t}{r}\Bigl(\frac{kr}{t}+kt\Bigr)-(r+1)t
=k+\frac{kt^2}{r}-(r+1)t\ \ge\ k,
\tag{5.3}
$$
where the last inequality is $k\ge r(r+1)/t$.

Consider the summand
$g=\operatorname{Tr}(\beta\prod_{i\in I_0}R_i^{s_i})$. For $i\in I_0\subseteq
S$ we have $d_i>r+1>e_i=r+1-i$ and $e_i+1=r+2-i\le r+1$, so
$$
\sum_{i\in I_0}\bigl(d_i-(e_i+1)\bigr)
=\sum_{i\in I_0}\bigl(d_i-(r+1)\bigr)+\sum_{i\in I_0}(i-1)
\ \ge\ k .
$$
Also $\max_{i\in I_0}d_i\le k$. Theorem 5 (applied to the $Y_i$: they are
linear subspaces of $F=\mathbb F_{2^k}$ with $k$ prime and $d_i>r+1\ge2$, so
they are not contained in a constant multiple of the only proper subfield
$\mathbb F_2$) gives a monomial
$M=\prod_{i\in I_0}X_i^{\sigma_i}$ with $\operatorname{wt}_2(\sigma_i)=e_i+1$
whose coefficient in
$g_0=\operatorname{Tr}(\beta\prod_{i\in I_0}Q_i(X_i)^{s_i})$ is nonzero.
Lemma 7 (with $P_i=Q_i$, $P_{ij}=Q_{ij}$) shows that the coefficient of $M$ in
$g$ is the same, hence nonzero.

It remains to show that no other summand of $h$, indexed by a $t$-set $J\neq
I_0$, contributes to $M$. Each $R_j$ is a polynomial in $X_1,\dots,X_j$ that is
a sum of monomials of total $\mathbb F_2$-degree at most $1$ (linearized
polynomials plus a constant). In the summand indexed by $J$,
$$
\operatorname{Tr}\Bigl(\beta\prod_{j\in J}R_j^{s_j}\Bigr)
=\sum_{w=0}^{k-1}\beta^{2^w}\prod_{j\in J}\ \prod_{l=0}^{e_j}R_j^{2^{l+w}}.
$$
A monomial $\prod_{i\in I_0}X_i^{\sigma_i}$ arising here is obtained by
choosing, for each pair $(j,l)$ with $j\in J$, $0\le l\le e_j$, one monomial
from $R_j^{2^{l+w}}$; each is a single variable $X_i$ (possibly raised to a
power of $2$) or a constant. No chosen non-constant monomial may be supported
on a variable $X_i$ with $i\notin I_0$, because the target $M$ has
$\mathbb F_2$-degree $0$ in those variables; constant choices are allowed and
contribute nothing to the degree. Hence, on any choice that contributes to
$M$, let $\varphi_j(l)\in I_0$ be the index of the chosen variable (if the
choice is constant, it contributes no degree and may be discarded from the
count). Then the $\mathbb F_2$-degree of the product in $X_i$ is at most
$|\varphi_j^{-1}(i)|$ summed over $j\in J$. Since the target has
$\mathbb F_2$-degree $e_i+1=r+2-i$ in each $X_i$, any nonzero contribution
requires
$$
\sum_{j\in J}\bigl|\varphi_j^{-1}(i)\bigr|\ \ge\ r+2-i
\qquad\text{for all }i\in I_0 .
\tag{5.5}
$$
The following combinatorial fact then forces $J=I_0$.

**Claim 8.1 (Claim 7.4).** Let $I_0,J\subseteq[r]$ with $|I_0|=|J|=t$. Suppose
there are functions $\varphi_j:\{0,1,\dots,r+1-j\}\to\{1,\dots,j\}$ for $j\in J$
satisfying (5.5) for every $i\in I_0$. Then $J=I_0$.

*Proof.* For $j\in J$,
$\sum_{i\in I_0}|\varphi_j^{-1}(i)|\le r+2-j$ (the size of the domain).
Summing (5.5) over $i\in I_0$ gives
$\sum_{i\in I_0}(r+2-i)\le\sum_{j\in J}(r+2-j)$, hence
$\sum_{i\in I_0}i\ge\sum_{j\in J}j$. Now fix $m\ge1$. Since
$\varphi_j(l)\le j$, the expression
$\sum_{j\in J}\sum_{l=0}^{r+1-j}\bigl(j^m-\varphi_j(l)^m\bigr)$ is nonnegative;
expanding it and using (5.5) gives
$$
\sum_{j\in J}(r+2-j)j^m\ \ge\ \sum_{i\in I_0}(r+2-i)i^m .
$$
Rearranging,
$\sum_j j^{m+1}-\sum_i i^{m+1}\le(r+2)\bigl(\sum_j j^m-\sum_i i^m\bigr)$;
by induction on $m$ this yields
$\sum_{i\in I_0}i^m\ge\sum_{j\in J}j^m$ for every $m\ge1$.

Let $i_1=\max I_0$ and $j_1=\max J$. Taking $m$ large gives $i_1\ge j_1$. If
$i_1>j_1$, then $\varphi_j(l)\le j\le j_1<i_1$ for all $j\in J$, $l$, so all
preimages of $i_1$ are empty, contradicting (5.5) (as $r+2-i_1\ge2>0$). Hence
$i_1=j_1$. Now (5.5) for $i=i_1$ reads
$\sum_j|\varphi_j^{-1}(i_1)|\ge r+2-j_1$, which equals the size of the domain of
$\varphi_{j_1}$; since all other $\varphi_j$ take values $\le j<j_1=i_1$,
necessarily $\varphi_{j_1}(l)=i_1$ for every $l$ in its domain. Removing
$i_1,j_1$ preserves (5.5) for the remaining $i$, and the power-sum inequalities
subtract $i_1^m=j_1^m$; the argument repeats for the second-largest elements,
then the third-largest, and so on, giving $I_0=J$. $\square$

By Claim 8.1, only the summand $I_0$ can contribute to $M$, and its
coefficient is nonzero as shown. Hence $h$ is nonconstant, so $f$ is
nonconstant on $A$. $\square$

**Corollary 9 (Theorem 2.2, $p=2$).** With the notation of Theorem 8, for any
constant $C>0$ with $D_0:=kr/t+kt+r(r+1)<C(kr)^{4/5}=CM^{4/5}$, the function
$f$ is an affine disperser for dimension $\lfloor CM^{4/5}\rfloor+1$ on
$M=rk$ bits.

*Proof.* By Theorem 8, $f$ is nonconstant on every affine subspace of
$F^r\cong\mathbb F_2^M$ of dimension $>D_0$. Since
$\lfloor CM^{4/5}\rfloor+1\ge CM^{4/5}>D_0$, every affine subspace of dimension
at least $\lfloor CM^{4/5}\rfloor+1$ has dimension $>D_0$. $\square$

---

## 6. The specific family used in `proof.md`

Fix an integer $n\ge1$ and let
$$
k=\text{least prime}>2n^{3/5},\quad r=\lceil n/k\rceil,\quad
t=\lceil\sqrt r\rceil,\quad m=rk .
$$
Pad $x\in\{0,1\}^n$ to $m$ bits by appending zeros. Split the padded string
into $r$ blocks of $k$ bits and, using the polynomial basis of $F=\mathbb
F_{2^k}$, read block $i$ as $x_i\in F$. The function of Definition 1.1 is
$$
f_n(x)=\operatorname{Tr}\Bigl(\beta\sum_{|I|=t}\prod_{i\in I}x_i^{s_i}\Bigr),
\qquad s_i=2^{i+1}-1,
\tag{6.1}
$$
where $\beta\ne0$ is the element representing the constant-coefficient map
$\pi$ (i.e. $\pi(z)=\operatorname{Tr}(\beta z)$). The exponent $s_i$ equals
$1+2+\dots+2^i=\sum_{j=0}^{i}2^j$, i.e. the sequence of Theorem 8 after the
reindexing $i\mapsto r+1-i$; since the symmetric sum is invariant under
permuting blocks, (6.1) is exactly the function of Theorem 8.

We verify the hypotheses of Theorem 8 for all sufficiently large $n$:

* $t<r$ and $k\ge r(r+1)/t$: see Lemma 4.3 of `proof.md`.
* $D_0=kr/t+kt+r(r+1)<12m^{4/5}$: see Lemma 4.3 of `proof.md`, where the
  elementary bounds $r\le 2n/k$, $t\ge\sqrt r$, $k<4n^{3/5}$ are used.

Therefore, for all sufficiently large $n$, $f_n$ is an affine disperser for
dimension $d(n)=\lfloor12m^{4/5}\rfloor+1$ on the $M=m$ padded bits. The
padding is an affine injection $\mathbb F_2^n\to\mathbb F_2^m$ (append a fixed
zero string), so an affine subspace of $\mathbb F_2^n$ of dimension $\ge d(n)$
maps to an affine subspace of $\mathbb F_2^m$ of the same dimension, on which
$f_n$ is nonconstant. This is the affine-disperser property at the *actual*
input length $n$.

---

## 7. Uniform polynomial-time field representation (Shoup)

To turn (6.1) into one fixed polynomial-time family for all $n$, we must
deterministically produce the field $F=\mathbb F_{2^k}$ for each prime $k$ in
time polynomial in $k$.

**Theorem 10 (Shoup 1988/1990, Section 3 and Theorem 3.2 of the paper).**
For each fixed prime $p$ there is a deterministic algorithm which, given $n$
in unary, outputs a monic irreducible polynomial of degree $n$ over
$\mathbb F_p$ in time polynomial in $n$. In particular, for $p=2$, a monic
irreducible $q_k$ of degree $k$ is computable deterministically in time
$\operatorname{poly}(k)$.

The algorithm reduces the problem deterministically in polynomial time to
factoring polynomials over the prime field; for $p=2$ this is Berlekamp's
deterministic factorization. We use Theorem 10 only to fix the representation
$F=\mathbb F_2[Z]/(q_k)$ and to make the evaluator of `proof.md` run in
polynomial time. It does not enter any lower-bound inequality, and the
lower bound of `gate-elimination.md` holds for any field representation.
This is the only external computational premise used; it is a standard
result in computational algebra, stated here with its exact scope.

---

## 8. What is used where

* Lemma 1 (image-subspace polynomial) and Lemma 2 (structural zero pattern)
  are used in Claims 5.1 and 5.2.
* Lemma 3 (permanent coefficient formula) and Remark 4 give the
  noncancellation criterion used in Theorem 5.
* Theorem 5 is applied inside Theorem 8 to the high-entropy block set
  $I_0$.
* Lemma 6 (echelon decomposition) and Lemma 7 (degree/concentration) reduce
  the general affine subspace $A$ to the independent-block setting and
  isolate the unique contributing summand.
* Claim 8.1 (power-sum rigidity) rules out contributions from $J\neq I_0$.
* Corollary 9 and §6 give the affine-disperser property at actual length
  $n$.

All statements are for the full binary field $\mathbb F_2$; no prime-field
generality beyond $p=2$ is needed.
