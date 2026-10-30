# 3.5n − o(n) gates for a fixed total family in P (joint coefficient/exponent proof)

## 0. Result and how to read this submission

We prove, for one fixed uniformly defined total Boolean family $(f_n)_{n\ge 0}$
that lies in P (hence in NP with empty witnesses),

$$
C_{B_2}(f_n)\;\ge\;\frac{7}{2}\,n-7\,d(n)
=\frac{7}{2}\,n-O\!\left(n^{4/5}\right)
\qquad\text{for every sufficiently large integer }n,
$$

where $d(n)=O(n^{4/5})$ is the affine-disperser dimension of $f_n$ defined
below (we use the explicit safe choice $d(n)=\lfloor12\,m(n)^{4/5}\rfloor+1$;
the constant is irrelevant to the coefficient). Consequently

$$
\liminf_{n\to\infty}\frac{C_{B_2}(f_n)}{n}\;\ge\;\frac{7}{2},
$$

in the full Li–Yang $B_2$ model at actual input length $n$. The declared
finite projection is therefore $(c,p)=(7/2,1)$, and the declared full
guarantee is $g(n)=\frac{7}{2}\,n$ (all logarithmic exponents $0$).

The bound $7/2=3.5$ is strictly stronger than the literature seed
$(31/10,1)=(3.1,1)$ of the same challenge; both results are lower bounds of
the form $cn-o(n)$ for one fixed family, and $7/2>31/10$.

The submission has one mathematical engine that is new relative to the seed:

1. **Gate elimination (reproduced in full in `gate-elimination.md`).**
   Theorem 4.1 of Mishra–Santhanam, *An Explicit 3.5n−o(n) Size Lower Bound
   against General Boolean Circuits*, ECCC TR26-260, 19 October 2026:
   *every $B_2$ circuit computing an affine disperser for dimension $d$ on
   $n$ bits has at least $3.5\,n-7d$ gates.* We reproduce the complete proof
   (definitions, invariants, the potential function, all operations and all
   ten cases) in `gate-elimination.md`, so the argument is checked from
   supplied evidence rather than from a citation.

2. **The fixed family (reproduced in `disperser-construction.md`).**
   The total family is the Ben-Sasson–Kopparty block construction
   (`BK` below), evaluated at every integer input length by zero-padding.
   Its affine-disperser dimension is $d(n)=O(n^{4/5})$. We state the
   construction and the supporting lemmas with proofs in
   `disperser-construction.md`. The only external computational ingredient is
   Shoup's deterministic irreducible-polynomial algorithm, stated precisely
   there; it is used only to make the family uniformly computable in
   polynomial time and does not enter the lower-bound accounting.

Everything else in this file is the assembly: membership in P/NP, the model
transfer, the parameter verification, the exact arithmetic, and the
all-length quantifiers.

---

## 1. The fixed total family

We use $\mathbb F_2=\{0,1\}$ and identify a bit string of length $k$ with an
element of the field $\mathbb F_{2^k}$ once a basis is fixed.

**Definition 1.1 (the family).**
Set $f_0(\varepsilon)=0$ on the empty input. For an integer $n\ge 1$ define

$$
k=k(n)=\text{the least prime strictly greater than }2n^{3/5},\qquad
r=r(n)=\left\lceil \frac{n}{k}\right\rceil,\qquad
t=t(n)=\left\lceil \sqrt{r}\right\rceil,\qquad
m=m(n)=r\,k .
$$

Fix once and for all a deterministic algorithm $\mathcal S$ with the property
of Theorem 10 in `disperser-construction.md` (Shoup). On input the prime $k$
(in unary), $\mathcal S$ outputs a monic irreducible polynomial
$q_k(Z)\in\mathbb F_2[Z]$ of degree $k$. Put $K=\mathbb F_2[Z]/(q_k)$ and use the polynomial basis
$1,Z,\dots,Z^{k-1}$; write $\varphi_i(x)\in K$ for the element whose
coefficient vector is the $i$-th consecutive block of $k$ bits of $x$, with
the first bit of the block the constant coefficient. Since $m\ge n$, append
$m-n$ zeros to $x$ before splitting, obtaining an injective
$\mathbb F_2$-linear map
$\Phi=(\varphi_1,\dots,\varphi_r):\{0,1\}^n\to K^r$. Let
$\pi:K\to\mathbb F_2$ be the constant-coefficient map, and set

$$
y_i=\varphi_i(x)^{\,2^{i+1}-1}\quad(1\le i\le r),\qquad
f_n(x)=\pi\!\left(\operatorname{Sym}_r^t(y_1,\dots,y_r)\right),
$$

where $\operatorname{Sym}_r^t$ is the sum of the products
$\prod_{i\in I}y_i$ over all $t$-element subsets $I\subseteq\{1,\dots,r\}$.
This is one fixed family: the prime $k$, the field representation, $\pi$,
and the exponents are all determined by $n$ through fixed deterministic
rules, and there is no per-length advice, no truth table, and no
nonconstructive choice. (The zero-padding is part of the definition of the
map $\Phi$; it is not a change of the input length, which remains $n$.)

The exponent sequence $2^{i+1}-1=1+2+\cdots+2^{i}$ is exactly the sequence
$1+p+\cdots+p^{i}$ of the BK construction at $p=2$. The order of the blocks is
irrelevant because $\operatorname{Sym}_r^t$ is symmetric; this is the
reindexing that identifies Definition 1.1 with BK Theorem 2.2 / Theorem 7.3.

---

## 2. The family is in P, hence in NP with empty witnesses

**Proposition 2.1.** There is one deterministic algorithm that, on input
$x\in\{0,1\}^n$, computes $f_n(x)$ in time polynomial in $n$.

*Proof.* (a) *Finding $k$.* Scan integers $k>2n^{3/5}$ upward and test
primality by trial division up to $\sqrt{k}$. By Bertrand's postulate there
is a prime in $(2n^{3/5},4n^{3/5}]$, so at most $2n^{3/5}$ candidates are
tested; each test costs $O(n^{3/10})$ divisions. This is polynomial in $n$.

(b) *Building $K$.* Run the fixed Shoup algorithm $\mathcal S$ on $k$; by
Theorem 10 of `disperser-construction.md` it runs in time polynomial in $k$,
hence in $n$, and outputs $q_k$. Field addition is bitwise XOR. Multiplication in $K$ is multiplication
of polynomials of degree $<k$ modulo $q_k$, costing $O(k^2)$ bit operations;
inversion, when needed, is obtained by the extended Euclidean algorithm in
$O(k^2)$.

(c) *The powers.* For each $i\le r$, the exponent $2^{i+1}-1$ has binary
expansion of length $i+1$, so $\varphi_i(x)^{2^{i+1}-1}$ is computed by
$i+1$ squarings and at most $i$ multiplications, i.e. $O(i)$ field
operations. Summed over $i\le r$ this is $O(r^2)$ field operations.

(d) *The symmetric polynomial.* Start with $E_0=1$ and $E_j=0$ for
$1\le j\le t$. Reading $y_1,\dots,y_r$ one at a time, update, in decreasing
$j$ from $\min(i,t)$ down to $1$,
$E_j\leftarrow E_j+y_iE_{j-1}$. Splitting the $t$-subsets of
$\{1,\dots,i\}$ according to whether they contain $i$ proves that after step
$i$, $E_j$ is the $j$-th elementary symmetric polynomial of
$y_1,\dots,y_i$; in particular the final $E_t$ equals
$\operatorname{Sym}_r^t(y_1,\dots,y_r)$. This is $O(rt)$ field operations.
(No division by $2$ is used, so characteristic two causes no difficulty.)

Since $k=O(n^{3/5})$, $r=O(n^{2/5})$, $t=O(n^{1/5})$, all the above are
polynomial in $n$. One deterministic algorithm therefore evaluates the whole
family at every length. $\square$

**Corollary 2.2 (NP membership and the verifier).**
Let $L=\{x: f_{|x|}(x)=1\}$. Define the polynomial-time verifier
$V(x,w)=1$ iff $w=\varepsilon$ (the empty string) and $f_{|x|}(x)=1$, with
witness bound $W(n)=0$. Then
$x\in L\iff \exists w,\ |w|\le 0,\ V(x,w)=1$: both directions are immediate
from Proposition 2.1. The witness is not an input to the target circuits.
The same language, evaluator, verifier and witness bound work at all lengths.

---

## 3. The circuit model and transfer

The challenge fixes Li–Yang's general $B_2$ model: deterministic acyclic
circuits with one output, two ordered input pins per gate, all sixteen binary
Boolean truth tables, unrestricted depth and fan-out, exact computation on
every input; only gates are counted, inputs/wires/output selection are free;
negations may be absorbed into the truth tables of descendant gates; the input
is $x$ alone, not an NP witness; repeated pins and shared subgraphs are
allowed.

The proof in `gate-elimination.md` is written for a circuit model with the
same gate count $|C|$: every gate of the input circuit is counted, while
input vertices may additionally be labelled by constants and gates that
depend on at most one input are first deleted (at no cost) during the
potential argument. This model is at least as permissive as the challenge's:
any challenge circuit is a circuit in the reproduction with exactly the same
gate count, because the reproduction never requires fewer gates for the input
circuit. The two facts used are:

* an input vertex labelled by a constant is free in the reproduction; the
  challenge's constant gates are ordinary counted gates there;
* a gate that depends on at most one input is deleted inside the argument,
  not assumed absent from the input circuit.

Formally, let $\mathcal C_{\text{ch}}$ be the set of challenge circuits for
$f_n$ and $\mathcal C_{\text{MS}}$ the set of circuits admitted by the
reproduction. Since every challenge circuit is also admitted there,
$C_{B_2}(f_n)=\min_{\mathcal C_{\text{ch}}}|\cdot|\ \ge\
\min_{\mathcal C_{\text{MS}}}|\cdot|$. The reproduced theorem lower-bounds the
right-hand side, hence the left-hand side. Nothing in the argument uses
per-length advice, truth-table oracles, nonconstructive existence, or a
restricted basis: the lower bound is for unrestricted $B_2$ circuits.

---

## 4. The affine-disperser property of $f_n$

**Definition 4.1.** A function $h:\{0,1\}^N\to\{0,1\}$ is an *affine disperser
for dimension $d$* if $h$ is nonconstant on every affine subspace of
$\mathbb F_2^N$ of dimension at least $d$ (equivalently, on every affine
subspace of dimension exactly $d$).

**Theorem 4.2 (BK block construction; Ben-Sasson–Kopparty 2010/2012).**
Let $N_{\!B}$ be prime, let $r_B,t_B$ be integers with $1\le t_B<r_B$,
$N_{\!B}\ge r_B(r_B+1)/t_B$, and put $M=r_BN_{\!B}$. Let $e_i=r_B+1-i$,
$s_i=\sum_{\ell=0}^{e_i}2^\ell$, and let
$h:K^{r_B}\to\mathbb F_2$ be
$h(x_1,\dots,x_{r_B})=\pi\bigl(\sum_{|I|=t_B}\prod_{i\in I}x_i^{s_i}\bigr)$.
Then $h$ is nonconstant on every affine subspace of
$\mathbb F_2^{M}\cong K^{r_B}$ of dimension $>D_0$, where
$D_0=N_{\!B}r_B/t_B+N_{\!B}t_B+r_B(r_B+1)$. Consequently, for any constant
$C>0$ with $D_0<CM^{4/5}$, $h$ is an affine disperser for dimension
$\lfloor CM^{4/5}\rfloor+1$.

The complete statement and proof are reproduced in
`disperser-construction.md` (key lemmas: the image-subspace polynomial, the
structural lemma on consecutive zero coefficients, the coefficient formula via
permanents, the independent-source theorem, Bourgain's decomposition, and the
final noncancellation argument). Note that $N_{\!B}$ is prime, so the only
proper subfield of $K=\mathbb F_{2^{N_{\!B}}}$ is $\mathbb F_2$; every subspace
of dimension at least $2$ is therefore not contained in an affine shift of a
proper subfield, which is the hypothesis used in the independent-source step.

**Lemma 4.3 (the family is an affine disperser at every large length).**
With $k,r,t,m$ as in Definition 1.1 and $d(n)=\lfloor 12m^{4/5}\rfloor+1$,
for all sufficiently large $n$, $f_n$ is an affine disperser for dimension
$d(n)$ on $N=n$ bits.

*Proof.* Fix $n$ large enough that $k\ge 3$, $r\ge 4$, $n/k\ge1$, and the
elementary bounds below hold. Set $N_{\!B}=k$, $r_B=r$, $t_B=t$, $M=m=rk$.

*Hypotheses of Theorem 4.2.*
(i) $k$ is prime by construction.
(ii) $t<r$: $t=\lceil\sqrt r\rceil\le r$ for $r\ge 1$, and $t<r$ for $r\ge3$.
(iii) $k\ge r(r+1)/t$: Since $t\ge\sqrt r$, $\sqrt r(r+1)\le 2r^{3/2}$ for
$r\ge1$, and $r=\lceil n/k\rceil\le n/k+1\le 2n/k$ (as $n/k\ge1$),
$$
\frac{r(r+1)}{t}\le\sqrt r\,(r+1)\le 2r^{3/2}\le 2\Bigl(\frac{2n}{k}\Bigr)^{3/2}
=4\sqrt2\Bigl(\frac{n}{k}\Bigr)^{3/2}
<4\sqrt2\Bigl(\frac{1}{2n^{2/5}}\Bigr)^{3/2}n^{3/5}
=2n^{3/5}<k,
$$
where we used $k>2n^{3/5}$, i.e. $n/k<\tfrac12 n^{2/5}$.

*Dimension accounting.* We verify $D_0<12m^{4/5}$. From $r\le 2n/k$ and
$t\ge\sqrt r$,
$$
\frac{kr}{t}\le k\sqrt r\le k\sqrt{2n/k}=\sqrt{2kn}.
$$
Also $kt\le k(\sqrt r+1)\le\sqrt{2kn}+k$, and
$$
r(r+1)\le \frac{2n}{k}\Bigl(\frac{2n}{k}+1\Bigr)\le
\frac{2n}{k}\cdot\frac{3n}{k}=\frac{6n^2}{k^2}.
$$
Bertrand's postulate gives $k<4n^{3/5}$, hence
$\sqrt{2kn}<\sqrt{8n^{8/5}}=2\sqrt2\,n^{4/5}$, and
$n^2/k^2<\frac14n^{4/5}$ because $k>2n^{3/5}$. Therefore, for all large $n$,
$$
D_0=\frac{kr}{t}+kt+r(r+1)
<2\sqrt2\,n^{4/5}+\bigl(2\sqrt2\,n^{4/5}+4n^{3/5}\bigr)
+\frac{6}{4}n^{4/5}
<8n^{4/5}\le 12m^{4/5},
$$
using $m\ge n$. Since the threshold in the definition of $d(n)$ is
$\lfloor12m^{4/5}\rfloor+1\ge 12m^{4/5}>D_0$, every affine subspace of
$\mathbb F_2^n$ of dimension at least $d(n)$ has dimension $>D_0$ after the
padding is applied, so Theorem 4.2 gives nonconstancy.

(The padding is an affine injection of $\mathbb F_2^n$ into
$\mathbb F_2^m=\mathbb F_2^{rk}\cong K^r$: appending a fixed zero string is a
linear map of rank $n$, so it carries an affine subspace of dimension $D$ to
an affine subspace of dimension $D$.)

Finally $d(n)\le n$ for large $n$, since $m\le n+k\le 2n$ and
$12(2n)^{4/5}<n$ once $n>(12\cdot2^{4/5})^{5}=12^{5}\cdot2^{4}$. $\square$

---

## 5. The gate-elimination lower bound

**Theorem 5.1 (Mishra–Santhanam, Theorem 4.1; proved in
`gate-elimination.md`).** Let $h:\{0,1\}^n\to\{0,1\}$ be an affine disperser
for dimension $d$. Then every $B_2$ circuit computing $h$ has at least
$$
\frac{7}{2}\,n-7d
$$
gates.

The proof is a potential-function gate-elimination argument. We summarise its
skeleton here so the assembly is self-contained; the complete case analysis is
in `gate-elimination.md`.

* State: an injective parameter map
  $\Phi(\vec u)=\vec b+L\vec u+\sum_{\{p,r\}\in P}\vec c_{p,r}\,pr$ from
  $\mathbb F_2^m$ to $\mathbb F_2^n$, a fair semi-circuit (an invertible
  affine system of equations with at most three nonconstant terms per row,
  followed by an acyclic $B_2$ part), and the set $P$ of disjoint parameter
  pairs. Invariant: $C(\vec u)=h(\Phi(\vec u))$ for all $\vec u$.
* Quantities: $s$ = circuit size (equations + acyclic gates), $m$ = number of
  parameters, $q=|P|$, $a=m-q$. Proposition 2.1 of the reproduction: the
  image of $\Phi$ contains an affine subspace of dimension $a$, so if $a\ge d$
  then $h\circ\Phi$ is nonconstant.
* Potential: $W=2s+3q+\sum_{x\text{ unpaired}}\operatorname{wt}(x)$, where
  $\operatorname{wt}(x)=0,6,7$ according as $\deg(x)=0,1,\ge2$.
* Lemma 9.1 of the reproduction (Lemma 4.2 of the source): while $a>d$, one
  polynomial-time iteration
  changes the state so that $m+s$ strictly decreases, $\Delta a\in\{0,-1,-2\}$,
  and $\Delta W\le 14\,\Delta a$. The iteration moves affine gates into the
  affine system, applies deletions $D1$–$D4$, affine reparameterisation $(AR)$,
  quadratic reparameterisation $(QR)$, and the ten-case analysis of the
  minimal AND gate.
* Start with $\Phi=\mathrm{id}$, $P=\varnothing$, empty affine system, and the
  input circuit as the acyclic part: $m=a=n$, $q=0$, $s=|C|$, so
  $W_0\le 2|C|+7n$.
* Iterate until $a_T\le d$. Summing $\Delta W\le14\Delta a$ and using
  $W_T\ge0$ gives $W_0\ge14(n-a_T)\ge14(n-d)$, hence
  $2|C|+7n\ge14(n-d)$ and $|C|\ge\frac72n-7d$. $\square$

(The reproduction also records why $7/2$ is the optimal leading coefficient
for the specific potential family used there; no reweighting improves it.)

---

## 6. Assembly and the declared guarantee

**Theorem 6.1 (main).** For the fixed family $(f_n)$ of Definition 1.1 there is
$N_0$ such that for every integer $n\ge N_0$,
$$
C_{B_2}(f_n)\ \ge\ \frac{7}{2}\,n-7\,d(n),
\qquad d(n)=\lfloor 12m(n)^{4/5}\rfloor+1 .
$$
Hence, for every $\varepsilon>0$ and all sufficiently large $n$,
$$
C_{B_2}(f_n)\ \ge\ \Bigl(\frac{7}{2}-\varepsilon\Bigr)n,
$$
and therefore
$$
\liminf_{n\to\infty}\frac{C_{B_2}(f_n)}{n}\ \ge\ \frac{7}{2}.
$$

*Proof.* Combine Lemma 4.3 with Theorem 5.1. For the second display, use
$m(n)\le n+k(n)\le 2n$ for large $n$ to get
$7d(n)\le 84(2n)^{4/5}+7=84\cdot2^{4/5}n^{4/5}+7=o(n)$, so
$C_{B_2}(f_n)/n\ge\frac72-o(1)$. $\square$

**Exact parameters.**
$c=7/2$, $p=1$, $q=r=0$ (the logarithms of the guarantee are absent), and
the finite projection is the same pair $(7/2,1)$. Both are reduced rationals:
$7/2$ and $1/1$. The declared guarantee $g(n)=\frac72 n$ is exactly the
leading term of Theorem 6.1, and the theorem gives
$\liminf C_{B_2}(f_n)/g(n)\ge1$ as required by the `polynomial_log` meaning.

**Comparison with the seed.** The seed declares $(31/10,1)$ for the same
BK/Shoup family using Li–Yang's Theorem 1.1. We do not dispute that
literature statement; we replace its use by the strictly stronger
Mishra–Santhanam bound $3.5n-7d$ (their Theorem 4.1) and keep the same kind of
family. Since $7/2=3.5>3.1=31/10$ and both are guarantees of the form
$cn-o(n)$ for fixed eligible families, the full-expression comparator
`expressions.compare` returns `stronger` for our guarantee against the seed's
guarantee (polynomial exponent equal, log exponents equal, coefficient
larger).

---

## 7. Quantifiers, all lengths, and conventions

* The family, the field representation, the prime $k(n)$, the padding and the
  evaluator are functions of the *actual* input length $n=|x|$ only. They do
  not depend on $\varepsilon$, on the circuit, on a witness, or on a chosen
  subsequence of lengths.
* The lower bound is uniform over every integer $n\ge N_0$; there is no sparse
  or limsup statement. The only loss is the additive $7d(n)=O(n^{4/5})$, which
  divided by $n$ tends to $0$. In particular the remainder is $o(n)$, which is
  exactly what the $c\,n^p-r(n)$ form of the rules permits.
* No practical threshold is asserted: $N_0$ is some absolute constant whose
  existence follows from the finitely many asymptotic inequalities in
  Lemma 4.3.
* The exponent $p=1$ is the power of the actual input length, not of a block
  length, a field size, or a witness length. The padding to $m=rk$ bits is
  internal to the definition of $f_n$; the lower bound is measured against the
  $n$ genuine input bits, and Lemma 4.3 accounts for the padding.
* The full guarantee is a liminf over all integer $n$; the projection is the
  finite pair $(7/2,1)$. Neither an endpoint nor a superpolynomial rate is
  claimed. We do not claim any exponent $p>1$; that would be a different and
  much stronger theorem.
* We do not claim that $3.5$ is optimal for unrestricted $B_2$ circuits; the
  reproduction notes only that no reweighting of the specific potential used
  there improves $7/2$.

---

## 8. Provenance and evidence map

| Evidence | Content |
| --- | --- |
| `bound.json` | Canonical declaration: version, convention, family, theorem, guarantee, projection. |
| `proof.md` (this file) | Result, family, P/NP membership, model transfer, disperser instantiation, assembly, quantifiers, provenance. |
| `gate-elimination.md` | Self-contained reproduction of ECCC TR26-260 §§2–4: circuits, parameters, fair semi-circuits, potential, operations, deletions, Lemma 5.2 (source Lemma 2.3), the ten cases, Lemma 9.1 (source Lemma 4.2), proof of Theorem 4.1. |
| `disperser-construction.md` | Ben-Sasson–Kopparty block affine disperser: subspace polynomials, structural lemma, permanent coefficient formula, independent-source theorem, Bourgain decomposition, Theorem 7.3, parameter verification; Shoup's deterministic irreducible-polynomial theorem. |

**Primary sources.**

* O. S. Mishra and R. Santhanam, *An Explicit 3.5n−o(n) Size Lower Bound
  against General Boolean Circuits*, ECCC TR26-260, 19 October 2026.
  Theorem 1.1, Theorem 4.1 (pp. 27–28 in the ECCC PDF). The gate-elimination
  proof is reproduced in `gate-elimination.md`.
* E. Ben-Sasson and S. Kopparty, *Affine dispersers from subspace polynomials*,
  preprint 12 March 2010; SIAM J. Comput. 41(4):880–914, 2012. Theorem 2.2 and
  Sections 3, 6, 7. The construction, lemmas and block theorem are reproduced
  in `disperser-construction.md`.
* V. Shoup, *New algorithms for finding irreducible polynomials over finite
  fields*, Math. Comp. 54:435–447, 1990; extended abstract FOCS 1988.
  Deterministic polynomial-time construction for fixed prime characteristic.
  Stated as Theorem 10 in `disperser-construction.md`.
* J. Li and T. Yang, *3.1n−o(n) circuit lower bounds for explicit functions*,
  STOC 2022. Only as the seed's comparator; not used in the proof.
* ECCC TR21-023 (Li–Yang source conventions for the $B_2$ model) as cited by
  `rules.md`.

**What is checked from supplied evidence.** The complete gate-elimination
argument, the complete disperser construction and block theorem (including
the finite-field lemmas), the padding and parameter verification, the P/NP
argument, and all exact arithmetic are in this submission. Shoup's theorem is
stated with its exact scope and is used only for uniform polynomial-time
computability of the field representation; it is not used in any lower-bound
inequality. No candidate-created dossier, alleged approval, or comparison hint
is used.
