# Certified finite path sum for the cubic vertex Gaussian wave

This note explains the numerical component implemented in `certify.py`. The integral bounds are analytic; floating point agreement is not used as a certificate. The output must be combined with the separately proved graph comparison and layout conversion.

## 1. The path process and the finite target

Let `(X_i)` be the centered stationary Gaussian process with

\[
\mathbb E X_iX_{i+d}=w_d=(1+d/3)2^{-d/2}.
\]

The covariance recurrence and the initial covariances show that this is exactly the stationary autoregression

\[
X_i=\sqrt2X_{i-1}-\tfrac12X_{i-2}+\epsilon_i/\sqrt{12},
\]

where the innovations are independent standard normal variables and independent of the past. Indeed `w_1=2 sqrt(2)/3`, `w_2=5/6`, the recurrence holds for every `d>=2`, and the innovation variance is `1/12`. Gaussian uncorrelatedness gives independence.

Writing `Y_i=(-1)^i X_i`, its transition density for positive states is

\[
k(z\mid x,y)=\sqrt{6/\pi}\exp[-6(z+\sqrt2y+x/2)^2].
\]

The initial pair density is

\[
f_2(x,y)=\frac3{2\pi}\exp[-\tfrac92(x^2+y^2)-6\sqrt2xy].
\]

The probability that a fixed `j`-edge path alternates is

\[
p_j=2\int_{(0,\infty)^{j+1}}f_2(x_0,x_1)
 \prod_{i=2}^j k(x_i\mid x_{i-2},x_{i-1})\,dx.
\]

For each fixed odd integer `J`, the graph proof uses the **finite** sum

\[
\gamma_J=\sum_{j=1}^J(-1)^{j+1}3\,2^{j-2}p_j.
\]

The two-step completion estimate below also proves convergence and encloses the infinite sum between consecutive even and odd truncations.

## 2. Truncation with positive path completions

Throughout, `n=j+1>=2`, `L>0`, and `N>=1`. Put

\[
c_2=\frac3{2\pi},\qquad c_1=\sqrt{\frac6\pi}.
\]

The positive-state pair and transition densities are

\[
f_2(x,y)=c_2\exp[-\tfrac92(x^2+y^2)-6\sqrt2xy],
\]
\[
k(z\mid x,y)=c_1\exp[-6(z+\sqrt2y+x/2)^2].
\]

The exact joint density is

\[
f_n(x_0,\ldots,x_{n-1})
 =f_2(x_0,x_1)\prod_{r=2}^{n-1}k(x_r\mid x_{r-2},x_{r-1}).
\]

### 2.1 Positive completions from an arbitrary anchor pair

The transition is reversible in the precise order-two sense

\[
f_2(x,y)k(z\mid x,y)=f_2(y,z)k(x\mid z,y).             \tag{A}
\]

Indeed both sides have prefactor `c_2*c_1` and exponent

\[
-6x^2-\frac{33}{2}y^2-6z^2
-12\sqrt2xy-12\sqrt2yz-6xz,
\]

which is symmetric under interchange of `x` and `z`. This is a direct
identity of the displayed densities, so no unproved assertion about a
reversed Markov process is needed.

For any adjacent anchor positions `a,a+1`, where `0<=a<=n-2`, repeated use
of (A) gives the exact factorization

\[
\begin{split}
f_n(x_0,\ldots,x_{n-1})
={}&f_2(x_a,x_{a+1})
 \prod_{r=0}^{a-1}k(x_r\mid x_{r+2},x_{r+1})\\
&\quad\cdot
 \prod_{r=a+2}^{n-1}k(x_r\mid x_{r-2},x_{r-1}).       \tag{B}
\end{split}
\]

An empty product is one. To verify the telescoping explicitly, move the
initial pair from `(x_0,x_1)` to `(x_1,x_2)` by (A), leaving the reverse
factor `k(x_0|x_2,x_1)`; then move it to `(x_2,x_3)`, leaving
`k(x_1|x_3,x_2)`, and continue until the requested anchor is reached.
All factors at indices greater than `a+1` stay in their original order.

For fixed `x,y>=0` and `z>=0`,

\[
k(z\mid x,y)\le c_1e^{-6z^2},\qquad
\int_0^\infty k(z\mid x,y)\,dz\le\frac12.           \tag{C}
\]

The first inequality follows because `z+sqrt(2)*y+x/2>=z`, and the second
uses the elementary Gaussian half-integral. With positive anchor values
fixed, integrate the leftmost coordinate `x_0`, then `x_1`, and so on
through the left arm of (B). Each integral removes one reverse factor and
costs at most `1/2`. Similarly integrate the rightmost coordinate
`x_(n-1)`, then `x_(n-2)`, and so on through the right arm. All remaining
coordinates are still positive when (C) is applied. Nonnegativity permits
these integrations by Tonelli's theorem. There are exactly `n-2` integrated
coordinates, so

\[
\int_{(0,\infty)^{n-2}} f_n\,dx_{\text{outside anchor}}
 \le 2^{-(n-2)}f_2(x_a,x_{a+1})                     \tag{D}
\]

on the positive quadrant of the anchor pair. This is a bound on the actual
conditional completion mass. It does not assume that the completion events
or the anchor coordinates are independent.

### 2.2 Omitted-tail bound

For `x>=L` and `y>=0`, retain the negative cross term in the pair density:

\[
\begin{split}
\int_L^\infty\int_0^\infty f_2(x,y)\,dy\,dx
&\le c_2\int_L^\infty e^{-9x^2/2}
       \left(\int_0^\infty e^{-6\sqrt2xy}\,dy\right)dx\\
&=\frac{c_2}{6\sqrt2}
       \int_L^\infty \frac{e^{-9x^2/2}}x\,dx\\
&\le\frac{c_2}{6\sqrt2 L}
       \int_L^\infty e^{-9x^2/2}\,dx\\
&\le\frac{c_2e^{-9L^2/2}}{54\sqrt2 L^2}.             \tag{E}
\end{split}
\]

The first inequality drops only the nonpositive exponent term `-9y^2/2`.
For the last inequality, `x/L>=1` on the integration domain, so

\[
\int_L^\infty e^{-9x^2/2}\,dx
 \le\frac1L\int_L^\infty x e^{-9x^2/2}\,dx
 =\frac{e^{-9L^2/2}}{9L}.
\]

For each possible outlier coordinate `x_i>L`, choose an adjacent pair
containing it: positions `(i,i+1)` if `i<n-1`, and `(n-2,n-1)` if `i=n-1`.
The pair density is symmetric, so (E) applies whether the outlier is its
first or second coordinate. Bound the other positive coordinates using (D),
and take the union bound over the `n` choices. The mass omitted on replacing
the positive orthant by `[0,L]^n` is at most

\[
\boxed{
T_n^{\mathrm{new}}
=\frac{n\,c_2\,e^{-9L^2/2}}
 {2^{n-2}\,54\sqrt2\,L^2}.
}                                                        \tag{F}
\]

The formula is valid also at `n=2`, when there are no completion coordinates.
It bounds one positive orthant; the alternating probability still multiplies
the quadrature and all its errors by two.

The factor `2^{-(n-2)}` is useful because the graph coefficient of `p_j`
grows as `2^{j-2}`. Those factors cancel in this tail contribution. The
result is nevertheless a finite union-bound computation for each `j`, not
an inference about an infinite alternating series.

### Positive completion in pairs of steps

For nonnegative initial states $x,y$, the probability that both next states are positive is at most the probability with $x=y=0$. Both positive transition densities in the two-step integral are pointwise decreasing in the starting coordinates. With independent standard normals $U,V$, the maximal event is $U>0,\ V>\sqrt2\,U$. Its probability is $\arctan(1/\sqrt2)/(2\pi)<1/10$. For the strict inequality, $\cos(2\arctan(1/\sqrt2))=1/3>(\sqrt5-1)/4=\cos(2\pi/5)$. A single positive extension has conditional probability at most $1/2$.

Put $H(k)=10^{-\lfloor k/2\rfloor}2^{-(k\bmod2)}$. Conditioning in pairs bounds the chance of $k$ positive steps by $H(k)$ uniformly over nonnegative starting pairs. Formula (B) factors the two arms with the anchor values held fixed. For a possible outlier at $i=0,\dots,n-1$, use the pair starting at $a=\min(i,n-2)$. Bound its one-coordinate tail by (E) and the two remaining arms by $H(a)H(n-2-a)$. The program uses the resulting single-orthant tail bound

    T_n = [c2 exp(-9 L²/2)/(54 sqrt(2) L²)]
          × sum over i = 0,...,n-1 of H(a) H(n-2-a).

For $a_j=3\,2^{j-2}p_j$, the one-step bound gives $a_{j+1}\le a_j$, and two steps give $a_{j+2}\le(2/5)a_j$. Thus $\sum_j a_j$ converges. The alternating sum $\Gamma$ lies between every even partial sum and the next odd partial sum, and $\Gamma=\inf_{J\ {\rm odd}}\Gamma_J$. Section 7 of the graph proof obtains a saving below this limit by changing the layout.

## 3. Gaussian moment bound for tensor quadrature

### 3.1. Conditional one-coordinate form and amplitude

Let

\[
f_n(x)=f_n(0)\exp(-x^TQ_nx/2),\qquad
f_n(0)=c_2c_1^{n-2},
\]

where `Q_n` is positive definite, `c_2=3/(2*pi)` and
`c_1=sqrt(6/pi)`. Fix every coordinate except `x_i=t`, and write

\[
x^TQ_nx=q_it^2+2bt+c,
\qquad q_i=(Q_n)_{ii}>0.
\]

Completing the square gives

\[
f_n(t,x_{-i})=a\exp[-q_i(t-\mu)^2/2],
\qquad \mu=-b/q_i,
\qquad a=f_n(0)\exp[-(c-b^2/q_i)/2].
\]

The residual `c-b^2/q_i` is nonnegative: it is the minimum of the
nonnegative quadratic form `x^TQ_nx` over real `t`, with the other
coordinates fixed. Hence

\[
0<a\le f_n(0).                                           \tag{1}
\]

This proves the needed amplitude bound directly; no unstated conditional
Gaussian or Schur-complement assertion is being used.

### 3.2. Gaussian characteristic function with normalization

For any `q>0`, let

\[
\phi_q(z)=\frac1{\sqrt{2\pi q}}e^{-z^2/(2q)},\qquad
F_q(s)=\int_{\mathbb R}e^{isz}\phi_q(z)\,dz.
\]

The Gaussian integral normalizes `phi_q` to one. For real `s`, dominated
differentiation is valid because every polynomial times `phi_q` is
integrable. Since `z*phi_q(z)=-q*phi_q'(z)`, integration by parts gives

\[
\begin{split}
F_q'(s)
 &=i\int_{\mathbb R} z e^{isz}\phi_q(z)\,dz\\
 &=-iq\int_{\mathbb R}e^{isz}\phi_q'(z)\,dz\\
 &=-qsF_q(s).
\end{split}
\]

The boundary terms vanish by Gaussian decay. With `F_q(0)=1`, the solution
of this elementary differential equation is

\[
F_q(s)=e^{-qs^2/2}.                                      \tag{2}
\]

No Fourier inversion theorem is required.

### 3.3. Exact even moments and derivative bound

Let `M_(2k)=integral_R z^(2k)*phi_q(z) dz`. A second integration by
parts, again with vanishing boundary terms, gives

\[
M_{2k}=q(2k-1)M_{2k-2},\qquad M_0=1.
\]

Thus

\[
M_{2N}=q^N(2N-1)!!=q^N\frac{(2N)!}{2^N N!}.              \tag{3}
\]

Differentiating the integral for `F_q` exactly `2N` times is justified by
the integrable majorant `|z|^(2N)*phi_q(z)`. Using (2), (3), and
`|exp(isz)|=1` on the real line yields

\[
\left|\frac{d^{2N}}{ds^{2N}}e^{-qs^2/2}\right|
 \le M_{2N}
 =q^N\frac{(2N)!}{2^N N!}.                              \tag{4}
\]

In fact equality in absolute value holds at `s=0`, but only the inequality
is needed. Applying (4) to the one-coordinate form and using (1) gives the
uniform real derivative bound

\[
\boxed{
|\partial_i^{2N}f_n(x)|
 \le f_n(0)q_i^N\frac{(2N)!}{2^N N!}.
}                                                        \tag{5}
\]

It holds for every real vector `x`, and in particular throughout the cube
used for quadrature.

### 3.4. Exact precision diagonals

The joint-density exponent is

\[
-\frac92(x_0^2+x_1^2)-6\sqrt2x_0x_1
-6\sum_{r=2}^{n-1}(x_r+\sqrt2x_{r-1}+x_{r-2}/2)^2.
\]

In the precision matrix, each transition square contributes `12` to its
current-coordinate diagonal, `24` to its previous-coordinate diagonal, and
`3` to its second-previous-coordinate diagonal. The initial pair contributes
`9` to each of its two diagonals. Therefore the exact list is

\[
(q_i)_{i=0}^{n-1}=\begin{cases}
(9,9),&n=2,\\
(12,33,12),&n=3,\\
(12,36,\underbrace{39,\ldots,39}_{n-4},36,12),&n\ge4.
\end{cases}                                              \tag{6}
\]

Consequently the exact integer sum needed below is

\[
S_{n,N}:=\sum_iq_i^N=\begin{cases}
2\cdot9^N,&n=2,\\
2\cdot12^N+33^N,&n=3,\\
2\cdot12^N+2\cdot36^N+(n-4)39^N,&n\ge4.
\end{cases}                                              \tag{7}
\]

There is no precision-matrix approximation in (6)--(7).

### 3.5. Tensor Gauss--Legendre error

The usual `N`-node Gauss--Legendre error on `[0,L]` has bound

\[
\frac{L^{2N+1}(N!)^4}{(2N+1)((2N)!)^3}
\sup|h^{(2N)}|.
\]

Apply (5) for coordinate `i`. The other `n-1` integration or positive
Gauss-quadrature operators each have supremum norm `L`. Telescoping the
`n`-fold tensor products and summing the coordinate bounds therefore gives

\[
\boxed{
E^{\mathrm{moment}}_{n,N}
=\frac{f_n(0)L^{n+2N}(N!)^3}
 {(2N+1)2^N((2N)!)^2}\,S_{n,N}.
}                                                        \tag{8}
\]

The cancellation of factorials is explicit: (5) contributes `(2N)!/N!`,
which reduces `(N!)^4/((2N)!)^3` to `(N!)^3/((2N)!)^2`.
The factor `2^N` remains in the denominator. The powers of `L` are
`L^(2N+1)` from the one-coordinate error and `L^(n-1)` from the remaining
operators, giving `L^(n+2N)`.

Write E^(moment)_(n,N) for the Gaussian-moment bound in (8). Its factors are exact integers or positive directed intervals; the density power is evaluated by repeated interval multiplication. The one-node Gauss normalization gives the midpoint bound L^3 sup|h''|/24.

### A product envelope for the tensor error

Write the density as $c_2c_1^{n-2}$ times the pair and transition exponentials. Fix a coordinate $i$ and allow it to vary on the real line, holding every other coordinate nonnegative. Pull out the factors independent of $x_i$. All quadratics in the remaining negative exponents are nonnegative at every real $x_i$, including the positive definite initial pair quadratic. Completing the square in their product yields a Gaussian with precision $q_i$ and amplitude at most one. The derivative estimate (5) applies to it.

For every $j\in\{2,\ldots,n-3\}$ at distance more than two from $i$ (zero-based indices), the full coefficient $39/2$ at $x_j^2$ survives in the pulled-out factors. All cross terms in their negative log density are nonnegative on the positive orthant. The derivative bound in (5) can therefore be multiplied by $\prod_j \exp(-39x_j^2/2)$ over these indices. There are at least $\max(0,n-9)$ of them. The corresponding integral has mass at most $\sqrt{\pi/78}<1/4$. For the positive quadrature rule the program encloses $\sum_k w_k\exp(-39x_k^2/2)$ and checks the upper endpoint is below $1/4$. The other integration or quadrature operators have mass $L$. Telescoping as in (8) gives the single-orthant error $E_{n,N}^{\rm moment}/(4L)^{\max(0,n-9)}$. This reasoning applies to every mixed tensor product in the telescoping sum.

## 4. Certified quadrature nodes, weights and arithmetic

Pi, root brackets, quadrature weights, density factors, analytic error bounds and the final alternating sum use Python's standard-library `decimal` with 100 decimal digits. Every Decimal lower endpoint uses `ROUND_FLOOR` and every upper endpoint uses `ROUND_CEILING`. Repeated positive tensor contractions use exact integers on the fixed-point grid proved in Section 5, with outward integer floor and ceiling. Square root and exponential are correctly rounded in Decimal; moving their rounded values one representable number outward encloses the exact transcendental value. Pi is enclosed by Machin's identity and 100-term alternating arctangent series, including an explicit next-term remainder.

The Legendre roots are enclosed by exact rational sign checks. At $x=a/b$ with $b>0$, set $t_k=k!b^kP_k(x)$. The Legendre recurrence gives $t_0=1,\ t_1=a$ and

    t_k = (2k-1) a t_(k-1) - (k-1)^2 b² t_(k-2).

Integer signs at the endpoints of each bracket are opposite. The $N$ disjoint brackets lie in $(-1,1)$ and isolate all roots. Evaluate the derivative at the rational center with that same integer recurrence and $(x^2-1)P_N'(x)=N(xP_N(x)-P_{N-1}(x))$. Directed division encloses the rational result. The error out to either end of the bracket is at most $N^4$ times its distance from the center. Indeed $\|P_n\|_\infty\le1$ on $[-1,1]$ by the standard Legendre integral representation $P_n(\cos t)=\pi^{-1}\int_0^\pi(\cos t+i\sin t\cos u)^n\,du$. Also $P_{k+1}'-P_{k-1}'=(2k+1)P_k$ by the Legendre generating function. Telescoping gives $\|P_n'\|_\infty\le n^2$, and differentiating that sum gives $\|P_n''\|_\infty\le n^4$. Evaluate the weight on the entire bracket using this derivative enclosure and the following exact formula.

\[
w=\frac2{(1-x^2)(P_N'(x))^2}
\]

Roots and weights are then mapped to `[0,L]`.

The tensor quadrature is evaluated by the exact Markov factorization, without changing its mathematical value. If `M[i,j]` already includes the quadrature weights in both live coordinates, initialize

\[
M[i,j]=w_iw_j f_2(x_i,x_j)
\]

and update

\[
M_{new}[j,k]=w_k\sum_i M[i,j]k(x_k\mid x_i,x_j).
\]

Every term is nonnegative. The exact fixed-point enclosure proved in Section 5 therefore gives lower and upper bounds on each contraction. Summing the final matrix encloses the complete tensor quadrature. Subtracting/adding `E_{n,N}`, adding the one-sided tail `T_n`, and multiplying by two encloses `p_j`. The finite signed sum is then evaluated with interval arithmetic, including the correct endpoint reversal for its negative terms.

For speed, the transition density is evaluated by the exact factorization

\[
k(z\mid x,y)=\sqrt{6/\pi}e^{-3x^2/2}e^{-12y^2}e^{-6z^2}
 e^{-6\sqrt2xy}e^{-6xz}e^{-12\sqrt2yz}.
\]

All exponential factors depend on at most two coordinates, so only O(N^2) exponential evaluations are required. The program precomputes the positive interval factors and combines each triple with directed multiplications. This is an algebraic identity, with every multiplication error still enclosed; it does not change the quadrature rule or its analytic error bound.

The implementation does not assume correctly rounded platform binary `libm` transcendental functions: the initial binary cosine guess is irrelevant after the root-bracket sign checks succeed.

## 5. Exact fixed-point contraction with quadratic storage

For the exact Gauss nodes \(x_i\) and weights \(w_i\), let

\[
\begin{aligned}
A_{ij}&=\exp(-3x_i^2/2-6\sqrt2x_ix_j),\\
B_{ik}&=\exp(-6x_ix_k),\\
C_{jk}&=\sqrt{6/\pi}\,w_k
 \exp(-12x_j^2-6x_k^2-12\sqrt2x_jx_k).
\end{aligned}
\]

These factors are nonnegative and satisfy the exact identity

\[
w_k k(x_k\mid x_i,x_j)=A_{ij}B_{ik}C_{jk}.
\]

Consequently

\[
T_{ij}=M_{ij}A_{ij},\qquad
M'_{jk}=C_{jk}\sum_i T_{ij}B_{ik}.
\]

Only two-index arrays are needed. In particular, the same Gauss tensor rule
can be evaluated with \(O(N^2)\) storage and \(O(JN^3)\) integer arithmetic.
Moving \(C_{jk}\) outside the sum is an algebraic identity and does not
approximate the transition density.

The implementation fixes the exact integer scale \(S=2^{224}\). An integer
pair \((u,v)\) represents the interval \([u/S,v/S]\). To convert a
nonnegative Decimal interval \([a,b]\), the code uses
`Decimal.as_integer_ratio()` and exact integer division to obtain

\[
u=\lfloor Sa\rfloor,\qquad v=\lceil Sb\rceil.
\]

No rounded Decimal multiplication or conversion to a binary floating point
number is used in this conversion. Since \(u/S\le a\le b\le v/S\), it is
outward. The initial state and the three factor arrays are converted in this
way from their already established Decimal enclosures.

Suppose \((m^-_{ij},m^+_{ij})\), \((a^-_{ij},a^+_{ij})\),
\((b^-_{ik},b^+_{ik})\), and \((c^-_{jk},c^+_{jk})\) enclose the exact
nonnegative quantities \(M_{ij},A_{ij},B_{ik},C_{jk}\), respectively,
with the common denominator \(S\). First compute

\[
t^-_{ij}=\left\lfloor\frac{m^-_{ij}a^-_{ij}}S\right\rfloor,
\qquad
t^+_{ij}=\left\lceil\frac{m^+_{ij}a^+_{ij}}S\right\rceil.
\]

Then \(t^-_{ij}/S\le T_{ij}\le t^+_{ij}/S\). Next compute the exact integer
sums

\[
d^-_{jk}=\sum_i t^-_{ij}b^-_{ik},\qquad
d^+_{jk}=\sum_i t^+_{ij}b^+_{ik}.
\]

Every product is nonnegative, so

\[
\frac{d^-_{jk}}{S^2}
\le\sum_i T_{ij}B_{ik}
\le\frac{d^+_{jk}}{S^2}.
\]

The next state is

\[
\begin{aligned}
(m')^-_{jk}
 &=\left\lfloor\frac{c^-_{jk}d^-_{jk}}{S^2}\right\rfloor,\\
(m')^+_{jk}
 &=\left\lceil\frac{c^+_{jk}d^+_{jk}}{S^2}\right\rceil.
\end{aligned}
\]

Dividing these last integers by \(S\) encloses \(M'_{jk}\): before the last
floor or ceiling the exact bounding products have denominator \(S^3\),
and the displayed division by \(S^2\) restores the represented denominator
\(S\). This proves by induction that every state array encloses the exact
tensor state. All the sums and products between the floor and ceiling
operations use Python's exact arbitrary-size integers. For a nonnegative
integer numerator \(u\) and positive integer denominator \(v\), the code
implements ceiling as \((u+v-1)//v\); floor is \(u//v\).

At path length \(j\), summing the lower and upper state integers gives exact
integer bounds \(s^-\) and \(s^+\). The interval for the quadrature is
\([s^-/S,s^+/S]\). Its endpoints are converted back to Decimal only by
explicit downward and upward division. Decimal construction from an integer
is exact. The subsequent analytic error and finite alternating-sum operations
therefore retain the same enclosure guarantee as Section 4.

This is a different directed-rounding implementation of the same finite
quadrature. Every quantization error is carried in its interval endpoints;
none is ignored or estimated from numerical agreement with another run.
The analytic truncation and Gauss-quadrature bounds proved in Sections 2--3
apply to this same finite tensor rule.

### Sparse positive contractions

The implementation checks that every scaled upper endpoint in B is at most S. For a scaled product M*A with rounded upper endpoint at most 8, it omits the nonnegative lower contribution and charges the whole upper endpoint times S to the dot product. If the sum of all scaled upper endpoints times the scaled upper endpoint of C is at most S, the exact output lies in [0,1/S]. This permits a short circuit. All other products use the floor and ceiling divisions above. Symmetric reuse of exponential factors is exact.

## 6. Parameters and target

`python -I -B certify.py` checks 75 Gauss nodes on [0,3] at each of 27 path lengths, using 100-digit directed Decimal enclosures and exact 224-bit fixed-point contractions. It verifies

    0.130647566419 < Gamma_26 <= Gamma <= Gamma_27 < 0.130647566425.

It also checks the preparatory pivotal-event floor 18079200775 at scale 10^12 from cutwidth.md Section 7.4, the all-center lower sum, the two covariance integrals, the central and outer rational margins, and the final score arithmetic. The public scripts and documentation of Galarnon42, DarkWingMcQuack and subsequent contributors supply the preceding path calculation. Each finite premise used here is reproduced by this executable, with the corresponding analytic error bounds above and below.

## 7. Conditional means and strip inequalities

Section 7.5 of cutwidth.md contains the good-patch integral and its numerical table. For convenience this section records how to obtain the covariances between its named scores and the three conditioning coordinates. For a pair of linear forms sum the products of coefficients times w_d. For an upper error add η times the sum of its positive products times k_d; for a lower error add η times the negative sum. Take k_d=d²(√2)^(d−1), k_0=0. The root regression coefficient for the residual of a linear form is its covariance with the root O. Since O has variance one the conditional covariance is obtained by subtracting the product of these two root covariances. The conditional centered means and variances in Section 7.5 are computed with these operations.

## 8. Integrated credit from adjacent and distance-two pairs

This note supplements the previously proved six-contrast tilted private-edge method. It yields the effective per-edge coefficient

    K >= 0.016117814324835210

per unit of total adjacent covariance excess, before the chosen finite-profile approximation allowance. The old tilted coefficient was 0.000672905379119700. This note analyzes the heavy-edge term only; a uniform improvement of the cubic cutwidth constant also requires the companion good-center and noncentral-threshold bounds.

### Preliminaries

Use the finite high-girth cover and base Hellinger fields of the cutwidth proof. The normalized critical tree wave X has w_d=(1+d/3)2^(-d/2) and lambda=2sqrt(2). Let Y be an independent limiting pullback of base fields and Z=sqrt(1-s)X+sqrt(s)Y. P denotes the two coordinates being differentiated and we condition on Z_P=0. In the standard Gaussian Price interpolation a pair contributes its covariance excess times the zero-pair density times the conditional mixed difference of the crossing-matching number F. All these mixed differences are nonpositive. Raw finite covariance differences are nonnegative by fiber Cauchy-Schwarz.

The preceding tilted proof gives two elementary bounds on Y. For each vertex t its wave defect R_t=sum_(q~t)Y_q-lambda Y_t has standard deviation at most c=3-lambda, and each centered contrast sqrt(3)(Y_a-Y_b) of distinct siblings has variance at most one. For completeness, the full normalized alpha-profile f on the universal tree satisfies sum_(q~t)f_q=(alpha^(-1)+2alpha)f_t except for a correction at t of norm tending to zero as alpha increases to 1/sqrt(2); also sum f_q^2<=(alpha^(-2)+2alpha^2) f_t^2 pointwise. Triangle and Cauchy-Schwarz on each fiber bound the sum of the pushed rows between (lambda-o(1))b_t minus a vector of norm o(1) and (3+o(1))b_t. Thus the defect norm is at most 3-lambda in the limit. For siblings, base covariance is at least w_2=5/6. Both are variance upper bounds after conditioning at zero. The covariances in any fixed finite patch form a compact family. The high-girth finite cover converges locally to the critical wave, after choosing alpha, finite truncation and positive tie noise successively before graph size as in the submitted proof.

We use the following Gaussian facts. If W is a nonnegative even function and g standard Gaussian, translating the density and averaging the opposite linear exponents gives E W(g+h)>=exp(-||h||^2/2) E W(g). Given a fixed linear contrast vector H of the centered conditional Y, tilting by exp(-r||H||^2/2), r=s/(1-s), gives total mass det(I+r Cov(H))^(-1/2) and a centered Gaussian law with every variance decreased. This follows by writing Y as a linear map of a standard normal vector and completing the square. Finally, for jointly centered Gaussian scalar forms V_i, the probability of the intersection of symmetric strips |V_i|<a_i is at least the product of the marginal strip probabilities (Sidak's Gaussian rectangle inequality). One proof conditions on the last coordinate. The probability of the earlier symmetric convex box under the resulting Gaussian translation is an even log-concave function of the conditioning coordinate by Prékopa's theorem. It therefore decreases with its absolute value and has nonnegative covariance with the indicator of the final symmetric interval. Induct, using a small independent-noise limit for singular covariances. Write H(t)=Pr(|N(0,1)|<t) for t>=0.

### Integrating the individual margins on adjacent pairs

At an adjacent pair the submitted tilted proof uses two sign-reflected private-edge events, whose indicators sum to at most -Delta F, and sixteen choices of the labeled private branches. For fixed labels the two-orientation witness is even in six independent critical contrasts g. The actual field, divided by sqrt(1-s), is the ideal recursion in g+sqrt(r)H_Y plus sqrt(r) times sibling-shared wave-defect errors. At a selected endpoint the three errors are R_i/2, (lambda/4)R_i+R_a/2, and (3/4)R_i+(lambda/4)R_a+R_z/2. The other side has the same formulas. Hence their standard deviations are bounded by tau=(.086,.208,.336) in the tilted law as well. Each H_Y contrast has variance at most one, so the tilt mass is at least (1-s)^3 by Hadamard's inequality.

Instead of fixed margins choose the actual absolute margins of the realized ideal wave. Conditional on Y, require sufficient symmetric margins to absorb its realized six errors. This condition is even in g, so the same translation bound applies. Now integrate over g and then the tilted Y. For each ideal witness sample the six actual margin thresholds are fixed; Sidak's inequality gives the product of their individual Gaussian strip probabilities. The adaptively labeled events are disjoint per orientation: if the specified signs persist, each specified-sign outer neighbor and opposite child is uniquely determined. Both orientations may be counted because their simultaneous occurrence forces -Delta F=2.

On one side x is the absolute outer-neighbor amplitude, k is the opposite-child amplitude, and z the absolute difference of its two children. Their three same-pair minimum margins are x, k, and (lambda k+x-z)/2. Let f_A be the density of |N(0,1/12)| and f_D that of |N(0,1/3)|. Given x the opposite child has amplitude k=(|D|-lambda x)/2 when positive, and its sibling has amplitude lambda x+k. Then one oriented surviving witness has probability at least (1-s)^3 J(s)^2, where

    J(s)=integral_(x,k>=0) f_A(x) 2f_D(2k+lambda x)
      H(x/(tau_1 sqrt(r))) H(k/(tau_2 sqrt(r)))
      integral_(0<=z<=lambda k+x) f_D(z)
          H((lambda k+x-z)/(2 tau_3 sqrt(r))) dz dk dx.

The square uses independence of the two sides in the centered ideal wave; no independence of their actual perturbations is assumed. The limiting adjacent zero-pair density is at least 3/(2pi)>.477 throughout interpolation since its correlation is at least w_1. The checker certifies the integral of one oriented lower probability is at least .003781166862163664. Counting two orientations and using .477 gives

    K_1 >= .003607233186504135.

### A different pivotal event at distance two

Suppose P=(u,v), the path between them is u-w-v, and the cover has girth sufficiently large that the displayed patch is a tree. Let t be the third neighbor of w. Require t to have the sign of w, and all four outer neighbors of u,v to have the opposite sign to w. No conditions on any more distant signs are needed.

On this event Delta_uv F<=-1. Remove u,v,w and let the exterior crossing matching size be M. Write A,B for the possible zero-or-one increases when u, respectively v, is individually restored with the color of w and can connect to its outer neighbors. If both centers have the other color, the only new edges form the two-edge star at w and the maximum is M+1. The two unequal assignments yield M+1+A and M+1+B: one center independently matches w. If both agree with w, their combined addition can increase exterior matching size by at most A+B. Indeed they add vertices in the same bipartition class. If A=B=0, a maximum matching of the exterior cannot gain an augmenting path starting at either new vertex individually, and an augmenting path cannot pass through the other new unmatched vertex. If exactly one indicator is zero, gaining two would contradict that zero on removing the other new vertex from the matching. The both-one bound is immediate. This proves the mixed-difference bound by the two equal values minus the two unequal values.

For the two endpoints conditioned to be zero in the ideal wave, w has variance 1/33 and t=lambda w by the wave equation. At each endpoint its two outer values are -w/2 plus and minus G/(2sqrt(3)) with two independent standard contrasts G, independent of the central value. More generally the regression of X_w on X_P is b(X_u+X_v), b=w_1/(1+w_2)=4sqrt(2)/11. Use the three ideal independent standard coordinates sqrt(33)(X_w-b(X_u+X_v)), G_u, G_v. For the base use their counterparts H_0=sqrt(33)(Y_w-b(Y_u+Y_v)), H_u,H_v. Before conditioning,

    Var H_0 <= 1 + 33·2 b^2(1-w_2) = 43/11 <3.91.

This follows by replacing the two negatively signed adjacent covariances by their lower bound w_1 and the positively signed endpoint covariance by one in the variance expansion. Both other variances are at most one. The tilt mass for these three coordinates is consequently at least (1-s)/sqrt(1+3.91r).

Given Z_P=0 and Y, the central score of Z/sqrt(1-s) is exactly the ideal reconstructed w using g+sqrt(r)H. The third neighbor's discrepancy from lambda w is sqrt(r)R_w; the two outer siblings at u, respectively v, have the shared discrepancy sqrt(r)R_u/2, respectively sqrt(r)R_v/2. For each fixed Y, ensure the event by requiring the symmetrically enlarged robust margins in the ideal reconstruction. This is an even event in the three driving coordinates, so translation, tilting, and then integrating over the ideal g give the same product-of-strips bound. Write f_0 for the |N(0,1/33)| density. At x=|w| the margin for the third neighbor is lambda x, and for each of the two outer pairs the allowed defect magnitude is x-z where z=|N(0,1/3)|<x. Thus

    E[-Delta_uv F | Z_P=0] >= (1-s)/sqrt(1+3.91r)
       integral_(x>=0) f_0(x) H(lambda x/(c sqrt(r)))
         [ integral_(0<=z<=x) f_D(z) H((x-z)/(c sqrt(r))) dz ]^2 dx.

This event already includes both signs of w, so there is no extra factor of two. The checker certifies its integrated probability at least .007684450176540460. Pair correlation is at least w_2 in the limit, so the zero density is at least 3/(pi sqrt(11))>.2878. Therefore each distance-two excess has coefficient

    K_2 >= .002211584760808344.

### Counting the distance-two excess without double counting

On the high-girth finite cover let b_y denote the unit key rows pulled back from the base (including the repeated base tie-noise coordinates if present). Temporarily measure the raw excesses from the exact critical values: let S_d be the sum of b-row covariance minus w_d over unordered pairs at cover distance d. Every distance-two pair has a unique common neighbor, since the cover girth exceeds four. Put E_y=sum_(x~y)b_x-lambda b_y. Since lambda=3w_1 and lambda^2=3+6w_2,

    sum_y <b_y,E_y> = 2 S_1,
    2 S_2 = sum_y (||sum_(x~y)b_x||^2-lambda^2)
          = 4 lambda S_1 + sum_y ||E_y||^2.

Consequently S_2>=2lambda S_1. For actual finite-profile base-minus-cover differences, replacing w_d by the finite high-girth covariances at d=1,2 costs at most an arbitrarily small constant times the cover order. All actual differences have the favorable sign. The total Price gain from the two disjoint pair sets therefore has effective coefficient K_1+2lambda K_2 for adjacent excess; exact matching-cut lifting then divides by the covering degree. Using 2lambda>5.65684 gives the coefficient stated at the start. If eta is the good-center cutoff, Section 7.5 of cutwidth.md charges the nine edges of the named patches, yielding the estimate min(6q,K eta/6) before arbitrarily small finite-profile errors.

### Exact reproduction and finite transfer

The submitted `certify.py` contains both integrals and their shared scalar helpers. Run:

    python -I -B certify.py

The shared exact scalar routines use a spatial step 1/200 and interpolation mesh s_j=(j/80)^2, j=1,...,79. Integrands and cell bounds decrease in s, so use the right endpoint and omit the untouched final tail. Each density is evaluated at a spatial cell upper endpoint and each strip at a lower amplitude endpoint. Inner remaining margins are rounded down before using a grid convolution. The program uses the conservative bounds 2.82842<lambda<2.82843 and c<.17158. For each absolute normal strip probability, round its argument down to .001 by integer square root, cap at 3.3, and integrate the degree-23 odd Taylor lower polynomial of exp(-t^2/2), with prefactor .7978. For an exponential density, use its seventh-degree odd Taylor lower polynomial at z/64, then square downwards six times. All positive fixed-point floors at scale 10^18 are downward bounds, and density prefactors are certified from pi<3.1416. The executable prints floor(K_1·10^18)=3607233186504135, floor(K_2·10^18)=2211584760808344 and floor((K_1+5.65684 K_2)·10^18)=16117814324835211.

The estimates derived for the ideal wave transfer to sufficiently accurate finite profiles up to any fixed error. Otherwise take a violating parameter/edge sequence, extract a convergent subsequence of the finite joint base-patch covariances and of s in a closed integration cell, and use convergence of the cover patch to the critical wave. The conditioning endpoint covariance stays nonsingular for s<=s_79<1. Each sign tested has positive conditional variance owing to the independent ideal cover residual (w and t each do as well in the distance-two case). Thus these finite unions of open sign cylinders have convergent probabilities, including the two-orientation count where relevant. This proves the required uniform limiting bounds without imposing any condition on exterior signs.
