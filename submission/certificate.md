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

## 6. Parameters, exact comparisons and reproduction

Run the standard-library program with Python 3 and its default parameters. The declared input is the empty JSON object. It checks the 93 Gauss nodes on [0,3] for each of 47 path lengths, with 100-decimal-digit directed enclosures and exact 224-bit fixed-point grid arithmetic. The output verifies

    0.1306475664231676642 < Gamma_46 <= Gamma
        <= Gamma_47 < 0.130647566423167665.

The even lower and the odd upper are separate comparisons in the code. The graph constant also uses the 0.000000047 saving proved in cutwidth.md Section 7.6. A separate exact rational Riemann sum in the same executable gives a floor of 18079200775 for 10^12 times its two-private-edge event lower sum. The corresponding probability is greater than .018 as used in the structural proof.

    python -I -B certify.py

With the default settings the program exits successfully only if both series comparisons and the pivotal check pass. It prints one JSON object with each alternating-path probability enclosure and the final comparisons. The proof of the uniform graph improvement and the circuit transfer is in the companion manuscripts.

The directed quadrature, Markov factorization, Gaussian moment bound and integer contraction derive from CircuitBreaker submissions fcf21648, 38e961ae, ac0dd6e1, 0514f3af and 3f2be7ea. The paired completion, coordinate envelope, root enclosure optimization, sparse arithmetic and pivotal sum are used here to certify the sharper bound.


## 7. Finite checks for the strengthened retiming argument

The same executable evaluates the good-center lower sum described in certificate.md Section 8 and the integrated tilted conditional probabilities in certificate.md Section 9. All entries of those two computations are either rational or rounded downward by integer divisions in a proved lower bound. The good sum must exceed 15 times 0.00000002. The integrated per-edge saving must have floor 672905379119700 at scale 10^18. The original .04-margin pivotal sum remains available for the preparatory Section 7.4. The program also checks the odd upper and even lower enclosures, then emits one JSON record.


## 8. A quantitative probability bound at good centers

Let $\eta=17/40000=.000425$. In the notation of the retiming rule, take an arbitrary fixed naming of neighbors $a,b,u$ at $o$, and of children $w,w'$ of $u$. Use $a_1,a_2$ and $w_1,w_2$ for the two outward children of $a$ and $w$. If all lifted edges in the radius-three tree ball have Hellinger excess at most $\eta$, then at every $0\le t\le .67$,

\[
 \Pr\{Z_o>t,\ Z_a,Z_b,Z_{a_1},Z_{a_2},Z_w,Z_{w_1},Z_{w_2}<t\}>3.6190\,10^{-7}.
\]

The assertion holds for the finite-profile, positive-tie-noise field once the parameters are sufficiently accurate, uniformly over the graph. Reflection gives the negative thresholds. The original conflict graph selects at least $g/15$ of the $g$ good vertices. Hence their simultaneous expected improvement exceeds $2.4126\,10^{-8}g$. One fixed orientation at each center suffices.

#### Row estimates

Work first with the complete normalized radial rows $f_x(z)\propto\alpha^{d(x,z)}$ on the universal tree and their pushed rows $b_x(p)=\|f_x|_{\pi^{-1}(p)}\|_2$. Write $\lambda=2\sqrt2$, $\lambda_\alpha=\alpha^{-1}+2\alpha$ and $w_d=(1+d/3)2^{-d/2}$. All $o(1)$ terms below are uniform in the base graph and tend to zero as $\alpha\uparrow1/\sqrt2$, followed by sufficiently accurate truncation and noise choices.

For vertices joined by a good tree path of length $d$ put $E_{xy}=\langle b_x,b_y\rangle-\langle f_x,f_y\rangle$. Then

\[
 0\le E_{xy}\le d^2\alpha^{-(d-1)}\eta=(k_d+o(1))\eta,
 \qquad k_d=d^2(\sqrt2)^{d-1},\quad k_0=0.                    \tag{1}
\]

Indeed, if $e_x(p)$ is the normalized fiber direction, $2E_{xy}=\sum_p b_x(p)b_y(p)\|e_x(p)-e_y(p)\|^2$. For successive vertices $i-1,i$ in the path the ratio $b_x(p)b_y(p)/(b_{i-1}(p)b_i(p))$ is at most $\alpha^{-(d-1)}$ by the adjacent norm-ratio bounds. Telescope the unit directions in the weighted Euclidean norm and use $E_{i-1,i}\le\eta$. As $w_d+k_d\eta<1$ for $1\le d\le6$, distinct tree vertices in the radius-three ball have distinct base projections. Equivalently the finite field covariances satisfy $c_{xy}\in[w_d-o(1),w_d+k_d\eta+o(1)]$. All local paths used here remain inside that ball.

For each of $v=o,a,u,w$ the centered Gaussian wave defect $D_v=\sum_{y\sim v}Z_y-\lambda Z_v$ satisfies

\[
   \operatorname{sd}(D_v)\le\sqrt{3(3-\lambda)\eta}+o(1)<.01480.       \tag{2}
\]

To see this let $R_v=\sum_{y\sim v}b_y-\lambda_\alpha b_v$. Except at $z=v$, the sum of the three neighbor entries at $z$ over the tree is $\lambda_\alpha f_v(z)$ and the sum of their squares is $(\alpha^{-2}+2\alpha^2)f_v(z)^2$. At $z=v$ the linear formula is short by $q_\alpha=(\alpha^{-1}-\alpha) f_v(v)=o(1)$ and the sum of squares is smaller. The fiber triangle inequality gives $R_v+q_\alpha\mathbf1_{\pi(v)}\ge0$ coordinatewise. The fiber sums of squares and Cauchy–Schwarz give

\[
 R_v\le [\sqrt{3(\alpha^{-2}+2\alpha^2)}-\lambda_\alpha] b_v.
\]

Also $\langle b_v,R_v\rangle\le3\eta+3r_1(\alpha)-\lambda_\alpha=3\eta+o(1)$. Multiply the upper bound for the nonnegative vector $R_v+q_\alpha\mathbf1_{\pi(v)}$ by that same vector coordinatewise and sum. This gives $\|R_v\|^2\le (3-\lambda)3\eta+o(1)$, proving (2). In passing to the actual field the finite truncation and tie noise perturb each of this bounded number of rows by an arbitrarily small amount. The finite root variance is one. All rounded inequalities used subsequently are strict at the limiting parameter values.

#### Three coordinates and their conditional means

Set $O=Z_o$, $R=2Z_u-Z_a-Z_b$, $Q=\sqrt3(Z_{w'}-Z_w)$. Subtract their respective regressions on $O$ from $R,Q$, and call the covariance of the resulting pair $K$. In the ideal critical wave the three starting coordinates have covariance $I_3$. Write $O=t+x$ and take the two residual coordinates to equal $K(y,z)^T$, with $x,y,z>0$. Their joint density in $dx\,dy\,dz$ is

\[
 (2\pi)^{-3/2}\sqrt{\det K}\exp[-((t+x)^2+(y,z)K(y,z)^T)/2].      \tag{3}
\]

The mean of $Z_v$ under these conditions is

\[
 c_{vo}(t+x)+[\operatorname{Cov}(Z_v,R)-c_{vo}\operatorname{Cov}(O,R)]y
             +[\operatorname{Cov}(Z_v,Q)-c_{vo}\operatorname{Cov}(O,Q)]z.
\]

The $t$ term is at most $t$. The coefficients on $x,y,z$ are bounded **above** as follows, including sufficiently small finite approximations:

| $v$ | $x$ | $y$ | $z$ |
|---|---:|---:|---:|
| $a$ or $b$ | .9433 | −.1610 | .0172 |
| $a_1$ or $a_2$ | .8358 | −.2196 | .0368 |
| $w$ | .8358 | .4730 | −.2810 |
| $w_1$ or $w_2$ | .7148 | .5055 | −.3920 |

Here are details to reproduce the table and density bound. Compute each covariance of two linear forms by summing products of the coefficients times $w_d$ and add error up to $\eta$ times the sum of the *positive* products times $k_d$, or down to the corresponding negative sum, by (1). Thus the error multipliers for $\operatorname{Cov}(O,R)$ and $\operatorname{Cov}(O,Q)$ are $\pm2$ and $\pm4\sqrt6$. For $\operatorname{Var}R,\operatorname{Var}Q,\operatorname{Cov}(R,Q)$ the ranges are $[-8k_2,2k_2],[-6k_2,0],\pm\sqrt3(2+2k_3)$. For the four rows in order the ideal $(c_{vo},\operatorname{Cov}(Z_v,R),\operatorname{Cov}(Z_v,Q))$ are

\[
 (w_1,-1/6,0),\ (w_2,-\lambda/12,0),\ (w_2,\lambda/6,-1/(2\sqrt3)),\ (w_3,1/2,-\lambda/(4\sqrt3)).
\]

Their one-sided error multipliers from (1), before the regression subtraction, are respectively $(k_1,2k_2,\sqrt3 k_3)$, $(k_2,2k_3,\sqrt3k_4)$, $(k_2,2k_1,\sqrt3k_2)$ and $(k_3,2k_2,\sqrt3k_3)$. For each $v$ use $(w_d+\eta k_d)2\eta$ or $(w_d+\eta k_d)4\sqrt6\eta$ to bound the regression subtraction above. This gives the four displayed upper bounds. The same estimates give $K_{11}\le1.005$, $K_{22}\le1+o(1)$, $|K_{12}|<.029$ and $(2\pi)^{-3/2}\sqrt{\det K}>.0615$. In particular, for all sufficiently accurate finite profiles,

\[
 (y,z)K(y,z)^T\le1.035y^2+1.030z^2\quad (y,z\ge0),
\]

by $2yz\le y^2+z^2$ with the strict slack just displayed.

#### The remaining centered Gaussian fluctuations

Condition on these three coordinates and subtract conditional means from every remaining coordinate. Use defect tolerances $|D_v-\mathbb E[D_v\mid O,R,Q]|\le r$ at the four stars. Also constrain the centered halves of the three differences $Z_a-Z_b,Z_{a_1}-Z_{a_2},Z_{w_1}-Z_{w_2}$ in absolute value by $d_1,d_2,d_3$ respectively. The root fluctuation is zero, as are those of $R,Q$. The star identities and the fixed sibling differences now bound the absolute fluctuations at the four rows of the table by

\[
 r/3+d_1,\quad 1.415(r/3+d_1)+r/2+d_2,\quad .972r,\quad 2.043r+d_3.
\]

Indeed the fluctuation at $u$ is the root defect divided by three. At $w$ its bound is $(\lambda/6+1/2)r$ since $Q$ is fixed. At a child of $w$ the midpoint has bound $(\lambda/2)(.972r)+(1/6+1/2)r<2.043r$.

Conditional variances are at most the unconditional variances. By (1) each half-difference has standard deviation at most $1/\sqrt{12}+o(1)<1/3.46$. The Gaussian rectangle inequality gives a product lower bound for the seven constraints, even with correlated or singular residuals. For completeness, given one centered Gaussian coordinate $X$ the other coordinates have conditional law $W_0+vX$ with $W_0$ an independent centered Gaussian. For a symmetric rectangle $C$, $h(s)=\Pr(W_0+vs\in C)$ is even and decreases with $|s|$. In the nondegenerate case write its integral by layer cake over centered ellipsoidal density level sets $L$. The nonempty intersections $C\cap(L+vs)$ have even volumes whose $1/d$ powers are concave in $s$ on their support by Brunn–Minkowski, so each decreases on $s\ge0$. Consequently $\mathbb E[h(X)\mid |X|\le c]\ge \mathbb E h(X)$. Induction over coordinates proves the product bound. Add independent Gaussian noise tending to zero to get the singular case; positive rectangle faces have zero mass for centered one-dimensional Gaussian marginals, including those of zero variance. Zero-width constraints can simply receive probability lower bound zero.

For $d_i>0$ the checker uses $H(3.46d_i)$, where $H(u)=.797v(1-v^2/6)$, $v=\min(1,\max(0,u))$. This follows from $e^{-s^2/2}\ge1-s^2/2$ and $\sqrt{2/\pi}>.797$. For defect tolerances $r=.01705,.02370,.03185$, (2) permits standardized thresholds $1.15,1.60,2.15$, giving lower masses $.748,.889,.967$ respectively. For example integrate $.797\sum_{j=0}^{11}(-s^2/2)^j/j!$ from zero to the corresponding threshold, obtaining lower bounds $.74902,.88941,.96736$. Raising the stated masses to the fourth power handles the four defects.

#### Finite rational lower sum

`certify.py (good_trigger_certificate)` uses boxes of side lengths $1/200,1/40,1/40$ in $x,y,z$ on $[0,.2]\times[0,3.5]\times[0,6.5]$. On each cell take the table's worst signed endpoints and negate them to obtain margins $A,B,C,D$. For each of the three tolerances with $C>.972r$ and $D>2.043r$ put

\[
 U=\min\{A-r/3,(B-r/2)/1.415-r/3\},\quad
 d_1\in\{.5U,.8U\},\quad d_2=B-r/2-1.415(d_1+r/3),\quad d_3=D-2.043r.
\]

Discard nonpositive cases. For every valid choice the centered rectangle fits inside the trigger, using its conditional mean bound. Equalities have probability zero in the unconditional field with positive tie noise. Take the highest lower bound among these alternatives, multiply by (3) at its conservative upper quadratic and upper box endpoints, and sum. The checker bounds each required exponential $e^{-v}$ from below by $(1-v/4096)^{4096}$ with downward fixed-point rounding after every squaring. Every margin and product is evaluated downwards or with the unfavorable endpoint. It returns

\[
 \left\lfloor 10^{15}\times(\text{lower sum})\right\rfloor=361900864.
\]

Thus each good center has the claimed probability at each central threshold. Finite-radius concentration, mesh interpolation, and key-tie handling are unchanged.


## 9. Integrated covariance gain at a heavy edge

This replaces the small-noise estimate in the adjacent pivotal-witness argument. The notation of that argument is as follows. X is the field on a sufficiently high-girth finite cover, Y is an independent pullback of the Gaussian base field, P=(i,j) is an adjacent pair, and Z=sqrt(1-s)X+sqrt(s)Y is conditioned on Z_P=0. The two sign-reversed private-edge witnesses each force the centered mixed difference of the crossing-cover function F to be at most -1. If both occur its value is -2, since all four outer neighbors can simultaneously be matched to their private leaves and the only additional central cost is the crossing indicator of ij.

The integrated gain per unit of raw adjacent covariance excess is at least

    K = 0.000672905379119700.

As usual arbitrary fixed finite approximation losses can be absorbed in the o(n) conclusion. Exact lifting divides by the degree of the cover. The same gain applies at every threshold by the unmodified F's zero-threshold maximum.

### Uniform control of the wave equation error

The ideal cover wave has lambda=2sqrt(2) and covariance w_d=(1+d/3)2^(-d/2). Write f_x for the full normalized positive alpha^distance tree row and b_x for its pushforward by fiber norms, with alpha increasing to 1/sqrt(2). For all source coordinates z except x exactly one neighbor is nearer to z. Therefore

    sum_{y~x} f_y = (alpha^(-1)+2alpha) f_x
                    -(alpha^(-1)-alpha) f_x(x) 1_x,
    sum_{y~x} f_y^2 <= (alpha^(-2)+2alpha^2) f_x^2

pointwise. Here f_x(x) tends to zero. On each fiber, the triangle inequality for the first expression and Cauchy-Schwarz on the three norms in the second imply

    (alpha^(-1)+2alpha)b_x - o(1) 1_{projection(x)}
       <= sum_{y~x} b_y <= sqrt(3(alpha^(-2)+2alpha^2)) b_x.

All row norms equal one. Thus in every limit of finite local base Gaussian covariances,

    sd(R_x) <= c=3-lambda,    R_x=sum_{y~x} Y_y-lambda Y_x.

The normalized finite truncation and independent tie noise approximate the full-profile rows uniformly, as in the submitted proof, so these bounds apply in the near-critical limiting argument. For any two siblings, their contrast H=sqrt(3)(Y_child1-Y_child2) has variance at most one in this limit: covariance domination gives Cov(Y_child1,Y_child2)>=w_2=5/6. Both variance bounds remain upper bounds after conditioning on Z_P=0.

### Symmetric translation of six independent contrasts

It suffices first to work with the exact ideal cover wave X and any limiting Gaussian pullback Y satisfying these bounds. A label choice specifies one of the two outer neighbors at each endpoint and one of the two outward children of each selected neighbor. Thus there are sixteen choices. For a fixed choice, the two orientations of the private-edge event require signs on the four outer neighbors, both children of each selected parent, and both children of each selected opposite child. Require margins m_d>0 at depths d=1,2,3. Under either orientation, success uniquely specifies the four label choices. Let W(g) be the sum of the two corresponding indicators in an ideal wave with P=0. Then W is even. It depends on six independent standard contrasts g: one at each endpoint, one at each of the two selected parents, and one at each of the two selected opposite children. All are independent of the original cover central pair X_P. A sibling pair in X has mean (lambda X_parent-X_grandparent)/2 and half-difference g/(2sqrt(3)).

Condition now on Z_P=0 and on the entire Y, so that the standard contrast vector g from X still has its original law. Put r=s/(1-s) and let H be the corresponding six contrasts of Y. Reconstruct the patch through the ideal recursion with zero endpoints and contrast vector g+sqrt(r)H, then multiply by sqrt(1-s). Its discrepancy from the actual Z values consists entirely of sqrt(s) times the wave defects of Y. For any sibling pair the two errors are identical. Up to the factor sqrt(s), at the three depths the discrepancies are respectively

    R_i/2,
    (lambda/4) R_i + R_a/2,
    (3/4)R_i+(lambda/4)R_a+R_z/2,

where a is the selected outer parent, z its selected opposite child; the other side uses j in place of i. These formulas follow by substituting the ideal and actual child recursions successively (lambda^2=8). Their standard deviations even before conditioning are bounded strictly by .086,.208,.336, respectively, by the triangle inequality and the bound c. For the fixed labels it suffices to control two absolute errors at each depth, since siblings share the error.

For any nonnegative even W and an independent standard Gaussian vector g,

    E W(g+h) >= exp(-||h||^2/2) E W(g).

Translate the Gaussian density and average the two signs of its linear exponent under g -> -g, using cosh>=1. If all six absolute wave errors (without sqrt(s)) are smaller than (m_d-.0001)/sqrt(r), whenever W(g+sqrt(r)H) is positive the corresponding actual signs are successful. Thus, averaging over Y conditional on Z_P=0, the probability of surviving labeled witnesses is at least E W(g) times

    E [ exp(-r||H||^2/2) 1[all six errors small] | Z_P=0 ].

The exponential tilt of this centered conditional Gaussian law is still Gaussian and decreases every variance. Indeed write the conditional Y as B times a standard normal vector, the six contrasts as A times that vector, and note that the tilted driving covariance is (I+r A^T A)^(-1)<=I. The total mass of the tilt is det(I+r Cov(H|Z_P=0))^(-1/2), at least (1+r)^(-3)=(1-s)^3 by Hadamard's inequality and the six unit upper variance bounds. Apply a union bound to the two absolute errors at each depth under the tilted law. This lower-bounds the last display by

    (1-s)^3 [1-4 sum_{d=1}^3 Phi(-(m_d-.0001)/(sigma_d sqrt(r)))],
    (sigma_1,sigma_2,sigma_3)=(.086,.208,.336).

The bound does not depend on the label choice. Sum the sixteen choices. The ideal total E W equals 2 Q(m_1,m_2,m_3), where Q is the one-orientation probability given by the submitted two-dimensional half-normal integral, with the lower amplitude endpoints x>=m_1, k>=m_2 and distribution argument lambda k+x-2m_3. At most one successful labeling per orientation can contribute on any actual realization. The two orientations' sum is dominated by -Delta F, including when both succeed.

This limiting lower bound transfers uniformly over the finite cover/base pairs for s in each closed cell in [0,.42], up to an arbitrarily small loss. Otherwise extract a violating sequence with finite patch covariance of X tending to the ideal wave, local Y covariances convergent by compactness, and s tending to a value in that cell. The two-by-two conditioning covariance stays invertible because s<=.42 and the cover's limiting pair is nonsingular. The limiting Y satisfies both structural variance bounds proved above. Conditional Gaussian convergence applies, and every exterior score in the patch has positive conditional variance, since the (1-s) times the independent cover residual variance given its own endpoints is positive. Therefore the probabilities of the finitely many strict sign cylinders converge, including their two-orientation sum. This contradicts the limiting uniform bound. The graph outside the patch is arbitrary throughout.

### Reproducible finite arithmetic

Use the following cells. The end coordinate is in units of .00001 and margins are in units of .01. The preceding cell ends at the left endpoint, starting at zero.

      50  2  2  2        100  2  2  4        200  2  4  4
     300  2  4  4        500  2  4  6        800  2  6  6
    1200  2  6  8       1800  2  8 10       2500  4  8 10
    3500  4 10 12       5000  4 10 14       7500  4 14 18
   10000  6 14 18      14000  6 16 22      18000  6 20 26
   24000  8 22 30      32000 10 24 34      42000 12 30 42

The companion exact checker evaluates Q from below on the region x<=.65 and k<=.85, tessellated by squares of side .01. Here f_A is the |N(0,1/12)| density, and f_D,H_D are the |N(0,1/3)| density and distribution. The formula evaluated is

    Q=[ integral_(x>=m1) f_A(x) integral_(k>=m2)
        2 f_D(2k+lambda x) H_D(lambda k+x-2m3) dk dx ]^2.

Indeed x is the amplitude at each pair of outer neighbors and k is the amplitude of the chosen neighbor's opposite child. The latter has children of its own sign with margins m3 precisely when their contrast has absolute value at most lambda k+x-2m3. The other child of the private parent has amplitude lambda x+k>=m2. The two ends are independent in the conditioned ideal wave. Use upper square endpoints for decreasing densities, lower square endpoints for the CDF, and the bounds 2.828<lambda<2.829. The half-normal coefficients exceed 11/4 and 11/8 by pi<22/7. A lower exponential is the odd Taylor polynomial through degree seven at z/64, followed by six downward-rounded squarings. For the CDF integrate the odd Taylor lower polynomial for exp(-3t^2/2) through degree nineteen in 3t^2/2, capping the argument at 5/4 by monotonicity; at nonpositive arguments the CDF is zero. The checker floor-rounds positive terms at scale 10^18. For a Gaussian tail upper bound, round the positive argument down to .001 using integer square root, cap it at three, then subtract from one half the integral of the odd Taylor lower polynomial for exp(-t^2/2) through degree nineteen in t^2/2, times .3988<1/sqrt(2pi).

For each cell the checker subtracts 10^(-9) from its Q lower bound, multiplies by the positive survival and tilt lower bound at the right endpoint and by the cell length, then sums. Reproduce with standard Python 3:

    python -I -B certify.py

The integers at scale 10^18 for the one-orientation integral and per-edge credit are respectively 715856786297553 and 672905379119700, rounded down. To pass from the former to the latter, double it and multiply by the adjacent zero-pair density lower bound .47 (the cover's limiting correlation rho gives density 3/(2pi)). No simulation is used.

An edge can be observed by at most 14 root centers (one of its endpoints must be within distance two, giving at most 2+4+8 preimages of those centers). Therefore at heavy-edge threshold eta a base graph with g good centers has sum_e D_e >= (n-g)eta/14. Uniformly approximate each heavy excess to retain .99eta and add the simultaneous good-center retiming credit qg/15. The asymptotic saving is thus at least min(q/15, 0.0000475840232377502 eta), up to the additional explicitly chosen numerical slack.

The executable also checks the two-pass rational comparisons used by Section 7.6 of cutwidth.md. The probability estimates of Sections 8–9 are the same ones used in the preceding submission.
