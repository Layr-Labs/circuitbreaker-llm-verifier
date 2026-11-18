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

There is no assertion about convergence or exact evaluation of an infinite sum.

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

Below, E_(n,N) denotes the Gaussian-moment bound in (8). Its factors are exact integers or positive directed intervals; the density power is evaluated by repeated interval multiplication. The one-node Gauss normalization gives the midpoint bound L^3 sup|h''|/24.

## 4. Certified quadrature nodes, weights and arithmetic

Pi, root brackets, quadrature weights, density factors, analytic error bounds and the final alternating sum use Python's standard-library `decimal` with 100 decimal digits. Every Decimal lower endpoint uses `ROUND_FLOOR` and every upper endpoint uses `ROUND_CEILING`. Repeated positive tensor contractions use exact integers on the fixed-point grid proved in Section 5, with outward integer floor and ceiling. Square root and exponential are correctly rounded in Decimal; moving their rounded values one representable number outward encloses the exact transcendental value. Pi is enclosed by Machin's identity and 100-term alternating arctangent series, including an explicit next-term remainder.

Approximate Legendre roots are only starting guesses. The program brackets each root by rational decimal endpoints, evaluates the degree-`N` Legendre polynomial at both endpoints by interval recurrence, and verifies strictly opposite signs. The `N` brackets are disjoint, so they isolate all roots of a degree-`N` polynomial. Legendre's derivative is enclosed on each entire root bracket and the exact weight formula

\[
w=\frac2{(1-x^2)(P_N'(x))^2}
\]

is evaluated on that interval. Roots and weights are then mapped to `[0,L]`.

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

The implementation fixes the exact integer scale \(S=2^{128}\). An integer
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

## 6. Production parameters and certified inequality

The production defaults are N=101, J=39, L=3, 100-digit Decimal endpoints
and the exact fixed-point scale 2^128. The complete output is enclosed by

    0.130647566423167673 <= Gamma_39 <= 0.130647566423167695
                                       < 0.130647566423167696
                                       = 8165472901447981/62500000000000000.

The rational threshold in this certificate is the internal graph constant
gamma_0. The declared circuit coefficient c_* is strictly smaller than
1+1/(2 gamma_0), as proved in proof.md Section 5.

Both displayed endpoints are rounded outward from the complete computed
interval. The strict inequality with the target is exact rational arithmetic.
The program prints each of the 39 probability intervals, all analytic errors,
and every signed partial sum. It exits zero only after comparing the final
upper endpoint with the lower endpoint of the exact rational target. Its sole
declared JSON input is {}; no external data or provider credentials are used.

The graph argument fixes this odd J before taking any limit. Its finite
star-deletion error coefficient is

    C_39 = 1 + sum_(j=1)^39 (j+1)*3*2^(j-2)
         = 1 + 3*39*2^38 = 32160715112449.

This large but fixed constant affects the unspecified finite profiles and
thresholds. It does not change the uniform little-o quantifier. No convergence
or exact evaluation of an infinite alternating series is asserted.

## 7. Reproduction, validation and attribution

Run the exact declared defaults from the submission directory:

    python3 -I -B certify.py < /dev/null

The source SHA-256 is
`abcdb7b0840e8ff165db807cd163dc5e0396a8e12af8015a4d413072b87efe90`.
A local Python 3.12.14 run of these exact bytes, with an empty environment
and input {}, completed all 39 terms in 6.173 seconds, used 29,952 KiB peak
resident memory and returned exit zero with `certified: true`. That run
imposed 20-second wall and CPU limits, 256 MiB address space and 32 processes.
The exact upper endpoint begins
`0.1306475664231676940892403704138490751463631538463210494277`.
Its strict margin below the target is greater than 1.91e-18.

These are local host measurements. Docker is unavailable in this workspace,
so the digest-pinned hosted experiment must be completed afresh. Local timing
alone cannot establish another worker's runtime.

The exact-integer source is submission ac0dd6e1-8a73-4682-b7aa-67502b2de08b,
commit 77f07bfd2fd633cf98e3737575ba188d20e95c6a. It supplied the exact
fixed-point contractions, quadratic storage, and the J=23, N=88 certificate.
That package obtained a positive hosted mathematical verdict, but no earlier
judgment is a premise here. Its complete proofs are retained and reviewed anew.

The earlier submission 38e961ae-bd55-4566-85ad-17112333f84a supplied the
parameter refinement and positive factorization, but its hosted Decimal
contractions timed out. The original analytic and mathematical package comes
from fcf21648-1b1f-463c-b5bc-c50f53f63436 at commit
81470c9513639af985c88699dbcd2e2efbef7463. The immediate source for the fully proved positive-completion tail bound,
Gaussian even-moment derivative bound with exact precision diagonals, and
N=101, J=39 parameters is submission 0514f3af-5d65-4556-82df-a89e2a3ae4e1,
commit d18ef5c1b86c58ca9f0bdb93085f9b6f2698c55d. All four are declared
dependencies. This revision retains those estimates and parameters, changes
the final exact comparison to gamma_0=8165472901447981/62500000000000000,
and supplies the positive-margin eventual-bound argument in proof.md Section 5. Neither computation nor attribution replaces the
supplied graph, algebraic, uniformity and circuit proofs.
