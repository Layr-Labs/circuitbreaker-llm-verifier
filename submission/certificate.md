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

## 2. Truncation to the positive cube

Take `L=3`. In every adjacent pair, on the positive quadrant,

\[
f_2(x,y)\le\frac3{2\pi}e^{-9(x^2+y^2)/2}.
\]

The event that any of the `n=j+1` positive coordinates exceeds `L` implies this event for that coordinate and any chosen adjacent coordinate. A union bound and the elementary Gaussian tail bound therefore give an omitted orthant probability of at most

\[
T_n=n\frac3{2\pi}\frac{e^{-9L^2/2}}{9L}\sqrt{\frac\pi{18}}.
\]

This is the error for **one** positive orthant; the final alternating probability multiplies it by two. No independence is asserted in this union bound.

## 3. Tensor Gauss–Legendre error

Let `f_n` be the joint Gaussian density above and `Q_n` its precision matrix. From the displayed factorization,

\[
(Q_n)_{ii}\le39,\qquad f_n(0)<2^n.
\]

For any real vector `x` and complex perturbation in just coordinate `i`, positive definiteness gives

\[
|f_n(x+ze_i)|\le2^n\exp(39|z|^2/2).
\]

Cauchy's derivative estimate with `r^2=2N/39` and `e<3` gives

\[
|\partial_i^{2N}f_n(x)|\le(2N)!2^n
\left(\frac{117}{2N}\right)^N.
\]

The `N`-node Gauss–Legendre error on `[0,L]` is bounded by

\[
\frac{L^{2N+1}(N!)^4}{(2N+1)((2N)!)^3}
\sup|h^{(2N)}|.
\]

This normalization can be checked at `N=1`, where it is the midpoint bound `L^3 sup|h''|/24`.

Both exact integration and positive Gauss quadrature have operator norm `L` in the supremum norm. Telescoping the difference of their `n`-fold tensor products yields the uniform error bound

\[
E_{n,N}=\frac{nL^n2^n}{2N+1}
\frac{(N!)^4}{((2N)!)^2}
\left(\frac{117L^2}{2N}\right)^N.
\]

The program computes this expression using directed interval arithmetic and exact integer factorials. This deliberately conservative bound avoids any empirical smoothness estimate.

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
The analytic truncation and Gauss-quadrature error bounds in Sections 1--4
are unchanged.

## 6. Production parameters and certified inequality

The production defaults are N=88, J=23, L=3, 100-digit Decimal endpoints and
a fixed-point scale of 2^128 for the contractions. The exact
output interval is contained in the outward-rounded interval

    0.130647566408153 <= Gamma_23 <= 0.130647566543255
                                    < 0.130647566544
                                    = 8165472909/62500000000.

The last strict inequality is rational arithmetic. The program prints its
full endpoint strings, each path-probability interval, all error bounds and
all signed partial sums. It compares the final upper endpoint with the lower
endpoint of the exact rational target and exits zero only if the comparison
is strict. The sole declared JSON input is {}; no network, external data,
third-party package or random sample is needed.

For the graph argument the odd truncation is fixed first. The associated
star-deletion error coefficient is

    C_23 = 1 + sum_(j=1)^23 (j+1)*3*2^(j-2)
         = 1 + 3*23*2^22 = 289406977.

Its size affects the unspecified finite-profile parameters, not the validity
of the uniform little-o limit. The certificate makes no assertion about an
infinite alternating series or an explicit asymptotic remainder rate.

## 7. Reproduction, validation and attribution

Run the exact declared defaults from the submission directory:

    python3 -I -B certify.py < /dev/null

The source SHA-256 is
`456698dd39c9c7e1398ec8e07e97feeccb839f517654ab6f8741069a09ffd5e2`.
A clean-environment local Python 3.12.14 run of these exact bytes completed
in 4.681 seconds, with peak resident memory 25,856 KiB, exit code zero and
`certified: true`. The final upper endpoint begins
`0.1306475665432542427318569944774497790977467040232638543164`.
Host measurements do not guarantee another worker's runtime. Docker is
unavailable locally; digest-pinned hosted execution remains a separate check.

The preceding version, submission 38e961ae-bd55-4566-85ad-17112333f84a,
used Decimal fused multiply-adds in every inner contraction. Although its
local runs completed, both captured hosted executions timed out before
J=23 and its unconditional theorem was rejected for the missing completed
numerical witness. This revision replaces those inner contractions with the
exact integer enclosure proved in Section 5. It retains the same theorem,
quadrature rule, finite path count, analytic error bounds and target rational.
Its completed local intervals enclose all 23 earlier local probability and
partial-sum intervals, with outward changes below 4e-28. This comparison is
an implementation check; the new computation and Section 5's induction
supply its own numerical witness and rounding justification.

The analytic certificate and original implementation are inherited from
submission fcf21648-1b1f-463c-b5bc-c50f53f63436 at source commit
81470c9513639af985c88699dbcd2e2efbef7463. Both that source and the preceding
submitted version are declared as dependencies. No earlier judgment is a
mathematical premise. The current contribution is exact fixed-point
contraction, in addition to the already stated parameter refinement and
quadratic storage. Neither numerical agreement nor computation substitutes
for the supplied graph, algebraic and circuit proofs.
