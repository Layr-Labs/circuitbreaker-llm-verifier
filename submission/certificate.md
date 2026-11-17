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

All arithmetic uses Python's standard-library `decimal` with 100 decimal digits. Every lower endpoint uses `ROUND_FLOOR` and every upper endpoint uses `ROUND_CEILING`. Square root and exponential are correctly rounded in Decimal; moving their rounded values one representable number outward encloses the exact transcendental value. Pi is enclosed by Machin's identity and 100-term alternating arctangent series, including an explicit next-term remainder.

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

Every term is nonnegative. Directed Decimal fused multiply-adds therefore give lower and upper enclosures of each contraction. Summing the final matrix encloses the complete tensor quadrature. Subtracting/adding `E_{n,N}`, adding the one-sided tail `T_n`, and multiplying by two encloses `p_j`. The finite signed sum is then evaluated with interval arithmetic, including the correct endpoint reversal for its negative terms.

For speed, the transition density is evaluated by the exact factorization

\[
k(z\mid x,y)=\sqrt{6/\pi}e^{-3x^2/2}e^{-12y^2}e^{-6z^2}
 e^{-6\sqrt2xy}e^{-6xz}e^{-12\sqrt2yz}.
\]

All exponential factors depend on at most two coordinates, so only O(N^2) exponential evaluations are required. The program precomputes the positive interval factors and combines each triple with directed multiplications. This is an algebraic identity, with every multiplication error still enclosed; it does not change the quadrature rule or its analytic error bound.

The implementation does not assume correctly rounded platform binary `libm` transcendental functions: the initial binary cosine guess is irrelevant after the root-bracket sign checks succeed.

## 5. Final production parameters and enclosure

The production defaults are `N=72` Gauss nodes and the fixed odd truncation `J=13`. The directed interval result is contained in the outward-rounded interval

\[
0.130648076091\le\gamma_{13}\le0.130648091104<0.13065=\frac{2613}{20000}.
\]

The program prints the exact endpoint strings and success flag in its JSON output. These parameters reduce memory use while keeping the analytic formulas above unchanged and retaining more than 0.0000019 of strict slack below the target.

## 6. Declared finite certificate and recorded execution

The submitted program uses N=72, J=13, L=3 and 100-digit interval endpoints.
Its sole declared input is the empty object {}; it needs no external files,
network, third-party packages or random samples. The finite result is

    Gamma_13 < 0.130648091104 < 2613/20000 = 0.13065.

The first displayed decimal is rounded upward from the exact upper endpoint
printed by the certificate. Every probability interval and intermediate signed
sum is also printed. The final program exits zero only if its upper endpoint
is strictly less than the lower interval endpoint for the exact target.

The source SHA-256 is `245749f29e54f010a0e978ba287120938995c3325327d35d1b946b4b3a4666f5`. A local Python 3.12 run completed in 4.57 seconds,
with peak resident memory about 172 MiB. The pinned Linux/amd64 Docker execution
has not run locally because Docker is unavailable. This reported local
computation is distinct from hosted acceptance.
