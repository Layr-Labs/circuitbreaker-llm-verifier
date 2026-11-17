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

## 5. Quadratic-storage contraction

For the nodes x_i and weights w_i define the nonnegative factors

    A_ij = exp(-3*x_i^2/2 - 6*sqrt(2)*x_i*x_j),
    B_ik = exp(-6*x_i*x_k),
    C_jk = sqrt(6/pi)*w_k*exp(-12*x_j^2-6*x_k^2-12*sqrt(2)*x_j*x_k).

The factorization in Section 4 gives exactly

    w_k*k(x_k | x_i,x_j) = A_ij*B_ik*C_jk.

Consequently the same tensor-quadrature update is

    T_ij = M_ij*A_ij,
    M'_jk = C_jk*sum_i T_ij*B_ik.

The code stores only these two-index arrays and the current state. Each lower
endpoint uses downward multiplication and fused multiply-add on lower
endpoints; each upper endpoint uses the corresponding upward operation on
upper endpoints. Every factor and state value is nonnegative, so induction
through the products and sums encloses the exact expression above. Pulling
C_jk outside the sum is an algebraic identity. There is no approximation of
the transition kernel. Storage is O(N^2), and arithmetic is O(J*N^3).
All bounds in Sections 1-4 continue to apply to this same finite tensor rule.

## 6. Production parameters and certified inequality

The production defaults are N=88, J=23, L=3 and 100-digit endpoints. The exact
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
`94c26fc316c7ef8cedc5513dc1361dedbe39662c4159dc207a472ad6c89c56da`.
A local Python 3.12.14 run of those exact bytes completed in 13.603 seconds,
with peak resident memory 28,272 KiB, exit code zero and `certified: true`.
Host measurements are not a guarantee of the hosted worker's runtime.
The digest-pinned Linux/amd64 Docker replay is unavailable locally because
Docker is absent; hosted execution remains a separate required check.

The analytic certificate and original implementation are inherited from
submission fcf21648-1b1f-463c-b5bc-c50f53f63436 at source commit
81470c9513639af985c88699dbcd2e2efbef7463. The changes are the finite parameters,
target rational and contraction order described above. At the original
N=72, J=13 settings, the factorized implementation and original implementation
produced the identical final upper endpoint to all 100 printed decimal digits;
individual probability endpoints differed by at most 2e-102. That cross-check
supports the implementation change; its mathematical justification is the
nonnegative directed-interval argument in Section 5. Neither the cross-check
nor the final computation substitutes for the written graph and circuit proofs.
