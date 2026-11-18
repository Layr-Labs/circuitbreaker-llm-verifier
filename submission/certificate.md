# Certified finite path sum with localized quadrature errors

This file proves the complete finite numerical premise used by cutwidth.md.
Its analytic reductions, quadrature validation, truncation errors and arithmetic
enclosures are all supplied. The standard-library computation certify.py is an
arithmetic witness, and is combined with the separately supplied graph,
algebraic, circuit and uniformity proofs.

## 1. The stationary path process and finite goal

Put c_2=3/(2*pi) and c_1=sqrt(6/pi).

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
\Gamma_J=\sum_{j=1}^J(-1)^{j+1}3\,2^{j-2}p_j.
\]

Only this finite sum is a premise of the submitted theorem; no infinite-series limit is used.


For `n>=2` put

\[
f_n(x_0,\ldots,x_{n-1})
=f_2(x_0,x_1)\prod_{r=2}^{n-1}k(x_r\mid x_{r-2},x_{r-1}).
\]

The retained Gaussian-path argument identifies

\[
p_j=2\int_{(0,\infty)^{j+1}}f_{j+1}(x)\,dx,
\qquad \Gamma_J=\sum_{j=1}^J(-1)^{j+1}3\,2^{j-2}p_j.
\]

Only the fixed odd value `J=65` is needed here. The proposed finite premise is

\[
\Gamma_{65}<\gamma_0,
\qquad
\gamma_0=\frac{130647566423167664341749381}{10^{27}}
=0.130647566423167664341749381.
\]

The declared program uses the exact 80-node Gauss--Legendre rule on `[0,3]`,
100-digit directed Decimal initialization, and exact integer contractions on
the `2^128` grid. Its analytic changes are proved next.

## 2. Reversible anchoring

Direct expansion gives

\[
f_2(x,y)k(z\mid x,y)=f_2(y,z)k(x\mid z,y)
\]

because either side is `c_2*c_1` times the exponential with exponent

\[
-6x^2-\frac{33}2y^2-6z^2
-12\sqrt2xy-12\sqrt2yz-6xz.
\]

Thus the joint density can be anchored at any adjacent pair `a,a+1`:

\[
f_n=f_2(x_a,x_{a+1})
\prod_{r=0}^{a-1}k(x_r\mid x_{r+2},x_{r+1})
\prod_{r=a+2}^{n-1}k(x_r\mid x_{r-2},x_{r-1}).\tag{1}
\]

To verify (1), successively move the initial pair one place to the right
using the displayed three-coordinate identity. Each move leaves behind the
corresponding reverse transition. All other factors are unchanged. Empty
products have value one. Taking the first forward transition into the anchor
also gives, for `0<=a<=n-3`,

\[
f_n=f_3(x_a,x_{a+1},x_{a+2})
\prod_{r=0}^{a-1}k(x_r\mid x_{r+2},x_{r+1})
\prod_{r=a+3}^{n-1}k(x_r\mid x_{r-2},x_{r-1}),\tag{2}
\]

where

\[
f_3(x,y,z)=c_2c_1\exp[-6x^2-\tfrac{33}2y^2-6z^2
-12\sqrt2(xy+yz)-6xz].\tag{3}
\]

These are exact density identities, not independence assumptions.

## 3. Positive completion bounds

For nonnegative old states and a nonnegative new state,

\[
k(z\mid x,y)\le c_1e^{-6z^2},
\qquad \int_0^\infty k(z\mid x,y)\,dz\le\frac12.\tag{4}
\]

For two successive positive new states, dropping only nonnegative shifts
inside squares gives

\[
k(z\mid x,y)k(w\mid y,z)
\le c_1^2e^{-6z^2-6(w+\sqrt2z)^2}.
\]

Set `u=sqrt(6)z`, `v=sqrt(6)(w+sqrt(2)z)`. The Jacobian of this change is
six and the domain becomes `v>=sqrt(2)u>=0`. Consequently its integral is

\[
\kappa=\frac1\pi\int_{v\ge\sqrt2u\ge0}e^{-(u^2+v^2)}\,du\,dv
=\frac{\arctan(1/\sqrt2)}{2\pi}<\frac14.\tag{5}
\]

The last equality follows directly by polar coordinates: the wedge has angle
`arctan(1/sqrt(2))` and `integral_0^infinity r exp(-r^2)dr=1/2`.

Fixing either an anchor pair or an anchor triple at positive values, strip
each exterior arm from its far endpoint. Integrating its last two coordinates
conditionally on the preceding two costs at most `kappa` by (5); a remaining
single coordinate costs at most `1/2` by (4). This works identically for the
reverse kernels in (1) and (2). All integrands are nonnegative, so Tonelli
justifies the integration order. For an arm of length `r` its positive
completion mass is therefore at most

\[
F(r)=\kappa^{\lfloor r/2\rfloor}2^{-(r\bmod2)}.\tag{6}
\]

The product `F(a)F(n-a-3)` bounds the two exterior arms of a triple at
`a,a+1,a+2`. The bound is uniform in its nonnegative anchor values.

## 4. Triple outlier tails

For an endpoint outlier `x>=L`, use (3) and discard only negative exponent
terms which do not contain `x`. The `y,z` integrations then separate:

\[
\begin{aligned}
\int_L^\infty\!\int_0^\infty\!\int_0^\infty f_3(x,y,z)\,dz\,dy\,dx
&\le\frac{c_2c_1}{72\sqrt2}
\int_L^\infty\frac{e^{-6x^2}}{x^2}\,dx\\
&\le\frac{c_2c_1e^{-6L^2}}{864\sqrt2 L^3}=:U_e.
\end{aligned}\tag{7}
\]

Indeed the separate exponential integrals are `1/(12sqrt(2)x)` and
`1/(6x)`. Also `x^(-2)<=x/L^3` for `x>=L`, and
`integral_L^infinity x exp(-6x^2)dx=exp(-6L^2)/12`.
Symmetry under `x,z` proves the identical bound for the other endpoint.

For a middle outlier `y>=L`, the `x,z` integrations each give
`1/(12sqrt(2)y)` after dropping their other nonpositive exponent terms.
Therefore

\[
\begin{aligned}
\int_L^\infty\!\int_0^\infty\!\int_0^\infty f_3(x,y,z)\,dx\,dz\,dy
&\le\frac{c_2c_1}{288}
\int_L^\infty\frac{e^{-(33/2)y^2}}{y^2}\,dy\\
&\le\frac{c_2c_1e^{-(33/2)L^2}}{9504L^3}=:U_m.
\end{aligned}\tag{8}
\]

Here `integral_L^infinity y exp(-(33/2)y^2)dy=exp(-(33/2)L^2)/33`.

For `n>=3`, associate an outlier position `i` with the triple beginning at

\[
a_i=\begin{cases}
0,&i=0,\\
n-3,&i=n-1,\\
i-1,&1\le i\le n-2.
\end{cases}
\]

At the two path endpoints the outlier is an endpoint of its triple; every
other outlier is the middle of its triple. Applying (2), (6), (7), (8), and
then the union bound gives the following bound on the mass outside `[0,L]^n`
in one positive orthant:

\[
T_n=\sum_{i=0}^{n-1}U_iF(a_i)F(n-a_i-3),
\quad
U_i=\begin{cases}U_e,&i=0\text{ or }n-1,\\U_m,&\text{otherwise.}\end{cases}\tag{9}
\]

This covers `n=3` too, with both arm lengths zero. Overlapping outlier events
only make this union bound conservative. The probability `p_(n-1)` is twice
the one-orthant integral, so its error allowance is twice (9).

### The two-coordinate case

No tail or quadrature approximation is used for `p_1`. Write the original
unit Gaussian pair as `(U,rU+sigma V)`, where `U,V` are independent unit
normals, `r=2sqrt(2)/3`, and `sigma=1/3`. In the rotationally symmetric
`(U,V)` plane each of the two opposite-sign events is a wedge of angle
`arctan(sigma/r)`. It follows that

\[
p_1=\frac{\arctan(1/(2\sqrt2))}{\pi}.\tag{10}
\]

The initial two-coordinate tensor state is still used to evaluate the
higher-dimensional tensor rules. Equation (10) replaces only the first
probability in the final signed sum.

## 5. Localized Gauss--Legendre error

Let `Q` be the exact positive `N`-node Gauss--Legendre operator on `[0,L]`,
and let `I` be ordinary integration on this interval. Both give the constant
function one the value `L`. Define

\[
q_\mathrm{half}=Q[c_1e^{-6z^2}],
\qquad h=\max\{1/2,q_\mathrm{half}\}.\tag{11}
\]

For nonnegative old states, either `I` or `Q`, applied in the new state of
`k`, has value at most `h`: use (4) for `I`, and its pointwise first
inequality plus positivity for `Q`. The program encloses
`q_half` using its already certified exact nodes and weights. For the final
parameters its upper endpoint is below `1/2`, so `h=1/2` exactly. More generally the formula below remains valid with any certified rational
upper bound on h. The program uses max(1/2,the computed upper endpoint of
q_half); for these production parameters this is exactly 1/2.

Fix a coordinate `i` of the `n`-coordinate density and put

\[
s_i=\min(2,i)+\min(2,n-1-i).\tag{12}
\]

Use the pair anchor `a=i` when `i<n-1`, and `a=n-2` for the last coordinate.
In (1), group the pair density and precisely the transition factors which
contain `x_i`. Call their product `G_i`, and the remaining product `R_i`.
The local product contains the pair plus exactly `s_i-1` transitions; its
other variables are exactly the coordinates at distances one or two from
`i` that exist within the path. Thus there are `s_i` such variables, and

\[
f_n=G_iR_i,
\qquad
G_i=c_2c_1^{s_i-1}\exp[-H_i/2],\tag{13}
\]

where `R_i` is independent of `x_i`. The quadratic form `H_i` is
nonnegative on all real vectors: it is the pair's positive definite
quadratic form plus the squares from the selected transitions. This also
checks the endpoints and `n=2,3` cases; no interior-only estimate is used.

Let `q_i` denote the coefficient of `x_i^2` in `H_i`. It equals the
corresponding precision diagonal of the full density, because the factors
in `R_i` do not contain `x_i`. Expansion gives

\[
(q_i)=\begin{cases}
(9,9),&n=2,\\
(12,33,12),&n=3,\\
(12,36,\underbrace{39,\ldots,39}_{n-4},36,12),&n\ge4.
\end{cases}\tag{14}
\]

Fix all variables other than `t=x_i`. Completing the square writes

\[
G_i=a\exp[-q_i(t-\mu)^2/2],
\qquad 0<a\le c_2c_1^{s_i-1}.\tag{15}
\]

The amplitude bound follows since the residual exponent is the minimum of
the nonnegative quadratic form `H_i` over the real variable `t`.

For completeness the Gaussian derivative estimate used here is

\[
\left|\frac{d^{2N}}{dt^{2N}}e^{-q(t-\mu)^2/2}\right|
\le q^N\frac{(2N)!}{2^NN!}.\tag{16}
\]

An elementary proof is to let `Z` be a centered normal of variance `q`.
Gaussian integration by parts gives
`F'(s)=-qsF(s)`, `F(0)=1` for `F(s)=E exp(isZ)`, hence
`F(s)=exp(-qs^2/2)`. Repeated differentiation under the integral is
justified by Gaussian polynomial integrability and gives an upper bound
`E Z^(2N)`. Integration by parts gives the recurrence
`E Z^(2N)=q(2N-1) E Z^(2N-2)`, hence exactly the right side of (16).
Translation by `mu` does not change the estimate.

Combining (13)--(16) yields, pointwise on the remaining positive coordinates,

\[
|\partial_i^{2N}f_n|
\le c_2c_1^{s_i-1}q_i^N\frac{(2N)!}{2^NN!}\,R_i.\tag{17}
\]

Apply any mixture of `I` and `Q` to the `n-1` remaining coordinates.
Strip each exterior arm from its far endpoint. There are exactly
`n-1-s_i` remaining transition kernels in `R_i`, each removed at a cost
at most `h` by (11). The `s_i` neighbors are then free and each costs
`L`. Positivity and commutation of finite sums and finite integrals justify
these eliminations for every mixture of operators. Thus

\[
\left(\prod_{r\ne i}\mathcal O_r\right)R_i
\le L^{s_i}h^{n-1-s_i},\qquad \mathcal O_r\in\{I,Q\}.\tag{18}
\]

This is the key improvement over multiplying the uniform derivative bound
by `L^(n-1)`: the unaffected kernels remain present and provide positive
completion bounds.

The usual one-coordinate Gauss--Legendre remainder on `[0,L]` is bounded by

\[
\frac{L^{2N+1}(N!)^4}{(2N+1)((2N)!)^3}\sup|g^{(2N)}|.
\]

Telescope `I^n-Q^n` one coordinate at a time. Apply this one-coordinate
bound with (17), and then the other positive operators using (18). Summing
the resulting `n` bounds proves

\[
\boxed{
E_{n,N}=
\frac{c_2L^{2N+1}(N!)^3}{(2N+1)2^N((2N)!)^2}
\sum_{i=0}^{n-1}q_i^N c_1^{s_i-1}L^{s_i}h^{n-1-s_i}.
}\tag{19}
\]

This holds for every `n>=2,N>=1,L>0`. At `n=2` it equals the retained
Gaussian-moment bound; nevertheless (10) is used instead of numerical
quadrature for `p_1`. In the final program every term in (19) is evaluated
by positive directed interval arithmetic.

## 6. Exact state rescaling and integer arithmetic

For exact Gauss nodes and weights let `M_j` be the usual two-coordinate
tensor state for the `j`-edge path. The unchanged factorization is

\[
w_k k(x_k\mid x_i,x_j)=A_{ij}B_{ik}C_{jk},
\]
\[
A_{ij}=e^{-3x_i^2/2-6\sqrt2x_ix_j},\quad
B_{ik}=e^{-6x_ix_k},\quad
C_{jk}=c_1w_ke^{-12x_j^2-6x_k^2-12\sqrt2x_jx_k}.
\]

The new state is the exact rescaling

\[
\widetilde M_j=4^{j-1}M_j.
\]

Thus its initial state is unchanged and its update is

\[
\widetilde M_{j+1}[j,k]
=4C_{jk}\sum_i\widetilde M_j[i,j]A_{ij}B_{ik}.\tag{20}
\]

The tensor integral is recovered as
`sum(tilde M_j)/4^(j-1)`. Equation (20) changes neither the quadrature
rule nor its analytic error. It keeps the state well resolved on a fixed
grid as the unscaled positive orthants shrink.

Let `S=2^128`. A nonnegative Decimal interval `[a,b]` is converted outward
to integers `[floor(Sa),ceil(Sb)]` using the exact fractions returned by
`Decimal.as_integer_ratio()`. Assume the integer pairs for the state and
factors enclose their exact nonnegative values divided by `S`.
First compute

\[
t^-_{ij}=\left\lfloor m^-_{ij}a^-_{ij}/S\right\rfloor,
\qquad
t^+_{ij}=\left\lceil m^+_{ij}a^+_{ij}/S\right\rceil.
\]

Then form exact nonnegative integer dot products

\[
d^-_{jk}=\sum_i t^-_{ij}b^-_{ik},\qquad
d^+_{jk}=\sum_i t^+_{ij}b^+_{ik}.
\]

The new endpoints are

\[
m'^-_{jk}=\left\lfloor 4c^-_{jk}d^-_{jk}/S^2\right\rfloor,
\qquad
m'^+_{jk}=\left\lceil 4c^+_{jk}d^+_{jk}/S^2\right\rceil.\tag{21}
\]

Before the last division, the represented products have denominator `S^3`;
division of the numerators by `S^2` restores denominator `S` and encloses
the rescaled update (20). Positivity proves this inductively at every step.
Python integer addition and multiplication are exact. `sum(map(mul,...))`
computes these same exact finite sums; it uses no floating point or matrix
library. For a nonnegative integer `u`, `u>>128` is exactly `floor(u/S)`
and `(u+S-1)>>128` is exactly `ceil(u/S)`. The corresponding 256-bit shifts
implement the divisions by `S^2` in (21).

The state sum is converted to Decimal using explicitly directed division
by `S`, then divided outward by the exact integer `4^(j-1)`. There is no
assumption that rescaling errors cancel: every rounding error is enclosed.

The exact tensor probability for `j>=2` is enclosed by

\[
2\bigl([s^-,s^+]+[-E_{j+1,N},E_{j+1,N}]+[0,T_{j+1}]\bigr).\tag{22}
\]

The signed sum reverses endpoints on its negative terms. For `j=1` use
(10) directly.

## 7. Nodes, weights, arithmetic and elementary function certificates

Pi, root brackets, quadrature weights, density factors, analytic error bounds and the final alternating sum use Python's standard-library `decimal` with 100 decimal digits. Every Decimal lower endpoint uses `ROUND_FLOOR` and every upper endpoint uses `ROUND_CEILING`. Repeated positive tensor contractions use exact integers on the fixed-point grid proved in Section 6, with outward integer floor and ceiling. Square root and exponential are correctly rounded in Decimal; moving their rounded values one representable number outward encloses the exact transcendental value. Pi is enclosed by Machin's identity and 100-term alternating arctangent series, including an explicit next-term remainder.

Approximate Legendre roots are only starting guesses. The program brackets each root by rational decimal endpoints, evaluates the degree-`N` Legendre polynomial at both endpoints by interval recurrence, and verifies strictly opposite signs. The `N` brackets are disjoint, so they isolate all roots of a degree-`N` polynomial. Legendre's derivative is enclosed on each entire root bracket and the exact weight formula

\[
w=\frac2{(1-x^2)(P_N'(x))^2}
\]

is evaluated on that interval. Roots and weights are then mapped to `[0,L]`.

The tensor quadrature is evaluated by the exact Markov factorization, without changing its mathematical value. For the unscaled tensor state, if `M[i,j]` already includes the quadrature weights in both live coordinates, initialize

\[
M[i,j]=w_iw_j f_2(x_i,x_j)
\]

and update

\[
M_{new}[j,k]=w_k\sum_i M[i,j]k(x_k\mid x_i,x_j).
\]

Every term is nonnegative. The exact fixed-point enclosure proved in Section 6 therefore gives lower and upper bounds on each contraction. For the production rescaled state, divide the final matrix sum by 4^(j-1), as proved in Section 6, to enclose the complete tensor quadrature. Subtracting/adding `E_{n,N}`, adding the one-sided tail `T_n`, and multiplying by two encloses `p_j` for j>=2; equation (10) supplies p_1 exactly. The finite signed sum is then evaluated with interval arithmetic, including the correct endpoint reversal for its negative terms.

For speed, the transition density is evaluated by the exact factorization

\[
k(z\mid x,y)=\sqrt{6/\pi}e^{-3x^2/2}e^{-12y^2}e^{-6z^2}
 e^{-6\sqrt2xy}e^{-6xz}e^{-12\sqrt2yz}.
\]

All exponential factors depend on at most two coordinates, so only O(N^2) exponential evaluations are required. The program precomputes the positive interval factors and combines each triple with directed multiplications. This is an algebraic identity, with every multiplication error still enclosed; it does not change the quadrature rule or its analytic error bound.

The implementation does not assume correctly rounded platform binary `libm` transcendental functions: the initial binary cosine guess is irrelevant after the root-bracket sign checks succeed.



The root and weight checks are retained from the frozen source: binary
cosines supply only starting guesses, Newton iterations choose centers,
and the program verifies strict opposite polynomial signs at both ends
of every bracket, disjointness of all `N` brackets, and nonvanishing of
the derivative enclosure. Since the Legendre polynomial has degree `N`,
these are all its roots. The exact positive weight formula is evaluated
outward on the brackets. The final bracket radius is the exact decimal
`10^-65`; initialization precision remains 100 digits. Failure of any
required inequality aborts the run.

Pi is enclosed by the same 100-term Machin-series construction. To enclose
`atan(x)` at `x=1/sqrt(2)` and `x=1/(2sqrt(2))`, use respectively 220 and
120 terms of

\[
\arctan x=\sum_{k=0}^{m-1}\frac{(-1)^kx^{2k+1}}{2k+1}+R_m,
\quad |R_m|\le\frac{x^{2m+1}}{2m+1},\qquad 0<x<1.
\]

The elementary alternating-series remainder follows by integrating the
finite geometric expansion of `1/(1+x^2)`. All series operations and the
next-term radii are directed intervals. Correctly rounded Decimal square
root and exponential evaluations are expanded one representable number
outward. The remaining arithmetic uses only exact integers and directed
Decimal intervals. No experimental agreement is a proof premise.

## 8. Production witness, reproducibility and attribution

The production parameters are N=80, J=65, L=3, 100-digit Decimal arithmetic,
fixed-point scale 2^128 and exact state rescaling 4^(j-1). The exact code is
reproduced from the submission directory with

    python3 -I -B certify.py < /dev/null

The sole declared JSON input is {}. The program uses only Python's standard
library and no external data or provider credentials. Its source SHA-256 is

    971824da89474b02868a113ae49191a436f51975b450388b44dd90f39032cf63

A resource-limited run of these exact bytes under Python 3.12.14, an empty
environment and input {} completed all 65 terms in 5.176 seconds, with
22,528 KiB peak RSS and exit zero. It imposed 20-second wall and CPU limits,
256 MiB address space and 32 processes. This is local host evidence; Docker
is unavailable locally, so the required hosted computation and mathematical
review must be run afresh. No previous judgment is a premise.

The complete computed interval is

    lower = 0.1306475664231676643417493777627387226446737087839483086429013850383116864845439257011675727543687882
    upper = 0.1306475664231676643417493806534603210307271447454970258635182301766950514285992790646531631160644148

Consequently the shorter displayed bounds are outward:

    0.1306475664231676643417493777 <= Gamma_65
       <= 0.1306475664231676643417493807
        < 0.130647566423167664341749381
        = 130647566423167664341749381/1000000000000000000000000000
        = gamma_0.

The strict margin from the complete upper endpoint to gamma_0 is greater
than 3.465e-28. The code exits zero only after this exact directed comparison.
It prints all 65 probability, quadrature-error, tail-error and signed-sum
intervals, plus q_half, the positive operator bound and completion constant.
The computed q_half upper endpoint is strictly below 1/2, justifying h=1/2.

The graph proof fixes J=65 before taking any limit. Its deletion coefficient is

    C_65 = 1+3*65*2^64 = 3597115094373362565121.

This is a large finite constant. It changes the unspecified profile and length
thresholds, while preserving the same uniform little-o graph quantifier.
The graph constant gamma_0 is internal: proof.md Section 5 uses its strictly
larger intermediate circuit coefficient to absorb the sublinear loss before
claiming the eventual coefficient

    315323783211583832170874691/65323783211583832170874691.

The immediate source is submission 3f2be7ea-aa22-4f7e-98e9-23b40557dbcf,
commit 737818f9911e00e9d4909c644c771e0724fb6c37. It supplied the complete
mathematical package and positive-slack eventual-bound assembly. Submission
0514f3af-5d65-4556-82df-a89e2a3ae4e1 supplied the earlier Gaussian-moment
and positive-completion error analysis with J=39; submission
ac0dd6e1-8a73-4682-b7aa-67502b2de08b supplied the exact integer contractions;
submission 38e961ae-bd55-4566-85ad-17112333f84a supplied the factored kernel
and parameter refinement, with its original hosted Decimal timeout disclosed.
The original complete mathematical package is fcf21648-1b1f-463c-b5bc-c50f53f63436,
commit 81470c9513639af985c88699dbcd2e2efbef7463. All five are declared
for attribution; their judgments are not reused.

This revision's new ingredients are the localized quadrature bound (19),
the reversible triple-outlier tails (9) with sharp two-step completion factors,
the exact first probability (10), and the rigorously rescaled contractions
(20)--(22). These support the new finite certificate. The complete graph and
circuit reductions remain necessary and are supplied in the other evidence.
