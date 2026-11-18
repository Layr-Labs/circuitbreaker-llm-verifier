# Finite path and buffered local Gaussian certificates

Sections 1--8 prove the finite path premise used by the baseline graph
argument. Sections 9--10 prove the quantitative local Gaussian estimates
used by the structural improvement in cutwidth.md Section 11. The analytic
reductions, quadrature validation, truncation errors, probability estimates,
covariance perturbation arguments and arithmetic enclosures are all supplied. The standard-library computation certify.py is an
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

The path argument uses only this finite sum; no infinite-series limit is used.


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
\Gamma_{65}<\gamma_base,
\qquad
\gamma_base=\frac{130647566423167664341749381}{10^{27}}
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

The path parameters remain N=80, J=65, L=3, 100-digit Decimal arithmetic,
fixed-point scale 2^128 and exact state rescaling 4^(j-1). Reproduce the
single combined certificate from the submission directory with

    python3 -I -B certify.py < /dev/null

The sole declared input is {}. The program uses only Python's standard
library and no external data or provider credentials. Its main routine
first checks the exact local-event and structural constants of Sections
9--10, then encloses all 65 path terms. All results are included in one
JSON object. The output's `target` is gamma_base, the finite path target;
`structural_certificate.gamma_new` is the strictly smaller graph constant.
The code does not compare Gamma_65 directly to gamma_new.

The path computation gives the complete interval

    lower = 0.1306475664231676643417493777627387226446737087839483086429013850383116864845439257011675727543687882
    upper = 0.1306475664231676643417493806534603210307271447454970258635182301766950514285992790646531631160644148

Consequently these shorter displayed bounds are outward:

    0.1306475664231676643417493777 <= Gamma_65
       <= 0.1306475664231676643417493807
        < 0.130647566423167664341749381
        = gamma_base.

The strict margin from the complete upper endpoint to gamma_base is greater
than 3.465e-28. Successful exit requires that directed comparison, all new
exact local-event assertions, and the exact structural and circuit-constant
comparisons. The program prints every path probability, quadrature-error,
tail-error and signed-sum interval, plus the positive quadrature/operator
and completion bounds. The q_half upper endpoint is below 1/2, allowing
h=1/2 in the analytic quadrature estimate.

The finite deletion coefficient is C_65=1+3*65*2^64=
3597115094373362565121. The graph proof fixes it before every profile or
graph-size limit. Its new Section 11 proves the additional uniform gain
delta_g=10^-18, giving

    gamma_new = gamma_base-delta_g
              = 130647566423167663341749381/1000000000000000000000000000.

Theorem 2 of cutwidth.md supplies H(gamma_new,h_new); proof.md Section 5
transfers it and leaves fixed positive slack to obtain the eventual
coefficient 315323783211583831670874691/65323783211583831670874691.
The finite arithmetic below supports those reductions, and does not
replace their supplied proofs. Execution records and immutable source
digests are kept outside the submission. Local Docker execution is
unavailable; hosted execution and mathematical evaluation must be fresh.

The immediate source is ae372f38-09d4-452a-af05-4954cdd2f168, commit
2a8b03fe0be18fdf4918dc5fc34da51534f1aca1. It supplied the complete
Sections 1--7 path proof and the retained 65-term certificate. The prior
numerical and mathematical packages are 3f2be7ea-aa22-4f7e-98e9-23b40557dbcf,
0514f3af-5d65-4556-82df-a89e2a3ae4e1,
ac0dd6e1-8a73-4682-b7aa-67502b2de08b,
38e961ae-bd55-4566-85ad-17112333f84a (its hosted Decimal timeout is
disclosed), and fcf21648-1b1f-463c-b5bc-c50f53f63436. All six are declared
for attribution, with exact source commits in proof.md. No prior judgment
is reused. Sections 9--10 and their corresponding exact arithmetic are
the present revision's new numerical contribution.

## 9. Quantitative buffered-curvature certificate

Set

\[
 q=\frac3{100000},\quad \theta=\frac1{100},\quad
 \zeta_c=\frac1{700000},\quad \eta=\frac{99}{70000000}.
\]

The curvature window has adjacent roots `u,v`, their four other neighbors `A`, the two children `B,C` of each `A`, and the two children `D,E` of each `B`: twenty non-root coordinates. At each root designate one positive and one negative `A`; prescribe `B` opposite its parent, `C` equal to its parent, and `D,E` equal to `B`. Call this event `E`. Its buffer isolates the relevant crossing components from arbitrary outside anchors. Their cover number is four for equal root bits and five for unequal root bits, so the mixed difference is `-2` on `E`.

**Quantitative lemma.** Let `C0` be the critical covariance on this window, with `w_d=(1+d/3)2^(-d/2)`. If `C` is the covariance of any centered unit-variance Gaussian vector satisfying `max|C-C0|<=zeta_c` entrywise, then its root pair is nonsingular and

\[
 \mathbb P_C(E\text{ relative to }t\mid X_u=X_v=t)>q/2
 \quad (|t|\leq\theta). \tag{C1}
\]

The full covariance may be singular. The critical probability at zero is separately greater than `1/20000>q`.

### 9.1 Exact critical probability lower bound

Write `sigma^2=1/12`, `r=2sqrt(2)/3`, and `phi(x)=sqrt(6/pi)exp(-6x^2)`. On the critical tree, a vertex `z` with parent `p` has children

\[
 X_{z_1}=\sqrt2X_z-X_p/2+\xi_z,
 \qquad X_{z_2}=\sqrt2X_z-X_p/2-\xi_z.
\]

The ten innovations at the roots and the `A,B` vertices are independent `N(0,1/12)`, independent of the root pair. Indeed the covariance values `w_1=r,w_2=5/6` make the variance of `sum_neighbors X-2sqrt(2)X` zero. The half-difference of two children has variance `(2-2w_2)/4=1/12` and covariance zero with all previously exposed coordinates, whose distances to the two children agree. Exposing pairs in parent-before-child order proves the claimed Gaussian independence.

Conditioning both roots to zero, the two `A` values at each root are `a,-a`. For `A=s*a`, writing `B=-s*b` gives innovation density `phi(b+sqrt(2)*a)` and `C=s(2sqrt(2)*a+b)`. The two children of `B` are `-s(sqrt(2)*b+a/2)±xi_B`. Consequently, if `E_M` requires magnitude at least `M` at every prescribed coordinate, then

\[
 \mathbb P(E_M\mid X_u=X_v=0)=I_M^2,
\]

\[
 I_M=\int_M^\infty\phi(a)
 \left[\int_M^\infty\phi(b+\sqrt2a)
 G(\sqrt2b+a/2-M)\,db\right]^2da, \tag{C2}
\]

where `G(h)=P(|N(0,1/12)|<=h)`. Both squares come from disjoint independent innovations.

Here is a fully specified rational lower sum for (C2). Put

\[
 c=69/50,\quad s_-=707/500,\quad s_+=283/200,
 \quad L(x)=(1-x/65536)^{65536}.
\]

Squaring proves `s_-<sqrt(2)<s_+`. The identity

\[
 22/7-\pi=\int_0^1x^4(1-x)^4/(1+x^2)\,dx>0
\]

and `c^2<21/11` imply `c<sqrt(6/pi)`. Also `L(x)<=exp(-x)` for `0<=x<65536`. For `h>=0`, let `z=min(h,1/2)` and put

\[
 g(h)=\frac{69}{25}
 \left(z-2z^3+\frac{18}{5}z^5-\frac{36}{7}z^7\right).
\]

Taylor's formula gives `exp(-y)>=1-y+y^2/2-y^3/6` for `y>=0`. Integrating it at `y=6x^2` proves `G(h)>=g(h)`. The polynomial is nonnegative and increasing on `[0,1/2]`: its derivative is the decreasing cubic `1-y+y^2/2-y^3/6` on `0<=y<=3/2`, whose minimum is `1/16`.

For `M=1/20`, take `Delta=1/400`, `a_j=M+j*Delta` for `0<=j<180`, and `b_k=M+k*Delta` for `0<=k<300`. Write `A_j=a_j+Delta`, `B_k=b_k+Delta`. Restricting (C2) to `a<=1/2,b<=4/5` gives

\[
 I_M\geq\Delta\sum_j cL(6A_j^2)
 \left[\Delta\sum_k
 cL\bigl(6(B_k+s_+A_j)^2\bigr)
 g(s_-b_k+a_j/2-M)\right]^2. \tag{C3}
\]

The densities decrease on the positive half-line and `G` increases. The upper `a` endpoint in the density and the lower `a` endpoint in the buffer therefore give a nonnegative inner bound uniform throughout each `a` cell, which can be squared. All arguments of `L` are below fourteen.

The right side is an explicit rational number. Downward integer arithmetic at denominator `10^30` gives

\[
 I_{1/20}\geq
 \frac{4024054722670646630717842689}{10^{30}}>\frac1{250}.
 \tag{C4}
\]

For `M=0,Delta=1/40`, twenty `a` cells and thirty-two `b` cells give

\[
 I_0\geq\frac{7239337205339735179200724088}{10^{30}},
 \qquad I_0^2>\frac1{20000}. \tag{C5}
\]

The function `certify_buffered_curvature()` checks these finite inequalities using only integers and rational numbers. It lower-bounds `L(x)` by rounding `1-x/65536` down, then squaring sixteen times with each product rounded down. It evaluates `g` as an exact rational before rounding down and rounds products and sums in (C3) down. All factors are nonnegative, so every rounding preserves the lower-bound direction. Thus (C2)–(C3) provide the analytic proof of the probability bound; the computation only verifies the displayed finite rational inequalities.

### 9.2 Conditional perturbations, including singular covariances

Temporarily let `0<zeta<=10^-5` and `max|C-C0|<=zeta`. The root covariance is `K=[[1,c_R],[c_R,1]]`, with `|c_R-r|<=zeta`, hence `0<c_R<19/20` and `||K^-1||_op<20`. The critical regression rows onto the roots are

\[
 (\sqrt2,-1/2),\quad (3/2,-\sqrt2/2),\quad
 (\sqrt2,-3/4), \tag{C6}
\]

at the three successive depths, or their swaps. Their sums `s_i` lie in `[3/5,1]`. Under conditioning `R=(X_u,X_v)=(t,t)`, the mean is

\[
 \mu_i=t(C_{iu}+C_{iv})/(1+c_R).
\]

Its coefficient differs from `s_i` by at most `(2+s_i)zeta/(1+c_R)<=3zeta`. Therefore

\[
 |\mu_i-t|\leq|t|/2. \tag{C7}
\]

Extract ten linear coordinates `I_z=(X_(z1)-X_(z2))/2` from the perturbed vector. Each extraction row has coefficient absolute sum one. At criticality their covariance is `sigma^2 I_10` and they are uncorrelated with `R`; under `C` those entries change by at most `zeta`. Writing `V=Cov(I,R)`, the conditional innovation covariance is `U=Cov(I)-VK^-1V^T`. Each entry of the subtracted term has magnitude at most `40zeta^2`. Hence

\[
 \|U-\sigma^2I_{10}\|_{\rm op}\leq10(\zeta+40\zeta^2)
 \leq(101/10)\zeta,
 \quad \lambda_{\min}(U)>(7/25)^2.
\]

Diagonalizing `U` and using `|sqrt(x)-sigma|=|x-sigma^2|/(sqrt(x)+sigma)` gives

\[
 \|\sqrt U-\sigma I_{10}\|_{\rm op}
 <(101/10)\zeta/(14/25)<19\zeta. \tag{C8}
\]

Let `W` denote the twenty non-root coordinates. At criticality `W=beta R+L I`. At the three depths, the nonzero absolute coefficients of a row of `L` are respectively `(1)`, `(sqrt(2),1)`, and `(3/2,sqrt(2),1)`. Thus `||L_i||_2^2<=21/4`.

For the perturbed law define `D_i=W_i-beta_i R-L_i I`. It vanishes identically under the critical law. As a linear form in the original twenty-two coordinates, its coefficient absolute sum is at most `3/2+sqrt(2)`, `5/2+3sqrt(2)/2`, or `13/4+2sqrt(2)` according to depth. To check the last bound, the coordinate and its last innovation combine to the average of its sibling pair, contributing one; earlier innovations contribute `3/2+sqrt(2)` and the root regression contributes `sqrt(2)+3/4`. The largest bound has square less than 37. Therefore `Var_C(D_i)<=37zeta`, since the critical variance is zero and each covariance entry changes by at most `zeta`. Gaussian conditioning only decreases this variance.

Under the conditional law set `g=U^-1/2(I-E[I|R])` and `W0=L sigma g`. Then `g` is standard Gaussian and `W0` has exactly the critical conditional law with both roots zero. Moreover

\[
 W-\mu-W^0=L(\sqrt U-\sigma I)g+(D-\mathbb E[D\mid R]).
\]

These terms need not be independent. The triangle inequality in `L^2`, (C8), and the residual variance bound imply, for each coordinate,

\[
 \operatorname{Var}(W_i-\mu_i-W_i^0)
 \leq\left(\sqrt{37\zeta}+\sqrt{1896}\,\zeta\right)^2
 \leq40\zeta. \tag{C9}
\]

The last inequality follows from `sqrt(37)<61/10`, `sqrt(1896)<44`, `sqrt(10^-5)<1/316`, and `2(61/10)44/316+1896/100000<3`. This constructs the coupling directly; no nonsingularity assumption on the other twenty coordinates or general singular square-root continuity theorem is used.

Take `zeta=zeta_c`. On `E_(1/20)`, (C7) leaves error allowance `h=1/20-theta/2=9/200`. A centered Gaussian of variance at most `v` satisfies `P(|Y|>=h)<=2exp(-h^2/(2v))`, by its moment generating function and Markov's inequality optimized at `h/v`. A union bound over (C9) gives failure probability at most

\[
 40e^{-h^2/(80\zeta_c)}=40e^{-567/32}<10^{-6}.
\]

For the final strict bound, `exp(1/2)>79/48`, `exp(7/32)>39/32`, and the exact integer inequality `39*79^35>40000000*32*48^35` imply `exp(567/32)>40000000`. Equations (C4) and the union bound now give

\[
 \mathbb P_C(E\mid R=(t,t))>
 \frac1{62500}-\frac1{1000000}
 =\frac3{200000}=q/2,
\]

proving (C1) without independence between the buffered event and the errors.

### 9.3 Interpolation and the constants used in the graph proof

For the nonnegative profile fields, cover comparison gives `0<=C_Y-C_X<=1` entrywise. Fix the initial local critical-window error at most `zeta_c/100`. Along `C_rho=C_X+rho(C_Y-C_X)`, the whole interval `0<=rho<=eta=99zeta_c/100` stays within `zeta_c`, so (C1) applies. The conditional mixed difference is at most `-q` and the adjacent pair density at `(t,t)` is at least

\[
 d_0=e^{-\theta^2}/(2\pi)>
 (1-1/10000)7/44>159/1000.
\]

Thus the strict covariance-comparison constant is

\[
 A=\eta qd_0>\frac{47223}{7000000000000000}.
\]

For the threshold inequality, `b0=(1-r)/2>1/35` follows from `sqrt(2)<99/70`, and `exp(-theta^2/2)>999/1000`. Hence its constant is

\[
 B_*=(3/2)qb_0d_0e^{-\theta^2/2}\theta^2
 >\frac{1429569}{70000000000000000}.
\]

With the graph proof's `alpha=1/6370000`, exact rational arithmetic gives

\[
 A\alpha>\frac{47223}{44590000000000000000000}>10^{-18},
 \qquad B_*>10^{-18}.
\]

These constants are fixed before the profile radius, positive coordinate noise, or input graph. All finite arithmetic in this subsection is checked by `certify_buffered_curvature()`.

## 10. Quantitative buffered majority-gain probability

The following finite Gaussian estimate supplies the gain-event constants
used in the graph argument:

\[
 p=10^{-17},\qquad \theta_g=\frac1{100},\qquad
 \zeta_g=\frac1{35000}. \tag{G1}
\]

Root the cubic tree at an edge `vw`. Let `a,b` be the other neighbors of
`v`, and `c,d` the other neighbors of `w`. Include both remaining neighbors
of each of `a,b,c,d`. These fourteen vertices form a tree window of
diameter five. The required pattern has positive coordinates at `v,w,d`
and the two children of `d`, and negative coordinates at `a,b,c` and their
six children. Signs at threshold `t` mean `b_t^+(z)=1[X_z>t]`.

**Claim.** Every centered Gaussian vector with unit variances and covariance
entrywise within `zeta_g` of this critical window's covariance has the
required pattern with probability at least `p/2`, for every
`|t|<=theta_g`. The covariance is allowed to be singular.

### 10.1 Critical innovation box

The critical tree covariance is
`w_j=(1+j/3)2^(-j/2)`. Write `r=w_1=2sqrt(2)/3`. The root pair has density

\[
 f(v,w)=\frac3{2\pi}
       \exp[-v^2/2-(9/2)(w-rv)^2]. \tag{G2}
\]

For each `z=v,w,a,b,c,d`, its children away from the root edge satisfy

\[
 X_{z,1}=\sqrt2X_z-X_{\rm parent}/2+\xi_z,\qquad
 X_{z,2}=\sqrt2X_z-X_{\rm parent}/2-\xi_z, \tag{G3}
\]

where the six innovations are independent `N(0,1/12)` variables and are
independent of the root pair. For `v,w`, their parent in (G3) is the other
root. To verify this representation, the critical covariance gives the
variance-zero relation `sum_neighbors X=2sqrt(2)X`. The half-difference of
two siblings has variance `(2-2w_2)/4=1/12`. Those siblings have equal
distances to every previously exposed vertex in a parent-before-child
exposure, so the difference is Gaussian-orthogonal to the past. This proves
both (G3) and independence. Each innovation has density
`sqrt(6/pi) exp(-6x^2)`.

Use the following product box, labeling `c` as the child receiving
`+xi_w`:

| Variable | Interval |
|---|---|
| `X_v` | `[1/5,21/100]` |
| `X_w` | `[103/100,53/50]` |
| `xi_v` | `[-1/100,1/100]` |
| `xi_w` | `[-33/20,-8/5]` |
| `xi_a,xi_b` | `[-9/50,9/50]` each |
| `xi_c,xi_d` | `[-7/25,7/25]` each |

The bounds `7/5<sqrt(2)<99/70`, verified by squaring, give

\[
 X_a,X_b\le-26/125,\qquad X_c\le-703/3500,\qquad
 X_d\ge2937/1000. \tag{G4}
\]

For example the three bounds are obtained by substituting, respectively,
`(sqrt(2),X_v,X_w,xi_v)=(99/70,21/100,103/100,1/100)`,
`(sqrt(2),X_w,X_v,xi_w)=(99/70,53/50,1/5,-8/5)`, and
`(sqrt(2),X_w,X_v,xi_w)=(7/5,103/100,21/100,-8/5)` into (G3).
Applying (G3) once more gives, for either child index `i`,

\[
 \begin{aligned}
 X_{a,i},X_{b,i}&\le-(7/5)(26/125)-1/10+9/50=-132/625,\\
 X_{c,i}&\le-(7/5)(703/3500)-103/200+7/25=-2581/5000,\\
 X_{d,i}&\ge(7/5)(2937/1000)-53/100-7/25=16509/5000.
 \end{aligned} \tag{G5}
\]

Together with the two root intervals, these prove every required sign with
absolute margin at least `1/5`.

The box probability has an elementary explicit lower bound. Integrating
`1/(1+x^2)<=1-x^2/2` on `[0,1]` gives `pi<10/3`, hence
`3/(2pi)>9/20` and `sqrt(6/pi)>4/3`. On the root rectangle,

\[
 v^2/2+(9/2)(w-rv)^2
 \le441/20000+(9/2)(131/150)^2
 =13817/4000<7/2.
\]

Here `r>14/15`, and `w-rv` is positive. Including the six innovation
endpoint exponents gives total exponent less than

\[
 \frac72+6\left[(1/100)^2+(33/20)^2
                    +2(9/50)^2+2(7/25)^2\right]
 =52913/2500<43/2.
\]

Multiplication of all interval lengths and density prefactors gives

\[
 P_0=\frac1{100}\frac3{100}\frac9{20}\frac1{50}\frac1{20}
         \left(\frac9{25}\right)^2
         \left(\frac{14}{25}\right)^2\left(\frac43\right)^6
     =\frac{4704}{152587890625}.
\]

Thus the box probability is greater than `P_0 exp(-43/2)`. The exponential
series gives `8/3<e<87/32<11/4`; for the upper bound, bound the tail from
degree four by `(1/24)/(1-1/5)=5/96`. Also `sqrt(e)<5/3`. Exact rational
comparisons therefore yield

\[
 e^{43/2}<(11/4)^{21}(5/3)<3\cdot10^9,
 \qquad P_0/(3\cdot10^9)>10^{-17}=p. \tag{G6}
\]

### 10.2 Stability under arbitrary positive semidefinite perturbations

If two positive semidefinite `m`-dimensional covariances have entrywise
difference at most `zeta`, their operator norm difference is at most
`m zeta`, by Cauchy--Schwarz. Moreover, for positive semidefinite `A,B`,

\[
 \|\sqrt A-\sqrt B\|_{\rm op}
 \le\sqrt{\|A-B\|_{\rm op}}. \tag{G7}
\]

Indeed, with `delta=||A-B||op`, `A<=B+delta I` implies
`sqrt(A)<=sqrt(B+delta I)<=sqrt(B)+sqrt(delta) I`; interchange `A,B` for
the opposite inequality. The operator monotonicity used here follows from

\[
 \sqrt A=\pi^{-1}\int_0^\infty A(A+sI)^{-1}s^{-1/2}\,ds:
\]

inversion reverses the positive definite order, so the integrand
`I-s(A+sI)^(-1)` preserves order. The integral follows by diagonalization
and the scalar substitution `s=lambda u^2`; singular matrices are included.

Couple `X=C_0^(1/2)G` and `Y=C^(1/2)G` with one standard normal vector.
By (G7), each coordinate error has variance at most `m zeta`. The Gaussian
Chernoff bound and a union bound give

\[
 \Pr\{\max_i|Y_i-X_i|\ge h\}
 \le2m\exp[-h^2/(2m\zeta)]. \tag{G8}
\]

For `m=14`, `zeta=1/35000`, and `h=19/100`, this is at most

\[
 28e^{-361/8}<28(3/8)^{45}<\frac1{2\cdot10^{17}}=p/2, \tag{G9}
\]

using `e>8/3` and an exact rational comparison. On the critical box, a
coupling error less than `19/100` leaves positive coordinates strictly
above `1/100` and negative coordinates strictly below `-1/100`. They have
the required signs for every threshold in `[-1/100,1/100]`. Subtract (G9)
from (G6) to prove the claim.

### 10.3 Combined rational structural contract

The graph proof's diameter-five estimate and oriented-window count give
`tau=zeta_g/21=1/735000` and `alpha=3tau/26=1/6370000`, with finite-profile
error at most `zeta_g/21`. Its combined covariance-and-gain loss is
`min(A alpha,p/8)` at central thresholds. The companion curvature proof
uses initial covariance error at most `zeta_c/100`, where
`zeta_c=1/700000`. Since `0<=(C_Y-C_X)_{ij}<=1`, it permits
`eta=99zeta_c/100=99/70000000`. With its proved `q=3/100000` and
`d0>159/1000`, the exact comparisons are

\[
 \begin{aligned}
 A&>47223/7000000000000000,\\
 A\alpha&>47223/44590000000000000000000>10^{-18},\\
 p/8&=1/(8\cdot10^{17})>10^{-18},\\
 B_*&>1429569/70000000000000000>10^{-18}.
 \end{aligned}
\]

Thus `epsilon0=10^-18` is a valid common structural loss. The function
`certify_buffered_gain()` checks all rational arithmetic above, including
the sign margins, probability bounds, coupling tail, and combined
constant. Its assertions supplement the analytic proofs given here.
