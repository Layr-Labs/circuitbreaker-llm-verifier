# Certificate for the alternating-run constant

The declared standard-library Python experiment `certify.py` proves

\[
\gamma_5=p_1-p_2+p_3-p_4+p_5<\frac{1341}{10000}=0.1341. \tag{C}
\]

Here \(p_k\) is the probability of \(k\) consecutive sign transitions in the centered, unit-variance Gaussian sequence whose correlations are

\[
c_0=1,\qquad c_d=\left(\frac14+\frac1{\sqrt2}+\frac d4\right)2^{-d/2}\quad(d\ge1).
\]

This is precisely the finite numerical claim used by the graph proof. Set \((A_n)_{ij}=(-1)^{i-j}c_{|i-j|}\) for \(1\le i,j\le n\), and let \(Q_n(A)\) denote the positive-orthant probability of a centered Gaussian with covariance \(A\). Central symmetry gives \(p_k=2Q_{k+1}(A_{k+1})\). Interval LDL decomposition verifies positive definiteness of \(A_6\), hence of its leading principal submatrices.

The certificate does not directly integrate a six-dimensional orthant. It proves \(\gamma_5\le U<0.1341\), where \(U\) is an explicit upper-bound expression involving only five one-dimensional integrals. The output field `gamma5_majorant` encloses \(U\), **not** \(\gamma_5\). Its lower endpoint is not asserted to be a lower bound on \(\gamma_5\).

## 1. A quadratic bound for the final Gaussian coordinate

Write \(A=A_5\) and partition \(A_6\) as
\(\begin{pmatrix}A&r^{\mathsf T}\\r&1\end{pmatrix}\), with \(r\) a row vector. Put \(\beta=rA^{-1}\) and \(\sigma^2=1-rA^{-1}r^{\mathsf T}>0\). If \(Y\sim N(0,A)\), the conditional distribution of the sixth coordinate, given \(Y\), has mean \(\beta Y\) and variance \(\sigma^2\). For \(\mathcal O=\{Y_i>0\text{ for every }i\}\), put \(q=\Pr(\mathcal O)\) and \(z=\beta Y/\sigma\). Consequently

\[
p_5=2\mathbb E[1_{\mathcal O}\Phi(z)]. \tag{1}
\]

Let \(a=-1/2\), let \(\phi\) be the standard normal density, and put
\(\kappa=\phi(0)e^{-1/2}/2=1/(2\sqrt{2\pi e})\). Since \(\Phi''(x)=-x\phi(x)\), differentiation of \(x e^{-x^2/2}\) shows \(|\Phi''(x)|\le 2\kappa\) for every real \(x\). Taylor's theorem with its integral remainder therefore gives the global bound

\[
\Phi(z)\le\Phi(a)+\phi(a)(z-a)+\kappa(z-a)^2. \tag{2}
\]

This inequality is valid for all real \(z\), without restricting the sign or size of the conditional mean.

## 2. Truncated first and second moments

For \(\mu\in\mathbb R^5\), define \(Q(\mu)=\Pr[N(\mu,A)\in(0,\infty)^5]\), and put \(b=\nabla Q(0)\), \(H=\nabla^2Q(0)\). Completing the square in the Gaussian density proves

\[
\mathbb E[e^{t^{\mathsf T}Y}1_{\mathcal O}]
=e^{t^{\mathsf T}At/2}Q(At).
\]

Differentiating at \(t=0\) yields

\[
\mathbb E[Y1_{\mathcal O}]=Ab,\qquad
\mathbb E[YY^{\mathsf T}1_{\mathcal O}]=Aq+AHA. \tag{3}
\]

All these derivatives are justified by Gaussian domination on a neighborhood of the origin. Writing \(S^{i}\) and \(S^{ij}\) for the Schur complements obtained by conditioning the indicated coordinates of \(A\) to zero, boundary differentiation gives

\[
b_i=\phi(0)Q_4(S^{i}),\qquad
H_{ij}=\frac{Q_3(S^{ij})}{2\pi\sqrt{1-A_{ij}^2}}\quad(i\ne j),\qquad
H_{ii}=-\sum_{j\ne i}A_{ij}H_{ij}. \tag{4}
\]

For the last identity, express \(\partial_iQ(\mu)\) as the marginal density at the boundary in coordinate \(i\) times the conditional orthant probability of the remaining coordinates. At \(\mu=0\), the derivative of the marginal density vanishes. Increasing \(\mu_i\) shifts the conditional mean in coordinate \(j\ne i\) by \(-A_{ji}\mu_i\), while increasing \(\mu_j\) shifts that mean by \(\mu_j\). The chain rule then gives the displayed identity. The first two identities follow from integrating the Gaussian density on one and two boundary hyperplanes.

Since \(A\) is invariant under coordinate reversal, \(b_i=b_{6-i}\) with one-based indices. Thus only three distinct four-dimensional conditional orthant probabilities are needed. From (3),

\[
m_1:=\mathbb E[z1_{\mathcal O}]=\frac{rb}{\sigma},\qquad
m_2:=\mathbb E[z^2 1_{\mathcal O}]
=\frac{(1-\sigma^2)q+rHr^{\mathsf T}}{\sigma^2}. \tag{5}
\]

Combining (1), (2) and (5) bounds \(p_5\) above by
\(2[\Phi(a)q+\phi(a)(m_1-aq)+\kappa(m_2-2am_1+a^2q)]\). Since \(p_4=2q\), the certificate evaluates the following algebraically identical combined expression, which avoids paying twice for the numerical uncertainty in \(q\):

\[
U=p_1-p_2+p_3+A_q q+B_r(rb)+C_H(rHr^{\mathsf T}), \tag{6}
\]

\[
\begin{split}
A_q&=2\left[\Phi(a)-a\phi(a)+\kappa a^2+\kappa\frac{1-\sigma^2}{\sigma^2}-1\right],\\
B_r&=\frac{2[\phi(a)-2\kappa a]}{\sigma},\qquad
C_H=\frac{2\kappa}{\sigma^2}.
\end{split}
\]

Equations (1)–(6) prove \(\gamma_5\le U\). Every expression is evaluated using outward-rounded interval operations; in particular the negative values of \(rb\) and \(rHr^{\mathsf T}\) are not treated as unsigned estimates.

## 3. Computing the needed orthant probabilities

For a positive definite covariance matrix \(B\), let \(D=\operatorname{diag}(B)\), \(B(t)=D+t(B-D)\), and let \(S^{ij}(t)\) be the Schur complement of its \(i,j\) principal block. Define

\[
K_{ij}(t;B)=\frac{B_{ij}}{2\pi\sqrt{B_{ii}B_{jj}-t^2B_{ij}^2}}.
\]

Then

\[
Q_n(B)=2^{-n}+\int_0^1\sum_{i<j}K_{ij}(t;B)Q_{n-2}(S^{ij}(t))\,dt. \tag{7}
\]

To derive (7), differentiate the Gaussian density along \(B(t)\). Its derivative is \(\sum_{i<j}B_{ij}\partial_i\partial_j\varphi_{B(t)}\), as follows by differentiating its Fourier transform, or by Gaussian integration by parts. Integrating twice over the positive orthant in coordinates \(i,j\) leaves their joint density at zero, namely \(1/(2\pi\sqrt{B_{ii}B_{jj}-t^2B_{ij}^2})\), times the conditional orthant probability of the remaining coordinates. Its covariance is the Schur complement. Integration in \(t\) gives (7), since \(B(0)\) has independent coordinates. Positive definiteness holds throughout this compact interpolation, justifying differentiation and integration by uniform Gaussian domination.

The explicit small-dimensional formulas are

\[
Q_0=1,\quad Q_1=\tfrac12,\quad
Q_2(B)=\tfrac14+\frac{\arcsin\rho_{12}}{2\pi},\quad
Q_3(B)=\tfrac18+\frac1{4\pi}\sum_{i<j}\arcsin\rho_{ij},
\]

where \(\rho_{ij}=B_{ij}/\sqrt{B_{ii}B_{jj}}\). The two-dimensional formula is the angular measure of a sector for two independent standard normals. Expanding \(\prod_i(1+\operatorname{sgn}Y_i)/8\) proves the three-dimensional formula: odd sign moments vanish by symmetry, and pair sign moments are \(2\arcsin\rho_{ij}/\pi\).

Thus (7) evaluates \(Q_4(A_4)\), \(Q_5(A_5)\), and the three distinct \(Q_4(S^i)\) using only five one-dimensional integrals. All codimension-two probabilities in (4) use the explicit \(Q_3\) formula. The unnormalized Schur complement is computed as

\[
S_{xy}=B_{xy}-\frac{B_{jj}B_{xi}B_{yi}+B_{ii}B_{xj}B_{yj}-B_{ij}(B_{xi}B_{yj}+B_{xj}B_{yi})}{B_{ii}B_{jj}-B_{ij}^2}.
\]

No individual coordinate normalization is performed until the small-dimensional arcsine formulas require it. Formula (7) allows arbitrary positive variances.

## 4. Certified elementary arithmetic

The script uses Python's standard library only. Its numerical model is IEEE-754 binary64 round-to-nearest addition, subtraction, multiplication and division, together with correctly implemented `math.nextafter`. Every interval constructor rejects nonfinite or reversed endpoints. Each elementary operation moves the lower rounded result one representable number toward negative infinity and the upper rounded result one toward positive infinity. Products inspect all four corner products and reject nonfinite products before taking minima and maxima. Divisors containing zero are rejected. Squaring uses monotonicity on each half-line and returns lower endpoint zero when its argument straddles zero.

`math.sqrt` is used only as an initial guess. A candidate lower endpoint is reduced until its outward-rounded square is at most the input lower endpoint; a candidate upper endpoint is increased until its downward-rounded square is at least the input upper endpoint. These comparisons certify the final square-root bounds without assuming the initial guess's accuracy.

No library exponential, trigonometric function or normal CDF is trusted. The script uses

\[
\arctan z=\sum_{k=0}^{24}\frac{(-1)^k z^{2k+1}}{2k+1}+R,
\qquad |R|\le |z|^{51}/51\quad(|z|\le1).
\]

Larger positive arguments are reduced by repeatedly applying
\(\arctan z=2\arctan(z/(1+\sqrt{1+z^2}))\) until their interval upper endpoint is at most the binary64 literal `0.2`. Oddness handles negative arguments; intervals straddling zero are bounded at their endpoints. Then \(\arcsin x=\arctan(x/\sqrt{1-x^2})\) for \(-1<x<1\).

The same series encloses \(\pi=16\arctan(1/5)-4\arctan(1/239)\). To verify this identity, put \(a_0=\arctan(1/5)\). Double-angle formulas give \(\tan(2a_0)=5/12\) and \(\tan(4a_0)=120/119\); hence \(\tan(4a_0-\arctan(1/239))=(28680-119)/(28441+120)=1\). This angle is between zero and \(\pi/2\), so it equals \(\pi/4\).

For \(|x|\le1\), the exponential is enclosed by its sum through degree 24 and the geometric tail bound

\[
\left|e^x-\sum_{k=0}^{24}\frac{x^k}{k!}\right|
\le\frac{|x|^{25}}{25!}\frac1{1-|x|/26}.
\]

For \(|x|\le1\),

\[
\Phi(x)=\tfrac12+\phi(0)\sum_{k=0}^{20}\frac{(-1)^k x^{2k+1}}{2^k k!(2k+1)}+R_\Phi,
\quad |R_\Phi|\le\phi(0)\frac{|x|^{43}}{2^{21}21!\,43}.
\]

This follows by integrating the exponential series; the terms decrease in absolute value, so the first omitted term bounds the alternating remainder. These formulas suffice for \(\Phi(-1/2)\), \(\phi(-1/2)=\phi(0)e^{-1/8}\), and \(\kappa=\phi(0)e^{-1/2}/2\). All remainder arithmetic is outward rounded. Rational constants use interval division of exact integers, and \(\phi(0)=1/\sqrt{2\pi}\).

## 5. Integration and coverage

The class `J` carries interval enclosures of an integrand and its first and second derivatives. Sums and products follow the usual derivative rules. For a composition \(h(f)\), the first derivative is \(h'(f)f'\), and the second is \(h''(f)(f')^2+h'(f)f''\). The code uses

\[
(x^{-1})'=-x^{-2},\quad (x^{-1})''=2x^{-3},\qquad
(\sqrt x)'=(2\sqrt x)^{-1},\quad (\sqrt x)''=-(4x\sqrt x)^{-1},
\]

\[
(\arcsin x)'=(1-x^2)^{-1/2},\qquad
(\arcsin x)''=x(1-x^2)^{-3/2}.
\]

All their operations use the interval primitives above. On a cell of length \(w\), midpoint quadrature has error at most \(w^3\sup|f''|/24\). Taylor expansion about its midpoint proves this: the linear term integrates to zero and the integral of \((x-m)^2/2\) is \(w^3/24\). The certificate intersects this enclosure with the natural interval range of the integrand times the cell length; both independently contain the integral.

The cover starts with \([0,1]\), and a cell is replaced only by its two midpoint children. **Before every cell evaluation**, depth greater than 50 is rejected. Thus all endpoints, lengths, and midpoints used in accepted cells are exact binary64 dyadic values. A separate strict interior-midpoint check also prevents degenerate splitting. These rules prove that the current leaves always cover \([0,1]\) with disjoint interiors. A failed interval domain check triggers subdivision and supplies no accepted enclosure.

Approximate totals and a priority queue choose subdivisions and a stopping point; they are not trusted evidence. Once the approximate stopping criterion is met, every accepted leaf interval is summed again using outward-rounded interval addition. This final sum is the certified integral. The script declares success only if the resulting upper endpoint of \(U\) is strictly below the lower endpoint of the interval quotient \(1341/10000\). Consequently adaptive choices or approximate stopping totals cannot cause a false success.

Printed enclosure endpoints are converted from their exact binary64 values to decimal strings with 18 digits after the point, rounding lower endpoints toward negative infinity and upper endpoints toward positive infinity using `Decimal.from_float` and explicit directed decimal rounding. The final inequality is checked before printing, using interval arithmetic on the exact rational target.

## 6. Reproduction and scope

The SHA-256 of `certify.py` is
`9d27d3aab5e085074477d529160c9e54ba2fa01fdb9d6f9d2313a6071db0ab09`.

```
python3 -I -B certify.py
```

The program reads no files, uses no randomness, and writes one JSON result to stdout. An input JSON value may be supplied on stdin and is ignored. `certified: true` proves (C) through (1)–(7) and the validated arithmetic above. The probability model, graph reduction and circuit lower-bound theorem remain separate mathematical obligations.

The local run used CPython 3.9.6 and completed in approximately 7.99 seconds. It reported

\[
U\in[0.134072072323763391,\ 0.134085456508304113],
\qquad\gamma_5\le U<0.1341.
\]

The five integrations used 144, 148, 83, 48 and 26 leaves respectively, with maximum depths 9, 9, 8, 7 and 5. The minimum LDL pivot was greater than `0.269752650521788150`. These finite checks used only the submitted Python source. Hosted execution and review statuses are reported separately; a local successful run is not official acceptance.
