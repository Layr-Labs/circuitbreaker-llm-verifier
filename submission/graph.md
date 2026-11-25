# Bounded-local fair cuts approaching the certified energy ceiling

A bounded-radius i.i.d. rule on the matched cubic tree is equivariant and fair, with expected crossing count per vertex below 12262902000000000000/10^20, counting edges once. Proof §1 transfers it to circuit graphs. Both phases in §5.1 are required.

## 1. A globally compatible Gaussian source, including branched neighborhoods

Contract the perfect matching to anchors of a loopless 4-regular multigraph Q. Pair its four labeled incidences into two straight continuations, one per physical vertex. The marked universal tree has two lines at each anchor. Set t=41/57, a=28/57, delta=t-a, so t^2+2a^2=1.

Give each anchor an independent standard Gaussian Z. For the oriented edge uv set B^0_uv=Z_u. For B^(j+1)_uv, take t times the depth-j inward message arriving at u from the straight opposite direction and a times each of the two others. The three predecessor half-trees are disjoint, so every B^j has variance one. The shells on uv use disjoint distance layers. Their normalized sum B^[J]=(B^0+...+B^(J-1))/sqrt(J) has variance one. The four inward sums at a root are iid, and inward families on disjoint half-trees are independent. The finite rules are equivariant and local. Applying the t,a,a recurrence to the normalized sums differs from the sum at uv only in the first and next omitted shell, giving mean-square norm sqrt(2/J). A diagonal subsequence of the bounded covariance arrays on the countable tree gives a Gaussian limit with the exact recurrence and the same half-tree independence.

The root score R is half the sum of its four incoming modes. Conditional on an entering B=t B_f+a(B_l+B_r), generate the three farther modes by independent standard normals Z,Y: B_f=tB+sqrt(2)a Z, W_0=sqrt(2)a B−t Z, B_l=(W_0+Y)/sqrt(2), B_r=(W_0−Y)/sqrt(2). The recurrence determines the backward mode. Farther half-trees are conditionally independent given their entering modes, because the inward families were independent and upstream data use each only through its entering mode.

Here is the two-real-state version used by the certificate; all transition coefficients are rational. Generate independent standard R,U shared by the two lines. Draw an independent V~N(0,1/2) for each line ell=+1,-1. For each initial ray on that line, set
```
D = a(R+ell U)/2 +/- aV,    E=((a+t)R+delta ell U)/2,    P=D+E.
```
At a fresh site draw W~N(0,1/2). Its anchor score is X=P-delta W. The straight continuing ray has new states
```
D_s=tD+2a^2 W,    P_s=tP+aD+2a^2 W.
```
For the perpendicular pair take another independent V~N(0,1/2); their two states are
```
D_o=a(D-tW) +/- aV,
P_o=tD+aP+(a^2-t^2-at)W +/- aV.
```
Beyond the displayed states the three half-trees use independent innovations. With entering mode B, D=aB. The new forward mode is tB+2aW and the backward mode is 2E−tB. The perpendicular pair sum is 2(D−tW), and the other pair sum is 2(E+aW). For either line E is half the sum of t times its own pair and a times the other pair. Substitution gives the recursions. The stationary score law and the energy and gain certificates use the covariances below. On a fixed finite ball the finite-shell thresholds converge in law, since each limiting equality has probability zero.

### 1.1 Exact line and crossing covariances
The root law and conditional recursion determine the limit, which inherits all marked-tree symmetries. For a straight ray at its first site, (D,P) has coefficients a/2 and (2a+t)/2 on R. The straight update sends these coefficients to (tD,tP+aD), and W is independent of R. Thus
    r(d) = (a+t/2) t^(d-1) + (d-1) a² t^(d-2)/2
for distance d>=1, omitting the second term at d=1, with r(0)=1. In any consecutive block the first R and the subsequent fresh W enter triangularly with nonzero diagonal, giving positive definiteness.

The coefficient of U at distance d is ell[r(d)-h(d)], with h(d)=a t^(d-1), by the same update starting at ell t/2. Lines with opposite ell share only R,U. The covariance between distances i,j on those lines is
    r(i)r(j)-(r(i)-h(i))(r(j)-h(j))
      = h(i)r(j)+h(j)r(i)-h(i)h(j).
Use r(d) for both same-ray and opposite-ray distances on one line.

## 2. Initial line-minimum policy and the first-pass energy

Each physical vertex has two C edges and an M pin. Its split indicator is NAE. The minimum split count on a line is the minimum vertex cover size of the raw disagreement paths, as proved in §4. It equals the maximum matching size by adjacent pairing. Greedy selection from left to right is optimal by exchanging the earliest edge into a maximum matching and inducting. For disagreement indicators d_i the greedy selection satisfies g_i=d_i(1−g_(i−1)). Nineteen substitutions leave a nonpositive residual. Averaging with A_j(T) the j-run probabilities gives
    c1 <= S2(T)+A3(T)-A4(T)+A5(T)
            -(A6-A7)(T)-(A8-A9)(T)-(A10-A11)(T)-(A12-A13)(T)
            -(A14-A15)(T)-(A16-A17)(T)-(A18-A19)(T).
Here S2(T)=P(X_0>T,X_2<=T)=A1-A2: the two adjacent disagreements have XOR equal to the endpoint disagreement; stationarity makes its two orientations equally likely. In a finite window greedy vs full-line counts differs only within its initial disagreement component; almost surely finite components suffice upon taking expectations and growing window size.

Use the least C-edge optimizer on finite windows, from submodularity. A ray context is the clamped cost at 0 minus cost at 1. A raw pin b sends 2b−1, so f(x,b)=sign(x+2b−1). Two repeats send 2b−1 irrespective of the boundary. Thus the context stabilizes by §4. For incoming x,y the two least C bits are [x+f(y,b)>0] and [y+f(x,b)>0]. Each edge uses the contexts at its two ends. Truncate rays at a fixed L and then use sufficiently many Gaussian shells to approximate that finite neighborhood.

### 2.1 Directed moment and energy interface

For centered consecutive Gaussian line scores, write s_i=(-1)^i. An alternating word plus its central reflection on m positions has total probability
    2^(1-m) * sum over even I⊂{0,...,m-1} of (product of s_i in I) E product of sign(X_i) on I.
Thus zero-threshold U5=S2(0)+A3(0)−A4(0)+A5(0) and the needed A_j(0) bounds follow from the centered sign intervals in §2.2. The exact outward cost follows by substituting the inequalities in §2.3; no frozen U5, moment, or gain outputs are read.


### 2.2 Directed Gaussian sign moments

Fix any principal block of the six consecutive Gaussian line scores, r(d) as in the gain proof. For an even sign moment, select the upper central index i (rather than the first) and multiply all covariances incident to it by v from zero to one. At zero the sign moment vanishes. For each other j the derivative contribution is
  (2/pi) * G_ij / sqrt(1-v²G_ij²) * M_remaining(v),
where M_remaining is the centered sign moment in the conditional covariance after pinning i,j at zero. This follows by differentiating a Gaussian density in a symmetric off-diagonal entry to obtain its mixed spatial partial, integrating twice by parts across the two zero hyperplanes, and multiplying by the pinned bivariate Gaussian density at zero. Compact uniform positive definiteness ensures Gaussian domination. At dimension two the angular measure of the intersection of two Gaussian halfspaces gives the sign moment 2 asin(rho)/pi. For the remaining indices k,l the Schur covariance equals
  G_kl-G_kj G_lj - (v²/(1-v²G_ij²))(G_ki-G_ij G_kj)(G_li-G_ij G_lj).
Normalize by the square roots of its diagonal entries. All true real-node covariance matrices are positive definite; the interval arithmetic encloses those true values. The program tests true rational positive pivots of the original six-point covariance.

For completeness, the same uniform complex bound applies to every selected central coordinate after conditioning. The diagonals of the inverse of the exact original six-point covariance are <=20/3, checked rationally. The bound is inherited by principal submatrices (inverse-diagonal variational formula sup_z(2z_i-z^T Gz)), by disconnecting the center (principal block plus singleton), and by convexity for real v. The inverse of a Schur complement is a principal block of the unconditioned inverse; multiplying by conditional variances <=1 when normalizing preserves the diagonal bound. Hence on disconnecting any selected center, whitening the other block and rotating sends the varying covariance to a diagonal matrix with eigenvalues 1+sv,1-sv and the remaining ones, where s²=1-1/(G^-1)_ii<=17/20; all coordinate changes can be chosen real and independent of v. The sign integrand after coordinate change is a fixed measurable function bounded by one.

Map the closed Bernstein ellipse of parameter 7/4 to [0,1]. The extreme real parts are -9/112 and 121/112. The absolute imaginary part is at most 33/112. Since s^2<=17/20, s<461/500. The two eigenvalues lambda_1, lambda_2 therefore have real parts >=1/256, at least one real part >=1, and absolute imaginary parts <=9/32.

Use the right-half-plane square root. The absolute integral of the complex Gaussian density is the product of sqrt(|lambda_i|/Re lambda_i) over the two varying directions. The two ratios |lambda_i|/Re lambda_i are bounded by 73 and 17/16 using the bounds on the real and imaginary parts above, so the absolute integral is less than 9. Under the normalized absolute density, the coordinate second moment is |lambda_i|^2/Re lambda_i. Differentiating the density therefore contributes at most (9s/2) sum_i (1/|lambda_i|+1/Re lambda_i) <= 9s*(256+1) < 4096. Strictly positive real parts of the inverse eigenvalues give analytic domination in a neighborhood. This bound applies to the combined derivative. The checker adds the quadrature error after it sums the contributions of all j. It encloses each nested four-sign moment before use.

#### Rule, nodes, weights and directed arithmetic

Gauss–Legendre quadrature is exact through degree 2N−1. Choose the roots of P_N by the Legendre recurrence. At a root x, the mapped weight on [0,1] is (1−x²)/(N² P_(N−1)(x)²) by the Christoffel–Darboux and derivative identities.

Opposite signs at 32 disjoint 2^-72 brackets cover all roots of P_32. Interval arithmetic verifies positive weights and their sum bracket. It rounds outward at 2^72 and uses integer square roots. Machin's identity gives pi. Two half-angle reductions express arcsine as arctangent, enclosed by alternating series and the next term. Clip correlations to [-1,1] and add the analytic integration error.
For an analytic derivative bounded by 4096 on the ellipse of parameter 7/4, the Laurent coefficient estimate via the Joukowski map bounds Chebyshev coefficient k>=1 by 2*4096*(4/7)^k. Truncating after degree 63 incurs uniform error at most 2*4096*(4/7)^64/(1-4/7). The integral and positive normalized GL32 rule each have norm one, so their difference is at most twice this outward error. Add it at each center-row integration, including each real-node inner four-sign moment. This proof works identically for the higher central interpolation coordinate chosen by the actual code.


### 2.3 Shift + finite-tail discounts

Let u=291/800 and kappa²=(1−r(2))/(1+r(2)). Bivariate conditioning and differentiation at both boundaries give −S2'(u)=phi(u)(2Phi(kappa*u)−1). Thus S2(0)−S2(u) is the integral of exp(−(s²+w²)/2)/pi over 0<=w<=kappa*s, 0<=s<=u. Substitute w=kappa*s*v and integrate the degree-five Taylor lower polynomial. Multiply the positive result by the lower square root of kappa² divided by the upper pi.

#### Moments of the Gaussian shift

For n=4,5,6 set C_ij=r(|i-j|), K=C^-1, l=K1, D=1'K1, Y~N(0,C), and Z=l'Y. Write s_i=(-1)^i and let O be the union of the two opposite alternating orthants. Use A(T)=P(Y-T1 in O) and M_(2k)=E[1_O Z^(2k)]. Gaussian density translation and central symmetry give
```
A(T)=exp(-DT^2/2) E[1_O cosh(TZ)].
```

For a set S of coordinates pinned at zero, let F(S,k) be the conditional integral of Z^k over the remaining signs s_i, times the marginal density of Y_S at zero. For empty S the density factor is one. Let J be the unpinned coordinates, C_S the conditional covariance on J, v=C_S l_J, and d=l_J'v. Since C l=1, the checker computes v_j=1-C_(j,S) C_(S,S)^(-1)1 and d=D-1'C_(S,S)^(-1)1, with zero correction for empty S.

The centered conditional density g satisfies (l_J'y)g(y)=-sum_j v_j partial_j g(y). Integrate this identity against (l_J'y)^(k-1) over the remaining signed orthant. Differentiation of the power gives (k-1)d F(S,k-2). The boundary at coordinate j=0 contributes s_j v_j F(S union {j},k-1), since restricting the joint density gives the additional pinned density factor. Gaussian decay removes the boundaries at infinity. Thus for k>=1,
```
F(S,k)=(k-1)d F(S,k-2) + sum_(j in J) s_j v_j F(S union {j},k-1),
```
where the first term is omitted when k=1. Central symmetry gives M_(2k)=2 F(empty,2k). The checker evaluates these moments for k=1,2,3 with directed intervals.

The recurrence reaches k=0 only with even |S|. For empty S, F(empty,0)=A(0)/2 by §2.1. Otherwise there are at most four unpinned coordinates. Expand the conditional orthant indicator as a product of signed factors 1+s_j sign(Y_j), divided by 2^|J|. Odd sign moments vanish and each two-sign moment is 2 asin(rho)/pi. If |J|=4, include the four-sign moment from §2.2. The standardized conditional covariance obeys the inverse-diagonal bound from that section. Divide the resulting probability by (2 pi)^(|S|/2) sqrt(det C_(S,S)). Reuse the directed square roots and two-sign routine of §2.2.

As a check, the second moment also equals D A(0)+A''(0). Differentiating each signed orthant probability twice with respect to the common threshold gives
```
A''(0)=sum_(i!=j) s_i s_j sqrt((1-C_ij)/(1+C_ij)) / pi
          * P(s_l Y_l>=0 for l outside {i,j} | Y_i=Y_j=0).
```
On the first derivative, the threshold on face i moves at speed one. Given Y_i=T, each remaining centered threshold moves at speed 1-C_ji. The derivative of the marginal density at zero vanishes, and adding the opposite orthants doubles the conditional contribution. The checker evaluates all 62 ordered pair contributions and verifies that each second-moment interval intersects its corresponding recurrence interval.

#### Upper and lower bounds for the shifted orthants

Whiten Y to R theta. On O, Z^2/R^2 <= beta by the exact face calculation below. The radial variable R^2 is chi-square with n degrees of freedom, independent of theta. Thus for k>=3,
```
M_(2k) <= M_6 beta^(k-3) product_(j=3 to k-1)(n+2j).
```
For n=4 and n=6 use upper M_2, M_4, and M_6 in the first three nonconstant terms of cosh(TZ). The subsequent bound terms have ratios T^2 beta (n+2k)/[(2k+1)(2k+2)]. This ratio decreases with k, so the checker bounds the tail by the geometric sum with the first ratio at k=3. If H is the resulting upper bound for E[1_O (cosh(TZ)-1)] and e+ is the positive degree-six upper Taylor bound for exp(-DT^2/2), with e+<1, then
```
A(T)-A(0) <= (e+-1) A(0)_lower + e+ H.
```
For n=5, A(T) enters with a minus sign. Use the positive degree-five lower exponential bound e-<1 and the lower moment bounds in the first three nonconstant terms. Their nonnegative tail gives
```
A(0)-A(T) <= (1-e-) A(0)_upper
              - e- [T^2 M_2_lower/2 + T^4 M_4_lower/24 + T^6 M_6_lower/720].
```

#### Exact active-face certificate

For m=4, 5, and 6 let C be the stationary line covariance matrix with entries r(i-j). Write K=C^-1, D=1^t K1, sigma_i=(-1)^i, M=diag(sigma)K diag(sigma), v=diag(sigma)K1. To maximize (v^t y)^2 over y>=0 with y^t M y<=1, enumerate both signs s and every nonempty support S. Set b=s v_S, z=M_SS^-1 b and retain the supports with z>=0. The maximum of b^t z and zero is the exact squared bound beta. Indeed a maximizing nonzero positive objective lies on the sphere and is interior to its positive support, where the Lagrange equation makes it a positive multiple of z. Conversely every retained nonzero z normalized to unit M_SS norm and extended by zeros is feasible and realizes that value. This verifies all faces by exact rational inversion without assuming which support maximizes the expression.

#### Longer alternating-run tails
For even m from 8 through 20, the discount (A_(m-2)−A_(m-1))(T) is the probability of the two opposite words with the first m-1 bits alternating and the last bit equal to its predecessor. Their centered signs are sigma_i=(-1)^i for i<m-1 and sigma_(m-1)=1. Sum the two §3.2 word bounds, using separate scales and shifts. Take the maximum with Jensen when available.

Let K=C_m^{-1}, D=1^t K 1, and Y_i independently halfnormal with positive rational scales s_i. Let M=diag(sigma)K diag(sigma), E=E[Y^t M Y]=sum K_ii s_i^2+(2/pi)sum_{i!=j} sigma_i sigma_j K_ij s_i s_j and L=sum_i sigma_i s_i sum_j K_ij. For the two opposite words translate by T*1 and write the Gaussian integral as an expectation under the product halfnormal density. Each orthant has common normalizing multiplier prod(s_i)/(2^m sqrt(det C_m)); the exponent inside the expectation has mean (m-E-D T^2)/2 +/- T sqrt(2/pi)L. Jensen on each orthant separately and summation give
A_{m-2}(T)-A_{m-1}(T) >= prod(s_i)/(2^{m-1}sqrt(det C_m)) exp((m-E-D T^2)/2) cosh(T sqrt(2/pi)L).
This directly uses the biased threshold instead of applying a separate worst-case translation lower. Use the explicit rational scale vectors in energy.py. They were selected numerically but all claims use only their explicit rational values, not optimality. Bound E from above by rational Machin bounds for pi, with sign-aware off-diagonal division. For any real z the odd Taylor polynomial through degree 11 lies below exp(z) (Taylor remainder nonnegative); verify positive at the chosen lower exponent. Bound cosh below by 1+T^2 L^2/pi_upper, and sqrt(det C_m) above by an authenticated integer square root at denominator 10^24. Exact SPD elimination validates the inverse and determinant.

## 3. Tight gain at the active threshold

Take T=-291/800, a=28/57, t=41/57. A raw pin is 1 when its score exceeds T. At one matched root, use the thresholded root R and four consecutive straight rays A,B (opposite directions of one line), C,D (opposite directions of the other). On a ray, truncate to six pins with largest external message +1. A pin b=1 sends 1 if its incoming message is not -1 and 0 otherwise; b=0 sends -1 if its incoming message is not +1 and 0 otherwise. This recursion is monotone in the incoming message. Write M, H, Z for the events that the resulting upper message is respectively -1, nonpositive, and exactly zero. Thus M and H imply the corresponding downward events for actual messages, though Z alone need not. We use Z for partitions of the *truncated* event only.

For two rays of the same endpoint at a positive root, the first-pass line optimization gives local improvement +1 when both messages are nonpositive and at least one is negative (Y), and nonnegative improvement when at least one is negative (N). These identities follow by substituting the 3 x 3 incoming trit pairs into the two adjoining edge bits [x+f(y,1)>0], [y+f(x,1)>0]. Since the two physical endpoints can change their shared matching bit only together, the nonnegative improvement for a pair equals 1[Y_1 N_2]+1[Y_2 N_1] if the root bit is 1. On dividing by two vertices and using line exchange, it suffices to lower-bound P(R>T, Y_AB, M_C or M_D). Starting with the largest boundary +1 makes this a lower bound for every deeper boundary. Other gains are not needed.

Use the covariance in §1.1 for the root and all ray positions. The straight recursion for (D,E) follows from E=P−D in §1.

### 3.1 Disjoint cylinders, not a sum of overlapping orthants

For p=P(R>T,Y_AB,M_C), split Y_AB into (M_A,H_B) and (Z_A,M_B), where H=M disjoint union Z. Evaluate all 64 six-pin strings and merge sibling prefixes with identical trits. Products of the indicated classes give 72 disjoint cylinders, 54 up to reflecting the first line. Record only coordinates in their prefixes.

Extending a ray from six to eight pins with boundary +1 can only lower its trit, by monotonicity. Partition each ray by its pair of six-pin and eight-pin trits using the same prefix merge. Add all 71 reflection classes of prefix triples passing the eight-pin condition and failing the six-pin one. Their cylinders are disjoint from each other and the six-pin event. The same 12 classes retain their rational witnesses. Seven coordinate sweeps propose positive scales for the other 59, rounded to denominator 10000. The exact lower bound only requires positive scales, and uses expansion center zero there. The overlap upper still applies because a negative trit requires a first zero and a zero among pins two and three, and a nonpositive trit requires a zero among its first two pins.

### 3.2 A Gaussian orthant lower bound with a positive remainder

For one selected cylinder take covariance G, precision K=G^{-1}, signs s_i (+ for the root and all 1 pins, - for 0), positive scales v_i, U=diag(s_i v_i), D=1^T K1. Under independent H_i=|standard normal|, set X=T1+UH. Change variables on the orthant to obtain the probability
    prod(v_i)/(2^n sqrt(det G)) E exp(L),
    L=H^T B H + l^T H - T^2 D/2,
    B=(I-UKU)/2,  l=-T U K1.
Put mu=E H_i=sqrt(2/pi), eta_i=H_i-mu, m_k=E eta_i^k, w=m_2. Raw halfnormal moments satisfy a_0=1,a_1=mu,a_k=(k-1)a_{k-2} by Gaussian integration by parts. Binomial expansion gives all m_k through k=6. In particular
    EL = w tr(B)+mu^2 sum_ij B_ij + mu sum_i l_i - T^2 D/2.
Set b_i=2mu sum_j B_ij + l_i, d_i=B_ii, e_ij=2B_ij for i<j, Hhat_i=b_i eta_i+d_i(eta_i^2-w). Then X0=L-EL is the sum of the centered Hhat_i and e_ij eta_i eta_j. Its variance v and third centered moment k are computed by:
    v = sum_i E Hhat_i^2 + w^2 sum_{i<j} e_ij^2,
    k = sum_i E Hhat_i^3
        + 6 sum_{i<j} e_ij E[Hhat_i eta_i] E[Hhat_j eta_j]
        + 3w sum_{i<j} e_ij^2 (E[Hhat_i eta_i^2]+E[Hhat_j eta_j^2])
        + m_3^2 sum_{i<j} e_ij^3
        + 6w^3 sum_{i<j<k} e_ij e_ik e_jk.
To see exhaustiveness, expand the square/cube and use independent coordinates. A lone eta at any index vanishes. Three univariate terms must belong to one index; a single edge needs one univariate term at each endpoint; two edges and one univariate term must repeat the same edge with the remaining term incident to it; three edges must be either identical or form a triangle. These alternatives give the displayed coefficients.

For any rational center s and any real y the Taylor cubic 1+y+y^2/2+y^3/6 <= exp(y), with nonnegative fourth-order remainder even for negative y. Thus
    E exp(L) >= exp(EL+s) [1-s+(v+s^2)/2+(k-3sv-s^3)/6].
The average polynomial is checked positive after adding a directed lower bound for the Taylor remainder. For any real y the integral remainder after the cubic is
```
e^y - (1+y+y^2/2+y^3/6) = (y^4/6) integral_0^1 (1-u)^3 e^(uy) du.
```
Let R(y)=exp(y)-1-y-y²/2-y³/6 and r(t)=R(-t) for t≥0. Choose a positive rational b≤sqrt(E[Y²]). Put f(t)=r'(t)/t and B(t)=r(t)-r(b)-f(b)(t²-b²)/2. Direct differentiation yields f'(t)=∫_0^1 u(1-exp(-tu))du, increasing and concave with f'(0)=0. For t≤b,
    B(t)=1/2 ∫_t^b (z²-t²)f'(z) dz ≥ f'(b)(b²-t²)²/(8b).
The inequality uses f'(z)≥zf'(b)/b. For t≥b, write B(t)=∫_b^t u∫_b^u f'(z)dz du≥ b f'(b)(t-b)²/2. Thus B(t)≥κ(t-b)² for κ=b f'(b)/8, or any directed lower bound on κ, and for t≤b the bound strengthens to κ(t-b)²(1+t/b)².

For y≤0 this immediately bounds R(y) below by r(b)+f(b)(y²-b²)/2+κ(y+b)². For y≥0, R(y)-r(y)≥y⁵/60 by the odd exponential series. When y≥b this gives a surplus at least (b³/240)(y+b)², since y⁵/(y+b)² is increasing. For 0≤y≤b, set x=y/b. Divide the surplus over r(b)+f(b)(y²-b²)/2 by (y+b)². The stronger B bound gives at least κ[(x-1)²+λx⁵] for λ=b³/(240κ). For any a≥0, x⁵≥5a⁴x-4a⁵. Completing the square gives (x-1)²+λx⁵≥1-4λa⁵-(5λa⁴-2)²/4. Maximize this lower bound over the finite rational anchors and zero, and cap it at one. For all y, the result times κ is a lower bound K on the surplus divided by (y+b)². The bound for y≥b follows since every anchor bound is at most λ, by testing x=1.

Therefore E R(Y)≥r(b)+K(E[Y²]-2bs+b²), dropping the nonnegative tangent contribution f(b)(E[Y²]-b²)/2. Compute b by a directed square root of the lower second moment and r(b) by the degree-15 Taylor minorant. Enclose f'(b) below by 1/2-(1-(1+b)P_15(-b))/b², or use f'(b)≥b/[3(1+b)] from its integral and 1-exp(-z)≥z/(1+z). All operations after the root are rational. Use zero for a negative interval lower bound on E[(Y+b)²]. If b=0 the unsigned Jensen bound suffices. Also retain the directed bound M2²/[6(4+sqrt(M2))] for M2=E[Y²]: the law with density 4(1-u)^3 has survival at most exp(-4u), so the integral is at least 1/(4+|Y|), and Cauchy–Schwarz gives that bound. Take the maximum. At expansion center zero the mean factor is also at least one by Jensen.

### 3.3 A quantitative overlap upper

By symmetry the desired probability equals 2p-q where q=P(R>T,Y_AB,M_C,M_D). Split which first-line ray is negative (cost a factor two by a union bound). Say A is negative, B nonpositive. Every negative ray must have pin 1 negative and at least one of pins 2 or 3 negative (here indices are positions starting at one); a nonpositive ray has at least one of pins 1 or 2 negative. This follows directly from the first three bits or the exhaustive trit control. Thus, for each orientation q is covered by 16 events with root positive and the seven indicated coordinates negative: A1,Aa,Bb,C1,Cc,D1,Dd with a,c,d in {2,3}, b in {1,2}. Some overlap; summing *upper* bounds remains safe.

For any of those eight-dimensional marginal events let Sg=diag(+,-,...,-), M=Sg K Sg, b=T Sg K1, C=T^2 1^T K1. All eight b_i>0 are checked exactly. Change variables to signed, unscaled deviations y_i>=0. Put r=b^T y and b'=b-T e_0. Then b'^T y>=r since T<0 and y_0>=0. From M^{-1}b=T Sg 1 and (M^{-1})_00=1, we get b'^T M^{-1}b'=C-T^2=V. The checker verifies V>0. Cauchy--Schwarz gives y^T M y>=r^2/V. Integrating a nonnegative function of r on the orthant uses density r^7/(7! prod b_i) by scaled simplex volume. Up to 5, partition r into 192 equal panels. On each panel [l,h] the decreasing factor exp(-r-r^2/(2V)) is at most its value at l, and the integral of r^7/7! is (h^8-l^8)/8!. Enclose exp(-x) from above for x>=0 by the reciprocal of sum_{j=0}^{12} x^j/j!, and round each panel sum upward at the probability scale in gain.py.

On the tail write L=5 and v=L^2/(2V). Bounding exp(-r) above by exp(-L) and substituting z=r^2/(2V) gives the upper expression
```
V^4 exp(-L-v) (1+v+v^2/2+v^3/6)/105.
```
The identity follows by three integrations by parts on the gamma integral. Apply the same exponential enclosure to the tail and the density factor exp(-C/2). Finally divide the sum, with its exp(-C/2) factor, by (2pi)^4 sqrt(det G) prod_i b_i using lower bounds in the denominator. Sum upward enclosures over the 16 events and the two orientations.
No independence across the roots' two lines is asserted.


Every finite matrix is scaled by a common denominator to an integer Gram matrix. Fraction-free Gauss-Jordan computes its adjugate with positive leading pivots (Sylvester determinant identity for the scaled row update); the checker multiplies it back to the computed determinant times identity, then scales the inverse and determinant appropriately. Gaussian square roots use integer squaring, pi is bracketed by a rational alternating Machin formula, and signed products round outward. Positive probability comparisons cannot disappear under optimized Python. The active rational energy-ceiling comparison in energy.py uses its freshly recomputed gain, not a claim copied from a different threshold.

## 4. The actual initial majority exceeds one half

It applies to the actual *least global* first-pass line minimizer. It does not insist that positive agreement C-edges stay 1 when both endpoints are split. The initial matching repair may turn a positive raw pin off (take 0 on ties) but never turns a negative raw pin on.

Let raw bit 1 mean that the Gaussian anchor score exceeds T. On one line, consider the graph whose edges join consecutive discordant raw pins. Its nontrivial components are alternating paths, almost surely finite in the limiting stationary Gaussian law. The number of split vertices in any C assignment is at least the minimum vertex-cover size of this disagreement graph: an edge with discordant endpoint pins forces a split at some endpoint. Conversely given a cover, unsplit vertices prescribe their C incidences consistently, because every disagreement edge is covered, so the minimum equals the sum of the component minimum cover sizes. In every minimizing assignment the split set is a minimum cover. Isolated vertices of the disagreement graph are never split.

An odd disagreement path on 2k+1 vertices has a unique minimum cover of size k, at alternating internal positions. Including an endpoint leaves a matching of size k in the remaining path. Excluding it forces its neighbor, so use induction. The C edges equal the common raw bit of uncovered vertices. On an even path the positive vertices form a minimum cover, so some optimum makes all its disagreement C edges zero. Use the common bits on agreement edges in that optimum. The least optimizer keeps those zero disagreement edges and may erase positive agreement edges.

Fix one physical endpoint at a root anchor with raw-positive pin and let A,B,C,D be the indicators that its four nearest raw neighbors are negative, with A,B on its own line and C,D on its partner's line. A raw-positive endpoint can lose majority after matching repair in only these cases:
1. Its original first-pass C bits are both zero. If both nearest neighbors are negative, charge (A B). If exactly one is negative, the other C edge is a positive agreement edge that took zero. Its other endpoint must also be split, so its farther C neighbor must be negative. Thus it is an occurrence of the four successive raw bits 0110, in one of two orientations. It cannot happen with neither nearest neighbor negative, as the endpoint then is isolated in the disagreement graph.
2. Its C bits are mixed and M was turned off. Since two disagreement edges of a run have the same bit and an isolated endpoint never splits, exactly one of A,B equals 1. For the shared matching bit to turn off, its other endpoint must not have both C bits one: if it did, keeping M would strictly lower its cost (the first endpoint is mixed). The other endpoint cannot be isolated, so C or D equals 1. Thus this loss is bounded by (A XOR B)(C OR D). For binary a,b,c,d one has
       (a+b-2ab)(c+d-cd) <= ac+ad+bc+bd-2ab(c+d).
Combining reflection and line-exchange symmetry gives an upper bound 4p_cross−4p_triple, where p_cross=P(R>T,A=C=1), p_triple=P(R>T,A=B=C=1).

There is a disjoint positive putback with raw-negative root. On its own line the raw five-pin pattern 11011 (center zero) puts the root at the unique covering position of its three-vertex odd disagreement path, so both its C bits are one irrespective of further pins. The one-sided matching repair cannot undo that majority. If p0=P(X_0>T) and y=P(R>T,A=B=1), a lower of the initial majority is
       p0 - y - 2*p_0110 -4*p_cross +4*p_triple +p_11011.       (1)
The consecutive-repeat event has positive probability. Covariance decay r(d)→0 gives asymptotically independent distant repeat indicators, so the variance of their ray average tends to zero. A ray with no repeat would have average zero, an event of probability tending to zero. Continuous finite-dimensional Gaussian laws and stabilization preserve all strict inequalities at sufficiently large finite radii.

### Directed calculation of y and p_cross

The root-neighbor covariance is r=97/114. The covariance of its two considered neighbors is either r_2=529/722 (opposite directions of one line) or 2 a r-a^2 (different lines). Condition on R=z. Each neighbor then has mean rz, variance v=1-r², and conditional correlation rho=(r_neighbor−r²)/v. Their centered mean and half-difference are independent normals of variance v(1+rho)/2 and v(1−rho)/2. The conditional probability both are below T is therefore the expectation over S=|N(0,1)| of

    Phi((T-rz)/sqrt(v(1+rho)/2)-sqrt((1-rho)/(1+rho))*S).

Every factor decreases with z and S on z>=T. The exact upper sum uses z∈[T,5], S∈[0,5], rectangle widths ≤1/64, upper Gaussian mass of each rectangle (factor 2 for S), and its lower z,S endpoints. The conditional argument is rounded *up* by an upper denominator in its negative first term and a lower coefficient for the subtractive term. Remaining tails are charged one, total ≤3[1-Phi(5)]. CDF arithmetic uses fixed midpoint h=1/1024 and Taylor brackets with outward fixed point. Derivative |phi''|≤1 follows by considering (x²−1) exp(-x²/2) on x²≥0, whose positive maximum is 2exp(-3/2). The full complete midpoint error is at most 5 h^2/24. For CDF requests between grid points use monotonicity and the neighboring grid brackets rather than fractional-panel interpolation. The Mills tail is checked via the outward upper at x=5, and negative arguments use exact CDF reflection. The other numerical branch computes p0=Phi(-T) with a downward bracket.

### The agreement-edge case: 0110

Write s=(-1,+1,+1,-1), C=(r(|i-j|)) for four consecutive line coordinates. At threshold zero the orthant probability of this word is

    (1+sum_{i<j}s_i s_j (2/pi)asin(C_ij)+M4)/16,

where M4 is the four-sign moment supplied by fresh validated quadrature on that covariance. Denote this probability as a function of the common threshold by p(t). Differentiation on the i-th face contributes at t=0
    −s_i phi(0) * P(other signs match | X_i=0)
and the conditional three-sign probability is (1+sum three pairs of signed (2/pi)asin(conditional rho))/8. Conditional signs of odd degree vanish at zero. Conditional covariance is C_jk−C_ji C_ik; all the arcsines and square roots are directed. Thus p'(0)>0 is checked from exact covariance and finite brackets.

The orthant probability p(t) is log-concave by Prékopa’s theorem. After translation by t·1 the shifted Gaussian density times the fixed convex orthant indicator is jointly log-concave. Its directed tangent at 0 has positive slope, so p(T)≤p(0)<.0121 for T<0. [Prékopa 1973](https://acta.bibl.u-szeged.hu/14411/)

### Remaining two positive witnesses

The Gaussian orthants defining P(R>T,A=B=C=1) and 11011 have exact crossed-line and straight-line marginals respectively. Use the same normalized tilted halfnormal cubic lower certificate as for the gain, with four variables in the first case and five variables in the second; unlike the gain, the latter root is negative, and the checker passes that sign explicitly. The witness scales and centers are merely rational proposed parameters, not stored output probabilities. Exact directed compilation gives

    y <= 25712328863/10^12; p_cross <= 11174938767/10^12;
    p_0110 <= 12083713264/10^12; p_triple >= 1278710786/10^12;
    p_11011 >= 14548396569/10^12;
    q_initial >= 566886136812/10^12 > 1/2.

This lower bound includes the matching repair.

## 5. Simultaneous local erasure cannot increase the same edge-split energy

On a finite graph let S be the set of positive edges of the initial assignment and E(S)=sum_v NAE on its incident bits. E is a submodular set function (each NAE is submodular). For a patch B of allowed edge variables choose the least minimizer X_B among subsets of S agreeing with S outside B. It exists by meet/join closure. Its removal set is C_B=S minus X_B. For any union D of previously selected removals, submodularity gives
```
E(S minus (D union C_B)) <= E(S minus D)+E(S minus C_B)-E(S minus (D intersection C_B)).
```
The last set is feasible on the B patch: it differs from S only by deletions inside B. Thus the last two costs have nonpositive difference. Induction proves that the sum of split indicators cannot increase for **any simultaneously activated subcollection**, with arbitrary overlap. All patches are minimized against the same original S, not sequentially recomputed.

At each anchor take a Q-ball of the radius below. Its patch may erase internal matching and C edge bits. Independent uniform marks activate patches with probability λ. Intersect the least minimizers and take the physical majorities. Transfer the finite-radius split bound through a finite marked 4-regular graph of sufficiently large girth. Extend the partial right translations of the reduced words of length at most L to permutations. The group has no nontrivial relation of length at most L by evaluation at the empty word. Open each Cayley vertex into two physical vertices joined by M.

---

### Rigorous upper bound on positive majority after all patches activate

The patch at a root anchor frees the internal edges of its ball. Extend each boundary C ray straight through 5 more physical vertices and free the C edges that lead to them. The extensions are disjoint on the matched cubic tree. Each matching edge at an extension vertex and the far C edge stay fixed.

Relax the patch by dropping the interior C ceilings, keeping the zero pin at each raw-negative interior matching edge, setting each extension matching edge to its raw bit, and setting the far C edges to one. These fixed bits are at least the actual bits because the matching repair and later operations only delete positive edges. Outside terms depend only on the fixed variables and are constant within each problem. The meet of the two least minimizers is feasible for the actual patch, and their join is feasible for the relaxed patch. Submodularity gives that the relaxed least minimizer is at least the actual one. The latter bounds the fully patched root bits from above.

Compute that relaxed tree solution exactly. A directed subtree reports its zero-cost minus one-cost difference on the incident edge to its parent, clamped to {-1,0,1}. At a cubic NAE vertex with two incoming messages x,y the outgoing difference is sign(x+y) after clamping the two incoming differences; this follows directly by minimizing four cases, equivalently paying NAE on the two child bits and parent bit. When a boundary edge is fixed 1 or 0 the corresponding infinite difference may safely be clamped: changing that bit affects at most one adjacent NAE by one, hence all parent fixed-bit minimum-value differences are unchanged. Combining the two physical endpoints at a matching anchor, a site with straight incoming x and perpendicular pair y,z sends
```
sign(x-1) if raw X<=T,          sign(x+sign(y+z)) otherwise.
```
A C extension vertex with fixed matching bit b sends sign(x+2b−1) toward the root. Start at the far fixed-one message +1. Write L_k for k straight steps. Each map is increasing and sends +1 to at most +1, so L_(k+1)<=L_k pointwise. With L_k on the child rays, the full update is at most one straight step because sign(y+z)<=1. Hence L_5 is a supersolution and the true coupled laws decrease at each full depth. With d full layers and 5 straight steps the root is computed inside a ball of radius at most d+5. At the root use its four inward messages z_0,...,z_3 (0,1 and 2,3 paired), b=1[R>T] meaning unpinned M, and compute
```
m = b AND [sign(z0+z1)+sign(z2+z3)>0],
e_i = [z_i+sign(z_(i xor 1)+(b?sign(z_(i xor 2)+z_(i xor 3)):-1))>0].
```
These are the simultaneous least optimum bits: an unfixed coordinate is 1 in the least minimizer iff its minimum with 1 is strictly cheaper than with 0. Each root C-edge inward difference is already in {-1,0,1}. Thus the target is the mean of (majority(e0,e1,m)+majority(e2,e3,m))/2, not a message density, an M density, or a number of split vertices. Both bits and the target are monotone. The independent terminal-root audit enumerates the 162 ternary root cases; the overlap argument is proved in §3.

### 5.1 Directed graph calculation

The trit law increases in D,P by the positive Gaussian transitions and monotone updates. Store lower negative tail m and upper positive tail p. Round D up, then round P up on the resulting grid line. Missing queries use the maximal trit.

Let gamma=(t−a)/t=13/41, beta=2045/1148, and Y=P−gamma D. Write x=tD+2a²W, z=aD−atW, u=tx, and v=2az. Then u+v=D, x²+2z²=D²+2a²W², the two child second coordinates are x+tY and beta z+aY, and the raw pin score is Y+(t−a)v/(2a²t). Auxiliary first coordinates sum to 2z. A bin endpoint x_i looks up contrast Y+x_i−gamma D_i after rounding its first coordinate to D_i.

The stored first-coordinate step is h=1/d before striding. Bins have width h/2. At main index k, their indices sum to 2k+1. For bin midpoints u_i,v_j the lower mass is
    (ell/sqrt(pi)) exp(D²/(2a²)−ell²/12) exp(−u_i²/(2a²t²)) exp(−v_j²/(4a⁴)),
where ell=(h/2)/(2a²t). Jensen gives this lower mass by integrating exp(−W²) over an interval of width ell and using the mean of W². For auxiliary index k, the bin indices sum to 4k+1. With midpoints b_i,b_j and ell=(h/2)/a, the same argument gives the mass
    (ell/sqrt(pi)) exp(z²/a²−ell²/12) exp(−b_i²/(2a²)) exp(−b_j²/(2a²)).
Use upper bin endpoints and assign omitted mass to the maximal output.

Use the degree-11 Taylor lower bound at b/256 and eight downward squarings. The exponential, normalization and tail scales are 10752, 2^40 and 16384. Split the bin containing the pin threshold by adding a lower bound on its pinned subinterval. Its two convolutions use the exact negative and nonpositive differences between fixed negative and unpinned matching preference. Check for carries before packing in base 2^64. Split multiplication by output residue.

Start each boundary with five straight sites, as in §5. At a pinned site the negative event is that its straight child is at most zero. At an unpinned site the positive event is that its child is at least zero. These formulas give the line transfer. Use 20 sweeps on d=212, with 373 contrast rows and stale rows as allowed by monotonicity in depth. The prefix checks both reconstructed seed tails pointwise against this result. Lift the same seed to d=306 by rounding up both physical coordinates as above and use two fine sweeps on 509 contrast rows. Earlier bounds remain valid at any greater depth by the collar comparison.

Join the start u and decoded endpoint v in lower-tail coordinates (m,1−p). Its auxiliary output has three lower quadratic Bernstein coefficients, obtained by the two same-input and one cross-input convolution. The main output combines these with u and v to form four cubic Bernstein coefficients. Subtract the chord and check strict positivity wherever v exceeds u. The other coordinates are at most their verified start bound. Exact midpoint subdivision to depth nine proves a residual of at least 1/(3 PATH_DEN 8^9) in table units, with PATH_DEN=441×2^72. Each slope is at most 16384, so 3 PATH_DEN 8^9×16384 physical levels suffice. The radius also includes five collar sites and at most 20×373+2×509 levels for the sweeps.

For the root, use independent unit Gaussians R,U and symmetry to restrict U to [0,4]. Use R in [−4,4]. Each line uses its largest coordinates in each rectangle and one auxiliary pass to bound the ordered categories --,-0 | -+,00 | 0+ | ++. If the ordered upper tail triples are (A,B,C) and (D,E,F), the positive endpoint expectation is A+D at a pinned root and 2(CE+BF−BE)+A(1−F)+(1−C)D otherwise. Check the payoff against the literal minima and nine simplex vertices. Assign omitted mass the maximal payoff. Directed Gaussian intervals and root steps 480/40000, doubled beyond magnitude 2, give 8161/16384, strictly below one half.

### 5.2 Exact fair marginal and cut density

The fixed terminal patch with all activations yields majority strictly below half by the relaxed outcome. Apply the combined repair calibration in the final addendum of proof.md, from the pre-motif state of §4. This yields a fair local rule with expected split cost below c. Strict margins persist at finite first-pass and shell truncations, and transitivity gives the marginal at each vertex. See proof.md for the edgewise transfer.

On each physical edge where its endpoints' majorities differ, its assigned bit agrees with at most one endpoint and is the unique minority edge at the other. Distinct crossing edges cannot charge the same split endpoint, since it has only one minority incidence. This proves that expected crossing count per vertex is bounded by expected split cost, hence strictly below c. This finishes the graph input.

The pair checker proposes factors and shifts by 32 floating-point coordinate sweeps. The exact Jensen check validates the rounded rational proposal. The cylinders and their 5,976 trit cases are unchanged.
