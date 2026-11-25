# Bounded-local fair cuts approaching the certified energy ceiling

On the cubic tree with a specified perfect matching, there is an equivariant bounded-radius rule with independent labels, fair at each vertex, with expected crossing count per vertex below 12282281766462836057/10^20, counting each edge once. Section 1 of proof.md transfers the edgewise bounds to finite graphs.

The graph parameters are a=28/57, t=41/57 and T=-23/64. Both directed replay phases in Section 5.1 are required.


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

Each physical vertex has two C edges and an M pin. Its split indicator is NAE. The minimum split count on a line is the minimum vertex cover size of the raw disagreement paths, as proved in §4. It equals the maximum matching size by adjacent pairing. Greedy selection from left to right is optimal by exchanging the earliest edge into a maximum matching and inducting. For disagreement indicators d_i the greedy selection satisfies g_i=d_i(1−g_(i−1)). Thirteen substitutions give an upper bound since the residual term is nonpositive. Averaging with A_j(T) the j-run probabilities gives
    c1 <= S2(T)+A3(T)-A4(T)+A5(T)
            -(A6-A7)(T)-(A8-A9)(T)-(A10-A11)(T)-(A12-A13)(T).
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

The checker proves disjoint 2^-72 root brackets for P_32 by opposite endpoint signs. There are 32, so all roots are covered. Interval arithmetic verifies positive weights and their sum bracket. It rounds outward at 2^72 and uses integer square roots. Machin's identity gives pi. Two half-angle reductions express arcsine as arctangent, enclosed by alternating series and the next term. Clip correlations to [-1,1] and add the analytic integration error.
For an analytic derivative bounded by 4096 on the ellipse of parameter 7/4, the Laurent coefficient estimate via the Joukowski map bounds Chebyshev coefficient k>=1 by 2*4096*(4/7)^k. Truncating after degree 63 incurs uniform error at most 2*4096*(4/7)^64/(1-4/7). The integral and positive normalized GL32 rule each have norm one, so their difference is at most twice this outward error. Add it at each center-row integration, including each real-node inner four-sign moment. This proof works identically for the higher central interpolation coordinate chosen by the actual code.


### 2.3 Shift + finite-tail discounts

Let u=23/64 and kappa²=(1−r(2))/(1+r(2)). Bivariate conditioning and differentiation at both boundaries give −S2'(u)=phi(u)(2Phi(kappa*u)−1). Thus S2(0)−S2(u) is the integral of exp(−(s²+w²)/2)/pi over 0<=w<=kappa*s, 0<=s<=u. Substitute w=kappa*s*v and integrate the degree-five Taylor lower polynomial. Multiply the positive result by the lower square root of kappa² divided by the upper pi.

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
For m=8,10,12,14, the discount (A_(m-2)−A_(m-1))(T) is the probability of the two opposite words with the first m-1 bits alternating and the last bit equal to its predecessor. Their centered signs are sigma_i=(-1)^i for i<m-1 and sigma_(m-1)=1.

Let K=C_m^{-1}, D=1^t K 1, and Y_i independently halfnormal with positive rational scales s_i. Let M=diag(sigma)K diag(sigma), E=E[Y^t M Y]=sum K_ii s_i^2+(2/pi)sum_{i!=j} sigma_i sigma_j K_ij s_i s_j and L=sum_i sigma_i s_i sum_j K_ij. For the two opposite words translate by T*1 and write the Gaussian integral as an expectation under the product halfnormal density. Each orthant has common normalizing multiplier prod(s_i)/(2^m sqrt(det C_m)); the exponent inside the expectation has mean (m-E-D T^2)/2 +/- T sqrt(2/pi)L. Jensen on each orthant separately and summation give
A_{m-2}(T)-A_{m-1}(T) >= prod(s_i)/(2^{m-1}sqrt(det C_m)) exp((m-E-D T^2)/2) cosh(T sqrt(2/pi)L).
This directly uses the biased threshold instead of applying a separate worst-case translation lower. Use the explicit rational scale vectors in energy.py. They were selected numerically but all claims use only their explicit rational values, not optimality. Bound E from above by rational Machin bounds for pi, with sign-aware off-diagonal division. For any real z the odd Taylor polynomial through degree 11 lies below exp(z) (Taylor remainder nonnegative); verify positive at the chosen lower exponent. Bound cosh below by 1+T^2 L^2/pi_upper, and sqrt(det C_m) above by an authenticated integer square root at denominator 10^24. Exact SPD elimination validates the inverse and determinant.

## 3. Tight gain at the active threshold

Take T=-23/64, a=28/57, t=41/57. A raw pin is 1 when its score exceeds T. At one matched root, use the thresholded root R and four consecutive straight rays A,B (opposite directions of one line), C,D (opposite directions of the other). On a ray, truncate to six pins with largest external message +1. A pin b=1 sends 1 if its incoming message is not -1 and 0 otherwise; b=0 sends -1 if its incoming message is not +1 and 0 otherwise. This recursion is monotone in the incoming message. Write M, H, Z for the events that the resulting upper message is respectively -1, nonpositive, and exactly zero. Thus M and H imply the corresponding downward events for actual messages, though Z alone need not. We use Z for partitions of the *truncated* event only.

For two rays of the same endpoint at a positive root, the first-pass line optimization gives local improvement +1 when both messages are nonpositive and at least one is negative (Y), and nonnegative improvement when at least one is negative (N). These identities follow by substituting the 3 x 3 incoming trit pairs into the two adjoining edge bits [x+f(y,1)>0], [y+f(x,1)>0]. Since the two physical endpoints can change their shared matching bit only together, the nonnegative improvement for a pair equals 1[Y_1 N_2]+1[Y_2 N_1] if the root bit is 1. On dividing by two vertices and using line exchange, it suffices to lower-bound P(R>T, Y_AB, M_C or M_D). Starting with the largest boundary +1 makes this a lower bound for every deeper boundary. Other gains are not needed.

Use the covariance in §1.1 for the root and all ray positions. The straight recursion for (D,E) follows from E=P−D in §1.

### 3.1 Disjoint cylinders, not a sum of overlapping orthants

Let p=P(R>T, Y_AB, M_C). The two disjoint alternatives for Y_AB are (M_A,H_B) and (Z_A,M_B), since H=M disjoint union Z. For each of M,H,Z compute the trit outcomes on all 64 binary words and combine pairs of sibling terminal prefixes. The compressed prefixes within each class partition exactly its set of six-bit strings. Taking the indicated product with the prefixes for M_C yields 72 disjoint cylinder events for p. Reflecting the first line pairs some identical Gaussian probabilities, leaving 54 classes with multiplicities verified from the exhaustive partition. Each cylinder is an orthant with shifted threshold on the root and only those pin coordinates its prefix specifies, never the unconstrained remainder. The tabulated rational scales and expansion centers are *untrusted witnesses*; no fitted probabilities, sample counts or quadrature outputs are read. Any positive scales and any center work if the directed final comparison passes.

The same event p is increasing when the three trits decrease. Extend the rays to eight pins with external +1. Partition each ray into disjoint prefixes classified by its pair of six-pin and eight-pin trits. The checker selects twelve classes of triples that qualify at eight and fail at six. It checks each prefix and its reflection multiplicity. These cylinders are disjoint from each other and from every six-pin cylinder. Add their directed lower masses to the one-ray bound. The old upper bound for overlap still applies. A negative trit needs a first zero and at least one zero among pins two and three, because two ones reset the incoming trit to +1. A nonpositive trit needs a zero among its first two pins. Both facts hold with any deeper boundary.

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
The probability density 4(1-u)^3 on [0,1] has survival (1-u)^4, at most e^(-4u). Its law is stochastically no larger than an exponential random variable of rate four. Use e^(uy) >= e^(-u|y|) and the decreasing function e^(-u|y|) under that smaller law. The displayed integral is at least 1/(4+|y|). Apply Cauchy–Schwarz to Y^2/sqrt(4+|Y|) and sqrt(4+|Y|) for Y=L-E[L]-s. Write M2=E[Y^2]=Var(L)+s^2. Then
```
E[Y^4/(4+|Y|)] >= M2^2/(4+E|Y|) >= M2^2/(4+sqrt(M2)).
```
Use the directed lower M2 in the numerator and the upper square root from the directed upper M2 in the denominator. The extra term is therefore a certified nonnegative addition to the average cubic. For the prefactor use exponential monotonicity at an outward lower endpoint followed by a degree-15 Taylor polynomial, globally below the exponential. Every scale and expansion center is an untrusted rational witness. No optimizer or fitted probability enters the computation.

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
No independence of raw pins, no physical graph finite-size heuristics, and no equality bounds on the finite-window exceptional runs are used. Finite-radius first-pass local approximations preserve strict slack: the positive consecutive-repeat event occurs in each direction almost surely. One direct justification is covariance decay r(d)->0, which gives variance of any size-k average of the repeat indicators tending to zero (Gaussian pairs of such indicators become independent at large lag, by convergence of their finite-dimensional covariances). The event of no repeat on a whole ray is bounded by the event that this average stays zero. The Gaussian joint distributions on every fixed neighborhood are continuous at all threshold tests.

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

    y <= 24279389015/10^12; p_cross <= 10568765751/10^12;
    p_0110 <= 12083713264/10^12; p_triple >= 1289323806/10^12;
    p_11011 >= 14532704536/10^12;
    q_initial >= 569310640684/10^12 > 1/2.

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

For a root edge or vertex, deleting under all active patches can only lower its bits from the center patch solution. That center patch still fixes every originally negative edge and has actual exterior boundary bits <=1. Relax it by **keeping only raw-negative matching pins** (these are truly zero by the one-sided initial repair), dropping all other interior upper bounds, and fixing each C-boundary bit to one except at edges whose first two exterior straight raw pins are zero. At those edges leave the bit free in the relaxed patch. Its actual initial C bit is zero, since its exterior context is −1 and its interior context is at most +1. Motifs only delete. Compare the objective on the physical vertices of the patch. Outside terms are constant under the actual fixed boundary, so do not include them after freeing those bits. Least minimizers can only increase: use submodularity for the two solutions, their meet feasible under the smaller upper bounds and smaller boundary, and join feasible under the larger constraints. Treat boundary as fixed variables in both comparisons. Consequently the relaxed least root majority bounds the fully patched physical majority from above. 

Compute that relaxed tree solution exactly. A directed subtree reports its zero-cost minus one-cost difference on the incident edge to its parent, clamped to {-1,0,1}. At a cubic NAE vertex with two incoming messages x,y the outgoing difference is sign(x+y) after clamping the two incoming differences; this follows directly by minimizing four cases, equivalently paying NAE on the two child bits and parent bit. When a boundary edge is fixed 1 or 0 the corresponding infinite difference may safely be clamped: changing that bit affects at most one adjacent NAE by one, hence all parent fixed-bit minimum-value differences are unchanged. Combining the two physical endpoints at a matching anchor, a site with straight incoming x and perpendicular pair y,z sends
```
sign(x-1) if raw X<=T,          sign(x+sign(y+z)) otherwise.
```
An exterior fixed-one bit gives +1, and a free boundary bit with no exterior cost gives 0. This boundary field is a pointwise supersolution of the directed recurrence. On a ray whose first two raw pins are zero its outgoing update is sign(x−1)≤0, and on every other ray the update is at most +1. The directed fields therefore decrease with depth under one coupling of the raw source. With d layers these incoming messages correspond to exactly the radius-d patch with C edges into the outside fixed. All updates increase in the child messages and in the raw X. At the root use its four inward messages z_0,...,z_3 (0,1 and 2,3 paired), b=1[R>T] meaning unpinned M, and compute
```
m = b AND [sign(z0+z1)+sign(z2+z3)>0],
e_i = [z_i+sign(z_(i xor 1)+(b?sign(z_(i xor 2)+z_(i xor 3)):-1))>0].
```
These are the simultaneous least optimum bits: an unfixed coordinate is 1 in the least minimizer iff its minimum with 1 is strictly cheaper than with 0. Each root C-edge inward difference is already in {-1,0,1}. Thus the target is the mean of (majority(e0,e1,m)+majority(e2,e3,m))/2, not a message density, an M density, or a number of split vertices. Both bits and the target are monotone. The independent terminal-root audit enumerates the 162 ternary root cases; the overlap argument is proved in §3.

### 5.1 Positive Gaussian convolutions from the certified boundary

At state D,P let m lower-bound negative mass and p upper-bound positive mass, with m+p≤1. For independent f,g and zero mass z=1−m−p, the sign sum has negative mass m_f z_g+(1−p_f)m_g and positive mass p_f z_g+(1−m_f)p_g. It is monotone in each child law. The true law increases in original D,P by the nonnegative Gaussian transitions and the monotone threshold update. Round up D first, then choose the smallest stored contrast with physical P at least the required value. Above the grid use the maximal law.

Write gamma=(t−a)/t=13/41 and beta=(a²−t²−at)/(−at)=2045/1148. The main grid stores states D=k h, P=Y+gamma D, with h=1/225, Y=j/225. The k and j axes have the symmetric integer sets specified below. At a parent node let x=tD+2a²W and z=aD−atW denote the straight and perpendicular first coordinates. Define u=t x, v=2a z. Since t²+2a²=1, substitution in §1 yields, identically,
    u+v=D,   x²+2z²=D²+2a² W²,
    P_s=x+tY,   P_o=beta z+aY,    X=Y+delta v/(2a²t).
Consequently for fixed Y the two conditional child distributions live on two fixed lines, and the raw threshold decision depends on v alone. The auxiliary grid tabulates the sign-sum across the two perpendicular rays at z=k h, P_o=E+z, E=j/225, using its own k axis and the same j axis. For a given (z,E), their first coordinates are x_1=z+aV and x_2=z−aV; their second coordinates are x_i+E. We integrate the whole real line, including both ray orientations.

Use symmetric integer axes with (full radius, last index) (117,270) for D, (126,286) for z and (144,356) for contrast. After each full radius use stride four and include the last index. Round the first coordinate up, then round the contrast up to dominate physical P at that stored first coordinate. If none dominates, use the top law. The fine integration mesh does not change these lookups.

In both integrations use the bin spacing H_b=h/2=1/450, while the stored first coordinate spacing remains h. Here is the parent integration rule for fixed D=k h. For each integer i use the bin (i−1)H_b <= u <= iH_b and j=2k+1−i, so (j−1)H_b <= v <= jH_b. Boundary points have zero Gaussian mass. Sample the straight child at the upper endpoint iH_b/t and its associated P_s; sample the auxiliary at the upper endpoint jH_b/(2a) and its associated P_o. Their second coordinates increase with their respective first coordinates. Classify the ENTIRE bin as potentially unpinned if Y+delta jH_b/(2a²t)>T; otherwise the second message is fixed −1. Thus all rounded laws and decisions stochastically dominate their true conditional counterparts. Only |i|,|j|<=580 are included in the parent convolution. The straight ray and the perpendicular aggregate are independent conditional on W, with the auxiliary V integral already inside the latter. Omitted domains and all deficits in the Gaussian mass bounds are assigned the top message.

Write u_i=(i−1/2)H_b, v_j=(j−1/2)H_b. The innovation bin width is ell=H_b/(2a²t). For uniform Z on [−ell/2,ell/2], Jensen gives E exp(−(W_mid+Z)²)>=exp(−W_mid²−ell²/12). Hence the bin mass (density exp(−w²)/sqrt(pi)) has the following positive separable lower bound:
    [ell/sqrt(pi) exp(D²/(2a²)−ell²/12)] exp(−u_i²/(2a²t²)) exp(−v_j²/(4a⁴)).
For an auxiliary node z=k h the two bins have i+j=4k+1, now in coordinates x_1,x_2, both sampled at their upper endpoints. Set |i|,|j|<=710, ell=H_b/a, with midpoints b_i=(i−1/2)H_b,b_j=(j−1/2)H_b. Since (x_1²+x_2²−2z²)/(2a²)=V², the corresponding lower bound is
    [ell/sqrt(pi) exp(z²/a²−ell²/12)] exp(−b_i²/(2a²)) exp(−b_j²/(2a²)).
These are lower bounds per disjoint bin, not Gaussian moment guesses, and missing mass is never renormalized.

Use Q=2^14 for tails, W_s=12288 for exponential factors and 2^40 for the bracketed weights. The degree-11 Taylor polynomial at b/256 lies below exp(b/256) for either sign of b. Check positivity, round down at 2^80 and square eight times with downward rounding. Squaring against an upper Machin pi certifies the lower 1/sqrt(pi). Check total lower mass <=1. Round both integrated lower CDF tails down, using 1−Pr(+) as the second tail. Check the probability simplex after each update.

On the straddling bin, initially treat the entire bin as unpinned. Put v_*=(T−Y)2a²t/δ. The pinned W interval runs from D/t−(T−Y)/δ to D/t−(j−1)H_b/(2a²t). Directed CDF differences lower-bound its mass. Use the same upper child endpoints for both laws. Fixing the second child to −1 adds nonnegative increments to both lower CDF tails. Multiply the unpinned tails by the full bin lower mass and the increments by the pinned lower mass. Add before rounding down. The increment numerators have scale Q² and the mass scale 2^26.

Pack coefficients in base 2^64. Check (2N+1)W_s²<2^64 and Q² max_j C_j<2^64 for weight-product coefficients C_j. With z=1−m−p, the negative and positive probabilities are (1−p_f)(1−p_g)−z_f z_g and (1−m_f)(1−m_g)−z_f z_g. For pinned g use (m,p,z)=(1,0,0). The checked simplex makes each subtraction nonnegative coefficientwise. Each term is at most Q² C_j, with the same weights for pinned bins. No digit carries or borrows. Add the straddling correction outside these products.

For odd index 2r+1 extract coefficient r of f_even g_odd + f_odd g_even. Auxiliary square coefficients at 4r+1 and 4r+3 are twice coefficient r of f_0 f_1 + X f_2 f_3 and f_0 f_3 + f_1 f_2, respectively (lanes collect residues modulo four). Mixed products use both orientations with the same shifts. GNU MP products have Python-integer checks and an exact fallback. One-sided CDF queries and cached packed operands give the same directed integers. Removing and restoring end zero limbs preserves each coefficient.

The depth-zero message is 0 on two exterior zeros and +1 otherwise. The scores are P−δW and tP+aD+2a²W−δW′, for independent variance-one-half W,W′. Intersect 1/16 bins in [−6,6] with W≥(P−T)/δ and use lower mass times the lower W′ tail at the right endpoint. At k/16 its threshold is [41(i+j)/225−57T+(1568/57)k/16]/13. Allocate missing mass to +1 and round up p. True laws decrease in depth. Initialize all auxiliary rows, then update auxiliary and main in increasing contrast order. Auxiliary tails tighten directly. Take the main meet (max m,min p). A row with every pinned correction bounded by the full lower bin mass has a monotone numerical update. Split the weight into nonnegative pinned and unpinned parts. Once such a row equals its meet, later updates tighten it directly. Use periods 1,2,4 as 100|j| passes 45·225 and 80·225. Update when the period divides step+|j|. Stale rows remain valid. Each sweep adds at most the row count to the depth bound.

The mixed backward differences of the endpoint table are stored as zigzag variable-length integers, then LZMA- and Base85-encoded. Invert by adding left and above and subtracting upper-left. Both witnesses store coarse lower-CDF values in two-tick units and dense corrections in whole ticks.

Use the lower CDF coordinates (m, 1−p), in units of Q. Let v1 be the decoded table. In lower CDF coordinates let v0 be the coordinatewise minimum of v1 and the table from 25 ordinary sweeps. This is a weaker upper law at that depth, preserves the simplex and gives v0<=v1. The affine table v(s)=(1−s)v0+s v1, 0<=s<=1, is a simplex law at every stored state. All lookups still use domination in the original physical coordinates D and P.

The exact auxiliary map on v(s) has quadratic Bernstein coefficients A00,A01,A11: apply the auxiliary convolution to the two child tables with the indicated endpoint indices. Either ordering of the mixed pair gives the same sum by interchanging the integration-bin indices. Round both lower tails of each coefficient down. Their Bernstein interpolation then gives a stochastic upper bound on the unrounded auxiliary law. The main update is bilinear in its two child laws. The pinned approximation, Gaussian deficits and threshold correction preserve this property (constant and linear terms can also be written bilinearly on probability simplexes). Use the resulting auxiliary upper law in the directed main integrand. Its cubic Bernstein coefficients, from unrounded main numerators, are: For a mixed packed product replace 2 f_r f_s by f_r g_s+f_s g_r with the same digit shift and carry bound.
    M(v0,A00), [2M(v0,A01)+M(v1,A00)]/3,
    [M(v0,A11)+2M(v1,A01)]/3, M(v1,A11).
For each moving coordinate subtract the cubic Bernstein coefficients of v(s). Check strict positivity of the residual by midpoint subdivision to depth at most 12. All tested coefficients have integer numerators and denominators dividing 3*PATH_DEN*8^12. The residual on the entire path is therefore at least 1/(3*PATH_DEN*8^12) in units of Q. One true iteration advances the parameter by at least 1/PATH_LAYERS, where PATH_LAYERS=3*PATH_DEN*Q*8^12. Indeed each coordinate of v(s) has slope at most Q. Unchanged coordinates retain their previous upper law because the true law decreases with depth. Thus v1 bounds the true law after PATH_LAYERS additional layers.

In the first phase interpolate a packed grid of coarse lower-CDF values with target knot spacing eight and subtract the packed integer corrections. A coordinatewise minimum with the regenerated 25-sweep table starts the exact path to this table. In the terminal phase check that digest and advance eight ordinary sweeps. Form a second candidate by adding coarse lower-CDF increments with target knot spacing 12. Always include the final endpoint of each axis in both grids. Both interpolations are separable cubic Lagrange interpolation with four consecutive knots and exact rational weights and integer rounding. Packed integer corrections finish this candidate. Each correction table stores a support bit mask followed by the selected integers. A fresh cubic verification proves the second path from the eight-sweep table. These constructions have no validity assumption. Evaluate the root at this endpoint. Patch radius is 2*PATH_LAYERS+33*(number of contrasts).

Condition on independent unit Gaussians R,U. Line exchange allows U>=0 with twice the mass. The two centers use D=a(R±U)/2 and P=(97R±41U)/114. One extra auxiliary pass bounds A=Pr(both positive), B=Pr(positive sum), C=Pr(nonnegative sum), with omitted mass assigned to ++. Use R in [−4,4] and U in [0,4], width at most 220/40000 up to magnitude 2 and 440/40000 outside. Use the largest R and each line's own largest signed contrast in a rectangle, then apply the dominating lookup. All tails increase in physical coordinates. The ordered categories are --,-0 | -+,00 | 0+ | ++. Write A,B,C and D,E,F for the two triples of tails at categories >=3,>=2,>=1. The expected positive endpoint count equals A+D when pinned and otherwise equals
    2(CE+BF−BE)+A(1−F)+(1−C)D.
The checker verifies the payoff over 72 unordered cases, 162 literal least root minima and nine vertex cases. Ordered clipping A<=B<=C and the monotone payoff justify inserting upper bounds. Include T as an exact R boundary. Lower Gaussian bin masses use variance 1 at scale 2^26. Allocate omitted mass payoff one after division by two. Sum integer products and round up to Q once.

For density exp(-x²/q)/sqrt(q*pi), q=1,2, use midpoint spacing 1/H, H=4096. The starting exponential is exp(-1/(4qH²)) and ratio i to i+1 is rho^(i+1), rho=exp(-2/(qH²)). Cubic and quartic alternating Taylor brackets initialize both. At scale 2^100 use directed multiplication and Machin pi with integer square roots. Sum at scale 2^50 with outward rounding. The complete and partial errors up to 8 are <=8/(12H²) and <=1/H² by |density''|<=2 and |density'|<=1. Reflect and use Mills tail <2^-32 beyond 8.


The terminal upper is 8187/16384, below one half.

### 5.2 Exact fair marginal and cut density

The fixed terminal patch with all activations yields majority strictly below half by the relaxed outcome. Apply the combined repair calibration in the final addendum of proof.md, from the pre-motif state of §4. This yields a fair local rule with expected split cost below c. Strict margins persist at finite first-pass and shell truncations, and transitivity gives the marginal at each vertex. See proof.md for the edgewise transfer.

On each physical edge where its endpoints' majorities differ, its assigned bit agrees with at most one endpoint and is the unique minority edge at the other. Distinct crossing edges cannot charge the same split endpoint, since it has only one minority incidence. This proves that expected crossing count per vertex is bounded by expected split cost, hence strictly below c. This finishes the graph input.
