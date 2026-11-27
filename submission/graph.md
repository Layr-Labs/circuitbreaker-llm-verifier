# Bounded-local fair cuts approaching the certified energy ceiling

On the cubic tree with a specified perfect matching, there is an equivariant bounded-radius rule with independent labels, fair at each vertex, with expected crossing count per vertex below 12233725047350865816/10^20, counting each edge once. Section 1 of proof.md transfers the edgewise bounds to finite graphs.

The graph parameters are a=28/57, t=41/57 and T=-293/800. Both directed replay phases in Section 5.1 are required.


## 1. A globally compatible Gaussian source, including branched neighborhoods

Contract the perfect matching to anchors of a loopless 4-regular multigraph Q. Pair its four labeled incidences into two straight continuations, one per physical vertex. The marked universal tree has two lines at each anchor. Set t=41/57, a=28/57, delta=t-a, so t^2+2a^2=1.


At each oriented edge update inward messages with weights t for the straight predecessor and a for the two perpendicular predecessors. Seed at the tail with iid standard normals. At depth j the message has variance one; shells and disjoint half-trees are independent. Sum depths j<J divided by sqrt(J). The recurrence error has L² norm sqrt(2/J) (two boundary shells). A covariance subsequence gives the exact equivariant Gaussian law. The anchor score is half its four independent inward modes. Conditional on entering B, fresh independent standard Z,Y give B_f=tB+sqrt(2)aZ, W_0=sqrt(2)aB−tZ, B_l,r=(W_0±Y)/sqrt(2); distinct future half-trees stay independent given their entering modes.

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
Beyond these states the three half-trees use independent innovations. With entering B, D=aB. The forward and backward modes are tB+2aW and 2E−tB. The perpendicular pair sums to 2(D−tW), and the other pair to 2(E+aW). The definition of E weights its own pair by t/2 and the other by a/2, giving the displayed updates. Finite-shell threshold comparisons converge on any fixed ball since each equality has probability zero.

### 1.1 Line and crossing covariances

The first straight ray has coefficients a/2 and (2a+t)/2 on R in D and P. The independent W and straight update send coefficients (x,y) to (tx,ty+ax). Thus the line covariance is r(d)=(a+t/2)t^(d−1)+(d−1)a²t^(d−2)/2 at distance d≥1, with the second term zero at 1 and r(0)=1. The first R and next fresh W enter triangularly with nonzero diagonal, proving positive definiteness on each line block. The U coefficient is ell[r(d)−h(d)], h(d)=a t^(d−1), by the same update starting at ell t/2 in P. Opposite lines share only R,U, so their distance-i,j covariance is h(i)r(j)+h(j)r(i)−h(i)h(j).

## 2. First-pass energy

A split vertex covers an edge of the raw disagreement graph if one of that edge's ends is split. The optimum equals the minimum vertex-cover size by the argument in §4. Each component is a path, so its minimum cover equals its maximum matching, obtained greedily from either end. With disagreement indicators d_i and selection g_i=d_i(1−g_(i−1)), expand 19 times and drop the nonpositive remainder. For j successive disagreements of probability A_j(T), stationarity gives
```
c1 ≤ S2(T)+A3(T)−A4(T)+A5(T)
      −(A6−A7)(T)−(A8−A9)(T)−(A10−A11)(T)−(A12−A13)(T)−(A14−A15)(T)−(A16−A17)(T)−(A18−A19)(T),
S2(T)=P(X0>T,X2≤T)=A1−A2.
```
The two consecutive disagreements have XOR equal to endpoint disagreement, and the latter has equally probable orientations. The boundary error is confined to the first disagreement component, which is almost surely finite by §4. The least optimizer has straight trit update f(x,b)=sign(x+2b−1), obtained by minimizing the cubic NAE cost. Two equal pins reset the trit. Incoming trits x,y give C bits [x+f(y,b)>0] and [y+f(x,b)>0]. Take finite ray and shell truncations using these resets.

### 2.1 Directed sign moments

For an alternating word on m sites plus its sign reversal, the centered probability is 2^(1−m) times the sum over even I of (product of alternating signs on I) E(product of Gaussian signs on I). Thus S2(0), A3(0), A4(0), and A5(0) follow from moments of at most six sites. For a centered Gaussian matrix G, scale the covariances incident to the upper central index i by v in [0,1]. The even moment vanishes at zero. Differentiation in the covariance and two spatial integrations by parts give the derivative contribution for j≠i
```
(2/pi) G_ij M_remaining(v) / sqrt(1−v²G_ij²).
```
Here M_remaining uses the covariance conditional on i,j at zero. Its entry k,l is
```
G_kl−G_kj G_lj−(v²/(1−v²G_ij²))(G_ki−G_ij G_kj)(G_li−G_ij G_lj).
```
Normalize its diagonals. The terminal two-sign moment is 2 asin(rho)/pi by angular measure.

For the integration error, exact rational elimination gives inverse diagonals ≤20/3 for the six-site covariance. Principal submatrices inherit this bound by the inverse-diagonal variational formula. It survives disconnecting a center, convex interpolation of that disconnect, and normalized Schur conditioning, since the inverse of a Schur complement is a principal inverse block and conditional variances are ≤1. Whiten the fixed block and rotate. The two varying eigenvalues are 1±sv, where s²=1−1/(G^-1)_ii≤17/20, and the other eigenvalues are one. The sign integrand remains a fixed function of absolute value ≤1. On the ellipse of parameter 7/4 mapped to [0,1], Re v lies in [−9/112,121/112] and |Im v|≤33/112. Thus each varying eigenvalue has real part ≥1/256, one has real part ≥1, and both have imaginary part at most 9/32 in absolute value. Their ratios |lambda|/Re lambda are at most 73 and 17/16. The complex Gaussian absolute integral is their square-root product, below 9. Its normalized second moment in each coordinate is |lambda|²/Re lambda. Differentiating the density bounds the combined derivative by (9s/2) sum(1/|lambda|+1/Re lambda)≤9s(256+1)<4096. Positive real inverses give analytic domination. The argument also covers the nested four-sign integrals.

Use GL32 on [0,1], exact through degree 63. The polynomial recurrence checks 32 disjoint opposite-sign root brackets. At an exact root x of P32 the mapped weight is (1−x²)/(32² P31(x)²) by the derivative identity. Enclose quadrature with positive interval weights at 2^-72. For the analytic derivative the kth Chebyshev coefficient is at most 2·4096·(4/7)^k by the ellipse Laurent bound. The quadrature error is at most twice the resulting uniform tail past degree 63, since both the integral and positive normalized rule have norm one. Machin's identity encloses pi, integer roots give outward square roots, and two half-angle reductions and alternating arctangent give the arcsine. At each level add the integration error after summing the contributions.

### 2.2 Shift and longer tails

Let u=293/800 and kappa²=(1−r(2))/(1+r(2)). Conditioning at both boundaries gives −S2'(s)=phi(s)(2Phi(kappa s)−1). Integrate exp(−(s²+w²)/2)/pi over 0≤s≤u, 0≤w≤kappa s. Set w=kappa s v and use the odd degree-five exponential lower polynomial, a lower kappa, and an upper pi.

For m=4,5,6 put C_ij=r(|i−j|), K=C^-1, l=K1, D=1'K1, Z=l'Y for Y~N(0,C), and O the two alternating orthants. Density translation gives A(T)=exp(−D T²/2) E[1_O cosh(T Z)]. For one alternating sign vector s, let F(S,k) integrate Z^k over the signs left after pinning S at zero, with the marginal pinned density included. With J its complement, set
```
v_j=1−C_(j,S) C_(S,S)^-1 1,   d=D−1'C_(S,S)^-1 1.
F(S,k)=(k−1)d F(S,k−2) + sum_(j in J) s_j v_j F(S union {j},k−1).
```
Use zero correction for empty S and omit the first term for k=1. This is integration by parts using (l_J'y)g(y)=−sum_j v_j partial_j g(y), with boundary density included exactly. Thus M_(2k)=E[1_O Z^(2k)]=2F(empty,2k). At k=0 expand conditional orthants into sign moments and divide by (2pi)^(|S|/2) sqrt(det C_SS); nonempty S then has at least two pins and leaves at most four variables. The second-moment control checks overlap with a separate interval for D A(0)+A''(0), calculated by conditional boundary derivatives:
```
A''(0)=sum_(i≠j) s_i s_j sqrt((1−C_ij)/(1+C_ij))/pi
         ·P(remaining signs | Y_i=Y_j=0).
```
On face i the other conditional thresholds move at rate 1−C_ji, and the marginal density derivative vanishes at zero. The checker counts 62 contributions.

Whiten Y to radius R and direction. On O, Z²/R²≤beta. To get beta exactly let M=diag(s)K diag(s), v=diag(s)K1. Maximize the two signed linear forms ±v over y≥0, y'My≤1. On any nonzero maximizing support S the multiplier gives a positive multiple of M_SS^-1(±v_S). Retain nonnegative solutions, normalize in M_SS norm, and maximize the squared values. Radius and direction are independent and R² is chi-square with m degrees of freedom. Consequently M_(2k)≤M_(2l) beta^(k−l) prod_(j=l to k−1)(m+2j). The replay uses l=6. For m=4,6 bound the nonconstant cosh sum by its first six upper terms and the decreasing geometric tail with ratio T² beta(m+2l)/((2l+1)(2l+2)). Write H for this upper sum and e+ for the degree-six upper exponential. Then A(T)−A(0)≤(e+−1)A(0)_lo+e+H. For m=5 use (1−e−)A(0)_hi−e− times the lower nonconstant sum, with positive degree-five e−.

For m=8,10,12,14,16,18,20 the discount is the two opposite words with m−1 alternating bits and the last pin repeated. Let s be the signs of one such word, v_i>0 the tabulated scales, K=C^-1, D=1'K1, and independent Y_i=v_i|N(0,1)|. Set E=E[Y' diag(s)K diag(s)Y], L=sum_i s_i v_i sum_j K_ij. Changing variables to Y in each orthant and applying Jensen gives the sum lower bound
```
prod(v_i)/(2^(m−1)sqrt(det C)) exp((m−E−D T²)/2) cosh(T sqrt(2/pi)L).
```
Evaluate E from the halfnormal first two moments with outward signed products; use the odd degree-11 exponential lower polynomial and cosh(x)≥1+x²/2. The input scales are rational witnesses. The exact determinant and inverse are checked.

For the first eight-site discount also use both tilted halfnormal charts of gain.py at centers −3/5000,−2767/10000. With L the exponent after the change of variables, coefficientwise floor L−z to denominator 2^44 (all monomials nonnegative) and expand the degree-five exponential minorant. Raw halfnormal moments are as in §3. Group by odd exponent count and choose the directed endpoint of sqrt(2/pi) according to coefficient sign. Multiply by a positive degree-15 lower bound for exp(z), divide by the upper determinant square root, and add the two charts. The sum exceeds 37311636/10^11. Retain the larger directed lower discount.

## 3. Tight gain at the active threshold

Take T=-293/800, a=28/57, t=41/57. A raw pin is 1 when its score exceeds T. At one matched root, use the thresholded root R and four consecutive straight rays A,B (opposite directions of one line), C,D (opposite directions of the other). On a ray, truncate to six pins with largest external message +1. A pin b=1 sends 1 if its incoming message is not -1 and 0 otherwise; b=0 sends -1 if its incoming message is not +1 and 0 otherwise. This recursion is monotone in the incoming message. Write M, H, Z for the events that the resulting upper message is respectively -1, nonpositive, and exactly zero. Thus M and H imply the corresponding downward events for actual messages, though Z alone need not. We use Z for partitions of the *truncated* event only.

For two rays of the same endpoint at a positive root, the first-pass line optimization gives local improvement +1 when both messages are nonpositive and at least one is negative (Y), and nonnegative improvement when at least one is negative (N). These identities follow by substituting the 3 x 3 incoming trit pairs into the two adjoining edge bits [x+f(y,1)>0], [y+f(x,1)>0]. Since the two physical endpoints can change their shared matching bit only together, the nonnegative improvement for a pair equals 1[Y_1 N_2]+1[Y_2 N_1] if the root bit is 1. On dividing by two vertices and using line exchange, it suffices to lower-bound P(R>T, Y_AB, M_C or M_D). Starting with the largest boundary +1 makes this a lower bound for every deeper boundary. Other gains are not needed.


Use §1.1 for the marked covariance. To bound p=P(R>T,Y_AB,M_C), split Y_AB into (M_A,H_B) and (Z_A,M_B). Merge six-pin prefixes with equal trits to get 72 disjoint products, 54 under reflection. Stored disjoint eight-pin charts add cases that fail the six-pin test. Upper-bound the overlap of the two reflected events by §3.3. A negative ray needs a first zero and a zero among pins two and three; a nonpositive ray needs a zero among its first two. Thus the overlap enclosure still holds for the extra classes.

### 3.2 Signed Gaussian cylinders

For covariance G, precision K, signs s and positive scale vector v, set U=diag(s_i v_i), z=T1. Under independent halfnormals H use X=z+UH. The cylinder mass is prod(v_i)/(2^n sqrt(det G)) times E exp(L), where L=H'BH+l'H−z'Kz/2, B=(I−U'KU)/2 and l=−U'Kz. One-edge cylinders may instead use z=T1+Db and U=DA for sign diagonal D, b≥0 and nonnegative lower triangular A with positive diagonal. The resulting suborthant uses det A for the Jacobian and the same exponent. Replay uses rational witness scales, offsets and centers.

Let mu=E H, eta=H−mu, m_j=E eta^j, w=m_2. The raw moments start at 1,sqrt(2/pi) and satisfy a_j=(j−1)a_(j−2). Thus EL=w tr(B)+mu² sum B+mu sum l−z'Kz/2. Write centered L as sum of J_i=(2mu sum_j B_ij+l_i)eta_i+B_ii(eta_i²−w) and e_ij eta_i eta_j with e_ij=2B_ij. Its variance is sum EJ_i²+w² sum e_ij². Its third moment is the sum of EJ_i³, 6e_ij E(J_i eta_i)E(J_j eta_j), 3w e_ij²[E(J_i eta_i²)+E(J_j eta_j²)], m_3² e_ij³, and 6w³ e_ij e_ik e_jk in their respective one-, two-, and three-index sums. Independence removes every other term (a lone eta).

For center s put Y=L−EL−s, variance v of L and third centered moment k. The cubic lower bound times exp(EL+s) uses mean 1−s+(v+s²)/2+(k−3sv−s³)/6. Add the nonnegative remainder R(Y)=exp(Y)−1−Y−Y²/2−Y³/6. For t≥0 write r(t)=R(−t), f(t)=r'(t)/t. Direct differentiation gives f'(t)=integral_0^1 u(1−exp(−tu))du, increasing and concave. For b≤sqrt(EY²) put B(t)=r(t)−r(b)−f(b)(t²−b²)/2. For t≤b, integrating f'(z)≥zf'(b)/b gives B(t)≥f'(b)(b²−t²)²/(8b); for t≥b, B(t)≥bf'(b)(t−b)²/2. Hence for y≤0 the tangent surplus in R(y) is ≥kappa(y+b)², kappa=bf'(b)/8. For y≥b, R(y)−r(y)≥y⁵/60 gives surplus ≥b³(y+b)²/240. For 0≤y≤b set x=y/b and lambda=b³/(240kappa), and use surplus ratio ≥kappa[(x−1)²+lambda x⁵]. The tangent x⁵≥5a⁴x−4a⁵ yields lower factor 1−4lambda a⁵−(5lambda a⁴−2)²/4. Maximize over the rational anchors and zero and cap at one. The same factor applies for y≥b by evaluating it at x=1. Thus if K is that factor times kappa, E R(Y)≥r(b)+K(EY²−2bs+b²), dropping the nonnegative mean tangent term.

Replay brackets b, r(b) and f'(b) with directed roots and odd degree-15 exponential lower polynomials. It can use f'(b)≥b/[3(1+b)] from 1−exp(−u)≥u/(1+u). Also the integral remainder and survival (1−u)^4≤exp(−4u) give E R(Y)≥(EY²)²/[6(4+sqrt(EY²))] by Cauchy–Schwarz. Use the stronger directed estimate and Jensen at zero center. Suggested centers, pivots and chart parameters are not input probability claims.

Five depth-six cylinders use the polynomial method of §2. Four eight-pin marginals use order five and one six-pin marginal uses order seven. Round the coefficients of the centered tilted log density downward. Its halfnormal monomials are nonnegative. The odd Taylor polynomial then bounds the exponential below. Group product moments by odd-coordinate count and round by coefficient sign. Degree at most fourteen in each coordinate permits four-bit exponents. Cylinder membership and multiplicities are unchanged.

### 3.3 A quantitative overlap upper

The desired probability is 2p-q. For q, take either orientation of a negative ray A and a nonpositive ray B on the first line. Section 3.1 bounds each orientation by the sum of the 16 marginals on the positive root and negative pins A1,Aa,Bb,C1,Cc,D1,Dd, where a,c,d are 2 or 3 and b is 1 or 2.

On each orthant use signed deviations y>=0, signed precision M, linear density coefficient b=T diag(sign) K1, and C=T² 1ᵀK1. Check b>0, put b'=b−T e_0, V=C−T², and r=b'ᵀy. The identity b'ᵀM^-1 b'=V and Cauchy–Schwarz give yᵀMy>=r²/V. On the uniform simplex at radius r, the root fraction Z=b'_0 y_0/r has beta(1,7) law. The linear exponential factor averages to exp(−r) E exp(gamma r Z) with gamma=−T/b'_0 in (0,1). Since E Z^j=j!/(8⋯(7+j)), the latter expectation is at most 1/(1−gamma r/8) for r<=5. On each of 192 panels bound the decreasing exponential at the left endpoint and this reciprocal at the right. The radial mass is (h^8−l^8)/8! and the density divisor is (2pi)^4 sqrt(det G) prod b'_i. For r>=5 use Z<=1 and the bound exp(−5(1−gamma)−v)V^4(1+v+v²/2+v³/6)/105, v=25/(2V), from the gamma integral. Multiply by exp(−C/2). Use reciprocal degree-12 exponential majorants and outward arithmetic throughout.

For each orientation order the 16 events parent-first along the replay tree. The marginal sum overcounts by at least the sum of the 15 intersections with parents. Each intersection pins the union of the two sets and uses the same signed chart bound. Reflect to get the second orientation.

Scale each rational Gram matrix to integers. Use fraction-free Gauss–Jordan with positive pivots and verify the adjugate identity. Integer square roots, Machin's formula for pi, and outward products bound the numerical factors.

## 4. Initial majority


Minimum split sets are minimum covers of the paths of discordant raw pins (§2). On an odd path the unique minimum cover is the alternating internal vertices: selecting an endpoint leaves a matching of the same minimum size. On an even path the positive vertices form a minimum cover and allow all discordant C bits zero, as in the least optimizer. Isolated vertices are unsplit. Repeats occur on both rays almost surely: their positive mean and covariance decay imply variance tending to zero for their empirical average. This justifies finite ray truncation.

At a raw-positive endpoint let A,B indicate negative nearest pins on its line and C,D on the partner line. If both its first-pass C bits are zero, charge AB or an adjacent positive agreement edge set to zero. In the latter case both ends are split, forcing pattern 0110 in one of two orientations. If its C bits are mixed, exactly one of A,B is one. Matching deletion then requires C or D, since two partner C ones make retaining M strictly better. The inequality (A+B−2AB)(C+D−CD)≤AC+AD+BC+BD−2AB(C+D) gives loss at most 4p_cross−4p_triple, where p_cross=P(R>T,A=C=1) and p_triple=P(R>T,A=B=C=1). The raw-negative pattern 11011 forces the two root C ones by the odd-path rule, giving disjoint recovery. Thus initial majority is at least
```
Phi(−T)−y−2p_0110−4p_cross+4p_triple+p_11011,
y=P(R>T,A=B=1).
```
For y and p_cross condition on R=z. The remaining two coordinates have mean rz, variance v=1−r², correlation rho=(r_neighbor−r²)/v, with r=97/114, r_neighbor=529/722 or 2ar−a². Their joint lower tail is the expectation over S=|N(0,1)| of Phi((T−rz)/sqrt(v(1+rho)/2)−sqrt((1−rho)/(1+rho)) S). This decreases in z,S on z≥T. Use rectangles of width ≤1/64 up to five, lower endpoints, directed upper masses, and a union-tail bound. Normal CDF uses midpoint 1/1024, |phi''|≤1 and monotonic directed brackets between grid points.

For 0110 at zero use the signed even moment expansion in §2.1. Its derivative is the sum over faces i of −s_i phi(0) times the conditional three-sign probability. Compute that probability by the three two-sign moments with covariance C_jk−C_ji C_ik. The derivative is positive. Prékopa log-concavity in the common threshold then gives p_0110(T)≤p_0110(0)<.0121. Use the halfnormal certificate in §3 for the two positive witnesses, with opposite root sign for 11011. The checker obtains:
```
y ≤ 25662444692/10^12, p_cross ≤ 11142921838/10^12,
p_0110 ≤ 12083713264/10^12, p_triple ≥ 1272776606/10^12,
p_11011 ≥ 14564912183/10^12,
initial majority after the matching repair ≥ 568150426130/10^12.
```
The additional motif debit leaves at least 542487981439/10^12. The repair calibration uses the pre-motif margin. The strict margins persist on finite neighborhoods by continuity of the Gaussian threshold distributions.

## 5. Simultaneous local erasure cannot increase the same edge-split energy

Let S be the positive edge set and E(S)=sum_v NAE. For a patch B choose the least minimizer among subsets of S agreeing outside B, which exists by submodular meet/join closure. Write C_B for its removal set. For any union D of earlier removals,
```
E(S − (D union C_B)) ≤ E(S−D)+E(S−C_B)−E(S−(D intersection C_B)).
```
The last set is feasible for patch B, so the last two costs sum to a nonpositive quantity. Thus independently activated patches, each minimized from the original S, can be intersected without raising cost. Transfer through a large-girth marked graph by extending the partial right translations on reduced words of length ≤L to permutations. Evaluation at the empty word excludes any nontrivial relation of length ≤L. Split each anchor along its marked M edge.

### Relaxed terminal majority

Free the internal edges of an anchor ball. Extend each boundary C ray through three straight physical vertices, freeing the joining C edges. Extensions are disjoint on the tree. Drop interior C ceilings, keep raw-negative interior M pinned zero, fix extension M edges to their raw bits and far C edges to one. The fixed bits dominate the actual ones. The meet of actual and relaxed minimizers is feasible for the actual patch and their join for the relaxed patch. Outside terms are constant in each problem. Submodularity and leastness then bound the actual root bits above by the relaxed least minimizer.

A ray reports the zero-cost minus one-cost difference on its parent edge, clamped to {-1,0,1}. Changing a child bit affects its adjacent NAE cost by at most one, which justifies clamping fixed-bit and other differences. The four child choices give sign(x+y) at one cubic vertex; joining the two ends of M gives sign(x−1) when pinned and sign(x+sign(y+z)) otherwise, for straight x and perpendicular y,z. On an extension use sign(x+2b−1). Let L_k start with +1 and take k extensions. Monotonicity gives L_(k+1)≤L_k and full update of L_k≤L_(k+1). Thus the true laws from L_3 decrease in each depth, and d full levels require radius at most d+3.

At the root use its four trits z_i, pairs (0,1) and (2,3), and b=[R>T]. The least optimum bits are
```
m = b AND [sign(z0+z1)+sign(z2+z3)>0],
e_i = [z_i+sign(z_(i xor 1)+(b?sign(z_(i xor 2)+z_(i xor 3)):-1))>0].
```
A free bit in the least optimum is one exactly when the cost at one is strictly smaller than at zero. The two physical majorities are those of (e0,e1,m) and (e2,e3,m), monotone in these inputs. The checker exhausts the 162 literal cases.

### 5.1 Directed graph calculation

The trit law increases in D,P by the positive Gaussian transitions and monotone updates. Store lower negative tail m and upper positive tail p. Round D up, then round P up on the resulting grid line. Missing queries use the maximal trit.

Let gamma=(t−a)/t=13/41, beta=2045/1148, and Y=P−gamma D. Write x=tD+2a²W, z=aD−atW, u=tx, and v=2az. Then u+v=D, x²+2z²=D²+2a²W², the two child second coordinates are x+tY and beta z+aY, and the raw pin score is Y+(t−a)v/(2a²t). Auxiliary first coordinates sum to 2z. A bin endpoint x_i looks up contrast Y+x_i−gamma D_i after rounding its first coordinate to D_i.

The stored first-coordinate step is h=1/d before striding. Bins have width h/2. At main index k, their indices sum to 2k+1. For bin midpoints u_i,v_j the lower mass is
    (ell/sqrt(pi)) exp(D²/(2a²)−ell²/12) exp(−u_i²/(2a²t²)) exp(−v_j²/(4a⁴)),
where ell=(h/2)/(2a²t). Jensen gives this lower mass by integrating exp(−W²) over an interval of width ell and using the mean of W². For auxiliary index k, the bin indices sum to 4k+1. With midpoints b_i,b_j and ell=(h/2)/a, the same argument gives the mass
    (ell/sqrt(pi)) exp(z²/a²−ell²/12) exp(−b_i²/(2a²)) exp(−b_j²/(2a²)).
Use upper bin endpoints and assign omitted mass to the maximal output.

Use the degree-11 Taylor lower bound at b/256 and eight downward squarings. The exponential, normalization and tail scales are 11264, 2^40 and 16384. Split the bin containing the pin threshold by adding a lower bound on its pinned subinterval. Its two convolutions use the exact negative and nonpositive differences between fixed negative and unpinned matching preference. Check for carries before packing in base 2^64. Split multiplication by output residue.

Start each boundary with 3 straight sites. With a fixed matching bit zero the negative event is that its straight child is nonpositive. With bit one the positive event is that its straight child is nonnegative. The code uses these formulas and the meet with the previous law. On the d=306 grid use three symmetric index axes with (full radius, last index) equal to (159,367), (171,389), and (196,484) for main, auxiliary and contrast. Use stride five past the full radius, retaining the endpoints, and limits 789 and 966 for the two bin integrations. Advance 10 sweeps in increasing contrast order. Outside the inner contrasts, stale rows only keep previous bounds. Each sweep advances at most 509 levels in the physical recurrence.

Join each consecutive pair of tables u,v in lower CDF coordinates (m,1−p), starting at the table after the sweeps. The auxiliary output has three lower quadratic Bernstein coefficients from both same-input products and one cross-input product, with downward rounding. The main output is bilinear in the straight and auxiliary laws, including the split pin adjustment. Its cubic coefficients are M(u,A00), [2M(u,A01)+M(v,A00)]/3, [M(u,A11)+2M(v,A01)]/3, M(v,A11). Subtract the chord and check positivity when the endpoint coordinate increases. Otherwise use the verified start and the monotone decrease of the true physical laws. Midpoint subdivision to depth 12 gives a residual at least 1/(3 PATH_DEN 8^12) per moving coordinate, with PATH_DEN=121·2^74. Because each positive slope is at most 16384, the path advances by at least one part in 3 PATH_DEN 8^12·16384 per physical level. 2 such lengths plus 10·509+3 bound the patch radius.

Each stored endpoint uses lower-CDF knots in units of two, mixed backward differences, and separable rational cubic interpolation. A mask locates corrections in units of two. The first table adds signed corrections. The second subtracts unsigned corrections. The checker tests both simplex inequalities and hashes the decoded integer tables. The prefix certifies the first path. The terminal shares its first endpoint by digest, verifies the next path, and checks both majority bounds, energy, and both entropy partitions.

For the root, use independent unit Gaussians R,U and symmetry to restrict U to [0,4]. Use R in [−4,4]. Each line uses its largest coordinates in each rectangle and one auxiliary pass to bound the ordered categories --,-0 | -+,00 | 0+ | ++. If the ordered upper tail triples are (A,B,C) and (D,E,F), the positive endpoint expectation is A+D at a pinned root and 2(CE+BF−BE)+A(1−F)+(1−C)D otherwise. Check the payoff against the literal minima and nine simplex vertices. Assign omitted mass the maximal payoff. Directed Gaussian intervals and root steps 1000/40000, doubled beyond magnitude 2, give 8098/16384, strictly below one half.

### 5.2 Exact fair marginal and cut density

The fixed terminal patch with all activations yields majority strictly below half by the relaxed outcome. Apply the combined repair calibration in the final addendum of proof.md, from the pre-motif state of §4. This yields a fair local rule with expected split cost below c. Strict margins persist at finite first-pass and shell truncations, and transitivity gives the marginal at each vertex. See proof.md for the edgewise transfer.

On each physical edge where its endpoints' majorities differ, its assigned bit agrees with at most one endpoint and is the unique minority edge at the other. Distinct crossing edges cannot charge the same split endpoint, since it has only one minority incidence. This proves that expected crossing count per vertex is bounded by expected split cost, hence strictly below c. This finishes the graph input.
