# Bounded-local fair cuts approaching the certified energy ceiling

On the cubic tree with a specified perfect matching, there is an equivariant bounded-radius rule with independent labels, fair at each vertex, with expected crossing count per vertex below 12252408923603907442/10^20, counting each edge once. Section 1 of proof.md transfers the edgewise bounds to finite graphs.

The graph parameters are a=28/57, t=41/57 and T=-293/800. Both directed replay phases in Section 5.1 are required.


## 1. A globally compatible Gaussian source, including branched neighborhoods

Contract the perfect matching to anchors of a loopless 4-regular multigraph Q. Pair its four labeled incidences into two straight continuations, one per physical vertex. The marked universal tree has two lines at each anchor. Set t=41/57, a=28/57, delta=t-a, so t^2+2a^2=1.

Give anchor u an independent standard Z_u. On uv put B^0_uv=Z_u and update by weights t on the straight predecessor, a on each perpendicular predecessor. Each B^j has variance one by disjoint predecessor half-trees. The shells on uv use disjoint distance layers, so B^[J]=sum_(j<J) B^j/sqrt(J) has variance one. The update differs from B^[J] by the first and next omitted shells, with L² norm sqrt(2/J). A diagonal covariance subsequence on the countable tree gives a Gaussian limit with the exact recurrence, marked-tree symmetries and independent inward families on disjoint half-trees.

The root score R sums its four independent inward modes and divides by two. Given entering B=t B_f+a(B_l+B_r), use independent standard Z,Y to sample B_f=tB+sqrt(2)aZ, W_0=sqrt(2)aB−tZ, B_l=(W_0+Y)/sqrt(2), B_r=(W_0−Y)/sqrt(2). The linear constraint fixes the backward mode. Farther half-trees remain independent conditional on their entering modes.

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

A split vertex covers an edge of the raw disagreement graph if one of that edge's ends is split. The optimum equals the minimum vertex-cover size by the argument in §4. Each component is a path, so its minimum cover equals its maximum matching, obtained greedily from either end. With disagreement indicators d_i and selection g_i=d_i(1−g_(i−1)), expand 13 times and drop the nonpositive remainder. For j successive disagreements of probability A_j(T), stationarity gives
```
c1 ≤ S2(T)+A3(T)−A4(T)+A5(T)
      −(A6−A7)(T)−(A8−A9)(T)−(A10−A11)(T)−(A12−A13)(T),
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

Whiten Y to radius R and direction. On O, Z²/R²≤beta. To get beta exactly let M=diag(s)K diag(s), v=diag(s)K1. Maximize the two signed linear forms ±v over y≥0, y'My≤1. On any nonzero maximizing support S the multiplier gives a positive multiple of M_SS^-1(±v_S). Retain nonnegative solutions, normalize in M_SS norm, and maximize the squared values. Radius and direction are independent and R² is chi-square with m degrees of freedom. Consequently M_(2k)≤M_6 beta^(k−3) prod_(j=3 to k−1)(m+2j). For m=4,6 bound cosh from above through k=3 and add a geometric tail with decreasing ratio T² beta(m+2k)/((2k+1)(2k+2)), starting at k=3. If H is this upper nonconstant sum, e+ is the positive degree-six exponential upper bound, use A(T)−A(0)≤(e+−1)A(0)_lo+e+H. For m=5, bound A(0)−A(T) above by (1−e−)A(0)_hi−e− sum of the first three nonconstant lower terms, using the positive degree-five e−.

For m=8,10,12,14 the discount is the two opposite words with m−1 alternating bits and the last pin repeated. Let s be the signs of one such word, v_i>0 the tabulated scales, K=C^-1, D=1'K1, and independent Y_i=v_i|N(0,1)|. Set E=E[Y' diag(s)K diag(s)Y], L=sum_i s_i v_i sum_j K_ij. Changing variables to Y in each orthant and applying Jensen gives the sum lower bound
```
prod(v_i)/(2^(m−1)sqrt(det C)) exp((m−E−D T²)/2) cosh(T sqrt(2/pi)L).
```
Evaluate E from the halfnormal first two moments with outward signed products; use the odd degree-11 exponential lower polynomial and cosh(x)≥1+x²/2. The input scales are rational witnesses. The exact determinant and inverse are checked.

## 3. Matching gain

For each of the four straight rays A,B,C,D, let M,H,Z mean its upper trit is negative, nonpositive, or zero after six pins with external +1. On a line at a positive root, putting Y=(both nonpositive, some negative), N=(some negative), the nine trit cases give gain one for Y and nonnegative gain for N. The two physical endpoints share M, so line exchange gives gain per vertex at least P(R>T,Y_AB,M_C or M_D). Exact downstream trits are no larger than the truncated ones, so these decreasing events are valid lower witnesses. The two disjoint Y alternatives are (M_A,H_B) and (Z_A,M_B) on the truncated trits. Sibling-prefix merging of all 64 words gives 72 disjoint cylinders for p=P(R>T,Y_AB,M_C), or 54 reflected classes. The twelve depth-eight extensions qualify at eight and fail at six. The checker classifies each ray by both trits, verifies disjointness and reflections, and adds the new cylinders to p.

For a selected covariance G and signs s_i, set K=G^-1, U=diag(s_i v_i) for positive scale witnesses v_i, and X=T1+UH for independent standard halfnormals H_i. The mass is prod(v_i)/(2^n sqrt(det G)) E exp(L), with
```
L=H' B H+l'H−T²(1'K1)/2, B=(I−UKU)/2, l=−TUK1.
```
Use mu=sqrt(2/pi), eta_i=H_i−mu, w=E eta_i², m3=E eta_i³, b_i=2mu sum_j B_ij+l_i, d_i=B_ii, e_ij=2B_ij, and F_i=b_i eta_i+d_i(eta_i²−w). Raw halfnormal moments obey a_0=1, a_1=mu, a_k=(k−1)a_(k−2). Thus
```
EL=w tr(B)+mu² sum_ij B_ij+mu sum_i l_i−T²(1'K1)/2.
v=Var(L)=sum_i E F_i²+w² sum_(i<j)e_ij².
k=E(L−EL)³=sum_i E F_i³
+6 sum_(i<j)e_ij E(F_i eta_i)E(F_j eta_j)
+3w sum_(i<j)e_ij² [E(F_i eta_i²)+E(F_j eta_j²)]
+m3² sum_(i<j)e_ij³+6w³ sum_(i<j<l)e_ij e_il e_jl.
```
Independence removes any lone centered factor. The remaining triples consist of one univariate index, one edge with its endpoints, a repeated edge with one endpoint, three equal edges, or a triangle. With center h and A=L−EL−h, the lower cubic expectation is 1−h+(v+h²)/2+(k−3hv−h³)/6. The integral Taylor remainder is A⁴/6 times integral_0^1(1−u)^3 exp(uA)du. The density 4(1−u)^3 is stochastically below Exp(4) by its survival, so this integral is at least 1/(4+|A|). Cauchy–Schwarz bounds E[A⁴/(4+|A|)] below by (v+h²)²/(4+sqrt(v+h²)). Use directed intervals and the odd degree-15 exponential lower polynomial for exp(EL+h). The result is a certified positive mass lower bound.

For the overlap q=P(R>T,Y_AB,M_C,M_D), one of A,B is negative; take A and double by symmetry. A negative trit forces a zero at pin 1 and another zero at pin 2 or 3. A nonpositive trit forces a zero among the first two. Thus sum upper bounds for the 16 events with the root positive and pins A1,Aa,Bb,C1,Cc,D1,Dd negative, where a,c,d∈{2,3} and b∈{1,2}. These requirements hold at any deeper boundary. For each marginal precision K, set S=diag(+,-,...,-), M=SKS, b=T SK1, C=T² 1'K1. Check b_i>0. For signed deviations y≥0, put r=b'y and b2=b−T e0. Then b2'y≥r and b2' M^-1 b2=C−T²=V>0. Cauchy–Schwarz gives y'My≥r²/V. The simplex density for r is r^7/(7! prod b_i). Integrate the decreasing exp(−r−r²/(2V)) up to five with 192 upper panels. For r≥L=5 substitute z=r²/(2V), use exp(−r)≤exp(−L), and integrate to get V⁴ exp(−L−v)(1+v+v²/2+v³/6)/105 with v=L²/(2V). Bound negative exponentials from above by the reciprocal degree-12 positive exponential series, and include exp(−C/2)/((2pi)^4 sqrt(det G) prod b_i). The gain lower bound is 2p−q. All inverses use verified exact elimination, roots use integer enclosure, and signed products are outward.

## 4. Initial majority

The split set of every line optimizer is a minimum cover of the discordant-pin paths. A discordant edge requires a split at an endpoint, and any cover leaves a consistent assignment to all unsplit vertices. On a 2k+1-vertex path the unique size k cover consists of alternating internal positions: selecting an endpoint still leaves a size-k matching, so exclude the endpoint, force its neighbor, and use induction. On an even path the positive vertices form a minimum cover, allowing zero on all disagreement edges there. The least optimizer keeps those zero values. Isolated positions are unsplit. Almost surely the two rays encounter repeats: the expected repeat indicator is positive and covariance decay makes the variance of its empirical average tend to zero. This proves stabilization for finite ray approximations.

At a raw-positive endpoint let A,B indicate negative nearest pins on its line and C,D on the partner line. If both its first-pass C bits are zero, charge AB or an adjacent positive agreement edge set to zero. In the latter case both ends are split, forcing pattern 0110 in one of two orientations. If its C bits are mixed, exactly one of A,B is one. Matching deletion then requires C or D, since two partner C ones make retaining M strictly better. The inequality (A+B−2AB)(C+D−CD)≤AC+AD+BC+BD−2AB(C+D) gives loss at most 4p_cross−4p_triple, where p_cross=P(R>T,A=C=1) and p_triple=P(R>T,A=B=C=1). The raw-negative pattern 11011 forces the two root C ones by the odd-path rule, giving disjoint recovery. Thus initial majority is at least
```
Phi(−T)−y−2p_0110−4p_cross+4p_triple+p_11011,
y=P(R>T,A=B=1).
```
For y and p_cross condition on R=z. The remaining two coordinates have mean rz, variance v=1−r², correlation rho=(r_neighbor−r²)/v, with r=97/114, r_neighbor=529/722 or 2ar−a². Their joint lower tail is the expectation over S=|N(0,1)| of Phi((T−rz)/sqrt(v(1+rho)/2)−sqrt((1−rho)/(1+rho)) S). This decreases in z,S on z≥T. Use rectangles of width ≤1/64 up to five, lower endpoints, directed upper masses, and a union-tail bound. Normal CDF uses midpoint 1/1024, |phi''|≤1 and monotonic directed brackets between grid points.

For 0110 at zero use the signed even moment expansion in §2.1. Its derivative is the sum over faces i of −s_i phi(0) times the conditional three-sign probability. Compute that probability by the three two-sign moments with covariance C_jk−C_ji C_ik. The derivative is positive. Prékopa log-concavity in the common threshold then gives p_0110(T)≤p_0110(0)<.0121. Use the halfnormal certificate in §3 for the two positive witnesses, with opposite root sign for 11011. The checker obtains:
```
y ≤ 25662444692/10^12, p_cross ≤ 11142921838/10^12,
p_0110 ≤ 12083713264/10^12, p_triple ≥ 1272669402/10^12,
p_11011 ≥ 14557214457/10^12,
initial majority after the matching repair ≥ 568142299591/10^12.
```
The repair calibration uses the pre-motif majority bound above. The strict margins persist on finite neighborhoods by continuity of the Gaussian threshold distributions.

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
