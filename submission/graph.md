# Bounded-local fair cuts approaching the certified energy ceiling

The marked cubic tree rule is equivariant, bounded-radius and fair. Independent labels give expected crossings per vertex below 12233148667589113194/10^20, one count per edge. Proof.md §1 transfers the directional bounds.

Parameters are a=28/57, t=41/57, T=-293/800. Use both replay phases.


## 1. A globally compatible Gaussian source, including branched neighborhoods

Contract the perfect matching to anchors of a loopless 4-regular multigraph Q. Pair its four labeled incidences into two straight continuations, one per physical vertex. The marked universal tree has two lines at each anchor. Set t=41/57, a=28/57, delta=t-a, so t^2+2a^2=1.


For every oriented inward message use predecessor weights t (straight), a,a (perpendicular) on independent standard normal seeds. Each depth has unit variance, and distinct shells and disjoint half-trees are independent. Sum the first J depths over sqrt(J). The two boundary shells give recurrence error sqrt(2/J) in L², so a covariance subsequence yields the exact equivariant law. The anchor score is half its four independent inward modes. Given entering B, independent standard Z,Y give B_f=tB+sqrt(2)aZ, W_0=sqrt(2)aB−tZ, B_l,r=(W_0±Y)/sqrt(2). Future half-trees separate given their entering modes.

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

Starting R coefficients in D,P are a/2,(2a+t)/2, updated by (x,y)→(tx,ty+ax). The line covariance is r(d)=(a+t/2)t^(d−1)+(d−1)a²t^(d−2)/2 at distance d≥1, with the second term zero at 1 and r(0)=1. The first R and next fresh W enter triangularly with nonzero diagonal, proving positive definiteness on each line block. The U coefficient is ell[r(d)−h(d)], h(d)=a t^(d−1), by the same update starting at ell t/2 in P. Opposite lines share only R,U, so their distance-i,j covariance is h(i)r(j)+h(j)r(i)−h(i)h(j).

## 2. First-pass energy

A split vertex covers an edge of the raw disagreement graph if one of that edge's ends is split. The optimum equals the minimum vertex-cover size by the argument in §4. Each component is a path, so its minimum cover equals its maximum matching, obtained greedily from either end. With disagreement indicators d_i and selection g_i=d_i(1−g_(i−1)), expand 19 times and drop the nonpositive remainder. For j successive disagreements of probability A_j(T), stationarity gives
```
c1 ≤ S2(T)+A3(T)−A4(T)+A5(T)
      −(A6−A7)(T)−(A8−A9)(T)−(A10−A11)(T)−(A12−A13)(T)−(A14−A15)(T)−(A16−A17)(T)−(A18−A19)(T),
S2(T)=P(X0>T,X2≤T)=A1−A2.
```
The two consecutive disagreements have XOR equal to endpoint disagreement, and the latter has equally probable orientations. The boundary error is confined to the first disagreement component, which is almost surely finite by §4. The least optimizer has straight trit update f(x,b)=sign(x+2b−1), obtained by minimizing the cubic NAE cost. Two equal pins reset the trit. Incoming trits x,y give C bits [x+f(y,b)>0] and [y+f(x,b)>0]. Take finite ray and shell truncations using these resets.

### 2.1 Directed sign moments

An alternating word and its reversal of signs have centered mass 2^(1−m) times the even sign-moment sum with alternating weights. This gives S2(0), A3(0), A4(0), A5(0) on ≤6 sites. Scale covariances at upper center i by v∈[0,1]. The isolated sign has mean zero. Covariance differentiation and two integrations by parts give each j≠i derivative (2/pi)G_ij M(v)/sqrt(1−v²G_ij²). The remaining conditional covariance given i,j at zero is
    G_kl−G_kj G_lj−[v²/(1−v²G_ij²)](G_ki−G_ij G_kj)(G_li−G_ij G_lj).
Normalize diagonals. Two signs give 2 asin(rho)/pi by angles.

Rational elimination bounds inverse diagonals by 20/3. The variational formula handles restriction, convexity handles interpolation to an isolated center, and Schur inversion plus variance normalization handle conditioning. Whiten the fixed block, then rotate to eigenvalues 1±sv, s²≤17/20. On the mapped 7/4 ellipse, Re v∈[−9/112,121/112], |Im v|≤33/112. Both Re eigenvalues are ≥1/256, one ≥1, and imaginary magnitudes ≤9/32. Ratios |λ|/Re λ ≤73,17/16 give absolute Gaussian integral <9 and normalized moments |λ|²/Re λ. Density differentiation gives magnitude ≤(9s/2)sum(1/|λ|+1/Re λ)<4096. Positive inverse real parts justify analytic domination including nested integrals.

GL32 is exact through 63. Opposite signs certify 32 disjoint Legendre brackets. Mapped weights (1−x²)/(32²P31(x)²) use 2^-72 brackets. Twice the Chebyshev tail of 8192(4/7)^k after degree 63 bounds integration error since both rules have norm one. Add it per level. Use Machin, integer roots, two half-angle reductions and alternating arctangent.

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

For the first eight-site discount sum the two signed orthant floors from the degree-five interpolation of §3.2. The sum is at least 373601877/10^12.

## 3. Tight gain at the active threshold

Take T=-293/800, a=28/57, t=41/57. A raw pin is 1 when its score exceeds T. At one matched root, use the thresholded root R and four consecutive straight rays A,B (opposite directions of one line), C,D (opposite directions of the other). On a ray, truncate to six pins with largest external message +1. A pin b=1 sends 1 if its incoming message is not -1 and 0 otherwise; b=0 sends -1 if its incoming message is not +1 and 0 otherwise. This recursion is monotone in the incoming message. Write M, H, Z for the events that the resulting upper message is respectively -1, nonpositive, and exactly zero. Thus M and H imply the corresponding downward events for actual messages, though Z alone need not. We use Z for partitions of the *truncated* event only.

For two rays of the same endpoint at a positive root, the first-pass line optimization gives local improvement +1 when both messages are nonpositive and at least one is negative (Y), and nonnegative improvement when at least one is negative (N). These identities follow by substituting the 3 x 3 incoming trit pairs into the two adjoining edge bits [x+f(y,1)>0], [y+f(x,1)>0]. Since the two physical endpoints can change their shared matching bit only together, the nonnegative improvement for a pair equals 1[Y_1 N_2]+1[Y_2 N_1] if the root bit is 1. On dividing by two vertices and using line exchange, it suffices to lower-bound P(R>T, Y_AB, M_C or M_D). Starting with the largest boundary +1 makes this a lower bound for every deeper boundary. Other gains are not needed.


Use §1.1 for the marked covariance. To bound p=P(R>T,Y_AB,M_C), split Y_AB into (M_A,H_B) and (Z_A,M_B). Merge six-pin prefixes with equal trits to get 72 disjoint products, 54 under reflection. Stored disjoint eight-pin charts add cases that fail the six-pin test. Upper-bound the overlap of the two reflected events by §3.3. A negative ray needs a first zero and a zero among pins two and three; a nonpositive ray needs a zero among its first two. Thus the overlap enclosure still holds for the extra classes.

### 3.2 Signed Gaussian cylinders

For covariance G, precision K, signs s, and positive scales v, use U=diag(s_i v_i), z=T1 and X=z+UH with independent halfnormals H. Mass is prod(v_i)/(2^n sqrt(det G)) times E exp(L), with L=H'BH+l'H−z'Kz/2, B=(I−U'KU)/2, l=−U'Kz. A one-edge chart can use z=T1+Db and U=DA for sign diagonal D, b≥0 and nonnegative lower triangular A with positive diagonal. Its suborthant uses det A for the Jacobian. Replay checks rational scales, offsets and centers.

Put mu=E H, eta=H−mu, m_j=E eta^j, w=m_2. Raw moments begin at 1,sqrt(2/pi) and recur a_j=(j−1)a_(j−2). Then EL=w tr(B)+mu² sum B+mu sum l−z'Kz/2. Centered L is the sum of J_i=(2mu sum_j B_ij+l_i)eta_i+B_ii(eta_i²−w) and e_ij eta_i eta_j, e_ij=2B_ij. Variance sums EJ_i² and w²e_ij². The third moment sums EJ_i³, 6e_ij E(J_i eta_i)E(J_j eta_j), 3w e_ij²[E(J_i eta_i²)+E(J_j eta_j²)], m_3² e_ij³, and 6w³ e_ij e_ik e_jk over the corresponding indices. Every omitted term has a lone centered factor.

At distinct nodes x_i and Lagrange bases l_i, set h_i=[1+(1−2l_i'(x_i))(x−x_i)]l_i². The Hermite floor is sum exp(x_i)h_i, by remainder exp(ξ)prod(x−x_i)²/(2k)!≥0. Two nodes use the centered variance and third moment. Jensen gives floor one. Multiply by a lower exp(EL) and choose exp(x_i) endpoints by signs of lower E h_i.

Five six-site and two eight-site orientations use three or four nodes. Lower L coefficients on the nonnegative orthant, expand powers and group moments by odd-coordinate count using directed sqrt(2/pi). Four-bit exponents cover degree ≤14. For |x|≤1/2 use degree-16 Horner with error ≤36/(35·2^17·17!), then outward square after range halving.

### 3.3 A quantitative overlap upper

For q take A negative and B nonpositive, then reflect. Cover the positive root and negative A1,Aa,Bb,C1,Cc,D1,Dd, a,c,d∈{2,3}, b∈{1,2}. For signed y≥0 set C=T²1'K1, b=T diag(sign)K1>0, v=b−Te_0, W=C−T². Put r=v'y, p_i=v_i y_i/r. From v'M^-1v=W,
    y'My=r²/W+2r²B(p), B=p'[diag(v)^-1 M diag(v)^-1−11'/W]p/2≥0.
The Jacobian is r^7/(7! prod v_i) with uniform simplex law. Integrate exp(−C/2−r−r²/(2W)+γrp_0−r²B), γ=−T/v_0, divided by (2pi)^4 sqrt(det G).

For I_n=∫r^(n−1)exp(−r−r²/(2W))dr on [0,∞),
    I_1=exp(W/2)sqrt(W)[sqrt(pi/2)−sqrt(W)sum_(k≥0)(−W/2)^k/(k!(2k+1))].
Use alternating 17/18. Parts give I_2=W(1−I_1), I_(n+1)=W[(n−1)I_(n−1)−I_n]. Round outwards.

Group coefficients of B^k times prod e_i! by h=e_0 into A_h. The simplex monomial mean is 7! prod e_i!/(7+|e|)!. The normalized radial tilt terms are γ^j I_(8+2k+j) sum_h A_h binom(j+h,h)/(7+2k+j)!. Beyond j=6 use geometric ratios ≤γ sqrt(W/(8+j)) at k=0 or γ sqrt((8+2k+j)W)/(j+1) otherwise, by p_0≤1 and I_(n+1)^2≤I_n I_(n+2)≤nW I_n². For b=2 use upper k=0, lower k=1, upper k=2 as u−ta+t²b/2, t≤min(1,a/b), from exp(−r²B)≤exp(−tr²B)≤1−tr²B+(tr²B)^2/2. Reflect C,D.

Subtract 15 parent intersections and seven disjoint excluded cylinders. The latter have 01011 on one of A,C,D, 00 on the other two and B=10, or B=1011 with A,C,D all 00 except at most one 0100. Those two long words send 0 and 1 for every exterior trit, violating the event. Subtract directed lower masses and reflect A,B.

## 4. Initial majority


Split sets are minimum covers of discordant paths (§2). The unique minimum cover of an odd path alternates inside. An endpoint leaves a same-size matching. Exclude it and induct. On even paths the positive cover permits all discordant C bits zero, as in the least optimizer. Isolated vertices are unsplit. Positive repeat mean and decaying covariance imply a repeat on either ray almost surely.

At a positive root let A,B and C,D indicate nearest negative pins on its line and partner line. Two zero root C bits require AB or a zeroed positive agreement edge. The latter splits both ends, forcing 0110 in one of two orientations. Mixed C bits require exactly one of A,B, and deletion of M requires C or D. Use (A+B−2AB)(C+D−CD)≤AC+AD+BC+BD−2AB(C+D). Pattern 11011 at a negative root forces two positive C bits by the odd-path rule. Thus majority is at least
```
Phi(−T)−y−2p_0110−4p_cross+4p_triple+p_11011,
y=P(R>T,A=B=1), p_cross=P(R>T,A=C=1), p_triple=P(R>T,A=B=C=1).
```

Given R=z, the pair has mean rz, variance v=1−r², rho=(r_neighbor−r²)/v, with r=97/114 and r_neighbor=529/722 or 2ar−a². Its lower tail averages Phi((T−rz)/sqrt(v(1+rho)/2)−sqrt((1−rho)/(1+rho))S) over S=|N(0,1)|. This decreases in z,S for z≥T. Use rectangles of width ≤1/64 up to five, lower endpoints, directed upper masses and union tails. CDF midpoint mesh 1/1024 uses |phi''|≤1 and directed monotone brackets.

For 0110 at zero use the signed even moment expansion in §2.1. Its derivative is the sum over faces i of −s_i phi(0) times the conditional three-sign probability. Compute that probability by the three two-sign moments with covariance C_jk−C_ji C_ik. The derivative is positive. Prékopa log-concavity in the common threshold then gives p_0110(T)≤p_0110(0)<.0121. Use the halfnormal certificate in §3 for the two positive witnesses, with opposite root sign for 11011. The checker obtains:
```
y ≤ 25662444692/10^12, p_cross ≤ 11142921838/10^12,
p_0110 ≤ 12083713264/10^12, p_triple ≥ 1272902052/10^12,
p_11011 ≥ 14572880678/10^12,
initial majority after the matching repair ≥ 568158896410/10^12.
```
The additional motif debit leaves at least 542496451718/10^12. The repair calibration uses the pre-motif margin. The strict margins persist on finite neighborhoods by continuity of the Gaussian threshold distributions.

## 5. Simultaneous local erasure cannot increase the same edge-split energy

Let S be the positive edge set and E(S)=sum_v NAE. For a patch B choose the least minimizer among subsets of S agreeing outside B, which exists by submodular meet/join closure. Write C_B for its removal set. For any union D of earlier removals,
```
E(S − (D union C_B)) ≤ E(S−D)+E(S−C_B)−E(S−(D intersection C_B)).
```
The last set is feasible for patch B, so the last two costs sum to a nonpositive quantity. Thus independently activated patches, each minimized from the original S, can be intersected without raising cost. Transfer through a large-girth marked graph by extending the partial right translations on reduced words of length ≤L to permutations. Evaluation at the empty word excludes any nontrivial relation of length ≤L. Split each anchor along its marked M edge.

### Relaxed terminal majority

Free the internal edges of an anchor ball. Extend each boundary C ray through three straight physical vertices, freeing the joining C edges. Extensions are disjoint on the tree. Drop interior C ceilings, keep raw-negative interior M pinned zero, fix extension M edges to their raw bits and far C edges to one. The fixed bits dominate the actual ones. The meet of actual and relaxed minimizers is feasible for the actual patch and their join for the relaxed patch. Outside terms are constant in each problem. Submodularity and leastness then bound the actual root bits above by the relaxed least minimizer.

Clamp the parent-edge cost difference to {-1,0,1}, since a child bit changes adjacent NAE cost by at most one. Four choices yield sign(x+y) at a cubic vertex. Joining M gives sign(x−1) when pinned or sign(x+sign(y+z)) otherwise for straight x and perpendicular y,z. Extend by sign(x+2b−1) starting at +1. The k-extension law L_k decreases, and its full update is ≤L_(k+1). Thus d full updates from L_3 need radius ≤d+3.

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

Start with three straight sites. With matching bit zero, a nonpositive child gives a negative output. With bit one, a nonnegative child gives a positive output. Meet with the preceding law. At $d=306$ the main, auxiliary and contrast axes have (full radius, endpoint) $(159,367),(171,389),(196,484)$ and stride five outside the full region. The integration limits are 789 and 966. Run ten increasing-contrast sweeps, retaining stale bounds on skipped rows. The code packs a fixed out-of-grid suffix separately. Skip auxiliary writes with no read before replacement. Likewise, the continuation omits auxiliary rows outside the main stencils. Each sweep advances at most 509 levels.

Join each consecutive pair of tables u,v in lower CDF coordinates (m,1−p), starting at the table after the sweeps. The auxiliary output has three lower quadratic Bernstein coefficients from both same-input products and their exact polarization, with downward rounding. The main output is bilinear in the straight and auxiliary laws, including the split pin adjustment. Its cubic coefficients are M(u,A00), [2M(u,A01)+M(v,A00)]/3, [M(u,A11)+2M(v,A01)]/3, M(v,A11). Subtract the chord and check positivity when the endpoint coordinate increases. Otherwise use the verified start and the monotone decrease of the true physical laws. Midpoint subdivision to depth 12 gives a residual at least 1/(3 PATH_DEN 8^12) per moving coordinate, with PATH_DEN=121·2^74. Because each positive slope is at most 16384, the path advances by at least one part in 3 PATH_DEN 8^12·16384 per physical level. 2 such lengths plus 10·509+3 bound the patch radius.

Decode units-of-two knots by mixed differences and separable rational cubic interpolation. Apply signed masked corrections to the first table and unsigned subtractions to the second. Check both simplexes and digests.

For the root, use independent unit Gaussians R,U and symmetry to restrict U to [0,4]. Use R in [−4,4]. Each line uses its largest coordinates in each rectangle and one auxiliary pass to bound the ordered categories --,-0 | -+,00 | 0+ | ++. If the ordered upper tail triples are (A,B,C) and (D,E,F), the positive endpoint expectation is A+D at a pinned root and 2(CE+BF−BE)+A(1−F)+(1−C)D otherwise. Check the payoff against the literal minima and nine simplex vertices. Assign omitted mass the maximal payoff. Round auxiliary row indices upward to multiples of six, capped at the last row. This increases the P coordinate at fixed D. Directed root steps 1250/40000, doubled beyond magnitude 2, give 8158/16384 < 1/2.

### 5.2 Exact fair marginal and cut density

Full activation gives majority below half. Use proof.md's calibration from the pre-motif law of §4. Strict margins persist under finite windows, and transitivity gives fair marginals.

If adjacent majorities differ, their shared edge is a minority at one endpoint. A split endpoint has only one minority edge, so expected crossings per vertex are at most expected splits, strictly below c.
