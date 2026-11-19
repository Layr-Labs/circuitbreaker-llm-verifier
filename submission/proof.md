# A marked Gaussian ordering with asymmetric signed margins

For every sufficiently large simple cubic H on n vertices

    cw(H) <= (gamma-epsilon0)n,
    gamma=500000000/3827915978,

for a fixed epsilon0>0. The separate `family.md` transfers this to one
fixed total uniform NP language with full-B2 circuit lower bound
`(2413957989/500000000) N=4.827915978 N` at every sufficiently large
actual integer length N, with zero subtracted remainder.

Set

```
eta=9/10000, theta=eta/2, T=49/2000,
U=1633094581/12500000000=.13064756648,
A=1900519160504672063862531/25000000000000000000000000000,
Csave=2818956503/10^14.
```

The manifest binds the replay; `nine-choices.txt` only selects boxes.
Only the ordinary double-short star is spent.

Use the complete positive radial/fibre source in
Puthranps/cubic-cutwidth-retiming@`6b299fde41307ff47eb4a699e140fd009e61a51d`,
cutwidth Sections4,7.2--7.5. Uniform errors include atoms and repeated
images and are fixed before graph size.

## 1. Crossing matchings and Gaussian comparison

For vertex bits `b` on a cubic graph, write `F(b)` for the minimum cover
of edges joining opposite bits. It equals their maximum matching size
by bipartiteness. For edge bits `x` put

```
Q(b,x)=sum_v NAE(b_v,(x_e : e incident to v)),
F(b)=min_x Q(b,x).
```

NAE indicates both bits occur. Nonzero-cost vertices cover all crossing
edges; conversely a crossing cover lets outside vertices consistently
fix their incident edge bits. NAE has nonpositive mixed differences.
Taking meet/join of joint minimizers proves F submodular. F is
complement-invariant and changing one bit costs at most one.

For centered unit Gaussian fields `X,Y` with covariance matrices `A,C`,
put `Z_s=sqrt(1-s)X+sqrt(s)Y` with the fields independent. At threshold
zero, Gaussian integration by parts gives the Price identity

```
(d/ds) E F(1[Z_s<0]) = sum_(i<j) (C_ij-A_ij) phi_(ij,s)(0,0)
                     * E[Delta_ij F | Z_(s,i)=Z_(s,j)=0].            (1)
```

Here phi is the bivariate density and Delta the four-corner difference
with other coordinates thresholded. Smooth the multilinear extension's
thresholds: the Gaussian Hessian identity pairs the off-diagonal terms;
the two boundary derivatives give (1), and diagonal coefficients vanish.
Positive-definite approximations extend by marginal continuity and
bounded F. The proof applies at any common threshold. Entrywise covariance
domination therefore decreases the expectation.

For a single field of covariance `C` let `f(t)=E F(1[X<t])` and
`W_j=X_j-C_ij X_i` under conditioning at `i`. Set `k_j=1-C_ij` for
that `i` and put `g_i(t)=E[D_i F(1[X_-i<t]) | X_i=t]`. The density
of the conditioned boundary gives `f'(t)=phi(t) sum_i g_i(t)`.
Symmetry gives `g_i(-t)=-g_i(t)`. On interpolating from `-t` to `t`
the conditioned tests are `W_j<k_j s`, so differentiating their
boundary integrals gives

```
2g_i(t)=integral_(-t)^t sum_(j!=i) k_j f_(W_j)(k_j s)
           E[Delta_ij F | X_i=X_j=s, other bits tested at s] ds.     (2)
```

Condition first on X_i and then on W_j=k_j*s; f_(W_j) is its density.
Thus f is even and decreases for t>=0. Independent coordinate noise
and marginal limits include singular fields.

## 2. From nested cuts to a layout

A vertex ordering with every prefix crossing cover <=M yields width
<=M+6 by the following edge interpolation.
[Cutwidth proof, Lemmas 2.2 and 3.1–3.2](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/737818f9911e00e9d4909c644c771e0724fb6c37/submission/cutwidth.md#L279).

Clamp same-anchor edge bits. Submodularity gives monotone least
minimizers along prefixes. Flipping one anchor costs at most M+1.
In the monotone interpolation let A gain and B lose nonzero list cost,
P retain it. Resetting any subset of A in the destination optimum
gives Hall's inequality: its changing B-neighborhood is at least its
size. Process P–P,B–P edges and three-edge A stars first. Remaining
A blocks have at most two B neighbors, hence form edges or half-edges
on B. Hall implies tree components with at most one half-edge, or
unicyclic components with attached trees. Peel leaves, then finish the
half-edge or cycle. At most one block is unfinished between blocks,
two within a block, so split lists cost <=M+3. The edge-position bag
adds at most two endpoints. In median-incident-edge vertex order,
a crossing edge before the median charges the outside vertex's lowest
edge, one after it the inside vertex's highest; both enter the bag.
At most one edge at the median remains. Thus width <=M+6.

## 3. Radial keys and uniform finite transfer

Every finite simple cubic graph has finite covers of arbitrarily large
girth: assign free generators outside a spanning tree; separate the
finitely many short nonempty reduced words by partial permutations on
0,...,length, extend to permutations, and take their product regular
action. For a degree-d cover, repeated bits have objective dF: lift a
base cover for <=, and average an upstairs cover to a fractional base
cover, summing its constraints on a maximum matching, for >=.

On the universal cubic tree use the unit positive row
f_v(z)=sqrt(c_alpha)*alpha^dist(v,z),
c_alpha=(1-2alpha^2)/(1+alpha^2), .7<=alpha<1/sqrt2.
Sphere counts give
r_j(alpha)=alpha^j[(j+1)-(j-1)alpha^2]/(1+alpha^2), tending to
w_j=(1+j/3)2^(-j/2). Fibre norms give physical unit rows b_v and
covariance domination by Cauchy--Schwarz. Mark with complete-source
D_uv=<b_u,b_v>-r_1(alpha)>=0 before truncation/noise.

Truncate at radius R and renormalize before fibre norms; add independent
physical tie noise. Geometric squared tails and contraction under fibre
norms give uniform covariance convergence. Upstairs repeat physical
noise along fibres. Choose alpha, then R/noise, then cover girth.
Corrected keys inherit original keys within radius3, at most22 times.
For fixed R, bounded differences gives variance O_R(n) for each fixed
threshold and key-interval count. Chebyshev on a fixed mesh, interval
expectation<=n*width/sqrt(2pi), and tail counts give a simultaneous
all-threshold bound. Ties inherited from a single original key cost<=22;
Section2 adds6. Choose these errors before graph size. No cover degree,
physical alias or graph-dependent parameter enters the limiting bound.

## 4. Certified reference expectation

For odd J, half the directed crossing nonbacktracking j-walk counts
N_j give cover upper bound sum_(j<=J)(-1)^(j+1)N_j on paths/cycles.
Deleting each degree3 crossing vertex changes this by O_J(1).
At the critical Gaussian wave, sum_(u~v)X_u=2sqrt2 X_v, so its three-edge
crossing star is impossible; its probability tends uniformly to0.
Thus the reference expectation per vertex is at most
Gamma_J=sum_(j<=J)(-1)^(j+1)3*2^(j-2)p_j+o(1), with p_j the alternating
j-path sign probability.

The checker proves Gamma_25<U using60-point Gauss--Legendre tensors on
[0,3],58-digit directed decimals, root brackets of radius1e-45 and
fixed scale2^104. The path AR(2)
recurrence is (sqrt2,-1/2), innovation variance1/12. The positive pair
density is3/(2pi)*exp[-9(x^2+y^2)/2-6sqrt2*xy], and extension density
sqrt(6/pi)*exp[-6(z+sqrt2*y+x/2)^2]. Anchor reversal permits each tail
coordinate to be charged with one/two-step completion bounds1/2,1/10.
The resulting tail at coordinate i, anchor h=min(i,m-2), is
c2*exp(-9L^2/2)/(54sqrt2*L^2)*H(h)H(m-2-h), L=3,
H(r)=10^(-floor(r/2))*2^(-(r mod2)); sum over i and both orthants.
The Gaussian Fourier moment bounds coordinate2N derivatives by
f_m(0)*q_i^N*(2N)!/(2^N*N!). Multiply by the Gauss remainder
L^(2N+1)(N!)^4/[(2N+1)((2N)!)^3], telescope coordinates, and use unchanged
interior exp(-39x^2/2) integral/quadrature masses<=1/4 for long paths.

Legendre roots have disjoint checked sign brackets; the derivative
identity and |P_N''|<=N^4 enclose weights. The positive transition
A_ij B_ik C_jk has A,B<=1,sum_k C_jk<=1. At scale S its per-step
rounding loss is at most(D_A+D_B+N*D_C+2N^2)/S. Contraction propagates
this explicit loss; the faster integer shifts are exact at S=2^104.
Apply lower or upper endpoints according to each alternating sign.
The mandatory replay checks every bound and returns Gamma_25<U.
This is the published path-certificate argument with independently
checked roots, tails, tensor error and rounding at the stated scale.

## 6. Safe advances and strict matching savings

Apply Gaussian comparison to the original field, then account for
advances at each threshold. Let B be its below-threshold set. All
named vertices below are physically distinct and listed stars complete;
lightness proves this faithfulness when counting events.

For the seven-key advance, name the root `o` with neighbors `a,b,u`.
Use complete stars at `a`, with further neighbors `a1,a2`, at `u`,
with further neighbors `w,w'`, and at `w`, with further neighbors
`w1,w2`. Put `a,b,a1,a2,w,w1,w2` below the threshold and `o,u,w'`
above. Delete `o` from the crossing graph. Then `a` is isolated, so
restoring the original `o` adds one to the maximum matching. After
advancing `o` into `B`, the only new crossing edge is `ou`. It can
be replaced by `uw` in any matching using it, since `w` has no
other crossing neighbor. The maximum matching has decreased by one.

For the nine-key advance take stars on the tree with edges

```
oa, ob, ou, aa1, aa2, bw, bw', uz, uz', ww1, ww2, zz1, zz2.
```

Put `a,a1,a2,w,w1,w2,z,z1,z2` below threshold and `o,b,u,w',z'`
above. After deleting `o`, the external crossing graph has the
isolated vertex `a` and closed pairs `b–w` and `u–z`. Originally
`oa` adds one to a maximum matching. After the advance the new
edges `ob,ou` can be covered by choosing `b,u` in minimum covers
of their closed pairs. Thus this advance also saves one.

The positive support alone is strict: starting from the guarded
configuration, turn other outside named bits inside. Submodularity
keeps the root marginal <=-1, its minimum possible value.
Guards only distinguish counted events. A seven-key guarded
success has two inside neighbors of its root, and determines `u,w`.
At most two seven names succeed, ordinarily weighted1/2 each. Nine
names are disjoint: the inside root neighbor and the inside child in
each other branch are unique. Their one-inside-neighbor roots are
disjoint from seven successes with two inside neighbors.

The weak-pair weight is zero.
At each root lower its key to the minimum of its original key and the
support maxima of eligible seven/nine names, with the graph/mark masks
fixed across thresholds. A moved root has new_key<t<=old_key, so its
certificate applies. Process counted strict roots first, then any other
admitted roots. Submodularity preserves all savings. Every inherited
key is an original key within distance three, with multiplicity at most22.

The seven-key count is strengthened without changing this event. Give
an eligible name weight1 when its uninspected inside neighbor b is
high-marked, and weight1/2 otherwise. If b is high-marked, the competing
private choice would inspect b and is ineligible. Thus total successful
weight at a root is still at most1. Its conditional probability floor
applies to each naming because b's uninspected edges enter none of its
covariance paths.

## 7. Source excess, defects and paired inequalities

For squared physical rows P_v(p)=b_v(p)^2, q=alpha^2,beta=q/(1+2q^2),
the exact nonbacktracking generating function gives
P=(1-3beta)(I-beta Adj)^(-1). Hence
sum_(u~v)P_u(p)=beta^(-1)P_v(p)-(beta^(-1)-3)1[p=v].
Adjacent fibre ratios h lie in[q,q^-1], and the symmetric nonnegative
summand is P_v(p)phi_q(h), phi_q(h)=sqrt(h)-sqrt(q)(1+h)/(1+q).
Its sum is D_vu. This includes every root atom and repeated image.

Telescope the fibre unit directions along a j-path. Row norm ratios
lie in[alpha^j,alpha^-j], giving
w_j-o(1)<=Cov(Y_v,Y_u)<=w_j+j^2*sqrt2^(j-1)*H+o(1)
if all path edges have excess<=H. The sharper per-path source bounds
in family Section4 are also used. At eta=.0009 these uppers are<1
through the diameter6 of each inspected stencil, proving all names
physically distinct, including relaxed terminal edges.

Let lambda0=2sqrt2,c=3-lambda0,R_v=sum_(u~v)Y_u-lambda0 Y_v.
Fibre triangle and the resolvent give a row differing by o(1) from
0<=r<=c*b_v. Put S=<r,b_v>=sum_(e~v)D_e+o(1). Then
||r||^2<=cS, and at distance j its covariance with a key lies between
sqrt2^(-j)S and sqrt2^j S, up to uniform errors. Conditioning a pivotal
interpolation endpoint at zero subtracts at least t*2^(-j)S^2.
Maximizing over S<=c gives c^2*h_j(t), where h0=1-t until1/2 then
1/(4t), and hj=1-t/2^j for j>=1. For a clean star with S<=3H the
bound is3cH-9tH^2/2^j; stratified first outward S<=.006+2*.007=.02.
Further conditioning/positive quadratic tilt decreases covariance.
At depth d the midpoint error has impulse coefficients
(i+1)2^(-i/2)/2,0<=i<d; use the triangle inequality for its SD.
All scaled deviations sqrt(t/(1-t))*SD are nondecreasing, justifying
right-time strips. No residual independence is asserted.

The paired off-edge law is sum_off D>=sigma*D_e-o(1),sigma=2021/1250.
For one fibre let x be mass in the branch toward the other endpoint,
h=alpha^2+(alpha^-2-alpha^2)x,x'=alpha^2(1-x)/h. The excess is
P_v f(x), f(x)=sqrt(h)-alpha-(alpha^-1-alpha)x. Concavity gives two
away terms>=P_v f(1-x) minus the endpoint root-atom Lipschitz debit;
include both endpoints' atoms. With y=sqrt(h)+1/sqrt(h),z=alpha^2+alpha^-2,
the paired numerator is R(y)-zy/(alpha+alpha^-1),
R(y)=sqrt(z^2-(sqrt((z+1)^2-z*y^2)-1)^2).
Its concavity follows from ((u-1)^2)''=2z((z+1)^2-u^3)/u^3>0.
The endpoint chord yields coefficient
[2(alpha+alpha^-1)sqrt(z-1)-2z]/(alpha+alpha^-1-2), tending to
9sqrt6+12sqrt3-15sqrt2-20>sigma. Root debits vanish uniformly.

Finally, at a vertex the weighted triangle inequality is
`phi(t)<=r phi(r)+s phi(s)` for `t,r,s` in `[1/2,2]` of sum three,
where `phi=phi_(1/2)`. Here is a scalar verification. Fix `t`; an
interior stationary point of the right side has `r=s` or
`sqrt(r)+sqrt(s)=9/(4sqrt(2))`. At the boundary set `r=1/2`,
`s=5/2-t`. In the rationalized `phi`, the numerators at `s,t`
agree. If `D(x)=sqrt(x)+sqrt(2)(1+x)/3`, then
`sD(t)-D(s)=sqrt(s)(sqrt(st)-1)+sqrt(2)(st-1)/3>=0`.
At `r=s` a minimum requires `r<=81/128` by the second derivative
of `x phi(x)`. Put `r=x^2/2`, `u=x-1` so `0<=u<=1/8`.
The difference after squaring the positives
`8+3x^3-4x^2-x^4` and `3sqrt(6-2x^2)` is
`u^2(15-6u-5u^2+8u^3+3u^4+2u^5+u^6)>=0`.
For the other stationary points put `L=r+s`, in `[81/64,41/32]`.
The difference multiplied by `sqrt(2)` is `Q_0(L)-sqrt(6-2L)`, with
`Q_0(L)=-729/256+27L/8+(2/3)(4-2L-L^2+(81/32-L)^2/2)`.
Here `Q_0'(L)=17/48-2L/3>=-1/2`, so the difference is increasing
and its value at `81/64` is the verified balanced case. Continuity
covers the boundary zeros.

By (3), the expectation with weights `P_v(p)` of the absolute
error in the sum of the three fibre ratios from three is at most
`2(beta^(-1)-3)`. Clamp the ratios to `[1/2,2]` and adjust them in
that box to sum three. Their weighted total displacement is `o(1)`.
The bounded derivatives of `phi_q` and `h phi_q(h)` on the common
compact interval show that `D_i<=sum_(j!=i) sum_p h_j d_j(p)+o(1)`.
In particular `D_i<=2 sum_(j!=i)D_j+o(1)`. When one low edge is
charged from both ends, the two weights of its symmetric summand
are `h` and `1/h` and their sum is at most `q+1/q`. This improves
the budget available from edges between unmarked vertices.

## 9. Outer derivative witnesses and Price transport

For endpoints at distance d>=2, alternate the middle colors and match
each middle third-neighbor color. At an endpoint require both outside
roots opposite its adjacent middle, or delete a same-colored closed
EVEN crossing arm with its boundary fixed. Its matching size is a
constant in all four endpoint flips. Let A,B be the separate exterior
matching increases and C the joint increase. For odd d the attachments
are in opposite bipartition classes; symmetric differences allow both
to be exposed when A=B=1, so C>=A+B. For even d they are in the same
class: C=2 forces A=B=1, while A=B=0 implies C=0 by the augmenting-path
criterion; hence C<=A+B. The alternating path contributes
floor((d-1+x+y)/2) in relative flips. Combining parities gives physical
mixed difference<=-1 for every exterior.

Under critical tree conditioning X0=Xd=s, write the centered middle
block J. Its score shifts are J_i+b_i*s,
b_i=(w_i+w_(d-i))/(1+w_d)-1. Third-neighbor shifts are
lambda0*J_i-J_(i-1)-J_(i+1)
 +(lambda0*b_i-b_(i-1)-b_(i+1)-(3-lambda0))*s.
Endpoint sibling midpoints are -(J_1+(b_1+3-lambda0)*s)/2 and its
reverse. Their contrasts are independent of J and each other by subtree
swaps. The AR(2) path law gives middle precision diagonals36 at the
ends,39 inside, first off-diagonal -18sqrt2, second6, determinant
9*12^(d-1)(1-w_d^2). For d=2 the exponent is33*x^2/2.
Both global signs are included. Alternating magnitudes make all density
cross terms positive. The checker contracts the JOINT middle density.

Even arms of sizes2,4,6 alternate along their closed component; each
complete star fixes all boundary signs. Independent ideal branch
contrasts give continuation density at b from
b+sqrt2*a+p/2 plus the directed threshold offset. Root/closing strips
use the exact midpoint shifts above. Freeze their sufficient corners,
then average the entire density exponent by Jensen with every shared
coordinate variance included. Arm lengths are distinguished by their
first closing sign; endpoint root signs distinguish baseline and arms.
Thus all retained arm/sign classes are disjoint before addition.
At distance2 its one-arm reflection is counted once. No-arm D4 keeps
its separately checked precision/stencil in family Section13.
The adjacent witness uses its two disjoint oriented private sign cases.
Every moving-mean offset is checked on its declared inner interval,
including the tighter adjacent .01 response. All omitted derivatives
are nonpositive.

For pair correlation w let r=(1-w)/(1+w). The density in(2) is
sqrt(r)/(2pi)*exp(-t^2/2-r*s^2/2). With w<=a rational upper evaluation w_d^+ of w_d plus1e-9, use
r_d=(1-w_d^+)/(1+w_d^+): sqrt(r)exp(-r*s^2/2) increases in r for
|s|<1. A response valid through b0 pays to T via d=min(T,b0) and

```
J(T,d)=d^2/2-(1/8+r_d/24)d^4+r_d*d^6/72
 +(T-d-(T^3-d^3)/6)*(d-r_d*d^3/6).
```

This integrates (1-t^2/2)(1-r_d*s^2/2) on0<=s<=min(t,d),0<=t<=T.
Multiply by directed multiplicity3*2^(distance-1), its checked response
floor and sqrt(r_d)/(2pi). Add a shorter-range improvement only as its
floor DIFFERENCE on that range. The exact total is
O(T)=7096191013529/(25*10^16)>Csave. Evenness/monotonicity extend it
to |t|>=T. Compare physical and reference fields at the SAME threshold,
then apply safe advances. This does not spend the central Price budget.

### Price transport

Set Z=sqrt(1-s)X+sqrt(s)Y, r=s/(1-s), and condition both pivots zero.
Ideal whitened middle/sibling drivers g remain independent standard,
independent of X_end,Y. Apply the same linear maps to Y to define H.
Drivers are g+sqrt(r)H(Y), with
Section 7's impulse errors. Conditional Y is centered. Fix sizes/labels
and both reflections; pairing g,-g uses cosh>=1, leaving
exp(-r||H||^2/2), of mass det(I+rV)^(-1/2), V=Cov(H|Z_end=0).
The centered tilt decreases error variances; Gaussian correlation
pays strip products. Only ideal coordinates factor; Y may repeat.
At zero the exact arm root/continuation densities are 2fD(x+2a) and
2phi_(1/sqrt(12))(b+sqrt(2)a+p/2), fD=density of |N(0,1/3)|.
Use Section 9's zero-offset closing margins and Section 7 impulses.
Pivot density is >=1/[2pi sqrt(1-w_d^2)].

For d>=3 let Q be full path precision, CE the endpoint covariance,
and S the conditional middle covariance with regression L. For
Bx=x_mid-L*x_end, M=B^T S^-1 B=Q-embed(CE^-1) is PSD and
tr(M Cov(X))=d-1.
Its adjacent entries are -12sqrt(2) at ends,
-18sqrt(2) inside; second neighbours have entry6; endpoint subtraction
adds w_d/(1-w_d^2). The source bound is
D2<=C*(D_left+D_right), C*=(sqrt(10)+3)/(2sqrt(2)). For fibre ratios
u,v>=1/2, u+v<=5/2 put A=sqrt(2u),B=sqrt(2v),S=A+B,P=AB,
z=(sqrt(10)+3)/2, psi=sqrt(uv)-(u+v+.5)/3. Then
6[C*(phi(u)+phi(v))-psi] equals
z(3S-4-S^2)+(2z-5)P+S^2+1. Min P gives
(z-1)(S-2)(3-S) on [2,3], or 3(S-3)[z-(S+3)/2] on [3,sqrt(10)],
both >=0. Use D2<=1/6 with alternating cap weights
0,theta,0,theta,..., theta=38-12sqrt(10). Adjacent terms cancel, giving
the actual trace <=(d-1)+2floor((d-1)/2)theta+2w_d/(1+w_d).
Schur/AM-GM/Hadamard bound the inverse determinant by
(1-s)^(1+(n+m)/2)/(1+mu_d*r)^((d-1)/2) for arm sizes n,m.
Use mu3,...,mu7=1.467,1.281,1.213,1.152,1.129.
For d2 H0=sqrt(33)[Y_w-(4sqrt(2)/11)(Y_u+Y_v)].
Var(H0)=1+(192/11)D2-24sqrt(2)(D_left+D_right)
<=571/11-16sqrt(10)<1.313. Clean means:

```
1+48H2(4sqrt(5)-5sqrt(2))/11 <1.04905,
1+24H2(C*-sqrt(2)) <1.111,
1+H2(24C*-2008sqrt(2)/95) <1.135,
1+H3(24C*-123sqrt(2)/7) <1.193,
1+(2/5)[4H3(12C*-18sqrt(2))+3/11] <1.117.
```

Joint root cases use Var(H)+3Var(R)/(2c^2)<=3/2, c=3-2sqrt(2),
H=sqrt(3)(Y_a-Y_b) at an endpoint.
Indeed for branch square roots A,B,W>=1/sqrt(2), sum squares3,
f=A+B+W-2sqrt(2) in [0,c], one has 2(A-B)^2+f^2/c^2<=1.
For f<=2c/5 use
(A-B)^2<=1/2-sqrt(2)f; otherwise use
(A-B)^2<=6-2(2sqrt(2)+f)^2/3 and convexity, with endpoint values
1 and (160-96sqrt(2))/25<1. Integrate. Choose h=Var(H) before tilt:
h<=p gives gain g=[1-(1-p)s]^(-1/2); otherwise use root variance
c^2 min(h0(s),1-2p/3). Baseline2 uses p=3/4 and min(gA,B)^2,
A,B the low/high side convolutions. Arms2 uses p=2/3 and the minimum
of LL,LH,HH EXTRA polynomials. D3 uses p=2/3 and minimum of four
complete middle integrals. Never add cases.

### D7–D10 and interfaces

Distance adjacency on a high-girth cubic cover obeys
`P1=t,P2=t^2-3,Pd=t*P_(d-1)-2P_(d-2)`. Put L=2sqrt(2).
For P7=t^7-13t^5+48t^3-44t its tangent residual at L is
(t-L)^2R(t), R=t^5+2Lt^4+11t^3+6Lt^2+56t+64L.
Here R'=5(t^2+4Lt/5)^2+(37/5)(t+30L/37)^2+632/37>0
and R(-3)=280L-708>0.
For z=t^2-8 the next polynomials are

```
P8=88+184z+94z^2+17z^3+z^4,
P9=t(48+120z+72z^2+15z^3+z^4),
P10=208+640z+508z^2+158z^3+21z^4+z^5.
```

Their tangent residuals are respectively

```
184(t-L)^2+z^2[(z+17/2)^2+87/4],
(t-L)^2[120(t+2L)+t(t+L)^2(z^2+15z+72)],
640(t-L)^2+z^2[76+t^2((t^2-3/2)^2+47/4)].
```

The odd expression is nonnegative for t>=0. For -3<=t<=0,
expanding P9/t at t^2=9/2 gives coefficients
`271/16,-17/4,-12,1,1`. Thus `|P9|<=9363/4<1920L`,
whereas its tangent `1968t-1920L<=-1920L`.
Let Sd be half the sum of excesses over directed distance-d pairs;
S1=sum_e D_e. All four residuals are nonnegative on [-3,3]. Pairing their
positive-semidefinite matrices with the physical row Gram matrix,
using unit diagonal and Pd(L)=3*2^(d-1)w_d, gives
`Sd>=multiplier*S1-o(n)` with multipliers 532,1040,1968,3620.
The last three are bounded by derivatives 368L,1968,1280L.
Finite-reference replacement has uniform o(n) error.

Use the alternating-path witness above, two endpoint contrasts and
optional size-two arms. All sign/arm classes and all distances are
disjoint. The common banded precision and trace formula give:

```
d   precision determinant   both-sign normalizer   pair density   mu
7   24541056                39.94                  .1665          1.129
8   305550144               56.2                   .1635          1.106
9   3748902912              78.3                   .1616          1.091
10  45586457856             109.2                  .1606          1.074
```

Endpoint determinants are 263/288,2183/2304,31/32,9047/9216.
The checker derives every entry from `9*12^(d-1)*(1-w_d^2)`,
the pi enclosure and the trace bound, without marginal factorization.

For D7--D10 use36 amplitude cells of width1/80 and a positive
banded two-coordinate state. Each insertion includes its new square,
adjacent and distance2 density products and strip; copy and close at
each retained distance, then continue the interior state. Nonnegative
sums are divided downward at scale10^18. Time endpoints are
(j/24)^2,j=2,3,4,5,6,7,8,10,12,14,17,20,23, with the tangent power
integration described below. All factors and determinants are derived
from the stated precision, pi enclosure and trace bound.
The current exact sums pay the retained global floors. Multiplication
by their spectral constants adds the checked global increment
65797/200000000 to K. No dirty-window debit is introduced.

The positive cell contractions use only ideal independent coordinates;
physical residuals are centered, correlated and paid by Gaussian
correlation after the common quadratic tilt. No independence between
physical fibre images is assumed.

On each spatial box freeze every strip at a sufficient corner. For the
remaining joint density exp(-Q), Jensen gives its integral at least the
box volume times exp(-E_uniform Q). For an affine form, E_uniform Q is
the squared form at the cell means plus its coefficient-squared cell
variances. This is applied to the WHOLE sum of density exponents even
when their factors share coordinates. Only then factor the scalar
exponential for contraction. The adjacent convolution has two independent
uniform coordinates; both variances are included. Endpoint reflection
and all old label multiplicities remain unchanged.

D2 no-arm uses640 cells of width1/800 on [0,4/5], time endpoints
(j/160)^2,j=1,...,159, and the original sign margins. On a time panel
[L,R], freeze its spatial row at R. For the global joint side minimum
min(g(t)A_t,B_t), g(t)=(1-t/4)^(-1/2), its squared row loses at most
(4-R)/(4-t) relative to R. Above1/2 the root cases coincide and no
compensation is needed. With Delta=R-L, a=3/[2(1-R)],
k=(mu-1)/[2(1+(mu-1)R)] and b=1/(4-R) below1/2 (otherwise0),
the exact lower panel multiplier over the old right rectangle is

    1+(a+k-b)Delta/2+(ak-ab-bk)Delta^2/3-abkDelta^3/4.

It integrates the product of the two positive convex tangents and
1-b(R-t). Clean D2 has b=0,mu=1.04905; global mu=1.313.
Every spatial floor and division precedes the joint min/square as
before. Carry-free integer convolution only evaluates those exact sums.

For D3--D10 and larger even arms, each 2j-arm supplies (1-t)^j.
Freeze strip rows at R; integrate the positive powers of1-t exactly.
The tangent of (1+(mu-1)t)^(-nu) at R is a lower bound. Under each
positive decreasing power weight, mean(R-t)>=(R-L)/2, giving the
checked positive tangent improvement. Integer/half-integer powers use
exact antiderivatives and directed square roots. Global/local arms of
size6 at D2/D3 use the same first-step clean window and global defect
bounds thereafter, hence unchanged dirty counts. Their six extra
contrast diagonals are <=1 and the same determinant payment applies.
The checker keeps baseline, size2, and at-least-one-size4/6 classes
disjoint; it never adds the baseline twice.

The D3--D6 contraction may first round strip/endpoint weights downward
at binary scale2^60; all sums stay positive and this only lowers the
certified integral. Its paid floors remain the same.

All time endpoints are below1. Omitted tails pay zero. Source errors,
atomless boundaries and strict floor gaps permit the uniform finite transfer.

## 10. Signed and nine-key probabilities on the actual field

Here eta=.0009,theta=eta/2,T=.0245. Family Section4 proves the
actual-source path coefficients k2=109/25,k3=561/50,k4=586/25 and
opposite-branch costs3.67,7.85,17.8,14.8,10.53,26,28.9 by concavity,
two-mass reduction and complete directed radical covers, including
atoms and zero masses. The checker recalibrates all means, covariances, densities and SDs at
these eta,T.

### Signed drivers and exact reconstruction

Use the seven-key stars in Section6, with every inspected edge<=eta.
Put O=Y_o, R=2Y_u-Y_a-Y_b, Q=sqrt(3)(Y_w'-Y_w). Condition on O=t+x
and use precision coordinates (R,Q)-gO=K(y,z), x,y,z>0, where
g=Cov((R,Q),O) and K=Cov((R,Q)|O). The density is

    sqrt(det K)/(2pi)^(3/2)
       *exp(-[(t+x)^2+(y,z)K(y,z)^T]/2).

The controls are

```
|g_R|<=2eta, |g_Q|<=sqrt(3)*3.67eta,
K_RR>=1-8k2*eta-(2eta)^2,
K_QQ>=1-6k2*eta-(sqrt(3)*3.67eta)^2,
|K_RQ|<=sqrt(3)*11.05eta+|g_R*g_Q|,
K_RR<=1+2k2*eta, K_QQ<=1.
```

The checked determinant and inverse gaps are strict. For the six score
groups a,b; a1,a2; w; w1,w2; u; w', use signs + for below-threshold
scores and - for above-threshold guards. The ideal y coefficients in
the corresponding upper signed means are
`(-1/6,-lambda/12,lambda/6,1/2,-1/3,-lambda/6)`; the ideal z
coefficients are `(0,0,-1/(2sqrt3),-lambda/(4sqrt3),0,-1/(2sqrt3))`.
Add the one-sided raw covariance errors

```
y: (2k2,17.8,2,2k2,2k2,2k3)*eta,
z: sqrt3*(7.85,14.8,k2,10.53,1,k2)*eta,
```

and the root regression costs c_v^upper |g_i|. The x coefficient is
the appropriate signed upper covariance with O; threshold bias is
(Cov(Y_v,O)-1)t with the same sign. Round upper coefficients up.
Precision coordinates require no inverse-K correction in these means.
Worst affine corners give six positive available margins f,g,m,e,u,v.

After subtracting these means let X,A_a,D_u,D_w be centered wave
residuals and D1,W,H the centered half-differences at a-b, a1-a2,
and w1-w2. Exact star reconstruction gives deviations

```
a,b:     X/3 +/- D1,
a1,a2:   sqrt2*X/3+sqrt2*D1+A_a/2 +/- W,
u:       X/3,
w,w':    sqrt2*X/3+D_u/2,
w1,w2:   M +/- H, M=X/2+sqrt2*D_u/2+D_w/2.
```

In the second line the selected a has +D1; the sufficient private
bound uses its absolute value. The centered middle score is
J=sqrt2*X/3+D_u/2 at w,w'. The root guard is unilateral
X>=-3u. All reconstructions are linear identities before absolute values.

For ordinary private repair choose a centered strip |A_a|<=r_a and
one convex trapezoid in (D1,W) contained in

    |D1|<=f-r/3, sqrt2|D1|+|W|<=g-sqrt2*r/3-r_a/2.

The conditional pair covariance/Schur lower bounds follow by projecting
first on O and then the two drivers, using their exact majorants. As in
`family.md` Section8, sum its four quadrant densities, use cosh>=1,
and the directed Schur/normalizer bounds. Integrating 1-z<=exp(-z)
over the trapezoid uses its area and second moments; negative lower
values are replaced by0. Both coordinates may be correlated.

The alternate shifted private variable P=D1+sqrt2*A_a/4 gives child
deviations sqrt2*X/3+sqrt2*P +/- W. Sufficient private constraints are

    |A_a|<=r_a, |P|<=f-r/3-r_a/(2sqrt2),
    sqrt2|P|+|W|<=g-sqrt2*r/3.

For the raw positive private defect row, S<=3eta,
Cov(A_a,D1) lies in[-S/2,S/4], and |Cov(A_a,W)|<=S/(2sqrt2).
Its root and residual-driver row norms are bounded by the exact radial
ratios. Residual driver l1 bounds for A_a,D1,P,W are respectively
25eta,25eta,34eta,59eta; the two-driver inverse norm is<50/47.
Before projection the P variance is between

```
1/12-(k2/2+3sqrt2/4)eta
and 1/12+(3sqrt2/8+3(3-2sqrt2)/8)eta,
```

and |Cov(P,W)|<=((1+k3)/4+3/8)eta. Subtract the root and driver
projection costs for lower variance bounds; add them for absolute
cross-covariance uppers. The checker verifies its Schur and density
constants at eta=.0009 directly. Each private recipe is one centered
convex symmetric event; no alternative probabilities are added.

For the terminal pair use the centered convex diamond D={|M|+|H|<=e}.
Let c=3-2sqrt2,z0=3eta. An upper SD for M is
(1+sqrt2/2)sqrt(c*z0), and Var(H)<=1/12. Raw row ratios give
|Cov(M,H)|<=(5sqrt2/8)z0. Enlarge the two raw diagonals to their directed
SD squares, keeping the raw off-diagonal; this dominates the conditional
covariance as well. Gaussian measure of a centered convex symmetric set
decreases under this enlargement (Anderson's inequality). If the resulting
correlation magnitude is at most r0, sum the four quadrant densities,
drop cosh and the reciprocal determinant factor, and integrate the
lower density

    2/(pi*a*b)*exp[-(x^2/a^2+y^2/b^2)/(2(1-r0^2))].

Its total mass is1-r0^2. Increasing both normalized independent folded
SDs to a,b decreases x+y<=e, so the checker pays
(1-r0^2) times the independent diamond probability. Bins of width.001
use Jensen mean-square folded masses and their right endpoints for the
remaining half-strip radius; e is rounded down on the .0005 grid.
The old ending recipe remains an alternative.

### Adaptive core and reflected asymmetric margins

For root radius r use the symmetric convex core
R(r,m)={|X|<=r, |D_u+(2sqrt2/3)X|<=2m}.
Projection on the root gives conditional variances at most
c*z0-z0^2 and c*z0-z0^2/2 for X,D_u, respectively. Enlarge the
diagonals to checked SD squares sx^2,sy^2, retaining the covariance.
This dominates the conditional centered covariance. To cover its unknown
correlation, divide[-1,1] into32 panels. In a panel use the maximal
absolute mean multiplier |2sqrt2/3+rho*sy/sx| at its endpoints and the
largest conditional SD sy*sqrt(1-rho^2), at the point nearest0.
On folded-root bins of width.001 use the upper root endpoint in that
mean and Jensen's exact mean-square marginal mass. Only bins with
2m larger than that mean contribute, so both normal CDF arguments
are positive and increasing either chosen upper loss lowers the mass.
Take the minimum over panels. CDF arguments may be rounded downward
to multiples of5 in their integer index. Singular endpoints follow by
continuity. Tabulate r=.010,.016,...,.046 and m in steps.001.
This covers every correlation sign.

Let a=min(r,3u), L=min(m,v), H=max(m,v), where m,v are the two middle
margins. The actual core is I=[-a,r]x[-v,m] in coordinates(X,J).
For any fixed symmetric private/terminal event E, central symmetry gives

    P(I intersect E)=1/2[P((I union -I) intersect E)
                            +P((I intersect -I) intersect E)].

The union contains R1=[-r,r]x[-L,L] and
R2=[-a,a]x[-H,H]; their intersection equals I intersect -I.
Consequently P(I intersect E)>=1/2[P(R1 intersect E)+P(R2 intersect E)].
Each term uses Gaussian correlation with the SAME private debits for
the full r. Round a down to an available root radius and L,H down to
the middle grid. If no available radius<=a, the second term pays0.
The downward integer average is a valid core floor. No covariance sign
or residual independence is assumed.

The following further sufficient recipes may accelerate the replay while
keeping the same required signed floors. In one reflection choose X<=0:
its root terms in both private lines are favorable and need no debit;
the root guard only needs r<=3u. In the opposite branch X>=0, any root
radius passes the root guard while the full private root debits remain.
Each pays half the probability of its fixed symmetric sufficient event.

There is also a two-sign combination. For r<=3u let R be its symmetric
core; A a private event paying the worst root debit; B a private event
omitting only the favorable X terms; D the terminal diamond above; and
S={|H|<=e}. On R intersect A intersect S, one reflection has M<=0 and
succeeds. On R intersect B intersect D, one has X<=0 and succeeds.
On their intersection, A and D make BOTH reflections succeed. Orbitwise
therefore the success indicator averaged with its negative is at least
half the sum of those two event indicators, even if A,B are not nested.
Gaussian correlation yields the lower bound

    (1/2)P(R)[P(A)P(S)+P(B)P(D)].

All these events are centered convex symmetric subsets of the same
Gaussian residual space. Each corner chooses one proved recipe for the same original success
cylinder; these alternatives change no name count.

The Gaussian success probability is log-concave in its affine mean
parameters by Prékopa's theorem. Within a driver box, interpolate its
corner lower probabilities log-concavely. Integrate the resulting
log-affine bound with the density, using its exact cell moments and
Jensen, or the checker’s equivalent geometric corner-product floor.
Private/core recipes may differ at corners because they only lower-bound
that SAME success cylinder. Every normalization, width and corner product is rounded downward. Fixed positive
integer convolutions/packed lanes only regroup exact sums; their carry
bounds and lane widths are checked. No cached search score is a mass.

At t=0,-T,+T the required signed floors, at scale10^14, are
225627942,197515123,207450756. The replay proves at least these floors
on disjoint driver boxes. Negative-margin boxes and all omissions pay0.

### Direct precision means for strict and relaxed nine names

Use the thirteen-edge faithful stencil of Section6. Strict names have
all edges<=theta; relaxed names allow eta only on the four terminal
child edges. Keep the same guarded sign cylinder. Put O=t+x and raw
drivers D=(Y_b+Y_u-2Y_a, sqrt3(Y_w'-Y_w),sqrt3(Y_z'-Y_z)).
For g=Cov(D,O), K=Cov(D|O), use D-gO=K*xi,xi=(y,z,h). Its density is

    sqrt(det K)/(2pi)^2 exp(-[(t+x)^2+xi^T K xi]/2).

The root, sibling, mixed-driver and opposite-ending inequalities in
`family.md` Section4 give lower diagonals
`(1-8k2*theta-g0^2,1-6k2*theta-g1^2,1-6k2*theta-g2^2)` and off-diagonal
magnitudes `(sqrt3*24.55theta+|g0g1|,sqrt3*24.55theta+|g0g2|,
35.3theta+|g1g2|)`, where
`|g|<=(2theta,sqrt3*3.67theta,sqrt3*3.67theta)`.
The comparison matrix is positive definite by strict diagonal dominance.
With lower diagonals r,q,q and magnitudes c,c,f, determinant monotonicity
under adding a positive diagonal gives

    det K >= r*q^2-r*f^2-2q*c^2-2c^2*f.

The checker pays normalizer247/10000. The upper quadratic retains actual
upper diagonals `(1+2k2*theta,1,1)` and these cross magnitudes.

For EACH score use its direct conditional mean

    E[Y_v|O,D]=c_v O+(Cov(Y_v,D)-c_v*g^T)xi.

Sign this expression for below/above guards. A positive raw covariance
coefficient uses its upper excess; a negative uses its lower excess.
Its opposite-branch contrast costs before sqrt3 are7.85theta at depth1,
14.8theta at depth2,26theta at depth3, or28.9theta at relaxed endings.
Charge the possible doubled terminal edge on its actual path. Regression
costs c_v^upper |g_i|. Threshold bias is bounded by(1-w_d)T.
Take the smallest margin over each score group. This avoids replacing
actual signed means by the older absolute inverse-K envelope.

The centered score identities and strips remain the original ones:
six wave-defect strips, four half-differences, and an optional centered
root/ending trapezoid. The two Q-driver endpoints have the SAME centered
loss after conditioning. An ending/root half-pair has covariance below
the recomputed bound and diagonal<=1/12, so the checked enlargement
.09 I dominates it; integrate one selected trapezoid with slope1.5.
Other halves/defects use centered strip probabilities and Gaussian
correlation. It is invalid to add both alternative ending-pair masses.
All SDs, Schur bounds and radius masses are regenerated at eta=.0009.

The code recomputes every supplied nonzero choice. Box density uses the
upper quadratic averaged on that box, then Jensen and downward integer
exponentiation. Reflection z<->h includes both branch labelings; the
envelope applies to each separately and does not assume exchangeability
of the actual covariance. The relevant reflection multiplicity and
root-coordinate volume are retained. Codes with any invalid margin
pay0. A reduced witness may omit any further codes: only nonnegative
mass is lost. In the submitted replay the strict and relaxed lower
floors are each at least1673/10^9. This exceeds respectively
u+v9+v11 and u in Section13 below. The unreduced witness gives the
stronger floors2336449/10^12 and59783/31250000000, which are unnecessary.

Nonsingular driver blocks, compact covariance controls and strict floor
gaps permit the common finite transfer for all names and |t|<=T.
## 11. Capped Price energy and simultaneous exceptional charges

Mark original complete-source edges high when D>eta and medium when
theta<D<=eta, and set D'=min(D,eta). The numerical global coefficient is
K0=1297660793724599979211/25000000000000000000000. The seven raw local
capacities, common multipliers and three exceptional-class dirty debits are

```
name       raw cap             multiplier   dirty b1,b2,b3
adj        61075651/10^10       1            29,29,29
D2         6354307/2500000000   2lambda      0,18lambda,18lambda
D3         993819/2500000000    19           0,327,327
D4         534053/10^10        18lambda     0,522lambda,522lambda
D5         5331/625000000      120          0,0,6040
D6         3957/2500000000     92lambda     0,0,7788lambda
stratD2    4288307/10^10        2lambda      0,18lambda,42lambda
```

Here lambda=141421/50000<2sqrt2. Retain fractions
(8353353717,10^10,10^10,5606863503,10^10,5468685896,10^10)/10^10
and call the spent raw amounts x1,...,x7. The SAME spends enter every
common credit and every dirty debit. Global floors are entered once;
local caps are local-minus-global replacement floors. Arms of length4/6
exclude baseline/size2 cases. The adjacent global-to-clean transfer
71931/10^10 is removed globally before being included locally.

For Ej={D>Hj}, H=(.0027,.006,.007), each dirty tree-name test is lifted
together to a high-girth cover and divided by covering degree. Projection
cannot increase its number of sites, even with aliases/short cycles.
The residual row R_v and S_v=sum_(e incident v)D_e obey Section7.
For branches across w--u put

```
U1=lambda0*b_u-b_w+R_u,
U2=6b_u-lambda0*b_w+lambda0*R_u+sum_(z~u,z!=w)R_z,
lambda0=2sqrt2.
```

Expanding branch inner products gives the positive source terms and
independent dirty-site counts below. A site q*S_v costs at most3q*eta.

```
D2: lambda*S_w                              6 centres
D3: 9D_ab,5S_a                              13 edges,14 sites
D4: 6lambda*S_w,3lambda*S_u                  14,30 sites
D5: 44D_ab,24S_u,14S_z                       29 edges,30,62 sites
D6: 24lambda*S_w,14lambda*S_u,8lambda*S_z     30,62,126 sites.
```

These entries follow by expanding U1/U2, subtracting ideal products,
and collecting each retained residual with its full nonnegative bracket.
For D3 that bracket is5b_a+lambda0 R_a; for D4 it is
3lambda0*b_u+6R_u+lambda0 sum R_w. For D5 the outward coefficients are
24,14 and for D6 they are14lambda0,8lambda0. Remaining products are
nonnegative up to uniform o(1), since b_u>=alpha*b_w. Summing endpoint
incidences yields exactly the table's common multipliers; the displayed
site counts yield its dirty debits. StratD2 uses6 centres at H2 and14 at H3.

Adding an incident exceptional edge introduces at most8 adjacent sites;
D2 introduces2, D3(4,4), D4(4,8), D5(8,8,16), D6(8,16,32).
Spanning-tree union bounds save per additional edge twice the resources

```
S=(21/2)eta*x1,
Tsrc=eta[lambda(6x2+189x4)+(231/2)x3],
t5=2220eta*x5, t6=2886lambda*eta*x6.
```

StratD2 saves4,10,6 centres on an incident middle/middle, top/top,
middle/top pair, respectively. Select disjoint pairs inside eligible
higher blocks, counting every unselected edge separately. This pays
12r7,30r7,18r7 with r7=lambda*eta*x7, even if dirty site sets intersect.
Retain D6 connected union savings on pure top components and all-top
higher blocks inside mixed components; these are different components.

Outside energy is only discarded K0(D-D'). Select isolated E2 edges,
isolated lower E1 edges, and one maximum-D edge from each pure higher
dyad; they form a matching. The paired sigma law supplies

```
Ij=K0(sigma*Hj-4eta)/2,
Jjk=(K0/2)[sigma*Hj-kH1-(4-k)eta]_+,
Jjstar=K0[(sigma-1)Hj-3eta]/2.
```

For Jjk transfer the unused K0(D-H1)/2 from each lower neighbor;
each has at most two accesses. Outside E1 each edge has at most one
access at either endpoint, so all half-charges fit its discarded energy.
Eligible isolated edges/pure higher dyads unshared with another access
retain the other half; shared components receive adjacent-window savings.
Indeed an access joining different exceptional components gives13 extra
overlapping adjacent names. Root a forest at ineligible components, or
use a spanning tree on p>=2 eligible components, yielding >=13p/2 extra
names. Thus the old first-class bonus min(I1,(13/2)eta*x1)=I1 survives,
and all higher singleton/dyad components can use
bsmall=(13/4)eta*x1, checked<=min(I2,I3,J2star,J3star). This is ONE
forest account, distinct from within-component savings.

Let m=10^-12,

```
Rj=K0(Hj-eta)-eta*sum_g b_jg*xg+min(S,Ij),
p0=R1+S-I1-m, P=p0+S,
tj=m-Rj-Tsrc (j=2,3),
v=max(t2-6r7,(t2+t3)/2-9r7,t3-t5-t6-15r7).
```

All branch/clamp signs are regenerated exactly. At this cutoff p0>0,
so pure lower dyads need NO additional dyad source payment. Pure lower
components are paid by2p0>=0,3p0+S>=0,P>=0; isolated lower edges by
R1+I1>=m. Higher singletons/dyads use respectively

```
Ij+bsmall-S-Tsrc>=tj,
J2star/2+6r7+bsmall/2>=t2,
J3star+18r7+bsmall>=t2+t3,
t5+t6+J3star/2+15r7+bsmall/2>=t3.
```

Larger all-middle, mixed, all-top components reduce to the size3 tests
S+Tsrc-3t2+12r7>=0,
S+Tsrc-t2-2t3+18r7>=0, and
S+Tsrc-3t3+4(t5+t6)+30r7>=0. Additional edges are nonnegative because
max(t2,t3)<=S+Tsrc. A mixed nonsingleton higher block has resource
>=2(S-v), by its three class compositions and those pair refunds.

For a mixed E1 component let l count lower edges, q_jk higher singleton
blocks of class j with k lower incidences, and r nonsingleton blocks.
Connectivity gives q+r>=1 and sum k*q_jk+r<=2l. Its remaining amount is

    F=sum(S+Jjk-Tsrc-tj)q_jk+2(S-v)r+P*l-2S.

For l=1 use the five (q,r) cases(1,0),(2,0),(0,1),(1,1),(0,2).
Their checked inequalities are p0+Jj1-tj-Tsrc>=0,
P/2+Jj1-tj-Tsrc>=0,P/2-v>=0,
p0+2S-Tsrc+Jj1-tj-2v>=0,(p0+3S)/4-v>=0.
For l=2 require both a_jk>=0 and a_jk+2p0>=0 for k=1,2,
a_jk=S+Jjk-Tsrc-tj, and both2(S-v)>=0,2(S-v)+2p0>=0.
This pays the case for either sign of p0.

For l>=3 set beta=p0/2+S(l-2)/(2l) and
mu_jk=S+Jjk+k*beta-Tsrc-tj. Then exactly

    F=sum mu_jk*q_jk+[2(S-v)+beta]r
                       +beta[2l-sum k*q_jk-r].

The lower exceptional edges are in no physical triangle or square:
the unicyclic covers give excess floors
(2+2sqrt6)/7-2sqrt2/3>.04 and
(2sqrt170+2sqrt80)/45-2sqrt2/3>.03. Further fibre coarsening increases
these affinities; fixed finite-source continuity preserves the strict
bounds, whereas lower edges have D<=.006. Thus two higher singleton
blocks share at most one lower link. Consequently at l=3 at most one
degree3 singleton occurs; at l=4,5 at most one degree4 occurs (at l=4
all other degrees are<=2); at l=6,7 at most1,2 degree4 singletons occur.
Two such degree4 blocks need at least7 distinct lower links; three need
at least9. These are source-driven cycle exclusions, not arbitrary
marked-graph assumptions.

For each l=3,...,7 the checker permits only those exceptional negative
coefficients. It bounds their loss per block by max(0,-mu_2d,-mu_3d),
and pays from the remaining2l-d*m slots at the minimum of beta,
2(S-v)+beta, and mu_jk/k for all eligible other degrees. It verifies
this remaining payment exceeds the total loss for the largest allowed
multiplicity m; every smaller multiplicity is easier. For l>=8 use
beta=p0/2+3S/8 and nonnegative coefficients; the remaining constant is
S*l/4-2S>=0. This exhausts all connected exceptional components.

All these tests are exact in the replay. After paying m per exceptional
edge and the vanishing source debits, they deliver

    K=sum_global K0+sum_g multiplier_g*xg+65797/200000000,
    K*sum_e D'_e-o(n), eta*K=A.

The D7--D10 increment needs no dirty charges. Source identities may
constrain several pools, but each global/local floor, discarded-energy
half-charge, selected pair/forest union refund and capped endpoint
payment has been allocated once. The finite source/mesh errors are
chosen beneath the final strict graph margin, uniformly in the graph.
## 12. Ordinary marked energy, direct runs and two medium payments

All amounts in this section are normalized by A. Put

```
sigma=2021/1250,a0=sigma-3/2,c=a0/4,h=(sigma-1)/4,
rho=2X/A, s=Z/A, g=(2/3)(s-c-rho/4), B=c+g/2,
k=33/100,d=(sigma-1)/8,nu=(1/4+rho)a0,nu2=s-1/8,
omega*=min(1/2,2s), kappa=59599/10^6.
```

The checked domain is .4<rho<.5,1/8+rho/4<s<1/4+rho/4.
The ordinary table follows by assigning high ends1/2, medium endsx/2,
shared low endsx/2+s and boundary low endsx+rho*min(x,1/4).
The weighted source triangle bounds pay their U2,U3,J propagation costs;
`family.md` Section6 gives the full derivation. Consequently

    K sum D' +2Z*sS+X*U2+(3X/2)*U3+(X-A/5)*J

pays A times the following vertex coefficients. Here n_ab counts vertices
with a high,b medium edges, sS counts low marked--marked links, Ui counts
unmarked vertices with i marked neighbors, and J counts U1--U1 edges.

```
01 1/2+rho/4+B   02 1/2-g       03 3/4-3B
10 3/4+h+s       11 3/4+h-B     12 1-h-2B
20 sigma/2       21 5/4-2h-B    30 3/2-3h.
```

For each H1--H1 high edge (H1=10,11), the old fixed paired bonus is
already removed before its off-edge residual is used below. Everything
uses capped high baseline1 and unchanged low/medium terms; no discarded
high surplus is allocated. The one/two-leaf mixed groups retain their
old fixed fees. Zero/one-H1 type20 centers have count reserves
V0=(2-sigma)/2=.1916,V1=(5-3sigma)/4=.0374.

For a two-H1-leaf type20 center let original high excesses1+se,1+sf,
center low l, leaf off-sums1/2+ye,1/2+yf, Y=ye+yf,S=se+sf.
Its ordinary mixed-star remainder is

    P=Y/2+f(l)-a0=R/2+[f(l)-l]+(sigma-1)S/2,
    f(l)=l/2+min(l/2+rho*min(l,1/4),s), R=Re+Rf.

The replay's complete recursive rectangle cover proves the critical
scalar inequality G(p)+G(q)+(4/5)phi(r)>=.0596(phi(p)+phi(q)),
p+q+r=3 in[1/2,2]. Here G is the paired residual at the sharp coefficient Cstar;
the actual residual at sigma is larger by (Cstar-sigma)*D>=0. It uses
G(t)>=(16/15)phi(t)(t-1)^2/t, rationalized phi interval bounds, and
F(t)=G(t)-.0596phi(t) whose F(t)-t^2 is concave. The lower chord
min(F(lo),F(hi))-(hi-lo)^2/4 and the quadratic bound are combined.
Recursive disjoint rectangles cover every feasible point; the checker
accepts1443 rectangles and rejects59 outside the domain. Finite ratio
clamping and the root-atom debit yield R+.8l>=kappa(xe+xf)-o(1).
For an actual boundary l<=1/4, f(l)-l=rho*l>=.4l, hence P>=kappa.
For boundary l>=1/4, P>=rho/4; for a shared port P>=s-a0.
Both exceed kappa. This uses the actual port branch, and pays one P.

The medium count retains B for every type11 stub not facing01.
Remove the saturated anchors02,03,12,20,21,30 from the high/medium
marked graph. The remaining type11 components are alternating runs
or cycles. In a run let N be its type11 count and e,t10,hA,mA count
its01,10,high-anchor,medium-anchor ends. Then e+t10+hA+mA=2.
Its original extra pool X_R consists of quarter-lows, opposite external
medium halves and H1--H1 residuals AFTER actual medium fees. It is
nonnegative. A mixed01--02--11 group retains G-c>=max(0,s11/2-c);
that single actual group belongs to its one run. The unchanged medium
count gives B(N-e). A high10 end retains nu BEFORE both its recovered
fee and its positive amount b. Thus its basic pool after b is

    B(N-e)+(nu-b)t10-o(N).

The direct row payment C and anchor debit p are valid under

```
0<C<min(B,p,nu-b,d,(B+omega* a0)/2,(2/3)(B+k/16)),
0<b<min(p,(s+(1/2+rho)a0)/2),
epsilon<min(.1916,.0374+p,kappa).
```

The last b condition is the disjoint10--10 concavity pool in
`family.md` Section12;10--11 uses its full nu before C+b.
For e=0 the basic count pays CN. For e=1 with a high end, the
shortest N=1 uses C<=p or C+b<=nu; longer runs have slope B-C>=0.
For e=1 with a medium end, N=2 pays2d+B-c; N>=4 pays
B(N-1)+k(N-2)/16. For e=2, N=2 pays2d; N>=6 pays
(B+k/16)(N-2). At N=4 only, reserve the unused full low remainder
above the original quarter share; g0(l)>=omega* l gives
2(B+omega* a0). These N4 high neighbors are all11, so this reserve
pays no10 fee. Cycles have B per vertex and nonnegative X_R.
This proves the direct C menu, replacing the old alpha mixture.
At anchors12,21,30, actual11 and10 high stubs are disjoint; -p per
high stub pays their run debit or b, since b<=p. At type20 split
-epsilon+2p by its actual11,10 leaves, reserving each run debit once.

For the new ordinary medium payments t at01,z02 at02, keep the full
B run count. The exact unused medium pool after old fees is
(g/2)Pm+gN0+(g/2)N1+(kappa/2)N2+B(Mmed-Rmed), where Nj counts02
centers with j neighboring01 leaves and Pm counts01--01 pairs.
The two-leaf kappa/2 follows from the SAME weighted scalar star with
medium cutoff1/2, after both old c fees; shared and large-boundary low
branches give larger amounts. These stars have no11 medium neighbor.
Subtracting the targets leaves coefficients
(g/2-2t),(g-z02),(g/2-z02-t),(kappa/2-z02-2t), respectively.
No B stub or mixed01--02--11 remainder is taken twice.

The remaining01 vertices are ends of the runs, or face saturated
medium anchors. For the short01--11==10 run its paired law retains
beta+l>=a0+1-m. The untouched fee pool is
F=g0(beta)+g0(l)-l/4>=min(nu+(1-m)/4,nu2).
Concavity reduces this piecewise-linear minimum to square/diagonal
vertices. Its own01 low surplus above its old f(1/4) baseline is
J01>=(m-1/2)/4. Thus F+J01>=min(nu+1/8,nu2), paying C+b+t.
This own-low surplus is used only here. Longer high10, medium-ended
and two01-ended runs use the preceding exact pools with one or two
t targets subtracted. Their tails have nonnegative slope.

For type20 high ends put r=p-[3C-2B]_+. A nonshort run with a far01
end pays its t first and releases only r-t; without a01 end it
releases at least r at each high end. A short end costs C+t.
The sufficient remaining guards used by the checker are exactly

```
t<=g/4, t+z02<=g/2, 2t+z02<=kappa/2,
t+C+b<=min(nu2,nu+1/8), t+3C<=2B+nu-b,
t+2C<=2d+B-c, t+4C<=3B+k/8,
t+C<=d, t+2C<=B+omega* a0, t+3C<=2B+k/8,
t<=r,
t+epsilon+C<=V1+2p,
t+epsilon+C+b<=kappa+2p,
2t+epsilon+C<=kappa+p+r,
2t+epsilon+2C<=kappa+2p.
```

The last line is two SHORT arms: use only ordinary P>=kappa and
nonnegative X_R. No .0782 or full double-Q pool, far-end reserve,
new dyad law, nonnegativity refund or extra low remainder is added.
The preceding line is short+long and counts BOTH near/far01 targets.
Other type20 cases use their old reserves and released actual debits.
If a run joins two centers it has no01 end and releases p at each;
there is no repeated far target. All component arms stay paid even
if an original high edge is clipped or exceptional: these arguments
spend only capped low/medium endpoint shares.

At saturated types03,12,21,30 charge t per actual medium01 neighbor
or per actual high-anchored run carrying a01 target. There are at
most a+b=3 such incidences at any such vertex. Thus row shifts are

```
01:+t,02:+z02,03:-3t,10:+b,11:+C,20:epsilon-2p,
12:-p-3t,21:-2p-3t,30:-3p-3t.
```

This partitions every01 mate/11 run and preserves the q01=2 fees.
Each residual, medium half, low share, B stub and anchor debit is spent
once, including shared low links. Bounded ownership/slopes give o(n)
error even on long runs.
## 13. One threshold quota and the exact marked ledger

Before any advance, the actual original field loses at least
b*min(t^2,r^2)*Pelig from its zero-threshold crossing-cover expectation,
where r=101/5000,b=5329214440081/10^16 and Pelig counts directed eligible
distance2 paths. This difference is ADDITIVE to the zero-threshold Price
comparison: it telescopes on the original base field, not on interpolation
pairs. The seven-edge path stencil is faithful by the covariance bounds.

Here is the actual-source argument for this deficit. Condition endpoints
u,v of u-w-v to s, and put Y_w=m*s+Z. The checker proves variance bounds
.0279<Var(Z)<.03328, endpoint residual SD<.0216, sibling-half SD<.29,
and centered endpoint-pair correlation magnitude<.191. For either outward
child a at u set L_a=Y_a+Y_w/2-sqrt2 Y_u and Z0=Y_w-t0(Y_u+Y_v),
t0=4sqrt2/11. The critical source gives
Cov(L_a,Z0)<=D_au+D_uw+D_wv. To verify this source-specific inequality,
partition a physical fibre by its nearest vertex of a-u-w-v. Normalize
P_u=1 and put P_a=1/2+3A/2,P_w=1/2+3W/2,P_v=1/4+3W/4+3V,
A,W>=0,A+W<=1,0<=V<=W. For d_ij=sqrt(P_i P_j)-I_ij, the desired
slack has coefficients(1+t0,1+sqrt2+t0/2,1+t0/2,t0,-1,-8/11)
on(au,uw,wv,av,aw,uv). Four times its negative Hessian is congruent to

```
[(1+t0+t0*q-s)/r, 1, -t0;
 1, (1+sqrt2+t0/2+(1+t0/2)q-r)/s, -(1+t0/2);
 -t0, -(1+t0/2), ((1+t0/2)s+t0*r-8/11)/q],
 r=sqrt(P_a),s=sqrt(P_w),q=sqrt(P_v).
```

The fixed69-box recursive cube cover(A=(1-W)p,V=Wz) checks all Sylvester
minors with directed radicals; its determinant lower exceeds.000709.
Thus the slack is concave and vanishes on the simplex vertices, hence
is nonnegative. Root atoms and finite-profile limits are uniform. Actual
endpoint regression differs from t0 by<.0045; the checker pays joint
covariance3.06eta and verifies the adverse Z slope<.099 at eta=.0009.

On a folded-Z amplitude cell [l,h] in[.04,.60], the two endpoint
centered convex diamonds |D'|+|R'|/2<H have sufficient margin
H=(.5-.099)l-(.101+.00443/2+.00218)r. A middle residual strip of
radius.1 pays its third-neighbor sign. Gaussian correlation bounds the
two endpoint diamonds and middle strip together by.999*p(H)^2.
Covariance enlargement to SDs.29,.0216 and |corr|<=.191 is valid.
Conditional on standardized folded defect z, pairing contrast signs
pays cosh>=1; its weighted density is at least
.7978*exp[-z^2/(2(1-.191^2))]. Bins of width1/40 through3.5 keep the
right endpoint in H-.0216*z/2 and Jensen's mean square in the density.
Amplitude bins of width1/400 use both variance endpoints for their
Jensen density lower bound. The resulting sign stencil forces mixed
difference<=-1 for every exterior, by the distance2 matching argument.
The exact derivative factors .3985*(4039/25000)*.71/2 give the stated b.
The scalar source guards are all checked at eta=.0009,r=.0202.

For each fixed signed naming, its ten original score inequalities form
a convex cylinder jointly in keys and threshold. Prékopa's theorem makes
its success probability log-concave in t, even for singular Gaussians by
restriction to their range. Let q0=225627942/10^14,qT=197515123/10^14,
L=3327/25000. The checker verifies qT*sum_(j=0)^6 L^j/j!>q0, so
q(t)>=q0 exp(-L|t|/T). For x<=r use the lower quadratic

    q0-q0*L*x/T+[b+q0*(1/2-L/6)*(L/T)^2]x^2,

obtained from the cubic Taylor lower bound, and complete the square.
For r<=x<=T use qT+b*r^2. Both exceed q=218978509/10^14.

Both eligible signed weight and Pelig are >=the SAME nonnegative count
M=max(0,6n-sum c7(a)n_ab),c7=(0,14,10,6). For the signed count this
follows exactly by writing i(v)=number of high-marked neighbors of an
unmarked v: L7=sum_(UU edges)(2-i(o))(2-i(u)), and the extra unique
high-partner weight is E01+E11. Their sum equals
6U-4rH+E12+4E22+6U3>=6U-4rH. For Pelig the exact count is6U0+2U1
>=6U-4rH. High boundary incidence substitution gives the displayed
marked bound for both. Thus the two successive expectation differences
can spend the SINGLE quota q with this marked debit, even if the affine
count is negative; taking M's positive part is legitimate.

The remaining strict/relaxed eligible counts are derived in `family.md`
Section5 and retain their physical-star and overlap multiplicities.
Write j=a+b,c9=(0,35,22,9),c11=(0,45,28,11), and

```
cR(a,b)=(0,36,24,12)[j]+8(3-j)1[a>0],
cA(a,b)=c9[j]-2(3-b)1[a=0<b].
```

The lower count constants are9,11,12,9 for L9,L11,LR,LA respectively;
L9,L11 are alternate lower bounds on the SAME strict count, LR/LA on
the relaxed count. Their(sS,U2,U3,J) coefficients are respectively
(26,12,30,4),(34,16,40,6),(24,12,24,0),(22,12,24,2).
Use the exact weights

```
(q,u,v9,v11,gap,X,E)=(218978509,167231461,100,100,100,
                         1778578981,394844746)/10^14,
e=u-gap, Z=13v9+17v11+11u+gap,
(A*C,A*p,A*b)=(335739986,335740086,205463189)/10^14,
A*t=A*z02=29115760/10^14, weak weight=0.
```

The score payment is qL7+v9L9+v11L11+(u-e)LR+eLA.
Require signed-plus-deficit>=q, strict>=u+v9+v11, relaxed>=u.
The two nine floors1.673e-6 suffice strictly. Its sS coefficient is2Z.
The checker verifies U2,U3,J coefficients cover X,3X/2,X-A/5 and every
ordinary table coefficient plus Section12's shifts covers
q*c7[a]+v9*c9[j]+v11*c11[j]+(u-e)cR(a,b)+e*cA(a,b).
All nonzero source/marked guards have strict rational slack. Adding them
cancels each count debit/overlap allowance exactly once and gives

    Csave=6q+9v9+11v11+12u-3e=2818956503/10^14.

Therefore the advanced original cut has expectation at most
(U-Csave)n+o(n) for |t|<=T. The outer bound pays more for |t|>=T.
The reference quadrature and final check prove the strict margin

    gamma-(U-Csave)=3539099480467/191395798900000000000000>0.

Choose one common radial parameter, finite radius, positive noise,
threshold mesh/tails and concentration accuracy inside that margin;
then choose graph size. Section2 gives the limiting cutwidth with fixed
positive slack below gamma. Disconnected components concatenate;
components below the common size threshold have bounded width.
The separate actual-length circuit/family transfer now applies.
