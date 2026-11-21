# A marked Gaussian ordering with paired signed and nine-key Q bounds

For every sufficiently large simple cubic H on n vertices

    cw(H) <= (gamma-epsilon0)n,
    gamma=1000000000/7656280233,

for a fixed epsilon0>0. The separate `family.md` transfers this to one
fixed total uniform NP language with full-B2 circuit lower bound
`(9656280233/2000000000) N=4.8281401165 N` at every sufficiently large
actual integer length N, with zero subtracted remainder.

Set

```
eta=57/50000, theta=57/100000, T=277/10000,
U=1633094581/12500000000=.13064756648,
A=11569604236818386923776429/125000000000000000000000000000,
Csave=1791867437/50000000000000.
```

The manifest binds the replay; `nine-choices.txt` only selects boxes.
The two central bands use different allocations on ONE nested ordering;
all marked owners and current source menus are proved in Sections11--13.
All displayed parameters are current; the old .0009/.0245 calibration
is superseded by the source and probability checks below.

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

The critical-wave alternating-path calculation gives Gamma_25<U,
with U=.13064756648. Its complete path law, tails, derivative error,
Legendre brackets, directed weights and integer-rounding proof are
in family.md AppendixB, included in this packet. In particular the
reference crossing-cover expectation is at most U*n+o(n). No prior
external numerical result is assumed.

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
in family Section4 are also used. At eta=.00114 these uppers are<1
through the diameter6 of each inspected stencil, proving all names
physically distinct, including relaxed terminal edges.

The reverse path law is exact before taking the critical limit. On one
fiber put d_xy=b_x*b_y-I_xy. If z extends x--...--y, dist(x,y)=k,
the positive radial rows give
I_xz=alpha*I_xy+alpha^k*I_yz-alpha^(k+1)*b_y^2. Therefore

    d_xz-alpha*d_xy-alpha^k*d_yz
      =(b_x-alpha^k*b_y)(b_z-alpha*b_y)>=0.

Both factors are nonnegative by the pointwise row ratios. Summing and
induction give D_(v0,vd)>=alpha^(d-1)*sum_path D, including root atoms.
Uniform finite-row/noise errors are absorbed as above. The stronger
13/20 dyad used by the current marked centers is proved in Section12.

Let lambda0=2sqrt2,c=3-lambda0,R_v=sum_(u~v)Y_u-lambda0 Y_v.
Fibre triangle and the resolvent give a row differing by o(1) from
0<=r<=c*b_v. Put S=<r,b_v>=sum_(e~v)D_e+o(1). Then
||r||^2<=cS, and at distance j its covariance with a key lies between
sqrt2^(-j)S and sqrt2^j S, up to uniform errors. Conditioning a pivotal
interpolation endpoint at zero subtracts at least t*2^(-j)S^2.
Maximizing over S<=c gives c^2*h_j(t), where h0=1-t until1/2 then
1/(4t), and hj=1-t/2^j for j>=1. For a clean star with S<=3H the
bound is3cH-9tH^2/2^j. The stratified first-outward star has
S<=H2+2H3=.0232 at the current windows; its variance must use this cap.
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

By the displayed resolvent identity, the expectation with weights `P_v(p)` of the absolute
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

Family AppendixC contains the complete current proof of the alternating
path pivotal witnesses, joint middle precision and arms, threshold
integration, conditional quadratic tilt, D2 tangent panels, D3--D10
power integration and spectral incidence. Every later reference to
Section9 includes that appendix. The only independent coordinates
factored there are ideal tree coordinates; physical residuals remain
correlated and use the stated common tilt and Gaussian correlation.
At the current join it pays outer credit
17930855452501/500000000000000000. Its global/local Price witnesses
and simultaneous dirty-window debits are allocated exactly in Section11.
No additional external comparison theorem is substituted.

## 10. Signed and nine-key probabilities on the actual field

Use eta=.00114,theta=eta/2,T=.0277. Family Section4 proves the
path coefficients k2=109/25,k3=561/50,k4=586/25 and opposite-branch
costs3.67,7.85,17.8,14.8,10.53,26,28.9 by complete scalar fiber
certificates. Family AppendixA proves the NEW joint Q/root comparison.
All means, covariances, density normalizers, defect SDs and grade
ratios below are recalibrated at these current parameters.

### Signed drivers, source means and exact reconstruction

On the seven-key stars let O=Y_o,R=2Y_u-Y_a-Y_b,
Q=sqrt3(Y_w'-Y_w), g=Cov((R,Q),O), K=Cov((R,Q)|O).
Condition on O=t+x and use precision coordinates (R,Q)-gO=K(y,z),
x,y,z>0. Its density is

    sqrt(det K)/(2pi)^(3/2)
      *exp(-[(t+x)^2+(y,z)K(y,z)^T]/2).

Write d_i=D_(o,i), J=k2/2, alpha0<1/sqrt2 the directed lower value.
The raw R variance excess is -4D_ua-4D_ub+2D_ab. Section7's reverse
and upper two-edge bounds give

    Var R>=1-(4J-2alpha0)(d_u+d_a+d_b)-(4J+2alpha0)d_u,
    K_RR<=1+2(k2-4alpha0)eta.

For a root-star cap f*eta, f=3,5/2,2, and d_u<=eta, subtract
g_R^2 to obtain its lower diagonal. The other controls are

```
|g_R|<=2eta, |g_Q|<=sqrt3*3.67eta,
K_QQ>=1-6k2*eta-(sqrt3*3.67eta)^2, K_QQ<=1,
|K_RQ|<=sqrt3*11.05eta+|g_R*g_Q|.
```

The least eigenvalue is>.94 even at f=3. The determinant checks
give lower density normalizers (6136901,6150156,6163383)/10^8
for f=3,5/2,2. The common upper quadratic has diagonals
(1003491987,1000000001)/10^9 and cross magnitude21835167/10^9.
There is no independent-driver substitution.

For a signed score s_v(Y_v-t), s_v=+ for inside and - for guards,
the exact conditional mean is

    s_v[(c_v-1)t+c_v*x+(Cov(Y_v,(R,Q))-c_v*g^T)(y,z)^T].

Here c_v=Cov(Y_v,O). In the R column retain the joint coefficients
of d_u,d_a,d_b before maximizing their nonnegative light caps.
For c_-<=c_v<=c_+, the six resulting error coefficients are

```
a,b: (2J-2c_-)_+ +(2J-alpha0+c_+)_+ +(c_+-alpha0)_+;
a1,a2: 17.8+2c_+;
w: 2-2alpha0^2+2(c_+-alpha0^2)_+;
w1,w2: 4J-4alpha0^3+2(c_+-alpha0^3)_+;
u: 2J+2c_++2(J-c_-)_+;
w': 4(k3/3)-2+2c_++2(k3/3-c_-)_+.
```

They follow by substituting Cov(R,O)=2d_u-d_a-d_b and charging each
raw covariance on its actual upper or reverse path. The replay takes
the smaller of these and the older valid upper bounds. For Q,
family AppendixA pays the six raw-plus-root costs
(7.02,13.38,6.92,11.19,3.66,6.92)eta before sqrt3, with the additional
3.67*k_d*eta^2 root-covariance correction. The ideal y coefficients
are (-1/6,-sqrt2/6,sqrt2/3,1/2,-1/3,-sqrt2/3), and the ideal z
coefficients are (0,0,-1/(2sqrt3),-sqrt2/(2sqrt3),0,-1/(2sqrt3)).

The checked upper (x,y,z) coefficients, in the same score-group order,
are the first three columns below divided by10^5:

```
 94395 -15833   1387   315
 83831 -21349   2646   917
 83831  47332 -27497   917
 71990  50917 -38606  1611
-94280 -32338    724   315
-83333 -44809 -27497   917
```

The last column bounds (1-w_d)*.055, at scale10^5. At t=sign*T,
an adverse threshold term scales that column by T/.055; a favorable
term retains -(1-c_v^upper)T. Every score covariance is<1.
Strict upward rounding gives six available margins f,g,m,e,u,v.

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

For the private variables use the entire centered convex regions

```
ordinary: |D1|<=f, sqrt2|D1|+|W|+|A_a|/2<=g;
shifted:  |P|+|A_a|/(2sqrt2)<=f, sqrt2|P|+|W|<=g,
          P=D1+sqrt2*A_a/4.
```

First reduce f,g by r/3,sqrt2*r/3 when the root radius r is debited.
The exact reconstruction above proves both sufficient. Use the better
probability bound, never their sum. Enlarge the RAW three-variable
covariance diagonals to SDs a,p,b, keeping its off-diagonal entries.
Here b=.289; p=.289 ordinary or .290 shifted; and a is respectively
.024224,.022113,.019779 for private caps3eta,5eta/2,2eta. These
enlargements also dominate the conditional residual covariance.
Anderson's inequality for centered convex symmetric sets applies.

Put S_a=cap*eta. The raw source gives
Cov(A_a,D1) in[-S_a/2,S_a/4], |Cov(A_a,W)|<=S_a/(2sqrt2),
and |Cov(D1,W)|<=(1+k3)eta/4. The shift adds the nonnegative
sqrt2*Var(A_a)/4 to Cov(A_a,D1), retaining its absolute boundS_a/2;
it adds at mostS_a/8 to the other half-pair covariance. Also
Var(P)<=1/12+(sqrt2+c)S_a/8. All enlarged diagonal guards are
checked strictly. Let x,y,z bound the three absolute correlations.
Their pairwise row sums are<1. The inverse-precision diagonal bounds
are the reciprocals of

```
a^2[1-(x^2+y^2)/(1-z)],
p^2[1-(x^2+z^2)/(1-y)],
b^2[1-(y^2+z^2)/(1-x)].
```

These follow by lower-bounding each two-coordinate correlation block
by its least eigenvalue1-|correlation|. The last two variances exceed
1543/20000. Fold all eight density signs. Jensen cancels every
off-diagonal exponent term, and Hadamard supplies prefactor
8/[(2pi)^(3/2)apb]. In h=|A_a|/2 coordinates, integrate its marginal
lower density in bins with endpoints
`.0005*(0,1,2,4,7,11,16,22,30,40,52,64)`.
Use the exact bin-mean square in its exponential. Its right endpoint
pays closing debit h in the ordinary recipe, or root debit h/sqrt2
in the shifted recipe. The remaining trapezoid has area and second
moments integrated exactly with `exp(-q)>=1-q`, using variance1543/20000.
Its pair prefactor is rounded down to38/5 or379/50; the equivalent
rescalings of84/11 are209/210 and4169/4200. Negative lower areas pay0.
This integrates each WHOLE convex private event, not Gaussian
correlation on nonconvex individual bins.

At the tight private cap replace each integrated recipe by the smaller
of its tight lower mass and(50/47) times its ordinary lower mass.
Consequently47/50 of that contribution remains at the ordinary cap.
Unrestricted private bins are compared to the original tight bins;
their minimum ratio applies after clipping as well. Both recipes use
the same bin edges. They keep one physical event and one name weight.

For the terminal pair D={|M|+|H|<=e}, let the three role caps be
z_o,z_u,z_w. An upper SD for M is
`a=(sqrt(c*z_o)+sqrt2*sqrt(c*z_u)+sqrt(c*z_w))/2`, and use b=.289
for H. Raw row ratios give
`|Cov(M,H)|<=sqrt2*(z_o+z_u+z_w/2)/4`.
Enlarge diagonals to these directed SDs, retaining the off-diagonal.
If s=rho^2, folding the four density signs and dropping cosh gives

    2/[pi*a*b*sqrt(1-s)]
       *exp[-(x^2/a^2+y^2/b^2)/(2(1-s))].

Thus its diamond mass is at least
`sqrt(1-s)*P0(e/sqrt(1-s))`, with P0 the INDEPENDENT enlarged
diamond mass. Divide the allowed s interval into eight panels;
round its prefactor downward at the right endpoint and its radius
downward using the left endpoint. Take the minimum over panels.

The independent diamond and the core below use folded-root bins of
width.001 and the mean-square Jensen density mass. In each retained
bin both standardized strip endpoints are nonnegative. The remaining
strip envelope g is decreasing and concave, hence above its secant;
the folded density f is decreasing. Their nonnegative covariance
under uniform measure proves
`integral_bin f*g >=(integral_bin f)*(g(left)+g(right))/2`.
Round both endpoint strips downward. Omit bins crossing the strip
feasibility boundary. This secant payment and the retained terminal
determinant normalization are part of the current replay.

### Adaptive core and reflected asymmetric margins

For root radius r use the symmetric convex core
R(r,m)={|X|<=r, |D_u+(2sqrt2/3)X|<=2m}.
Projection on the root gives conditional variances at most
c*z0-z0^2 and c*z0-z0^2/2 for X,D_u, respectively. Enlarge the
diagonals to checked SD squares sx^2,sy^2, retaining the covariance.
This dominates the conditional centered covariance. To cover its unknown
correlation, divide[-1,1] into48 panels. In a panel use the maximal
absolute mean multiplier |2sqrt2/3+rho*sy/sx| at its endpoints and the
largest conditional SD sy*sqrt(1-rho^2), at the point nearest0.
On folded-root bins of width.001 use that panel mean multiplier and
Jensen's mean-square marginal mass, and the secant-bin bound above. Only bins with
2m larger than that mean contribute, so both normal CDF arguments
are positive and increasing either chosen upper loss lowers the mass.
Take the minimum over panels. CDF arguments are rounded downward to their integer strip indices. Singular endpoints follow by
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

At the two endpoints -T,+T the strongest grade (root,private,middle,
terminal caps)=(2,2,5/2,5/2)eta pays316980084/10^14 and342931255/10^14.
The meshes have y,z index widths(8,4),(30,16), corresponding to physical
widths(.02,.04),(.075,.16), and use813458,15113 disjoint boxes. Root
edges are the fixed listed partitions of[0,.3], in units1/2000. All
negative-margin boxes and omissions pay0. No zero-threshold replay
or original-cut quadratic deficit is spent in this composition.

For a FIXED physical name and its grade, the original ten signed
score inequalities are affine jointly in the Gaussian field and t.
Their cylinder is convex. Gaussian integration (Prékopa on the range
of the covariance for a singular limit) makes its probability log-concave
in t. The smaller endpoint floor therefore pays all t in[-T,T].
Different endpoint cells/recipes are allowed because they bound the
SAME physical cylinder, not different success events.

There are36 grades: root/private caps each2,5/2,3 timeseta and
middle/terminal caps each5/2,3 timeseta. For each take the product of
(1) the minimum positive-entry core-table ratio at its root/middle caps,
(2) the analogous midpoint-table ratio at its three caps,
(3) private ratio1,47/50,or the minimum unrestricted/tight bin ratio,
and(4) its root density-normalizer ratio. These positive comparisons
hold for every same corner recipe, pass through maxima and geometric
corner products, and multiply the strongest endpoint floor. Core or
terminal reference entries zero contribute no mass. The checker
verifies that each of the36 scaled common floors strictly exceeds
the selected signed quota minus its exact role charges in Section13.
The minimum grade slack is positive. Eligibility and grades are fixed
by original marks, so this comparison is uniform across thresholds.

### Direct precision means for strict and relaxed nine names

Use the thirteen-edge faithful stencil of Section6. Strict names have
all edges<=theta; relaxed names allow eta only on the four terminal
child edges. Keep the same guarded sign cylinder. Put O=t+x and raw
drivers D=(Y_b+Y_u-2Y_a, sqrt3(Y_w'-Y_w),sqrt3(Y_z'-Y_z)).
For g=Cov(D,O), K=Cov(D|O), use D-gO=K*xi,xi=(y,z,h). Its density is

    sqrt(det K)/(2pi)^2 exp(-[(t+x)^2+xi^T K xi]/2).

The root, sibling, mixed-driver and opposite-ending inequalities in
`family.md` Section4 give lower diagonals
`(1-(8k2-2sqrt2)*theta-g0^2,1-6k2*theta-g1^2,1-6k2*theta-g2^2)` and off-diagonal
magnitudes `(sqrt3*24.55theta+|g0g1|,sqrt3*24.55theta+|g0g2|,
35.3theta+|g1g2|)`, where
`|g|<=(2theta,sqrt3*3.67theta,sqrt3*3.67theta)`.
The comparison matrix is positive definite by strict diagonal dominance.
With lower diagonals r,q,q and magnitudes c,c,f, determinant monotonicity
under adding a positive diagonal gives

    det K >= r*q^2-r*f^2-2q*c^2-2c^2*f.

The checker pays normalizer123513/5000000. Its lower comparison
matrix has least eigenvalue>.93. The upper quadratic retains actual
upper diagonals `(1+2(k2-2sqrt2)*theta,1,1)` and these cross magnitudes.

For EACH score use its direct conditional mean

    E[Y_v|O,D]=c_v O+(Cov(Y_v,D)-c_v*g^T)xi.

Sign this expression for below/above guards. For a signed raw covariance
sum, assign each positive distance-d excess the upper per-edge factor
k_d/d and each negative one the reverse factor alpha0^(d-1).
Accumulate coefficients on the ACTUAL path incidences before maximizing
their nonnegative caps. Strict edges have captheta; relaxed terminal
edges have cap2theta. The opposite-branch contrasts also have separate
upper costs7.85,14.8,26 times theta by depth, with28.9 at relaxed endings.
Keep the smaller valid bounds. Regression costs c_v^upper|g_i|.

In either Q column improve the full raw-plus-root expression by
family AppendixA on each matching seven-key subtree. The sole doubled
edge in the reused selected-ending-child row costs1.88theta extra;
its root-covariance error is multiplied by4/3 in the quadratic term.
The far opposite terminal scores keep26theta or28.9theta and the old
regression bound. These comparisons use the SAME physical rows as K.
For t=sign*x on a half-band[a,b], an adverse signed threshold term
uses b(1-c_lower); a favorable term retains -a(1-c_upper).
Take the minimum available margin in each score group.

The centered score identities and strips remain the original ones:
six wave-defect strips, four half-differences, and an optional centered
root/ending trapezoid. The two Q-driver endpoints have the SAME centered
loss after conditioning. An ending/root half-pair has covariance below
the recomputed bound and diagonal<=1/12, so the checked enlargement
.09 I dominates it; integrate one selected trapezoid with slope1.5.
Other halves/defects use centered strip probabilities and Gaussian
correlation. It is invalid to add both alternative ending-pair masses.
For the correlated root/ending half-pair the raw variances are<=1/12.
Let p_j=k_j strict, or (2,3k2/2,4k3/3,5k4/4) relaxed. Its absolute
conditional covariance is at most
`(p2+p4)theta/4+p3*theta^2/4+26*113*theta^2/.93 <13/2000`.
The26,113 bound the raw regression-row norms and.93 the driver inverse
gap. Thus .09 I dominates the pair covariance. A single selected
root/ending trapezoid with slope1.5 is integrated by its area and
second moments, using exp(-q)>=1-q and prefactor700/99. Other centered
strips are paid jointly by Gaussian correlation. All calibrations
are regenerated at eta=.00114. The core/terminal SD uppers are
.017130,.022115, from3c*theta and5c*theta. At radii
(.018,.02504,.03365,.028), their respective strip lower masses are
(706639,856187,950508,897852)/10^6 and
(584309,742470,871882,794520)/10^6; the Gaussian prefactor19947/25000
and the integrated degree31 Taylor lower polynomial certify them.

The code recomputes every supplied nonzero choice. Box density uses the
upper quadratic averaged on that box, then Jensen and downward integer
exponentiation. Reflection z<->h includes both branch labelings; the
envelope applies to each separately and does not assume exchangeability
of the actual covariance. The relevant reflection multiplicity and
root-coordinate volume are retained. Root bins have widths.005 on
[0,.1] then.01 on[.1,.2]; breadth bins have width.05 on[0,3], and
both ending coordinates width.1 on[0,4], with820 unordered pairs.
For the root density take the larger box-mean square at BOTH band
endpoints; convexity covers all intervening thresholds. The other
driver density is shared across bands. Every choice code is checked
on the actual margins; invalid or omitted choices pay0. No search
score is accepted as mass. The resulting floors at scale10^12 are

```
band          strict     relaxed
[-T,0]       1790351     1389839
[0,T]        1974557     1537579.
```

These are compared with u+v9+v11 and u in the TWO allocations of
Section13. Neither the old .0009 floors nor its clipped choices are
used. Different allocations do not change the nested keys or masks.

Nonsingular driver blocks, compact covariance controls and strict floor
gaps permit the common finite transfer for all names and |t|<=T.

## 11. Capped Price energy and exceptional components

Use eta=57/50000,theta=eta/2,T=277/10000,D'=min(D,eta),
lambda=141421/50000. Section9's global/local replacements give
K0=1297660793724599979211/25000000000000000000000 and

```
group    raw cap                  multiplier   dirty b1,b2,b3
adj      29379553/5000000000       1            29,29,29
D2       23969081/10^10            2lambda      0,18lambda,18lambda
D3       372789/10^9               19           0,327,327
D4       502461/10^10              18lambda     0,522lambda,522lambda
D5       19651/2500000000          120          0,0,6040
D6       3019/2000000000           92lambda     0,0,7788lambda
stratD2  1999423/5000000000          2lambda      0,18lambda,42lambda.
```

Spend x_g=cap_g*f_g with
f=(6445176870,10^10,10^10,3854668533,10^10,8512011682,10^10)/10^10.
The same spends enter credits and debits. Local replaces global;
size4/6 excludes baseline/size2. Adjacent71931/10^10 is removed
globally before entering its clean replacement.

Set H=(.0031,.0069,.00815),Ej={D>Hj}. On a common high-girth cover
expand the branch identities U1=lambda0*b_u-b_w+R_u and
U2=6b_u-lambda0*b_w+lambda0*R_u+sum_(z~u,z!=w)R_z,
lambda0=2sqrt2. Retain every nonnegative residual bracket, using
b_u>=alpha*b_w. D2 retains lambda*S_w with6 sites; D3 retains
9D_ab,5S_a with13 edges,14 sites; D4 retains6lambda*S_w,3lambda*S_u
with14,30 sites; D5 retains44D_ab,24S_u,14S_z with29 edges,30,62
sites; D6 retains24lambda*S_w,14lambda*S_u,8lambda*S_z with30,62,126
sites. A q*S site costs3q*eta. Endpoint incidence sums give exactly
the common/dirty table. StratD2 has6 sites at H2,14 at H3.
Projection and division by cover degree preserve these upper counts.

Incident exceptional edges introduce at most8 adjacent sites,
2 for D2,(4,4)D3,(4,8)D4,(8,8,16)D5,(8,16,32)D6. Tree unions
refund twice S,L,t5,t6 per added edge, where

```
S=(21/2)eta*x1, L=eta[lambda(6x2+189x4)+(231/2)x3],
t5=2220eta*x5,t6=2886lambda*eta*x6,r7=lambda*eta*x7.
```

Disjoint stratD2 pairs refund12r7,30r7,18r7 for middle/middle,
top/top,middle/top; unpaired edges remain single. D6 refunds on pure
top components and top blocks of mixed components have distinct owners.
Only discarded K0(D-D') is collateral. Select isolated E2 edges,
isolated lower E1 edges and maximum-D edges in pure higher dyads:
they form a matching. The paired source law gives

```
Ij=K0(sigma*Hj-4eta)/2,
Jjk=(K0/2)[sigma*Hj-kH1-(4-k)eta]_+,
Jj*=K0[(sigma-1)Hj-3eta]/2, sigma=2021/1250.
```

Each lower neighbor transfers K0(D-H1)/2 at most twice; outside E1
there is at most one access per endpoint. Unshared eligible owners
keep the other half. Each shared connection between components adds13
overlapping adjacent names. A forest rooted at ineligible components,
or a tree on p>=2 eligible components, gives>=13p/2 names. Therefore
min(I1,(13/2)eta*x1)=I1 survives and higher singletons/dyads keep
b0=(13/4)eta*x1<=min(I2,I3,J2*,J3*). This forest is separate from
the within-component refunds. Put m=10^-12 and

```
Rj=K0(Hj-eta)-eta sum_g b_jg*xg+min(S,Ij),
p0=R1+S-I1-m,P=p0+S,tj=m-Rj-L (j=2,3),
v=max(t2-6r7,(t2+t3)/2-9r7,t3-t5-t6-15r7).
```

All clamp signs are checked (only J24 has negative raw argument).
Here p0>0: lower components use2p0,3p0+S,P; isolated edges use
R1+I1>=m. Higher singleton/dyad tests are
Ij+b0-S-L>=tj, J2*/2+6r7+b0/2>=t2,
J3*+18r7+b0>=t2+t3, t5+t6+J3*/2+15r7+b0/2>=t3.
Larger pure cases use S+L-3t2+12r7>=0,
S+L-t2-2t3+18r7>=0,S+L-3t3+4(t5+t6)+30r7>=0;
additional edges use max(t2,t3)<=S+L. Mixed nonsingleton higher
blocks retain2(S-v).

For l lower edges,q_jk higher singleton blocks,r nonsingleton blocks,
sum k*q_jk+r<=2l. The remaining pool is
F=sum a_jk*q_jk+2(S-v)r+Pl-2S,a_jk=S+Jjk-L-tj.
At l=1 the five patterns(1,0),(2,0),(0,1),(1,1),(0,2) are paid by
p0+Jj1-tj-L,P/2+Jj1-tj-L,P/2-v,
p0+2S-L+Jj1-tj-2v,(p0+3S)/4-v.
At l=2 use min(a_jk,a_jk+2p0)>=0 for k=1,2 and
min(2(S-v),2(S-v)+2p0)>=0. For l>=3,beta=p0/2+S(l-2)/(2l),
F=sum mu_jk*q_jk+[2(S-v)+beta]r+beta[2l-sum k*q_jk-r],
mu_jk=a_jk+k*beta.
Lower edges are in no triangle/square: their unicyclic-cover excess
floors(2+2sqrt6)/7-2sqrt2/3>.04 and
(2sqrt170+2sqrt80)/45-2sqrt2/3>.03 exceed H2; coarsening increases
affinity. Thus singleton blocks share<=1 lower link. At l=3 there
is<=1 degree3 block; at l=4,5 <=1 degree4 (all others<=2 for l=4);
at l=6,7 <=1,2 degree4 blocks. Two degree4 blocks need7 lower links,
three need9. For each stratum with z exceptional degree-d blocks,
the checker verifies that(2l-dz)min(beta,2(S-v)+beta,mu_jk/k for
other allowed degrees)>z max(0,-mu_2d,-mu_3d). At l>=8 use
beta=p0/2+3S/8: all coefficients and Sl/4-2S are nonnegative.
This exhausts components, reserving m per exceptional edge.

Hence K sum D'-o(n) is paid, with
K=K0+sum multiplier_g*xg+65797/200000000 and
A=eta*K=11569604236818386923776429/125000000000000000000000000000.
The D7--D10 increment has no dirty debit. Each floor, half-charge,
pair/forest refund and capped endpoint is allocated once.

## 12. Marked endpoint, run and center payments

Divide all source payments in this section by A. Put

```
a0=sigma-3/2, c=a0/4, h=(sigma-1)/4, k=33/100,
rho=2X/A, s=Z/A, g=(2/3)(s-c-rho/4), B=c+g/2,
d=(sigma-1)/8, nu=(1/4+rho)a0, nu2=s-1/8,
g0(l)=min(l/2+rho min(l,1/4),s), f(l)=l/2+g0(l),
kappa=59599/10^6, omega=9/20.
```

The exact domain is2/5<rho<1/2 and1/8+rho/4<s<1/4+rho/4.
In particular g0(l)>=min(.9l,.225)>=omega*l for0<=l<=1/2.
Let C,p,b,t,z,epsilon be respectively the normalized type11 payment,
high-anchor debit, type10 coupon, type01/type02 payments and type20 bonus.
They are specified in Section13. Write delta=(C-B)_+, p*=p+(C-p)_+.

The unchanged endpoint/stub proof in family Section6 assigns high ends1/2,
medium x/2, shared low x/2+s and boundary low x+rho min(x,1/4).
It pays XU2+(3X/2)U3+(X-A/5)J and one2Z per shared low. Put
H1={10,11}, R_e=sum_off x-sigma*x_e. Its high-stub identity leaves
V0=.1916,V1=.0374 at20 with zero/one H1 neighbors. H1--H1 residual
halves are used only after their fixed2h bonuses. The ordinary table is

```
01 1/2+rho/4+B   02 1/2-g       03 3/4-3B
10 3/4+h+s       11 3/4+h-B     12 1-h-2B
20 sigma/2       21 5/4-2h-B    30 3/2-3h.
```

Use family Section7's phi(t)=sqrt(t)-sqrt2(1+t)/3,G(t),Cstar and
bilateral k=.33, with both root atoms. Its bound
G(t)>=(16/15)phi(t)(t-1)^2/t and concavity of
G(t)-.0596phi(t)-t^2 give the weighted-star chord
min(F(lo),F(hi))-(hi-lo)^2/4. Widest-axis bisection proves
G(p)+G(q)+.8phi(r)>=.0596(phi(p)+phi(q)),p+q+r=3 in[1/2,2]:
1443 leaves pass,59 are outside. The smaller kappa pays source errors.

Let Pm,Rm count01--01 and other01 medium links,M=3n03+n11+2n12+n21,
and Nj count02 centers with j01 leaves. Beyond the old medium table
the exact count leaves(g/2)Pm+gN0+(g/2)N1+(kappa/2)N2+B(M-Rm).
For N2 the pool after both c fees is R/2+f(l)-l+(sigma-1)sum s_i/2;
the star gives R+.8l>=kappa. Small boundary ports use f(l)-l>=.4l;
larger/shared ports pay more. These stars have no11 medium neighbor.
All actual fees and the full B remainder survive.

### Runs and their disjoint terminal pools

Remove02,03,12,20,21,30 anchors. The11 components alternate high/medium;
a run with N11 vertices,e01 ends,m medium ends has count B(N-e) and
N-e-m internal medium incidences. Opposite medium halves, residual
halves after bonuses and quarter-lows form its old pool AFTER actual
fees. The10 fee reserve is separate; a20 arm has no H1--H1 residual.

A good H1 connector of excess m uses R_i+R_j+L>=km, hence km/4
from quarter-residual/quarter-low shares. Designate an adjacent all11
pair when present. Its other medium<=1 gives L'>=(sigma-1-m)_+;
after quarter-low draws it keeps>=l/5. Pay another kL'/4: a pair
is designated<=2 times and k/2<=1/5. Thus each good incidence pays
d*=k(sigma-1)/8; require delta<=d*. The only good connector between
two10--11 pairs (N=2) uses the two10 delta reserves instead. Exclude
all connectors meeting separately selected terminal pairs.

For all11 terminal pairs let a=sigma-1/2. Their paired law is
S+L-R>=a-m and connector law R_i+R_j+L_c>=km. For small lows,
g0(l)>=.9l. Weights(.5 connector,.4 outer law) give end group
G=.4a0+k/2=.21172. For closed N=4 use(.7 connector,.2 each outer)
to get F4=.4a0+.7k=.27772; for N=6 use(.25 each connector,.5 each
outer) to get F6=a0+k/2=.2818. Residual,low,medium coefficients are
<=.5,.9,.5 respectively. A large low gives.225; add the opposite
paired a0/2 at N=4 or the unaffected end group at N=6. Two large
lows give.45. At N>=8 disjoint end groups pay2G, excluding<=8
medium incidences; >=(N-10)_+ good incidences remain. The missing
tail beyond N=8 costs<=2delta, giving closed_N8 below. Medium
anchors lose at most c per actual fee.

On01--11==10, paired law beta+l>=a0+1-m gives whole available
pool>=g0(min(a0+1-m,1/2))+.75(m-.5)>=s by concavity; this includes
the10 fee, external half and01 own low surplus once. For N=3, let
y=beta+l at10. Pooling both high residual halves and far low/medium
shares gives g0(y)+omega max(a-m,km-y). For y<1/4 its minimum is
at y=a-m,m=1, giving omega(k+a0)=.20106; y>=1/4 pays more.
At N>=5 retain nu at10 and far reserve r*=ka/2; >=N-5 good
incidences pay the tail except2delta. For closed N=2 full lows and
opposite halves pay omega(sigma-1); the old restriction C<=d is absent.

With no01 end, cycles have every incidence good; high anchors lose<=2
incidences, charged2delta. A10 end keeps min(nu,nu2)-b>=delta.
One medium end of length3 retains c from its paired low/opposite
half after fees; from length5 retain(133/180)G-c>5d*, excluding
terminal pairs and the high-end connector (>=N-7 good incidences).
Two medium ends use2d-2c>2d* at N=2 or G-2c>6d* with>=N-6 good
incidences. One01/medium uses N=2,F4,F6 and then2G-c at N>=8;
2G-F6>4d* pays its entire tail. One01/high uses the short charge
C+t; N=3 retains omega*a0>2d*, and N>=5 retains G-4delta
beyond length3, greater than ka/[4(1+k)]. None pays a10 fee.

A12 high anchor keeps debit p: an unused medium stub gives B, or,
if both face01, its two half-surpluses, the near11 medium half
(opposite half at N=1) and near full low pay omega*a0=.05256 by
the paired law. The first connector is excluded;01 targets are separate.
Thus(C-p)_+<min(B,.05256) suffices. Types21,30 pay p*=max(p,C)
per high stub. Require p>=2delta and b<p. A10--10 pair separately
requires2b<s+(1/2+rho)a0, its post-bonus concavity pool.

### Type20 source, outside energy and recovered connectors

We need the following stronger dyad law on a-b-c:
sum_(four outward incidences at a,c)D+D_g/10>=(13/20)(D_ab+D_bc)-o(1).
Here g is the third edge at b, and a closing edge ac is counted twice.
For a scalar proof define A(t)=sqrt(t(5t/2-1))-sqrt2(7t/2-1)/3.
A and phi have the common rationalized numerator2(t-1/2)(2-t)/9.
Thus A>=v phi follows from
sqrt(t)[1-v sqrt(5t/2-1)]+(sqrt2/3)[(1+v)+(1-7v/2)t]>=0.
The directed upper pairs(t,sqrt(5t/2-1)) for
v=(13/20,7/10,3/5,3/4,11/20,53/100,1/2) are respectively
(21/16,121/80),(19/16,45/32),(23/16,13/8),(17/16,31/24),
(5/3,9/5),(7/4,15/8),(2,2). Substitution proves the inequality;
for the first two rows split at1, use sqrt(5t/2-1)<=5/4 below1,
and sqrt(t)>=1,sqrt2/3<1/2 above1. Their positive margins are
33/6400,3/640; the other rows have nonnegative brackets.
For z=5/2-t, rationalized numerators agree and increasing denominators
give phi(z)>=phi(t) for t>=5/4, >=(7/5)phi(t) at t>=5/3,
and >=(3/2)phi(t) at t>=7/4. Assume p>=q,p+q+r=3.
If p<=21/16 apply the first row twice; if p<=23/16 use the3/5,7/10
rows plus shifted concavity phi(q)+phi(r)>=phi(5/2-p)>=phi(p);
otherwise q<=17/16 and the last three rows together with the
3/4 row prove A(p)+A(q)+phi(r)/10>=(13/20)(phi(p)+phi(q)).
The fiber formula and uniform atom/clamping argument of Section7
give the claimed physical dyad law, including repeated images.

At a two-H1 center put high excesses1+s_i, low l, leafoff-sums
1/2+y_i,Y=sum y_i. After fixed fees its pool P0=Y/2+f(l)-a0>=kappa,
with P0=(R1+R2)/2+f(l)-l+(sigma-1)sum s_i/2. The dyad gives also
P0>=.0332+1.35l on small boundary lows; larger boundary/shared lows
pay.35-a0 or s+.0332.

Dyad and paired laws give h_i<=108750/42523+o(1): use leaf
off-sums A_i<=3/2,A1+A2+l/10>=.65(h1+h2),h_j+l+A_i>=sigma*h_i.
Thus eta*h_i<H1=.0031, and every arm/off-edge is nonexceptional.
Its discarded global K0 high surplus has no Section11 owner.
The exact5/8<K0/K<1 in particular pays rate omega.

A shallow leaf (10 or length-one11--01) retains omega*y_i from
unused low/opposite-medium shares after old baselines. Two shallow
leaves have Y>=.3-l/10>=1/4, hence S0=kappa+9/80. One shallow leaf
and the free arm surpluses pay omega(y_i+sum s_i)>=omega(a0-l).
Combine P0's two small-low bounds to cancel l, obtaining
M0=(2kappa+.0332)/3+omega*a0; other ports pay more.

On a nonfar11 arm of length>=3 its first connector was excluded.
Its draw X_i=zeta(R_i'+l_i+l_i'),zeta=47/400, satisfies
X_i>=zeta*l_i and X_i>=zeta(km-R_i). Residuals fit since zeta<=1/8;
lows fit since zeta+k/4=1/5. Two recoveries use distinct lows and
<=2zeta of a pair. A medium-terminal group retains133/180 of each
resource: low fractions<=zeta/omega, residual fractions<=2zeta.

With no far01 ends and j=1,2 long arms, P0+sum X_i has capacity
N0+jN1, N0=(1-2zeta)kappa,N1=zeta*k/2, and costs(2+j)delta.
With one shallow leaf keep m+l_i=1/2+y_i in the first connector:
X_i>=tau[k(1/2+y_i)-R_i], tau=zeta/(1+k). Since
R_i<=2P0-2rho*l, the combined pool is at least
[1-2tau(1-k)]P0+[tau((2-2k)rho-k)-omega]l
+(omega+tau*k)a0+tau*k/2. Use rho>=2/5 and put weight
w=[omega-(4/5-9k/5)tau]/(27/20) on.0332+1.35l and the remaining
positive weight on kappa. This gives
M2=82335849341/718200000000. Shorter nonfar arms retain M0.

A far01 arm with N>=5 retains F=G-4delta at its far end. For
j=1,2 length-three far arms, combine P0>=sum(m_i-.5)/2+f(l)-a0
with P0>=kappa in weights.9,.1; add omega(sigma-.5-m_i) per far
pair. Cancellation gives K_j+.9f(l),K_j=kappa/10+omega*j*(sigma-1)-.9a0.
Adding the shallow omega(a0-l) gives
MF=min(K_1+omega*a0,M0+F). Near/far resources have distinct owners.

A one-H1 center receives O=7/200 from its other high edge.
A saturated neighbor has an unused stub h>O; a zero-H1 neighbor
reserves2O. Two one-H1 centers pair on their unique other high edge,
a matching. If combined leaf surplus U and center lows L give
U+L>=2a0, their half shares pay both bonuses. Otherwise the dyad and
paired laws imply,chi=13/20,
(2chi-2(1-chi)/sigma)h<=(1+(1-chi)/sigma)(1+U)
+(11/10+(1-chi)/sigma)L. Substituting U+L<2a0 gives h<2 and leaf
arms<2, hence no exceptional accesses. Each end's half-surplus/low
plus free high surplus pays(K0/K)a0/2>O, with no two-H1 owner.

The complete normalized menu below is checked strictly in both bands,
with delta<=d*,(C-p)_+<min(B,.05256),b<p and2b<s+(1/2+rho)a0:

```
t<=g/4; t+z<=g/2; 2t+z<=kappa/2;
p>=2delta; C+t>=2delta; min(nu,nu2)>=b+delta;
s>=C+b+t;
2B+omega(k+a0)>=b+3C+t;
2B+nu+ka/2-2delta>=b+3C+t;
2d+B-c>=2C+t;
3B+F4-c>=4C+t; 3B+F6-c>=4C+t+2delta;
omega(sigma-1)>=C+t;
B+F4/2>=2C+t; 2B+F6/2>=3C+t;
3B+G>=4C+t+delta;
V0-2O+2p>=epsilon;
V1+O+2p>=epsilon+b, epsilon+C+t;
S0+2p>=epsilon+2b, epsilon+b+C+t, epsilon+2C+2t;
min(M0-delta,M2-2delta)+2p>=epsilon+C+t, epsilon+b;
MF+2p+2B>=epsilon+4C+2t, epsilon+b+3C+t;
min(K_2,K_1+F,kappa+2F)+2p+4B>=epsilon+6C+2t;
min(K_1,kappa+F)+2p+2B>=epsilon+3C+t+2delta;
min(kappa-2delta,N0+N1-3delta,N0+2N1-4delta)+2p>=epsilon.
```

The center guards exhaust zero/one H1,two shallow,shallow/nonfar,
shallow/far,10/far,two far,one far,two nonfar. Nonfar lengths1,2 cost
delta per end; longer arms cost2delta with first-connector recovery.
Two ends on one run use its B count once; the good-incidence count
and terminal exclusions pay all additional length.

### Type02 role charges and final row shifts

After t,z, the medium remainder is(g/2-2t)Pm+(g-z)N0+
(g/2-t-z)N1+(kappa/2-2t-z)N2+B(M-Rm). Keep B for runs/type12.
A02 role charge D02/A occurs only at incident sum>2eta. Write its
mediums1/2+s_i; then s1+s2+l>1. With two01 leaves the old post-fee
pool is>=1/4-2c=.1916, paying D02 above kappa/2. With one01 leaf,
s1+l>1/2 and the unused own half s1/2+f(l)>=1/4 pays c+D02/A;
the opposite half remains for11. With no01 leaf use g. Exact guards:
g-z-D02/A>0,1/4-c-D02/A>0,.1916-kappa/2-D02/A>0.

The mixed01--02--11 actual remainder G_actual-c remains its run's.
All other01 targets are run ends or saturated medium neighbors.
Charge t per such medium incidence or high-anchored01 run; a saturated
vertex has<=3 incidences. Relative to the ordinary table the shifts are

```
01:+At, 02:+Az, 03:-3At-D03, 10:+Ab, 11:+AC,
20:A(epsilon-2p), 12:-A(p+3t),
21:-A(2p*+3t), 30:-A(3p*+3t).
```

D02 is paid internally, with no second row debit. Each fixed bonus,
actual fee, medium half, low share, connector/end group, center and
outside-cap resource has one specified owner. Errors have bounded
multiplicity even on long runs and sum to o(n).

## 13. Role grades, two signs and the exact marked ledger

Section10 supplies36 grade ratios R_grade against(2,2,5/2,5/2),
using caps(2,5/2,3)eta at o,a and(5/2,3)eta at u,w. Its physical
cylinder and endpoint log-concavity give Q*R_grade on[-T,T],
Q=316980084/10^14. In units10^-14, intermediate role charges are
(40842978,19018805,0,0), high charges(71228409,58127297,82536630,9109886).
Each grade passes Q*R_grade+sum charges>q. Only03 needs3eta; only02
with incident sum>2eta needs an intermediate root/private cap.
Each vertex has<=6 weighted names per role:3 incident choices and2
child choices, with private/partner choices of total weight<=1.
This retains weight1 for a high uninspected partner,1/2 otherwise.
Thus D03=1326013332/10^14,D02=359170698/10^14; Section12 pays D02
internally. No actual-base deficit or curvature credit is spent.

For high-marked H,let Ui count unmarked vertices with i high neighbors
and Eij their edge classes. The weighted seven identity is
L7=sum_(UU edges)(2-i(v))(2-i(w))+E01+E11
=6U-4rH+E12+4E22+6U3>=6U-4rH.
Thus L7>=6n-sum c7(a)n_ab, c7=(0,14,10,6).

The unchanged branch enumeration of family Section5 applies at theta.
Its strict/relaxed populations yield L9,L11,LR,LA with constants9,11,12,9
and the following marked/overlap coefficients (r=sum(3-a-b)n_ab-2sS,
U0=U-r+U2+2U3,rM>=2n01+n02-2sS):

```
j=a+b, c9=(0,35,22,9), c11=(0,45,28,11),
cR(a,b)=(0,36,24,12)[j]+8(3-j)1[a>0],
cA(a,b)=c9[j]-2(3-b)1[a=0<b],
           sS U2 U3 J
L9         26 12 30 4
L11        34 16 40 6
LR         24 12 24 0
LA         22 12 24 2.
```

L9,L11 bound the same strict count,LR,LA the same relaxed count;
negative affine bounds are allowed, with no duplicate event payment.

Use these integers divided by10^14;Cabs=AC,pabs=Ap,babs=Ab,E=Aepsilon,
e=u-gap:

```
band      q         u         v9       v11      gap X           E
[-T,0]    316979984 138983800   4781365 35269835 100 2067212632 832074027
[0,T]     316979984 149025542  39628377      100 100 1994982270 872240994
band      Cabs      pabs      babs
[-T,0]    784951776 497089265 497089165
[0,T]     784951776 517172748 517172648
At=Az=100/10^14 in both bands; weak weight=0.
Z=13v9+17v11+11u+gap, rho=2X/A, s=Z/A.
```

The combination is qL7+v9L9+v11L11+gap*LR+eLA. Section6 applies:
strict floors(1790351,1974557)/10^12 exceed u+v9+v11; relaxed
floors(1389839,1537579)/10^12 exceed u, strictly. Its shared
coefficient is2Z; U2,U3,J exceed X,3X/2,X-A/5. The ordinary table
plus Section12's shifts covers
q*c7[a]+v9*c9[j]+v11*c11[j]+gap*cR(a,b)+e*cA(a,b).
All nonzero rows and source/overlap/role guards are strict; cancellation
gives

```
Cminus=6q+9v9+11v11+12u-3e=1791867437/50000000000000,
Cplus =6q+9v9+11v11+12u-3e=143990663/4000000000000.
```

Both bands analyze one ordering: Section6's masks,names and key lowering
are fixed; only scalar allocations depend on sign(t). Original-field
threshold monotonicity retains zero-threshold Price before advancing.
Thus for |t|<=T the expected advanced cover is
(U-min(Cminus,Cplus))n+o(n),U=.13064756648.

For |t|>=T, Section9's same-threshold outer comparison plus safe
advances pays O=17930855452501/500000000000000000>Cminus.
Adjacent and no-arm even D6,D8,D10 baselines extend to T; the tighter
.01 adjacent contribution is only its floor difference. Other arms
retain their horizons. All sign/distance/arm classes are disjoint;
outer credit spends no central Price resource. Exactly

```
gamma=1000000000/7656280233,
gamma-(U-Cminus)=2719189980821/382814011650000000000000>0.
```

Fix one common source, radius, positive tie noise, mesh/tails and
concentration accuracy below this margin and every strict row/source
gap, then graph size. Section2 gives cutwidth below gamma*n with
fixed slack. Small components have bounded width; disconnected ones
concatenate. The unchanged actual-length full-B2 family transfer gives
1+1/(2gamma)=9656280233/2000000000=4.8281401165 for all sufficiently large
actual lengths, with zero explicit remainder.
