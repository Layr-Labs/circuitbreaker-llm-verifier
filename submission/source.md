# Local probability bounds for the actual radial source

The source and notation are defined in `graph.md`. The original finite
base Gaussian keys are denoted Y. Its full-alpha adjacent excesses D_e
are used for marking before finite truncation and jitter. Every limiting
inequality below has a uniform arbitrarily small finite-source allowance:
first choose alpha close to 1/sqrt(2), then a fixed finite radius and
positive jitter, before choosing the graph size. No repeated base noise
is replaced by fresh tree noise. The rational checks described here are
implemented independently in `certify.py`.

## Fiber rigidity, defect bounds and incidence constraints

Let u_x(p) be the unit direction of a full radial row on fiber p, and
b_x(p) its norm. Then, for neighboring lifted vertices,

    2D_xy=sum_p b_x(p)b_y(p)||u_x(p)-u_y(p)||^2.

Adjacent row entries, and hence adjacent fiber norms, have ratios between
alpha and alpha^-1. Telescope the directions along a d-edge tree path,
apply Cauchy--Schwarz to its d summands, and use
`b_x b_y/(b_(i-1)b_i)<=alpha^(-(d-1))`. This proves that a path whose
edges all have excess at most eta satisfies

    w_d-o(1)<=Cov(Y_x,Y_y)<=w_d+k_d eta+o(1),
    k_d=d^2(sqrt(2))^(d-1).                          (S1)

The positivity calculation in `curvature.md` gives, at a vertex with a
light complete star, for `R_v=sum_neighbors Y-lambda Y_v`,

    Var(R_v)<=3c eta,
    0<=Cov(R_v,Y_x)<=3eta (sqrt(2))^dist(v,x)          (S2)

in the limiting patch laws, where lambda=2sqrt(2), c=3-lambda.
Indeed the corrected nonnegative defect row is at most c b_v and has
inner product at most 3eta with b_v. This proves the variance bound.
The coordinatewise path ratio `b_x<=alpha^-d b_v` proves the second
bound. The root correction has norm tending uniformly to zero.

Two source-specific constraints will be needed in `proof.md`:

    D_e<=2 sum_(other edges at v)D_f+arbitrarily_small_error,
    sum_(four off-edge incidences at uv)D_f
          >=(8/5)D_uv-arbitrarily_small_error.        (S3)

Here is a derivation including finite-alpha root masses. Put a=alpha,
b=1/alpha, and
`f(x)=sqrt(a^2+(b^2-a^2)x)-a-(b-a)x`. In one fiber about v let T be
its total squared mass, C its root atom, and A_1,A_2,A_3 its three branch
masses. The excess at incident edge i is exactly `T f(A_i/T)`.
Rationalizing gives

    f(x)=(b-a)^2 x(1-x)
       /[sqrt(a^2+(b^2-a^2)x)+a+(b-a)x].

The denominator lies between 2a and 2b, so
`f(x)<=a^-2 f(1-x)`. Concavity gives
`f(y)+f(z)>=f(y+z)` when y+z<=1. Its Lipschitz constant is
`L_alpha=(b-a)^2/(2a)`. Since branch fractions sum to 1-C/T, summing
these inequalities over fibers proves

    D_i<=a^-2(D_j+D_k)+a^-2 L_alpha c_alpha.           (S4)

For a paired edge use x=A/T in the branch from v towards u and
`h=a^2+(b^2-a^2)x`. In the reverse direction the fiber mass is T'=Th
and its branch fraction is `x'=a^2(1-x)/h`; both endpoint atoms are
included correctly by this complement identity. The four off-edge
terms are at least
`T[f(1-x)+h f(1-x')]-L_alpha(C_v+C_u)`.
Write s=a+b, z=a^2+b^2, r=sqrt(h), y=r+1/r. Then 2<=y<=s,
`f(x)/r=1-y/s`, and the bracket divided by r equals

    R(y)-zy/s,
    R(y)=sqrt(z^2-(w-1)^2),
    w=sqrt((z+1)^2-zy^2), 1<=w<=z-1.

For `7/10<=alpha<1/sqrt(2)`, one has 5/2<=z<3, and
`((w-1)^2)''=2z((z+1)^2-w^3)/w^3>0`.
Thus R is concave. The chord between R(2)=2sqrt(z-1) and R(s)=z gives
`f(1-x)+h f(1-x')>=C_alpha f(x)`, where

    C_alpha=[2s sqrt(z-1)-2z]/(s-2).

Summing the endpoint atoms gives precisely `2L_alpha c_alpha` of debit.
The limit of C_alpha is
`9sqrt(6)+12sqrt(3)-15sqrt(2)-20>8/5`; the strict inequality reduces
by squaring positive sides to 289>288. Equation (S4), its root debit,
and this paired inequality prove (S3) uniformly before graph size.

## The two strict ordering cylinders

For the seven-key event name root o, its neighbors a,b,u, the other
neighbors a1,a2 of a, the other neighbors w,w' of u, and the other
neighbors w1,w2 of w. Its nine inspected edges are the complete stars
of o,a,u,w. The positive support is
`{a,b,a1,a2,w,w1,w2}`; the guarded event puts o,u,w' above threshold t
and all seven support keys below t. At a root there are twelve names,
each assigned weight 1/2. A successful guarded event fixes u and w,
so its total weight is at most the indicator that the root moves.

For a nine-key event use the fourteen-vertex tree with edges

    oa,ob,ou,aa1,aa2,bw,bw',uz,uz',ww1,ww2,zz1,zz2.

Its positive support is `{a,a1,a2,w,w1,w2,z,z1,z2}`. The guards
o,b,u,w',z' are above t. Its thirteen inspected edges are the complete
stars at o,a,b,u,w,z. Choose a in three ways and w,z in two ways each:
the twelve guarded names are disjoint, because a is the unique original
below-threshold neighbor of o and w,z are unique below-threshold children.
They are disjoint from the seven-key guarded events, which have two
below-threshold neighbors of o. The light cutoffs are

    eta7=19/40000, eta9=19/80000.

Equation (S1) ensures faithfulness of each light projected patch: for
its distinct lifted vertices at distances up to six, `w_d+k_d eta<1`.
Equal physical projections would have covariance one, a contradiction.

Both positive supports certify a strict one-unit decrease of F, even
without guards. For seven keys, remove o. Its neighbor a is isolated
in the original crossing graph, so oa adds one to the exterior optimum.
After flipping o the only possible new crossing edge is ou; if present
it may be replaced by uw, since w has no other crossing neighbor.
For nine keys, after removing o the two branches b-w and u-z are
closed two-vertex components, while a is isolated. Originally oa costs
one; after flipping o both new attachments are covered at no extra
cost by the b,w and u,z covers. This first proves the guarded case.
Restoring any omitted guard decreases the root marginal by submodularity,
and one-bit Lipschitzness prevents a marginal below -1, proving the
guard-free support assertion.

## Seven-key Gaussian envelope

On its actual light patch use `O=Y_o=t+x`,
`R=2Y_u-Y_a-Y_b`, `Q=sqrt(3)(Y_w'-Y_w)`. After regressing R,Q on O,
write their residual vector as `K_2(y,z)^T`, with x,y,z positive.
The density in these coordinates is
`(2pi)^(-3/2)sqrt(det K_2) exp[-((t+x)^2+(y,z)K_2(y,z)^T)/2]`.
Using (S1), its lower envelope is

    .0621 exp[-((T+x)^2+1.038y^2+1.033z^2)/2].       (S5)

The six upper signed conditional means, at T=.055, are the following
integer rows divided by 10^5. The first four mean score minus t;
the last two mean t minus score.

| Tested scores | x | y | z | constant |
|---|---:|---:|---:|---:|
| a,b | 94329 | -16039 | 1921 | 315 |
| a1,a2 | 83603 | -21780 | 4113 | 917 |
| w | 83603 | 47316 | -28012 | 917 |
| w1,w2 | 71567 | 50606 | -39010 | 1611 |
| u | -94280 | -32705 | 522 | 315 |
| w' | -83332 | -45350 | -28012 | 917 |

Here are reproducible covariance bounds for that table. The signed ideal
R coefficients are `(-1/6,-lambda/12,lambda/6,1/2,-1/3,-lambda/6)`.
Their upper error multipliers are
`(2k2,2k3,2k1,2k2,2k2,2k3)eta`, plus the O-regression cost
`2eta(w_d+k_d eta)`. The signed ideal Q coefficients are
`(0,0,-1/(2sqrt(3)),-lambda/(4sqrt(3)),0,-1/(2sqrt(3)))`.
Their error multipliers are `sqrt(3)(k3,k4,k2,k3,k1,k2)eta`, plus
`sqrt(3)k2 eta(w_d+k_d eta)`. The root coefficients use the upper
`w_d+k_d eta` for the first four rows and lower w_d for the last two.
The constant is `(1-w_d)T`. The corresponding lower residual variances
are `1-8k2 eta-(2eta)^2` and `1-6k2 eta-(sqrt(3)k2 eta)^2`; the
absolute covariance is bounded by
`sqrt(3)(2+2k3)eta+2eta*sqrt(3)k2 eta`. These yield (S5) by the
determinant and quadratic-form bounds. `certify.py` checks every rational
enclosure with strict slack.

For T=.03 or .02 replace each constant h by `ceil(h*T/.055)`; all other
coefficients stay fixed. The exact rows of constants are respectively
`(172,501,501,879,172,501)` and `(115,334,334,586,115,334)`.
Centered defects at o,a,u,w have deviation at most .01565 by (S2).
For common tolerances r=.018,.02504,.03365 their centered strip masses
exceed .748,.889,.967. The three centered sibling half-differences have
deviation at most 1/sqrt(12)<1/3.46. The star equations bound the six
score fluctuations by

    r/3+d1,
    1.415(r/3+d1)+r/2+d2,
    .972r,
    2.043r+d3,
    r/3,
    .972r.

Subtract these from the negative upper means. Inequality (C1) in
`curvature.md` pays the product of four defect and three half-difference
strips, regardless of dependence. Use the lower strip polynomial
`.797(v-v^3/6+v^5/40-v^7/336)` capped at v=1.6, with standardized
arguments rounded down to 1/10000. On each box try the two choices
d1=.5U,.8U, where `U=min(A-r/3,(B-r/2)/1.415-r/3)` for its first two
lower margins A,B. Take the largest compatible d2,d3. Maximizing among
these valid rectangles is legitimate.

The finite lower sum uses widths `(1/250,1/50,1/50)` on
`[0,.2]x[0,3.6]x[0,6.8]`. Densities use upper endpoints; each signed
affine mean uses its worst box corner. The density exponential is
bounded below by `(1-u/4096)^4096` with downward squarings. Nonpositive
margin boxes and the entire exterior are simply omitted. The original
implementation in this submission independently computes these positive
sums. It may stop only after a complete x-slice when its exact partial
mass already exceeds the corresponding used floor; every remaining
box has nonnegative contribution. This proves the strict usable floors

    q=.000000273249 for |t|<=.03,
    q=.000000283 for |t|<=.02.                       (S6)

## Nine-key Gaussian envelope

Use actual drivers

    O=Y_o=t+x, R=Y_b+Y_u-2Y_a=y,
    Q_b=sqrt(3)(Y_w'-Y_w)=z,
    Q_u=sqrt(3)(Y_z'-Y_z)=h,

with x,y,z,h positive. These four are independent standard normals
only in the ideal critical law. For actual conditioning put
`gamma=(2eta,sqrt(3)k2 eta,sqrt(3)k2 eta)` with eta=eta9. The covariance
K of the last three drivers conditional on O is bounded below in its
diagonal by
`(1-8k2 eta-gamma0^2,1-6k2 eta-gamma1^2,1-6k2 eta-gamma2^2)`.
Its off-diagonal absolute bounds are
`r=sqrt(3)(1+3k3)eta+gamma0 gamma1` and
`q=6k4 eta+gamma1 gamma2`. Let M have those lower diagonals and
negative off-diagonals (-r,-r,-q). It is strictly diagonally dominant.
Its inverse W is entrywise positive, and a Neumann-series comparison
proves `|K^-1|<=W`. With `s=x+T`, `U=(y,z,h)+gamma s`, the driver
density is bounded below by

    .0252 exp[-(s^2+U^T W U)/2].                    (S7)

Hadamard uses `det K<=1+2k2 eta`; the prefactor and all matrix bounds
are checked with rational radical uppers sqrt(2)<1.415,sqrt(3)<1.733,
`k2<5.66,k3=18,k4<45.28`.

The conditional means of the six star defects, the root half-difference
`(Y_b-Y_u)/2`, and three terminal half-differences must all be retained.
Their following formulas specify the entire envelope. For a defect at v,
put `p_j=1.415^dist(v,j)`, `c0=p_o`, and

    c_v=(max(p_b+p_u,2p_a),1.733 max(p_w,p_w'),
                                     1.733 max(p_z,p_z')).

Then its absolute mean is bounded by
`B_v=3eta[c0 s+(c_v+c0 gamma)^T W U]`, by (S2).
For a half-difference with root covariance bound h0 and raw driver
bounds hvec, the absolute mean is at most

    B_half=h0 s+(hvec+h0 gamma)^T W U.               (S8)

For the root half-difference use
`h0=eta/2`, `hvec=(k2 eta,1.733(1+k3)eta/2,1.733(1+k3)eta/2)`.
For a terminal pair l1,l2 use
`h0=eta*k_dist(l1,o)/2` and

    hvec=(eta/2)(k_dist(l1,b)+k_dist(l1,u)+2k_dist(l1,a),
                1.733(k_dist(l1,w)+k_dist(l1,w')),
                1.733(k_dist(l1,z)+k_dist(l1,z'))),

with `k0=0,k1=1,k2=5.66,k3=18,k4=45.28,k5=100`. Equal ideal distances
to each sibling cancel the ideal covariance, which proves (S8).

Before adding these mean corrections, the ideal midpoint gaps are

| Scores | midpoint minus t |
|---|---|
| a | rho*x-y/3+(rho-1)t |
| b,u | rho*x+y/6+(rho-1)t |
| a1,a2 | (5/6)x-lambda*y/6-t/6 |
| w,z | (5/6)x+lambda*y/12-(z or h)/(2sqrt(3))-t/6 |
| w',z' | same with the contrast sign reversed |
| w1,w2 or z1,z2 | x/sqrt(2)+y/4-lambda*(z or h)/(4sqrt(3))+(1/sqrt(2)-1)t |

These are algebraic star reconstructions, not a Markov assumption about
Y. Their defect-mean upper corrections are respectively

    B_o/3;
    B_o/3+B_root_half;
    lambda B_o/6+B_a/2;
    lambda B_o/6+lambda B_root_half/2+B_b/2 (or B_u/2);
    B_o/2+3B_root_half/2+lambda B_b/4+B_w/2
           (or the reflected expression at z).

For each terminal pair also subtract its own B_half from (S8).
For centered defects with common radius r and centered half-differences
d0,da,dw,dz, the required centered score deviations are bounded by
`r/3`, `r/3+d0`, `.972r+da`, `.972r+1.415d0`, and
`1.708r+1.5d0+(dw or dz)` in the corresponding rows.
Use the same three r choices and masses as above, now for six defects.
Use the same lower strip polynomial for the four half-differences.
Inequality (C1) pays their product, including singular laws.

The fixed rational lower sum has widths `(.02,.1,.1,.1)` on
`[0,.2]x[0,3]x[0,4]x[0,4]`, 480,000 boxes. Upward-rounded mean
coefficients and inverse entries have scale 10^9. For each box choose
da as its available maximum, and d0 as 1/4,1/2 or3/4 of its admissible
maximum; dw,dz are then their available maxima. Discard any infeasible
choice. (S7) uses the upper box endpoints and the same exact exponential
lower rule. This gives the strict floor

    p=247011/10^12 for every |t|<=.03.               (S9)

All values in (S5)-(S9) concern actual finite-base laws. Strict scalar
gaps and covariance compactness transfer the limiting envelopes uniformly:
the driver blocks stay invertible, each score threshold has a continuous
Gaussian marginal, and all named key comparisons have positive variance.
A violating sequence would converge to one of the already paid limiting
patch laws, contradicting the strict lower sum. Thus one set of finite
source resources works for all light names and all central thresholds.

## Relaxing only the four terminal child edges

Keep cutoff theta=19/80000 on the complete core stars o,a,b,u of the
nine-key patch, but allow eta=2theta on the two outward child edges
at each of w,z. Every path among the eight scores in the four drivers
uses core edges only. Thus gamma,M,W, the driver density, and the root
half-difference regression in (S7)-(S8) are unchanged.

For a star with unequal edge bounds, the proof of (S2) retains their
sum. At w,z that sum is at most 5theta instead of 3theta. Multiply
their defect mean bounds by 5/3. Their variance upper bound is now
5c theta, still below .01565^2, which is checked exactly. All other
defect means are unchanged.

Every path from an ending child to a driver vertex contains at most one
edge with the larger cutoff. Keeping the individual terms in the direction
telescope proving (S1) gives covariance error at most

    d(sqrt(2))^(d-1)sum_path D_e
       <=(1+1/d)theta k_d.

Consequently multiply k_d by 1+1/d only in the two ending-pair mean
rows of (S8); keep the a1,a2 row unchanged. Equal ideal distances still
cancel the ideal covariance in each half-difference. Their centered
variances stay at most 1/12, and every named projection remains faithful
by the larger eta bound through diameter six.

All sign and density formulas are otherwise identical. The independent
standalone program reruns the same positive box lower sum with these
two changes and stops as soon as its exact partial mass exceeds

    u=301/10^10.                                   (S9R)

It returns a strict rational partial bound above u; every unvisited box
has nonnegative mass. This stopping rule does not infer anything from
a timeout or sample. It is an exact finite certificate with a proved
nonnegative remainder. The same covariance compactness and strict sign
boundary argument makes this an actual-source floor, uniformly at
|t|<=.03 with fixed finite resources before graph size.

Let LR be its eligible name count. Since LR>=L9, pay
`(p-u)L9+u LR` for this same collection of disjoint guarded nine-key
events: on old-eligible names the probability floor is p, and on the
extra names the floor is u. This uses a convex split of p on the old
names, not two copies of their gain. The strict support certificate is
unchanged, and its negative original marginal keeps it disjoint from
the neutral groups. The ordering still inherits keys within radius three.

The program pairs nine-key boxes under the z/h branch swap. It checks
that gamma and the inverse majorant are invariant, the o,a and root-half
rows are fixed, and the b/u, w/z and terminal-half rows exchange exactly
after rational rounding. The density quadratic, minimum margins, strip
products and tolerance maximization are invariant. Thus summing only
iz<=ih, with weight two off the diagonal and one on it, reproduces the
full original lower sum exactly. For the relaxed early stop each retained
one- or two-box orbit is nonnegative, so the same stopping proof applies.

## The actual-source weak-pair cone

Use the twelve-vertex patch with m adjacent to r,s,h, r adjacent to a,a',
s adjacent to b,b', and a,b each adjacent to their two terminal children.
The five complete stars are m,r,s,a,b, totaling eleven edges. All must
be light at theta=19/80000. Require m,r,s,h above t and all other
vertices below t, with a<a',b<b' selecting the private children.

Take the seven actual driver forms

    Y_m,
    2lambda Y_m-3(Y_r+Y_s),
    sqrt(3)(Y_r-Y_s), sqrt(3)(Y_a-Y_a'), sqrt(3)(Y_b-Y_b'),
    sqrt(3)(Y_a1-Y_a2), sqrt(3)(Y_b1-Y_b2).

Their ideal critical covariance is the identity, by the orthogonal
sibling-contrast decomposition. For each two signed coefficient rows,
bound its covariance error by summing separately positive and negative
products against `theta*k_distance`, taking the larger absolute endpoint.
This gives an entrywise majorant P with row sums below one. The inverse
majorant `(I-P)^-1` bounds the actual inverse by its Neumann series;
the product of upper diagonal variances bounds its determinant. All
these rational matrix calculations are reproduced by the program.

At threshold t parameterize the drivers through

    Y_m=t+z, Y_r=t+xz, Y_s=t+yz,
    .23<=z<=.8, 0<x,y<=.12,
    -.65z<G_r,G_s<0, |J_a|,|J_b|<.75z.

The Jacobian is `6sqrt(3)z^6`. Uniformly for |t|<=.03, absolute driver
amplitudes divided by z are bounded by
`(1.131,5.703,.209,.65,.65,.75,.75)`. The inverse/determinant controls
bound the driver density below by `exp(-22z^2)/650`.
Regression of the five actual defects using (S2) gives absolute means
bounded by `(.065,.070,.070,.094,.094)z`. Restrict all their centered
parts to .1z; each mass exceeds .95 since z>=.23 and (S2) places it
beyond two standard deviations. Their joint mass is at least .95^5
by (C1), without requiring independent defects.

The h sign follows from `Y_h=lambda Y_m-Y_r-Y_s+R_m`. At r, the
a,a' midpoint is `(lambda Y_r-Y_m+R_r)/2` and the negative G_r makes
a the lower child. Its terminal midpoint is
`(lambda Y_a-Y_r+R_a)/2`. Using the displayed driver, mean and strip
bounds makes both child gaps and both terminal gaps strictly negative;
the s side is identical. The adverse threshold terms are at most
`c|t|/2` per midpoint. The program verifies the strict rational sign
inequalities with radical enclosures.

Integrate z^6 exactly on the 57 cells [.23,.24],...,[.79,.80], evaluating
the decreasing exponential at each upper endpoint and using
`exp(-u)>=(1-u/256)^256`. Multiply by the volume
`6sqrt(3)(.12)^2(.65)^2(1.5)^2`, density factor 1/650 and .95^5.
The result strictly exceeds

    w=4695/10^12.                                  (S10)

The seven-driver covariance stays uniformly invertible. Scores and the
two selected-child differences have nonzero variance. The same compactness
argument transfers this strict cone lower bound to fixed finite source
resources uniformly over every light canonical weak-pair name.

## Threshold deficits on the cover and on the actual base

First consider the high-girth reference cover. In (G3) retain adjacent
and distance-two ordered pairs only. For adjacent pairs a zero-conditioned
private-leaf witness with margins .04 has probability above .018:
with f_A,f_D as in `curvature.md`, its probability is the square of

    integral_(x,k>=.04) f_A(x)2f_D(2k+lambda x)
                               H_D(lambda k+x-.08) dk dx.

Here H_D is the absolute N(0,1/3) distribution. The program sums .01
boxes with x=.04,...,.59 and k=.04,...,.79. Density lower functions
are `(11/4)(1-6x^2/64)^64` and
`(11/8)(1-(3/2)y^2/64)^64`; H_D uses the integrated seventh-order odd
Taylor lower polynomial, capped at 5/4. Upper density endpoints and
lower margins give the stated strict .018 bound. The two sign-reflected
orientations count with curvature two on overlap. At moving endpoint
value s their depth-one/two/three relative threshold shifts have absolute
coefficients below .086,.208,.336. For |s|<=.03 this preserves the
zero-pair signs with positive margin. Finite conditional convergence
therefore gives response above .0358.

For distance two, given endpoints equal to s, write
`Y_w=(4lambda/11)s+Z`, with Z~N(0,1/33) in the critical reference law.
The third-neighbor gap is `lambda Z-s/11`; the two endpoint sibling
midpoints have the opposite sign on the event
`|Z|=x>=.01`, with both absolute N(0,1/3) contrasts below x-.008.
The critical probability is bounded below by

    integral_(.01)^(.60) sqrt(66/pi) exp(-33x^2/2)H_D(x-.008)^2 dx.

Using .005 cells, upper density endpoints, lower strip endpoints,
`sqrt(66/pi)>4.58`, `H_D(z)>=(69/50)(z-z^3/2)` for z<=.6, and
`exp(-u)>=(1-u/128)^128`, the exact sum is strictly above .045.
The type-(0,0) distance-two matching certificate gives mixed response
at least one on this event, uniformly for |s|<=.03. Uniform finite
conditional convergence preserves the strict floor.

For adjacent pairs use `k_j>.056`, residual density >1.1 and three
ordered choices per root. For distance two use `k_j>.165`, residual
density >.71 and six choices. All bands have fixed slack near the
critical covariance. With phi(t)>.39 and T=.03, integrating (G3) gives

    M(0)-M(T)
      >[3*.39*.056*1.1*.0358+6*.39*.165*.71*.045]T^2/2*n
      =671223267/10^14*n >67/10^7*n.                (S11)

Evenness and monotonicity extend this to all |t|>=T. Comparing cover
and original base fields at that SAME threshold and dividing by cover
degree transfers (S11). A safe nonlinear ordering cannot increase F.

An additional, distinct deficit holds directly in the actual base field.
Take a light seven-edge patch consisting of the complete stars of a
distance-two path's endpoints u,v and middle w. Write h for w's third
neighbor, and a,b,c,d for the four outer neighbors. Its vertices are
faithful by (S1). Given `Y_u=Y_v=s`, |s|<=.03, let
`rho=Cov(Y_u,Y_v)`, `r1=Cov(Y_u,Y_w)`, `r2=Cov(Y_v,Y_w)`,
`m=(r1+r2)/(1+rho)`, `Z=Y_w-ms`. Sufficiently accurate finite choices
give

    .8332<rho<.8361, .9428<r1,r2<.9433,
    1.026<m<1.030, .029<Var(Z)<.032.                 (S12)

The variance identity is
`1-(r1+r2)^2/[2(1+rho)]-(r1-r2)^2/[2(1-rho)]`.
For `D_u=(Y_a-Y_b)/2`, (S1) gives covariance bounds against u,v,w of
.00024,.00428,.00135. The entries of
`Gram(u,v)^-1(r1,r2)^T` are below .52 in magnitude. After further
conditioning on Z=z,

    |E D_u|,|E D_v|<.13|z|+.003|s|,
    |E R_x|<.095|z|+.017|s|, x=u,v,w.

The first is Gaussian regression using Var(Z)>.029; the second follows
from regression Cauchy--Schwarz and sd(R)<.016. Centered half-differences
have sd<.29 and centered defects have sd<.016. On .2<=|z|<=.3 restrict
the two centered half-differences to .05, two endpoint defects to .02
and the middle defect to .1. The outer-pair sign margin is at least
`.3225|z|-.1125*.03-.06>=.001125`; the h margin is positive as well,
using `|lambda m-3|<.1`. Thus the type-(0,0) distance-two certificate
holds for every exterior.

The half-normal density of Z exceeds .93 on that interval; use
`4.45(1-(45/29)/128)^128>.93`. The two half-difference strip masses
exceed .1365, the two small defect masses exceed .78, and the last
exceeds .999. These follow from the same lower Gaussian Taylor
polynomials; for the last tail use `2exp(-6.25^2/2)<.001`.
By (C1), the conditional expected negative mixed difference is at least
`.093(.1365)^2(.78)^2(.999)>.001`.

There are 6n directed distance-two paths. Each physical edge is inspected
28 times by their seven-edge patches: 8 endpoint-path, 4 middle-third
and 16 outer positions. If h7 counts edges with D>eta7, at least
`max(0,6n-28h7)` directed patches survive. Their endpoint pairs are
distinct: a second common neighbor would violate the patch's faithful
distance-four bound. For these pairs (S12) gives k>.163 and residual
density >.71 on the interval. Equation (G3), phi>.39, and integration
therefore give for every real t

    f_Y(0)-f_Y(t)>=b min(t^2,.03^2) max(0,6n-28h7),
    b=451347/20000000000.                            (S13)

For |t|>.03 use monotonicity. This is an ACTUAL-base deficit and can
be added to the zero-threshold covariance gain by telescoping
`reference -> f_Y(0) -> f_Y(t)`. It is not a second charge to a Price
pair, and no comparison is applied to the nonlinear corrected objective.
