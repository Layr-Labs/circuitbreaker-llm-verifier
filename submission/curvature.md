# Three disjoint pair sets in one covariance budget

Notation is that of `graph.md`: F is the full crossing matching, Y is
the actual positive radial base field repeated on a high-girth cover,
and X is an independent reference field. All limiting formulas use
`lambda=2sqrt(2)`, `c=3-lambda`, `w_d=(1+d/3)2^(-d/2)` and
`Z_s=sqrt(1-s)X+sqrt(s)Y`. The actual finite parameters are chosen after
the finitely many windows, meshes and strict scalar slacks in this file.
The graph size is chosen last. No independent-branch assertion is made
about Y. Its coordinates may coincide on different covering fibers.

## Source bounds and the Gaussian tools

For complete alpha profiles write `mu=alpha^-1+2alpha`. Outside its
root the tree row satisfies `sum_neighbors f=mu f`; the correction at
the root has norm tending to zero. On every coordinate,
`sum_neighbors f_neighbor^2<=(alpha^-2+2alpha^2)f^2`.
Taking norms within a fiber gives, up to that root correction,

    0<=sum_neighbors b_neighbor-mu b
      <=[sqrt(3(alpha^-2+2alpha^2))-mu]b.

In the critical limit the coefficient on the right is c. Consequently
the actual defect `R_v=sum_neighbors Y-lambda Y_v` has standard deviation
at most c. A sibling contrast `sqrt(3)(Y_a-Y_b)` has variance at most one,
since its covariance is at least w_2=5/6. Conditioning a centered Gaussian
at zero decreases its covariance matrix.

We use three elementary Gaussian facts. If g is standard normal and
W(g)>=0 is even, then
`E W(g+h)>=exp(-||h||^2/2) E W(g)`: write the translated density and
average its opposite linear exponents, using cosh>=1. For a centered
Gaussian vector H, weighting by `exp(-r||H||^2/2)` has mass
`det(I+r Cov(H))^-1/2`; completing the square also shows that all
conditional variances decrease. Hadamard bounds this mass below using
the individual variance upper bounds. Finally centered Gaussian strips
obey

    Pr[all |L_i|<a_i] >= product_i Pr[|L_i|<a_i].     (C1)

This is the centered Gaussian rectangle inequality, including singular
laws by adding vanishing independent noise. To see the direction,
condition on the last scalar: the earlier symmetric convex box has an
even log-concave Gaussian-translation probability and thus decreases in
the absolute value of that scalar. Its covariance with a centered-interval
indicator is nonnegative; induction gives the product bound. Log concavity
here follows by integrating the Gaussian density on translated convex
sets (equivalently the Gaussian form of the Prékopa integral inequality).
Only this classical centered-rectangle form is used.

## Finite matching certificates

A closed endpoint-rooted even path can be removed at a constant matching
cost in every assignment of its attachment vertex: replace any use of
the attaching edge by the path's perfect matching. Every boundary edge
of a certified path is required noncrossing, so its component is closed
in the full exterior, not merely in the displayed window.

For adjacent u,v, require opposite colors on the two outer neighbors
at each endpoint. For one orientation certify a closed even path on
the zero neighbor at u and the one neighbor at v. If C is the remaining
exterior cover cost, let a,b in {0,1} be its extra costs for covering the
other attachment roots. After deleting the two paths, the four costs
are `C+a,C+b,C+1,C+1+min(a,b)`, with the two equal-color corners first.
Thus `-Delta_uv F=2-max(a,b)>=1`. If the reverse orientation is certified
too, delete all four paths; only the central crossing indicator remains,
giving -Delta=2. For each orientation the path sizes and physical labels
are mutually exclusive. This proves an additive two-orientation bound
without an overlap penalty.

For a distance-two path u-w-v, require w's third neighbor to share its
color. At either endpoint choose type 0 (both outer roots opposite w)
or type 2 (one root same as w on a closed two-vertex path, the other
opposite). Remove all certified paths. Let M be the remaining exterior
matching number and A,B the individual exterior gains on restoring u,v
in w's color. Their joint gain C is at most A+B because these new
vertices belong to the same bipartition side. If both individual gains
vanish, an augmenting path cannot pass through the other newly unmatched
same-side vertex; if exactly one vanishes, a joint gain two contradicts
that zero upon removing the other new vertex. The both-one case is
immediate. The four corner costs are `M+1,M+1+A,M+1+B,M+C`, hence
`Delta=C-A-B-1<=-1`. Types (0,0),(2,0),(0,2),(2,2) and their physical
labels are disjoint. The first is the baseline event; the other three
can be credited separately.

For distance three u-a-b-v require opposite middle colors and matching
third-neighbor colors. Choose the same type 0/2 outward alternatives
at u and v relative to their own middle colors. Again remove closed
paths. Now the two exterior additions lie on opposite bipartition
sides, so C>=A+B. When A=B=1, maximum exterior matchings exposing the
two required neighbors can be combined along their alternating components:
desired opposite-color exposed vertices cannot be the two ends of one
even alternating path. If either gain is zero the other gain alone
suffices. With flip coordinates relative to the two middle colors, the
corner values are `M+1+C,M+1+B,M+1+A,M+2`. Their mixed difference is
`1+C-A-B>=1`. Exactly one flip coordinate reverses the ordinary bit,
so the ordinary mixed difference is at most -1. All these arguments
permit arbitrary connections in the rest of the crossing graph.

## Positive strip integrals

Write `H(x)=Pr[|N(0,1)|<x]`, `r=s/(1-s)`. Let f_A,f_D,f_0 be the
absolute-normal densities with variances 1/12,1/3,1/33 respectively, and
let `phi_sigma` be the signed density with variance sigma^2=1/12.
Define the shared defect bounds

    ell_d=(c/2) sum_(j=0)^(d-1)(j+1)2^(-j/2).

They follow by solving the linear sibling midpoint recurrence. The
positive impulse coefficients are `(j+1)2^(-j/2-1)`, so the triangle
inequality for standard deviations bounds each actual shared pair error.

For a closed even n-path, with a_0=0 and positive alternating magnitudes
a_1,...,a_n, put

    D_n(a)=2phi_sigma(a_1)
       product_(j=1)^(n-1) 2phi_sigma(a_(j+1)+sqrt(2)a_j+a_(j-1)/2),
    M_n=lambda a_n+a_(n-1),
    P_n(s)=integral_(a>=0) D_n(a)
       product_(d=1)^n H(a_d/(ell_d sqrt(r)))
       integral_0^(M_n) f_D(z)H((M_n-z)/(2ell_(n+1)sqrt(r))) dz da.

Each factor two sums the two possible physical continuation labels. The
terminal integral tests the two closing children. Every required sibling
on an internal step has larger magnitude than the selected child, so
the same strip suffices for both. Independent contrasts here belong
only to the ideal critical tree field: Gaussian orthogonality under
sibling swaps proves their independence and the recursion.

Fix physical labels and sizes before conditioning or tilting. Given
`Z_u=Z_v=0` and Y, the independent ideal zero-pair contrasts retain their
standard laws. The normalized actual scores have the ideal reconstruction
in `g+sqrt(r)H_Y` plus shared defect errors. The robust two-orientation
event is even in g. Translation gives the centered quadratic tilt in
H_Y; Hadamard and the variance-one bounds give factor
`(1-s)^((n+m+2)/2)`. In that tilted law each correlated defect strip is
paid by (C1) using the realized ideal margins. The ideal two half-fields
factor only after this step. Set

    B_2=(1-s)P_2(s),
    E_8=sum_(n=4,6,8)(1-s)^(n/2)P_n(s).

The extra adjacent response excluding size pair (2,2) is at least
`2(1-s)(2B_2 E_8+E_8^2)`. This is a positive expansion in disjoint
size classes; it does not subtract two separately bounded probabilities.
The pair density is at least `3/(2pi)>.477`. Therefore the extra
adjacent credit is bounded below by

    K1_extra=.954 integral_0^1 (1-s)(2B_2 E_8+E_8^2) ds. (C2)

For a sharper size-(2,2) baseline use the larger error bounds
`tau=(.086,.208,.336)`, which strictly exceed ell_1,ell_2,ell_3. Put

    J(s)=integral_(x,k>=0) f_A(x)2f_D(2k+lambda x)
      H(x/(tau_1 sqrt(r))) H(k/(tau_2 sqrt(r)))
      integral_0^(lambda k+x) f_D(z)
        H((lambda k+x-z)/(2tau_3 sqrt(r))) dz dk dx.

The same six-driver tilt gives the baseline credit

    K1_base=.954 integral_0^1 (1-s)^3 J(s)^2 ds.       (C3)

For distance two, the critical conditional middle is N(0,1/33), its
third neighbor equals lambda times it, and the two outer pairs have
midpoint -X_w/2 and independent half-contrast `G/(2sqrt(3))`.
The actual normalized central residual coordinate
`sqrt(33)[Y_w-(4sqrt(2)/11)(Y_u+Y_v)]` has variance at most 43/11<3.91,
by expanding and using adjacent covariance floors and ceiling one.
The outer contrasts have variance at most one. Define

    L0(s,x)=integral_0^x f_D(z)H((x-z)/(c sqrt(r))) dz,
    L2(s,x)=integral_(a,b>=0) 2f_D(x+2a)2phi_sigma(b+sqrt(2)a)
      H(a/(ell_1 sqrt(r)))H(b/(ell_2 sqrt(r)))
      integral_0^(lambda b+a) f_D(z)
          H((lambda b+a-z)/(2ell_3 sqrt(r))) dz db da,
    D2(s)=(1-s)/sqrt(1+3.91r).

The root density `2f_D(x+2a)` includes both root-label choices and the
Jacobian. It is not multiplied by a second independent root density.
The old no-attachment and new attachment credits are respectively

    K2_base=.2878 integral_0^1 D2(s) integral_(x>=0)
       f_0(x)H(lambda x/(c sqrt(r))) L0(s,x)^2 dx ds,
    K2_extra=.2878 integral_0^1 D2(s) integral_(x>=0)
       f_0(x)H(lambda x/(c sqrt(r)))
         [2(1-s)L0(s,x)L2(s,x)+(1-s)^2 L2(s,x)^2] dx ds. (C4)

The half-normal f_0 already includes both signs, so no extra factor two
appears. The extra terms correspond exactly to the three physical
attachment classes above. The three core translated drivers give D2;
each closed two-path contributes two more variance-one contrasts and
thus a factor 1-s. The Gaussian tools and defect impulses already proved
give every strip in (C4), without assumptions on correlations in Y.

For distance three the critical middle pair conditioned on the endpoints
has variances 1/18 and covariance sqrt(2)/36. Its independent sum/difference
modes have variances `v_+=(2+sqrt(2))/36`, `v_-=(2-sqrt(2))/36`.
The regression rows on u,v are `(1/sqrt(2),1/3)` and its reversal.
Put `d=1/sqrt(2)+1/3`, `e=1/sqrt(2)-1/3`. Expanding the actual
residual sum and difference gives upper variances

    4+4d^2-4d(w_1+w_2),
    2(1-w_1)+2e^2(1-w_3)+4e(1-w_1).

Dividing by 2v_+,2v_- bounds the whitened actual mode variances by 5,9.
All signs of expansion coefficients are retained when replacing a
covariance by its floor or ceiling. In the two opposite-sign middle
orientations, the joint magnitude density is

    f(x,y)=(18sqrt(2)/pi)exp[-18(x^2+y^2+sqrt(2)xy)].

It already sums both orientations. With `V(s,x)=L0(s,x)+(1-s)L2(s,x)`,
the resulting credit is

    K3=.225 integral_0^1 (1-s)/sqrt((1+5r)(1+9r))
      integral_(x,y>=0) f(x,y)
       H((lambda x+y)/(c sqrt(r))) H((lambda y+x)/(c sqrt(r)))
        V(s,x)V(s,y) dx dy ds.                       (C5)

Here .225 is below the zero-pair density sqrt(2)/(2pi). The four core
drivers consist of the two whitened middle modes and two endpoint
contrasts; Hadamard gives the displayed determinant factor. Extra branch
choices contribute their factors 1-s. The two actual middle third-neighbor
errors are R_a,R_b and give the two explicit strips. Fixing each finite
branch class before applying its Gaussian tilt justifies summing their
disjoint contributions afterwards.

## Exact lower sums and finite transfer

`certify.py` implements only positive lower sums for (C2)-(C5). Every
decreasing spatial density uses an upper cell endpoint; every increasing
margin/strip uses a lower endpoint. Negative omitted regions never occur.
For all integrals the scalar lower rules are:

* `1.41421<sqrt(2)<1.41422`, `c<.17158`, `1/sqrt(2)<.70711`;
* `exp(-u)>=(sum_(j=0)^7(-u/64)^j/j!)^64` on the used range,
  computed with a downward fixed-point floor before six squarings;
* for H, round the argument down to .001 and cap it down at 3.3,
  then use `.7978 sum_(j=0)^23(-1/2)^j t^(2j+1)/(j!(2j+1))`;
* use density normalization lower constants 2.7638,1.3819,4.5834,
  and 8.1027 for f_A,f_D,f_0,f respectively. Machin's arctangent
  bounds for pi and rational radical bounds check these constants;
* all positive arrays use integer denominator 10^18 with downward
  products and divisions. Final cell integration uses exact Fractions.

For (C3) use spatial width 1/200 with `0<=x<=.7,0<=k<=.9`, and time
endpoints `(j/80)^2`, j=1,...,79. The inner z density is sampled at its
RIGHT endpoint. Its remaining margin is rounded down. The same fine
time/spatial grid evaluates K2_base on `0<=x<=.8`. The upper density
endpoint is essential: sampling the left endpoint would not be a lower
bound. The certified usable floors are

    K1_base >=3607233186504135/10^18,
    K2_base >=2203451254088009/10^18.                 (C6)

For (C2) use path cells j/80, j=0,...,47, maximum even size eight, and
time endpoints `(j/24)^2`, j=1,...,23. The positive AR2 contraction
keeps the last two amplitude indices; the transition has density
`2phi_sigma(a_new+sqrt(2)a_current+a_previous/2)`. For the terminal
convolution use lower `m=floor(2.82842*a_index+previous_index)`; the
z-cell ending at q/80 contributes density there and strip margin m-q.
The first transition uses a_0=0 exactly. For (C4) retain the same
48-cell arm arrays and 64 central x-cells. Its type-2 root density is
`2f_D((i+1+2(j+1))/80)/80`. For (C5) reuse these exact arm arrays,
and use 64 cells in each middle coordinate. Its density is sampled at
both upper endpoints; its two strip margins use
`floor(2.82842*i+j)` and its transpose. Integer square roots bound every
determinant downward. The fixed exact checks prove strict inequalities

    K1_extra>143/10^6,
    K2_extra>14/10000,
    K3>216/10^6.                                   (C7)

All retained integrand factors decrease with s, so right time endpoints
give lower rectangles. Amplitudes outside the stated boxes, larger
branches and the last interpolation tail have nonnegative coefficients.
For finite source transfer, fix these windows, classes, meshes and strict
slacks first. A putative sequence violating a uniform bound has a
subsequence of its fixed-dimensional actual covariance matrices and
interpolation values converging to a law covered above. The conditioned
pair stays nonsingular on the retained compact s interval, and the
independent ideal field gives every tested exterior score positive
conditional variance. Strict sign cylinders therefore converge in
probability. Their limiting lower bound contradicts the violation.
This chooses one sufficiently accurate alpha,R,jitter and cover girth
uniformly before graph size, with arbitrarily small loss per vertex.

## Summing the pair sets

Let S_d be total covariance excess above w_d on unordered distance-d
pairs of the high-girth cover; its order is N. The repeated base key
rows form a positive semidefinite Gram matrix C, whether or not their
covariances commute with adjacency A. With lambda=2sqrt(2),

    S_2>=2lambda S_1,
    S_3>=19 S_1.                                   (C8)

For the first identity put `L_v=sum_(u~v)b_u-lambda b_v`. Expansion gives
`sum_v<b_v,L_v>=2S_1` and
`2S_2=4lambda S_1+sum_v||L_v||^2`, since
`lambda^2=3+6w_2`. For the second use the distance-three polynomial
`P_3(A)=A^3-5A`, valid when girth exceeds six. There are five
non-geodesic length-three walks per adjacent ordered pair. On the
physical spectrum [-3,3],

    P_3(t)-P_3(lambda)-19(t-lambda)
       =(t-lambda)^2(t+2lambda)>=0.

Taking its trace against the positive Gram matrix C proves the result,
because `3w_1=lambda` and `12w_3=P_3(lambda)`. Replacing the critical
covariances in these formulas by finite reference covariances changes
them by an arbitrarily small uniform multiple of N. All actual
base-minus-reference differences remain nonnegative by fiber comparison.

Apply (G2) once, retaining these three DISJOINT sets of unordered pairs,
and discard its other nonpositive terms. Equations (C6)-(C8), with
`2lambda>5.65684`, give the total adjacent-excess coefficient

    3607233186504135/10^18+143/10^6
       +5.65684*(2203451254088009/10^18+14/10000)
       +19*(216/10^6)
      >K=28238380378679347/10^18.                   (C9)

For each adjacent base edge there are exactly d lifted edges in a
degree-d cover. Exact matching lifting divides the integrated expectation
gain by d. The resulting original-base estimate is
`U*n-f_Y(0)>=K E-arbitrarily_small_error*n`. All finite-source debits
are paid by the strict scalar slacks before graph size. No term from
a particular distance set is reused as a separate gain.
