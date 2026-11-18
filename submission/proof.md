# Paired radial source constraints improve a full-B2 lower bound

For the single total norm-one language defined below, there exists N_0 such
that EVERY integer N>=N_0 satisfies C_B2(L_N)>=c_* N, where

    c_*=482709692616067094773114047344423/100000000000000000000000000000000,
    c_*=4.82709692616067094773114047344423.

All sixteen Boolean tables, shared DAGs, repeated pins and actual input
length are covered. The eventual subtracted remainder is zero; N_0 is
unspecified. The proof strengthens a Gaussian crossing-cover construction
by constraints on the actual radial source: excess on a heavy edge forces
additional excess on neighboring light edges. Closed even-path witnesses
also supply additional curvature. These gains are kept in one Price budget.
The graph coefficient is gamma=358496277066174326467971/2744000000000000000000000, and the unrounded circuit
coefficient is 1+1/(2gamma). No numerical screen is used as a theorem.

The graph/model/field baseline is reproduced from public source a901153f78d426225407c02dec82995e51192f6d,
building on galarnon42, Arthur Ostapenko, mooselumph, DarkWingMcQuack and
McGuire–Schlesinger–Williams. The two-pass scheduling argument reproduced in
cutwidth.md is public in Puthranps/cubic-cutwidth-retiming, research/two-pass-bound.md
at 6b299fde41307ff47eb4a699e140fd009e61a51d. Only that elementary scheduling
argument is used from the newer source, not its rejected distance-two
numerical floor. All mathematical premises used here are supplied below
or in the other declared evidence; no prior judgment is a premise.

## 1. The family: one evaluator on every bit string

### 1.1 Length parameters and input layout

For an integer N >= 16 define

    k  = the largest integer with 4^k <= N,      t = 2^k,
    m0 = ⌊(N − t)/(t + 1)⌋,   m = m0 if m0 is odd, m = m0 − 1 if m0 is even,
    ℓ  = N − t − m − t·m.

Thus t is a power of two, m is odd and gcd(m, t) = 1. An input x = x_1 … x_N is split, in this order, into the
**R-prefix** x_1..x_t, the **Q-prefix** x_{t+1}..x_{t+m}, the **payload** of t·m bits c_{u,v} (0 <= u < m,
0 <= v < t; c_{u,v} is bit number t + m + v·m + u + 1), and ℓ **leftover** bits. So N = t + m + t·m + ℓ exactly.

### 1.2 The evaluator A

On x ∈ {0,1}^N:

1. If N < 16, output 0.
2. Compute k, t, m, ℓ as above. If m < 1, output 0.
3. Let R(W) = W^t + Σ_{v<t} x_{v+1} W^v ∈ F_2[W] and Q(Z) = Z^m + Σ_{u<m} x_{t+u+1} Z^u ∈ F_2[Z].
   If IrrTest(R) or IrrTest(Q) rejects, output 0.
4. Let F = F_2[Z]/(Q) and E = F[W]/(R). Form z = Σ_{v<t} (Σ_{u<m} c_{u,v} Z^u) W^v ∈ E.
5. Let e = Σ_{j=0}^{t−1} 2^{m·j}. Compute z^e in E by square-and-multiply. Output 1 iff z^e = 1.

The leftover bits are never read. **IrrTest(P)**, for a monic P ∈ F_2[X] of degree d >= 1, accepts iff
X^{2^d} ≡ X (mod P) and, for every divisor d' of d with 1 <= d' < d, gcd(P, (X^{2^{d'}} − X) mod P) = 1, where
X^{2^j} mod P is computed by j successive squarings modulo P.

L = {x : A(x) = 1} and L_N = L ∩ {0,1}^N. This defines L_N at **every** length, including N < 16, the guard
m < 1 (which never fires for N >= 16, `circuit.md` Lemma 6.1(a)) and every invalid prefix. The polynomials R and
Q are input bits checked by the evaluator; they are not advice and not the output of a construction algorithm.
When both prefixes pass, E is a field with 2^{tm} elements (`circuit.md` Lemmas 2.2 and 2.4).

## 2. Uniform P and NP with empty witnesses

A runs in O(N^4) bit operations on every input (`circuit.md` §2.2), so L ∈ P. Let V(x, w) accept exactly when w
is the empty string and A(x) = 1, with witness bound W(N) = 0. If x ∈ L, the empty witness is accepted. If
V(x, w) = 1 for some permitted w, then A(x) = 1 and x ∈ L. This proves both directions with one fixed
polynomial-time verifier for all lengths. Witnesses are not circuit inputs.

## 3. Circuit model and coverage

The bound is for the pinned full B2 model (Li–Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3):
deterministic acyclic circuits with one output, two ordered input pins per gate, any of the 16 binary truth
tables per gate (constants, projections and negated projections included), unrestricted fan-out and depth,
repeated pins and shared subcircuits allowed, exact computation on every input of length N, nonuniform circuits.
Every gate is counted once; inputs, wires and output selection cost nothing. The formal model is `circuit.md` §1.

* **All gate types, repeated pins, sharing.** The only gate-dependent step is the local rule
  β(e_g) = φ_g(β(p_1), β(p_2)) of `circuit.md` Theorem 4.4, valid for every φ_g; degenerate gates are counted like
  all others. A gate with both pins fed by one node has two distinct pin edges (`circuit.md` Definition 4.1),
  counted twice in Lemma 4.2; their equal values are forced by the copy rules, not assumed.
* **No substitution into circuits.** The proof never fixes inputs of a circuit. It compiles the whole circuit
  (with its R-, Q-prefix and leftover inputs) into a read-once program and restricts the *program*
  (`circuit.md` Lemma 3.1, Corollary 4.5). Hence constant nodes created by substitutions never arise.
* **Constant nodes.** If the model admits free constant nodes at all, `circuit.md` Lemma 1.1(i) removes them by
  absorbing each constant into the truth table of the gate it feeds; a gate with two constant pins becomes a
  constant-type gate reading x_1 twice. No gate is added, no gate is removed except non-ancestors of the output
  (Lemma 1.1(ii)), and no uncharged basis element is introduced. If the model has no constant nodes this step is
  void.
* **Terminal cases.** L_N is neither constant nor a projection for N >= 256 (`circuit.md` Proposition 6.2(1)),
  so after Lemma 1.1 the output is a gate and every other gate has fan-out >= 1. Inputs of fan-out 0 get no
  vertex and are never read (Theorem 4.4); unused payload variables are counted by r (Proposition 6.2(2)). Wire
  graphs that are trees or unicyclic are handled in Proposition 6.2(4).
* **Every length, invalid prefixes.** A is defined on all strings and the lower-bound transfer holds at every integer N >= 256.
  Every circuit must compute L_N on all inputs, including invalid prefixes; the lower bound uses only the
  restriction of the compiled program to one valid prefix pair (`circuit.md` §2.4, Corollary 4.5), which exists
  at every N >= 16 (`circuit.md` Lemma 2.3).

## 4. Actual-source local excess and paired endpoints

Fix eta=17/40000 and L=77/5000 for the heavy and superheavy cutoffs. For the following exact local inequality put e_alpha=alpha^(-2)*L_alpha*c_alpha. It asserts D_i<=alpha^(-2)(D_j+D_k)+e_alpha at every vertex; L_alpha and c_alpha are computed explicitly below.

### 1. Exact three-branch fiber partition

Let pi:T3->H be the universal cover of any finite simple cubic H. Fix
`0<alpha<1/sqrt(2)` and use the complete normalized rows

    f_V(z)=sqrt(c_alpha) alpha^dist(V,z),
    c_alpha=(1-2alpha^2)/(1+alpha^2).

The cubic sphere counts give unit norm. For the three neighbors U_i of V,
partition T3 into `{V}` and the three components B_i of T3-V. On branch
B_i the ratio `f_U_i/f_V` is alpha^(-1); on its complement it is alpha.
This includes the root V, whose ratio is alpha.

On one projection fiber pi^(-1)(p), let

    A_i=sum_(z in B_i, pi(z)=p) f_V(z)^2,
    C=c_alpha*1[p=pi(V)], T=C+A_1+A_2+A_3.

Every T is positive in a connected graph. A zero-mass fiber, if using
the same argument in another setting, simply contributes zero. The
pushed squared rows satisfy

    b_v(p)^2=T,
    b_u_i(p)^2=alpha^2 T+(alpha^(-2)-alpha^2)A_i.

Meanwhile the IDEAL fiber inner product is
`alpha T+(alpha^(-1)-alpha)A_i`. Therefore that fiber's nonnegative
affinity excess is

    d_i(p)=T f_alpha(A_i/T),
    f_alpha(x)=sqrt(alpha^2+(alpha^(-2)-alpha^2)x)
                  -alpha-(alpha^(-1)-alpha)x.          (2)

Summing d_i(p) over p gives exactly
`D_vu_i=<b_v,b_u_i>-<f_V,f_U_i>`. No copied-noise independence or tree
factor-of-IID hypothesis enters this identity. It permits every collision
among lifts in a fiber.

### 2. Scalar inequality and explicit root debit

Write a=alpha and b=alpha^(-1). The function f in (2) is concave,
nonnegative and zero at 0,1. Rationalizing gives

    f(x)=(b-a)^2 x(1-x)
          /[sqrt(a^2+(b^2-a^2)x)+a+(b-a)x].             (3)

The denominator lies in `[2a,2b]`. The numerator is unchanged by
`x -> 1-x`, so `f(x)<= (b/a) f(1-x)=alpha^(-2)f(1-x)`.
This includes the zero endpoints. Concavity and f(0)=0 imply
`f(y)+f(z)>=f(y+z)` whenever y,z>=0 and y+z<=1.

Also f is Lipschitz with constant

    L_alpha=f'(0)=(b-a)^2/(2a),

because its derivative decreases to `f'(1)=-(b-a)^2/(2b)`.
On a fiber, put `x_i=A_i/T` and c=C/T. Their sum is `1-c`, so

    f(x_i) <= alpha^(-2) f(1-x_i)
           <= alpha^(-2)[f(x_j+x_k)+L_alpha c]
           <= alpha^(-2)[f(x_j)+f(x_k)+L_alpha c].

Multiply by T and sum over fibers. Exactly `sum_p C=c_alpha`, giving
(1). All root, constant and boundary mass in this full-profile identity
has been counted. There is no truncation shell because these are the
complete rows used by the deterministic marking.

Both `alpha^(-2)->2` and `e_alpha->0` as alpha approaches critical.
The latter is a uniform scalar bound, independent of H and graph size.
Actual finite truncation and tie noise are applied AFTER this marking,
as in the source; their separate uniform row-approximation debits are
unchanged.

The factor 2 in the limiting LOCAL inequality cannot be reduced merely
using the three-branch simplex relaxation. In (3),
`f(x)/f(1-x)->b/a` as x decreases to zero. Take
`(x_i,x_j,x_k)=(x,1-x,0)` and zero root mass. This sharpness is about
the local fiber-mass relaxation, not construction of an extremal graph.
### 3. Heavy vertices must carry light excess

Fix the SAME heavy cutoff eta as in the canonical-trigger lemma and
suppose alpha is close enough that `e_alpha<eta`. Define

    tau=alpha^2(eta-e_alpha)>0, theta=tau/eta<1/2,
    E=sum_(edges e) D_e.

Every vertex incident to exactly ONE heavy edge has the other two
edges light. Equation (1) forces their total excess to exceed tau.
Thus this is a feasibility constraint on actual edge-excess weights,
not just the binary heavy marking.

Let s_i be the number of vertices incident to exactly i heavy edges
(i=1,2,3), and let a count the s_1 vertices whose two light edges both
exit S0. Then

    m=(s_1+2s_2+3s_3)/2,
    r<=s_1+a+s_2.

The earlier induced-U calculation gives

    6n-H0 <=6|S0|+4r
           <=10s_1+10s_2+6s_3+4a.                   (4)

The light edges at the a special vertices are boundary edges, so each
is counted once and is incident to no other s_1 vertex. Their total
weight is at least tau*a. The required light-edge sums of the remaining
s_1-a vertices count each other light edge at most twice. Hence

    E_light >= tau*a+tau*(s_1-a)/2=tau*(s_1+a)/2,
    E >= eta*m+tau*(s_1+a)/2.                        (5)

For the paired-edge statement set a=alpha, b=alpha^(-1), s=a+b,
c=a^2+b^2=s^2-2 and
C_alpha=[2s sqrt(c-1)-2c]/(s-2). For alpha in [7/10,1/sqrt(2)),
the claim is that the sum of the four off-edge excesses at the two
endpoints of e is at least C_alpha D_e-2L_alpha c_alpha.
The proof retains both root atoms, as follows.
### 1. The reverse branch, including both root atoms

Fix adjacent lifts V,U in the cubic universal cover and one base fiber.
Let B be the component of `T3-V` containing U. For
`f_V(z)=sqrt(c_alpha)*alpha^dist(V,z)`, set

    T=sum_fiber f_V(z)^2, A=sum_(fiber intersect B) f_V(z)^2,
    x=A/T, h=a^2+(b^2-a^2)x.

Zero T contributes zero and can be omitted. Along B the ratio f_U/f_V
is b and on its complement it is a. Consequently the U-centered fiber
mass is `T'=Th`. The branch from U towards V is exactly `T3\B`, not
that set with either endpoint removed. Its U-centered mass is therefore

    A'=a^2(T-A), x'=A'/T'=a^2(1-x)/h.                 (3)

In particular V's root atom is included in A' and U's is excluded.
This identity holds with both root atoms present. It gives
`h(x')=1/h(x)` and, for the same edge excess,

    d_e=T f(x)=T' f(x'),
    f(x)=sqrt(a^2+(b^2-a^2)x)-a-(b-a)x.              (4)

One can also obtain equality of the two expressions by symmetry of
the fiber affinity and ideal fiber inner product. No independence or
tree limit is used.

At V, write `C_V=c_alpha*1[p=pi(V)]` for its root atom. The two branches
away from U have fractions summing to `1-x-C_V/T`. Concavity of f,
`f(0)=0`, and its L_alpha Lipschitz bound therefore give

    sum_(two off-edge terms at V) d_g
        >= T f(1-x)-L_alpha C_V.                     (5)

The identical argument at U gives `T'f(1-x')-L_alpha C_U`.
There are no additional endpoint debits. Summing over base fibers has
`sum C_V=sum C_U=c_alpha`, regardless of any projection collisions.
It remains to prove the scalar inequality

    f(1-x)+h f(1-x') >= C_alpha f(x).                (6)

### 2. Exact scalar reduction and concavity

Put `r=sqrt(h)` and `y=r+1/r`. Since `a<=r<=b` and `ab=1`,

    2<=y<=s, f(x)/r=1-y/s.

Indeed the linear term in f is `(1+h)/s`. The left side of (6), divided
by r, is

    R(y)-cy/s,
    R(y)=sqrt(c/r^2-1)+sqrt(cr^2-1)
        =sqrt(c^2-(w-1)^2),
    w=sqrt((c+1)^2-cy^2), 1<=w<=c-1.                (7)

All square roots in these equalities are nonnegative. They follow by
squaring the sum defining R once; its radicand is
`cy^2-2c-2+2w=c^2-(w-1)^2`.

Here `5/2<=c<3`. For example alpha^2 is in `[49/100,1/2]`,
so `c<=1/2+100/49<3`, and monotonicity of `z+1/z` on this interval
gives its lower bound 5/2. Differentiating at fixed c gives

    ((w-1)^2)''=2c*((c+1)^2-w^3)/w^3 >0.            (8)

The strict sign follows from `w<=c-1<2` and `(c+1)^2>=49/4>8`.
Thus `c^2-(w-1)^2` is strictly concave and positive. The square root
is increasing and concave, so R is strictly concave on `[2,s]`.
Its endpoint values are

    R(2)=2sqrt(c-1), R(s)=c.

It consequently lies above their chord, which is exactly

    R(y) >= C_alpha+(c-C_alpha)y/s.                  (9)

Subtract `cy/s` to obtain (6), including its zero endpoints. In the
interior `2<y<s` the inequality is strict. Equality at `y=2` means
`h=1` and `x=x'=a^2/(1+a^2)`. The other two branch masses at each
endpoint may each concentrate on one branch in the zero-root fiber
relaxation, so both subadditivity steps can be sharp there. Sharpness
of that relaxation does not assert realization by a finite graph.

For completeness, as `y` increases to s the ratio in (6) tends to c:
`R'(s)=0`, while the denominator is `1-y/s`. At critical alpha the
endpoint ratio is therefore 5/2, strictly above the symmetric minimum.
There is no unexamined endpoint singularity.

The C_alpha and root debit just derived are exact on the complete profiles.
Their limiting coefficient C_*=9sqrt(6)+12sqrt(3)-15sqrt(2)-20 exceeds 8/5:
this is equivalent to 5sqrt(6)>8+3sqrt(2), then 17>12sqrt(2), then 289>288.
Both L_alpha and C_alpha have finite limits while c_alpha tends to zero.
One sufficiently critical alpha therefore pays this coefficient with any
prescribed uniform positive allowance before truncation and graph size.
Fix sigma=8/5. At a vertex with one heavy edge D>eta, its light sum is
at least tau=alpha^2(eta-e_alpha), tending to eta/2. Its heavy edge obeys
D<=2alpha^(-2)eta+e_alpha, tending to 4eta. Thus L>4eta excludes a
superheavy edge at every heavy-degree-one endpoint uniformly.

## 5. The single clean-window Price budget

### 1. A local small-excess wave bound

Write lambda=2sqrt(2), c=3-lambda, and eta=17/40000. As in the source,
f_x is the complete normalized positive alpha-radial cover row,
alpha<1/sqrt(2), and b_x is its fiber-norm pushforward. Define

    D_xy=<b_x,b_y>-<f_x,f_y> >= 0

on an adjacent lifted pair. It depends only on the physical edge.
Put lambda_alpha=alpha^(-1)+2alpha and
R_x=sum_(y~x)b_y-lambda_alpha*b_x. The exact tree identities and fiber
triangle/Cauchy--Schwarz inequalities give, coordinatewise,

    -q_alpha*1_(pi(x)) <= R_x <= c_alpha*b_x,
    q_alpha=(alpha^(-1)-alpha)f_x(x) -> 0,
    c_alpha=sqrt(3(alpha^(-2)+2alpha^2))-lambda_alpha -> c.

All errors here are uniform in the graph. The negative part of R_x
has norm at most q_alpha. On its nonnegative coordinates,
R_x^2<=c_alpha*b_x*R_x. Since ||b_x||=1, it follows that

    ||R_x||^2 <= c_alpha*<b_x,R_x>+o(1).

Moreover

    <b_x,R_x> = sum_(y~x)D_xy+3r_1(alpha)-lambda_alpha.

The final two terms tend to zero. Consequently if all three incident
excesses are at most L, the actual base Gaussian wave error satisfies

    Var(sum_(y~x)Y_y-lambda*Y_x) <= 3cL+o(1).           (1)

Replacing lambda_alpha by lambda, truncating at a sufficiently large
fixed radius, normalizing and adding sufficiently small positive tie
noise changes (1) only by another uniform o(1). This uses the same
fixed-alpha geometric-tail bound and fiber-norm contraction as the
source; it does not presume a graph-specific rate or an infinite graph
Gaussian limit.

A clean central-edge window comprises its central edge and the four,
eight and sixteen outward edge positions: 29 in total. They are the full
stars at the two endpoints, their four outward neighbors and the eight
next vertices, exactly all defect centers over the sixteen private-branch
label choices. The three sibling-shared reconstruction errors are
R_i/2, (lambda/4)R_i+R_a/2,
(3/4)R_i+(lambda/4)R_a+R_z/2, on each side. Their global deviation
bounds .086,.208,.336 are linear in the common defect deviation bound c.
If all window excesses are <=L, the proved variance <=3cL+o(1)
therefore divides ALL three error bounds by any d with 3d^2L<c.
Conditioning the JOINT Gaussian on the actual mixed pair at zero and
applying the centered quadratic tilt only decrease each error variance.
No independence of the errors or base-field-only conditioning is used.
The contrast variances, determinant factor (1-s)^3 and pair-density
floor .47 remain unchanged. Compactness on the same eighteen closed
cells transfers these strict limiting bounds uniformly: the reference
pair stays invertible and every tested score has positive conditional
variance from the independent reference component. This is the same
finite-source passage proved in certificate.md Section 9.

For general d>=1, the same proof permits all three error deviations to be
divided by d when 3d^2 L<c=3-2sqrt(2). Any ideal margins m'_k>a=1/10000
are allowed provided the actual tail arguments (m'_k-a)/(sigma_k/d) are
recomputed. Take d=77/40, L=77/5000, and exactly half each original margin
in certificate.md Section 9. These are not coupled margins and their
standardized tails are not invariant. The standalone certify.py recomputes
all 54 tails, 18 cells and 36 half-grid starts, giving
K1_old=1555486141231646/10^18. The baseline is
K0_old=672905379119700/10^18. Exact strict gaps are
6/35-3d^2L=12807/56000000>0, c>6/35, and L-4eta=137/10000>0.

For m heavy edges and t edges with D>L, a central edge is clean if its
29-position witness window has all D<=L. Reversing a nonbacktracking path
of length zero through three shows that one superheavy physical edge
contaminates at most 1+4+8+16=29 centers, even with projected repetitions.
Thus at least (m-29t)_+ heavy edges are clean. The ORIGINAL field Price sum
pays K0_old*E+(K1_old-K0_old)*sum_cleanheavy D, E=sum_all_edges D.
Only the clean increment is added; overlapping witness windows do not
prevent distinct central derivative terms adding linearly.

## 6. Larger closed even paths supply additional curvature


For THIS larger-event calculation fix epsilon=1/10000 and set
x=sqrt(s/(1-s)), s=x^2/(1+x^2), 0<=x<=1/16 (equivalently 0<=s<=1/257).
At physical depth d=1,...,5 require ideal margin
m_d(x)=epsilon+(d/3)x. These inputs are fixed throughout the analytic
integration and are unrelated to the eighteen-cell margins used for W22.
The positive epsilon leaves the actual reconstructed sign margin
sqrt(1-s)*epsilon after the source perturbation. All references to the
fixed margin schedule in this section mean this displayed schedule.

Write `lambda=2sqrt(2)`, `sigma=1/sqrt(12)`, and
`rho=2sqrt(2)/3`. X is the ideal critical cubic-tree wave, with covariance
`w_d=(1+d/3)2^(-d/2)`. Condition its adjacent root pair P=(u,v) to zero
using the canonical Gaussian conditional law. Its two half fields are
independent. At a vertex with parent/vertex values a,b its children are

    (lambda*b-a)/2 +/- sigma*g,    g~N(0,1).              (3)

The contrasts g at distinct vertices are independent and independent of
X_P. This is a fact about this critical law, not a Markov assertion for
an arbitrary finite-base Gaussian source.

Y is instead the limiting pullback of the ACTUAL base profile rows. It
need not have independent branch coordinates and repeated base vertices
may coincide. The public row calculation gives, in every local limit,

    sd(R_z)<=c=3-lambda,
    R_z=sum_(w~z)Y_w-lambda*Y_z,
    Var(H_z)<=1,    H_z=sqrt(3)(Y_child1-Y_child2).         (4)

For clarity, the residual estimate retains the finite source law. For a
full normalized positive alpha^distance tree row f_z, alpha increasing
to 1/sqrt(2), push forward each fiber by its Euclidean norm to obtain
b_z. The tree row identity and fiber triangle/Cauchy--Schwarz inequalities
give coordinatewise

    (alpha^-1+2alpha)b_z-o(1)1_projection(z)
      <=sum_(w~z)b_w
      <=sqrt(3(alpha^-2+2alpha^2))b_z.

The lower and upper coefficients tend to lambda and 3, respectively;
the row norms are one. Hence the limiting residual row has norm at most
3-lambda. Sibling covariance is at least w_2=5/6 by profile domination,
which gives the second bound in (4). Normalized finite truncation and
independent tie noise approximate these full rows as in the retained
source. This is the needed source boundary: covariance domination alone
does not imply the residual bound for arbitrary Gaussian fields.

### Exact two-vertex probability needed by the analytic floor

Let p2 denote the selected branch's probability of exactly two vertices. Then p2=(2/pi)[acos(3/sqrt(21))-acos(5sqrt(2)/sqrt(63))]>6/25. The following derivation pays this premise.

### 1. Actual conditional law and the factor of four

Write lambda=2sqrt(2), sigma=1/sqrt(12). The critical covariance is
\(w_d=(1+d/3)2^{-d/2}\). Given an edge's parent and child keys x,y,
the two outward children are

\[
 (\lambda y-x)/2\ \pm\sigma G,\qquad G\sim N(0,1).             \tag{2}
\]

The sibling-difference variables at subsequent vertices are independent
standard normals, independent of the root-edge pair. For completeness,
a swap of the two subtrees makes the corresponding difference odd.
At two different vertices, the deeper swap fixes the earlier difference;
it also fixes the root-edge keys. The relevant covariances vanish, so
joint Gaussianity gives independence. Harmonic sums reconstruct the keys.

The edge covariance block is positive definite since w_1=2sqrt(2)/3<1.
Conditioning at (0,0) means the canonical centered Gaussian conditional
law given by its Schur complement. It is valid at this probability-zero
value and does not require full conditional-law nonsingularity.

At one endpoint u, the two outward keys are ±sigma X, with X standard
normal. The positive selected one is A=sigma|X|. Its children are
sigma(sqrt(2)|X|±Y), with Y independent standard normal. For X>0,
choose C=sigma(sqrt(2)X−Y). If Y>sqrt(2)X then C<0 and the other
child of A is positive. The two outward children of C are

\[
 \sigma(3X/2-\sqrt2Y)\ \pm\sigma Z,                           \tag{3}
\]

with Z independent standard normal. Both are negative iff
|Z|<sqrt(2)Y−3X/2. Thus

\[
 \mathcal C=\{X>0,\ Y>\sqrt2X,
                  |Z|<\sqrt2Y-3X/2\}                         \tag{4}
\]

makes A--C the COMPLETE crossing component after u is removed. No deeper
key is constrained: every edge leaving these two vertices is noncrossing.
Conversely, a two-vertex component requires one negative child of A,
and no positive child of that child. The two choices are disjoint and
have equal probability by Y-sign symmetry. The selected |X| has twice
the standard-normal density on X>0. Consequently

\[
 p_2=4\Pr(\mathcal C).                                        \tag{5}
\]

The factor four comes from these TWO symmetries, not from four outside
branches. The deeper Z may be selected from the chosen child's independent
difference: selection uses X,Y only, so Z remains independent standard.

### 2. Exact solid angle without quadrature

For centered nondegenerate Gaussian forms L1,L2,L3, sign expansion gives

\[
 \Pr(L_1,L_2,L_3>0)=\frac18+
 \frac{\arcsin\rho_{12}+\arcsin\rho_{13}+\arcsin\rho_{23}}{4\pi}.
                                                                  \tag{6}
\]

Indeed odd sign products have zero expectation by global sign reversal.
The pair identity E[sgn(Li)sgn(Lj)]=2arcsin(rho_ij)/pi follows by
measuring angular wedges of a rotationally invariant planar normal.

Apply (6) to L1=X, L2=Y−sqrt(2)X,
L3=sqrt(2)Y−3X/2−Z. These forms are linearly independent, have variances
1,3,21/4, and correlations

\[
 \rho_{12}=-\sqrt{2/3},\qquad
 \rho_{13}=-3/\sqrt{21},\qquad
 \rho_{23}=5\sqrt2/\sqrt{63}.                                 \tag{7}
\]

On L1,L2>0, the bound b=sqrt(2)Y−3X/2>X/2>0 holds. Therefore Z<b
and Z>−b cover the entire first-two-form event; their intersection is
(4). Z-sign symmetry and the bivariate sign formula yield

\[
 \Pr(\mathcal C)
 =2\Pr(L_1,L_2,L_3>0)-\Pr(L_1,L_2>0)
 =\frac{\arcsin\rho_{13}+\arcsin\rho_{23}}{2\pi}.              \tag{8}
\]

The rho_12 terms cancel. Arcsin-to-arccos conversion gives (1) with (5).
Equivalently, (4) is the spherical quadrilateral with extreme rays
(1,sqrt(2),±1/2), (0,1,±sqrt(2)); its solid angle is twice the
arccosine difference in (1), and its Gaussian mass is that angle/(4pi).

### 3. A conservative rational lower bound

Let delta be the positive arccosine difference in (1). Both angles are
in (0,pi/2), and tangent subtraction gives

\[
 \tan\delta=\frac{20-\sqrt{78}}{10\sqrt3+2\sqrt{26}}>\frac25.    \tag{9}
\]

It suffices for the last inequality that
5sqrt(78)+20sqrt(3)+4sqrt(26)<100. Squaring positive rational endpoints
gives sqrt(78)<177/20, sqrt(3)<7/4, sqrt(26)<51/10. Their weighted sum
is 1993/20<100. Since 0<delta<pi/2,

\[
 \delta>\arctan(2/5)
 >\frac25-\frac{(2/5)^3}{3}
 =\frac{142}{375}>\frac{66}{175}>\frac{3\pi}{25}.              \tag{10}
\]

The strict arctangent bound follows by integrating
1/(1+t²)>1−t² for t>0; the last inequality uses pi<22/7.
Equations (5), (8), and (10) prove (1).

### The finite closed-component certificate with arbitrary outside

The assertion labelled (1) in this local argument is -Delta F>=W_L for the two-orientation finite witness W_L defined below. Only this pointwise assertion is needed, not an infinite-branch limit.


Here is the four-corner calculation used in that assertion. Remove u,v
and the two certified closed even paths. Every such path contributes its
fixed minimum cover cost independent of u,v, since its attaching endpoint
can be covered in a minimum cover at no extra cost. Let C be the cover
cost of the remaining outside graph, and a,b in {0,1} the added costs of
the two remaining single-root attachments in the equal-color corners.
After subtracting the constant closed-path costs, these corners cost C+a
and C+b. One unequal corner has only the central crossing edge, so costs
C+1. The other has both outside attachments and the central edge, so
costs C+max(1,a+b): it costs at least one, and if a=b=1 a cover costing
C+1 would have exactly one central vertex and an optimal outside cover
covering the other attachment, contradicting its single-root cost1.
If either single-root cost is zero, include the opposite central vertex
and that optimal outside cover to attain C+1. If both are one, the two
central vertices attain C+2. Thus the mixed difference is
a+b-1-max(1,a+b)=max(a,b)-2<=-1. If both orientations are certified, all
four attaching paths are closed and even, so a=b=0 and the difference
is -2. This calculation is uniform over every connection in the outside.
### 1. The new arbitrary-outside finite certificate

In a finite simple cubic graph fix adjacent u,v and all other bits.
Assume the outward tree window through depth L+1 is faithful, with
all its displayed incidences and physical labels distinct. Delete
u,v when forming outside crossing components.

For endpoint z in {u,v} and color i in {0,1}, let Q_(z,i,L) mean:

* the two outward neighbors of z have opposite colors;
* its color-i neighbor lies in a COMPLETE endpoint-rooted crossing path
  of even size at most L;
* that whole component and all its incident boundary edges are seen
  in the faithful window, with every edge leaving the path noncrossing.

This is a finite sign predicate. It does not infer closure from an
unseen boundary or from a finite-depth zero-activity message choice.
Put

    W_L=1[Q_(u,0,L) Q_(v,1,L)]
         +1[Q_(u,1,L) Q_(v,0,L)].                       (2)

Then (1)'s pointwise inequality holds for the FULL crossing-cover
objective, regardless of every outside bit and connection.

Indeed an even endpoint-rooted path has a minimum cover containing its
attaching root at no extra cost. A certified closed path is an entire
outside component and is disjoint from all other attaching roots.
Remove the two components certified in one summand of (2); each leaves
only its constant optimum cost. The remaining two attaching roots can
interact arbitrarily through the remaining outside graph. The explicit
four-corner optimization proved immediately above gives
`Delta F<=-1`. If BOTH summands occur, all four disjoint certified
components have free covered roots. The equal-color corners have no
extra cost, while each unequal-color corner has exactly the one unit
for uv. Thus `Delta F=-2`. Elsewhere global submodularity gives
`Delta F<=0`. This proves (1) without assuming independent outside
components that were not certified closed.

For L=2, W_L is the familiar weighted K2 certificate. For L>=4 it is
a strict extension, not another calculation of the already known
four-corner identity.
### 2. Positive arm integrals, including all finite cones

Let `phi_sigma` be the centered N(0,sigma^2) density and define

    H_sigma(t)=0 if t<=0,
               2 Phi(t/sigma)-1 if t>0.

Select the POSITIVE outward root at a zero endpoint, then consider its
crossing component after deleting the central pair. Let p2(m1,m2,m3)
be the probability that this component is a closed two-vertex path with
all its prescribed physical signs at depth d having magnitude >=m_d.
Let p4(m1,...,m5) be the analogous closed four-vertex path probability.
The sign-prescribed sibling at each internal vertex has the same sign
as that vertex. All boundary edges of the chosen component are required
to be noncrossing; no unseen key is constrained.

Use alternating magnitudes a,b,c,d>0 along the four-vertex path, so its
keys are a,-b,c,-d. The other keys forced by (3) are

    -a at depth 1;
    lambda*a+b at depth 2;
    -(lambda*b+a+c) at depth 3;
    lambda*c+b+d at depth 4;
    -(sqrt(2)d+c/2) +/- sigma*g at depth 5.

The side key at each of depths 2--4 automatically exceeds the selected
key's margin at that depth. Therefore EXACTLY

    p4(m)=integral_(a>=m1,b>=m2,c>=m3,d>=m4)
      [2 phi_sigma(a)]
      [2 phi_sigma(b+sqrt(2)a)]
      [2 phi_sigma(c+sqrt(2)b+a/2)]
      [2 phi_sigma(d+sqrt(2)c+b/2)]
      H_sigma(sqrt(2)d+c/2-m5) da db dc dd.              (5)

Similarly,

    p2(m)=integral_(a>=m1,b>=m2)
      [2 phi_sigma(a)] [2 phi_sigma(b+sqrt(2)a)]
      H_sigma(sqrt(2)b+a/2-m3) da db.                    (6)

The first factor 2 chooses the actual sign-selected root. Each later
factor 2 chooses which child continues the unique crossing path. Thus
(5) is the sum of SIXTEEN fixed long-arm cones, (6) of FOUR short-arm
cones. A long arm uses five independent standard contrasts and a short
arm three. Before integrating the final contrast, each fixed event is
a strict finite intersection of linear halfspaces in these drivers.
The two halves are independent under the critical pair-zero law.
Equations (5)--(6) give effective positive integrals, not a numerical
approximation or an unconditioned tree substitution for the base source.

### 3. A useful analytic lower mass

All inequalities in this section hold for nonnegative margins whenever
their displayed lower factors are nonnegative. The exponent in (5),
with sigma^2 removed, is

    Q=a^2+(b+sqrt(2)a)^2+(c+sqrt(2)b+a/2)^2
                +(d+sqrt(2)c+b/2)^2.

Apply `2ij<=i^2+j^2` to each nonnegative cross term. Then

    Q<=D_a a^2+D_b b^2+D_c c^2+D_d d^2,
    D=(15/4+3sqrt(2)/2, 15/4+3sqrt(2),
       7/2+5sqrt(2)/2, 3/2+sqrt(2)).                    (7)

Using `7/5<sqrt(2)<10/7`, `25/8<pi<22/7`, and `sqrt(12)<7/2`,

    sqrt(D_a)<=5/2, sqrt(D_b)<=3,
    sqrt(D_c)<=8/3, sqrt(D_d)<=7/4,
    sqrt(D_a D_b D_c D_d)<32.

Diagonalizing only in this UPPER bound for Q makes the lower density a
product of half-normal densities with scales sigma/sqrt(D_i), times
`1/sqrt(product D_i)`. Drop the positive c/2 from the terminal H.
For D a half-normal with scale sigma/sqrt(D_d), the zero-margin terminal
factor has the exact angular identity

    E H_sigma(sqrt(2)D)
       =(2/pi) atan(sqrt(2)/sqrt(D_d))>14/33.            (8)

Indeed the ratio is >4/5, `atan(4/5)>2/3`, and `2/pi>7/11`.
The Lipschitz constant of H_sigma is <14/5. Also
`H_sigma(sqrt(2)d)<=4d` and the half-normal D density is <49/10.
Consequently removing d<m4 costs at most `(49/5)m4^2`, and shifting the
terminal threshold by m5 costs at most `(14/5)m5`. The first three
half-normal tail probabilities are at least `1-sqrt(2/pi)m_i sqrt(D_i)/sigma`.
This proves the explicit lower polynomial

    p4(m)>=P4(m):=(1/32)
      (1-7m1)(1-(42/5)m2)(1-(112/15)m3)
      [14/33-(49/5)m4^2-(14/5)m5].                      (9)

For example, at ZERO margins this already gives `p4(0)>7/528`, much
larger than the earlier intentionally coarse driver-box lower bound.

For the short arm it is better to retain its known exact zero-margin
mass than to diagonalize away correlations. The three margin losses in
(6) are bounded as follows. Removing a<m1 costs <(14/5)m1. The b-density
before the terminal factor is bounded by its value at zero,

    integral_(a>=0)4 phi_sigma(a)phi_sigma(b+sqrt(2)a)da
       <=sqrt(8/pi)<8/5.

Thus removing b<m2 costs <(8/5)m2. Finally the total a,b density before
the terminal factor, the probability that the branch has at least two
vertices, is

    (2/pi)atan(1/sqrt(2))<2/5.

For the last inequality use `1/sqrt(2)<5/7`, `atan(5/7)<5/8` and
`pi>25/8`. Multiplying this mass by the Lipschitz shift `(14/5)m3`
pays the remaining margin loss. The retained p2(0)>6/25 therefore gives

    p2(m)>=L2(m):=6/25-(14/5)m1-(8/5)m2-(28/25)m3.      (10)

No Gaussian integration is required for (7)--(10). The arctangent and
pi inequalities have short rational alternating-series replays.

### 4. Overlap weights and arbitrary outside persistence

Use the original four-corner mixed difference Delta F of the minimum
crossing-cover objective, with the physical central vertices u,v varied
and all other bits fixed. On the faithful outward depth-five tree window
define W22 as the old two-orientation size-two/size-two witness count,
with ANY prescribed public margins. Define W42 to count the two opposite
orientations where u's attaching branch is a CLOSED four-vertex path and
v's is a CLOSED two-vertex path. Define W24 by exchanging the endpoints.
These are actual finite-graph sign cylinders, with no magnitude margin
required. A superscript m below denotes their additional ideal margin
restriction.

For ANY graph beyond this window and ANY assignment of its bits,

    -Delta F>=W22+W42+W24.                              (11)

For each orientation the size pairs (2,2),(4,2),(2,4) are mutually
exclusive. Every attaching root of a closed even path can be covered in
a minimum cover without increasing its cost. Removing those root-edge
components leaves a common outside optimization. One counted orientation
then gives Delta<=-1; if both orientations count, all four attaching
roots have free minimum covers and the central edge gives Delta=-2.
This is the frozen arbitrary-outside cover-removal argument. It does
not require exterior components to be disjoint, finite, or tree-like.

For W42 there are 16*4=64 fixed choices of physical labels. Under the
ideal critical zero-pair law their probabilities sum to p4*p2 in one
orientation. Each successful actual sign cylinder uniquely specifies its
labels: both root siblings have opposite signs and at every internal
step exactly one child is opposite. For fixed labels the sum W of the
two orientations is a nonnegative EVEN function of EIGHT independent
standard ideal contrasts g. Consequently the ideal total over labels
and both endpoint orders is

    E(W42^m+W24^m)=4 p4(m1,...,m5)p2(m1,m2,m3).         (12)

The factor four is two orientations times two endpoint orders, not a
double count of a single curvature unit. If two orientations overlap,
(11) supplies both units. This also explains why the new events may be
added to the existing public W22 with its own, different margin schedule.

### 5. Eight-contrast reconstruction and its five error centers

Set `Z=sqrt(1-s)X+sqrt(s)Y`, with X,Y independent, and condition Z_P=0.
Then condition on the entire Y. The selected ideal contrasts g retain
their independent standard law: they are independent of X_P and Y.
Write r=s/(1-s) and let H be the eight corresponding Y sibling contrasts
from (4). Reconstruct an ideal patch with zero endpoints and contrast
vector `g+sqrt(r)H`, and multiply the reconstructed keys by sqrt(1-s).

The difference between actual and reconstructed Z at each sibling pair
is sqrt(s) times a common error center. With ancestors of the long arm
u,A,C,E,G, the depth-one through depth-five centers are

    e1=R_u/2,
    e2=(lambda/4)R_u+R_A/2,
    e3=(3/4)R_u+(lambda/4)R_A+R_C/2,
    e4=(lambda/4)R_u+(3/4)R_A+(lambda/4)R_C+R_E/2,
    e5=(5/8)R_u+(lambda/4)R_A+(3/4)R_C
                         +(lambda/4)R_E+R_G/2.          (13)

They follow from `e_child=(lambda e_parent-e_grandparent+R_parent)/2`
with both central errors zero. The short arm uses the first three
centers with its own ancestors. Thus only EIGHT absolute errors need be
controlled: two at each depth 1,2,3, and one at each depth 4,5. The
triangle inequality and (4) give the following exact and rational upper
standard deviations, also after Gaussian conditioning:

    depth  exact upper                         sigma_d
      1    3/2-sqrt(2)                            .086
      2    (sqrt(2)-1)/2                          .208
      3    7/4-sqrt(2)                            .336
      4    sqrt(2)/2-1/4                          .458
      5    13/8-3sqrt(2)/4                        .565.  (14)

For any nonnegative even W and standard Gaussian vector g,

    E W(g+h)>=exp(-||h||^2/2) E W(g).

This follows by translating its density and averaging the two signs of
the linear exponent. A margin-m reconstructed witness survives whenever
each |e_d|<(m_d-epsilon)/sqrt(r). Its actual key margins are then at least
sqrt(1-s)epsilon. Averaging over conditional Y bounds its success by
`E W(g)` times the mass of the small-error event under the tilt
`exp(-r||H||^2/2)`.

This tilt remains a centered Gaussian and decreases all error variances.
Its normalizing mass is
`det(I+r Cov(H|Z_P=0))^(-1/2)>=(1+r)^(-4)=(1-s)^4`, by Hadamard and
the eight variance bounds Var(H_i)<=1. No independence of Y's contrasts
or its residuals is used. A union bound for the eight absolute errors
therefore gives the source-bound survival factor

    T_s(m)>=(1-s)^4 [1
       -4 sum_(d=1)^3 Phi(-(m_d-epsilon)/(sigma_d sqrt(r)))
       -2 sum_(d=4)^5 Phi(-(m_d-epsilon)/(sigma_d sqrt(r)))]_+.
                                                               (15)

Under the fixed margin schedule of Section 6, every ratio in (15) is `d/(3sigma_d)>14/5`.
Mills' inequality, `1/sqrt(2pi)<2/5`, and the positive Taylor bound
`exp(98/25)>50` give `Phi(-14/5)<1/350`. There are sixteen one-sided
tails, hence UNIFORMLY on the entire interval 0<=x<=1/16,

    T_s(m)>=(167/175)(1-s)^4.                           (16)

At s=0 use its limit; no division by zero is necessary. Combining
(9)--(12) and (16) yields

    E[-Delta F | Z_P=0]
      >= the public W22 lower bound
         +4 P4(m(x)) L2(m(x)) (167/175)(1-s)^4.          (17)

All factors P4,L2 are positive throughout 0<=x<=1/16. Their positive
factorizations make this a direct analytic assertion; it is not inferred
from a sampling grid.
### 6. Exact integral, finite source transfer, and consumer boundary

At an adjacent pair the raw excess D_e is indexed by an UNORDERED physical
edge, in the same Price normalization as public Section 9. Its zero-pair
density in the limiting interpolation is at least
`3/(2pi)>.47`, since its correlation is at least rho. At sufficiently
accurate finite reference it still exceeds .47 uniformly on this compact
s interval. Thus the additional credit per unit D_e is at least

    .47*(167/175) integral_(0<=x<=1/16)
          8x P4(m(x))L2(m(x))/(1+x^2)^6 dx.             (18)

Here `ds=2x/(1+x^2)^2 dx` and the tilt is `(1+x^2)^-4`.
Bound the denominator from above by `(257/256)^6`. The numerator is a
degree-seven polynomial with rational coefficients, whose integral is
exact. The replay returns

    integral_0^(1/16) 8x P4(m(x))L2(m(x)) dx
      = 14189070234185433397678189
        /1856250000000000000000000000000,

    K_extra_lower
      =33365691995445067173331275188535296
        /9962431135604594950564205646514892578125
      >3349151/10^12>3/10^6.                           (19)

The first integral is recorded exactly by the replay; the second fraction
and strict comparisons are the paid scalar discriminator. There is no
numerical integration error or omitted Gaussian tail in (18)--(19):
the deliberately coarse analytic lower densities and tails were already
paid in (9),(10),(16), and polynomial integration is exact.

For completeness, transfer to the actual finite base/cover pairs is the
same finite-patch compactness argument as public Section 9, now on the
fixed depth-five patch. Otherwise take a violating sequence with ideal
reference patch covariances tending to w_d, convergent local base-pullback
covariances, and s tending to a point in [0,1/257]. The limiting Y obeys
(4). The conditioned pair block is uniformly invertible, since 1-s is
bounded away from zero and the reference pair block tends to the positive
definite critical block. Each exterior score retains positive conditional
variance from the independent critical reference residual. Hence strict
finite sign-cylinder probabilities converge. Labels and the fixed margins m_d(x)
vary continuously in s. The positive actual safety margin
sqrt(1-s)epsilon and the preceding limiting bound contradict the
violation. This transfer is uniform on this compact interval up to any
chosen fixed approximation loss. It uses the actual finite profile and
base covariance; no base Gaussian branching law is assumed.

The outside graph remains arbitrary, by (11), throughout this argument.
One then uses the same finite Price interpolation and exact cover lifting
as the retained consumer. In particular K_extra multiplies the same RAW
unordered adjacent excess as K_public; a subsequent clean/raw conversion
or good-center tradeoff is a separate consumer calculation. The fixed
slack between (19) and `3e-6` absorbs all finite approximation losses
without requiring a new profile radius in this lemma.

For each orientation the closed size pairs (2,2), (4,2), (2,4) are mutually
exclusive. If both orientations occur, the mixed difference is -2. Hence
the additional k=3/10^6 can be added both to K0_old and K1_old without
reusing a curvature unit: set K0=K0_old+k and K1=K1_old+k. Their clean
increment stays unchanged. This global k needs no larger CLEAN window.

### From the local mixed differences to the physical Price budget

On a degree-k sufficiently high-girth finite cover let X be its finite
reference field and Y the repeated base field, independent, with unit
covariances A,C and C>=A entrywise by fiber Cauchy–Schwarz. With
Z_s=sqrt(1-s)X+sqrt(s)Y, smoothing the threshold indicators in cutwidth
Lemma1.1, applying Gaussian integration by parts, and integrating the two
derivative coordinates across zero gives the unordered-pair identity

 d/ds E F({Z_s<0})=sum_(i<j)(C_ij-A_ij) phi_(ij,s)(0,0)
                         E[Delta_ij F | Z_(s,i)=Z_(s,j)=0].

The two ordered off-diagonal terms cancel the factor1/2 in that lemma.
This identity holds for s<1, where independent reference tie noise makes
the covariance positive definite. Integrate first on compact subintervals;
bounded convergence of F passes to the endpoint. Every mixed difference
is nonpositive. Retain only unordered adjacent pairs on the short compact
intervals used by certificate.md and by the larger-event proof; discard
the remaining nonpositive contributions. On each retained interval the
conditioning pair is uniformly nonsingular and the density bounded.

Each physical adjacent edge has k lifts. Its finite base-minus-cover
covariance differs from the deterministic complete-alpha D_e by an
arbitrarily small uniform error after fixed truncation/jitter. Local
curvature integrals similarly approach their stated floors uniformly,
including clean classifications from that SAME complete-alpha source.
There are 3kn/2 adjacent pairs, so their total error is at most a fixed
constant times kn times the chosen accuracy. Exact repeated-bit cover
lifting (cutwidth Lemma4.2) and division by k therefore give

 E F_base({Z<0}) <= (Gamma_J+accuracy)n
                   -K0*E-(K1-K0)*sum_cleanheavy D_e+O(accuracy*n).

The size-pair argument above permits both local estimates in this ONE
sum even on overlapping interpolation intervals. Letting accuracy be
arbitrarily small pays the complete-profile budget. Lemma1.1(b) applies
to the ORIGINAL Gaussian field to give the same upper bound at every
threshold. Only then is the pointwise two-pass debit added; no Gaussian
comparison is applied to the nonlinear final keys.

## 7. Paired population accounting and the fixed parameter


For the following population proof define r=K1/K0, lambda=L/eta,
sigma=8/5 and
rho=min((r+sigma)/28,(lambda+29+27sigma)/794,(4r+1+2sigma)/96,
        (4lambda+143+54sigma)/2744,(2r+1)/38,(2lambda+85)/1102).
These definitions supply the six ratios called (1) below.
Let S_i consist of vertices of heavy degree i, with size s_i. Let
`S0=S_1 union S_2 union S_3`, and let A be those S_1 vertices whose TWO
light edges both go outside S0. Set `a=|A|, b=s_1-a`. As before,

    m=(s_1+2s_2+3s_3)/2,
    6n-H0<=P, P=10s_1+10s_2+6s_3+4a.                 (3)

Let j count heavy edges with both endpoints in S_1, and j_AA count those
with both endpoints in A. These j edges are vertex-disjoint, and

    j >= (s_1-2s_2-3s_3)/2, j>=0,
    j_AA >= j-b, j_AA>=0.                            (4)

Indeed every S_1 heavy stub not paired with another S_1 consumes a stub
at S_2 or S_3; there are at most `2s_2+3s_3` such stubs. Every one of
the j edges not counted by j_AA uses at least one of the b remaining
S_1 vertices. These facts require no assumption that induced-U two-walks
already provide arbitrary star alternatives.

For a vertex v in S_1, write ell_v for its two incident light excesses.
At finite alpha, the earlier local law gives `ell_v>=tau`, with

    tau=alpha^2(eta-e_alpha),
    e_alpha=alpha^(-2)L_alpha*c_alpha.

For each of the j heavy pairs the new exact law gives

    ell_v+ell_u >= chi,
    chi=C_alpha*eta-2L_alpha*c_alpha.                 (5)

Take alpha sufficiently critical that `chi>=2tau>0`. In the sum of
light-edge excesses, one may charge weight 1 to ell_v at v in A and
weight 1/2 at the other S_1 vertices. Every light edge is charged at
most once: the A edges go outside S0 and cannot meet another S_1;
all remaining charged endpoints have weight 1/2. On a paired heavy edge,
minimizing this weighted sum under (5) and ell_v,ell_u>=tau gives an
additional `(chi-2tau)/2`, or `chi-2tau` if both endpoints lie in A.
Because heavy pairs are vertex-disjoint, these increments add without
reusing a vertex's lower bound. Hence the EXACT finite source inequality is

    E_light >= tau*(s_1+a)/2
                    +(chi-2tau)*(j+j_AA)/2.          (6)

There is no fresh-IID assumption on cover fibers. This is an inequality
for the actual copied-finite-base source and its deterministic D values.

As alpha tends to critical, tau tends to eta/2 and chi to C_*eta.
Consequently, for any fixed `1<sigma<C_*`, uniformly over finite graphs,

    E_light >= eta*(s_1+a)/4
                     +(sigma-1)*eta*(j+j_AA)/2-o(n). (7)

The o(n) here denotes an arbitrarily small prescribed coefficient times
n after choosing one alpha. All counts in (6) are O(n), so the debit is
uniform and precedes finite radius, jitter and graph size.

### 2. Keep the existing two-level budget only once

Let t count superheavy edges D>L. The strict gap L>4eta and the earlier
single-endpoint source law ensure, for sufficiently critical alpha, that
no superheavy edge has an S_1 endpoint. Thus

    t<=s_2+(3/2)s_3.

The original-field Price credit remains exactly

    K0*E+(K1-K0)*sum_(clean heavy e)D_e.

There are at least `(m-29t)_+` clean heavy edges. Heavy and light excess
are disjoint. Combining these facts with (7), normalized by K0*eta,
gives the SINGLE lower bound

    B=m+(lambda-1)t+(s_1+a)/4
                +(r-1)(m-29t)_+
                +(sigma-1)(j+j_AA)/2,               (8)

up to the uniform o(n). No full K1 credit is added a second time to E.

As in the earlier joint proof, replace t by min(t,m/29) if needed.
For t>m/29 this weakens (8) and preserves the superheavy incidence
restriction. Set

    p=m-29t>=0,
    d=s_2+(3/2)s_3-t>=0,
    k=sigma-1>0,
    z=27t+p-2d=(s_1-2s_2-3s_3)/2.

The population and old portion of (8) become

    P=794t+28p-18d-9s_3-4b,
    B_old=(lambda+56)t+(r+1)p-d-b/4.                 (9)

By (4) and monotonicity, (8) is at least

    B_old+(k/2)[z_+ +(z_+-b)_+].                    (10)

This keeps partial availability b>0 rather than replacing all S_1
vertices by the all-boundary case.

### 3. A complete six-ray algebraic certificate

Write rho for (1). First rho>=1/14: each of its six entries is at least
1/14 whenever r>=1, lambda>4, sigma>1. In particular
`18rho-1>=0` and `4rho-1/4>=0`.

We prove (10) is at least rho*P. The term `9rho*s_3` in the difference
is nonnegative and may be dropped throughout.

If z<=0, write `d=(27t+p)/2+f`, f>=0. Then

    B_old=(lambda+85/2)t+(r+1/2)p-f-b/4,
    P=551t+19p-18f-9s_3-4b.

The fifth and sixth entries of (1), together with the two nonnegative
coefficients just noted, prove the claim in this entire region.

Suppose z>=0, so `0<=d<=(27t+p)/2`. For b>=z, the difference of
(10) and rho*P increases with b, with coefficient `4rho-1/4>=0`.
Its minimum occurs at b=z. For `0<=b<=z` the difference is affine in b,
so it suffices to check b=0 and b=z. On either choice it is affine in d,
so it suffices to check d=0 and d=(27t+p)/2. The latter has z=b=0 and
is already the preceding case. At the two remaining corners the exact
coefficient pairs are

    d=0,b=0:
       B=(lambda+29+27sigma)t+(r+sigma)p,
       P=794t+28p-9s_3;

    d=0,b=z=27t+p:
       B=(lambda+143/4+27sigma/2)t
                                  +(r+1/4+sigma/2)p,
       P=686t+24p-9s_3.                             (11)

The first four entries of (1) are exactly the t and p ratios in (11).
All their coefficients in B-rho*P are nonnegative. This exhausts the
regions, including no heavy edges, no superheavy edges, shared light
edges and partial boundary availability. No optimization oracle or
numerical rule grid is used. The six ratios are the six rays checked
by this explicit piecewise-affine certificate, not a claim that every
ray is realizable by a finite graph/source.

Let rho be the minimum of the six ratios displayed at the start of this
population argument, equivalently
min((r+sigma)/28,(lambda+29+27sigma)/794,(4r+1+2sigma)/96,
    (4lambda+143+54sigma)/2744,(2r+1)/38,(2lambda+85)/1102).
Here r=K1/K0, lambda=L/eta and sigma=8/5. Each ratio is >=1/14;
the preceding piecewise affine proof establishes Credit>=K0*eta*rho*P-o(n).
Cutwidth.md supplies H0>=max(0,6n-P) and a TWO-PASS retiming debit
at least (q/5)H0 on |t|<=.67. Therefore the total saving is at least
min(6q/5,6eta*K0*rho)n-o(n). The existing public two-pass rule is used
instead of a packing mask; its weaker trigger floor still exceeds the
binding edge endpoint. This changes no declared circuit coefficient.

For q=361900864/10^15, eta=17/40000, L=77/5000,
K0=675905379119700/10^18 and K1=1558486141231646/10^18,
the unique active ratio is rho=31819/233240. Thus
delta=645198997746292029/2744000000000000000000000,
while 6q/5>delta. The certifier checks every ratio and strict surplus.

For every prescribed error choose ONE alpha sufficiently critical to pay
the paired root debit, sigma<C_*, the superheavy exclusion and all strict
variance gaps. Then fix radius, positive jitter and threshold mesh for all
retained bounds, before graph size. All counts are O(n) uniformly over the
marking induced by this same source. No finite profile is asserted to
attain the limiting delta exactly. Cutwidth.md supplies the uniform
all-threshold and multigraph conclusion at gamma=U-delta.

## 8. All actual integer lengths and the declared coefficient


Use U=26129513284633533/200000000000000000, delta=645198997746292029/2744000000000000000000000, gamma=358496277066174326467971/2744000000000000000000000. Theorem 1 of circuit.md
applies to this graph constant and the single nonnegative nondecreasing
h_gamma=o(D) reconstructed in cutwidth.md. Put c0=1+1/(2gamma) and
s0=ceil(c0*N). It gives at EVERY integer N>=256 the common remainder

 R(N)=5sqrt(N)+(1/2)sqrt(N)log2(N)+2
   +[7sqrt(N)+sqrt(N)log2(N)+6+3log2(2s0)+h_gamma(2s0)]/(2gamma),
 C_B2(L_N)>=c0*N-R(N).

Since s0=O(N), h_gamma(2s0)/N=[h_gamma(2s0)/(2s0)]*(2s0/N) tends to
zero, as do all other terms divided by N. The family and the remainder
are uniform over all circuits and every actual integer length. Exact
arithmetic gives c0=1730496277066174326467971/358496277066174326467971 and
c0-c_*=138374610764266285024267/
35849627706617432646797100000000000000000000000000000000>0.
The strict gap absorbs R(N) at EVERY sufficiently large integer N.
This proves the entire eventual zero-remainder declaration and its
exponent-one finite projection for the same fixed P, hence NP, language.
