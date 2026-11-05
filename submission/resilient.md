# A polynomial-graph resilient disperser from elementary moments

This component is a deterministic monotone function of M bits. We prove that both
output values have positive constant probability even when at most
M^(1-delta) coordinates are arbitrary functions of the others, provided
the other coordinates are sufficiently close to bounded-wise uniform.
We only need dispersal. No claim of negligible bias is made.

This has the same purpose as the resilient output function in the
non-oblivious bit-fixing paradigm, but the construction and proof below
are supplied in full. Its ingredients are polynomial graph intersections,
a direct negative-correlation calculation, Cauchy--Schwarz and truncated
inclusion-exclusion. No correlation inequality for general Boolean
circuits or non-oblivious-source extractor is assumed as a black box.

## 1. Exact construction

Fix delta=1/2^20 and the positive integer a=4/delta=2^22. All constants
are fixed independently of M. At an exceptional small M for which the
following shapes are unavailable, define the output to be zero.
For sufficiently large M put

    L_M=ceil(log_2(M+2)),
    q=the largest power of two <=floor(M/L_M), z=log_2 q,
    w=floor(log_2(floor(q/(8a(z+1))))), U=w*q.

In particular q is asymptotic to M/log M within fixed factors, w is
asymptotic to log_2 q, U<=M, and eventually U>=M/4. Construct F_q by
the deterministic irreducible search proved in seeded.md. Regard the
first U input coordinates as points (x,y), where x ranges over the first
w field elements in their fixed binary order and y ranges over F_q.
Ignore all remaining coordinates. Guard w>=a+1 and w<=q.

Let the row set consist of the L=q^a polynomials f over F_q of degree at
most a with constant term zero. For each row and field element b define
its graph block

    B_(f,b)={(x,f(x)+b): x among the chosen w elements}.

The q blocks in one row are disjoint and cover the U coordinates. Blocks
in two different rows intersect in at most a points: equality is a root
of a nonzero polynomial of degree at most a, because the nonconstant
parts of the two row polynomials differ.

We use only q_0 of the blocks per row. This permits accurate balancing
without relying on a power of two to approximate an exponential mean.
Let T=4z+16 and use the rational approximation

    lambda_z=sum_(j=1)^T 1/(j*2^j).

The series for -log(1-t) at t=1/2 follows by integrating the geometric
series, and gives 0<ln 2-lambda_z<2^-T. Define

    q_0=floor(2^w*(az+1)*lambda_z),

and take the first q_0 offsets b in the fixed field order. Guard
1<=q_0<=q. Eventually q/32<=q_0<=q/8. Each selected graph block is a
tribe: its value is the AND of its w input bits. The row value is the OR
of its q_0 tribes, and the final function F_M is the AND of all L row
values. This completely specifies its value and order at every length.

Since a is fixed, enumerating q^a rows, their q_0 graph blocks and w
coordinates takes M^O(1) time. The field search also takes polynomial
time in q, and the rational approximation uses O(log M) terms with
polynomial bit length. The degree is very large but fixed. No good-set
or coalition information is supplied to the evaluator.

Write mu=q_0*2^-w. Exact rounding and the approximation error show

    mu=ln(2L)+o(1).

Indeed their absolute errors are bounded by 2^-w and (az+1)2^-T,
respectively. Consequently exp(-mu)=(1+o(1))/(2L). The use of the
natural logarithm here is only a probability calculation; all sizes and
the actual computation were specified by integers and rationals.

## 2. Two elementary negative-correlation facts

A finite binary vector is called negatively associated here if, for any
two increasing real functions on disjoint sets of its coordinates, their
covariance is nonpositive. The same property then holds for two decreasing
functions, by negating both. For nonnegative functions with the same
monotonicity on pairwise disjoint sets, their product expectation is at
most the product of expectations, by induction: the product of previous
nonnegative functions has the same monotonicity on their union.

Consider independent uniform bits on a block of b>=1 coordinates,
conditioned not to be all one. This distribution is negatively associated.
For increasing functions g,h on disjoint subsets, let g_1,h_1 be their
values when all their inputs are one, and let G,H be their unconditional
expectations. Set c=2^-b. Before conditioning their expectation product
is GH, by independence. After removing the all-one outcome, direct
substitution gives the covariance

    -c*(g_1-G)*(h_1-H)/(1-c)^2 <=0.

Both differences are nonnegative. This proves the claim for arbitrary
increasing functions, not just for individual coordinates. The same
calculation applies to any block size left after other bits are fixed.

Independent products of negatively associated vectors are negatively
associated. For two such vectors X,Y, condition on Y, apply the property
to the disjoint X coordinate sets of g,h, then average. Their conditional
expectations are increasing functions of disjoint Y coordinate sets;
apply the property for Y. Iterating proves the general product claim.
Unconditioned independent bits are a special case. These observations
will control the joint failure probability of two rows without invoking
a general-purpose correlation theorem.

## 3. Uniform good bits with a non-oblivious coalition

Let B be any fixed set of at most Q=M^(1-delta) bad positions. All other
positions are independent uniform bits for this section. The bad values
may be arbitrary functions of all those good bits. Since F_M is monotone,
its value for every sample lies between the values obtained by fixing
all bad bits to zero and all bad bits to one. Denote these two functions
of the good bits by F_0 and F_1. It is enough to show

    Pr[F_0=1]>=c and Pr[F_1=0]>=c

for one fixed c>0 eventually, uniformly over the positions in B.
Because M is bounded above and below by q times a fixed power of log q,
we may use Q<=q^(1-delta)*poly(log q) in all following estimates. All
constants and powers in those polylogarithmic factors are fixed.

### 3.1 Fixing the coalition to zero

A bad-zero bit disables a tribe containing it. In a single row at most
Q of its q_0 disjoint tribes are disabled. The remaining tribes are
independent and each is one with probability 2^-w. The probability that
this row fails is therefore at most

    (1-2^-w)^(q_0-Q) = exp(-mu+o(1)).

Here Q*2^-w=q^(-delta)*poly(log q)=o(1), and q_0*2^(-2w)=o(1).
Thus each row fails with probability at most (1+o(1))/(2L), uniformly
in B. A union bound over the L rows gives

    Pr[F_0=1]>=1/2-o(1).

No independence between different rows is required for this direction.

### 3.2 Few rows contain a heavily controlled tribe

Put s=ceil(delta*w/2). Call a row heavy if one of its selected tribes
contains at least s bad positions. Eventually s>=a+1. Any a+1 distinct
points on a graph have distinct x coordinates and determine at most one
polynomial of degree at most a, hence at most one global graph block
including its offset. The collections of (a+1)-subsets inside different
graph blocks are consequently disjoint. Counting such subsets of B gives

    number of heavy rows
    <=number of heavily intersected graph blocks
    <=binom(Q,a+1)/binom(s,a+1).

Relative to L=q^a this is at most

    q^(1-delta*(a+1))*poly(log q)=o(1),

because delta*(a+1)>1. Thus L-o(L) rows are light. We may ignore the
heavy rows when proving that SOME row fails with the bad bits fixed to
one.

### 3.3 One light row

In a light row each tribe has fewer than s bad positions. With bad bits
one, a tribe with j bad positions is one exactly when its remaining
w-j good bits are all one, with probability p_e=2^(-w+j). Its good-bit
sets remain disjoint within the row. Its total trigger mean mu_i differs
from mu by at most

    sum_e 2^-w*(2^j-1) <=Q*2^(-w+s)
                           =q^(-delta/2)*poly(log q)=o(1).

The sum of the j values is at most Q and 2^j-1<=j*2^s was used. Also
max p_e<=2^(-w+s), which tends to zero. It follows that the row-failure
probability p_i is

    p_i=prod_e(1-p_e)=exp(-mu_i+o(1))
                       =(1+o(1))/(2L).

The logarithmic approximation error is bounded by a fixed multiple of
sum p_e^2<=max p_e*mu_i=o(1). Every estimate is uniform over light rows
and over the coalition.

### 3.4 Two light rows

Let E_i be failure of light row i. Conditional on E_i under uniform good
bits, its disjoint remaining tribes are independently conditioned not
to be all one. Bits not in a selected tribe of row i remain independent
uniform. Section 2 makes this whole conditional distribution negatively
associated.

Take a remaining tribe A from a different light row j. It has at most
w good bits and intersects each row-i tribe in at most a coordinates.
If a row-i remaining tribe has b good bits and u of them lie in A, then

    Pr[those u bits are one | row-i tribe is not all one]
      =2^-u*(1-2^(-(b-u)))/(1-2^-b)
      >=2^-u*(1-2^(-(b-a))).

Every such b is at least w-s, and at most w row-i tribes intersect A.
Independence of the conditioned row-i blocks therefore shows

    Pr[A is one | E_i] >=Pr[A is one]*(1-tau),
    tau=w*2^a*2^(-w+s)=q^(-1+delta/2)*poly(log q)=o(1).

The row-j failure events for its separate tribes are decreasing functions
on disjoint good-bit sets. Negative association yields

    Pr[E_j | E_i]
      <=prod_(tribes A in row j)(1-Pr[A is one | E_i])
      <=exp(-(1-tau)*mu_j)=exp(-mu_j+o(1)),

because tau*mu_j=o(1). Combining with p_i gives the uniform pair bound

    Pr[E_i and E_j]<=(1+o(1))/(4L^2).

This proof uses only the two-row intersection bound and the explicit
negative-association calculation. No conditional independence of the
entire row values was presumed.

### 3.5 Both forced extremes have positive mass

Let Z count failing light rows with bad bits one. Sections 3.2--3.4 give

    E Z=1/2+o(1), E binom(Z,2)<=1/8+o(1).

For every nonnegative integer Z, 1_(Z>=1)>=Z-binom(Z,2). Thus

    Pr[F_1=0]>=Pr[Z>=1]>=3/8-o(1).

Together with Section 3.1, both outputs of F_M have probability bounded
below, for example by 1/4 at all sufficiently large lengths, regardless
of the bad coordinates' functional dependence on the good ones.

## 4. Bounded-wise uniformity suffices

We now replace uniform good bits by a distribution whose marginal on
any at most t good coordinates has variation distance at most gamma
from uniform. Neither conditional independence on the bad values nor
uniformity of a marginal including bad positions is assumed.

A row failure, with bad bits fixed to one or zero, is the event that
none of its remaining AND triggers occurs. The joint failure of two
rows is likewise the absence of any trigger from those two rows. Each
trigger depends on at most w good bits. For a collection of triggers
write X for their count and let

    S_r=sum_(|J|=r) Pr[all triggers indexed by J occur]=E binom(X,r).

Truncating inclusion-exclusion at an even and the preceding odd order
sandwiches Pr[X=0], with gap S_R at the chosen terminal order R. This
follows pointwise from the binomial identities for a fixed integer X;
then take expectations. The two adjacent orders may be interchanged to
match parity without affecting any estimates.

We need a small UNIFORM remainder S_R. For one row its triggers are
independent. For two rows, decompose X=X_i+X_j. For u>0, Cauchy--Schwarz
and the independence WITHIN each row give

    E(1+u)^X
      <=sqrt(E(1+u)^(2X_i)*E(1+u)^(2X_j))
      <=exp(mu_total*(u+u^2/2)).

Here each row contributes products 1+p_e*((1+u)^2-1); their logarithms
are bounded by the corresponding p_e*((1+u)^2-1). A one-row collection
obeys the same upper bound. This argument requires no independence
between the rows and no general read-k inequality.

Since S_R is the nonnegative coefficient of u^R in E(1+u)^X,

    S_R <=exp(mu_total*(u+u^2/2))/u^R.

Every relevant mu_total is at most 3ln(2L) eventually: this covers one
or two light rows in the one-extreme, as well as zero-extreme rows.
Choose a fixed sufficiently large integer C_R, then
R=ceil(C_R ln(2L)) adjusted to either parity, and u=sqrt(C_R/6).
The logarithm of the displayed bound is at most

    3ln(2L)*(sqrt(C_R/6)+C_R/12)
      -(C_R/2)*ln(2L)*ln(C_R/6)+O(log C_R).

The negative term dominates as C_R increases. Fix, for definiteness,
C_R=2^20. The coefficient is less than -20, so S_R<=L^-10 eventually.
This completely explicit coarse choice suffices. The same conclusion
holds at both consecutive truncation orders.

Each product of at most R triggers depends on at most Rw good bits.
If t>=Rw, the bounded-wise assumption changes its expectation by at
most gamma. There are at most (2q_0)^r products at order r. Thus the
entire truncated probability changes by at most

    (R+1)*(2q_0)^R*gamma.

The uniform Bonferroni gap adds at most L^-10. Therefore, if

    t>=(ceil(2^20 ln(2L))+2)*w,
    (ceil(2^20 ln(2L))+3)*(2q_0)^(ceil(2^20 ln(2L))+2)*gamma<=L^-10,

each individual or pair failure probability used in Section 3 changes
by at most 2L^-10. These are conservative sufficient conditions; the
slight extra two handles both parities.

Summing these errors over at most L individual terms and L^2 pair
terms still gives o(1). Consequently the conclusions

    Pr[F_0=1]>=1/2-o(1), Pr[F_1=0]>=3/8-o(1)

remain valid. Monotonicity once again sandwiches arbitrary non-oblivious
bad functions between those extremes pointwise. Both output values
therefore occur with positive probability on the original source.

In the main application write L_N=floor(log_2 N). The actual seed
list has log_2 M<=2^21*L_N. The exact integer marginal-order and
error-charge bounds in the main proof are at most 2^85*L_N^2 and
2^87*L_N^2 eventually. The supplied h>=2^90*L_N^2 and
E=2^90*L_N^2 pay them with fixed factors. Thus the borderline
logarithmic exponent two is justified by explicit constants, rather
than by claiming that equal powers have an asymptotic strict margin.
All metadata guard branches still define the small-length extension.

## 5. Interpretation of the component

The good set is fixed for each source, and its bits may be correlated
beyond the bounded marginals. The bad coordinates may depend arbitrarily
on those good bits and on additional randomness. The proof did not
condition on their realized values. This is exactly why the monotone
forced extremes were used.

The ignored suffix of M-U coordinates does not change the argument:
projection retains the same bounded-wise marginal guarantees, and the
number of bad used positions is at most Q. The estimates throughout
used M=Theta(q log q), so no equality between the full length and a
special field-compatible length is assumed. The component is one total
uniform polynomial-time family at all lengths, with documented zero
extensions on finitely many exceptional shapes.
