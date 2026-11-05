# Fixed logarithmic-loss family and native full-B2 theorem

This file and its five companion proofs support the exact declaration in
bound.json. Mathematical arguments are supplied below.

## 0. Component map

The component map and provenance are retained in sampler-and-disperser.md.

# A full-B2 circuit-to-affine-flat bridge

All algebra is over F_2. A circuit is a finite acyclic, single-output DAG;
each gate is any one of the sixteen two-input Boolean tables and costs one.
Wires, reuse, depth, sharing and fan-out are unrestricted. Free constants
are allowed, a stronger model for the lower bound.

Attribution: the upstream source named in the supplied proof is GordoAR,
submission 6622943f-7b2d-4f8f-b73c-979c5232484a, version
226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e, proof.md, Chapters II–III;
its public antecedents are 5e9fcf7f-8338-4729-9faa-a18994661af6 and
bf313f35-d97b-4836-b1e8-30afd1157ac2. The supplied proof attributes the
affine-system framework to Mishra–Santhanam, ECCC TR26-260 (2026-10-19).
This reconstruction uses the supplied circuit-reduction.md; upstream sources
were not independently inspected. No upstream theorem is assumed below.

## Lemma

Let d >= 0 be an integer. If f:F_2^n -> F_2 is nonconstant on every affine
subspace of dimension greater than d, then

    C_fullB2(f) >= 4n - 8(d+1).                         (1)

Inputs need not all be essential.

### Representation and potential

During the proof we permit signed references to signals. Signs entering a
gate are absorbed in its truth table. A constant, unary or repeated-parent
gate is replaced by a constant or signed earlier signal. These proof
representations do not assert free native input negation: we start with the
native count and need no conversion of the terminal representation back
to a circuit. The output can be a constant or signed parameter, auxiliary
or gate output.

The six inessential tables are two constants and four unary tables. The
remaining ten are XOR, XNOR and the eight tables

    c + (x+b)(y+e),     b,e,c in F_2.                   (2)

Call the latter product gates. These exhaust the basis.

Use M parameter coordinates u, of which 2q are in disjoint pairs. Each pair
ranges over three of its four corners, and every unpaired (free) coordinate
ranges over F_2. Write D for this product domain. Maintain an affine map
Phi:F_2^M -> F_2^n injective on the entire cube, and a representation of
f(Phi(u)) for every u in D. Equality outside D is not required.

Auxiliaries v satisfy a square system

    A v = N u + c,     A invertible,                   (3)

with at most three distinct nonconstant variables in each row, after
cancellation. A row has an auxiliary, since A is invertible. Remaining
binary gates read parameters, auxiliaries or earlier gates. Gate outputs
never enter (3). Let s be rows plus gates. A parameter's children are the
distinct rows and gates containing it, not pins and not the output reference.
Let I count free parameters having a child. Set

    a=M-q,             V=s+4I+2q.                     (4)

Initially M=n, q=0, Phi is the identity, the system is empty and s=|C|.
Thus a=n and V<=|C|+4n, including circuits omitting inputs.

For each pair, fixing one endpoint opposite its forbidden coordinate
leaves the other endpoint unrestricted. Such safe fixings give a flat
inside D of dimension a; Phi maps it to an a-flat in the original inputs.
Consequently a constant represented output is impossible if a>d. If
a>=d+2, an affine output is also impossible: on a safe a-flat its constant
fiber has dimension at least a-1. A product gate with affine parents cannot
be the output either: on that flat a nonconstant affine factor has a zero
fiber of dimension a-1, and constant factors give a constant output. This covers
constants, parameter/auxiliary literals and terminal gates with affine parents.

Until a<=d+1, lose eight units of V per lost unit of a and never increase
V at fixed a.

### Algebra and cost-free preparations

All coordinate changes below compose Phi with an affine injection. Replacing
a free x by w=x+ell(other parameters), with ell affine, is a bijection of
the whole coordinate cube and of the domain's free fibers; ell may involve
paired coordinates. Fixing a free coordinate lowers a by one. A safe fixing
of a paired endpoint, or restriction of its pair to any two allowed points
parametrized as a line, lowers M,q by one and preserves a. Restricting two
free coordinates to three corners increases q by one and lowers a by one.

Move an affine gate with no product-gate ancestor into (3), processing its
affine predecessors first. Its output becomes a new auxiliary with its
one defining row. The new auxiliary matrix is block triangular with new
diagonal entry one. The row has width at most three; s and parameter use
do not change. Degenerate gates are deleted instead.

The following simplifications decrease s by one, preserve a, and do not
increase I:

* Delete a constant or literal gate, substituting its output. This includes
  constant pins, equal/complementary pins, unary tables, and a gate with
  affine parents whose linear parts are dependent. In that last case a
  parent is constant or the parents are equal up to a bit, so its result is
  a constant or parent literal.
* Delete a gate with no child unless it supplies the output.
* For a row of width at most two, solve for one auxiliary v as a constant
  or signed other variable, substitute, and delete that row and column.
  Other widths do not increase. Add this row to every other row containing
  v, then expand the determinant along v's column: the remaining auxiliary
  matrix is invertible.

A substitution cannot give an unused free parameter its first child: that
parameter already occurred in the object supplying the substituted alias.
Thus the stated use accounting also holds through iterated propagation.

We need one further elimination. An auxiliary Q with coefficient one on a
free x has the form Q=x+ell(other parameters). Exchange their roles: Q is
free and x auxiliary, without expanding any equation. If b=A^(-1)N_x,
then b_Q=1. Replacing A's Q-column with N_x gives determinant
det(A)b_Q=1, preserving invertibility and row widths. Update Phi by
x=Q+ell. Both x before and Q after have a child, so the exchange preserves
I, s, M, q and V.

If free x has exactly one child, a row E, N_x is nonzero. Some Q therefore
has b_Q=1. Exchange x,Q. Now auxiliary x occurs only in E and, during the
loop a>=d+2, cannot be the output by the output exclusion. Delete
E and x, which E determines but nothing else uses. Expanding along
x's one-entry column proves invertibility of the remaining system. The
old free x is replaced by at most one used free Q, so I does not increase.
This deletion loses one row charge and preserves a.

### Clean all pairs

Perform available affine transfers and the preceding eliminations. Then
apply the first applicable pair rule below, propagate constant/literal
substitutions and short rows, and restart. A fixed parameter deletes each
of its old gate children, since a binary table with one fixed pin is unary,
and each old row child, since that row shortens to width at most two.
Count each removed object once. All four rules preserve a.

1. If endpoint p is unused, safely fix its partner h. The new free p is
   unused; q falls by one, I does not increase, and V falls by at least two.
2. If an endpoint has at least two children, safely fix it. At least two
   distinct objects disappear, and only its partner can become newly used
   and free. Hence Delta s<=-2, Delta I<=1, Delta q=-1, so Delta V<=0.
3. Otherwise each endpoint has exactly one child. If partners share that
   child B, its gate output or sole auxiliary is a function of the pair
   (a row B has just the endpoints and that auxiliary). Two of the three allowed points
   give the same bit. Restrict the pair to their affine line. Delete B
   after substituting that constant bit. No other object used the pair,
   so the new line coordinate is unused: Delta V<=-1-2=-3.
4. Otherwise suppose p's sole child is a row E. Its partner h has a distinct
   sole child B. Safely fix h and delete B. If propagation deletes E too,
   two distinct charges have been consumed. If not, now-free p still occurs
   only in E: its only original occurrence could not spread without E's
   deletion. Apply the one-child-row elimination. At most one new used
   free parameter arises, allowing the role exchange. As in rule 2,
   Delta V<=-2+4-2=0.

When these fail, every paired endpoint has exactly one gate child; partners
have different children, and paired parameters do not occur in (3).
In particular every auxiliary depends only on free parameters.

### Prepare a product source

Continue only while a>=d+2. If all gates transfer away, the output is affine,
a contradiction. Choose a product gate G with no product ancestor; its
predecessors are affine and have been transferred. After the simplifications
its parents have nonzero independent linear parts. Turn any auxiliary
parent Q into a free parameter by a role exchange. Such a free pivot exists
because Q is nonconstant and all paired columns of N vanish. If the other
parent is free, independence supplies a pivot different from it; if paired,
Q is independent of that coordinate. With two auxiliary parents, exchange
the first and then use this same argument for the second.

After each exchange, check for simplifications and pair cleanup. If any
is possible, perform a completed action and restart. Otherwise keep this
same G: at most two exchanges make its parents distinct parameters. G is
not the output, by the product-output exclusion, and has a gate child by
dead-gate cleanup. It cannot have a row child.

### Restrictions and the sharing cases

Each restriction lowers a by one. The new safe dimension
is at least d+1, so the restricted output is nonconstant. Freeze old gates
and incidences immediately before the restriction. Propagate topologically,
recording old outputs as constants or signed aliases. Every old gate hit
by a constant pin is removed. If an old output becomes constant, all its
old children are removed too, even if it first became an alias. Charges
refer to distinct removed old objects, not to propagation events.

**Two free parents x,y.** Restrict to the three corners on which G has its
majority value, creating a pair. G and its k>=1 distinct children disappear,
even with shared descendants or children that also read x or y. Both x,y
leave I; no unused old free parameter acquires a first use. Thus

    Delta q=1, Delta M=0, Delta I<=-2, Delta s<=-(1+k),
    Delta a=-1,        Delta V<=-(1+k)-8+2<=-8.        (5)

**Paired p and free x.** Write h for p's partner and B for h's sole gate
child; B!=G. Safely fix h and fix x to kill one factor in (2) at G. G is
constant, B disappears, and newly free p loses its sole use. If G has a
child other than B, it is a third distinct deletion. Otherwise B is G's
sole child and reads the two constants h,G, so B is constant too. It cannot
be the output, since the restricted output is nonconstant. Its old child
exists by cleanup, differs from B by acyclicity and from G because G reads
only parameters, and is a third distinct deletion. Hence

    Delta M=-2, Delta q=-1, Delta I<=-1, Delta s<=-3,
    Delta a=-1,        Delta V<=-3-4-2=-9.            (6)

**Paired p and paired r.** By pair rule 3 their pairs differ. Let h,k be
their partners. Safely fix h,k and then fix
now-free p to kill G's product. G is constant and r is unused. Write B_h,B_k
for the partners' sole gate children, neither equal to G. Four distinct
old gates disappear, as follows.

If B_h!=B_k, already G,B_h,B_k disappear. A child of G outside this set is
a fourth. Otherwise choose a child B_h of G. It reads h,G and is constant.
If it has a child outside the set, that is fourth; otherwise its only
possible child is B_k, which then reads k,B_h and is also constant. None
of these constant gates is the output. Cleanup supplies a further child,
and acyclicity and G's parameter-only inputs put that child outside the
three-gate set.

If B_h=B_k=B, B reads h,k and is constant. G and B read only parameters,
so neither reads the other. Neither is the output; each has a child.
If different children can be chosen, these provide two further deletions.
Otherwise both child sets are the same singleton {J}; J reads G,B and is
constant. J is not the output, so it has a child, distinct from J by
acyclicity and from G,B by their parameter-only inputs. This is fourth.

Only p,r become free; p is then fixed and r unused. Therefore

    Delta M=-3, Delta q=-2, Delta I<=0, Delta s<=-4,
    Delta a=-1,        Delta V<=-4-4=-8.              (7)

These cases exhaust the overlaps of partner children and descendants.
Every charge belongs to one restriction and one distinct old gate.

### Termination and conclusion

Every completed action preserves representation on the current D and
strictly decreases the nonnegative integer M+s: eliminations decrease s,
pair rules decrease M, and (5)–(7) decrease s. Neither M nor s increases.
Preparatory transfers decrease the number of gates; there are finitely
many. A fixed-source preparation uses at most two role exchanges before
an action. Thus preparation cannot cycle. The geometric exclusions prevent
termination inside a>=d+2, so the process reaches a_T<=d+1. Each completed
action satisfies

    V_before - V_after >= 8(a_before-a_after).

As V_T>=0, summing along the one history gives

    |C|+4n >= V_start >= V_T+8(n-a_T) >= 8(n-d-1),

which proves (1). If n<=d+1 the right side of (1) is nonpositive and no
process is needed. Constants, literals, unused inputs and repeated pins are
covered. Each unit of s is an original gate charge, moved or consumed once.

# Exact finite assembly of the two affine breakers

The dependency order is noncircular:

1. Finite seeded-linear extraction and the constant-advice standard breaker.
2. A constant-degree generalized affine breaker.
3. The weak-seed-hashing branch of the equal-length Lemma 5.1 wrapper.
4. Near-full extraction and the general-input standard breaker.
5. An arbitrary-degree generalized affine breaker.

Every constant below is specified by integer/rational operations or a search
whose termination is proved. No search tests an extractor property or a source
distribution. Logarithms are binary. Let lg(x) be the least nonnegative integer
j with 2^j >= x, for positive rational x. All real-looking powers appearing in
error bounds are analytical bounds, not instructions to compute floating-point
numbers. All chosen parameters can be computed with exact integer arithmetic.

## 1. Inputs and their precise scope

All primitive contracts invoked below are proved in the five companion
files. Their exact component map and provenance are in sampler-and-disperser.md.

Advice condensation uses eta=1/10000000 and the supplied ring expansion
gamma=1/10000 above prime threshold2^20000. The CZ finite guards imply
that threshold; substituting this expansion closes its conditional premise.
Canonical prime setup and full ordinary-condenser proofs are in the
companion files, with the component map giving their precise provenance.

The proved outer bound is polynomial in degree and linear in advice,
error exponent, output length and log of ambient source length.

## 2. A reusable exact seeded-extractor envelope

For integer height h>=11, B=2^h, source length u>=1 and error exponent e>=1,
define

    b=e+h+4,
    rho_u=b+2 lg(4u)+2,
    S_u=least power of two >=(4h+1)rho_u,
    G_u=2B(4b+S_u+5),
    P_u=(h+7)S_u+18b+18.

Instantiate the finite recursive primitive for output z>=1 with leaf output
ceil(z/B), then truncate. Its exact threshold K(u,z) and seed length D(u,z)
are computed by its own algorithm and satisfy

    K(u,z)<=2z+G_u,       D(u,z)<=6z/B+P_u.                 (2.1)

The call is legal whenever its computed K(u,z)<=u. The inequality
128 P_u <= G_u holds if B>=64(h+25), since
P_u<=(h+25)(S_u+b+1) and G_u>=2B(S_u+b+1).

A useful envelope, valid with the definitions above, is

    G_u <= C(e+h+lg(u+1)+6),    C=64B(h+1).                (2.2)

Indeed rho_u=e+h+10+2lg u <=2(e+h+lg(u+1)+6),
S_u<2(4h+1)rho_u, and
4b+S_u+5 <= (16h+9)(e+h+lg(u+1)+6).
Both G and P are nondecreasing in u. No unspecified extractor constant
appears in these envelopes.

## 3. A polynomial-degree affine assembly from any linear standard interface

### 3.1 The standard interface required here

Suppose a standard one-tampering breaker C is supplied with an integer c>=1.
On source length u, advice length a, output z and requested exponent e, its
actual seed length d_C and source threshold k_C satisfy

    d_C,k_C <= c[z+a+e+lg(u+1)+1].                          (3.1)

It accepts arbitrary correlated source pairs independent of arbitrary
correlated seed pairs, honest seed uniform, and retains both seeds with
error <=2^-e. Its only remaining domain guard is k_C<=u; any other finite
guard must already be built into the interface. The fixed-advice interface
below is used only at its one fixed advice length. The final standard
interface in Section 7 works at every a>=1.

### 3.2 Exact affine parameters

Given n,t,a,m,E>=1, set

    e=E+lg(26t),                 j=lg t,
    h=20+2 lg(c+t+16),          B=2^h,
    C=64B(h+1),                 R=t+5,
    T_c=(t+1)(6c+2+6j),
    Q=T_c+(t+2)(3c+1)+t+4,
    K_c=T_c+(t+2)R+6c+3.

Here h is a seeded-extractor height, while j is the number of affine-merging
rounds. They are different parameters. Select w as the least power of two
satisfying

    w>=m+a+e+1,
    w>=lg(Qw+1),
    w>=G_U,       U=max(n,Qw),                              (3.2)

where G is computed at the just-selected h,e. Define

    d=Qw,          k=K_c w,
    s0=dx=dy=w,    d0=3cw,    r=Rw.                        (3.3)

Invoke C at source length d, output dx=w, advice a, error 2^-e. Its real
seed and threshold are at most 3cw, because the bracket in (3.1) is at most
w+(a+e+1)+lg(d+1)<=3w. Pad its seed to d0 bits.

Use (2.1), with the common h,e, to instantiate

    L0: input n, output d0, padded seed s0;
    Lr: input n, output r,  padded seed dx;
    E0: input d, output dx, padded seed dy;
    Lm: input r, output dy, padded seed dx.

Every seed here has exactly the dimension required in the affine proof.
Return the first m bits of its r-bit final output. The extension branch for
m>r is never needed because w>=m.

### 3.3 Every finite affine guard

Our choice of h gives B>=64(h+25) and
B>=12 max(3c,R). For example put s=c+t+16>=18. Then
B>=2^20 s^2, h+25=45+2lg s<=47s, and
max(3c,R)<=3s. These imply both inequalities directly.
Consequently (3.2), (2.1) and 128P_U<=G_U give

    D(u,z)<=w/2+w/128<w

for every listed target z in {3cw,Rw,w}, and

    K0<= (6c+1)w,    Kr<=(2R+1)w,
    KE<=3w,         Km<=3w.                               (3.4)

All relevant ambient lengths are <=U: in particular Q>=R and Q>=3c.
The L0,Lr guards follow from k<=n. The E0,Lm guards follow from
d=Qw>=3w and r=Rw>=6w. The standard C guard follows from d>=3cw.
Thus the sole external entropy guard is the explicit integer test k<=n.

The affine transcript budget is exactly

    T=(t+1)(s0+2d0+dx)+j(t+1)(4dy+2dx)=T_c w.

Since e<=w, all sufficient inequalities of the generalized affine proof hold:

    r >= (t+1)dy+Km+e,
    k >= max(K0,T+Kr+t*r+e),
    d >= (t+1)(s0+d0)+k_C+e,
    d >= T+KE+t*dx+e.

For the first line, its right side is at most (t+5)w=r. The second is
bounded by max((6c+1)w,[T_c+(t+2)R+2]w)<=K_c w. The third right side is at
most [(t+1)(1+3c)+3c+1]w, and the fourth at most (T_c+t+4)w; Q dominates
both. Also s0<=d and dy<=r. The resulting error is <=26t*2^-e<=2^-E.
This checks actual dimensions rather than using uncalibrated O notation.

### 3.4 The selector halts, with a polynomial bound

Put

    V=e+h+m+a+lg(n+1)+lg(Q+1)+10.

Let w* be the least power of two >=16 C^2 V. Then
w*<32 C^2 V. As n,Q,w*>=1,

 lg(max(n,Qw*)+1)<=lg(n+1)+lg(Q+1)+lg(w*+1),
 lg(w*+1)<=6+2lg C+lg V<=6+2C+V.

It follows that the bracket in (2.2) is at most 2V+2C+6<=4CV,
so G_max(n,Qw*)<=4C^2V<=w*/4. The other two tests in (3.2) also hold:
the first is bounded by V, and the second by 2V+2C+6<=4CV<=w*.
Thus w<=w*<32C^2V. The selector consists only of finite dyadic trials.

A convenient fully numerical envelope is

    d,k <= A(c)(t+1)^9 [m+a+E+lg(n+1)+1],
    A(c)=2^90(c+17)^9.                                    (3.5)

Here is explicit accounting. Set s=c+t+16 and H=m+a+E+lg(n+1)+1.
Then B<2^22 s^2, h+1<=23s, C<2^33s^3,
Q,K_c<=20s^2, and V<=47sH. The C bound uses64*2^22*23<2^33. Therefore
w<2^77s^7H, d,k<2^82s^9H, and
s<=(c+17)(t+1). The slack in 2^90 covers every endpoint.
The constants do not depend on t; no 2^t occurs. This is the essential
outer linear-logarithm, polynomial-degree bound.

For parameter choices with k>n, define the affine map as the all-zero m-bit
string and make no extraction claim. This is a total computable family.

## 4. Fix the finite condenser constants first

Set

    eta=1/10000000,     theta=eta/100000,
    beta=theta/2,      tau=9/10,     nu=1/10.

Compute the ordinary-condenser constants literally as follows. Use its
scalar-growth eta_s=1/1000000, theta_s=eta_s,
kappa_s=beta*eta_s/8 and P0=ceil(512/(beta*eta_s)). Let J be the least
nonnegative integer satisfying

    (1+theta_s)^J * (beta/4) >= 9/10.

This rational search terminates: (1+theta_s)^J>=1+J theta_s. Put

    ell=8*3^J,
    c0=min(kappa_s,3/10,1/2),
    C0=4(3^J-1)+2,
    gamma_o=min(c0/2,beta/100,(1-tau)/100,1/100),
    m_R=max(16,P0,ceil(2 lg(C0)/c0)),
    a_o=gamma_o/(16ell),    g_o=a_o/4,
    T_o=ceil((tau-beta)/g_o),
    A_o=ell^T_o,            D=(ell+1)^T_o,
    m_*=max(m_R,ceil(4/a_o),ceil(2(lg(ell+1)+1)/a_o),
            ceil(4/gamma_o)),
    q_*=max(ceil((8/gamma_o)lg(2T_o)),ceil((2/nu)lg(2D))),
    R_o=A_o*max(m_*,q_*),
    c_o=min(gamma_o/8,nu/2)/A_o.

For every integer r>=R_o, the fixed ordinary map has exactly D rows,
q=ceil(r/A_o) bits per row, and error <=2^-c_o r to an exact convex mixture
with a fixed 0.8q-entropy row in each component. These assertions are the
finite module's (12)--(14), not a new existential choice.

For the advice condenser use 10 blocks of r bits and advice capacity
b=10s, where s=floor(theta*r). Pad an actual advice of length a<=b with zeros
to b bits before every CZ call. Thus the CZ input length is u=r+b and the
honest coordinate entropy is r-2b. Its deficit is
rho=3b/(r+b)<=3eta/10000<eta/2. If eta*r>=64, the prime guard, the rounded
flat-source entropy guard and the additional 2^20000 prime threshold all
hold. The nonmalleable output length is h_N=floor(eta*(r+b)), and error is
at most 91*2^-eta(r+b), by the stated CZ input.

Set delta=2^-s and B_A=19b=190s in the advice-condenser lemma. If
r>=max(ceil(64/eta),ceil(28/theta)), then s>=4, 2^s>=10, h_N<=r, and

    h_N >=192s = B_A+s+log(1/delta).

Indeed h_N>=eta*r-1, 192s<=192theta*r, and eta*r>=64. Its precise error
therefore obeys

    e_A<=14*2^-s+91*2^(-eta*r+190s)
        <=105*2^(-theta*r/2)<=2^(-theta*r/4).               (4.1)

The last step uses theta*r>=28 and 105<2^7. It supplies good-row entropy
s>=beta*r. These constants are fixed before any breaker parameter is chosen.

## 5. Fix the inner constant-degree affine breaker

Let

    L_rows=11D,     t0=2L_rows-1,     a0=lg L_rows,
    M0=2^a0,
    j0=lg(8M0+2),
    cbase=max(128,3M0+8),
    c_in=2^32 [3M0+8+128(2M0+4)] [j0+lg cbase+2].

Section 8 of `constant-advice-breaker.md` gives interface (3.1) with this
c_in at the fixed advice length a0. To see the error conversion, use its
internal exponent e+j0 for requested exponent e. Its dyadic variable is
less than 2^31[z+e+j0+lg u+lg cbase+1]. Its actual seed is (3M0+8) times
that variable, and its source threshold is 128(2M0+4) times that variable.
The displayed c_in dominates both (3.1) bounds, with every additive
constant absorbed by j0+lg cbase+2. Its source guard is exactly the one
in (3.1). Its exponential advice cost depends only on a0, already fixed
by the condensers, and NEVER on the later arbitrary degree t.

Use Section 3 at c=c_in, degree t0, advice a0. Define the fixed coefficient

    F=A(c_in)(t0+1)^9,
    kappa0=1/(25600 A_o),
    mu=kappa0/(16F).                                      (5.1)

For equal source length d, the inner affine output and requested exponent
will both be

    m=floor(mu*d),       E_F=floor(mu*d).

Its computed seed and source requirements are at most
F[2mu*d+a0+lg(d+1)+1]. If

    F[a0+lg(d+1)+1] <= kappa0*d/8,                         (5.2)

both requirements are at most kappa0*d/4. We will pad its seed to the
available m0-bit string. Padding retains the entire good and tampered seeds;
it is not a new seed-independence hypothesis.

## 6. An exact equal-length Lemma 5.1 theorem

### 6.1 Parameters, error rate, and output rate

For d>=1 use the weak-seed branch exclusively. Set

    h_IP=floor(7d/20),    r=floor(d/100),    v=10r,
    q=ceil(r/A_o),        m0=floor(q/64),
    eH=ceil(7q/32),       BH=lg(d*m0+1),    rH=BH+m0+2eH,
    sY=d/20, sW=d/25, sX=sZ=d/10,
    m=E_F=floor(mu*d).

Use IP output length v, the advice and ordinary condensers from Section 4,
the explicit weak-seed linear hash on (right source Y, row S), and the
inner breaker of Section 5. Row labels are fixed a0-bit encodings.
The two honest d-bit sources need entropy >=9d/10. Their copies are
arbitrarily correlated within their respective independent pairs.

For d>=200, r>=d/200, hence q>=d/(200A_o). If d>=25600A_o, q>=128 and

    m0>=q/128>=kappa0*d.

Directly from q=ceil(r/A_o)<=r+1 and d>=200, q<=d/8. Thus m0<=d/512,
and the primitive entropy guard
m0+2eH<=29q/64+2<=floor(d/10) holds. For d>=100 the last inequality follows
from 29d/512+2<=d/10-1. If

    d/(200A_o) >=64 lg(d+1)+64,                            (6.1)

then q>=32BH+64, because BH<=2lg(d+1). This implies 2rH<=q, exactly the
hash's seed-fit guard. Its error is <=2^(-3q/160).

The residual entropy ledger is, without any rounding loss in its hypotheses,

    KY>=9d/10-2h_IP-d/20-m0-d/25 >= d/10,
    KX>=9d/10-2v-d/10-m-d/10 >=2d/5.

Here m<=d/10, since mu<1/10. By (5.2), the inner affine seed requirement
is <=kappa0*d/4<=m0 and its source requirement is <=kappa0*d/4<=KX. Its
own source guard is therefore met, not merely assumed. The actual seed
is padded to m0 and its output is m bits.

The advice capacity is b=10floor(theta*r)>=5theta*r>=theta*d/40, so

    a <= theta*d/40                                      (6.2)

suffices. At all lengths satisfying the finite guards, the error terms obey

    eIP <=2^(-d/80),
    eA  <=2^(-theta*d/800),
    eC  <=2^(-c_o*d/200),
    eH  <=2^(-3d/(32000A_o)),
    eF  <=2^(-mu*d/2).

For the first, use the integration module's safe IP exponent:
h_IP-d/5-v-1>=d/20-2; d>=80 then gives the claim.
For the last, require mu*d>=2. Tail errors sum to at most 4*2^(-d/25).
Consequently define rational constants

    lambda=min(1/80,theta/800,c_o/200,3/(32000A_o),mu/2,1/25),
    sigma=lambda/2,       omega=mu/2.

The exact Lemma 5.1 ledger has error <=12*2^-lambda*d. If lambda*d>=8,
this is <=2^-sigma*d. Its output floor(mu*d)>=omega*d when mu*d>=2.
The complete original right pair (Y,Y') is retained by that wrapper's
probability/mixture proof. Neither the weak-seed hash alone nor a claim of
joint two-source strongness is being substituted for that proof.

### 6.2 One computable threshold that works for ALL larger integers

We must not test one length and assume its rounded guards persist.
The following elementary threshold lemma handles the two logarithmic tests.
For positive rationals A and nonnegative rational b, find the least integer
j>=2 with

    A(j+2+b) <=2^j.                                      (6.3)

This exact dyadic search terminates, and for every integer d>=2^j,
A[lg(d+1)+b]<=d. To prove the all-length statement, on
2^k<=d<2^(k+1), lg(d+1)<=k+1, and
2^k/(k+1+b) is nondecreasing for k>=2. Equation (6.3) leaves one extra
unit of slack. Exponential growth also proves termination.

Apply this lemma to

    (A,b)=(8F/kappa0,a0+1)       for (5.2),
    (A,b)=(12800A_o,1)          for (6.1).

Let N_log be the maximum of the two resulting powers of two. Define

    R_adv=max(ceil(64/eta),ceil(28/theta)),
    N51=max(200,25600A_o,200R_adv,200R_o,
            ceil(2/mu),ceil(8/lambda),N_log).

For EVERY integer d>=N51, all guards just listed hold. Thus we have obtained
one explicit deterministic equal-length breaker with output floor(mu*d),
error <=2^-sigma*d, and advice range (6.2), under the stated finite inputs.
There is no additional unspecified C_51 or asymptotic minimum length.
For smaller d, or an advice length violating (6.2), define zero output and
assert no extractor guarantee. This is only a totality convention, not a
small-input extraction claim.

## 7. The general-input standard breaker, with a fixed linear coefficient

Let n,a,m,E>=1. Set H=m+a+E+lg(n+1)+1 and choose the fixed integer

    c_std=ceil(2^31+2N51+100/theta+4/sigma+4/mu+100),
    k_seed=c_std*H.                                      (7.1)

The genuine domain guard is k_seed<=n. A source X of n bits and entropy
at least k_seed, and an independent uniform k_seed-bit Y, suffice.
Their tampered pairs are arbitrary. We construct m output bits and error
at most 2^-E, retaining Y,Y'.

Use the near-full finite extractor at entropy k_seed and exponent E+3;
let its computed complete seed length be h_seed and its output length be
L=floor(9k_seed/10). Its internal primitive has e=E+7 and a_near=E+22. At input n its S<2880F_n,
G<2^24F_n, H_near+1<=2^25F_n and P<2^16F_n, where
F_n=E+lg(n+1)+1. These follow by direct substitution in near-full-seeded.md:
rho=E+28+2lg n<=32F_n. Thus k_seed>=2^30F_n implies both exact guards
k_seed>=16(H_near+1) and k_seed>=1000P. Our choice is larger.
The theorem supplies h_seed<k_seed/100 and error <=2^(-E-3).

Take disjoint first h_seed and last L slices of Y. The first slice seeds
extraction of X to A; the last is B. Do the identical operations on X',Y'.
Apply the Section 6 equal-length breaker to (A,B), then truncate to m bits.

Here are every remaining wrapper guard. Since k_seed>=10,
L>=4k_seed/5. Hence L>=N51, a<=theta*L/40,
sigma*L>=E+2, and floor(mu*L)>=m, by (7.1). The slices are disjoint because
h_seed+L<0.91k_seed<k_seed. For the entropy tail tau=2^(-E-2),

    h_seed+E+2 <=L/10,

since h_seed<k_seed/100 and k_seed>=100(E+3). No exposed slice is omitted
from the final guarantee. The finite general-input proof, with
eta=2^(-E-3), tau=2^(-E-2), and epsilon_L=2^(-E-2), gives total error

    2eta+tau+epsilon_L =3*2^(-E-2)<2^-E.

Its conditional-kernel step restores the full original Y,Y' after retaining
the used slices; it permits all correlations within those seed copies.
The actual seed length and source threshold are both exactly k_seed.
Therefore Section 3.1 holds with c=c_std, at ALL a,m,E>=1, with the single
external guard k_seed<=n. This exports an actual finite interface, not an
uncalibrated statement that all its lengths are O(log n).

## 8. Final arbitrary-degree generalized affine breaker

Apply Section 3 using the standard breaker of Section 7, c=c_std. On
arbitrary degree t>=1, advice length a>=1, desired output m>=1, exponent
E>=1 and ambient source length n, run the literal selector (3.2) and the
literal parameter algorithm (3.3). If k<=n, the resulting map satisfies

 TV((F(X0,Y0,alpha0),F(X1,Y1,alpha1),...,F(Xt,Yt,alphat),Y0,...,Yt),
    (U_m,F(X1,Y1,alpha1),...,F(Xt,Yt,alphat),Y0,...,Yt)) <=2^-E,

whenever Xj=Aj+Bj, H_inf(A0)>=k, Y0 is uniform on d bits,
(A0,...,At) is independent of (B0,...,Bt,Y0,...,Yt), and alpha0 differs
from every alphaj. All copies within either tuple can be arbitrarily
correlated, including additional randomness not determined by the honest
copy. No affine-function or deterministic-tampering restriction is added.
The full seed Y0 is retained as well as every tampered seed.

Its finite resource bound is

    d,k <= C_final (t+1)^9 [m+a+E+lg(n+1)+1],
    C_final=2^90(c_std+17)^9.                              (8.1)

C_final is one fixed computable integer depending only on the supplied
finite modules. It is independent of n,t,a,m,E. In particular for m=1,
seed length is linear in log n, advice length and error exponent, with
polynomial dependence on degree. This avoids the two prohibited losses:
there is no (log n)^2 term and no exponential-in-t coefficient.

The sharper normalized F involving t*lg(t)^3 is NOT established here and
must not be quoted as a consequence. Equation (8.1) is the replacement
proved from these particular finite primitives. If a final sampler schedule
uses (F), it must be adjusted to (8.1); exponential output/seed tradeoff in
its fixed height still dominates this fixed polynomial degree.

## 9. Uniformity, exactness, and remaining obligations

All parameter searches are now proved terminating: the two dyadic seeded
selectors, the rational ordinary-growth J search, finite constant formation,
and the log-domination thresholds. Exact seed lengths and entropy thresholds
are computed using the original primitive recurrence, even though the
conservative envelopes decide a legal pad length. The evaluator never
receives entropy witnesses, good-row labels, source decompositions, coupling
coins, proof histories, or a per-length truth table.

For fixed t and outer-logarithmic m,a,E, auxiliary dimensions are O(log n).
The seeded-linear algorithms are polynomial in their dimensions and B;
the finite ordinary/CZ/weak-seed/bilinear setup costs 2^O(local dimension),
so their complete composition is polynomial in the actual outer n. Fixed
constants can be enormous. This does not assert polynomial bit-time jointly
in arbitrary t and n: the exponent of n may depend on t. The intended
sampler eventually fixes t once, after solving its degree feedback.

For ordinary seeded primitives or standard breakers with the complete
source tuple independent of the complete seed tuple, unused suffixes may be
retained by their conditional kernel given all used prefixes. The outputs
depend only on the source tuple and those prefixes. This applies to the
primitive pads in Section 3 and to the zero-additive-shift independent-side
uses of the inner affine map in Sections 5--6. It does not assert a generic
kernel principle when arbitrary additive B-tuples correlate with suffixes.

For that more general affine case, use this separate consequence of the
universal contract just proved. Fix an admissible A tuple and advice. For
honest used seed y and a deterministic finite profile
r=(B0,...,Bt,Y1,...,Yt), let Delta(y,r) be the output-tuple TV from an
independent uniform honest output with unchanged tampered tuple. Choose a
maximizing profile r_y separately for each y. Uniform y with this profile
is an admissible right-side law independent of A. Since the contract retains
y, it gives the uniform-y average of max_r Delta(y,r) at most epsilon.
For any actual complete right-side variable R, including ignored suffixes,
conditioning on (y,R) fixes one profile. Its average conditional output TV
is bounded by the same average of maxima. Thus the whole right side and
all ignored suffixes may be retained even with B-suffix correlation, at
no additional error. The finite maximizing profiles are proof devices only;
the evaluator never searches for them. The proof does not assume uniformity
or independence of tampered suffixes.

# Cubic-certificate sumset disperser: one fixed all-length family

All logarithms are binary; lg(x)=ceil(log2 x) for positive x is computed
by bit length and power-of-two comparison. The finite primitive proofs
and their exact guards are supplied above and in the five companions.
The sumset architecture is due to Chattopadhyay--Liao,
https://arxiv.org/abs/2110.12652v1 , with the ordinary breaker following
Xin Li, https://eccc.weizmann.ac.il/report/2023/023/revision/1/ .
The new outer argument uses a cubic sign certificate. No previous verdict
or accepted theorem is a mathematical premise.

## 1. Finite interfaces and sharper seed padding

(L) The finite source-linear strong extractor in seeded-linear.md has the
following specialization. Fix integer h>=11 and put u=h+1. For actual input
length n, let L=lg(n+1), requiring L>=(8h+2)(h+10). With primitive error
exponent L and leaf output 128uL, its output length and padded seed are

  d=128*2^h*u*L, s=B0(h)L, B0(h)=24u^2+895u.

Indeed a_E=L+h+4 and rho=a_E+2lg(4n)+2<=3L+h+10. The exact
power-of-two selector satisfies S_E<=2(4h+1)rho, so

  S_E/L<=24h+6+(8h+2)(h+10)/L<=24h+7=24u-17.

The general seed bound D<=6ell+(h+7)S_E+18a_E+12 therefore gives

  D/L<=768u+(u+6)(24u-17)+18+(18h+84)/L
     <24u^2+895u-66<B0(h).

Here L>=h+6 makes the last fractional term <18. Also a_E<=2L,
S_E<=24uL, and 4<=4L. Thus
K0<=2*2^h(128uL+4a_E+S_E+4)
  <=2*2^h(152uL+12L)<=306*2^h*uL<3d, since u>=12.
The actual guard K0<=n is required. Pad the seed to s by ignoring its suffix.
The sampler at failure1/100 uses entropy at most3d+7. This extractor
is source-linear and has the supplied uniform polynomial evaluator.

(G) Select the optional outer ordinary constructor with k_f=2048,
eta_f=1/465,sigma=1/1024,rho_f=1/2+1/(2k_f)+eta_f/(8k_f),beta=1/2.
Use R_G=3900702720ell_O,nu=1/R_G,tau=1/2+2/R_G. Its fixed integers
M,D0,t0 give row length v(t)=ceil(t/M) and mixture row entropy
(1/2+1/R_G)v, with error1/2. Put b=(33R_G/64)M and
kappa_G=33(R_G-2)/128. Here b is an integer multiple of1024.
For left t=b*s bits and right s bits, use the graph at z=L under guards

 t>=t0, v(t)>=64, L>=200R_G+2, 100R_G*lg(v(t))<=L.

It applies O directly to its t-bit index, with no padding enumeration.
Thus v=(33R_G/64)s. The walk has m>=197v/(100R_G)>s bits;
its exact-degree field setup and strict support proof are supplied.
Every left set of size at least2^(t/2) has more than2^(1-L)*2^s neighbors.
The degree is at most D0*(1+floor((kappa_G+1)s/L)). For outer s=B0L,
L>=2ceil(kappa_G*B0)+3 gives at most1+ceil(kappa_G*B0) walk positions.
Write Samp(x,alpha,j) for the resulting source-linear sampler.
At bad-set density A^(-1/2), outside exceptional mass1/100, its all-bad
row count is STRICTLY LESS than sqrt(A), by the strict support proof.
Since sqrt(A) is integral, there are at most sqrt(A)-1 bad rows.

(F) The supplied degree-nine contract exports one fixed computable integer
C_F=2^90*(c_std+17)^9 and a total parameter/evaluation algorithm. For
n,t,a,m,E>=1 it returns a generalized affine breaker with d_F,k_F satisfying

  d_F,k_F<=C_F*(t+1)^9*(m+a+E+lg(n+1)+1).

When k_F<=n, for arbitrary correlated left tuples (A0,...,At) independent
of the full right tuple (B0,...,Bt,Y0,...,Yt), Hinf(A0)>=k_F and Y0 uniform,
its honest output on Aj+Bj is 2^-E-close to uniform jointly with all
tampered outputs and every seed, provided honest advice differs from each
tampered advice. The copies need not be deterministic functions of their
honest copy. The supplied proof includes all seed-retention and finite
parameter details. For fixed t and a,m,E=O(log n) evaluation is polynomial
in actual n.

## 2. Bad seeds for source-dependent seed maps

Fix a correlated left tuple (A_0,...,A_r), fixed advice differing from
the honest advice, and deterministic maps L_i(A_i). Their joint range has
at most2^ell values. Assume Hinf(A_0)>=K+ell+lg(1/epsilon).
Call w bad if some b_0,...,b_r,w_1,...,w_r make the honest bit in

  F(A_0+b_0,w+L_0(A_0),advice_0),
  (F(A_i+b_i,w_i+L_i(A_i),advice_i))_(i=1)^r

more than theta from independent uniform jointly with the other outputs.
The bad-seed density is at most2epsilon/theta.

Choose a witnessing right tuple for each bad w, and take W uniform
independent of the left tuple. Set B_i=b_i(W),Y_0=W+L_0(A_0),
Y_i=w_i(W)+L_i(A_i). Condition on all L_i(A_i). The left and right tuples
then are independent and Y_0 is uniform. The supplied bounded-range entropy
inequality loses at most ell bits of average conditional entropy; Markov
leaves Hinf(A_0)>=K outside probability epsilon. On each good branch (F)
has error at most epsilon. Retain the leakage values and seeds; they
recover W=Y_0+L_0(A_0). Averaging and dropping the leakage gives joint TV
at most2epsilon. On bad W the chosen witness has conditional TV>theta,
so their density is at most2epsilon/theta. All witnesses are proof devices,
not computed by the evaluator. Fewer than t tamperings can be padded and
their extra outputs discarded.

## 3. Three parities and a robust cubic certificate

Let A=2^(pL)>=2^128. Set theta=1/(16A^2),
epsilon=1/(128C*A^(11/2)). The row advice alpha||(j-1) has exactly
a=pL+lg(C) bits; there are NO input-dependent tags. Define

  R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),alpha||(j-1)),
  D(x)=1 iff 2 sum_alpha R_alpha(x)>A.

Take independent sources X_1,X_2 of min-entropy at least k. For every
nonempty Q of at most three rows choose its least row as honest, and
enumerate its C possible honest neighbors. Expand all row XORs into at
most3C calls; their alpha,j labels are distinct. Apply Section2 with
every A_i=X_1 and seed maps Samp(X_1,alpha,j), of joint leakage at most3Cd.
It suffices that

  k>=max(k_S,K+3Cd+lg(1/epsilon)).                       (D1)

There are at most3CA^3 instances, without advice-profile conditioning.
Their bad-seed union B has density at most

  3CA^3*2epsilon/theta=(3/4)A^(-1/2)<=A^(-1/2).                  (D2)

B depends on the law of X_1, not on the second-source point, and is never
constructed. Apply the sampler to X_2. With probability at least99/100,
a fixed x_2 has a row set G omitting q<=sqrt(A)-1 rows, g=|G|>=A-sqrt(A)+1, with
a neighbor per row outside B. For any Q subset G use its honest row's good
neighbor, shifts B_i=x_2 and seed offsets Samp(x_2,alpha,j).
Source-linearity gives the exact seed identities. The selected parity is
theta-close to fair, since XOR of all other calls with the honest bit
preserves the independent uniform bit. Thus its absolute sign expectation
is at most eta=2theta=1/(8A^2). No independence of the full vector is
asserted. The seed-map entropy loss is already included in Section2.

For the ACTUAL good-row sign sum S, cancelling repeated indices gives

  |E S|<=g eta, g-g^2 eta<=E S^2<=g+g^2 eta,
  |E S^3|<=g^3 eta.                                    (D3)

Odd tuples leave a nonempty set of one or three indices. Put
r=sqrt(A),delta=r-1, proof-only integers since pL is even. Let Y be
the actual total sign sum, so |Y-S|<=q<=delta. If D is constantly zero,
Y>=0 implies S>=-delta; then Q_+(S)=(S+delta)(S-r)^2>=0.
If D is constantly one, Y<0 implies S<delta; then
Q_-(S)=(delta-S)(S+r)^2>=0. Yet (D3) makes both expectations negative:

 E Q_+(S),E Q_-(S)
 <=-(r+1)g+r^3-r^2
     +eta*(g^3+(r+1)g^2+(r^2-2r)g)
 <=-A-1 + A/8+(r+2)/8 <=-3A/4-1<0.                       (D4)

The first line expands Q exactly; the absolute linear coefficient is
r^2-2r for r>=2. The base uses g>=r^2-r+1, giving -(r+1)g+r^3-r^2<=-A-1.
The error bound uses g<=A and r+2<=r^2=A for r>=2. This contradiction
includes ties and needs no perturbation or distributional approximation.
Thus D is nonconstant on the support of X_1+X_2.

## 4. Fix the constants and prove finite feedback

For each h>=11 put u=h+1 and

 B0=24u^2+895u, T=128*2^h*u, p=b*B0,
 C=D0*(1+ceil(kappa_G*B0)),
 Fbound=C_F*(3C+1)^9*(13p/2+2).

Select the least h with Fbound<=T. This is a fixed exact integer search.
Here b is a multiple of1024, so pL is even. For u>=12, B0<=100u^2.
Put cP=100b,cC=D0*(2+ceil(100kappa_G)),J=C_F*(3cC+1)^9*(13cP/2+2).
Then Fbound<=J*u^20. At h_*=max(512,2lg(J+1)), the elementary bound
19log2(h_*+1)<=h_*/2 gives Fbound<T/128, so the search halts.
At512 the left side is below190<256; its difference from h/2 decreases
thereafter, since ln2>1/2. Every constant is input-independent.

For actual n>=1 take L=lg(n+1), A=2^(pL),s=B0L,d=TL. Request
t=3C,m=1,a=pL+lg(C),E=(11p/2)L+lg(128C).
Then 2^-E<=epsilon and

  m+a+E+L+1=(13p/2+1)L+lg(C)+lg(128C)+2<=(13p/2+2)L,

so d_F,k_F<=Fbound(h)L<=TL=d. Use the first d_F seed bits in the padded
breaker, retaining the full d-bit seeds. For completeness, fix the left
tuple and, at each uniform honest raw seed, maximize output TV over all
right profiles. This is a legal (F) instance since that whole right side
is independent of the left tuple; retention of the raw seed bounds the
average of those maxima by2^-E. Any actual right variable, including padded
suffixes, fixes one such profile when conditioned with the raw seed. Thus
retaining that entire right variable incurs no extra error. No independence
among suffixes is assumed, and maxima are never evaluator searches.

Set

  H=(3C+1)T+11p/2+1.                                      (4)

For K=k_F<=TL, K+3Cd+lg(1/epsilon)<=HL and k_S<=3TL+7<=HL.
Thus k=HL supplies(D1). Require all exact primitive legality checks and

 n>=1, L>=(8h+2)(h+10), HL<=n,
 L>=lg(C)+lg(128C)+2, pL>=t0, v(pL)>=64,
 L>=200R_G+2, 100R_G*lg(v(pL))<=L,
 L>=2ceil(kappa_G*B0)+3.                                       (5)

These guards depend on no unknown source or witness. The walk degree is then
exactly bounded by C. All guards eventually hold at EVERY integer length:
v(pL)=(33R_G/64)*B0*L<=2pL, and log L/L tends to zero. An explicit threshold
follows from log2 L<=sqrt(L) for L>=16: if D=100R_G and
L>=max(16,4D*lg(2p),4D^2), then D*lg(v(pL))<=L. Fixed constant and
linear guards also hold eventually. Finally H*lg(n+1)<=n eventually,
since2^j/(j+2) tends to infinity for j=floor(log2 n).

## 5. The fixed total evaluator and native bound

Precompute the fixed constants in Section4. On input x of actual length n,
output0 if n=0 or any guard fails. Otherwise build the canonical extractor
in (L), seed-padded to s, and the (G) graph on pL-bit row indices, with its
s-bit projection. Pad its degree by repeating an existing neighbor to C.
For each alpha in lexicographic {0,1}^{pL} and j=1,...,C compute

  F(x,E_n(x,Gamma(alpha,j)),alpha||(j-1)).

XOR within each row and output its specified strict majority over rows.
There is no seed, advice, tag or source input besides x. All field and graph
setup is the supplied canonical finite setup, and all (F) parameters are
computed by its supplied algorithm. Polynomially many rows A<=2^p(n+1)^p
and calls of polynomial actual-n cost give a uniform polynomial evaluator;
all loop exponents and setup constants are fixed. The fixed height search
is finite constant setup independent of n. Its cost may be enormous.

For an affine flat S=a+V of dimension at least k=HL, take X_1 uniform on V
and X_2 independently uniform on a+V. Both have entropy at least k, and
their sum is uniform on S. Section3 proves nonconstancy on S. The ordinary
native bridge proved at the start applies with integer d=k-1, yielding

  C_fullB2(f_n)>=4n-8k=4n-8H*lg(n+1)
                  >=4n-16H*log2(n)

eventually at every integer actual n>=2. No directional or affine-tampering
property is needed or claimed for this family. The native model permits all
sixteen binary tables, free constants, unrestricted sharing and fan-out.
Set f_0=0. The NP verifier accepts exactly the empty witness when this fixed
total evaluator outputs1 and rejects every nonempty witness. Tags or advice
are not circuit inputs. The exact eventual guarantee uses K=16H and one
common threshold. Its projection is (4,1).
