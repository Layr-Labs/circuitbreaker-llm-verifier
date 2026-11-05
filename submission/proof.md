# Fixed logarithmic-loss family and native full-B2 theorem

This file and its five companion proofs support the exact declaration in
bound.json. Mathematical arguments are supplied below.

## 0. Component map

Material source: kaankolcu, submission0165533a-47ba-4e5b-a0c8-3cc7605b3da0,
immutable56cb98a8efa90be26a2804e41f37e9ec143e48f5. Its complete family,
NB, graph and native bridge are retained here. Used-prefix entropy and
bad-set lifting build on submission865cdc91-9181-45f7-a233-ad62151d36f7.
No independent rediscovery or prior-judgment premise is claimed. The new
analysis and explicit coefficients are stated in outer SectionsD1--D3,4,7.

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
2. Lemma A (Section 3) at degree 21 around it: NB's inner breaker F.
3. The CZ/IP advice condenser, the AGHP extractor and F: the standard
   breaker NB (affine-and-standard-wrappers.md Section 4).
4. Lemma A at the sampler's degree with C=NB (Section 5): the outer
   generalized affine breaker.

Every constant below is specified by integer/rational operations or a search
whose termination is proved. No search tests an extractor property or a source
distribution. Logarithms are binary. Let lg(x) be the least nonnegative integer
j with 2^j >= x, for positive rational x. All real-looking powers appearing in
error bounds are analytical bounds, not instructions to compute floating-point
numbers. All chosen parameters can be computed with exact integer arithmetic.

## 1. Inputs and their precise scope

The finite theorem in `../circuit-log-primitives-20261105T0130/seeded-linear.md`
is used with its exact K recurrence and exact seed length, not an imported
GUV/CGL asymptotic theorem. Its output may be truncated, and its seed may be
padded by ignoring a suffix. `constant-advice-breaker.md` Section 8 supplies
a one-tampering standard breaker for arbitrary correlated source and seed
pairs, retaining both complete seeds.

The deterministic bilinear primitive and the AGHP extractor with its weak-seed
lift are the sibling modules `bilinear-extraction.md` and `weak-seed-hashing.md`
in that same primitives directory. The affine composition is Sections 1--4 of
`../circuit-log-extractor-research-20261105T0106/notes/affine-transformation.md`.
The standard breaker is Theorem NB (affine-and-standard-wrappers.md Section 4);
it replaces the earlier equal-length and general-input wrappers.

For the advice condenser we use `cz-conditional-entropy.md` Section 8 with
eta=1/10000000, and its expansion input from `finite-approximate-ring.md`
Section 5: gamma=1/10000 and prime threshold 2^20000. The CZ source and prime
size guards below imply that threshold. Its canonical prime setup is the one
in `prime-fields-and-encoding.md`. From `ordinary-condenser.md` only
its scalar growth theorem G (Sections 5--7) is used, by the sampler. These
are declared mathematical inputs, not proved again or strengthened here.
The ring theorem has exactly the expansion contract left conditional in the
CZ file. Substitution with the displayed gamma and eta closes that premise.

This note does NOT assume normalized contract (F) in
`recursive-sampler-schedule.md`, nor any closed-form degree envelope: Lemma A
gives the outer breaker's seed and threshold as Qw and K_c*w, Lemma B bounds
NB's two coefficients at the outer call, and outer Lemma F1 bounds w.

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
guard must already be built into the interface. The constant-advice breaker
is used only at its one fixed advice length 4, inside NB's inner breaker.
The standard breaker NB (affine-and-standard-wrappers.md Section 4) works at
every a>=1.

### 3.2 Lemma A: a two-parameter affine assembly

Let n,t,a,m,E>=1 and integers delta0,kappa>=1 be given. Put

    e=E+lg(26t),  j=lg t,  R=t+5,
    h_A = least h>=11 with 2^h>=64(h+25) and 2^h>=12 max(delta0,R),
    B=2^h_A,
    T_c=(t+1)(2+2delta0+6j),
    K_c=max(2delta0+1, T_c+(t+2)R+2),
    Q=max((t+1)(1+delta0)+kappa+1, T_c+t+4).                        (A1)

Here h_A is a seeded-extractor height, while j is the number of
affine-merging rounds. They are different parameters. Let w be a power of
two passing the three tests

    w>=m+a+e+1,   w>=lg(Qw+1),   w>=G_U, U=max(n,Qw),                (A2)

with G_U from Section 2 at height h_A and exponent e. Let C be a standard
breaker with the contract of Section 3.1 (arbitrary correlated source pair
independent of arbitrary correlated seed pair, honest seed uniform, distinct
advice, both complete seeds retained, error <=2^-e), and suppose that C,
called at source length d=Qw, output w, advice a and exponent e, has actual
seed length <=delta0*w and actual threshold k_C<=kappa*w. Set

    d=Qw, k=K_c*w, s0=dx=dy=w, d0=delta0*w, r=Rw,                   (A3)

pad C's seed to d0 bits by ignoring a suffix, and use (2.1), at the common
h_A and e, to instantiate

    L0: input n, output d0, padded seed s0;
    Lr: input n, output r,  padded seed dx;
    E0: input d, output dx, padded seed dy;
    Lm: input r, output dy, padded seed dx.

Return the first m bits of the r-bit final output; the extension branch for
m>r is never needed because w>=m. If k<=n, the evaluator of
affine-and-standard-wrappers.md Section 2.1 with these maps is a degree-t
generalized affine breaker with output m, seed d, threshold k and error
<=2^-E, retaining every seed (contract (4) there).

Proof. The h_A rule gives B>=64(h_A+25), hence 128P_u<=G_u (Section 2), and
B>=12max(delta0,R).
Ambient lengths. L0 and Lr read n; E0 reads d=Qw; Lm reads r=Rw<=Qw because
Q>=T_c+t+4>t+5=R. All are <=U, and G_u, P_u are nondecreasing in u, so every
call has G_u<=G_U<=w and P_u<=w/128.
Seeds and thresholds. The targets are z in {delta0*w, Rw, w}, so
6z/B<=6max(delta0,R)w/B<=w/2 and (2.1) gives

    D(u,z)<=w/2+w/128<w,
    K0<=2delta0*w+w,  Kr<=2Rw+w,  KE<=2w+w,  Km<=2w+w.               (A4)

Each primitive's seed fits its padded seed (s0, dx, dy, dx are w). Legality:
K0,Kr<=k<=n since K_c>=2delta0+1 and K_c>=T_c+(t+2)R+2>2R+1; KE<=3w<=d and
Km<=3w<=r since Q>=3 and R>=6. C's guard: k_C<=kappa*w<=d since Q>=kappa+1.
C's seed is padded from at most delta0*w to exactly d0.
Contract (3) of affine-and-standard-wrappers.md Section 2.1, with eps=2^-e and
q=e. The transcript budget is

    T=(t+1)(s0+2d0+dx)+j(t+1)(4dy+2dx)=(t+1)(2+2delta0)w+6j(t+1)w=T_c*w.

The first test of (A2) gives m<=w and e<=w. Then s0=w<=d, dy=w<=r and m<=r;
r=(t+5)w>=(t+1)w+3w+e>=(t+1)dy+Km+q; k=K_c*w>=max(K0,T+Kr+t*r+e) because
K0<=(2delta0+1)w and T+Kr+tr+e<=(T_c+2R+1+tR+1)w=(T_c+(t+2)R+2)w;
d=Qw>=(t+1)(s0+d0)+k_C+e because the right side is at most
((t+1)(1+delta0)+kappa+1)w; and d>=T+KE+t*dx+e because the right side is at
most (T_c+3+t+1)w. Every map has error <=2^-e: the primitives by (2.1) at
exponent e, C by hypothesis. Contract (4) there then gives error
(13*2^j-8)2^-e<26t*2^-e<=2^-E, since 2^j<2t. Seed retention is part of (4).
This checks actual dimensions rather than using uncalibrated O notation.

Remarks. (i) Lemma A only needs w to pass (A2), so a caller may add tests to
the selector; NB's inner breaker (affine-and-standard-wrappers.md Section
4.2) does so. (ii) The least power of two passing (A2) exists and is bounded
in Section 3.3; at the outer call, outer Lemma F1 gives an explicit bound.

### 3.3 The selector halts, with a polynomial bound

With h=h_A, B=2^h_A and C=64B(h_A+1) as in (2.2), put

    V=e+h+m+a+lg(n+1)+lg(Q+1)+10.

Let w* be the least power of two >=16 C^2 V. Then
w*<32 C^2 V. As n,Q,w*>=1,

 lg(max(n,Qw*)+1)<=lg(n+1)+lg(Q+1)+lg(w*+1),
 lg(w*+1)<=6+2lg C+lg V<=6+2C+V.

It follows that the bracket in (2.2) is at most 2V+2C+6<=4CV,
so G_max(n,Qw*)<=4C^2V<=w*/4. The other two tests in (A2) also hold:
the first is bounded by V, and the second by 2V+2C+6<=4CV<=w*.
Thus w<=w*<32C^2V. The selector consists only of finite dyadic trials.

For parameter choices with k>n, define the affine map as the all-zero m-bit
string and make no extraction claim. This is a total computable family.

## 4. Fix the advice-condenser constants first

Set

    eta=1/10000000,     theta=eta/100000,     beta=theta/2.

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

## 5. The standard breaker NB and Lemma B

Theorem NB (affine-and-standard-wrappers.md Section 4) is a standard breaker
with the contract of Section 3.1 at c=c_std=480594400667724074115063104, with
the single guard k_C<=n, at every a,z,E>=1. In the dependency order above it
uses the constant-advice breaker inside Lemma A at degree 21 (its inner F),
the advice condenser with the constants of Section 4, and the AGHP extractor.
The outer assembly uses NB only through Lemma A, and Lemma B gives the two
coefficients that Lemma A needs at the outer call.

**Lemma B (NB at the outer call).** In Lemma A let C=NB, with the outer
parameters of outer Section 2 (t=4C_S-1, m=1, a and E from (O2), where C_S is
the sampler degree), and let w pass (A2). If Qw<=n, then NB's actual seed
length d_C<2^18*w and its threshold k_C<2^59*w.

Proof. NB is called with output z=w, exponent e and source length n'=Qw; put
Ee=e+4 and L'=lg(n'+1). By (A2), w>=G_U=2B(4b+S_U+5)>=8Bb with
b=e+h_A+4>Ee, so Ee<w/(8B). The outer advice a=pL+2+lg C_S is at most
(6p+1)L+33+lg C_S<=e, so a<w/(8B). Since Qw<=n, L'<=L=lg(n+1), and
L<=e/(6p+1)<w/(8B(6p+1)). Also w>8Be>=8B(6p+1)>48Bp, and B>=12R by the h_A
rule.
(i) By (kC) of affine-and-standard-wrappers.md Section 4.7, at z=w and L=L',

    k_C < kz*w+kE*Ee+kL*L'+ka*a+k0
        < [kz+(kE+ka)/(8B)+kL/(8B(6p+1))+k0/(48Bp)]*w.

With 8B>=96R>4.758*10^14 and 48Bp>=576Rp>5.312*10^27 the bracket is below
3.29871*10^17+(1.6534*10^25)/(4.758*10^14)+(5.297*10^23)/(4.758*10^14)
+(2.2543*10^27)/(5.312*10^27)<3.29874*10^17<2^59=576460752303423488.
(ii) By (dC) of affine-and-standard-wrappers.md Section 4.7, at z=w and
L=L',

    d_C < dz*w+dE*Ee+dL*L'+130a+d0
        < [dz+(dE+130)/(8B)+dL/(8B(6p+1))+d0/(48Bp)]*w.

With the same lower bounds for 8B and 48Bp the bracket is below
201015.1+(1.0075*10^13)/(4.758*10^14)+(3.228*10^11)/(4.758*10^14)
+(1.3737*10^15)/(5.312*10^27)<201016<2^18=262144.

NB's own guard holds: k_C<2^59*w<=Qw=n' because Q>=kappa+1. So Lemma A
applies at the outer call with (delta0,kappa)=(2^18,2^59).

## 6. Uniformity, exactness, and remaining obligations

All parameter searches are now proved terminating: the dyadic selectors of
Lemma A (Section 3.3), of the constant-advice breaker and of NB's inner
breaker (Lemma F there), and finite constant formation. Exact seed lengths and entropy thresholds
are computed using the original primitive recurrence, even though the
conservative envelopes decide a legal pad length. The evaluator never
receives entropy witnesses, good-row labels, source decompositions, coupling
coins, proof histories, or a per-length truth table.

For fixed t and outer-logarithmic m,a,E, auxiliary dimensions are O(log n).
The seeded-linear algorithms are polynomial in their dimensions and B;
the finite CZ/AGHP/bilinear setup costs 2^O(local dimension),
so their complete composition is polynomial in the actual outer n. Fixed
constants can be enormous. This does not assert polynomial bit-time jointly
in arbitrary t and n: the exponent of n may depend on t. The outer
assembly fixes t=4C-1 once, with its constants.

For ordinary seeded primitives or standard breakers with the complete
source tuple independent of the complete seed tuple, unused suffixes may be
retained by their conditional kernel given all used prefixes. The outputs
depend only on the source tuple and those prefixes. This applies to the
primitive pads in Section 3 and to the zero-additive-shift independent-side
use of NB's inner breaker (affine-and-standard-wrappers.md Section 4.6, step 7). It does not assert a generic
kernel principle when arbitrary additive B-tuples correlate with suffixes.

# Tagged affine-permutation logarithmic-loss family: outer assembly

All constants below are fixed once; none is length-dependent advice.
The affine-breaker contract (F') is proved below from Lemmas A, B and F1.

This retained peer composition supplies the complete finite seeded-linear,
scalar-growth, graph and affine-breaker constructions
credited in their companion modules. The sumset framework is due to
Chattopadhyay--Liao, https://arxiv.org/abs/2110.12652v1 , with the
advice-condenser breaker architecture following Xin Li,
https://eccc.weizmann.ac.il/report/2023/023/revision/1/ . The native full-B2
bridge is credited in minimal-circuit-bridge.md. No earlier accepted result
or model judgment is a mathematical premise.

All logarithms are binary. Write lg(x)=ceil(log2 x) for positive x, computed
by integer bit length and a power-of-two comparison. All string lengths,
indices and loop counts are integers. All arithmetic parameter algorithms
use exact integers and rational comparisons, not floating point.

## 1. The exact inputs to this assembly

(L) The finite source-linear strong extractor in seeded-linear.md has the
following specialization. Fix integer h>=11 and put u=h+1. For actual input
length n, let L=lg(n+1), requiring L>=(8h+2)(h+10). With primitive error
exponent L and leaf output 14uL, its output length and padded seed are

  d=14*2^h*u*L, s=B0(h)L, B0(h)=24u^2+211u.

Indeed a_E=L+h+4 and rho=a_E+2lg(4n)+2<=3L+h+10. The exact
power-of-two selector satisfies S_E<=2(4h+1)rho, so

  S_E/L<=24h+6+(8h+2)(h+10)/L<=24h+7=24u-17.

The general seed bound D<=6ell+(h+7)S_E+18a_E+12 therefore gives

  D/L<=84u+(u+6)(24u-17)+18+(18h+84)/L
     <24u^2+211u-66<B0(h).

Here L>=h+6 makes the last fractional term <18. Also a_E<=2L,
S_E<=24uL, and 4<=4L. Thus
K0<=2*2^h(14uL+4a_E+S_E+4)
  <=2*2^h(38uL+12L)<=78*2^h*uL<6d, since u>=12.
The actual guard K0<=n is required. Pad the seed to s by ignoring its suffix.
The sampler at failure1/100 uses entropy at most 6d+7. This extractor
is source-linear and has the supplied uniform polynomial evaluator.

(G') Theorem G' of sampler-and-disperser.md Section 4, with P=2437500,
p=P*B0(h), C=2p/3, A=2^(pL) and Lambda_S=A^(1/2)/2^15 (Lambda there): a
deterministic graph Gamma from left pL-bit strings to right s-bit strings, of
degree exactly C, such that whenever L passes the walk guard (W) there, every
set of at least Lambda_S left vertices has more than 2^(1-L)*2^s distinct
neighbors. Gamma is
the four rows of one AB+C step followed by a short-start Cayley walk and an
s-bit prefix; its uniform evaluator costs 2^O(pL). With the base sampler of
(L) it gives (23') there: for every k_S-source and every test of density at
most epsilon0=2^-L, fewer than Lambda_S rows are all-bad, except on source
probability 1/100.

(F') For every n passing the guards (5') of outer Section 4, Lemma A with
C=NB, at the request (O2) of outer Section 2, returns an exact seed length
d_F, source entropy requirement k_F and a one-bit generalized affine
breaker, where (outer Lemma F1)

  d_F=Qw_A<=TL,   k_F=K_c*w_A<=TL.                                      (F')

For arbitrary correlated tuples (A0,...,At) independent of
(B0,...,Bt,Y0,...,Yt), with H_inf(A0)>=k_F, Y0 uniform on d_F bits, and
alpha0 different from every alphaj, its honest output on Xj=Aj+Bj is
2^-E-close to uniform jointly with every tampered output and EVERY seed
Y0,...,Yt. The field is F_2, so addition is XOR. Neither source copies nor
seed copies are required to be deterministic functions of their honest copy.
For the fixed t and a,E=O(log n), its complete evaluator is polynomial in
the actual n (Lemma A's evaluator with NB, whose cost is stated in
affine-and-standard-wrappers.md Section 4.7).

## D1. Algorithm and finite hypotheses

Write L=lg(n+1), s_tag=log2 A=pL and w=2, with p>=128 and L>=1024.
Build the canonical field F_(2^L), and split x into L-bit coefficients of
P_x, padding its final block with zeros. For a full pL-bit row alpha,
let beta,gamma be its first and second L-bit blocks. Define
H_alpha(x)=Prefix_2(gamma*P_x(beta)), a two-bit linear map in x.
Each (beta,gamma) has exactly A/2^(2L) row preimages. For Delta!=0,
P_Delta is nonzero of degree ceil(n/L)-1<n/L, so its root fraction is
at most n/(L*2^L)<=1/L. Off those roots, multiplication by P_Delta(beta)
permutes gamma and exactly 1/4 of gamma values have zero two-bit
prefix. Thus at most A*(1/L+1/4)<=257A/1024 rows have H_alpha(Delta)=0.
The root bound follows by factoring X-beta at distinct roots. The supplied
finite-field proof makes the canonical search halt in 2^O(L) time.
Full alpha remains in the advice label despite repeated (beta,gamma).
No independence of hash values across different rows is asserted.

Use the existing source-linear sampler Samp(x,alpha,j), j in [C], seed
length d, tolerating test density epsilon0=2^-L with fewer than Lambda_S
all-bad rows ((23') in (G')), failure probability 1/100, at source
entropy k_S. The arbitrary enumerated field index alpha replaces [A].
For t=4C-1 and a=s_tag+w+lg(C), request the padded all-seeds-retained affine
breaker F of output length one, entropy K, and error at most

    epsilon = 2^(-L)/(2^33 C A^6).

The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1, and w-bit H_alpha(x). Define

    R_alpha(x)=XOR_(j=1)^C F_use(x,Prefix_dF(Samp(x,alpha,j)),advice(x,alpha,j)),
    D(x)=1 iff 2 sum_alpha R_alpha(x)>A.

Assume the input source entropy satisfies

    k >= max(k_S, K+4w+4C*d_F+2lg(1/epsilon)).                 (D1)

Here F_use is the raw d_F-seed map in(F'); the padded F evaluates precisely
this prefix expression. Set K=k_F=K_c*w_A, the actual raw threshold in(F').
The evaluator is unchanged, including all guards.
All concatenations, zero padding and majority ties are specified. There
are no seed or advice inputs to D. The tags are internal computations on
the actual n-bit x.

## D2. Bad seeds with arbitrary correlated left copies

All seeds in this lemma have length d_F and F_use denotes the raw(F')
map. Each seed map is the d_F-bit prefix of its corresponding sampler map.
Fix an arbitrary correlated tuple (A_0,...,A_r), fixed advice differing
from the honest advice, and deterministic seed maps L_i of A_i. Assume
Hinf(A_0)>=K+ell+lg(1/epsilon), where the joint seed leakage has at most
2^ell values. Call w bad if some b_0,...,b_r,w_1,...,w_r make the honest
bit in the tuple

  F_use(A_0+b_0,w+L_0(A_0),advice_0),
  (F_use(A_i+b_i,w_i+L_i(A_i),advice_i))_(i=1)^r

more than theta away from a fresh independent uniform bit jointly with
the other outputs. Its bad-seed density is at most2epsilon/theta.

Choose one witnessing tuple for each bad w. With W uniform independent
of the full left tuple, set B_i=b_i(W), Y_0=W+L_0(A_0), and
Y_i=w_i(W)+L_i(A_i). Conditional on all L_i(A_i), the complete left and
right tuples are independent, and Y_0 is uniform. The bounded-range
inequality in the retained sampler proof loses at most ell bits of average
conditional entropy of A_0, even when the leakage involves the other left
copies. Apply its conditional lift with breaker error epsilon and
exceptional probability epsilon. The retained honest seed and leakage
determine W. Averaging the local distances over W, whose independence
leaves the original left-tuple distribution unchanged, gives the claimed
density. The arbitrary-copy quantifiers and full seed retention are those
of the supplied actual contract (F'); no copies need entropy except A_0.
All witnesses are finite proof objects, never evaluator searches.

## D3. Actual parities and pointwise collision counts

Fix independent sources X_1,X_2 of entropy at least k, and write
A_0=X_1, A_1=Lambda X_1, B_0=X_2, B_1=Lambda X_2+b, with Lambda
invertible. Suppose tau(x)=Lambda x+b has no fixed point on the sum
support. Row alpha collides when H_alpha(x_1+x_2)=H_alpha(tau(x_1+x_2)).
The tag root bound gives at most N_col=257A/1024 collisions at EVERY point
(x_1,x_2). We use actual bits throughout; no random replacement is made.

For each nonempty set Q of at most four distinct row-copy indices, expand
its at most4C calls, hence at most t=4C-1 tampered calls. If some alpha occurs in only one copy in Q, choose
such an index as honest. Otherwise every alpha in Q occurs in BOTH
copies; call such Q fully doubled and choose its least index as honest.
Enumerate the C honest neighbor choices. Condition on just the SELECTED
left tags (H_alpha(A_sigma))_((sigma,alpha) in Q), at most 2^(4w) cells.
Discard cells of mass <epsilon/2^(4w), of total probability at most epsilon.
On every remaining cell X_1, and also Lambda X_1, has honest entropy at
least k-4w-lg(1/epsilon). Enumerate just the selected right tags,
at most 2^(4w) choices. All selected advice strings are then fixed.

If Q has a singleton alpha, its honest advice label alpha,j differs from
all other selected calls' labels, regardless of collisions at other rows.
If Q is fully doubled, the selected profiles determine all its collisions.
On profiles with no selected collision the honest advice is again distinct
from all tampered advice; two copies of the same label have distinct tags.
Tampered advice strings may coincide with one another in either case.
Section D2 therefore applies in all these cases, with joint USED-prefix seed leakage
at most4C*d_F. There are at most sum_(i=1)^4 binom(2A,i)<=64A^4 sets Q,
and hence at most64C A^4*2^(8w) instances including both profiles and honest
choices. It suffices that

  k>=max(k_S,K+4w+4C*d_F+2lg(1/epsilon)).              (D1')

Choose theta=1/(1024A^2), epsilon=2^(-L)/(2^33*C A^6), and let B_used be the
union of their bad-seed sets in the d_F-bit seed space. As 64*2^16*2*1024=2^33, Section D2 gives

  |B_used|/2^d_F<=64C A^4*2^(8w)*2epsilon/theta=2^33 C A^6 epsilon=2^(-L).    (D2')

Lift B_used to B=B_used x {0,1}^(d-d_F). Its full-d density is EXACTLY
the used-prefix density, and a full sampler seed avoids B precisely when
its used prefix avoids B_used. This lift is a proof object. No ignored
source-dependent suffix is conditioned on or assumed independent.
B depends on the left source law and tau, not the right source point;
the program never constructs B. Apply the sampler to both B_0 and B_1,
each of entropy k. Except on probability2/100 both conclusions hold.
Fix such an x_2. Intersect the good row sets to obtain G, omitting
q<=2(Lambda_S-1)<2^(-14)A^(1/2) rows. Since 2^(-14)A^(-1/2)<=1/256,
g=|G|>=255A/256>=A/2. Each sigma,alpha in {0,1}xG has its own
neighbor j_(sigma,alpha) whose honest seed lies outside B. The two copies
need not choose the same j. Source-linearity gives every seed identity
with shifts B_sigma=x_2 or Lambda x_2+b and offsets
Prefix_dF(Samp(B_sigma,alpha,j)). Prefixes preserve source-linearity.
The second sampler is essential; its seeds need not be fixed translates
of the first sampler's seeds.

Put eta=2(theta+epsilon)<=3/(1024A^2). For each Q on G, using its honest
index's good j, the actual parity has absolute expectation at most eta
if Q is not fully doubled. If fully doubled, the bound is instead

  eta+Pr[at least one selected alpha collides].       (D3a)

Indeed good cells/profiles without collisions have bias at most2theta;
exceptional entropy cells have total mass at most epsilon and collision
cells cost at most their probability. These costs sum to at most
2theta+epsilon+Pr[collision]<=eta+Pr[collision].
Bad-cell bounds apply individually to each Q; simultaneous good cells
are not assumed. All probabilities here are over X_1 for the fixed x_2.
The pointwise collision bound implies sum_(alpha in G) Pr[collision at alpha]
<=N_col. No independence of rows, copies, or collision events is asserted.

## D4. A robust quartic sign certificate

Let X_alpha,Y_alpha be the actual signs in the two copies on G, with
S_0=sum X_alpha and S_1=sum Y_alpha. Put Z_alpha=X_alpha Y_alpha,
Z=sum Z_alpha and zeta=E Z. Section D3 gives

  |zeta|<=g*eta+N_col.

In S_0^3 S_1, the triples reducing to one row have coefficient 3g-2.
Every other surviving parity has a singleton row label, hence bias <=eta.
The same applies to S_0 S_1^3. In S_0^2 S_1^2, all-even tuples contribute
g^2. Fully doubled two-row sets contribute
4 sum_(alpha<beta) E Z_alpha Z_beta=2(E Z^2-g).
Every other surviving parity again has a singleton. Thus, with
|e31|,|e13|,|e22|<=g^4 eta,

  E S_0^3 S_1=(3g-2)zeta+e31,
  E S_0 S_1^3=(3g-2)zeta+e13,
  E S_0^2 S_1^2=g^2-2g+2E Z^2+e22.                 (D3')

Marginal fourth moments have only single-copy surviving parities, so
E S_sigma^4<=3g^2+g^4 eta<=4g^2. No collision error is needed there.
The elementary variance identity gives E Z^2>=(E Z)^2=zeta^2.

Define Q_+(x,y)=xy(x-y)^2 and Q_-(x,y)=-xy(x+y)^2.
Both have central term -2x^2y^2. Set r=|zeta|/g. Since
N_col/g<=(257/1024)*(256/255)=257/1020 and eta<=3/(1024A^2),
we have 0<=r<13/50. Both quartic expectations are at most

  2(3g-2)|zeta|-2g^2+4g-4zeta^2+4g^4 eta
  <=[6r-4r^2-2+4/g+12/1024]g^2<=-g^2/2.          (D4')

For the final inequality, 6r-4r^2 increases on [0,13/50],
and 6r-4r^2-2 equals -444/625 at13/50; also 4/g<1/1000.
This uses the exact matched-pair second moment rather than bounding each
collision term separately. All identities concern actual deterministic signs.

Let Y_sigma be the ACTUAL sum over ALL rows in that copy. Pointwise
|Y_sigma-S_sigma|<=delta=q<2^(-14)A^(1/2). Since g>=A/2,

  delta/sqrt(g)<2^(-14)A^(1/2)/(A/2)^(1/2)=2^(-27/2)<1/10000.

The L4 norms of S_sigma are at most2sqrt(g), and those of Y_sigma at
most3sqrt(g). The norm triangle inequality follows from Holder applied
to (|X|+|Z|)|X+Z|^3 and division by the cube of the L4 norm (zero norm
is trivial). For four factors, Holder itself follows by Cauchy-Schwarz
twice: E|abcd|<=(E a^4 E b^4 E c^4 E d^4)^(1/4).

Telescoping the four factors in any degree-four monomial shows its
expectation changes by at most4delta(3sqrt(g))^3. Each Q has sum of
absolute coefficients four. Therefore

  |E Q_+(Y)-E Q_+(S)|, |E Q_-(Y)-E Q_-(S)|
    <=432delta g^(3/2)<(432/10000)g^2<g^2/8.            (D5')

Both actual expectations are strictly negative. If the actual majority
XOR were identically zero, Y_0Y_1>=0 pointwise and Q_+(Y)>=0, impossible.
If that XOR were identically one, Y_0Y_1<=0 and Q_-(Y)>=0, also impossible.
The weak product inequalities include all ties: the evaluator outputs one
exactly for a negative sign sum, and zero for a nonnegative sum.
Thus both actual majority-XOR values occur. No sixth or eighth moments,
small-tail inference, or independence of the full sign vector is used.

All sums and collision events above are actual deterministic quantities.
For an affine S=a+V of dimension at least k, take X_1 uniform on V and
X_2 independently uniform on a+V. Their sum lies in S. Good x_2 have
probability at least98/100. Fix any such x_2; the preceding contradiction
proves both values of D(x)+D(tau(x)) occur on S. Thus the property holds
for every invertible affine tau without fixed points on S, and includes
every nonzero translation.


## 2. Fix every constant before seeing n

Fix h=173, a constant of the construction (there is no search), and put

  u=h+1=174,   B0=24u^2+211u=726624+36714=763338,
  T=14*2^h*u=609*2^175
   =29165305762103947335632291860454321650744128919056678912,
  P=2437500,  p=P*B0=1860636375000,
  C=2p/3=1240424250000 (the sampler degree; written C_S in Lemma B),
  t=4C-1=4961696999999,

where 2^175=47890485652059026823698344598447161988085597568237568.
Lemma A at C=NB uses delta0=2^18, kappa=2^59 (Lemma B) and

  R=t+5=4961697000004,  j=lg t=43 (2^42<t<=2^43),
  h_A=46, since 12max(2^18,R)=12R=59540364000048 lies in (2^45,2^46]
     and 2^46>=64(46+25), while 2^45<12R;  B=2^46,
  T_c=(t+1)(2+2*2^18+6*43)=4961697000000*524548=2602648237956000000,
  (t+2)R=24618437119833808485000004,
  K_c=max(2^19+1,T_c+(t+2)R+2)=24618439722482046441000006,
  Q=max((t+1)(1+2^18)+2^59+1,T_c+t+4)=2602653199653000003,

so max(Q,K_c)=K_c. These constants depend on no input length, source,
experiment or oracle. The division by three is integral because 6 divides P.

For each actual n>=1 let L=lg(n+1), A=2^(pL), s=B0L and d=TL. Request

  t=4C-1, m=1,
  a=pL+w+lg(C), w=2,
  E=(6p+1)L+lg(2^33*C).                                       (O2)

Here lg C=41 (2^40<C<=2^41) and lg(26t)=47 (26t=129004121999974 lies in
(2^46,2^47]), so Lemma A's exponent is e=E+lg(26t)=(6p+1)L+e0 with
e0=33+41+47=121. Then 2^-E<=2^(-L)/(2^33*C A^6)=epsilon. The advice
concatenates the pL-bit alpha, lg(C)-bit j-1, and w-bit polynomial-evaluation
tag from Section D1. It has exactly a bits and is evaluated internally.
Below, w_A denotes Lemma A's selected power of two; w=2 is the tag length.

## 3. Termination and the exact feedback inequality

Put

  g1=4(6p+1)+2(4h_A+1)(6p+3)=2244p+1114=4175268025501114,
  g0=4(e0+h_A+4)+2(4h_A+1)(e0+h_A+10)+5=684+65490+5=66179.

The constants of Section 2 satisfy

  4B*max(Q,K_c)*(g1+1) <= T,                                      (F1)

because this reads K_c(g1+1)<=609*2^127, and
K_c(g1+1)=102788584211005831563203964778889868006690<1.028*10^41, while
609*2^127=103615980727425762124597567962973420388352>1.036*10^41.

**Lemma F1.** If L>=g0 and HL<=n (H from (O3) in Section 4), the least power of
two w_A passing (A2) at (O2) satisfies w_A<4B(g1+1)L. Hence Qw_A<TL<=n,
d_F=Qw_A<=TL and k_F=K_c*w_A<=TL.

Proof. For n>=1, lg(4n)=lg n+2<=L+2. With b=e+h_A+4 as in assembly
Section 2, rho_n=b+2lg(4n)+2<=e+h_A+10+2L and
S_n=pow2ceil((4h_A+1)rho_n)<2(4h_A+1)rho_n. Substituting e=(6p+1)L+e0,

  G_n=2B(4b+S_n+5)<2B[4(e+h_A+4)+2(4h_A+1)(e+h_A+10+2L)+5]=2B(g1*L+g0).

Let w*=pow2ceil(2B(g1*L+g0)); then G_n<w*<4B(g1*L+g0)<=4B(g1+1)L, using
L>=g0. By (F1), Qw*<TL<=HL<=n, so U=max(n,Qw*)=n and the third test of (A2)
passes at w*. The first passes since m+a+e+1=(7p+1)L+4+lg C+e0<2B*g1*L. The
second: lg(Qw*+1)<=lg(n+1)=L<w*. So w_A<=w*, and
Qw_A, K_c*w_A<4B*max(Q,K_c)(g1+1)L<=TL<=HL<=n.

So Lemma B's hypothesis holds at every n passing (5'), Lemma A applies there
with C=NB, and (F') holds. Pad the breaker to d-bit seeds by taking its first
d_F bits. All other finite guards are still checked explicitly by the
evaluator.

All original d-bit seeds remain retained by the following consequence of
(F')'s universal right-side quantifiers. A simple kernel argument given only
the used seeds would not suffice when additive shifts correlate with suffixes.

Fix the admissible left tuple A. For an honest raw seed y and a deterministic
profile r=(B0,...,Bt,Y1,...,Yt), let Delta(y,r) be the TV between the output
tuple over A and an independent uniform honest output together with the
unchanged tampered-output tuple. For each y choose a maximizing profile r_y
from this finite domain. Uniform y with profile r_y is a legal (F') instance,
since the whole right side is independent of A and every retained seed is a
function of y. Retention of the honest seed gives

  2^(-d_F)*sum_y max_r Delta(y,r) <= 2^-E.

For an actual full right-side variable R independent of A jointly with the
uniform honest raw seed y, conditioning on (y,R) fixes a profile. The joint
output distance retaining (y,R) is the average of its Delta(y,r) values,
and is therefore bounded by the displayed average of maxima. Thus the
WHOLE right side, including all padded seed suffixes, may be retained with
no extra error. No independence among those suffixes is assumed. The finite
maxima are proof devices, never evaluator searches. This proves the padded
d-seed guarantee even when B-tuples correlate arbitrarily with suffixes.

## 4. The original guard and a sharper used-prefix entropy requirement

Keep the evaluator's ORIGINAL guard constant

 H=(4C+1)T+12p+3=144709410103943034489126839574476258663392646546947510738880357178915.                                  (O3)

All actual maps, selected dimensions, h=173, P=2437500 and zero guard
branches remain identical to the published predecessor. Only the proof's
required source entropy and the declared explicit remainder change.
The proof of Lemma F1 gives the sharper literal bound
w_A<4B(g1 L+g0). Since w_A is an integer, put

 D_*=K_c+4C Q=12913601191198413637931532000006,
 W_1=4B g1=1175233470238672725053783670784, W_0=4B g0-1=18627732483734503423,
 H_1=D_*W_1+12p+2=15176496341210369504762710248528659224211697328193411426524706,
 H_0=D_*W_0+74+2lg C=240551108391279267699559511855696315288536443020694.

Then w_A<=W_1 L+W_0. The tagged-parity conditioning costs 4w=8 bits,
and 2lg(1/epsilon)=(12p+2)L+66+2lg C. Thus the actual requirement(D1)
with USED-prefix leakage is at most

 K+4C*d_F+4w+2lg(1/epsilon)
 =D_*w_A+(12p+2)L+74+2lg C <=H_1 L+H_0.

Define the proof-only integer k(n)=H_1 L+H_0. Also k_S<=6TL+7<=k(n):
H_1>10^61,6T<10^57,H_0>=7. Finally H_1+H_0<2*10^61<H, so for L>=1
k(n)<HL<=n whenever the ORIGINAL guard holds. All raw entropy/source
legality checks remain valid; no sampler seed or auxiliary map changes.
Ignored sampler suffixes never entered the left entropy conditioning.
The conditional-copy bad-seed proof and lifted set supplied(D1) directly;
no claim of suffix independence or prior-verdict premise is used.

The evaluator tests these literal finite guards:

  n>=1;  L>=max(1024,(8h+2)(h+10),g0,74+2lg C)=253638;
  6000000000*lg((B0+2)L)<=L;  HL<=n.                              (5')

Here (8h+2)(h+10)=1386*183=253638 dominates 1024, g0=66179 and
74+2*41=156. It also computes and checks the supplied algorithms' exact
parameters and guards. The already proved L-specialization, (G') and (F')
guarantee that these checks pass when (5') holds. No test of an unknown
extractor/graph property, source distribution, or existence witness occurs.

There is one fixed threshold N beyond which (5') holds for EVERY integer n.
Indeed L tends to infinity with n, the walk guard (W) holds at every
L>=2^39 (sampler-and-disperser.md Section 4.1), and the fixed lower bound on
L then holds too. Finally H*lg(n+1)<=n eventually: choose j=floor(log2 n),
so 2^j<=n<2^(j+1) and L<=j+2. The elementary exponential sequence 2^j/(j+2)
is increasing (consecutive ratio 2(j+2)/(j+3)>1) and unbounded, so it tends to
infinity: H(j+2)<=2^j for every j from some fixed j0 on. This also proves
eventual k_F<=n.

This is an all-integer-length argument, not a favorable subsequence or a
liminf used to invent a remainder. A small-input zero branch is allowed
because only an eventual lower bound is claimed.

## 5. The fixed evaluator

Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:

1. If n=0, output 0. Compute L and the parameters above. If any finite
   guard in (5') or its exact supplied parameter check fails, output 0.
2. Build the exact base extractor E_n from (L), seed-padded to s bits.
   Its data argument is x and its output length is d.
3. Build the (G') graph on left pL-bit strings: the four rows of the AB+C
   step, the short-start Cayley walk and the s-bit prefix (sampler-and-
   disperser.md Section 4). Its degree is exactly C.
4. Build the canonical field of degree L and the tags in Section D1.
   For each alpha in lexicographic {0,1}^{pL} and j=1,...,C compute

     y_(alpha,j)=E_n(x,Gamma(alpha,j)),
     b_(alpha,j)=F_n(x,Prefix(y_(alpha,j),d_F),
                       alpha || encode(j-1) || H_alpha(x)),

   using the ONE breaker parameter setting (O2). Let R_alpha be the XOR
   of these C bits.
5. Output 1 iff strictly more than half of the A row bits are one;
   otherwise output zero. This is exactly Section D1's evaluator.

F_n in this pseudocode denotes the raw d_F-seed breaker; the preceding
proof uses its padded d-seed form. The two evaluations are identical.
No arbitrary seed, affine basis, irreducible polynomial advice, entropy
witness, good row, tampering, or source decomposition is supplied to this
algorithm. Constants, field bases, finite prime searches, orderings and
zero branches are fixed by the component algorithms.

The base sampler tests density epsilon0=2^-L. By (23'), for every test of
density at most epsilon0 fewer than Lambda_S rows are wholly bad, except on
source probability 1/100. Source
linearity survives graph index selection. The degree and advice lengths match Section D1, and (D1) was checked
above. Sections D2--D4 prove the affine-permutation property at k=k(n)=H_1 L+H_0,
and in particular the directional property.

## 6. Uniform polynomial time in actual input length

All one-time constant computations terminate, independent of n. Their
finite cost is part of the fixed machine constant, not length-dependent
advice. Parameter arithmetic has O(log n) bits beyond these fixed constant
lengths. There are A*C <= C*[2(n+1)]^p breaker calls, a fixed polynomial.
The base linear extractor is polynomial in n at fixed h. The graph's
prime p_m has m=pL/3 bits and is found in 2^O(pL) bit operations, so its
evaluator over all A rows is also polynomial in n. The walk fields have fixed-algorithm arithmetic.

By (F'), each breaker call is polynomial in n for the fixed t=4C-1 and
a,E=O(log n). Its local prime searches and all finite field setup costs
2^O(log n) with fixed coefficients. The tag field search and arithmetic cost 2^O(L), and all A polynomial
evaluations cost a fixed polynomial in n. Counting the majority uses an O(log A)
bit counter. Therefore the whole fixed total evaluator is polynomial in
the ACTUAL n. The polynomial degree and multiplier can be extremely large
but do not depend on n. No enumeration of extractor truth tables or
verification over all source distributions is hidden in this cost.

## 7. Native circuit bound, NP eligibility and expression

The assembly is nonconstant under every fixed-point-free affine permutation
on a source flat of dimension at least k=k(n)=H_1 L+H_0. In particular, for every
nonzero v and every affine k-flat S its translation derivative is
nonconstant. We prove a finite bridge under the WEAKER hypothesis

  for every nonzero v and affine k-flat S, some x in S has f(x)!=f(x+v).

For integers k>=1 this implies C_fullB2(f)>=4n-4k. For n<=k this is trivial.
Otherwise use the retained affine-map, three-corner-pair bridge with d=k-1,
stopping at a=M-q=k. Separation implies nonconstancy on every flat of
dimension at least k, by choosing a nonzero tangent translation of any
hypothetically constant flat.

At every state with a>=k an unused parameter is impossible: for a free
parameter keep its coordinate direction in a safe a-flat; for a paired
parameter fix its partner safely. Independence gives zero derivative
there, and injective Phi preserves its nonzero direction. Absence includes
the output reference. Pair rule1 is therefore impossible; rule3 creates an
unused free line coordinate at the same a and is impossible as well.
Mixed and paired/paired productive actions create the unused coordinate
identified in the retained ledgers at new a>=k, so only FF actions occur.

Put W=s-2q. Each FF action removes at least two objects and adds one pair,
dropping W by at least 4 as a drops by one. Basic deletions, transfers,
role exchanges and one-child-row elimination do not increase W. Each
remaining pair rule2 or4 removes at least two objects and one pair, also
not increasing W. Termination thus gives exactly n-k FF actions and

  W_T<=W_0-4(n-k)<=|C|-4(n-k).

All parameters lie in the output's undirected incidence component;
otherwise the invertible auxiliary equations split by components and
the output is independent of a parameter, already excluded. In each
component the row/auxiliary block is square: full row and column rank
of invertible A force both inequalities for that block. Each gate-output
variable belongs with its defining object. Thus the output component has
M+s_0 variable vertices and s_0 hyperedges of size at most three. A connected
hypergraph of s_0 such edges has at most 2s_0+1 vertices, since adjoining
an intersecting edge adds at most two. Hence s>=s_0>=M-1. At a=k,

  W_T=s-2q>=k-q-1>=-1.

The only integer case with W_T<0 is q=k, M=2k, s=M-1. The tight incidence
tree proof and paired read-once zero-period lemma, proved in
sampler-and-disperser.md, exclude exactly this case under our weaker
separation hypothesis. Thus W_T>=0, proving the finite bound. Its uniform
additive constant zero is optimal for this weaker hypothesis: n=k=1,
f(x)=x satisfies separation and has zero gates. This does not assert
optimality under the stronger nonconstant-derivative hypothesis.

For this SAME fixed evaluator, at every integer n beyond the original
common threshold, C_fullB2(f_n)>=4n-4k(n)=4n-4H_1*lg(n+1)-4H_0.
Since lg(n+1)<=log2(n)+1 for every integer n>=1, set

 K=4H_1=60705985364841478019050840994114636896846789312773645706098824,
 K0=4(H_1+H_0)=60705985365803682452615958064912874944269574573927791478181600.

Then C_fullB2(f_n)>=4n-K*log2(n)-K0 for EVERY integer actual n>=N. (O4)
Both coefficients are exact62-digit positive integers, each smaller than
the predecessor's69-digit coefficient578837640415772137956507358297905034653570586187790042955521428715660.
This is a strict improvement of the declared explicit remainder, for the
same total family, not merely an unencoded Big-O constant. No optimality
of these conservative coefficients is asserted.
The NP verifier accepts exactly the empty witness when this fixed total
polynomial-time evaluator outputs1 and rejects every nonempty witness.
Tags, fields and seeds are internal, not circuit inputs. All original
all-integer eventual guards and uniformity arguments remain unchanged.
The explicit remainder is exactly K*log2(n)+K0, with projection(4,1).
