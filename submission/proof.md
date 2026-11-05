# Fixed logarithmic-loss family and native full-B2 theorem

This file and its five companion proofs support the exact declaration in
bound.json. Mathematical arguments are supplied below, independently of
any local check, prior judgment, external review or submission status.

## 0. Atomic proof map and selected construction

All mathematical premises are supplied in these six files. A legacy component
filename below means its included titled component, not an external file.
Historical scratch-directory prefixes in component prose have no operative
meaning. Source URLs identify attribution; no cited paper or prior judgment
replaces a proof. Component hashes identify fuller independently authored
source versions and are not external assumptions. Unused alternatives and
historical scope discussions were omitted as recorded in the preparation log.

- linear-and-bilinear-primitives.md contains seeded-linear.md (Sections 1--6,8
  and attribution 10), near-full-seeded.md(1--4), weak-seed-hashing.md(1--4),
  and the complete bilinear-extraction.md. The exact finite weaker linear
  primitive is used; no stronger imported GUV/CGL formula is assumed.
- sum-product-and-ordinary-condensation.md contains finite-approximate-ring.md
  Sections 1--5,7 and ordinary-condenser.md. Its scalar sum-product proof is
  the latter's Sections6--7; the unused alternative ring corollary is omitted.
- restricted-nonmalleability-and-breaking.md contains prime-fields-and-encoding.md
  Sections 1--3, cz-conditional-entropy.md Sections 1--8, and
  constant-advice-breaker.md Sections 1--8. The ring theorem supplies precisely
  the CZ expansion premise at gamma=1/10000; eta=1/10000000 satisfies its
  prime guard. Shared CZ Section 3 is elementary and has no deep assumption,
  so its use by ordinary condensation is noncircular.
- affine-and-standard-wrappers.md is compact-affine-wrappers.md. Its Section 2
  replaces affine-transformation.md for the used m<=r branch. Sections 3,5
  replace condenser-integration.md's advice/equal-length reductions; its Section 3.4
  replaces condenser-integration6.2's restricted CZ/IP reduction; Section 6
  replaces raz-order.md6's general-input proof. Only the weak-seed hash branch
  is selected, not Raz's extractor. The constant-advice breaker is used ONLY
  inside the fixed-degree inner affine breaker.
- sampler-and-disperser.md is compact-sampler-disperser.md. Sections 2--4
  supply the corrected real-entropy deterministic graph;5 supplies
  sampler-conversion-support.md;6 supplies seed-retention-and-leakage.md's
  used all-seeds branch;7--9 supply four-moment-disperser.md. No unused
  point-line incidence theorem, dummy-seed adapter or b=8 graph is assumed.
- This file includes minimal-circuit-bridge.md, breaker-parameter-assembly-v2.md
  and fixed-logarithmic-family-v2.md. References to their unversioned names
  mean these v2 bodies. The finite assembly supplies the exact degree-nine
  contract(F) used by the final family; this discharges the last conditional
  construction input, subject to the complete proofs supplied here. The old
  recursive-sampler-schedule.md is historical and not a premise.

Component status sentences report when each module was written. They neither
introduce a missing external premise nor establish current review approval.
The selected acyclic chain is finite linear extraction -> fixed-advice
breaker -> inner affine breaker; finite ring/CZ + ordinary condensation +
weak hash -> equal-length breaker; near-full extraction -> standard breaker;
finite affine transformation -> outer degree-nine breaker; deterministic
graph/sampler + four moments -> one fixed unseeded family; native B2 bridge
-> the eventual full 4n-O(log n) inequality. All parameter and evaluator
steps for that exact chain are specified in the final two components.

Precise legacy equation aliases: condenser-integration.md Section 5.10's P1
bilinear-prefix calculation is affine-and-standard-wrappers.md Section 5,
equation (17). Its Section 6.2 and interface (N1) are that file's Section 3.4,
equations (11)--(12), with the joint smoothed conditional-entropy law spelled
out there and in the retained CZ Sections 2,8. They are not generic ordinary
extractor assumptions. References to original preparation/audit files in the
sampler source list are provenance only; all their mathematical interfaces
are reproved in the included sampler Sections 1--9. The broader original
circuit-reduction.md is provenance; the complete used bridge follows below.

The direct circuit proof source is GordoAR, submission 6622943f-7b2d-4f8f-b73c-979c5232484a:
https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e/submission
Attribution record:
https://api.yukon.org/api/submissions/6622943f-7b2d-4f8f-b73c-979c5232484a/attribution
Its upstream affine-system framework is Mishra--Santhanam:
https://eccc.weizmann.ac.il/report/2026/260/revision/1/download .
The circuit proof is independently written and proved below, not adopted
through the source's status or an assumed license to republish its prose.


<!-- COMPONENT minimal-circuit-bridge.md SHA256 6f0e38dc61dc1bcb41e48ce4da31b18bcde6a4fe10e424ffc2b8b7181b304ca0 -->

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

## Consequence for a fixed uniform polynomial-time family

Suppose one fixed total deterministic algorithm computes f_n(x) on every
x in {0,1}^n, at every positive integer length n, in time (n+1)^O(1).
Fix an integer H>=1, independent of n. Assume f_n is nonconstant on every
affine flat of dimension at least

    h(n)=H ceil(log_2(n+1)).

At each actual input length n use d=h(n)-1 in (1). Then

    C_fullB2(f_n) >= 4n-8H ceil(log_2(n+1))
                   = 4n-O(log(n+1)).                 (8)

If h(n)>n the disperser condition is vacuous and (8) remains valid, possibly
negative. There is no reindexing, padding length, or subsequence in (8).
A construction only at special lengths does not suffice: its uniform
extension must satisfy the flat property at every other positive length.
H and the algorithm are fixed, with no length-dependent advice.

Set f_0 arbitrarily to make the language total on all finite binary strings.
The stated disperser condition is for n>=1: demanding it literally at n=0
with threshold h(0)=0 would be impossible, since the singleton 0-dimensional
flat is constant. Equivalently use max(1,h(n)) at n=0. The n=0 lower-bound
inequality itself is trivial. The resulting language L={x:f_|x|(x)=1} is
in P and hence in NP with zero witness bits: a verifier runs the algorithm
on (x,empty) and rejects nonempty witnesses; the witness bound is zero.
This proves an implication, not a disperser construction or candidate claim.


<!-- COMPONENT breaker-parameter-assembly-v2.md SHA256 939e500da50d7f28c936b3c7c189ba79e08cd0a7aaef7f7d85783576ada6106e -->

# Exact finite assembly of the two affine breakers

Date: 2026-11-05 UTC. Status: **parameter closure from the stated finite
primitive contracts, not an independent audit of those contracts**. The
ordinary-condenser, CZ conditional entropy and finite approximate-ring inputs
were independently audited before this assembly was frozen; their local PASS
reports do not substitute for an audit of this composition. It does not certify a circuit candidate, evidence-budget fit,
an external review, or a submission. No frozen source module was edited.

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

The finite theorem in `../circuit-log-primitives-20261105T0130/seeded-linear.md`
is used with its exact K recurrence and exact seed length, not an imported
GUV/CGL asymptotic theorem. Its output may be truncated, and its seed may be
padded by ignoring a suffix. `constant-advice-breaker.md` Section 8 supplies
a one-tampering standard breaker for arbitrary correlated source and seed
pairs, retaining both complete seeds. Its independent local audit passed; this assembly still states its use
explicitly.

The deterministic bilinear primitive, weak-seed hash, near-full extractor,
and finite Lemma 5.1 probability/mixture proof are the sibling modules
`bilinear-extraction.md`, `weak-seed-hashing.md`, `near-full-seeded.md`, and
`condenser-integration.md` in that same primitives directory. The affine
composition is Sections 1--4 of
`../circuit-log-extractor-research-20261105T0106/notes/affine-transformation.md`.
The standard general-input wrapper is Section 6 of that directory's
`raz-order.md`, with its near-full primitive replaced by the finite module.

For the advice condenser we use `cz-conditional-entropy.md` Section 8 with
eta=1/10000000, and its expansion input from `finite-approximate-ring.md`
Section 5: gamma=1/10000 and prime threshold 2^20000. The CZ source and prime
size guards below imply that threshold. Its canonical prime setup is the one
in `prime-fields-and-encoding.md`. For the ordinary condenser we use all of
`ordinary-condenser.md`, including its finite scalar-growth constants. These
are declared mathematical inputs, not proved again or strengthened here.
The supplied local audits passed the following exact input versions:

- constant-advice breaker: 60b045160a2ce5c8758794d865162d0c497ac800f1f21da9fd33c4fb7eb7d129
- ordinary condenser: 7a7240e059dd33418502cee33510d9db1e03c5f433cc7493e6b1c2729827eede
- CZ conditional entropy: 2360bd602dddac5017ecb58a25f439322bb46e03807c309ee4f41a9983106355
- finite approximate ring: f92e2b27b39874fd96b6f31aac5bd7c87a167c3e85a90d1bdf1bc1658c528a84

The ring theorem has exactly the expansion contract left conditional in the
CZ file. Substitution with the displayed gamma and eta closes that premise.
A later rejection of an input would make this result conditional on repair;
arithmetic cannot repair a false primitive.

This note does NOT assume normalized contract (F) in
`recursive-sampler-schedule.md`. Its own proved outer bound is coarser in the
degree, but still **polynomial in degree and linear in advice, error exponent,
output length and log of ambient source length**. That is sufficient for an
exponential-versus-polynomial final sampler feedback calculation, once that
calculation is redone with this bound.

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
Q,K_c<=20s^2, and V<=47sH. Therefore
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

The outstanding work is not a missing parameter inequality in this note:
it is independent audit of this assembly against the pinned inputs, integration with the actual finite graph/sampler using
(8.1), and the final whole-family/circuit/evidence-size certificate. If an
audit changes an exported module or guard, the corresponding constants here
must be recomputed; no frozen candidate acquires approval by association.


<!-- COMPONENT fixed-logarithmic-family-v2.md SHA256 56e41083fc186eb66a16c6be7f6c6246770031d834f091acd95b912ea78cbd9b -->

# A fixed actual-length logarithmic-loss family: finite outer assembly

Status: complete outer arithmetic and evaluator specification CONDITIONAL
ON the exact finite affine-breaker contract (F) below. Its proposed derivation
in breaker-parameter-assembly.md is undergoing separate checking. This file
does not certify that derivation or claim external review, final evidence fit,
challenge intake, acceptance or submission. All constants below are fixed
once; none is advice depending on the input length.

This is an independently written composition of the finite seeded-linear,
ordinary-condenser, graph, affine-breaker and four-moment constructions
credited in their companion modules. The sumset framework is due to
Chattopadhyay--Liao, https://arxiv.org/abs/2110.12652v1 , with the improved
ordinary breaker following Xin Li,
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
length n, let L=lg(n+1), requiring L>=h+6. With primitive error exponent L
and leaf output 128uL, its output length is

  d=128*2^h*u*L,

its seed length is at most s=1200u^2 L, and its entropy requirement is at
most 3d. Its actual primitive guard K0<=n is required. Pad its seed to s by
ignoring an unused suffix. The sampler conversion at failure probability
1/100 uses entropy at most 3d+7. The construction is linear in its n-bit
data argument for every fixed seed and has a uniform polynomial evaluator.

(G) Use the real-entropy graph revision and the ordinary-condenser constants
at beta=1/4,tau=127/128,nu=1/128. These define positive integers M,D0,t0,
with row length v(t)=ceil(t/M), mixture entropy at least 63v/64 and error
at most 1/2 for all t>=t0, and kappa=1/M. Set b=4M. The graph on left
t=b*s bits, with right s bits and z=L, is proved under its literal guards

  t>=t0, v(t)>=64, L>=4, L>=lg(v(t)).

It has degree at most D0*(1+b*s/(4L)), and every set of at least 2^(t/4)
left vertices has strictly more than 2^(1-L)*2^s distinct neighbors.
The graph is a deterministic row-condenser followed by the finite field
Cayley walk and s-bit projection, not a random existence choice.

(F) There is one fixed positive computable integer C_F and a total
deterministic parameter/evaluation algorithm with this contract. Given
n,t,a,m,E>=1 it returns an exact seed length d_F, source entropy requirement
k_F and an m-bit generalized affine breaker, where

  d_F,k_F <= C_F*(t+1)^9*(m+a+E+lg(n+1)+1).              (F)

When k_F<=n, for arbitrary correlated tuples (A0,...,At) independent of
(B0,...,Bt,Y0,...,Yt), with H_inf(A0)>=k_F, Y0 uniform on d_F bits, and
alpha0 different from every alphaj, its honest output on Xj=Aj+Bj is
2^-E-close to uniform jointly with every tampered output and EVERY seed
Y0,...,Yt. The field is F_2, so addition is XOR. Neither source copies nor
seed copies are required to be deterministic functions of their honest copy.

For each fixed t and a,m,E=O(log n), the complete evaluator of (F) is
polynomial in the actual n. The candidate finite construction exports
C_F=2^90*(c_std+17)^9; the whole derivation, not just this number, must pass
its independent audit before (F) is regarded as discharged. No sharper
t*log(t)^3 estimate from an earlier conditional schedule is assumed.

(D) The finite fourth-moment lemma states the following. For a linear
somewhere sampler with A left indices and C choices, rare-test density
A^(-3/4), failure probability 1/100, output d and input entropy k_S, use
a degree-4C all-seeds-retained generalized affine breaker with raw error
at most e/2, where

  e=1/(800*C*A^9).

With source entropy

  k >= max(k_S, k_F+4C*d+lg(2/e)),                     (D)

XOR its C outputs per row using distinct pair advice, then take majority
over all A rows. Provided A>=2^20, each output bit occurs with probability
at least 99/5000 on the sum of any two independent sources of entropy k.
The companion proof includes the conditional Markov lift, consistent
g(W),f_i(W) bad-seed witnesses, all four-moment inequalities and constants.
None of the proof's distributions or witnesses is an evaluator input.

## 2. Fix every constant before seeing n

The upstream constants b,D0,t0,C_F are already fixed. Define

  cP=1200*b,
  cC=D0*(1+300*b),
  cT=4*cC+1,
  c0=lg(cC), c1600=lg(1600*cC),
  cQ=10*cP+c0+c1600+7,
  J=C_F*cT^9*cQ,
  h=max(512,2*lg(J+1)), u=h+1,
  B0=1200*u^2,
  T=128*2^h*u,
  p=b*B0,
  C=D0*(1+b*B0/4).

The division by four is integral because b is a multiple of four. Every
constant has an unambiguous finite algorithm. Astronomical constants need
not be written as decimal strings in the theorem: they are fixed integers
defined by terminating exponentiation and finite upstream parameter searches.
This is not an assertion that the evaluator is practical at small lengths.

For each actual n>=1 let L=lg(n+1), A=2^(pL), s=B0L and d=TL. Request
the breaker parameters

  t=4C, m=1,
  a=pL+lg(C),
  E=9pL+lg(1600C).                                       (1)

Then 2^-E<=1/(1600C A^9)=e/2, exactly as (D) requires. Since A is a power
of two and C is a fixed integer, encode a pair (alpha,j) by the pL-bit
alpha followed by the lg(C)-bit binary representation of j-1. This is
injective for j=1,...,C. If C=1 the suffix is empty; a is still positive.

## 3. The degree feedback inequality has a finite solution

For u>=2,

  C<=cC*u^2, t+1<=cT*u^2, p=cP*u^2,
  lg(C)<=c0+2lg(u), lg(1600C)<=c1600+2lg(u).

Using L>=1, equations (1) imply

  m+a+E+L+1
    <= [10p+lg(C)+lg(1600C)+3]*L
    <= cQ*u^2*L.                                       (2)

For the last inequality use lg(u)<=u and hence 4lg(u)<=4u^2; the
constant terms c0+c1600+3 can also be multiplied by u^2. Thus (F) gives

  d_F,k_F <= J*u^20*L.                                 (3)

For h>=512, 19log2(h+1)<=h/2. At 512 it follows from log2(513)<10,
giving the stronger bound 190<256. The derivative of the left side minus
the right side is 19/((h+1)ln 2)-1/2<0 thereafter; for example ln2>1/2
follows from its integral over [1,2]. The chosen h therefore gives

  u^19<=2^(h/2),   J<2^(h/2),
  J*u^20 < 2^h*u < T.

Consequently d_F,k_F<d for every L>=1. The exact computed d_F may depend
on n through its legitimate parameter selectors, but never exceeds the
sampler output length. Define the padded breaker on d-bit seeds by using
their first d_F bits.

All original d-bit seeds remain retained by the following consequence of
(F)'s universal right-side quantifiers. A simple kernel argument given only
the used seeds would not suffice when additive shifts correlate with suffixes.

Fix the admissible left tuple A. For an honest raw seed y and a deterministic
profile r=(B0,...,Bt,Y1,...,Yt), let Delta(y,r) be the TV between the output
tuple over A and an independent uniform honest output together with the
unchanged tampered-output tuple. For each y choose a maximizing profile r_y
from this finite domain. Uniform y with profile r_y is a legal (F) instance,
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

## 4. One fixed all-length guard and entropy coefficient

Put

  H=max(3T+7,(4C+1)*T+9p+lg(1600C)).                      (4)

H is a fixed positive integer. Since k_F<=TL, k_S<=3TL+7 and
lg(2/e)=9pL+lg(1600C), equations (3)--(4) imply that k=HL satisfies (D)
for every L>=1. Also H>=3T+7 covers the exact primitive entropy guard
K0<=n whenever HL<=n.

The evaluator tests these literal finite guards:

  n>=1; L>=h+6; HL<=n;
  pL>=t0; v(pL)>=64; L>=4; L>=lg(v(pL)).                 (5)

It also computes and checks the supplied algorithms' exact parameters and
guards. The already proved L-specialization, graph interface and (F)
guarantee that these checks pass when (5) holds. No test of an unknown
extractor/graph property, source distribution, or existence witness occurs.

There is one fixed threshold beyond which (5) holds for EVERY integer n.
Indeed L tends to infinity with n. Since v(pL)<=pL, it suffices to have

  L>=max(4,2lg(p)+2),

for then log2 L<=L/2 and log2 p<=L/2-1, so lg(pL)<=L-1<=L.
The guards pL>=t0 and v(pL)>=64 follow, for example, from
L>=ceil(t0/p) and L>=ceil(64M/p). The remaining fixed lower bounds on L
are automatic eventually. Finally H*lg(n+1)<=n eventually: choose j=floor(log2 n), so
2^j<=n<2^(j+1) and L<=j+2. The elementary exponential sequence 2^j/(j+2)
tends to infinity. This also proves eventual k_F<=n.

This is an all-integer-length argument, not a favorable subsequence or a
liminf used to invent a remainder. A small-input zero branch is allowed
because only an eventual lower bound is claimed.

## 5. The fixed evaluator

Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:

1. If n=0, output 0. Compute L and the parameters above. If any finite
   guard in (5) or its exact supplied parameter check fails, output 0.
2. Build the exact base extractor E_n from (L), seed-padded to s bits.
   Its data argument is x and its output length is d.
3. Build the (G) graph on left pL-bit strings, using z=L, the ordinary
   row-condenser and the canonical finite field walk. Project to s bits.
   Pad its degree by repeating an existing neighbor to exactly C.
4. For each alpha in lexicographic {0,1}^{pL} and j=1,...,C, compute

     y_(alpha,j)=E_n(x,Gamma(alpha,j)),
     b_(alpha,j)=F_n(x,Prefix(y_(alpha,j),d_F),encode(alpha,j)),

   using the ONE breaker parameter setting (1). Let R_alpha be the XOR
   of these C output bits.
5. Output the majority of all A bits R_alpha, breaking ties in favor of 1.

F_n in this pseudocode denotes the raw d_F-seed breaker; the preceding
proof uses its padded d-seed form. The two evaluations are identical.
No arbitrary seed, affine basis, irreducible polynomial advice, entropy
witness, good row, tampering, or source decomposition is supplied to this
algorithm. Constants, field bases, finite prime searches, orderings and
zero branches are fixed by the component algorithms.

The base sampler tests density epsilon0=2^-L. The graph yields rare-test
density epsilon1=A^(-3/4)<=epsilon0, because p>=4800, and at most
2A^(1/4) wholly bad rows except on source probability 1/100. Source
linearity survives graph index selection. The degree and all advice
strings match (D), so the fourth-moment proof applies literally.

## 6. Uniform polynomial time in actual input length

All one-time constant computations terminate, independent of n. Their
finite cost is part of the fixed machine constant, not length-dependent
advice. Parameter arithmetic has O(log n) bits beyond these fixed constant
lengths. There are A*C <= C*[2(n+1)]^p breaker calls, a fixed polynomial.
The base linear extractor is polynomial in n at fixed h. The graph's
ordinary-condenser argument is pL bits, so its 2^O(pL) evaluator is also
polynomial in n. The walk fields have fixed-algorithm arithmetic.

By (F), each breaker call is polynomial in n for the fixed t=4C and
a,E=O(log n). Its local prime searches and all finite field setup costs
2^O(log n) with fixed coefficients. Counting the majority uses an O(log A)
bit counter. Therefore the whole fixed total evaluator is polynomial in
the ACTUAL n. The polynomial degree and multiplier can be extremely large
but do not depend on n. No enumeration of extractor truth tables or
verification over all source distributions is hidden in this cost.

## 7. Affine flats, native B2 loss and expression meaning

For every guarded n, the output takes both values on the sum of any two
independent n-bit sources of entropy HL. Let a+V be any affine flat of
dimension at least HL. Choose X1 uniform on V and X2 independently
uniform on a+V. Both have that entropy, and X1+X2 is uniform on a+V.
There is no need to split the dimension in half. Thus the evaluator is
nonconstant on every such affine flat.

The self-contained native full-B2 bridge with d_flat=HL-1 gives

  C_fullB2(f_n) >= 4n-8H*lg(n+1)                         (6)

for every sufficiently large integer actual n. For n>=3,
2(n+1)<=n^2 implies lg(n+1)<=2log2 n; n=2 satisfies this directly.
Consequently (6) implies

  C_fullB2(f_n) >= 4n-16H*log2 n                        (7)

eventually at all integer lengths. This is a full 4n-O(log n) inequality
for one fixed P family, hence one fixed NP language with witness bound
zero: accept (x,empty) exactly when the evaluator outputs 1 and reject
every nonempty witness. Header and witness bits are not hidden inputs.

Once every premise is genuinely supported, expression-v2 can encode (7)
as an eventual lower bound with leading powers (1,0,0), coefficient 4,
and a `big_o` remainder of powers (0,1,0), with no coefficient field. The
fixed existential K can be 16H. The finite projection is (4,1), since
log n/n tends to zero. No declaration file or scoring claim is made here.

The outstanding premise in THIS outer derivation is the full exact (F)
construction and its independent audit; compact rewrite integrity, complete
offline dependency closure and final atomic evidence constraints also
remain to be checked before any final-candidate claim.
