# CL-v1: at four explicit witness states, degree-class weights cannot lift the ternary-pair step accounting above 4

This is an **unscored barrier lemma**. It asserts no circuit lower bound and
no upper bound on the circuit size of any function. It is a finite statement
about one per-step accounting scheme for one gate-elimination procedure,
evaluated at four explicit states.

## 0. Statement

Let R be the ternary-pair procedure defined in Section 1.

Let R' be any variant of R that differs only at R's minimal-AND step. That
step chooses G, the first gate in list order. When both parents of G are
unpaired parameters, R' may, instead of GP, fix one of G's two parents to 0
or 1 and then exhaust D1 and D2. That step ends there. D3, D4 and M belong to
the next iteration.

Consider every potential of the form

    W = s + sum over unpaired parameters u of k(min(deg u, 3))
          + sum over pairs {p,r} of pi(min(deg p, 3), min(deg r, 3)),

where k(0) = 0, the values k(1), k(2), k(3) are arbitrary reals, pi is an
arbitrary real symmetric function on {0,1,2,3}^2, and kmax = max(k(1), k(2), k(3)).

**Lemma.** Let S1, S2, S3 be the states of Section 2, and S1' the state that
R reaches from S1 (Section 3). Suppose a real delta satisfies

    Delta W <= delta * Delta a

for the step R takes at each of S1, S1', S2 and S3. Then delta - kmax <= 4.

The same conclusion holds for the steps R' takes at those four states,
whatever R' chooses at S1 and S2:

- If R' fixes a parent of G at S2, to either value, the four inequalities
  give delta - kmax <= 11/3, whatever R' does at S1. This also holds if the
  M and D4 performed in the next iteration after a 1-fixing at S2 are
  counted as part of that step.
- If R' fixes a parent of G at S1 to its non-controlling value 1, the S1
  inequality alone gives delta - kmax <= 3.

For k(1) = k(2) = k(3) = 5, pi = 3 and delta = 9, all four inequalities for R
hold with equality and delta - kmax = 4.

**Motivation, not part of the lemma.** In the summation argument of the
ternary-pair method, the per-step inequality is summed over a run that ends when
a <= d+1. With nonnegative weights and delta >= 0, W_0 <= |C| + kmax*n and
W >= 0, so the
sum gives |C| >= (delta - kmax)*n - delta*(d+1). Thus delta - kmax is the
coefficient such an accounting certifies. A potential g*s + ... with g > 0
reduces to g = 1 after dividing every weight by g.

## 1. The framework and the procedure R

These definitions follow the ternary-pair method in CircuitBreaker
submission a6a363c8, file elimination.md Sections 1-6; provenance is in
Section 6. Only what the witness states need is stated, and everything is
checkable from this file alone. At each witness state the first applicable
row is among those stated here, so R's step there is determined. At S1'
the only freedom is P2's symmetric choice between x and y, which Section 3
treats.

**States.** A state consists of:

- parameters u_1, ..., u_m;
- a set of disjoint marked pairs. Each pair {p,r} has an excluded corner
  (e_p, e_r), and its allowed set is the other three points of F_2^2;
- a fair affine system A v = N u + c of auxiliaries v, with A square and
  invertible and at most three distinct nonconstant terms per equation;
- an acyclic list of B2 gates, each reading parameters, auxiliaries, earlier
  gates or constants. No gate feeds an equation;
- an output, which is a literal of an internal value.

We use the identity map from parameters to inputs; it plays no role below.
Other notation:

- s is the number of equations plus gates.
- deg w is the number of distinct gates and equations that read w. An
  equation reads every value occurring in it, including its auxiliary. Pins
  and the output reference are not counted.
- a = m - q, where q is the number of pairs.
- The safe value of a paired p is 1 - e_p. Fixing p to it leaves both values
  of its partner allowed.

Write AND(u,v) = u*v and XOR(u,v) = u+v over F_2. Each AND below has
exceptional corner (1,1), majority value 0 and controlling value 0 in each
argument.

**Deletions and preparation**, exhausted at the start of every iteration:

- **D1.** Delete a gate that has a constant pin, a repeated pin, a constant or
  unary table, or is AND-type on two affine parents with dependent linear
  parts. Reroute every use of it to the resulting constant or literal. A
  negation is absorbed into the reader's table.
- **D2.** Delete a non-output gate with no reader.
- **D3.** Delete an equation with at most two nonconstant terms, substituting
  the literal it defines.
- **D4.** If an unpaired parameter x is read only by an equation E, exchange x
  with an auxiliary Q that depends on x: Q becomes a parameter and x an
  auxiliary. Then delete x and E.
- **M.** A gate of XOR type whose two inputs are parameters or auxiliaries
  becomes one equation, with a new auxiliary for its value. This changes
  neither s nor the degree of any parameter.

**Paired rows.** The first applicable row is applied, in this order:

- **P0.** A paired parameter of degree 0: fix its partner safely.
- **P2.** A paired p of degree at least 2: fix p safely. D1 then deletes its
  readers.
- **PC.** Partners with a common reader: restrict the pair to an allowed line
  on which that reader is constant.
- **PE.** A paired p of degree 1 whose reader is an equation E: fix its partner
  r safely, delete r's reader H, then apply D4 to p and E.

**Minimal-AND step.** Otherwise let G be the first gate in list order. If a
parent is an auxiliary, exchange it with a parameter. If both parents x, y
are unpaired, apply GP:

- **GP.** Add the pair {x,y}, excluding the exceptional corner of G's table.
  On the new domain G is constant. Delete G, and then, by D1, every reader of G.

The rows for paired parents (AP1, AP2) are not reached by the witness states.

**Change of W.** A step's Delta W counts:

- every deletion the step makes, including its D1 and D2 closure;
- every change of degree class of an unpaired parameter;
- every new or removed pair, and every change of a pair's profile;
- the charge k of a parameter that leaves or enters the unpaired set.

## 2. The witness states

For an opaque gate h_i we write h_i = AND(t, t'), with two fresh parameters
t, t' that occur nowhere else. No step below touches them, and their degrees
never change. Gates are listed in their list order.

**S1.** The parameters are x, y and the t's. The gates are:

    G = AND(x,y);
    h1, h2, h3, h4, h5 (opaque);
    A = XOR(x,h1), B = XOR(x,h2), C = XOR(y,h3), D = XOR(y,h4), J = XOR(G,h5);
    O1 = AND(A,B), O2 = AND(C,D), O3 = AND(O1,O2), OUT = AND(O3,J).

The output is OUT. There are no pairs and no equations. Degrees: x is read by
G, A, B and y by G, C, D, so deg x = deg y = 3.

**S2.** The parameters are x, y, z and the t's. The gates are:

    G = AND(x,y);
    h1, h3, h6, h7, h8 (opaque);
    A = XOR(x,h1), C = XOR(y,h3), J = XOR(G,z);
    J1 = AND(J,h6), J2 = AND(J,h7), J3 = AND(J,h8);
    O1 = AND(A,C), O2 = AND(J1,J2), O3 = AND(O2,J3), OUT = AND(O1,O3).

The output is OUT, with no pairs and no equations. Degrees: deg x = 2 (G, A),
deg y = 2 (G, C) and deg z = 1 (J).

**S3.** The parameters are p, r, z and the t's. There is one pair {p,r}, with
excluded corner (1,1). There is one auxiliary v and one equation

    E: v + p + z = 0,

so A = (1) is invertible. The gates are:

    ha, hb, hc, hd, he, hf, hg (opaque);
    H = XOR(r,ha);
    V1 = AND(v,hb), V2 = AND(v,hc), V3 = AND(v,hd);
    Z1 = AND(z,he), Z2 = AND(z,hf), Z3 = AND(z,hg);
    O1 = AND(H,V1), O2 = AND(V2,V3), O3 = AND(Z1,Z2);
    O4 = AND(O1,O2), O5 = AND(O3,Z3), OUT = AND(O4,O5).

The output is OUT. Degrees: deg p = 1 (E), deg r = 1 (H), deg z = 4
(E, Z1, Z2, Z3) and deg v = 4.

In each state all invariants hold: the system is fair and invertible or
empty, the gates are acyclic, no gate feeds an equation, and the output is a
gate. They are syntactic states. They are not claimed to arise from a circuit
for an affine disperser.

**No earlier row applies.**

- In S1 and S2 no gate has a constant or repeated pin, and every table
  depends on both inputs. Every AND on two affine parents (G and the opaque
  gates) reads two distinct parameters, so its parents are independent.
  Every non-output gate has a reader. There is no equation, so D3 and D4 do
  not apply. Every XOR gate has a gate input, so M does not apply. There are
  no pairs.
- In S3 the same checks pass. V1-V3 and Z1-Z3 each have a gate parent, so
  D1's affine-parents clause does not reach them. E has three terms. D4 does
  not apply: z has degree 4, p and r are paired, and the t's are read only by
  gates. M does not apply because H has a gate input. P0 and P2 fail because
  p and r have degree 1. PC fails because p's reader E differs from r's reader H. So PE is
  the first applicable row.

In S1 and S2 the minimal AND is G, whose parents x, y are unpaired
parameters, so R applies GP.

## 3. The steps and their exact changes

**S1, GP.** Pair {x,y}, excluding (1,1). G is 0 on the domain, so G is
deleted. J = XOR(0,h5) is a literal of h5, so J is deleted and OUT reads h5.
No other gate gets a constant pin or loses its last reader, so Delta s = -2.

Before the step, x and y were unpaired with degree 3, contributing 2*k(3).
Afterwards they are paired, with readers A, B and C, D: profile (2,2), charge
pi(2,2). Delta q = 1, so

    Delta a = -1,   Delta W = -2 - 2k(3) + pi(2,2).

The resulting state is S1'.

**S1', P2.** S1' has no deletion, M or P0, and x has degree 2. Fixing x to
its safe value 0 makes A = h1 and B = h2, so D1 deletes both and
Delta s = -2. Now y is unpaired with degree 2, and the pair is gone.
Delta m = Delta q = -1, so

    Delta a = 0,    Delta W = -2 + k(2) - pi(2,2).

If R picks y, the swap x<->y, A<->C, B<->D, h1<->h3, h2<->h4, O1<->O2 gives
the same values.

**S1, the fixings of R'.**

- x := 0: G = 0, and then J = h5, A = h1, B = h2; all four are deleted.
  x leaves, and y drops to degree 2. Delta a = -1 and
  Delta W = -4 - 2k(3) + k(2).
- x := 1: G = y, so G is deleted and J reads y. A and B become negated
  literals of h1 and h2 and are deleted. J still has the gate input h5, so
  no M follows. y keeps degree 3 (C, D, J). Delta a = -1 and
  Delta W = -3 - k(3).

Fixing y is symmetric.

**S2, GP.** Pair {x,y}. G is deleted, and J = XOR(0,z) = z is deleted, so
J1, J2, J3 now read z. Delta s = -2. Afterwards x and y have degree 1 each
(A and C), so the profile is (1,1). z rises from degree 1 to degree 3. So

    Delta a = -1,   Delta W = -2 - 2k(2) + pi(1,1) + k(3) - k(1).

**S2, the fixings of R'.**

- x := 0: G = 0, then J = z and A = h1; G, J and A are deleted. x leaves,
  y drops to degree 1 and z rises to degree 3. Delta a = -1 and
  Delta W = -3 - k(2) + (k(1) - k(2)) + (k(3) - k(1)) = -3 - 2k(2) + k(3).
- x := 1: G = y, so G is deleted and J = XOR(y,z). A is deleted. x leaves.
  y keeps degree 2 (C, J) and z keeps degree 1. With D1 and D2 closure
  only, Delta a = -1 and Delta W = -2 - k(2).

  The next iteration first applies M: J becomes the equation
  v_J + y + z = 0. Then z is unpaired and read only by that equation, so D4
  exchanges z with v_J, the only auxiliary depending on z, and deletes z and
  the equation. M changes nothing in W. D4 gives:
  - -1 to s, for the equation;
  - k(1) - k(2), as y drops to degree 1;
  - -k(1), as z (class 1) leaves the unpaired set;
  - +k(3), as v_J enters it with degree 3 (J1, J2, J3).

  That totals -1 - k(2) + k(3). Counting this M and D4 inside the step gives
  Delta a = -1 and Delta W = -3 - 2k(2) + k(3), the same as the 0-fixing.

Fixing y is symmetric.

**S3, PE.** Fix r to its safe value 0. Then H = XOR(0,ha) = ha, so H is
deleted and O1 reads ha. Now p is unpaired and read only by E, so D4 applies.
The only auxiliary is v = p + z, which depends on p. It becomes a parameter
w, and p becomes an auxiliary; then p and E are deleted. So Delta s = -2
(H and E).

- The pair {p,r}, of profile (1,1), disappears.
- w is a new unpaired parameter read by V1, V2, V3, so it has degree 3.
- z drops from degree 4 to degree 3, both in class 3.
- m falls by one (r leaves, w replaces p) and q falls by one.

So

    Delta a = 0,    Delta W = -2 - pi(1,1) + k(3).

These values are exactly those recorded by the audit experiment.

## 4. Proof of the lemma

Put c = delta - kmax. Each Delta a = -1 row reads Delta W <= -delta, and each
Delta a = 0 row reads Delta W <= 0. Each inequality below uses kmax >= k(i)
for the stated i.

**(1) S1, using GP or a 0-fixing.** For GP, adding the S1 inequality and the
S1' inequality gives -4 - 2k(3) + k(2) <= -delta. The x := 0 and y := 0
fixings give this directly. With kmax >= k(3),

    c <= 4 + k(3) - k(2).

**(1') S1, using a 1-fixing.** delta <= 3 + k(3). With kmax >= k(3), c <= 3.

**(2) S2, using GP.** delta <= 2 + 2k(2) - pi(1,1) - k(3) + k(1). With
kmax >= k(1),

    c <= 2 + 2k(2) - k(3) - pi(1,1).

**(2') S2, using a 0-fixing, or a 1-fixing counted with the next M and D4.**
delta <= 3 + 2k(2) - k(3). With kmax >= k(3), c <= 3 + 2k(2) - 2k(3).

**(2'') S2, using a 1-fixing with D1 and D2 closure only.** delta <= 2 + k(2).
With kmax >= k(2), c <= 2.

**(3) S3, PE.**

    pi(1,1) >= k(3) - 2.

Each choice at S1 and S2, under either accounting of a 1-fixing at S2, is
covered by one of four cases:

- **S1 uses a 1-fixing.** Then c <= 3 by (1').
- **S1 uses GP or a 0-fixing, and S2 uses a D1/D2-closed 1-fixing.** Then
  c <= 2 by (2'').
- **S1 uses GP or a 0-fixing, and S2 uses a 0-fixing or a 1-fixing with the
  next M and D4.** Twice (1) plus (2') gives
  3c <= 8 + 2k(3) - 2k(2) + 3 + 2k(2) - 2k(3) = 11, so c <= 11/3.
- **S1 uses GP or a 0-fixing, and S2 uses GP.** This is the case of R
  itself. Twice (1), plus (2) and (3), gives

      3c <= 8 + 2k(3) - 2k(2) + 2 + 2k(2) - k(3) - pi(1,1)
         =  10 + k(3) - pi(1,1)
         <= 12,

  so c <= 4.

So c <= 4 for R and for every choice of R'. Any fixing at S2 gives at most
11/3, and a 1-fixing at S1 gives at most 3.

**Attainment.** For k(1) = k(2) = k(3) = 5, pi = 3 and delta = 9, R's four
values of Delta W are -9, 0, -9 and 0, at Delta a = -1, 0, -1 and 0. Every
inequality is an equality, and c = 9 - 5 = 4. QED.

No step of this proof uses a sign condition on the weights, the audit, or any
external theorem.

## 5. Scope and relevance

**Relevance.** rules.md allows a relevant partial lemma as an unscored string
claim. scientific-cases.md lists "The unchanged LP optimum" as a valid
unscored lemma: a finite ceiling certificate for the seed's accounting. This
lemma is of the same kind, for the ternary-pair accounting of Section 1,
whose flat weights (k = 5, pi = 3, delta = 9) meet these rows with equality.

Any degree-class reweighting of the per-step form above that certifies more
than 4 must do one of two things:

- change one of the four steps (GP at S1 and S2, P2 at S1', or PE at S3)
  other than by the fixings of R', which Section 4 shows also stay at or
  below 4;
- or exclude these states.

**Context, not used in the proof.** Two earlier remarks consider related
accountings with one flat charge for every used unpaired parameter:

- the table lemma in a6a363c8 Section 8, for W = g*s + i*I + p*q and this
  procedure's row bounds;
- Remark 4.3 of the Mishra-Santhanam paper, about their closely related
  procedure, where the potential s + 4q + 6I is stated to be optimal.

The present lemma lets the parameter charge depend on its degree class and
the pair charge on both members' degree classes. It also allows the fixing
variant R'. Nothing here relies on those remarks.

The exact trade-off for R, with t = k(3) - k(2):

- (1) gives c <= 4 + t;
- (2) and (3) give c <= 2 + 2k(2) - k(3) - (k(3) - 2) = 4 - 2t.

Raising t gains t at S1 but loses 2t at S2. Twice the first bound plus the
second gives 3c <= 12 for every t.

The lemma does **not** exclude:

- analyses that use properties of reachable states beyond the syntactic
  invariants;
- other rules for selecting the minimal AND. For example, choosing an
  opaque gate first at S1 avoids this witness;
- amortized or multi-step accountings not of the per-step form above;
- potentials with other terms, for example ones that look at the gate inputs
  of XOR readers;
- other procedures, other proofs, or a coefficient above 4 by any other
  method.

It makes no claim about the correctness of the 4n theorem of a6a363c8 or
5e9fcf7f. It only uses the procedure's definitions as stated in Section 1.

The audit experiment rebuilds S1, S2 and S3 from their gate lists. It:

- checks the syntactic conditions showing that no earlier row applies;
- applies each stated step with constant and literal propagation, including
  the M and D4 after the S2 1-fixing;
- compares every change of a potential term with Section 3;
- re-derives each choice's bound in exact rational arithmetic, and checks the
  attainment equalities.

It is supplementary: the proof above does not depend on it.

## 6. Provenance

- **Procedure and framework:** CircuitBreaker submission a6a363c8 (0xmoei),
  elimination.md Sections 1-6. It reproduces, with audit clarifications,
  submission 5e9fcf7f (ArthurOstapenko et al.).
- **Parameter method:** Om Swostik Mishra and Rahul Santhanam, "An Explicit
  3.5n - o(n) Size Lower Bound against General Boolean Circuits", ECCC Report
  TR26-260 (2026). Its revision 1, retitled "An Explicit 4n - O(log(n)) Size
  Lower Bound ...", contains Remark 4.3.
- **Contributed here:** the degree-class potential family, the witness
  states S1, S1', S2, S3, the exact step computations and the case
  certificates.

Written with Claude Opus 5.5 in Claude Code.
