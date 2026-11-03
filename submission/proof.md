# A complete supplied proof of C_B2(F_N) >= 4N - o(N) for one fixed total NP family

This document is self-contained evidence for the declaration in `bound.json`. Part A
defines one total language L in P (hence in NP with empty witnesses) and gives a
single polynomial-time evaluator for all lengths. Part B records the affine-disperser
premise for the payload, with its exact citation and the parameter arithmetic
verified here. Part C proves from scratch the class theorem that turns that premise
into a lower bound. Part D assembles the numerical remainder, the quantifiers and the
explicit thresholds. Part E states exactly which premises are supplied here and which
are taken from published sources.

The score claim is the finite pair `[4, 1]`: for the one family below, at actual input
length N and in the full B2 model of the challenge rules,

    C_B2(F_N) >= 4N - 20*N^(3/5) - 160*N^(4/5) for every integer N >= 2^20,

so liminf_{N->infinity} C_B2(F_N)/N >= 4. No endpoint, no exponent above one, no
separation of complexity classes and no optimality of the constant is claimed.

---

## 1. Model, notation and bookkeeping conventions

**1.1 Circuits.** We use the full B2 model fixed by `rules.md`: deterministic acyclic
circuits with one output, gates of fan-in two with ordered input pins, all sixteen
binary Boolean operations available, unrestricted depth and fan-out, repeated pins
allowed and shared subgraphs allowed, exact computation required on every input, and
inputs x alone (never NP witnesses). Size is the number of retained gates; input
vertices, wires and output selection are free. Constant nodes introduced during
substitutions are free. A negation on a wire is absorbed into the operation of the
gate that reads it, so a negated projection costs nothing, while a *retained*
negated-projection gate (a gate labelled with a unary or constant operation) counts
as one.

**1.2 Field conventions.** F_2 = {0,1} with addition XOR and multiplication AND, so
negation is addition of 1. Every operation h: F_2^2 -> F_2 has a unique multilinear
form h(t1,t2) = c + b1 t1 + b2 t2 + b12 t1 t2. The XOR-type operations (two of the
sixteen) have the form t1+t2+c; the AND-type operations (eight of the sixteen) have
the form (t1+a)(t2+b)+c. The two constant operations and the four unary operations
complete the sixteen. An operation is *affine* if b12 = 0.

Throughout, log means log base 2, ceil/floor are written explicitly, and all rational
constants are exact.

---

## 2. Part A: the fixed total family

**2.1 Length parameters.** Fix N >= 0. If N < 2^20 set F_N identically 0. Otherwise
let

    k = the least prime with k^5 > 32 N^3,      m = N - k,
    r = ceil(m/k),                              t = ceil(sqrt(r)).

*Lemma A1 (parameter bounds).* For every N >= 2^20: 2N^(3/5) < k < 5N^(3/5),
N/2 < m < N, 3 <= r <= N^(2/5), 2 <= t <= 2N^(1/5), t < r, and k >= 2 r(r+1)/t.

*Proof.* k^5 > 32N^3 iff k > 2N^(3/5) (both sides positive: 32^(1/5) = 2). By
Bertrand's postulate, applied to the integer n_0 = ceil(2N^(3/5)), there is a prime in
(n_0, 2n_0], so

    k <= 2 n_0 = 2 ceil(2N^(3/5)) <= 4N^(3/5) + 2 < 5N^(3/5),

the last inequality holding for N^(3/5) >= 2, i.e. N >= 3. Hence 2N^(3/5) < k < 5N^(3/5).
Since 5N^(3/5) <= N/2 iff N^(2/5) >= 10 iff N >= 10^(5/2) < 317, we get
m = N - k > N/2 for N >= 2^20. Then m/k < N/(2N^(3/5)) = N^(2/5)/2, so
r = ceil(m/k) <= N^(2/5)/2 + 1 <= N^(2/5) whenever N^(2/5) >= 2, and
r >= m/k >= (N/2)/(5N^(3/5)) = N^(2/5)/10 >= 25 for N >= 2^20, so 3 <= r. Next
sqrt(r) >= 5 and t = ceil(sqrt(r)) <= sqrt(r) + 1 <= 2 sqrt(r) <= 2N^(1/5), and t < r
because ceil(sqrt(r)) < r for r >= 3. Finally put u = N^(2/5) >= 256. From r <= u/2+1
and r+1 <= u/2+2 we get r(r+1) <= (u/2+1)(u/2+2) <= (u/2+2)^2 <= u^3 for u >= 8, and
t >= sqrt(r), so

    r(r+1)/t <= r(r+1)/sqrt(r) = sqrt(r)(r+1) <= (u/2+2)^(3/2) <= u^(3/2) = N^(3/5),

while k > 2N^(3/5) = 2u^(3/2). Hence k >= 2 r(r+1)/t. QED

**2.2 Definition of F_N.** Let x be an N-bit string, N >= 2^20, written as
x = (q_0,...,q_{k-1}, p_0,...,p_{m-1}) with the first k bits first. Let

    Q(X) = X^k + q_{k-1}X^(k-1) + ... + q_1 X + q_0  in F_2[X]

be the monic degree-k polynomial whose lower coefficients are read off the first k
bits, and let K = F_2[Z]/(Q) with the polynomial basis 1, Z, ..., Z^(k-1).

If Q(1) = 0, or if X^(2^k) is not congruent to X modulo Q, output F_N(x) = 0.

Otherwise Q divides X^(2^k) - X, so by Lemma A2 (with k >= 5 here, since
k > 2N^(3/5) >= 2^13 for N >= 2^20) Q is irreducible, K is the field with 2^k
elements, and
we proceed. Split the payload bits p_0,...,p_{m-1} in order into r consecutive blocks
of k bits, zero-padding the final block if m < rk (this happens only in the last
block, since r = ceil(m/k)); read block i as the element

    z_i = p_{(i-1)k} + p_{(i-1)k+1} Z + ... + p_{ik-1} Z^(k-1)  in K,

where a padded bit is read as 0. Put

    y_i = z_i ^ (2^(r+2-i) - 1)  in K,   for 1 <= i <= r,

let e_t(y_1,...,y_r) = sum over the t-element subsets I of {1,...,r} of prod_{i in I} y_i
be the degree-t elementary symmetric polynomial in K, and output

    F_N(x) = Tr_{K/F_2}( e_t(y_1,...,y_r) )  in {0,1},

where Tr is the field trace. This defines F_N on every N-bit string, and L = {x :
F_|x|(x) = 1} is a total language. The family is fixed: no per-length advice, no
oracle, no nonconstructive choice. All lengths, including N < 2^20 and the case
k >= N, are covered by the definition.

*Remark on the exponent order.* The exponents 2^(r+2-i)-1 for i = 1..r are the
descending order of the geometric sums 1+2+...+2^(e_i), e_i = r+1-i, used in the
cited disperser theorem of Part B; reversing the block index converts them to the
ascending order 2^(i+1)-1, and Proposition 3.3 of that source identifies the two
presentations. Because our r, t are the source's r, t (r = ceil(m/k), t = ceil(sqrt r))
and the map from payload bits to (z_1,...,z_r) is a fixed injective F_2-linear map,
the two conventions define the same function up to relabelling of blocks.

*Lemma A2 (irreducibility test).* Let k >= 3 be prime and let Q in F_2[X] be monic of
degree k. Then Q divides X^(2^k) - X if and only if Q is irreducible.

*Proof.* (=) If Q is irreducible of degree k, then F_2[X]/(Q) is the field with 2^k
elements, whose multiplicative group has order 2^k - 1; every element a satisfies
a^(2^k) = a, so X^(2^k) = X (mod Q).
(=) The polynomial Y^(2^k) - Y is squarefree (its derivative is 1) and its roots are
exactly the elements of F_(2^k); equivalently, it is the product of all monic
irreducibles over F_2 whose degree divides k, each occurring once. Since k is prime,
those degrees are 1 and k. Write Q = prod_j Q_j, the factorization into monic
irreducibles, with deg Q_j = d_j in {1,k} and sum_j d_j = k. The monic irreducibles of
degree 1 are X and X+1, so at most two indices have d_j = 1. If some d_j = 1 and some
other index has d_j = k, then k = sum_j d_j >= k + 1, a contradiction; if exactly one
d_j = 1 and all others are k, then k = 1 + (multiple of k) is impossible for k >= 3; if
exactly two have d_j = 1 and there are no degree-k factors, then k = sum_j d_j <= 2 < k,
a contradiction. Hence every d_j = k, and sum_j d_j = k forces a single factor: Q is
irreducible. QED

*Remark.* No hypothesis on Q(0) or Q(1) is needed. Consequently the evaluator's extra
rejection of Q(1) = 0 changes nothing: a monic irreducible of degree k >= 2 has
Q(0) = 1 and Q(1) = 1 (otherwise X or X+1, both of degree 1, would divide it), so it
passes both tests, while X^(2^k) = X (mod Q) implies irreducibility by the lemma.
Accepted prefixes are therefore exactly the irreducible ones, and rejected prefixes
simply output 0; no field structure is ever assumed on a rejected prefix.

**2.3 One polynomial-time evaluator.**

*Lemma A3.* There is a single deterministic algorithm that, on input x in {0,1}^N,
computes F_N(x), and runs in time polynomial in N.

*Proof.* Compute k by scanning candidates k = k_0, k_0+1, ... from k_0 = floor(2N^(3/5))+1
and testing k^5 > 32N^3 by integer arithmetic and k by trial division up to sqrt(k).
Only O(N^(3/5)) candidates are scanned (Lemma A1), each test costs O(N^(2/5)) bit
operations, so this is polynomial. Irreducibility: compute X^(2^k) mod Q by k
successive squarings, each a squaring of a polynomial of degree < k followed by
reduction modulo Q, costing O(k^2) bit operations, hence O(k^3) = O(N^(9/5)) in total;
compare with X. Powers y_i = z_i^(2^(r+2-i)-1): the exponent has bit length r+2-i <= r,
so repeated squaring uses at most O(r) multiplications in K, i.e. O(r k^2) bit
operations per block; summed over r blocks this is O(r^2 k^2) = O(N^2) bit operations.
Elementary symmetric polynomial without division by integers: start E_0 = 1 and
E_j = 0 for 1 <= j <= t, and for i = 1..r update E_j <- E_j + y_i E_(j-1) in descending
order of j from min(i,t) down to 1. Splitting the t-subsets by whether they contain i
proves the recurrence; it uses O(r t) = O(N^(3/5)) field multiplications, each
O(k^2) = O(N^(6/5)) bit operations. The trace is Tr(a) = a + a^2 + ... + a^(2^(k-1)),
computed in O(k) squarings in K. All steps are deterministic and polynomial in N;
the field operations use only the polynomial basis and the reduction rule Z^k = Q + X^k
(i.e. Z^k = q_(k-1)Z^(k-1)+...+q_0). QED

**2.4 NP membership, both directions.** Let A be the algorithm of Lemma A3. Define the
verifier V(x,w) = 1 iff w is the empty string and A(x) = 1, with witness bound
W(N) = 0. Then
  (i) if x in L then w = empty satisfies |w| <= 0 and V(x,w) = 1;
  (ii) if V(x,w) = 1 for some w with |w| <= 0 then w is empty and A(x) = 1, i.e. x in L.
Both directions hold, V runs in time polynomial in |x| (Lemma A3 plus a constant-size
check), so L is in P and hence in NP. Witness bits are not circuit inputs; the target
circuits of the lower bound receive x alone. This also discharges the eligibility
requirement: one fixed polynomial-time evaluator and one fixed polynomial witness
bound work for every length, and the language and algorithms do not change across the
asymptotic quantifiers.

**2.5 The hard restriction.** Fix once and for all an irreducible monic polynomial
Q* of degree k with Q*(0) = Q*(1) = 1; such a polynomial exists for every k >= 1
(e.g. take any irreducible of degree k over F_2, which exists by the standard counting
argument, and note that an irreducible of degree k >= 2 has both Q*(0) = 1 and
Q*(1) = 1, since X and X+1 have degree 1). Let x* be the N-bit string whose first k
bits are the lower coefficients of Q* and whose payload bits are zero, and define the
*payload function*

    g_N : {0,1}^m -> {0,1},  g_N(p) = F_N(x* with payload p).

Then g_N is exactly the function of Part B (the prefix is valid, so the evaluator
never returns 0 through the prefix test), and g_N is obtained from F_N by fixing k
input bits to constants.

*Lemma A4 (restriction).* If a B2 circuit C computes F_N on {0,1}^N, then there is a
B2 circuit C' with |C'| <= |C| that computes g_N on {0,1}^m.

*Proof.* Substitute the k fixed bits into C: every input vertex carrying a fixed bit
becomes the corresponding constant, and every gate that reads a constant is retained
but now computes a unary or constant operation (or, if both pins are the same
constant, an operation on one value); no gate is created and no gate needs to be
deleted, so |C'| <= |C|. QED

---

## 3. Part B: the payload is an affine disperser

**3.1 Definition and the cited theorem.** A function f: F_2^n -> F_2 is an *affine
disperser for dimension d* if it is non-constant on every affine subspace of F_2^n of
dimension at least d. The cited theorem below delivers exactly this monotone form (it
gives non-constancy on every affine subspace whose dimension exceeds a stated
threshold). Under the weaker non-monotone reading, in which only subspaces of dimension
exactly d are required, every step of Part C still goes through: an affine subspace S of
dimension at least d contains an affine subspace of dimension exactly d, and
non-constancy on that subspace already contradicts C_0 being constant on all of S.

*Theorem B1* (Ben-Sasson and Kopparty, "Affine Dispersers from Subspace Polynomials",
author preprint dated 2010-03-12, Theorem 7.3, printed p. 34; special case p = 2).
Let t < r be integers and let n be a prime with n >= r(r+1)/t. For i in [r] put
e_i = r+1-i and theta_i = 1 + 2 + 4 + ... + 2^(e_i) = 2^(e_i+1) - 1. Let
gamma in F_(2^n) be nonzero, and let f: (F_(2^n))^r -> F_2 be

    f(x_1,...,x_r) = Tr( gamma * sum_{I subset [r], |I| = t} prod_{i in I} x_i^(theta_i) ).

Let A be any F_2-affine subspace of (F_(2^n))^r with

    dim A > n r / t + n t + r (r+1).

Then |f(A)| > 1.

The same source's Theorem 2.2 states the resulting disperser with the canonical
parameter choice ("n the smallest prime bigger than 2 m^(3/5), r = ceil(m/n),
t = ceil(sqrt r)") and concludes non-constancy on every affine space of dimension
greater than 6 m^(4/5). We use the more general Theorem B1 with our own explicit
parameters, because our k is the least prime above 2N^(3/5) rather than the least
prime above 2m^(3/5), and we want the bound to be a statement about the literal
parameters of Part A.

**3.2 Verification of the hypotheses for our parameters.** Take n = k, r, t as in
section 2.1, gamma = 1, and identifiers K = F_(2^k) = F_2[Z]/(Q). By Lemma A1, k is
prime, t < r, k >= 2 r(r+1)/t (so in particular k >= r(r+1)/t), and with
u = N^(2/5) >= 256,

    nr/t  <= k * sqrt(r)      <= 5N^(3/5) * N^(1/5)   =  5 N^(4/5),
    nt    <= k * 2 sqrt(r)    <= 5N^(3/5) * 2N^(1/5)  = 10 N^(4/5),
    r(r+1)<= (u/2+1)(u/2+2)   <= u^2                  =        N^(4/5),

where we used r/t <= r/sqrt(r) = sqrt(r) <= u/2 + 1 <= u = N^(2/5) and
t = ceil(sqrt r) <= 2 sqrt(r) <= 2N^(1/5). Hence the theorem's dimension threshold
satisfies

    d(N) := ceil( n r / t + n t + r (r+1) ) <= 16 N^(4/5).

Therefore, for every N >= 2^20, the function f of Theorem B1 with these parameters is
non-constant on every affine subspace of (F_(2^k))^r of dimension at least d(N), i.e.
it is an affine disperser for dimension d(N).

**3.3 Transfer to the payload.** Let phi: {0,1}^m -> (F_(2^k))^r be the map that reads
the payload as r consecutive k-bit blocks in the polynomial basis (zero-padding the
last block). phi is F_2-affine and injective, and g_N = f o phi by construction of
section 2.2. If A is an affine subspace of {0,1}^m of dimension at least d(N), then
phi(A) is an affine subspace of (F_(2^k))^r of the same dimension, so f is
non-constant on phi(A), hence g_N is non-constant on A. Thus

    g_N is an affine disperser for dimension d(N) <= 16 N^(4/5), for every N >= 2^20.

**3.4 What this part does and does not supply.** The statements of Theorem B1 and of
the source's Theorem 2.2 are quoted exactly; the parameter verification, the
identification of our function with the source's (Remark 2.2, Proposition 3.3 of the
same source) and the transfer to the payload are proved here. The internal proof of
Theorem B1 (the source's Sections 3, 6 and 7: subspace polynomials, the main
structural lemma, and the monomial-isolation argument for the composed polynomial) is
not reproduced in this file; it is a published, refereed argument of the cited source,
and it is not a computational claim. Section E states this scope limitation again.

---

## 4. Part C: the class theorem (complete proof)

This part proves, from scratch:

*Theorem C0 (class theorem).* Let f: F_2^n -> F_2 be an affine disperser for dimension
d, computable in polynomial time. Then every B2 circuit C_0 on n inputs that computes
f has size |C_0| >= 4n - 10d.

The proof is a gate-elimination argument carried out on an underlying space of
independent parameters, following the architecture of Mishra and Santhanam, "An
Explicit 4n - O(log n) Size Lower Bound against General Boolean Circuits", ECCC
TR26-260, revision 1 (2026-11-03), Sections 2-4. Every case below is checked here with
explicit bookkeeping, so the argument is complete as written; only its provenance is
external.

**4.1 Parameters.** Throughout gate elimination we maintain

  - an injective F_2-affine map phi: F_2^m -> F_2^n, written phi(u) = b + L u with
    b in F_2^n and L in F_2^(n x m); its variables u = (u_1,...,u_m) are the
    *parameters*;
  - a set P of disjoint pairs of parameters;
  - a fair semi-circuit C that reads parameters (and auxiliaries) and computes
    C(u) = C_0(phi(u)) for every u in Z(P), where

        Z(P) = { u in F_2^m : p r = 0 for every {p,r} in P }.

Initially phi is the identity (b = 0, L = I, m = n), P is empty, Z(P) = F_2^n, and C
is the input circuit C_0, so s = |C_0|.

A parameter in a pair of P is *paired*, and its partner is the other member; the
remaining parameters are *unpaired*. For a pair {p,r} the constraint pr = 0 admits
exactly the three lines p = 0, r = 0 and p = r + 1 of F_2^2, which is why the
three-point domains below are legitimate.

**4.2 Fair semi-circuits.** A *fair semi-circuit* on the parameters u consists of

  (i) an affine system Av = N u + c of r equations in auxiliary variables
      v = (v_1,...,v_r), where A in F_2^(r x r) is invertible (*fairness*: a unique
      solution v for every u) and every row of (A | N) has at most three nonzero
      entries;
  (ii) an acyclic part: an ordinary circuit, topologically ordered, whose gates read
       constants, parameters, auxiliary variables or earlier gates;
  (iii) an output: a constant or a literal of a parameter, an auxiliary variable or a
       gate. No variable of an equation is designated as the output.

The size is s = (number of equations) + (number of gates in the acyclic part).
Invertibility of A makes every auxiliary variable an affine function of the
parameters, v = A^(-1) N u + A^(-1) c, so the semi-circuit computes a total function
of u.

*Moving an affine gate into the system.* If a gate J of the acyclic part has an affine
operation h(t1,t2) = b1 t1 + b2 t2 + c and both its inputs are constants, parameters or
auxiliary variables, replace J by the equation J + b1 t1 + b2 t2 = c and treat J as a
new auxiliary variable. The new equation has at most three nonconstant terms, the new
auxiliary matrix is block triangular with the old A and a 1 on the diagonal, hence
invertible, and the computed function is unchanged. One gate is replaced by one
equation, so s is preserved. Iterating in topological order moves every affine gate
with no AND-type ancestor into the system and leaves m, q, a and s unchanged while not
increasing the potential.

**4.3 States and the four kinds of change.** The state is (C, phi, P), summarised by

    s = |C|,  m = (number of parameters),  q = |P|,  a = m - q.

Note a = dim of the subspace of Z(P) obtained by fixing one parameter of each pair,
which is the dimension that the disperser premise constrains. The four changes are:

  * **(AR) affine reparameterisation.** For an unpaired parameter x and an affine map
    alpha, substitute x -> w + alpha(u_{-x}) with w fresh; implemented in the circuit
    by exchanging the roles of x and an auxiliary variable Q that depends on x (swap
    the column of x in N with the column of Q in A), so that no equation changes,
    A stays invertible, and only the map phi changes. Then
    Delta m = Delta s = Delta q = Delta a = 0 and Delta W = 0 (both x and w have
    degree >= 1, see 4.6).
  * **(UF) fixing an unpaired parameter.** Impose x = alpha(u_{-x}) (possibly a
    constant). Delta m = -1, Delta q = 0, so Delta a = -1.
  * **(PF) fixing a paired parameter.** For a pair {p,r} substitute p -> 0, or r -> 0,
    or p -> r+1, and remove the substituted parameter and the pair. Delta m = -1,
    Delta q = -1, so Delta a = 0.
  * **(PC) adding a pair.** For distinct unpaired x,y add {x,y} to P. Delta m = 0,
    Delta q = +1, so Delta a = -1; phi is unchanged.

Each change replaces phi by phi o psi with psi injective, maps Z(P) into Z(P), and
preserves the invariant C(u) = C_0(phi(u)) on Z(P). (For (PF) this uses that the three
choices are exactly the lines inside pr = 0.)

**4.4 Degrees and the potential.** The *parents* of a gate are its two input pins; the
gate is a *child* of each internal value it reads (reading the same value twice counts
as one child); the parents of an equation are the internal values occurring in it with
nonzero coefficient, and the equation is a child of each of them. The *degree* deg(t)
of an internal value is its number of distinct children. A parameter is *unused* if
deg = 0 and *isolated* if deg = 1. Degree is syntactic and counts direct occurrences
only.

*Rerouting* means replacing every occurrence of a gate or auxiliary variable by the
constant or literal it has become. Whenever we reroute to a literal of a parameter in a
deletion below, that parameter already occurs in the deleted gate or equation; hence a
parameter that is unused before a deletion is still unused after it.

The weight of an unpaired parameter is

    wt(x) = 0 if deg(x) = 0,   wt(x) = 6 if deg(x) >= 1,

and the potential is

    W = s + 4q + sum over unpaired x of wt(x)   >= 0.

Paired parameters contribute only through 4q. For any quantity X, Delta X denotes
(after) - (before) across one step; a step may group several operations.

**4.5 Extracting a constant subspace.** *Lemma C1.* Suppose a > d and the output of C
is either (i) an affine function of the parameters, or (ii) a literal of an AND-type
gate whose two parents are affine functions of the parameters. Then, given those
affine functions, one can compute in polynomial time an affine subspace
U subset Z(P) of dimension at least a - 1 >= d, together with b in F_2, such that
C(u) = b for every u in U.

*Proof.* Represent an affine function by its coefficient vector, psi(u) = sum_j
alpha_j u_j + alpha_0. Fix one parameter of each pair to 0; this gives an affine
subspace U_0 subset Z(P) of dimension m - q = a, with point 0 and basis the standard
vectors e_j of the a unpaired parameters. An affine function psi is constant on U_0 iff
alpha_j = 0 for all basis vectors e_j of U_0, and then its value is alpha_0. Otherwise
pick a basis vector e_j0 with alpha_j0 = 1; then {u in U_0 : psi(u) = 0} is the affine
subspace with point e_j0 and basis {e_j + alpha_j e_j0 : j != j0}, of dimension a - 1.

Case (i): if the output function psi is constant on U_0, take U = U_0 and b = alpha_0;
otherwise take U = {psi = 0} (dimension a - 1) and b = 0.
Case (ii): the output is (psi_1 + s)(psi_2 + s') + s'' for bits s, s', s''. If both
psi_1 + s and psi_2 + s' are constant on U_0 with values c_1, c_2, take U = U_0 and
b = c_1 c_2 + s''. Otherwise some psi_i + s is non-constant on U_0; take U = {psi_i + s = 0}
(dimension a - 1) and b = s''. In every case C is constant, equal to b, on U.
Since a > d we have dim U >= a - 1 >= d. The construction inspects at most two
coefficient vectors against a basis of U_0, so it is polynomial time. QED

**4.6 Deletions.** Each deletion below satisfies Delta a = 0 and Delta W <= -1.

  * **(D1) delete a gate that depends on at most one input.** This covers constant and
    unary gates, gates with a constant input, gates whose two pins carry the same
    value, an AND-type gate whose two inputs are parameters or auxiliary variables
    with dependent linear parts, and any gate that has become constant on all of Z(P).
    Reroute all its uses to the resulting constant or literal and delete the gate.
    Delta m = Delta q = Delta a = 0, Delta s = -1, and no unused parameter becomes used,
    so Delta W <= -1.
  * **(D2) delete an unused gate** (no children, not the output): Delta s = -1,
    Delta m = Delta q = Delta a = 0, Delta W <= -1.
  * **(D3) delete an equation with at most two nonconstant terms.** Such an equation
    must contain an auxiliary variable (otherwise its row of A would vanish); write it
    as v = U + c and substitute everywhere, then delete v and the equation. Delta s = -1,
    Delta m = Delta q = Delta a = 0, Delta W <= -1.
  * **(D4) unit parameter deletion.** Let x be an unpaired isolated parameter whose
    unique child is an equation E. Since the column of x in N has a single nonzero
    entry (in the row of E), the vector A^(-1) N_x is zero; choose an auxiliary
    variable Q with coordinate 1 in A^(-1) N_x so that Q = x + alpha(u_{-x}), apply (AR)
    to exchange x and Q, and delete the (now auxiliary) x together with E. Delta
    m = Delta q = Delta a = 0, Delta s = -1; the removed parameter had weight 6 and the
    new parameter Q has weight at most 6, and no other unpaired parameter gains weight,
    so Delta W <= -1.

Invariants under all operations: every equation keeps at most three nonconstant terms
and the auxiliary matrix stays invertible (for (AR) via det A' = det A * t_j with
t_j = 1; for the deletions via a determinant expansion along a column with a single
nonzero entry).

**4.7 Preliminary reductions.** At the start of every step, test whether the output has
form (i) or (ii) of Lemma C1; if so, apply Lemma C1 and stop, outputting
S = phi(U) and b. Otherwise repeatedly apply (D1)-(D4) and move every affine gate with
no AND-type ancestor into the affine system, until none applies, and then perform the
first applicable case among the following three. (Throughout, "Cases 1-3 do not apply"
is available as a hypothesis for later cases.)

  * **Case 1: a paired parameter p has degree 0 or at least 2** (r its partner).
    If deg(p) = 0: use (PF) to fix r to 0. The pair is removed and p becomes unpaired
    but unused, so it contributes 0. Delta m = -1, Delta s = 0, Delta q = -1,
    Delta a = 0, and no initially unused parameter becomes used, hence
    Delta W <= -4.
    If deg(p) >= 2: use (PF) to fix p to 0 and delete two of its children by (D1)/(D3).
    Only r (newly unpaired) can gain weight, by at most 6. Delta m = -1,
    Delta s <= -2, Delta q = -1, Delta a = 0, hence Delta W <= -4 - 2 + 6 = 0.
    If Case 1 does not apply, every paired parameter is isolated.
  * **Case 2: two partners have the same child.** Let p, r be partners with common
    child H, which is then their unique child. If H is a gate, its two pins are p and r,
    and H takes the same value at two of the three points (0,0), (0,1), (1,0) of
    pr = 0; the line through those two points is one of p = 0, r = 0, p + r = 1. Use
    (PF) to restrict to that line, so H becomes constant, and delete H by (D1). If H is
    an equation it has exactly three nonconstant terms, hence the form v + p + r = c
    with v auxiliary; impose p = r + 1 by (PF), so H becomes v = c + 1, deleted by
    (D3). In both cases the surviving member of the pair becomes unused. Delta m = -1,
    Delta s <= -1, Delta q = -1, Delta a = 0, and no initially unpaired parameter gains
    weight (we reroute only to constants), hence Delta W <= -4 - 1 = -5.
    If Cases 1-2 do not apply, every paired parameter is isolated and has a child
    different from its partner's child.
  * **Case 3: a paired parameter p occurs in the affine system.** Since Case 1 does not
    apply, p is isolated and its unique child H_p is an equation; likewise for r with
    unique child H_r, and since Case 2 does not apply, H_r != H_p. Use (PF) to fix r to
    0 and delete H_r (by (D1) if it is a gate, by (D3) if it is an equation); p does
    not occur in H_r, so p is now unpaired and isolated with unique child H_p. If H_p
    has at most two nonconstant terms, delete it by (D3); if it has exactly three, apply
    (D4) to p and H_p. Delta m = -1, Delta s <= -2, Delta q = -1, Delta a = 0, the new
    parameter replacing p has weight at most 6, and no initially unpaired parameter
    gains weight, hence Delta W <= -4 - 2 + 6 = 0.
    If Cases 1-3 do not apply, every paired parameter is isolated, occurs only in the
    acyclic part, and has a child different from its partner's.

**4.8 Restructuring the minimal AND.** Assume a > d and that no preliminary reduction
applies. The acyclic part is nonempty, because otherwise the output would be an affine
function of the parameters and Lemma C1 would have terminated the algorithm. Let G be
the first gate of the acyclic part in the topological order; G is AND-type (all affine
gates without AND-type ancestors were moved into the system), and, since (D1) does not
apply, both pins of G carry parameters or auxiliary variables with independent linear
parts. Call G the *minimal AND*.

If a pin of G is an auxiliary variable Q, then Q is a nonconstant affine function of
the parameters, so it depends on some parameter x; that x must occur in an equation
(otherwise the column of x in A^(-1)N would vanish), and since no Case 1-3 applies,
every parameter occurring in the affine system is unpaired, so x is unpaired with
deg(x) >= 1. Apply (AR) to exchange Q and x: the new parameter (formerly Q) has degree
at least 2 because Q occurs in an equation and is a pin of G, so
Delta m = Delta s = Delta q = Delta a = Delta W = 0. Then re-run the preliminary
reductions; if one of Cases 1-3 applies, perform it and end the iteration. Otherwise
repeat for the second pin, which exists and is dependent on a further parameter because
the two pins were linearly independent before the (AR) and an invertible affine change
preserves linear independence. After at most two applications of (AR) both pins of G
are parameters x, y. Write G = (x + a)(y + b) + c with bits a, b, c. If the output is a
literal of G, apply Lemma C1 and stop. Otherwise G is not the output and, since (D2)
does not apply, G has at least one child, and every child of G is a gate of the acyclic
part.

**4.9 Eliminating the minimal AND.** Both cases below have Delta a = -1 and
Delta W <= -10.

  * **Case 4: both pins unpaired.** Complement x and y if necessary (substituting
    x -> w + 1 in phi and in the circuit, which absorbs the negation into constants and
    operation tables and leaves s, m, q, a and W unchanged, the fresh parameter w having
    the same degree as x) so that G = x y + c. Add the pair {x,y} by (PC); now every
    parameter assignment in Z(P) satisfies xy = 0, so G = c on Z(P). Delete G by (D1),
    rerouting it to the constant c, and delete one child of G by (D1) (every child of G
    now has a constant input). x and y lose their weights 6 and the new pair contributes
    +4; no other unpaired parameter gains weight because G is rerouted to a constant and
    the deleted child is rerouted to a constant or a literal of its other input.
    Delta m = 0, Delta s <= -2, Delta q = +1, Delta a = -1, hence
    Delta W <= -2 + 4 - 6 - 6 = -10.
  * **Case 5: at least one pin is paired.** By Cases 1-3 each paired pin of G has G as
    its unique child, and its partner has a unique child, a gate of the acyclic part
    different from G. Let k in {1,2} be the number of paired pins.
    *k = 1:* say x unpaired, y paired with partner r and r's unique child H_r. Use (PF) to
    fix r to 0 (removing the pair), then (UF) to fix x to the constant c so that
    G computes c; delete G by (D1) and reroute. y becomes unused (G was its only child).
    H_r now has the constant input r, so delete it by (D1). x loses its weight 6. No
    initially unpaired parameter gains weight: G is rerouted to a constant and H_r to a
    constant or a literal of its other input. Delta m = -2, Delta s <= -2, Delta q = -1,
    Delta a = -1, hence Delta W <= -2 - 4 - 6 = -12.
    *k = 2:* let r_x, r_y be the partners (distinct, since Case 2 does not apply) with
    unique children H_{r_x}. Use (PF) twice to fix r_x and r_y to 0 (removing both
    pairs), then (UF) to fix x to c so that G computes c; delete G by (D1) and reroute;
    y becomes unused; delete H_{r_x} by (D1). Delta m = -3, Delta s <= -2, Delta q = -2,
    Delta a = -1, and x, y contribute nothing before (paired) or after (x fixed, y
    unused), hence Delta W <= -2 - 8 = -10.

**4.10 The iteration lemma.** *Lemma C2.* Suppose a > d. Then, in deterministic
polynomial time, the algorithm either outputs an affine subspace S of F_2^n of
dimension at least d and a bit b with C_0(x) = b for all x in S, or transforms the state
so that m + s strictly decreases, Delta a is in {0,-1} and Delta W <= 10 Delta a.

*Proof.* If the output has form (i) or (ii), Lemma C1 gives U subset Z(P) of dimension
>= a - 1 >= d with C = b on U; set S = phi(U). Since phi is injective and affine,
dim S >= d, and C_0(phi(u)) = C(u) = b for u in U, so C_0 = b on S.
Otherwise: the preparatory operations (moving affine gates into the system, (AR))
preserve m, s and a and do not increase W; each deletion has Delta a = 0 and
Delta W <= -1; Case 4 and Case 5 have Delta a = -1 and Delta W <= -10; Cases 1-3 have
Delta a = 0 and Delta W <= 0 (computed above: -4, <= 0, -5, <= 0). So the whole
iteration has Delta a in {0,-1} and Delta W <= 10 Delta a. Deletions (D1)-(D4) and
Case 4 decrease s and leave m unchanged, Cases 1-3 and 5 decrease m, and no step
increases s, so m + s strictly decreases. Each iteration performs at most s moves of
affine gates into the system and at most two exchanges, each test being a
Gaussian-elimination computation over F_2 or a constant-size table lookup, hence
polynomial time. QED

**4.11 The class theorem.** *Proof of Theorem C0.* Let C_0 compute the affine disperser
f for dimension d, and suppose |C_0| < 4n - 10d. Start with phi the identity, P empty,
m = a = n, q = 0, s = |C_0|. All parameters are unpaired with weight at most 6, so

    W_0 <= |C_0| + 6n < (4n - 10d) + 6n = 10(n - d).

Apply Lemma C2 repeatedly, stopping only when it outputs a subspace. If the first t
iterations all return new states, then summing Delta W <= 10 Delta a over them gives

    0 <= W_t <= W_0 + 10(a_t - n) < 10(n - d) + 10(a_t - n) = 10(a_t - d),

so a_t > d, and Lemma C2 applies again. Each returning iteration strictly decreases
m + s, whose initial value is n + |C_0|, so after at most n + |C_0| iterations Lemma C2
outputs S of dimension >= d and b in F_2 with C_0 identically b on S. Since f is an
affine disperser for dimension d, f is non-constant on S, so some x in S has
f(x) != b = C_0(x), contradicting that C_0 computes f. Hence |C_0| >= 4n - 10d. QED

**4.12 The constant 4 is not an artefact of a bad potential choice.** Consider any
potential of the form V = s + alpha q + beta I, where I counts used unpaired
parameters, and let D be the drop in V that the analysis can guarantee per unit
decrease of a. The five estimates of section 4.6-4.9 impose

    D <= 2 - alpha + 2 beta      (Case 4:        Delta V = -2 + alpha - 2 beta)
    D <= 2 + alpha + beta        (Case 5, k = 1: Delta V = -2 - alpha - beta)
    D <= 2 + 2 alpha             (Case 5, k = 2: Delta V = -2 - 2 alpha)
    0  >= -alpha - 2 + beta      (Cases 1 and 3: the newly unpaired parameter may gain
                                  the full weight beta while s drops by 2 and q by 1)

while the initial accounting gives W_0 <= |C_0| + beta n and therefore the coefficient
c = D - beta. The last inequality reads beta <= alpha + 2, so

    c = D - beta <= (2 - alpha + 2 beta) - beta = 2 - alpha + beta <= 2 - alpha + (alpha + 2) = 4,

with equality throughout precisely when beta = alpha + 2 and the Case 4 estimate is the
binding one. Our choice alpha = 4, beta = 6 satisfies beta = alpha + 2 and gives
D = 10 (Cases 4 and 5 with k = 2 give exactly -10, Case 5 with k = 1 gives -12, and
Cases 1 and 3 give exactly 0), so c = 10 - 6 = 4, the maximum of this family. Hence no
reweighting of this potential improves the constant, and a larger coefficient would need
a structurally different argument. None is claimed here; this paragraph is only a
limitation statement about the technique.

---

## 5. Part D: assembly, remainder and quantifiers

**5.1 The lower bound.** Fix N >= 2^20. By section 3.3 the payload function g_N on
m = N - k bits is an affine disperser for dimension d(N) <= 16 N^(4/5), and g_N is
computable in polynomial time (it is a restriction of the polynomial-time evaluator of
Lemma A3: the prefix is fixed, and the remaining computation is the same). Applying
Theorem C0 with n = m and d = d(N) gives

    C_B2(g_N) >= 4m - 10 d(N).

By Lemma A4, any circuit for F_N restricts to a circuit for g_N of no larger size, so
C_B2(F_N) >= C_B2(g_N). Using m = N - k > N - 5N^(3/5) and d(N) <= 16N^(4/5):

    C_B2(F_N) >= 4(N - 5N^(3/5)) - 160 N^(4/5)
               =  4N - 20 N^(3/5) - 160 N^(4/5),
    R(N) := 20 N^(3/5) + 160 N^(4/5) <= 180 N^(4/5).

**5.2 All lengths, not a subsequence.** The definition of F_N and of k, r, t is a
function of the integer N alone and is defined for every N; the bound above is proved
for *every* integer N >= 2^20, with no divisibility or sparsity condition, since the
zero padding of the last block is inside the injective map phi and the counting is
against the actual input length N, which includes the k prefix bits. For
0 <= N < 2^20 the bound is not asserted; it is not needed, because the liminf
quantifier only sees sufficiently large N.

**5.3 The coefficient and the explicit threshold.** R(N)/N = 20 N^(-2/5) + 160 N^(-1/5)
<= 180 N^(-1/5) for N >= 1, so R(N)/N -> 0 and

    liminf_{N -> infinity} C_B2(F_N)/N >= 4.

Explicitly, for eps > 0 and every integer N >= N_eps := max(2^20, ceil((180/eps)^5)),
we have 180 N^(-1/5) <= eps, hence C_B2(F_N) >= 4N - eps N = (4 - eps) N. The
requirement is exactly liminf >= c for c = 4 and p = 1, with a proved bound
c N - R(N) where R(N)/N -> 0, and the same single family serves the coefficient and
the exponent, since both coordinates describe this one theorem.

**5.4 Deterministic arithmetic spot-checks (bounded computation, not a proof).** A
script (`param_check.py`, exact integer arithmetic, no floating point in any test)
recomputed k, m, r, t for 6,245 lengths: every N in [2^20, 2^20 + 4000), the powers
N = 2^e for 20 <= e <= 53, N = 10^e for 6 <= e <= 16, the 200 values 10^16 - j for
0 <= j < 200, and 2,000 pseudo-random N in [2^20, 10^16] (fixed seed). For each it
checked, literally:

    k prime; k^5 > 32N^3; no prime lies in (2N^(3/5), k)   -- k is the least such prime
    k < 5N^(3/5); m > N/2; 3 <= r <= N^(2/5); 2 <= t <= 2N^(1/5); t < r
    k >= 2 r(r+1)/t
    ceil(nr/t + nt + r(r+1)) <= 16 N^(4/5)
    20 N^(3/5) + 160 N^(4/5) <= 180 N^(4/5);  4N - R(N) > 0

All 6,245 lengths pass every test with no failure. The largest observed ratios were
max k/N^(3/5) = 2.0042 against the proved bound 5, and max d(N)/N^(4/5) = 3.0861 against
the proved bound 16, so the margins are comfortable rather than tight. Example rows:
N = 2^20 gives k = 8209, m = 1,040,367, r = 127, t = 12, d(N) = 201,643;
N = 10^12 gives k = 31,697,867, m = 999,968,302,133, r = 31,547, t = 178.
This finite check supports, but does not replace, Lemma A1 and the estimates of
section 3.2, which are proved for every N >= 2^20; it establishes no circuit statement.

---

## 6. Part E: evidence scope, locators and limitations

**6.1 What is proved inside this file.** Lemma A1 (parameter bounds, all N >= 2^20),
Lemma A2 (irreducibility criterion), Lemma A3 (a single polynomial-time evaluator),
section 2.4 (both directions of NP membership with witness bound 0), section 2.2 and
section 2.5 (the family and the hard restriction), Lemma A4 (restriction of circuits),
section 3.2 (parameter verification for the cited disperser theorem), section 3.3
(transfer of the disperser property to the payload), all of Part C (Lemma C1, the
deletions, Cases 1-5, the iteration lemma and Theorem C0 with its full accounting), and
section 5 (remainder, all-length quantifiers, explicit threshold).

**6.2 Which premises come from published sources.** Two published theorems are used as
premises, with exact locators:

 1. Ben-Sasson, Kopparty, *Affine Dispersers from Subspace Polynomials*, author
    preprint dated 2010-03-12, Theorem 7.3 (p. 34) and Theorem 2.2 (p. 6),
    https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf (SHA-256
    922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf). Used for
    Theorem B1 only, i.e. for the affine-disperser property of the payload. Its
    internal proof (Sections 3, 6, 7 of that paper) is not reproduced here; the
    hypotheses are verified here in section 3.2 and the conclusion is applied to the
    literal function of section 2.2 via the injective map phi.
 2. Mishra, Santhanam, *An Explicit 4n - O(log n) Size Lower Bound against General
    Boolean Circuits*, ECCC TR26-260, revision 1 (2026-11-03),
    https://eccc.weizmann.ac.il/report/2026/260/revision/1/download/ (SHA-256
    a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921). Used for the
    architecture of Part C (Theorem 4.1 there). Part C is nevertheless written out
    here in full: every state, operation, deletion, case, inequality and the final
    summation are stated and derived, so the reader does not need that paper to check
    the argument.
 3. The model conventions of Li and Yang, ECCC TR21-023 (2021-02-20),
    https://eccc.weizmann.ac.il/report/2021/023/download/ (SHA-256
    4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062), Sections 1.1,
    2.1, 2.4 (footnotes 3-4) and 3.3, are the B2 conventions fixed by `rules.md`; they
    are used only for the model, in Part C's definitions and in the size accounting.

**6.3 Honest limitations.**
 * The affine-disperser premise is external (item 1 above). This file verifies its
   hypotheses for our parameters and transfers it to the payload, but does not
   re-derive the subspace-polynomial argument. An adversary could therefore press on
   that citation rather than on Parts A, C and D.
 * The class theorem's provenance is external (item 2), although its proof is supplied
   here in full. Section 4.12 records that no reweighting of the potential improves the
   constant 4, so this submission deliberately claims exactly 4 and no more.
 * The bound is proved for every integer N >= 2^20; shorter lengths are defined (as
   zero) but are not claimed to satisfy anything.
 * Nothing here says that affine dispersers need superlinear circuits, that the constant
   4 is optimal for the family, or that any exponent above 1 holds. No class separation
   follows from a fixed coefficient at p = 1.
 * The reported spot-checks in section 5.4 are a bounded computation over a stated
   range of lengths; they are arithmetic verification, not a proof of the analytic
   estimates, and they establish no circuit statement.