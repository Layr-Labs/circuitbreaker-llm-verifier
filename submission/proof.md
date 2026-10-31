# A fixed full-B2 family with error O(N^(1/100))

This file and supporting_proofs.md together are the complete mathematical
evidence. Both are required. The standalone reader manuscript contains the
same theorem with additional exposition and an earlier parameter example;
it is not needed to supply a missing premise of these two files.

## Joint theorem and scope

One fixed total P language, defined by the permanent choice M=100 in
Chapter I, has N-bit characteristic function F_N satisfying

    C_B2(F_N) >= 4N-8N^(1/100)-6N^(3/50000)-23

for EVERY sufficiently large integer actual input length N. In particular
C_B2(F_N)>=4N-9N^(1/100) eventually. The model has all sixteen binary truth
tables, arbitrary sharing and fan-out, one output, and exact computation.
Every header and padding bit is counted. The exact NP verifier in Chapter I
accepts iff its witness is empty and the evaluator returns one, with witness
bound zero. The family and evaluator do not change within the asymptotic
quantifiers.

This implies liminf C_B2(F_N)/N>=4, with leading expression 4N and finite
projection (4,1). It substantially improves round six's N^(4/5) error but
does not assert coefficient above four. No practical numerical threshold
is claimed: the fixed polynomial degree and eventual threshold are enormous.
The same proof permits each separately fixed M>=2; it does not let M grow
with N in the one language or infer an N^o(1) error for one family.

Chapter I defines and analyzes the total evaluator; II and III prove all
gate accounting; IV proves field existence and validation. In the second
file, V constructs the seeded extractors, VI proves cubic bias-to-flat,
VII proves the cubic trace theorem, and VIII gives a direct reciprocal-
parameter proof of the full recursive composition. Section numbers inside
chapters are local. All logarithms are base two unless marked ln. Affine
sources are uniform on nonempty F_2-affine flats, with entropy equal to
dimension. No review, finite check, hosted status, or citation is a premise.

# Chapter I. Concrete total evaluator and actual-length wrapper

This chapter fixes one language, gives finite evaluation rules at every input length, and proves the actual-length accounting. The atomic component and recursive composition proofs accompanying this chapter are part of the theorem. The earlier reviewers' conclusions are not premises.

## 1. Fix the language constants once

Fix the integer `M=100` once and for all. In the definitions below it is a constant in the program, not an input parameter. Write

    kappa=1/M, beta=kappa/1000, d=1000*M, C=d+5,
    gamma=kappa/(1000*C), Aexp=222000*M,
    B=ceil(2*Aexp/gamma)+2,
    tau=kappa^2/(4000000000*(B+1)),
    q_tree=200000000*M*(B+1), J=20*M-1, H=J*q_tree.

These are explicit integers or rational constants. For a payload length n define

    K=ceil(n^kappa), K0=ceil(n^(kappa/2)),
    L0=floor(n^(kappa/6)), Rerr=ceil(n^beta),
    w=ceil(n^(6*beta)), h=2*w+1,
    Pconf=floor(n^(beta/4)), G=floor(n^gamma),
    t=ceil(n^tau), F=t^B, ell=H*F,
    p_u^(a)=ell^(6*(J-u)+3*a), 0<=u<J, a in {1,2}.

Here h is the field-header length, ell is the number of padded challenge rows, and Pconf is the full confidence length. These quantities are distinct. The symbol q_tree is the fixed number of ordinary tree levels between successive active levels; it is not the number of bits fixed in the final circuit restriction.

For exact rational powers, `floor(n^(a/b))` means the largest integer z satisfying `z^b<=n^a`. The ceiling is defined similarly. Binary search with integer power comparisons implements these operations; a real-number oracle is not used. All exponents and denominators in this chapter are fixed constants.

## 2. One finite, total component guard

Only well-definedness is needed at small lengths. The following explicit guard is a convenient, stronger choice. It checks every finite component precondition that the evaluator uses. It does not need to test a source distribution, an affine flat, or the asymptotic probability calculations in the proof.

First, if `floor(log2 n)<H^2`, set the guard to false without constructing the large tree constants as expanded unary objects. H itself is an ordinary fixed integer represented in binary, and H^2 can be computed once by integer multiplication. This harmless cutoff makes small inputs return promptly. It is not asserted to imply the later tests; perform all later tests if the cutoff passes.

For the short-seed component A, use the atomic theorem with entropy exponent `kappa/2`, error `1/8`, and output length L0. Its deterministic parameter calculation is:

    Ln=ceil(log2 n), ML=ceil(log2 L0),
    Bcode_A=9+Ln+2*ML,
    s_A=least power of two at least Bcode_A,
    a_A=2*s_A,
    d_A=ceil(kappa*Ln/6),
    D_A=ceil(128*a_A^2/d_A), Rseed=2^D_A.

As usual the logarithmic ceiling of one is zero. Require `n>=2^(18*M)` and L0>=1 before using these formulas. The stated threshold ensures d_A>=1. The atomic greedy set-system proof then guarantees construction of all L0 sets of size a_A. To remove any possible implementation ambiguity, the evaluator enumerates a_A-subsets in lexicographic order and chooses the first eligible one at every step. The full seed set is all D_A-bit strings in lexicographic order.

Define a routine `Bpars(Kb,mb)` for the longer-seed component. It uses the fixed global input length n, error exponent Rerr, and overlap d. Compute

    Bcode=3+ceil(log2 n)+2*ceil(log2 mb)+2*Rerr,
    s=least power of two at least Bcode,
    a0=2*s,
    q0=least power of two at least a0,
    D=q0^2.

The routine accepts only if

    1<=mb<=n, 1<=Kb<=n,
    d<a0, q0^(d+1)>=mb,
    Kb>=mb*2^d+1+3*ceil(log2 mb)+3*Rerr.       (Finite-Graph)

These are conservative integer versions of the exact reconstruction inequalities proved in the atomic component chapter. Acceptance guarantees a strong input-linear extractor of output mb, seed length D, and error at most `2^(-Rerr)`. Its graph sets are the first mb degree-at-most-d polynomials in the prescribed coefficient-vector order. Every field in this component is constructed by the atomic tower or polynomial-size field search; no extra header is used.

For the merger, **pad every chunk to exactly G rows with zero rows**. This is permitted by the affine merger proof and makes all merger output lengths the same. Set

    bmerge=ceil(log2 G), T_j=ceil(G/2^j),
    L_0=L0,
    L_(j+1)=floor(L_j/(16*2^d)), 0<=j<bmerge.

The guard requires all of the following:

1. `G>=2`, `t>=B`, `t^H<=n`, `ell>=2`, `ell*Pconf<=L0`, and `p_0^(2)<=Pconf`.
2. For every `0<=j<bmerge`, require `T_j*L_j<=n`, `T_j*w<=floor(L_j/2)`, and `L_(j+1)>=1`; require `Bpars(floor(L_j/2),L_(j+1))` to accept. Denote its returned seed length by D_j.
3. Require `Bpars(K0,L0)` to accept; denote its seed length by D_reg. Require `L_bmerge>=D_reg` and `G^(B-1)>=Rseed`.
4. Set `M_Q=floor(sqrt(h)/64)`. Require `h>=2^16`, `h>=(64*Rerr)^2`, and `M_Q>=D_j` for every merger step j. The seed length D_reg is supplied by the final merger output; it need not also be supplied directly by Q.

If any test fails, the guard is false. Otherwise it is true. Every test is an exact finite integer comparison. The square-root floor in M_Q can be obtained from the integer square root of h. This list is intentionally redundant: for example, the finite reconstruction inequality is stronger than is needed simply to evaluate the component, but it gives a transparent finite correctness guarantee whenever evaluation proceeds.

The bounds on probabilities and the source-entropy reserves in the composition proof are not runtime tests. They are properties proved to hold for all sufficiently large n. Defining arbitrary outputs at smaller n does not affect the eventual lower bound.

## 3. Why this guard eventually passes at every integer length

The reciprocal-parameter proof supplies the following strict margins; they are recorded here to verify all finite guards, without selecting a sparse length subsequence.

The short-seed theorem is valid once `n>=2^(18M)`. Its seed count satisfies `Rseed<=n^(221200M)<=n^Aexp`. Also `G>=n^(gamma/2)` eventually, and `(B-1)*gamma/2>Aexp`, so `G^(B-1)>=Rseed` eventually.

The tree satisfies `tau*H=1-kappa/20<1`, so `t^H<=2^H*n^(1-kappa/20)<=n` eventually. We have `t>=B` and G>=2 eventually. Moreover

    ell<=n^e eventually, where e=kappa^2/2000000000,
    e+beta/4<kappa/6,
    (120M)*e<beta/4.

These show `ell*Pconf<=L0` and `p_0^(2)=ell^(120M)<=Pconf` eventually.

The halving merger has at most `gamma*log2 n+1` steps. There is no circular use of a minimum row length: with the integer `Ddiv=16*2^d=2^(C-1)`, repeated integer division gives the exact identity `L_j=floor(L0/Ddiv^j)`. Consequently, uniformly for all `j<=ceil(log2 G)`, we have `L_j>=2^(-C)*n^(kappa/6-(C-1)*gamma)-1`. For sufficiently large n this is at least the following weaker bound, both for every current row length and for the final row length:

    L_j>=2^(-C-1)*n^(kappa/6-C*gamma)
       =2^(-C-1)*n^(497*kappa/3000)
       >=n^(kappa/8)

eventually. Since

    gamma+6*beta<kappa/8,
    gamma+kappa/6<1,

we obtain every slice-cost and input-length guard. For the merger component, `Kb=floor(L_j/2)>=n^(kappa/8)/3`, and `mb=floor(L_j/(16*2^d))`. The advice term `mb*2^d` is at most L_j/16, whereas Kb is asymptotically L_j/2. The remaining reconstruction term is `O(log n+n^beta)=o(L_j)`, because beta<kappa/8. For regeneration, `Kb=K0` and `mb=L0`; the exponent inequality kappa/6<kappa/2 makes its advice term negligible compared with K0, and beta<kappa/2 handles the rest. The design conditions d<a0 and q0^(d+1)>=mb hold eventually by the atomic graph-set argument, since beta(d+1)>1 and mb<=n.

The long-seed components all have seed length at most `2304*n^(2*beta)` eventually. The cubic Q produces `M_Q>=sqrt(h)/128=Omega(n^(3*beta))` bits, so every D_j fits. Also `sqrt(h)/64>>Rerr` since 3beta>beta, giving `h>=(64*Rerr)^2`. The merger's final output exponent kappa/8 exceeds 2beta, so it contains D_reg bits. All remaining fixed lower thresholds, including `h>=2^16` and the preliminary `floor(log2 n)>=H^2`, hold eventually.

These comparisons prove that there exists an integer n_guard such that the guard is true for **every** n>=n_guard. Each comparison is a fixed positive power gap with explicit constants. The elementary facts `n^a/n^b ->0` for a<b and `log n=o(n^c)` for c>0 follow, for example, by substituting n=2^v and comparing a polynomial in v with the exponential series. Thus the threshold is effective and no number-theoretic length-selection premise is hidden here.

## 4. The total Boolean evaluator

For an actual input string z of length N<4, output zero. For N>=4 choose the largest integer n>=1 satisfying

    n+h(n)<=N,  h(n)=2*ceil(n^(3/(500M)))+1.

Such an n exists because 1+h(1)=4. Interpret the first h(n) bits as coefficients `a_0,...,a_(h-1)` of

    Q(T)=T^h+sum_(i=0)^(h-1) a_i*T^i.

The next n bits are the payload x; the remaining rho bits are padding y. All these positions and lengths depend only on N.

Let base(z) be zero if the finite guard above is false or if Q fails the polynomial-time irreducibility test proved in the field-support chapter. Otherwise let base(z) be the output of the fully specified recursive challenge-response evaluator, with the fixed reciprocal parameters above and the valid field representation `F_2[T]/Q(T)`. It uses one such field for every cubic call. Its merger chunks are padded to exactly G as specified above. Define

    F(z)=base(z) xor parity(y).

In particular, invalid-header and failed-guard inputs still take part in the final padding parity. The empty parity is zero. This defines a Boolean value on every finite binary input string.

Fix all ordering conventions as follows. Bits in field coefficient vectors are listed from lowest to highest basis index. All finite enumerations use increasing binary lexicographic order; a set's elements are listed in increasing universe index. Pad short inputs by appending zeros, list matrix bits in row-major order, partition lists into consecutive chunks, and enumerate increasing child tuples lexicographically. An interval of length v is split into t consecutive children by giving the first v mod t children length floor(v/t)+1 and the rest length floor(v/t). These conventions specify all ties and all coordinate embeddings.

For clarity, the recursion can be implemented as a deterministic tree computation. Split each coordinate interval as just specified for H levels. Compute child challenges before parent challenges. At each parent choose the first child whose full-Pconf challenge matches none of the enumerated response matrices, using child 1 if none does. Prefix the chosen challenge, then append the BE outputs on all increasing B-tuples of children and pad this appended part to F rows. Starting at the root, follow these stored child choices. At active depths u*q_tree, compare the challenge at the next active node with the responses on the active node, at confidence lengths p_u^(1) and p_u^(2). Define FAIL to be the Boolean value 0 and PASS to be 1. Return the second comparison bit at the first active node whose first comparison says FAIL; if no such node occurs, return zero. This is exactly the recursive evaluator proved in the composition chapter with the present parameters; it makes no source-dependent choices.

## 5. Total polynomial time and the NP verifier

For fixed M all constants, polynomial degrees, tree depth H, and block count B are fixed. The integer-power parameters, header search for n, and finite guards are polynomial time in N. Binary search finds the maximal n because n+h(n) is strictly increasing; even a linear scan would be polynomial. All intermediate integers have polynomially many bits in N, with constants depending on M. The preliminary cutoff avoids expanding enormous fixed objects on small inputs, but total polynomial time would hold in the formal asymptotic sense even without this convenience.

When the guard succeeds, the tree has at most O(H*t^H)=O(Hn) nodes. Each node appends at most t^B rows, which is bounded by a fixed power of n. The response list has `Rseed<=n^(221200M)` entries. The short-seed set system is obtained by enumerating at most `2^D_A=Rseed` subsets. All code coordinates are computed individually over fields whose **bit dimensions** are polynomial in n; the exponentially long codeword is never materialized. All graph-design fields have polynomial cardinality and can be found by the proved small-field search. The only larger odd-degree field is supplied in the h-bit header and validated in polynomial time. Every gate of the evaluator's arithmetic is a uniform bit operation; no advice or uncounted seed occurs.

Consequently there is one deterministic Turing machine deciding F in time N^c for some fixed constant c. This is a P language. An NP verifier is the following exact relation:

    V(z,u)=1 if and only if u is the empty string and F(z)=1.

It rejects every nonempty witness. Its witness bound can be the constant zero. In both directions, an input has an accepting witness exactly when its evaluator value is one.

## 6. Payload theorem used by the length wrapper

The atomic extractor proofs and the recursive composition proof, with the fixed reciprocal parameters above, establish:

> For every sufficiently large integer n and **every** irreducible degree-h(n) modulus Q, the payload function D_(n,Q) computed by the recursive evaluator is nonconstant on every affine flat of dimension at least `K=ceil(n^(1/M))`.

The finite guard passes eventually by Section 3. The additional proof-only inequalities for the recursive selection also hold at all sufficiently large n: selected block entropy has exponent at least 3kappa/4, exceeding the response threshold exponent kappa/2, while all imposed linear equations have total length exponent less than beta; the largest rare-event exponent is `ell*p_0^(2)=o(Pconf)`, and the smallest confidence ell^9 exceeds any fixed multiple of log n. These facts are derived with exact constants in the composition chapter. This wrapper does not replace that proof by a guard or a finite test.

For large n, K<=n-1. Thus D_(n,Q) is nonconstant and is not a literal or a complemented literal: such a literal is constant on an (n-1)-flat, contradicting the displayed property. Applying the proved ordinary-affine-disperser full-B2 gate theorem with d=K-1 gives

    C_B2(D_(n,Q))>=4n-8K.                         (Core)

## 7. Every overhead bit earns a distinct restriction deletion

Put `q_fix=N-n=h+rho`. There is at least one valid monic degree-h binary irreducible polynomial, by the field-existence proof. Fix one such Q. Its degree is h>=3. Since an irreducible polynomial of degree greater than one has no root at 1, Q(1)=1. Flipping any one of its h coefficient bits changes this value to zero, so the resulting monic degree-h polynomial has the factor T+1 and is invalid.

For sufficiently large n, D_(n,Q) takes value one on at least one payload. Choose that payload and padding zero. The valid header makes F=1, whereas changing any single header coefficient makes the base bit zero and hence F=0. Therefore **every header bit is an essential input of F**. Every padding bit is essential as well, because it participates in the final parity for every setting of the other inputs.

Start from any full-B2 circuit for F. Normalize constant, unary and repeated-parent gates and prune dead gates, retaining original labels on survivors. Fix all h header bits to Q and all rho padding bits to zero. The restricted function is exactly D_(n,Q), which is nonconstant and not a literal. The proved frozen-port restriction lemma removes at least one distinct retained gate for each essential fixed root, hence at least q_fix gates in total. Its proof uses one actual restriction with the original pin incidences frozen; it does not add losses from different hypothetical restrictions. The deleted labels, any earlier normalization labels, and the retained core labels are disjoint.

Every survivor remains a legal B2 table after literal and constant propagation; sign changes can be absorbed into a retained table because the final output is not a literal. Thus (Core) applies to the retained circuit. Counting deleted and retained labels once each yields

    C_B2(F_N)>=q_fix+4n-8K
               =4N-3*q_fix-8K.                   (Counted)

Chapter III proves the restriction lemma and Chapter II proves the full-B2 gate theorem; neither is an appeal to an audit judgment.

## 8. All-length arithmetic and the final bound

Let `a=3/(500M)`. Because 0<a<1, the function x^a is increasing and has increments `(n+1)^a-n^a<=1` for n>=1. For example this follows from the derivative a*x^(a-1)<=1; the same inequality also follows directly from concavity. Hence

    ceil((n+1)^a)-ceil(n^a) is either 0 or 1.

It follows that `n+h(n)` increases by either 1 or 3 when n increases by one. Maximality of the selected payload length therefore gives

    0<=rho=N-n-h(n)<=2,
    q_fix=h+rho<=2*n^a+5.

This proves actual-length coverage without assuming N is prime, a power, or a member of any special subsequence. Moreover n(N) tends to infinity. In fact every N>=n0+h(n0) has n(N)>=n0. Once the payload, guard and proof thresholds all hold beyond a fixed n0, the theorem holds at every integer N beyond this explicit monotone transfer.

Finally K=ceil(n^kappa)<=n^kappa+1, so (Counted) implies

    C_B2(F_N)
       >=4N-3*(2*n^a+5)-8*(n^kappa+1)
       =4N-8*n^kappa-6*n^a-23
       >=4N-8*N^kappa-6*N^a-23.

For the one fixed choice M=100 this is

    C_B2(F_N)>=4N-8*N^(1/100)-6*N^(3/50000)-23

for every sufficiently large integer N. Since a<kappa, it also gives

    C_B2(F_N)>=4N-9*N^(1/100)

eventually. All input bits, including the h-bit modulus and at most two padding bits, are counted in N. The circuit proof covers all sixteen binary Boolean tables, one output, and arbitrary sharing and fan-out.



# Chapter II. Full-B2 gate bound from ordinary affine dispersal


Theorem. If f:F2^n -> F2 is nonconstant on every affine flat of dimension
greater than an integer d>=0, then C_B2(f)>=4n-8(d+1).

The potential is V=s+4I+2q. Its surplus is different from the surplus of
the earlier W=s+5I+3q histories; no earlier token or AND-mass inequality is
transferred to it.

## 1. State and domain

Keep M independent binary parameters. There are q disjoint parameter pairs; each pair {p,h} ranges over F_2^2 minus one excluded corner (e_p,e_h). All other parameters range freely. Let D be this product domain.

Keep an affine injection Phi:F_2^M -> F_2^n, defined and injective on the whole parameter cube. The represented computation equals f(Phi(u)) for u in D. Equality outside D is not required.

The computation consists of:

* a square invertible affine system A v = N u + c, where each row has at most three distinct nonconstant terms;
* acyclic B2 gates, with no gate output entering the affine system;
* a selected output that is a constant or a literal of a parameter, auxiliary, or gate value.

Let s count all retained rows and gates, one charge each. A parameter's children are the distinct rows and gates containing it; the output reference is not a child. Let I count used free (unpaired) parameters. Define

    a = M-q,     V = s+4I+2q.

Initially M=a=n, q=0, Phi is the identity, s=|C|, and V<=|C|+4n.

Fixing one endpoint of each pair to the complement of its excluded coordinate leaves its partner free. Such a safe fixing picks an affine a-flat inside D. Its image under Phi is an affine a-flat of the same dimension. Therefore a constant represented output is impossible whenever a>d.

While a>=d+2, an affine represented output is also impossible: on a safe a-flat it has a constant subflat of dimension at least a-1>d. The same is true of an AND-type output whose parents are affine values. On the safe flat, fix a nonconstant affine factor to its blocking value; if there is no such factor, the output is affine or constant. Thus a minimal AND cannot be the selected output during the active loop.

## 2. Legal coordinate operations

The following preserve the domain form and affine injectivity, with the corresponding affine substitution also made in Phi.

1. Replace a free coordinate x by w+ell(other coordinates), where ell is affine. This is a bijection on each free two-point fiber, even if ell uses paired coordinates. It adds no gate.
2. Fix a free coordinate to a constant: M and a each fall by one; q is unchanged.
3. Safely fix one paired endpoint: M and q each fall by one, the partner becomes free, and a is unchanged.
4. Restrict a pair to any affine line contained in its three-point domain: parametrize the line by one free coordinate; M and q fall by one, with a unchanged.
5. Restrict two free coordinates to three corners of their square by declaring them a new pair: q rises by one and a falls by one; M and Phi do not change.

No nonlinear change of parameter is used.

## 3. Affine movement and elementary deletions

After absorbing negations into literal edges/tables, a nondegenerate binary gate is XOR-type or AND-type. All other tables are constant or depend on at most one parent.

Move an affine gate having no AND ancestor into the affine system, in topological order. Add its output as an auxiliary and its defining equation as a row. Its old gate charge becomes its row charge; nothing is added to s. The enlarged auxiliary matrix is block triangular with new diagonal entry one. Its row has at most three terms. The move preserves M,q,a,I,V.

Use the following elementary deletions.

**D1.** Delete a constant or one-parent-literal gate and reroute its uses. This includes fixed pins and repeated pins. If an AND has affine parents with dependent linear parts, those parents are constant or equal up to a constant, so its output is likewise constant or a literal and can be deleted.

**D2.** Delete a gate with no children unless it is the selected output.

**D3.** A row with at most two nonconstant terms contains an auxiliary, since A is invertible. Solve it as v=U+c, where U is absent or one other value, substitute, and remove its row and auxiliary column. Substitution does not increase row width. Adding the solved row to each row containing v and expanding on its column proves that the remaining auxiliary matrix is invertible.

For D1 and D3, any parameter receiving a rerouted use was already used in the deleted object; thus I cannot increase. D2 also cannot increase I. Each deletion removes a distinct retained charge and decreases V by at least one.

**D4.** Suppose a free parameter x has sole child row E. The column N_x is nonzero, so some auxiliary Q has x-coefficient one in A^{-1}N_x. Write Q=x+ell(other parameters). Exchange x and Q in role: Q becomes the new free parameter, and x becomes auxiliary. Substitute x=Q+ell in Phi. The new auxiliary matrix replaces A's Q-column by N_x and is invertible, since

    det(A with Q-column replaced by N_x)
      = det(A) * (A^{-1}N_x)_Q = 1.

No equation's terms change, so row width is unchanged. The old x now occurs as an auxiliary only in E, and it is not the selected output by the active terminal test. Remove x and E; expanding along its single-occurrence column again proves invertibility. The old used free x is replaced by at most one used free Q, so I does not increase. One charge is removed, and V drops by at least one.

The same role exchange can prepare an auxiliary AND parent Q without deleting a row. In that case the old free x is used in the system, and the new free Q is used by the AND; I is unchanged. The other statistics are unchanged. This is an invertible free-coordinate substitution, not an encoder.

## 4. Pair cleanup

At the start of each iteration, first move all available affine gates into the system. Then apply elementary deletions whenever possible, followed by the first applicable pair rule below. After any completed action, restart the iteration. In particular, finish the required pair cleanup after affine movement and before preparing an AND source. Propagate all forced constants and literals using D1/D3. When counting, count distinct old objects: an object already deleted supplies its one charge, never a second charge. Fixing a pin or shortening a row eventually removes the affected object, even if other propagation removes it first.

**P0: an unused paired endpoint.** If p is unused, safely fix its partner h and remove their pair. The surviving free p is unused. Then

    Delta a=0, Delta q=-1, Delta I<=0,
    V_before-V_after >= 2.

Any further deletions improve the bound.

**P2: an endpoint with at least two children.** Safely fix that endpoint and remove its pair. Its two distinct affected children are removed. Only its partner can become a new used free parameter. Therefore

    Delta a=0, Delta s<=-2, Delta I<=1, Delta q=-1,
    Delta V <= -2+4-2=0.

After P0/P2 fail, every paired endpoint has exactly one child.

**PC: partners with the same child B.** If B is a row, its third term is an auxiliary and its output is an affine function of the two endpoints. If B is a gate, its output is a Boolean function of them. Among three allowed pair points, two give the same output value. Those two points form an affine line inside the allowed domain. Restrict to it and remove the constant B. Both old endpoints occurred only in B, so the surviving line coordinate is unused. Hence

    Delta a=0, Delta s<=-1, Delta I<=0, Delta q=-1,
    Delta V<=-3.

**PE: a paired endpoint p has sole child row E.** Its partner h has a different sole child B, since PC is absent. Safely fix h and remove B and the pair. The free p still occurs only in E. If E has already disappeared in propagation, it supplies the second distinct deletion. Otherwise apply D4 to remove E. At most one newly free parameter is used, including a replacement parameter introduced by D4. Thus

    Delta a=0, Delta s<=-2, Delta I<=1, Delta q=-1,
    Delta V<=0.

When these rules all fail, every paired endpoint has one gate child, partners have distinct children, and paired parameters have zero columns in N.

## 5. Prepare a minimal AND

After affine movement and the required elementary and pair cleanup, if no acyclic gate remains, the selected output is affine, contrary to the active terminal test. Hence some AND G has no AND ancestor. Its parents are affine values. Elementary cleanup ensures their linear parts are nonzero and independent.

Paired columns of N are zero, since PE is absent, so all auxiliaries depend only on free parameters. For an auxiliary parent Q of G, exchange it with a free parameter on which it depends. If the other parent is already a free parameter, independence supplies a nonzero coefficient outside that parameter; if it is paired, Q has no dependence on it. If both parents are auxiliary, after the first exchange independence supplies a free pivot distinct from the first new parent. Thus at most two exchanges make both parents parameters.

Recheck elementary and pair rules after each exchange and complete any resulting action. If none arises, continue preparing the fixed G. This needs no endless source reselection. G has a child: it is not the selected output by the terminal test and it is not dead by D2.

## 6. Source restrictions and their charges

The parents of G are distinct parameters. Every AND-type binary table has one exceptional corner and one majority value at its other three points.

### GP: both source inputs free

Let G read x,y and have k>=1 distinct gate children. Restrict (x,y) to its three majority points and record them as a new pair. G is constant. Remove G and all its k children by fixed-pin propagation. These are 1+k distinct old gates, even if a child also reads x or y or one child is an ancestor of another. Each of x,y was used and now leaves I; old free parameters cannot acquire a first use through rerouting.

    Delta M=0, Delta q=1, Delta a=-1,
    Delta I<=-2, Delta s<=-(1+k).

Consequently

    Delta V <= -(1+k)-8+2 = -7-k <= -8.

### AP1: one source input paired

Let G read paired p and free x; h is p's partner. Safely fix h and remove its pair. Then fix x to block G. The now-free p has no remaining use, since G was its sole child. Let B be h's old sole child. PC guarantees B!=G. Remove G,B.

If G has a child J!=B, delete J as the third distinct object. Otherwise B is G's sole child and reads both G,h. Both inputs become constant, so B becomes constant and cannot be the selected output: the resulting safe dimension is a-1>d. A child of B supplies the third distinct old gate. Acyclicity makes it different from G,B.

This count is independent of propagation order. If B first becomes a literal of G, its later rerouted uses are still fixed when G becomes constant.

    Delta M=-2, Delta q=-1, Delta a=-1,
    Delta I<=-1, Delta s<=-3.

Hence

    Delta V<=-3-4-2=-9<=-8.

### AP2: both source inputs paired; four deletions

Let G read p,r, whose respective partners are h,k. The pairs differ; otherwise PC would apply to G. Safely fix h,k and remove both pairs. Fix the newly free p to block G. Then r becomes unused, since its sole child was G. Let B_h,B_k be the old sole children of h,k. Neither equals G, since PC is absent.

**Case 1: B_h!=B_k.** Removing G,B_h,B_k gives three distinct charges. If G has a child outside {B_h,B_k}, it is a fourth. Otherwise choose a child inside that set, say B_h. Its parents are G,h, which both become constant, so B_h becomes constant. If it has a child outside {G,B_h,B_k}, that child is a fourth deletion. It cannot feed G, whose parents are p,r, or itself, by acyclicity. If its only remaining child is B_k, then B_k reads B_h and k, both constant, and becomes constant too. Nonconstant output and acyclicity force a child outside the three-object set, supplying the fourth deletion. A selected output inside this constant chain is impossible after restriction, since a-1>d.

**Case 2: B_h=B_k=B.** The gate B reads h,k, so after fixing them B becomes constant and does not read G. G reads p,r, so neither G nor B reads the other. Both constant nodes have a child, since the restricted output is nonconstant. If they have different children, those children together with G,B give four distinct deletions. Otherwise their common child reads both constants, becomes constant, and has a child distinct from these three by acyclicity and nonconstant output. This is again a fourth deletion.

Every counted object is an old charged gate that is actually removed by the selected restrictions and propagation. If a previous substitution removed it already, its existing deletion counts once. There is no addition of two alternative restrictions or of hypothetical deletion sets.

The state changes are

    Delta M=-3, Delta q=-2, Delta a=-1,
    Delta I<=0, Delta s<=-4.

Therefore

    Delta V<=-4-4=-8.

All source restrictions end at safe dimension a-1>=d+1>d. Thus each use of nonconstant restricted output is valid. Further constant, literal, or dead propagation removes charges and never introduces a first use of an old free parameter, so it only improves all bounds.

## 7. Coverage, termination, and summation

Run while a>=d+2. Every active iteration performs finite affine movement and elementary/pair cleanup, or prepares a fixed minimal AND using at most two legal exchanges and applies GP/AP1/AP2. These exhaust its possible parameter-parent types.

Every completed action preserves the representation, never increases M or s, strictly decreases the nonnegative integer M+s, and obeys

    V_before-V_after >= 8(a_before-a_after).

Preparatory moves have zero changes in these quantities and cannot create an infinite preparation loop by the fixed-source protocol above. The process therefore terminates at a_T<=d+1. Since V_T>=0,

    |C|+4n >= V_0
             >= V_T+8(n-a_T)
             >= 8(n-d-1),

which yields

    |C| >= 4n-8(d+1).

If the loop never starts, n<=d+1 and the claimed bound is nonpositive. This completes the proof for every n.

Equivalently, let U_0=n-I_0 and let S_8 be the actual sum of V-surpluses. The exact identity is

    |C| = 4n+4U_0+S_8+V_T-8a_T.

No retained gate or row is counted twice, and no charged object is created. All sixteen binary tables are covered by constant/literal cleanup and the XOR/AND classification. Sharing and fan-out are unrestricted throughout.


# Chapter III. Fixed-root restriction with disjoint original gate labels


Lemma. Normalize an acyclic single-output full-B2 circuit, deleting dead,
constant, unary, and repeated-pin gates and propagating constants/literals.
Let K_0 be a set of primary inputs fixed to constants. Suppose the final
restricted output is nonconstant. Let F be the number of incidences from
K_0 to gate input pins in the normalized circuit. Completing constant and
literal propagation removes at least F distinct normalized old gates.
If every member of K_0 is essential, F>=|K_0|.

Proof. Give every original gate a permanent label. Initial normalization
and pruning delete a set R_0 of labels; surviving gates keep their own
labels even when substitutions change their input names or truth tables.
Freeze incidences AFTER that normalization, which is essential: projection
contraction can increase a root's normalized fan-out. The savings R_0
will not be used in the following estimate.

Fix the roots and process the normalized DAG in topological order. Record
for every old signal its final canonical constant or signed signal alias.
A gate with a constant input, or two equal/complementary input signals,
is constant or unary and is removed. Retain an essential binary gate
otherwise. Do not dead-prune during this count; later pruning is harmless.

Let T_0 be the removed normalized gate labels. Let Z be the set of old
gate-output symbols whose recorded final aliases are constants. These are
syntactic constants from the stated substitutions, not merely functions
that happen to be semantically constant. Every v in Z has a unique old
producer g_v in T_0. For g in T_0, let r(g) count its two frozen input
pins whose source is either in K_0 or in Z. Let d(v) be the frozen
outdegree, counted by pins, of an old gate-output symbol v.

Each incidence leaving K_0 or Z ends at a removed gate: its final input
pin is constant. Conversely these are exactly the incidences counted by
r(g). Hence the exact identity is

    F + sum_(v in Z) d(v) = sum_(g in T_0) r(g).       (CP)

If out(g) is not in Z then r(g)<=1; two recorded constant parents would
make its output a recorded constant. Always r(g)<=2. Each v in Z has
d(v)>=1: the initially pruned circuit had no dead gate, and the selected
output cannot be in Z because its restricted value is nonconstant.
Consequently

    |T_0|-F
      = sum_(g in T_0, out(g) notin Z) [1-r(g)]
        +sum_(v in Z) [1-r(g_v)+d(v)]
      >=0.

All summands are nonnegative, and every gate appears once. Arbitrary
fan-out and sharing are already included in the incidence sums. An
essential root must have at least one use in the normalized circuit, so
F>=|K_0| when all roots in K_0 are essential. This proves the lemma.

Native-basis detail. Signed aliases are a bookkeeping device. A sign on
a gate input is absorbed into its binary table. A sign on the selected
output can be absorbed into its last retained gate. If the restricted
function is a raw input literal this final absorption may instead require
special treatment, so we do NOT infer an additive native-circuit bound
in that case. All applications below leave a function that is neither
constant nor an input literal. Its last gate therefore exists and the
retained circuit is native B2 with no extra charge. The core gate proof
in Chapter II also permits these signed aliases throughout.

The sets R_0, T_0, and the retained old labels are disjoint. Thus a lower
bound L on the final retained native circuit gives a lower bound L+F
on the original circuit. This is one actual restriction history, not an
addition of hypothetical branches, current-state ranks, or reused charges.


# Chapter IV. Elementary irreducible-modulus support

This supports odd degrees ell=2w+1 in the cubic extractor. It is not restricted to prime extension degrees.

## 1. Existence for every positive integer ell

Take a splitting field over F_2 of X^(2^ell)-X, obtained by successively adjoining roots of its irreducible factors. Its derivative is one, so it has 2^ell distinct roots. These roots are closed under addition and multiplication by Frobenius, and under inversion of nonzero roots, so they form a field K with exactly 2^ell elements.

Frobenius T:x->x^2 has T^ell=identity on K. Thus every orbit length divides ell. For every proper divisor d of ell, the points whose orbit length divides d are roots of X^(2^d)-X, and there are at most 2^d of them. If ell>=2, the union of all shorter orbits therefore has size at most

    sum_(d=1)^floor(ell/2) 2^d
       =2^(floor(ell/2)+1)-2 < 2^ell.

So some alpha has Frobenius orbit length exactly ell. Its orbit polynomial

    product_(j=0)^(ell-1) (X-alpha^(2^j))

has coefficients fixed by Frobenius, hence in F_2. Every nonzero F_2 polynomial vanishing at alpha vanishes at all ell distinct orbit elements, and therefore has degree at least ell by the root bound. The orbit polynomial is its monic irreducible minimal polynomial of degree ell. Degree one is immediate separately.

## 2. A polynomial-time certificate check

For a monic binary polynomial P of degree ell, the following are equivalent to irreducibility:

1. X^(2^ell) == X mod P.
2. For every prime divisor p of ell,

       gcd(P, X^(2^(ell/p))-X)=1.

If P is irreducible, its quotient is a field of size 2^ell, proving the first condition. Its root orbit has length ell and so cannot have length dividing ell/p, proving each gcd condition.

Conversely, the first condition makes P a divisor of the squarefree polynomial X^(2^ell)-X. Each irreducible factor of P has degree e dividing ell, because its roots have Frobenius orbit length e and lie in the root field above. If e<ell, choose a prime p dividing ell/e. Then e divides ell/p, so this factor also divides X^(2^(ell/p))-X, contrary to condition 2. Every factor therefore has degree ell; since deg P=ell, P is itself irreducible.

Factor the integer ell by trial division. Compute X^(2^j) mod P for j=1,...,ell by successive squaring, retaining the needed checkpoints ell/p, and run polynomial Euclidean gcd at these checkpoints. Naive polynomial arithmetic uses polynomially many bit operations in ell. No search for a valid modulus is needed when its coefficients are provided as a prefix.


## Contribution sources and submission scope

The complete gate argument in Chapter II develops the three-point pair
framework from submission 5e9fcf7f-8338-4729-9faa-a18994661af6 and our
coefficient-four proof submitted as bf313f35-d97b-4836-b1e8-30afd1157ac2.
The input-modulus method and algebraic starting point were developed in
submissions f4f2c464-3045-4c32-a6e4-32478c05f343 and
2633389c-f681-4e93-b3c8-a86793dcc031. The affine-system method follows
Mishra and Santhanam, ECCC TR26-260, 2026-10-19. All arguments needed here
are reproduced. Source judgments are not proof premises.

The two structured claims assert the same complete M=100 family theorem.
The guarantee is g(N)=4N, with zero log and loglog exponents. Its finite
projection is (4,1). No checker or experiment is declared: this is a
proof-only submission. The two proof files supply all mathematical evidence.
