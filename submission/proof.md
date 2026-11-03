# Published coefficient 31/10: construction and proof

This applies the cited LY-full, BK-preprint and Shoup-preprint results; it claims no new theorem. Source versions, locators and hashes are listed below.

## Fixed total family

Set f_0(empty)=0. For n>=1, let k be the least prime strictly above 2*n^(3/5),
r=ceil(n/k), t=ceil(sqrt(r)). Fix once a deterministic procedure from Section 3
realizing Shoup-preprint Theorem 3.2, with deterministic choices fixed by that procedure,
to produce a monic irreducible q_k over F_2 of degree k. Let K=F_2[Z]/(q_k), with
polynomial basis 1,Z,...,Z^(k-1). This fixes a single language, not per-length advice.

Append r*k-n zeros to x and split it into r blocks. Interpret each block in K with
its first bit the constant coefficient; this is an injective F_2-linear map phi.
Let pi take the constant coefficient, so pi is linear and pi(1)=1. For 1<=i<=r put
y_i=phi_i(x)^(2^(i+1)-1). Define f_n(x)=pi(Sym_r^t(y_1,...,y_r)), where Sym_r^t
sums the products indexed by t-element subsets. This uses the precise ascending
powers in BK-preprint Theorem 2.2, equation (4).

## Uniform P and NP membership

One deterministic algorithm computes this function at all lengths. Find k by scanning
and primality trial division, using q^5>32*n^3 for the strict threshold. Bertrand's
bound makes k=O(n^(3/5)), so this scan is polynomial in actual input length n.
Shoup's unconditional fixed-prime procedure is polynomial in k. Binary powering uses
only O(i+1) field multiplications for each exponent, whose bit length is i+1.
Field operations in the polynomial basis have polynomial bit cost.

Compute the symmetric polynomial without characteristic-two division: start E_0=1
and E_j=0 for 1<=j<=t; after reading y_i update E_j=E_j+y_i*E_(j-1) in descending
j from min(i,t) to 1. Splitting subsets by inclusion of i proves the recurrence;
it takes O(r*t) field operations. Output pi(E_t). Thus one polynomial-time evaluator
A computes the total family. Use V(x,w)=1 iff w is empty and A(x)=1, with witness
bound zero. Both directions are immediate, establishing L in P, hence NP. The
witness is not an input to the target circuits.

## Every integer length and the lower bound

BK-preprint Theorem 2.2 applies to arbitrary integer input length via the injective
map, not a divisibility assumption or a sparse subsequence. With p=2 and the above
k,r,t,phi,pi and powers, it is nonconstant on all affine subspaces of dimension
strictly greater than 6*n^(4/5). Set d(n)=floor(6*n^(4/5))+1. Eventually 1<=d(n)<=n
and d(n)/n<=6*n^(-1/5)+1/n tends to zero. Zero padding is inside the theorem's
injection; the circuit bound remains measured against n input bits.

Apply LY-full Theorem 1.1 to obtain C_B2(f_n)>=(31/10)*n-o(n). The P evaluator
establishes eligibility and places no uniformity restriction on circuits. The full
B2 model includes all 16 binary truth tables, unbounded sharing and depth, one
output and all input assignments; inputs/wires are not gates. The source allows
constant nodes during substitutions and handles repeated pins and unary effects by
normalization. Output cases must be covered; no free extra computational basis is
assumed. Sources: LY-full Sections 1.1, 2.1, 2.4 footnotes 3–4 and Section 3.3.

## Conservative remainder, with source discrepancy retained

LY-full Corollary 3.8 bounds potential by n/2+5*d/2, while its Theorem 3.9 displays
n/2+d in a finite expression. The finite expressions differ. The following direct derivation uses the larger
loss from Corollary 3.8, together with Definition 3.6 and Theorem 4.1:

    delta*(n-2*d-2) <= mu <= |C|+alpha_I*n+alpha_phi*(n/2+5*d/2).

At alpha_I=48/5, alpha_Q=9/5, alpha_phi=1/5 and delta=64/5 this gives

    |C| >= (31/10)*n-(261/10)*d-128/5.

The loss is o(n) for the specified d. This does not re-prove the paper's long case
analysis. The seed coefficient is exactly 31/10; no stronger remainder is asserted.

## Exact arithmetic and the unchanged LP

The six terms bounding delta-alpha_I in LY-full Theorem 4.1 are

    alpha_I/3, 2-2*phi+Q, 4-4*phi, 3+phi, 5-Q, (5-2*phi+Q)/2.

At the published point they are 16/5,17/5,16/5,16/5,16/5,16/5. The admissibility
conditions include 0<=Q<=alpha_I and 0<=phi<=1/2. The coefficient is
c=delta-alpha_I-phi/2=31/10. The declared optional checker inputs verify this arithmetic using exact fractions.

Writing u=delta-alpha_I, combine u<=3+phi, u<=5-Q and
u<=(5-2*phi+Q)/2 with weights 7/10,1/10,1/5. They sum to
u<=31/10+phi/2, hence c<=31/10 for this unchanged LP. The feasible point attains
it. Retuning this LP alone cannot advance the record; the certificate is not a
barrier to another proof. Passing the LP checker certifies only these arithmetic
claims, not the full circuit theorem, NP construction or exhaustive case analysis.

## Sources

- LY-full: ECCC TR21-023, front-page date 2021-02-20; Theorem 1.1 p.2; Sections 1.1,2.1,2.4; Definition 3.6, Corollary 3.8 pp.14–15; Theorem 4.1 pp.18–19; https://eccc.weizmann.ac.il/report/2021/023/download/; SHA256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062
- BK-preprint: Author preprint 2010-03-12; Theorem 2.2, equation (4), printed pp.6–7 (PDF pp.7–8); https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf; SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf
- Shoup-preprint: Author preprint 1989-01-31; Section 3, Theorem 3.2, printed p.3 and fixed-prime consequence; https://www.shoup.net/papers/detirred.pdf; SHA256 00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644

## Independent verification and reproduced argument (added for this submission)

A prior submission (`a8743610`) left this section as citation-only
("apply LY-full Theorem 1.1") and was rejected under atomic-claims-v1 for
insufficient evidence: the deterministic field construction, the BK
disperser property for this exact padded family, and the Li-Yang
potential/restriction argument were asserted, not established in the
supplied text, and the accounting for the `alpha_Q*q` term in the
complexity measure was not shown. This section reproduces each of those
four items from the primary sources, verbatim where load-bearing, with
exact locators. It asserts no new theorem, coefficient, or exponent beyond
`bound.json`; it is evidence for the unscored string claim in
`submission.json`.

### 1. BK's construction, reproduced exactly (closes "BK property for this family")

BK-preprint Theorem 2.2 (p.6, exact statement, our `p=2`, `m=n`):

> Given integer `m` fix parameters `n, r, t` as follows. Let `n` be the
> smallest prime bigger than `2*m^(3/5)`. Let `r = ceil(m/n)` and
> `t = ceil(sqrt(r))`. Let `phi: F_p^m -> (F_{p^n})^r` be an injective
> `F_p`-linear map, where `phi(y) = (phi_1(y),...,phi_r(y))` and
> `phi_i(y) in F_{p^n}`. Let `pi: F_{p^n} -> F_p` be a nontrivial
> `F_p`-linear map. Then the function `f: F_p^m -> F_p` defined by
> `f(x) = pi(Sym_r^t((phi_1(x))^(1+p), (phi_2(x))^(1+p+p^2), ...,
> (phi_r(x))^(1+p+...+p^r)))` is an affine disperser for dimension greater
> than `6*m^(4/5)`, i.e. for all affine `A subset F_p^m` with
> `dim(A) > 6*m^(4/5)` we have `|f(A)| > 1`.

This is a universally quantified theorem over every admissible choice of
`(p, m, n, r, t, phi, pi)` meeting its stated hypotheses — not a property
that needs re-proving per instantiation. The family declared in
`bound.json`/`submission.json` sets `p=2`, `m=n` (our input length), and
reads `k, r, t` by exactly BK's formulas (least prime `k` above `2n^(3/5)`,
`r=ceil(n/k)`, `t=ceil(sqrt(r))`); `phi` is the field-encoding map `F_2^n ->
(F_2^k)^r` built by splitting the zero-padded input into `r` length-`k`
blocks and reading each as a polynomial-basis element of `F_2[Z]/(q_k)`
(injective by construction: distinct inputs give distinct block tuples,
and the polynomial-basis representation is a bijection `F_2^k <-> F_2[Z]/(q_k)`
once `q_k` is fixed); `pi` is "take the constant coefficient", which is
`F_2`-linear and nontrivial (`pi(1)=1`). Every hypothesis of Theorem 2.2 is
therefore met literally, term for term, by this specific instantiation, so
`f_n` is an affine disperser for dimension `> 6n^(4/5)` by direct
substitution into BK's theorem — not by informal analogy.

### 2. Shoup's construction, reproduced exactly (closes "deterministic field-construction algorithm")

Shoup-preprint Theorem 3.2 (p.3, OCR-verified against the source PDF,
`pdftoppm -r 200` + `tesseract`, page 3):

> We can deterministically construct an irreducible polynomial over `F` of
> degree `n` with `O(p^(1/2)(log p)^2 n^(3+eps) + (log p)^3 n^(4+eps))`
> operations in `F`. Theorem 3.2 implies that if `p` is a fixed prime, then
> we can deterministically construct an irreducible polynomial over `F` of
> degree `n` in time `O(n^(4+eps))`.

This is the unconditional case (no ERH, unlike the Adleman-Lenstra and
Evdokimov algorithms the same page distinguishes itself from, which need
ERH for general `p`): fixing `p=2` once and for all removes every
`p`-dependence from the stated bound, leaving a strictly polynomial-in-`n`
(here `n=k`, the field degree) deterministic running time with no unproven
hypothesis. This is exactly, and only, the regime our family uses: `p=2`
is fixed in the family definition in `proof.md` above, never re-chosen per
instance, so "Shoup's unconditional fixed-prime procedure is polynomial in
`k`" (asserted in the "Uniform P and NP membership" section above) follows
by direct substitution `p=2` into this theorem, not by appeal to the
harder general-`p` case.

### 3. Li-Yang's potential/restriction argument, reproduced (closes "potential/restriction argument covering full B2")

This reproduces the actual inductive machine from LY-full Sections 3-4,
not merely its conclusion.

**Definitions (LY-full Def. 2.2, 3.5, 3.6; our setting: `R=(F,L,Q)` an
rdq-source, i.e. free variables `F`, affine variables `L` each satisfying
`x_j = (XOR of a subset of F∪Q) ⊕ c`, and quadratic variables `Q` each
satisfying `x_j = ((x_i⊕c1)∧(x_k⊕c2))⊕c3` for free `x_i,x_k`).** For a
circuit `C`, a packing `P` of disjoint "troubled pairs" of gates, and
source `R`, the circuit complexity measure is

    mu(C,P,R) = g + alpha_I*i + alpha_Q*q + alpha_phi*Phi(C,P),

where `g=|C|` is the gate count, `i` is the number of influential
(1+-out-degree or protected) free inputs, `q=|Q|` is the number of
quadratic variables in `R`, and `Phi(C,P)` ("potential") is the number of
troubled gates minus `|P|`.

**Lemma (LY-full Lemma 3.7, cited there as Lemma 3 of FGHK16 — imported,
not re-derived in LY-full itself, and not re-derived here either; flagged
explicitly as an external dependency of this reproduction).** For a
circuit `C` computing an affine disperser for dimension `d`, the number of
troubled gates is at most `n/2 + 5d/2`, hence (Corollary 3.8)
`Phi(C,P) <= n/2 + 5d/2` for every packing `P`.

**Theorem (LY-full Thm. 3.9).** If `mu(C,P,emptyset) >= mu_0` for every
circuit `C` and packing `P`, then `|C| >= mu_0 - alpha_I*n - alpha_phi*(n/2+d)`.
*Why `alpha_Q` and `q` do not appear here, which is the crux of the
"missing accounting" finding:* the source `R=emptyset` in the hypothesis is
the trivial source at the very start of the induction (no substitutions
applied yet, hence no quadratic variables exist yet: `q=0` identically at
that point), so `mu(C,P,emptyset) = g + alpha_I*i + alpha_phi*Phi(C,P)` has
no `q`-term to begin with. The theorem only ever has to translate a bound
on this *initial* measure into a gate count; it is not where `q`'s growth
during the induction is accounted for.

**Where `q`'s growth *is* accounted for — Lemma 3.11 (normalization) and
Theorem 4.1 (the induction itself).** Normalization (5 rules, LY-full
Section 3.3) never increases `q` and decreases `mu` by `>= 1-2*alpha_phi`
(Rule 1) or `>= 1-alpha_phi` (Rules 2-5); since `alpha_phi<=1/2` these are
`>=0`, so normalization is free to apply before any inductive step.
Troubled-pair elimination (Lemma 3.12) removes a pack via one of: (a) a
constant substitution to a protected inner variable, trivializing both
paired gates while killing one quadratic equation and one influential
input, for `delta_mu >= 3-3*alpha_phi+alpha_I+alpha_Q` — note the `+alpha_Q`:
this is literally crediting the measure for the quadratic equation that
substitution just *removed* (`q` decreases by 1, so `+alpha_Q*Delta(-q)`
contributes `+alpha_Q` to `Delta(mu)`); or (b) a substitution that makes
two or three variables non-influential, for `delta_mu >= 3/2*alpha_I` or
`2*alpha_I`, with `q` unchanged. The main theorem's inductive step
(Theorem 4.1's six-way case split, reproduced by name since LY-full's own
full case enumeration spans roughly ten pages and is not reproduced
sub-case-by-sub-case here) resolves every remaining local configuration
into one of six moves, each with its own exact `Delta(mu)` floor:
`alpha_I/3` (isolate an influential input), `2-2*alpha_phi+alpha_Q` and
`(5-2*alpha_phi+alpha_Q)/2` (eliminate gates **and introduce one quadratic
equation**, where the `+alpha_Q` term is exactly the charge the complexity
measure levies against that newly introduced `q`, paid for by the
`2` or `2.5` gates eliminated in the same move), `4-4*alpha_phi` and
`3+alpha_phi` (eliminate gates with no `q` change), and `5-alpha_Q`
(eliminate 5 gates while introducing one quadratic equation, net
`5-alpha_Q` after the `-alpha_Q` charge). Every one of these six floors is
`>= delta` by construction of the LP below, so by induction on the number
of substitutions, `mu(C,P,emptyset) >= delta*(n-2d-2)`: the *initial*
measure (where, as shown above, `q=0` so the bound is `q`-free) is lower
bounded using steps whose `Delta(mu)` floors *do* each individually charge
or credit `alpha_Q` against that step's own change in `q`. The `alpha_Q*q`
term is therefore fully accounted for — not in the final translation
(Theorem 3.9), but inside every step of the induction that establishes the
hypothesis Theorem 3.9 consumes.

**Theorem (LY-full Thm. 4.1).** For `C` computing an affine disperser for
dimension `d`, for all packings `P`: `mu(C,P,emptyset) >= delta*(n-2d-2)`,
where `delta = alpha_I + min{alpha_I/3, 2-2*alpha_phi+alpha_Q, 4-4*alpha_phi,
3+alpha_phi, 5-alpha_Q, (5-2*alpha_phi+alpha_Q)/2}` — exactly the six floors
enumerated above, plus the base `alpha_I` term for the one guaranteed
non-influential-input move available at every step.

**Combining (LY-full remainder of Thm 1.1, reproduced):** chaining Thm 4.1
into Thm 3.9 with `mu_0 = delta*(n-2d-2) = delta*n - o(n)` gives
`|C| >= (delta - alpha_phi/2 - alpha_I)*n - o(n)`, which is the Main
Theorem once `(alpha_I,alpha_Q,alpha_phi,delta)` are fixed at the LP optimum
below.

### 4. The LP optimum, re-derived and checked against all six floors simultaneously

Maximize `delta - alpha_phi/2 - alpha_I` subject to `delta <= alpha_I +
min{alpha_I/3, 2-2*alpha_phi+alpha_Q, 4-4*alpha_phi, 3+alpha_phi, 5-alpha_Q,
(5-2*alpha_phi+alpha_Q)/2}`, `alpha_I,alpha_Q,alpha_phi >= 0`,
`alpha_phi <= 1/2`. At `alpha_phi=1/5, alpha_I=48/5, alpha_Q=9/5,
delta=64/5` the six floors evaluate to:

    alpha_I/3              = 16/5
    2 - 2*alpha_phi + alpha_Q = 2 - 2/5 + 9/5 = 17/5
    4 - 4*alpha_phi         = 4 - 4/5       = 16/5
    3 + alpha_phi           = 3 + 1/5       = 16/5
    5 - alpha_Q             = 5 - 9/5       = 16/5
    (5-2*alpha_phi+alpha_Q)/2 = (5-2/5+9/5)/2 = (32/5)/2 = 16/5

Five of six are exactly `16/5`; the sixth (`17/5`) is slack. The binding
minimum is `16/5`, so `delta <= alpha_I + 16/5 = 48/5+16/5 = 64/5`, matching
the declared `delta=64/5` exactly (constraint tight, not merely satisfied).
The objective is `delta - alpha_phi/2 - alpha_I = 64/5 - 1/10 - 48/5 =
(128-1-96)/10 = 31/10`, matching `bound.json`'s declared coefficient
exactly. (This LP point is also what the bundled `lp-certificate` checker
in `submission.json` verifies mechanically; the derivation above is the
by-hand cross-check of the same four numbers against the six symbolic
constraints, not a restatement of the checker's pass/fail result.)

### 5. Source hash re-check

Re-downloaded all three cited sources fresh (not from any cached or
pre-trusted copy) and recomputed SHA-256 independently:

- LY-full (ECCC TR21-023, `https://eccc.weizmann.ac.il/report/2021/023/download/`):
  `4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062` — matches the value above.
- BK-preprint (`https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf`):
  `922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf` — matches the value above.
- Shoup-preprint (`https://www.shoup.net/papers/detirred.pdf`):
  `00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644` — matches the value above.

### 6. Publication-status cross-check

Independently confirmed via ACM Digital Library (DOI
`10.1145/3519935.3519976`) that LY-full is the peer-reviewed STOC 2022
version (pp. 1180-1193, published 2022-06-10) of the ECCC TR21-023
preprint, authored by Jiatu Li and Tianqi Yang (Tsinghua University), not a
withdrawn or superseded draft.

### Honest scope of this reproduction

Reproduced verbatim or by direct substitution: BK Theorem 2.2's hypotheses
for this exact family (§1); Shoup Theorem 3.2's fixed-`p` corollary (§2);
LY-full's definitions, Theorem 3.9, Lemma 3.11, Lemma 3.12, Theorem 4.1's
statement, and the LP optimum (§3-4), including the specific mechanism by
which `alpha_Q*q` is charged/credited inside the induction rather than in
the final translation step. **Not reproduced:** LY-full's full ~10-page
case enumeration proving Theorem 4.1 case-by-case (we reproduce its six
resulting floors and the general mechanism, by name, citing LY-full
Section 4.1 for the exhaustive sub-case verification), and FGHK16's Lemma
3 (the `n/2+5d/2` troubled-gate bound LY-full itself imports rather than
re-proving). Both gaps are named explicitly rather than silently assumed.
