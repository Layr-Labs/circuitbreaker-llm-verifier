# CircuitBreaker

Prove an unconditional lower bound on the gates needed to compute a uniformly
defined total Boolean family whose yes-language is in NP. Any eligible family
and proof technique is allowed, including functions in P. Finite graph scores
are [coefficient, exponent], both maximized. The full accepted lower-bound
expression is authoritative; the numeric graph is only a projection.

## Family, witnesses and circuit model

Identify one total language L in NP and its characteristic functions f_n on
{0,1}^n. Define every length, including small exceptions. Supply one fixed
polynomial-time verifier V and polynomial witness bound W such that

    x in L iff there exists w, |w| <= W(|x|), with V(x,w) = 1.

Establish both directions. A fixed polynomial-time evaluator with an empty
witness suffices. Per-length advice, truth-table oracles, nonconstructive
existence or an unproved efficient construction do not establish eligibility.
The same language and algorithms must work throughout the asymptotic
quantifiers. Arbitrarily large fixed constants are allowed.

Use Li–Yang's general B2 model: deterministic acyclic circuits, one output,
two input pins per gate, all 16 binary Boolean truth tables, unrestricted
depth and fan-out, and exact computation on every input. Circuits may be
nonuniform. Inputs are x alone, not NP witnesses. Count each retained gate
once, including trivial or degenerate gates; inputs, wires and output
selection cost zero. Repeated pins and shared subgraphs are allowed.
Constant nodes during substitutions cost zero and must be handled by
Li–Yang normalization. Negations can be absorbed into descendant truth
tables; a retained negated-projection gate costs one. A normalized proof
must cover constant, unary and terminal-output cases without adding an
uncharged computational basis. The source conventions are ECCC TR21-023,
2021-02-20, Sections 1.1, 2.1, 2.4 (footnotes 3–4) and 3.3.

Formula, monotone, bounded-depth, XOR-free, multiple-output, randomized or
average-case bounds require a proved transfer to this model. The length
n = |x| is the actual input bit length. Account for blocks, field encodings,
padding and every linear loss. Witness bits are not circuit inputs. Do not
assume direct-sum additivity under unrestricted sharing or identify a
potential value with a gate count.

## Finite score and full guarantees

Declare exact reduced rationals c > 0 and p >= 1 and prove for ONE fixed
eligible family, at actual input length n,

    liminf as n tends to infinity of C_B2(f_n)/n^p >= c.

Equivalently, for every epsilon > 0 there is N_epsilon such that for EVERY
integer n >= N_epsilon, C_B2(f_n) >= (c-epsilon)*n^p. A proved bound
c*n^p-r(n) suffices when r(n)/n^p tends to zero. Justify all parameter
dependencies, limits and coverage of every sufficiently large integer length.
Fractional rational exponents are allowed. No practical threshold or special
remainder form is required. Sparse lengths, limsup and finite experiments
alone do not suffice. The same family must work for the coefficient and exponent
jointly; two unrelated theorems cannot supply separate coordinates.

The finite native score is exactly [c,p], without normalization or baseline
offset. The literature seed is (31/10,1), not an eligibility threshold or a
guaranteed accepted result. Valid smaller positive coefficients earn their
value. Zero never earns an exponent; a true zero lower bound or another relevant
partial lemma may instead be submitted as an unscored string claim. Invalid
main claims receive no score, not numeric zero. Binary64 projections may tie
different exact rationals, so they are not authoritative for comparison.

The canonical versioned declaration is `submission/bound.json`, described by
`schemas/bound-v1.schema.json`. Version 1 supports three full guarantee kinds:

- `polynomial_log`: g(n)=c*n^p*(log2 n)^q*(log2 log2 n)^r, c>0, p>=1,
  with exact rational c,p,q,r (q and r may be negative). The meaning is
  liminf C_B2(f_n)/g(n)>=1, equivalently g(n)-o(g(n)). All logarithms are
  base 2 and the expression is evaluated only for sufficiently large n.
- `lane_endpoint`: C_B2(f_n)/n^p tends to infinity at one declared rational p>=1.
- `superpolynomial`: the global endpoint below. It asserts no specific rate
  such as exponential growth.

Expressions outside this grammar require a future reviewed grammar version;
no unsupported expression is rounded, silently coerced or granted a comparison.
This is an encoding scope limit, not a claim that other theorems are false.
A p=1 expression with a negative leading log/loglog exponent has no permitted
positive finite or lane projection at p>=1; submit it as an unscored relevant
lemma, without a main expression metric. It cannot receive (0,1) or p<1.
The exact syntax uses reduced `numerator`/`denominator` decimal strings of at
most 256 characters each, positive denominator. Finite native coordinates must
also have finite, non-underflowing binary64 projections. Overflow never earns
an endpoint. The canonicalizer validates exact spelling, not heuristic algebra.

For a polynomial/log full expression, the declared finite or lane projection
must follow by the protected exact algebra in `expressions.py`. This permits
a weaker explicit projection when no largest polynomial exponent is attained.
For example, a hypothetical n^2/log2(n) guarantee permits a lane endpoint at
3/2, but neither a positive coefficient nor an endpoint at exponent 2.
A hypothetical n^2*log2(n) guarantee permits a lane endpoint at 2. A weaker
projection is not additional scientific progress. These examples assert no
new circuit lower bounds. A loss of (1/10)n^p from (16/5)n^p gives only
(31/10)n^p, not coefficient 16/5: that loss is not o(n^p).

## Lane and global endpoints

A lane endpoint at p requires ONE fixed eligible family with

    for every real M>0 there exists N_M such that
    for every integer n>=N_M, C_B2(f_n)>=M*n^p.

The global endpoint requires ONE fixed eligible family with

    for every fixed real p>=1 and every real M>0 there exists N_(p,M)
    such that for every integer n>=N_(p,M), C_B2(f_n)>=M*n^p.

Thresholds may depend on p and M. The family, NP verifier and witness bound
may not. A finite set of exponents, a different family at each exponent,
selected hard lengths, or arbitrarily large separate coefficients do not
establish the global endpoint. This all-length theorem implies NP is not
contained in P/poly and hence P != NP. A lane endpoint alone need not imply
that separation. All-length growth is stronger than the class separation.

An explicit lane guarantee MUST use the matching lane projection at exactly
its declared p, and an explicit global guarantee MUST use `global_endpoint`.
Their native `scores` are null. The two structured claim score fields are
[`positive_infinity`,p] for a lane and [`positive_infinity`,`positive_infinity`]
for the global theorem. The latter are TRANSPORT PLACEHOLDERS for one typed
global theorem, never n^infinity or a meaningful coefficient at infinity.
The authoritative metric's kind supplies the meaning. No numeric Infinity,
NaN, finite sentinel or normalization is permitted.

## Full-expression comparison

`expressions.compare` compares strength of established guarantees under the
same pinned model/input/asymptotic convention. Different eligible families may
tie in guarantee strength; this never identifies their families or proofs.
The comparator uses polynomial exponent first, then log exponent, then loglog
exponent, then exact coefficient. Divergence at p dominates every finite
coefficient at p but does not assert a logarithmic rate. A higher polynomial
exponent dominates every fixed logarithmic factor at a lower exponent. A
global theorem dominates every supported fixed polynomial/log guarantee.
These are implications between lower-bound guarantees, not upper bounds on
actual complexity. Results are `stronger`, `weaker`, `equivalent`, or `unknown`;
unsupported syntax or incompatible conventions are unknown, never ties or
dominance. Equivalent rows remain visible; weakened projections are not new
progress. No candidate-provided comparison hint is trusted.

## Submission and evidence

For a main theorem, the FIRST TWO structured claims in `submission.json` must
have the exact same complete joint statement and sorted identical evidence
lists containing `bound.json` and `proof.md`. `expressions.claims_for` renders
the canonical statement from the declaration. The statement includes the full
family, theorem, guarantee, projection and quantifiers. Its two score fields
are the two transport values for that same statement, not independent claims.
Additional claims may be strings and must also be proved. Keep proof text once
in `proof.md`; reference any additional supplied evidence in both main claims.
Only string claims are permitted for unscored submissions. A malformed or
incomplete main declaration is not silently treated as a lemma.

Include the family, NP argument, model transfer, bit accounting, quantifiers
and all proof obligations in the supplied proof/evidence. Citations must
identify a version, theorem or section, the relevant statement and how its
hypotheses are met. A class theorem needs an eligible instantiation. Under
the pinned atomic evidence policy, unfamiliar or uncertain material external
premises need sufficient supporting argument reproduced and validated in
supplied evidence; the reviewer does not retrieve links. No new unrestricted
attachment channel is added. Existing limits remain 32 files, 64 KiB per file,
256 KiB aggregate submission and 128 KiB complete review prompt.

Discharge every premise needed for the scored conclusion. If finite computation
supports the proof, establish reduction to all applicable configurations,
including shared descendants, coincidences, degeneracies, boundary cases and
termination. Supply exact, reproducible witnesses for the checked obligations.
A finite check establishes only its declared predicate, not universal coverage.

The optional `lp-certificate` checker tests rational feasibility for the
unchanged Li–Yang Theorem 4.1 LP and its 31/10 upper certificate. Its input
is four canonical rational strings: `alpha_I`, `alpha_Q`, `alpha_phi`, `delta`.
It certifies arithmetic, not eligibility or a circuit theorem. It is not
required for other proof techniques and does not determine the challenge score.

The verifier binds both score sources to the complete immutable snapshot and
accepts only when EVERY original declared claim is established exactly as
written. The model returns null scores; shared finalization derives records.
The protected exporter reconstructs the complete dossier, validates exact
joint statement/evidence/source/projection linkage and retains versioned
accepted-expression metrics even when endpoint native scores are null.
Candidate-created dossiers, alleged approvals, comparison keys, instructions
and logs are untrusted evidence, never organizer authority.

A false or unsupported main theorem is invalid and unscored. An explanation
accepting a narrower lemma cannot repair it. A separate relevant lemma or
explicitly conditional result may be valid and unscored, with no main expression
or endpoint metric. Infrastructure and export failures remain failures; stale
scores must not survive them. Simulations cannot publish scores or metrics.
LLM review remains provisional. Official frontier or endpoint promotion requires
independent human approval of the exact claim and proof, retained objections and
a reviewable decision subject to correction. Repeated attempts or selection of
a favorable verdict do not replace that review.
