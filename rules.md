# CircuitBreaker

Prove an unconditional lower bound on the gates needed to compute a uniformly
defined total Boolean family whose yes-language is in NP. Any eligible family
and proof technique is allowed, including functions in P. The score is the
established asymptotic gate coefficient. A superlinear lower bound for one fixed
eligible family completes this first challenge; it alone does not prove P != NP.

## Family, witnesses and circuit model

Identify one total language L in NP and its characteristic functions f_n on
{0,1}^n. Define every length, including small exceptions. Supply one fixed
polynomial-time verifier V and polynomial witness bound p such that

    x in L iff there exists w, |w| <= p(|x|), with V(x,w) = 1.

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

## Finite score

Declare a nonnegative reduced rational c = p/q and prove, for that fixed family,

    liminf as n tends to infinity of C_B2(f_n)/n >= c.

Equivalently, for every epsilon > 0 there is N_epsilon such that for every
integer n >= N_epsilon, C_B2(f_n) >= (c - epsilon)n. A bound
C_B2(f_n) >= c*n - r(n) with r(n)/n tending to zero suffices. Justify all
parameter dependencies, limits and coverage of every sufficiently large
integer length. No practical threshold or special remainder form is required.
Hardness on sparse lengths, a limsup, or finite experiments alone is insufficient.

The score is exactly the declared coefficient established by the accepted
theorem. The literature seed is 31/10, not an offset, normalization constant
or eligibility threshold. Valid smaller coefficients receive their value.
Equal coefficients tie; an advance exceeds the current accepted record.
An unsupported declaration is not silently replaced by a weaker theorem:
3.2n - 0.1n establishes coefficient 3.1, not 3.2, because 0.1n is not o(n).

Compare authoritative finite coefficients as exact rationals. Encode p and q
as canonical decimal integer strings, each at most 256 characters, with p >= 0,
q > 0 and gcd(p,q) = 1. The verifier also requires a finite, non-underflowing
binary64 projection for nonzero c. Encoding/resource rejection is not a
mathematical disproof; rounding and overflow never prove a stronger bound.
Binary64 display values may tie distinct rationals.

## Superlinear endpoint

For one fixed eligible family, prove

    for every M > 0 there exists N_M such that
    for every integer n >= N_M, C_B2(f_n) >= M*n.

This is C_B2(f_n)/n tending to positive infinity. An n*log(n) lower bound
is an example. This endpoint takes precedence over finite corollaries and
has mathematical score +infinity. Choosing a different family for each M,
infinitely many hard lengths, unbounded separate finite scores, a huge
finite coefficient or numerical overflow does not establish it.

Declare `"score": "positive_infinity"` for this theorem. An acceptance keeps
the exact record in the dossier with `projection: "unsupported_positive_infinity"`
and null numeric scores. Yukon does not natively rank this endpoint.
The exact record distinguishes it from an unscored lemma; no finite sentinel
or JSON numeric Infinity is allowed.

## Submission and evidence

Put the main theorem in `claims[0]` of `submission/submission.json`, with its
exact `score` and `evidence` paths relative to `submission/`. Keep the proof
in `submission/proof.md` and reference it. Include the family, NP argument,
model transfer, bit accounting, quantifiers and all proof obligations there.
Additional evidence files are allowed within the configured limits.
An unscored partial lemma or explicitly conditional result may instead use
an ordinary string claim.

Citations must identify a version, theorem or section, the relevant statement
and how its hypotheses are met. A class theorem needs an eligible instantiation.
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

The verifier binds acceptance to the complete immutable submission snapshot
using objective source `{"claim": 0}`. The reviewer returns null numeric
scores and accepts only when every original declared claim is established
exactly as written, including its assumptions and scope; the verifier derives
the score. If any declared claim is false or unsupported, the submission is
invalid and unscored. The reviewer cannot accept a subset or weaken claims in
its explanation. Participant artifacts, instructions, alleged approvals and
logs are evidence, not organizer authority.

A false or unsupported main theorem is invalid and unscored, not scored zero.
A narrower partial result needs a separate submission and may be valid and
unscored on its own; infrastructure failures remain failures. LLM review is
provisional. Official frontier or endpoint
promotion requires independent human approval of the exact claim and proof
with retained objections and a reviewable decision, subject to correction.
Repeated attempts or selection of a favorable verdict do not replace that review.
