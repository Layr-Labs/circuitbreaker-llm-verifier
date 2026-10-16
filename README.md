# CircuitBreaker

Prove a stronger unconditional circuit lower bound for an explicit Boolean
family in NP. The [rules](rules.md) define eligibility and scoring.
The starting [proof](submission/proof.md) instantiates the published 31/10
coefficient using a uniform polynomial-time family.

Edit `submission/`: state your theorem and exact coefficient in the first
claim of [submission.json](submission/submission.json), and reference your
proof file in its `evidence` list. Define the family and establish its NP
membership, every theorem hypothesis and the claimed asymptotic bound.
Every original declared claim must be established exactly as written for
acceptance. A narrower lemma needs a separate submission.

Check the submission structure with the pinned llm-verifier installation:

```sh
llm-verifier check --root .
```

If using the optional exact LP checker, declare its four rational parameters
in `checks`. Its [known-answer cases](checker-tests.json) run with the existing
digest-pinned Docker image in a process without model credentials:

```sh
llm-verifier test-checkers --root .
```

[Scientific cases](scientific-cases.md) illustrate accepted and rejected
mathematical claims. Local checks do not supply a mathematical review.
