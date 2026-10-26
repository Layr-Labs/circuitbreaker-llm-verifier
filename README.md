# CircuitBreaker

Prove a stronger unconditional circuit lower bound for one fixed, uniformly
defined total Boolean family in NP. The [rules](rules.md) define scientific
eligibility, joint coefficient/exponent scoring and explicit endpoints.
Native finite scores are `[c,p]`, both maximized. The full accepted expression
is authoritative; the graph is only a projection. Endpoint scores are null
and their typed expression metrics remain available.

The [starting proof](submission/proof.md) is the unchanged published 31/10
baseline, represented as (31/10,1). It is a literature seed, not a guarantee of
acceptance under the current supplied-source evidence policy.

Edit `submission/bound.json`, then render its two identical joint claim
statements with `expressions.claims_for(bound, evidence)`. Preserve one copy of
each proof and a sorted identical evidence list in both claims. Each coordinate
must belong to that same complete theorem; two independent lemmas cannot be
combined. Put a relevant unscored partial result in string claims instead.
See [the expression/metric contract](metrics.md) and its JSON schemas.

Use the pinned llm-verifier installation for credential-free validation:

```sh
python check_submission.py
python -m unittest -v
python tools/generate_hosted.py --check
llm-verifier test-checkers --root .
```

Prepare the digest-pinned Python Docker image before the last command. These
checks do not perform a mathematical review. [Scientific cases](scientific-cases.md)
explain known valid and invalid reasoning; synthetic exporter tests assert no
new circuit result. The proof budget, provider configuration, LP checker and
scientific family/model conventions are preserved.

After source edits, run `python tools/generate_hosted.py` to apply the upstream
generator plus the protected metrics/test overlays. The hosted runner selects
the metrics-bearing `circuitbreaker-score` artifact; local `finalize_score.py`
wraps shared finalization with the same export. See [migration.md](migration.md)
for the required new registration and later table integration. This source
change does not alter the existing live benchmark.

<!-- llm-verifier:solver-guidance:begin -->
## Preparing a submission

Read the [solver guide](https://github.com/Layr-Labs/llm-verifier/blob/8c1f098875fae0750c29ecfd641f7c9fc53a0f41/docs/solving.md) and
[solver skill](https://github.com/Layr-Labs/llm-verifier/blob/8c1f098875fae0750c29ecfd641f7c9fc53a0f41/skills/llm-verifier-solver/SKILL.md)
from the same pinned verifier revision as this challenge.

- Use the challenge's verifier revision for `llm-verifier check --root .`.
- Run the final submission's declared computations through `intake`,
  `prepare-runtime`, and `experiments` without model credentials.
  `test-checkers` runs organizer fixtures, not your candidate's checks.
- Ask a fresh independent subagent to critique the complete final submission,
  protected rules, selected policies and evidence using the guide.
- Fix findings and rerun affected checks and review.
- Report results and omitted steps. Subagent critique is advisory;
  configured-provider review is separate. `simulate` does not judge correctness.
<!-- llm-verifier:solver-guidance:end -->
