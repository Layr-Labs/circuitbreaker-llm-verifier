# CircuitBreaker

Prove a stronger unconditional circuit lower bound for one fixed, uniformly
defined total Boolean family in NP. The [rules](rules.md) define scientific
eligibility, joint coefficient/exponent scoring and explicit endpoints.
Native finite scores are `[c,p]`, both maximized. The full accepted expression
is authoritative; the graph is only a projection. Endpoint scores are null
and their typed expression metrics remain available. The additive
`circuitbreaker_frontier` envelope supports the registered Complexity table.
[The complete presentation template](examples/circuitbreaker-presentation.json)
combines that table with coefficient history grouped by the selected power of n.
It uses the merged Yukon contract from `f340c97dad9f170941e03c233a8a263ae76b8998`;
compatible API and UI releases must be running before it is selected.

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
python tools/generate_frontier_fixtures.py --check
llm-verifier test-checkers --root .
```

Prepare the digest-pinned Python Docker image before the last command. These
checks do not perform a mathematical review. [Scientific cases](scientific-cases.md)
explain known valid and invalid reasoning; synthetic exporter tests assert no
new circuit result. The proof byte budget, per-call provider configuration, LP
checker and scientific family/model conventions are preserved.

Review uses `builtin:math-decomposed-v1`: a planner partitions all submitted
source into at most eight regions, a whole-proof reviewer and independent region
workers assess the original claims, and a defender/adjudicator handles findings.
This profile adapts the offline review pattern from the UCLA Moonshot Math
Verifier; it does not add external retrieval or change atomic claim acceptance.
Every file needs at least one region, so **submit at most eight files**, including
`submission.json`, `bound.json`, proof files and any experiment code. This profile
limit applies even though the generic `max_files` ceiling remains 32. The baseline
uses three files. Keep supporting material together without dropping evidence.

After source edits, run `python tools/generate_hosted.py` to apply the upstream
generator plus the protected metrics/test overlays and validate the companion
presentation against the generated objectives and metric contract. Page settings
remain in the separate template, not in `benchmark.json`. The hosted runner selects
the metrics-bearing `circuitbreaker-score` artifact; local `finalize_score.py`
wraps shared finalization with the same export. See [migration.md](migration.md)
for replacement-registration inputs, the configured review budget, and binding
both views with Yukon's existing page-preparation helper. Registration queues
baseline review immediately; there is no model-free create-paused mode.
The verifier pin is `5c3e66acffc95aa5a9cfc31183dbefc7ed24b9bf` (version 0.7.0); the optional
structured-claims protocol remains disabled. This source change does not alter
the existing live benchmark or page.

<!-- llm-verifier:solver-guidance:begin -->
## Preparing a submission

Read the [solver guide](https://github.com/Layr-Labs/llm-verifier/blob/5c3e66acffc95aa5a9cfc31183dbefc7ed24b9bf/docs/solving.md) and
[solver skill](https://github.com/Layr-Labs/llm-verifier/blob/5c3e66acffc95aa5a9cfc31183dbefc7ed24b9bf/skills/llm-verifier-solver/SKILL.md)
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
