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
The planner response schema enforces the same one-to-eight region bound as
runtime validation. Failed reviews retain separate bounded diagnostics in
`llmv-initial-review-failure` or `llmv-final-review-failure`; those artifacts are
not judgments and cannot produce a score.

Each review prompt has a 2 MiB (2,097,152-byte) ceiling, matching the pinned
verifier's default. Planning includes the original evidence and numbered source;
defense and adjudication also include accumulated reports. The former 128 KiB
ceiling rejected an existing seven-file submission before its first model call.
The ceiling remains enforced without truncation. It does not guarantee that every
possible combination of source text, reports and responses fits; the separate
response, request, token and time limits still apply.

After source edits, run `python tools/generate_hosted.py` to apply the upstream
generator plus the protected metrics/test overlays and validate the companion
presentation against the generated objectives and metric contract. Page settings
remain in the separate template, not in `benchmark.json`. The hosted runner selects
the metrics-bearing `circuitbreaker-score` artifact; local `finalize_score.py`
wraps shared finalization with the same export. See [migration.md](migration.md)
for replacement-registration inputs, the configured review budget, and binding
both views with Yukon's existing page-preparation helper. Registration queues
baseline review immediately; there is no model-free create-paused mode.
Hosted artifact downloads retry transient network/HTTP failures at most four
attempts per operation within a shared three-minute transfer deadline. Inputs
are selected by immutable artifact ID and verified archive/file hashes. The
metrics job uses the original review's wheel, intake and dossier and emits a
separate `circuitbreaker-export-lineage` artifact binding those inputs, the source
and config, the score bytes, and distinct review/delivery attempts.

The runner declares deterministic baseline recovery for this metrics job. A
compatible Yukon deployment can validate a later metrics-only delivery without
repeating the review. Only the two derived score/lineage artifacts use replacement
uploads after a partial delivery; original review inputs remain immutable.
Recovery checkpoints freeze the actual delivered artifact IDs and hashes. Older
runs without this declaration are not eligible for mixed-attempt recovery.

The verifier pin is `1f61aa60c420e65c28786b56500b4e9276b5f53b` (version 0.7.1); the optional
structured-claims protocol remains disabled. This source change does not alter
the existing live benchmark or page.

<!-- llm-verifier:solver-guidance:begin -->
## Preparing a submission

Read the [solver guide](https://github.com/Layr-Labs/llm-verifier/blob/1f61aa60c420e65c28786b56500b4e9276b5f53b/docs/solving.md) and
[solver skill](https://github.com/Layr-Labs/llm-verifier/blob/1f61aa60c420e65c28786b56500b4e9276b5f53b/skills/llm-verifier-solver/SKILL.md)
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
