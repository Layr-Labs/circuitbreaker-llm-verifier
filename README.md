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
new circuit result. The exact seven-file candidate fixture under `tests/fixtures`
is retained only for transport regression and is not an accepted proof. The proof byte budget, per-call provider configuration, LP
checker and scientific family/model conventions are preserved.

Review uses `builtin:probabilistic-escalation-v1` with probability **1 (100%)**
and `override` mode. Its explicit children are `builtin:single-v1` (primary)
and `builtin:math-decomposed-v1` (Moonshot audit), with no child provider overrides.
Both run for every submission. **The Moonshot decision controls acceptance and
scores**, including when the two judges disagree. Their decisions and comparison
remain in the dossier. An audit failure retains diagnostics but cannot yield a
score from a successful primary. A failed primary also prevents a judgment.

The audit adapts the offline UCLA Moonshot Math Verifier pattern: a planner
partitions all source into at most eight regions, a whole-proof reviewer and
region workers assess the claims, and a defender/adjudicator handles findings.
No external retrieval is added. Both children inherit the same atomic claim and
scientific policies. All dependency-ready calls can run concurrently: the primary,
planner and whole-proof review can start together, and region reviews can start as
soon as the plan is ready. The optional `limits.max_concurrent_calls` is omitted;
there is no fixed four-worker limit. They share 23 requests and 753,664 reserved
completion tokens across at most one experiment follow-up. Identical checker
requests are deduplicated; distinct requests must fit
the existing one-follow-up limit. Each phase retains a 1,200-second wall deadline
and the pipeline retains 2,400 review seconds. These ceilings do not guarantee
completion; see [the budget derivation](migration.md#configured-review-budget).

Every file needs an audit region, so **submit at most eight files**, including
`submission.json`, `bound.json`, proof files and experiment code. This effective
limit applies despite the generic `max_files` ceiling of 32. The baseline uses
three files; the existing candidate uses seven. The planner schema enforces the
same one-to-eight range as runtime validation. Separate failure artifacts
(`llmv-initial-review-failure` or `llmv-final-review-failure`) are diagnostics,
not judgments. Failed-audit diagnostics can retain a successful primary, but the
configured override requires the Moonshot audit to complete before scoring.
Overdeep malformed audit output is retained as bounded diagnostic text with its
attempts charged.

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

The verifier pin is `c3196332cbd0ea75a8f478673e04a78ae03e581d` (version 0.8.2); the optional
structured-claims protocol remains disabled. This source change does not alter
the existing live benchmark or page.

## Reusing accepted judgments

Submission reuse is temporarily disabled by default while it needs further work.
The `main` track declares `submissionAdmission: {"version": 1,
"submissionReferencesEnabled": false}`. After explicit owner opt-in, a solver
can cite up to 16 immutable submission IDs in the original
`submission/submission.json`:

```json
"prior_results": [{"submission_id": "SOURCE_SUBMISSION_ID"}]
```

Explain how the source's original declared claims support the current proof.
Hosted intake resolves references through the platform and freezes the accepted
source bytes; copied dossiers or discovery responses cannot authorize reuse.
Source visibility, acceptance, dependency checks and any required human review
still apply. Prior scores and checker witnesses do not replace current scoring
or required checks, and incidental source prose is not an accepted claim.

For an **already imported benchmark**, changing or reimporting this source does
not update live admission. Its owner must separately PATCH that benchmark with
[the admission body](examples/submission-admission-update.json), which now disables
references. Explicit opt-in remains possible by changing its
`submissionReferencesEnabled` value to `true`. An admission update requires no
reorg and applies only to **new candidates**. Retries and rejudgments
keep their original admission snapshots; an originally disabled candidate must
be submitted anew to use references. Local intake likewise needs explicit enabled
`YUKON_EVALUATION.admission` plus trusted resolution; `benchmark.json` is not an
intake authorization fallback. Reuse is separate from attribution declarations;
this change does not alter DAG settings or enable `structured_claims`.

See the pinned verifier's [admission contract](https://github.com/Layr-Labs/llm-verifier/blob/c3196332cbd0ea75a8f478673e04a78ae03e581d/docs/submission-admission.md)
and [prior-result format](https://github.com/Layr-Labs/llm-verifier/blob/c3196332cbd0ea75a8f478673e04a78ae03e581d/docs/prior-results.md).

<!-- llm-verifier:solver-guidance:begin -->
## Preparing a submission

Read the [solver guide](https://github.com/Layr-Labs/llm-verifier/blob/c3196332cbd0ea75a8f478673e04a78ae03e581d/docs/solving.md) and
[solver skill](https://github.com/Layr-Labs/llm-verifier/blob/c3196332cbd0ea75a8f478673e04a78ae03e581d/skills/llm-verifier-solver/SKILL.md)
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
