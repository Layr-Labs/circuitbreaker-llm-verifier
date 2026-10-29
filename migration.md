# Coefficient/exponent registration and legacy retirement

This repository declares only the current two-objective `main` track. The legacy
single-objective registration is not part of its manifest or current page. Retire
that legacy registration through the supported close and soft-archive lifecycle,
retaining historical submissions, judgments and policy versions. Remove its page
bindings and selectors; do not recreate a legacy UI track during page preparation.

Do not reimport or reorg the legacy registration with this manifest: a normal
reorg cannot change its objectives. These source instructions do not execute
registration, archival or page changes.

For a new registration, independently review and use the frozen source containing:

- Objectives `gate_coefficient` and `gate_exponent`, both maximized, in that order.
- The new scoring/review versions, current fixed-family validity policy,
  atomic acceptance, and the decomposed review profile/budget below.
- The pinned shared verifier and regenerated hosted caller with the bounded
  metrics overlay. The runner must select `circuitbreaker-score`.
- The published (31/10,1) seed declaration and its unchanged proof evidence,
  reviewed under the current supplied-source policy. Passing structural tests
  or accepting an older baseline does not guarantee scientific acceptance.

Keep the closed, soft-archived legacy history outside the current page and
repository track list. Do not rewrite old single-coordinate results into accepted
two-coordinate results without review of the exact joint theorem and source.
Results tied under full-expression comparison retain both proofs; a weaker
projection of an existing expression is not fresh progress.

## Prepared page and source inputs

[The complete presentation template](examples/circuitbreaker-presentation.json)
is copied from Yukon commit `f340c97dad9f170941e03c233a8a263ae76b8998` (merged
PR1283). It combines `results.scoreTracks` grouped by `gate_exponent` with the
`circuitbreaker-v1` Complexity tab. The default numeric view is coefficient over
time at the highest eligible finite exponent. Accepted null-score endpoints
remain in the full-expression table. The older [table-only fragment](examples/complexity-result-tab.json)
and its conformance corpus are unchanged.

The template is a presentation, not an import request: it has no benchmark UUID
or publication flag and uses `unbound/circuitbreaker/main`. `benchmark.json`
contains no page configuration; Yukon stores it separately in
`challenge_pages.presentation`. `tools/generate_hosted.py` validates this companion
against the generated objectives and existing table contract, so regeneration
cannot silently leave it bound to a different score coordinate or metric namespace.
Full platform schema validation and replacement-source binding still apply.

Prepare activation with these explicit inputs:

- The target API origin and compatible API/Challenges UI releases containing
  both registered `circuitbreaker-v1` and `results.scoreTracks`. A merged PR alone
  does not establish that the selected deployment supports either field.
- The approved published challenge source SHA, repository
  `Layr-Labs/circuitbreaker-llm-verifier`, root directory `null`, `main` track,
  exact two-objective manifest and policy versions above. Keep verifier pin
  `5c3e66acffc95aa5a9cfc31183dbefc7ed24b9bf` and artifact `circuitbreaker-score`.
- A separately selected replacement registration name/source branch, actual
  replacement UUID, access/publication choices and any existing page revision
  used by the authorized rebind. Legacy UUID
  `015f2bee-4f6d-4db6-8e51-8616071c0454` belongs only to closed, soft-archived history;
  exclude it from current page bindings and track selectors. Preserve the current
  two-objective registration when resolving source-branch or name conflicts.
- The replacement benchmark JSON read from that API and the exact reviewed
  source SHA. Retain its baseline job, workflow run, selected evaluation and
  dossier identities for audit. Do not fabricate a replacement UUID or relabel
  the old single-coordinate results.

Use a real source branch pointing at the reviewed commit, then verify the
registration's resolved `sourceRef`. A bare 40-hex commit in `sourceBranch` passes
the current name validator and can resolve during import; baseline dispatch also
uses a managed branch at that commit. It is nevertheless unsuitable for a live
registration: candidate validation calls `getBranchHead`, resolves
`heads/<sourceBranch>`, and creates a PR against that branch even for pinned Pareto
evaluations. Without an actual branch, submissions fail. The recorded evaluation
SHA provides the immutable pin; a bare SHA is not a replacement for the branch.

Use the merged Yukon's existing pure preparation helper, from its checkout:

```sh
bun scripts/prepare-circuitbreaker-page.ts \
  --benchmark-json /path/to/replacement-benchmark.json \
  --benchmark-id "$REPLACEMENT_BENCHMARK_ID" \
  --expected-source-ref "$REVIEWED_SOURCE_SHA" \
  --api-url "$TARGET_API_ORIGIN" \
  --out /path/to/circuitbreaker-unpublished.json
```

The helper checks UUID/source equality, repository/root, ordered objectives and
directions, manual promotion and policy versions. It binds the clone, source
check, submission commands and rules link, validates the full import, and writes
`published: false` without network calls or overwriting an output. Its embedded
presentation must match the reviewed challenge template before use. It rejects
the legacy UUID. Import/rebind, page publication and opening submissions remain
separate authorized writes; the preparation command performs none of them.

## Registration and qualification boundary

There is **no supported model-free create-paused registration** in the inspected
Yukon contract. `POST /api/benchmarks` resolves the source and queues baseline
validation. For this schema-v2 manifest, `createChallenge` creates the track in
`validating` state and inserts its baseline job atomically. The GitHub Actions
runner creates/uses a baseline source branch and dispatches `benchmark.yml`, which
can invoke the configured model. `hidden` controls listing, not execution. Neither
`paused` nor a no-evaluation flag is accepted by the registration contract.

The lifecycle is `validating -> ready -> open`; pausing is `open -> paused`, not
a way to suppress the import job. Leaving a successful registration in `ready`
keeps submissions closed, but its baseline evaluation has already run. Offline
source/template preparation is the no-model alternative until registration and
its review spend are authorized. These details follow the pinned platform
`src/api/routes.ts`, `src/store/postgres.ts`, `src/platform/state.ts` and
`src/execution/github-actions.ts`; check the deployed contract before activation.

Technical readiness is distinct from scientific qualification. A successful
Pareto baseline workflow can make the evaluator ready with null scores; `ready`
does not certify that the baseline theorem was accepted. Before claiming an
accepted baseline, inspect its current trusted dossier, complete joint claims,
evidence binding and both exported namespaces. Before scientific promotion,
retain independent human review under the manual-promotion policy. The unchanged
(31/10,1) proof is a literature seed and has not been newly qualified by these
source, schema, synthetic or CI tests.

## Configured review budget

`verifier.json` selects `builtin:math-decomposed-v1` at verifier version 0.7.0.
All members inherit the unchanged `gpt-6-astra`, high reasoning, JSON-schema
responses and background mode, with one attempt, at most 32,768 completion
tokens and a 1,200-second request timeout. The profile permits at most eight
source regions. A planner, whole-proof reviewer, one worker call per region,
conditional defender and adjudicator share the pipeline budget.

The shared limits are **21 provider requests, 688,128 reserved completion tokens,
and 2,400 review seconds**. The largest follow-up path uses ten initial calls
(planner, whole-proof reviewer, eight region workers) and eleven final calls
(whole-proof reviewer, eight region workers, defender, adjudicator). The initial
round returns its experiment request before defense/adjudication, and the final
round reuses the frozen plan. A round without a follow-up uses at most twelve
calls. Each round also has a hard 1,200-second wall deadline; these limits do not
guarantee every permitted call will finish.

There is at most one experiment follow-up; initial and follow-up computation
each allow 20 seconds, with 40 computation seconds total. Generic source limits
remain 32 files, 65,536 bytes per file, 262,144 submission bytes and 131,072 prompt
bytes per call. The selected profile requires at least one region per submitted
file, so its **effective file limit is eight**, including `submission.json` and
`bound.json`; excess files fail preflight before model calls.

Those are request/token/time limits, **not a dollar-spend ceiling**. Completion
tokens are reserved at the full per-call cap even for short responses; input
billing is additional. Registration retries or further evaluations start
additional pipelines. This source repin does not activate a registration, run a
model, change credentials, or establish scientific qualification. Approve the
exact activation/evaluation scope separately.
