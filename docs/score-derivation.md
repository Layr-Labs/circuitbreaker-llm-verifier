# Retained-validity score derivation

This adapter uses Yukon's ordinary full-history score-only reorg. Yukon retains the selected validity result and human decisions, issues the exact prior-result descriptor, and rejects a derived score if those selected identities have changed. The verifier judges whether a new structured metric faithfully represents the already accepted full claims. It does not re-prove them or emit a new validity decision.

The default `verifier.json`, scoring policy v2 and hosted workflow remain unchanged. Derivation uses an explicitly materialized source profile with scoring policy `coefficient-exponent-expression-v3` and capability `yukon-retained-result-v1`. A complete protected evidence plan and a reviewed immutable verifier commit are required before that source can be used.

## Preparation and routing

1. Capture the complete current Yukon submission snapshot with selected score and validity evidence. Stage an organizer-controlled artifact index and exact original GitHub artifact ZIPs. The index format is `{"version":"circuitbreaker-derivation-artifacts-v1","records":[...]}`. Each record contains `submission_id`, `source_files_sha256`, and `artifacts`. The latter maps `intake` and `accepted_results` to exact repository/run/attempt/artifact IDs and archive/member SHA-256 values; it is null for retained-invalid results. The generic verifier documentation defines each artifact reference. Local archive filenames are `<artifact_id>.zip`.
2. With the new verifier installed, run `python tools/derivation_profile.py --snapshot SNAPSHOT --artifact-index INDEX --archive-root ARCHIVES --verifier-revision COMMIT --out NEW_DIRECTORY`. It validates every retained record, rejects missing or extra index entries, and creates a separate source directory. No source is produced on an incomplete plan. Original archives are not copied into that source or executed. Inspect `reorg/score-derivation-report.json` and the generated source before publishing it.
3. The generated workflow first checks the protected configuration and complete plan hashes. Baselines and candidate jobs explicitly requesting validity use the ordinary verifier with unchanged validity, review standard, judging profile and budgets. Only score-only jobs carrying Yukon's retained-result descriptor use `derive-score.yml`. Missing descriptors, stale selections, missing artifacts and failed derivations never fall back to full review.
4. Use the normal authorized Yukon full-history reorg after source review, a representative nonpublishing pilot, artifact-availability checks and budget review. Materialization and unit tests do not perform this step. Do not substitute a future-only activation when the scoring policy changes.

The captured 100-row inventory motivated this design, but neither generator nor adapter hardcodes three records or those 100 IDs. For that snapshot the platform would request 80 retained score-only candidates, 18 candidates needing both fields, and one baseline; two cancelled rows were excluded. These are execution counts, not a prediction of 99 proof reviews or an approved model budget. The 18 missing-validity jobs and baseline remain separate ordinary work.

## Result behavior

- A rich valid claim can yield a v2 eventual lower bound at exactly the same leading expression, family, theorem text, convention and numeric projection. The derivation reviewer must cite the literal original full joint claim and judge the remainder transcription. Structural checks alone do not establish it.
- A remainder cannot be invented from a leading coefficient, reviewer prose, a liminf statement, or unspecified little-o loss. Already represented eventual bounds cannot be tightened by this adapter.
- An unsupported transcription retains the original declared bound and its actual v1/v2 version. It does not promise preservation of a different metric created by a previous derivation.
- Retained-invalid records and authenticated unscored supporting lemmas produce explicit null scores without a provider call. A main claim cannot be silently discarded as an unscored lemma.
- Original `metrics.circuitbreaker` values and full joint statements retain their original provenance. `metrics.circuitbreaker_frontier` contains the selected representation with new scoring provenance; `metrics.circuitbreaker_derivation` records the separate derivation and origin. Top-level `derived_from` echoes Yukon's exact descriptor. There is no `valid` field.

The shared workflow uses authenticated empty accepted-result resolution and emits only `accepted-results-receipt.json` for the new run. It creates no accepted bundle or replacement acceptance anchor. Original accepted claims remain associated with the retained validity job.

The final export installs the verified current-run wheel, verifies exact intake/result artifact descriptors and transfer receipts, requires the current candidate/repository/verifier/run identities, and rejects simulated results. The original artifacts are inert data. Ordinary exporter and baseline-recovery checks remain intact.

## Offline validation and cost

`python -m unittest tests.test_score_derivation -v` covers domain transitions and routing with the ordinary pinned verifier. `tests.test_derivation_integration` requires the new verifier and covers the complete sealed derivation, profile generation, full-plan completeness, null-score branch, and hosted export. The dedicated feature CI job requires the new module; ordinary legacy CI reports these feature integration tests as skips while retaining historical byte/hash checks under the original verifier.

For a nonpublishing pilot, use the generic `simulate-derivation` or a separately budgeted real derivation and `python export_score_derivation.py --derivation RESULT --intake INTAKE --preview --out PREVIEW`. Preview output has `publishable:false` and no top-level `scores`. Synthetic decisions exercise plumbing, not transcription quality or theorem correctness.

The proposed derivation configuration permits one provider request and at most 8,192 completion tokens per applicable claim. This is a configured ceiling, not measured spend or demonstrated extraction accuracy. Baseline and missing-validity work retain their ordinary budgets. No paid calls or production reorg are performed by these preparation tools or tests.
