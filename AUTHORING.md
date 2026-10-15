# Generic llm-verifier challenge

This starter adapts the approved-definition/edit-map process from
Layr-Labs/yukon-template@72a6b96b4b5385f199a5a29ecdd742aeb45ad55f.
It uses the generic llm-verifier phase protocol. Stock AutoYukon init, Moonshot,
accept/revision/reject and its sealed runtime are not being invoked.

The controller keeps the approved definition and generated receipt outside this
export. Use `llm-verifier context --root . --task build --definition PATH
--receipt PATH`, then edit only the returned paths. `llm-verifier check --root .
--definition PATH --receipt PATH` validates protected bytes, the bounded config,
candidate snapshots and scripted public cases. Managed context/check reject all
unlisted files/directories and symlinks, skipping only root .git metadata without
traversal. Keep receipts, environments, caches and evaluation output external
when validating an authoring export. It executes no Python candidate
or checker and makes no model calls. Scripts only test orchestration.

`python3 tools/public_conformance.py --root .` separately executes declared
public checker cases with Docker and no model credentials. Prepare the image
explicitly first with `docker pull IMAGE_FROM_VERIFIER_CONFIG`. Execution never
pulls implicitly. The public-conformance workflow does that setup and needs no
cross-repository access or model secrets. These known-answer tests are not an
independent mathematical assessment and do not emit official scores.

Standalone controllers synchronize approved edits back into their external
definition, then initialize a fresh directory/receipt to finalize the inventory.
Managed builders do not change rules, workflows, configuration, instruction
files, scaffold inventory or the external receipt. Baseline/example content and
independently labeled public cases remain visible and are not private tests.

## Live integration

Install llm-verifier at the reviewed full commit pinned in yukon-scaffold.json.
Use separate intake, prepare-runtime, experiments, review-initial, optional
followup/review-final, and finalize commands from the verifier README. Only the
review subprocess/jobs receive OPENROUTER_API_KEY. There is no combined local
execute-with-secrets command. `yukon run` does not orchestrate these phases.

The benchmark workflow uses a pinned reusable workflow and a matching source
revision. Before hosted use, the verifier repository must permit reusable
workflow access for this caller. A private/internal verifier also needs read-only
checkout access via VERIFIER_CHECKOUT_TOKEN when the caller token is insufficient.
Configure OPENROUTER_API_KEY only for authorized review. A generated workflow is
not evidence that access, hosted execution, import or qualification has occurred.
No settings are changed by this starter.

Finite scores use the raw coefficient in binary64 transport; exact reduced
rational values remain in checker evidence and dossier score_records. Precision
ties remain a platform limitation. An explicitly accepted terminal record is
retained as positive_infinity with projection status unsupported, and its Yukon
scores are null. It cannot enter native Yukon frontier ordering; manual review
and promotion must inspect the typed dossier. Ordinary unscored lemmas have no
terminal record. No normalization, sentinel or numeric Infinity is used.
