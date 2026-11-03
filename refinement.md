# Historical representation refinement

`export_refinements.py` prepares a finite owner-page configuration overlay over
existing selected results. It does not change a score, claim, review, intake,
dossier, accepted-results bundle, workflow or scientific decision. Ordinary v1
and v2 exports, credit and accepted-result reuse remain unchanged.

This path handles authentic provisional v1 results whose **entire accepted
joint declaration** already states an explicit remainder. Both original atomic
claims must cover that exact declaration. The organizer supplies and checks the
full v2 transcription; the tool cannot infer a remainder from a liminf statement
or prove that natural-language theorem text implies a chosen AST. Structural
validation alone does not authorize publication.

## Input and execution

Use [the synthetic plan](examples/refinement-plan.json) as a shape example. A
plan has version `circuitbreaker-refinement-plan-v1` and 1–32 records, each with:

- `id`: a distinct 1–160 character ASCII identifier.
- `bundle`: the authentic extracted original artifact directory described below.
- `downloadReceiptSha256`: the raw-byte SHA256 of its `download-receipt.json`,
  pinned through the organizer's trusted artifact acquisition channel.
- `platform` and `platformSha256`: a captured selected public submission JSON
  (`{"submission": ...}`) and its independently pinned raw-byte SHA256.
- `originalSource`: the frozen original producer checkout, with the original
  submitted files restored under `submission/` from the authenticated intake.
- `originalPython`: an already prepared trusted interpreter environment containing
  that run's exact original wheel. Retain its virtualenv path; do not resolve the
  executable symlink to a system Python.
- `targetBound`: the full organizer-supplied v2 declaration.

Relative paths resolve against the plan directory. The supplied receipt and
snapshot pins are trust inputs, not proof of authentic origin by themselves.
Never accept a candidate-created receipt/dossier or let candidate content select
an interpreter or source checkout. The exporter downloads and installs nothing.

```sh
python3 export_refinements.py --manifest trusted-plan.json --out new-mapping-output
```

The output directory must not exist and must be outside the original source and
artifact directories. All records must pass before final output is created.
Previous output and all input artifacts are preserved. Outputs are:

- `manifest.json`: exact input plan bytes.
- `receipt.json`: source/content identity, target declarations, original wheel,
  archive and source digests, selected review and full claim bindings, unchanged
  status, and hashes of four byte-identical reconstructed outputs.
- `representation-refinements.json`: the prepared value of a page tab's optional
  `representationRefinements` field, following
  [owner-refinements-v1.schema.json](schemas/owner-refinements-v1.schema.json).

`audit.manifestSha256` and `audit.receiptSha256` hash the exact bytes of the first
two outputs. The receipt precedes record assembly and never contains final
record hashes. It binds exporter source by a file manifest and its content digest;
it does not assert an unpublished Git commit identity. Temporary full replays are
removed after use; submitted proofs are not copied to final outputs.

## Original replay boundary

The bounded implementation supports original producer
`2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851`, verifier
`c3196332cbd0ea75a8f478673e04a78ae03e581d` (0.8.2), and its three-phase,
first-attempt successful GitHub dispatch artifact layout. Other histories fail
closed. `refinement/original-source-v1.json` pins all 49 original protected files.
Each of the nine original ZIPs must match its GitHub metadata digest; every
extracted member must match the archive. Required directories are `llmv-package`,
`llmv-intake`, `llmv-initial-experiments`, `llmv-initial-review`, `verifier-dossier`,
`circuitbreaker-score`, `benchmark-score`, `verifier-accepted-results`, and
`circuitbreaker-export-lineage`. `run.json`, `artifacts.json` and
`download-receipt.json` bind their run and artifact identities.

Replay runs in a separate process under the original interpreter, with `-I -B`
and a fresh empty `pycache_prefix`. Every installed verifier package file must
match the original wheel. The frozen original exporter and algebra execute from
verified source bytes. Credentials are not forwarded. Network, review and
experiment execution boundaries raise errors; there is no fallback to new review.
Original config/rules, complete atomic claims, submitted files, selected review,
lineage and platform evidence are checked. Dossier, rich score, native score and
accepted-results output must all reproduce their original bytes.

The platform's `candidateKey` and the verifier's submitted-files digest are
different identities. The former is pinned verbatim from the selected platform
packet; the latter remains bound by the original intake and audit receipt.
`origin.sourceSha256` is the original first atomic claim's score-source digest,
not a Git SHA. Both score and validity evidence objects are retained in full.

The separate current v2 parser requires identical family, theorem, convention,
finite projection and leading guarantee. Only the expression version and the
explicit eventual remainder representation may differ. The current bounded
exporter requires original `review`, `human: null`, and `accepted_provisional`
status; it does not turn these into human acceptance.

## Owner selection and hashing

Only an explicitly configured `circuitbreaker-v2` metrics-frontier tab consumes
`representationRefinements`, with version `circuitbreaker-owner-refinements-v1`.
Owner publication of the page revision is the mapping authority. There is no
`approved` flag or fabricated signature. Page assembly, compare-and-swap and
publication happen separately; this command makes no live writes. The owner may
set `useForFrontier: true` for consistent table/history use.

The consumer checks the exact current benchmark/submission, candidate, packet,
both selected evidence objects, and raw rich/frontier namespace digests. It also
checks current public validity, access/revocation and human rejection gates.
Changed identities or invalid mappings retain the original valid metric with a
diagnostic; no stale target is applied. Future changes to a benchmark's active
source alone do not invalidate an unchanged selected packet. The original full
joint statement remains visible alongside the derived representation and audit
identity. This mechanism grants no accepted-result reuse or credit.

Shared digest bytes are Python
`json.dumps(value, sort_keys=True, ensure_ascii=True, separators=(",", ":"), allow_nan=False)`.
Object keys sort by Unicode code point. Escape code units **at or above U+007F**
(including DEL), using lowercase `\uXXXX`; astral characters use surrogate pairs.
Use standard JSON escapes for controls, quote and backslash. Reject lone
surrogates, all JSON numbers and non-JSON values. Only strings, booleans, null,
arrays and objects are admitted to this hash contract; rational values use
canonical decimal strings. `recordSha256` hashes the complete record except its
own top-level field. Original namespaces are hashed whole, without normalization.
The synthetic canonical vectors exercise DEL, Unicode key order, controls,
astral text, booleans and null.

Ordering compares real-valued guaranteed envelopes by sufficient implication.
It does not establish semantic inequivalence for integer circuit sizes and has
no integer-rounding equivalence policy. This representation change requires no
new benchmark registration or scoring execution.
