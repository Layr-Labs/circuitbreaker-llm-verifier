# Accepted expression interface, versions 1 and 2

`export_metrics.py` rebuilds a trusted shared-verifier dossier and appends
`metrics.circuitbreaker` plus the additive table envelope
`metrics.circuitbreaker_frontier` to its unchanged native projection. The rich
`circuitbreaker` record keeps its original fields and meaning and follows
the schema selected by the validated declaration version:

| Declaration and rich metric version | Rich schema | Declaration schema | Table version / comparator |
| --- | --- | --- | --- |
| `circuitbreaker-expression-v1` | [v1 metric](schemas/accepted-metric-v1.schema.json) | [v1 bound](schemas/bound-v1.schema.json) | `circuitbreaker-frontier-v1` / `circuitbreaker-v1` |
| `circuitbreaker-expression-v2` | [v2 metric](schemas/accepted-metric-v2.schema.json) | [v2 bound](schemas/bound-v2.schema.json) | `circuitbreaker-frontier-v2` / `circuitbreaker-v2` |

Legacy schemas, fixtures and statement semantics remain unchanged. Authentic
legacy dossiers under review v2/scoring v1 remain exportable as v1 only; the
legacy review path rejects every v2 declaration. Current review v3/scoring v2
accepts either declaration version. Mismatched review/policy pairs are rejected. A newly
reviewed v1 declaration still exports the v1 tuple, even under the new scoring
policy. Version 2 can also carry any legacy guarantee kind. `expressions.py` supplies
additional exact rational, domain and implication validation. Both schema shape
and runtime validation matter. The full statement and its two source fields
must be jointly accepted, not just present in a snapshot.

The record includes the full reviewed declaration, both identical joint
statements, exact score records, hashes of all referenced evidence, and the
intake, dossier, selected review, complete submission and verifier identities.
The `accepted_provisional` status denotes the configured atomic LLM assessment;
it does not denote independent human approval. Store the original dossier too.

| Accepted declaration | Native `scores` | Full metric |
| --- | --- | --- |
| Finite pair | `[float(c), float(p)]` | Exact full expression and exact pair |
| Lane endpoint | `null` | Explicit lane or stronger full polynomial/log guarantee, exact p |
| Global endpoint | `null` | Explicit `superpolynomial` theorem |
| Relevant string-only lemma | `null` | Absent |
| Invalid claim | `null` | Absent |
| Simulation or failed exporter | No publishable output | Absent |

For explicit endpoints, score fields containing `positive_infinity` are
transport conventions only. In particular the two such fields for a global
endpoint are not an infinite exponent and coefficient. No JSON numeric Infinity
is emitted. Exact scientific ordering never reads native floats.

## Registered Complexity table

The prepared generic metrics-frontier table selects the `circuitbreaker-v2`
comparator; deployment requires a compatible consumer. [The complete presentation](examples/circuitbreaker-presentation.json)
combines the table with coefficient history grouped by exponent using the merged
Yukon `results.scoreTracks` contract. It selects `gate_exponent`, numeric descending
order and the highest eligible finite group by default; `results.defaultView`
selects coefficient over time. Each group retains its complete eligible history,
including results dominated on the global Pareto frontier. These are display
groups, not new benchmark tracks or different scoring rules. Group keys use native
binary64 exponent values; exact rational distinctions remain in full expressions.
Null-score endpoints never create numeric groups. Eventual bounds use their
leading term for the same finite or endpoint projection rules; that projection
does not discard the remainder from either accepted metric.

[The table-only fragment](examples/complexity-result-tab.json) remains available
for existing integrations. Both templates add `resultTabs` at the chosen page or
track presentation level. The tab selects
`circuitbreaker_frontier`, version `circuitbreaker-frontier-v2`, scope
`uniform-np-full-b2-actual-length-all-lengths-v1`, model `full-b2-v1`, and notation
`lower-bound`. The single `circuit_size` column is maximized and initially
sorted descending. The fragment requires that registered comparator extension;
it does not itself publish source, update a page, or register a benchmark.

The v2 consumer accepts both authenticated v1 and v2 envelopes from that same
namespace, validating each exact version/comparator/scope/model/notation/resource
tuple. V1 envelopes must contain only v1 ASTs. The v1 consumer/config remains
unchanged and rejects v2; changing just a label never makes an envelope valid.

The [v1](schemas/frontier-v1.schema.json) and [v2 table schema](schemas/frontier-v2.schema.json)
describe the exact version-specific envelope:
`version`, `comparator`, `scope`, `model`, `notation`, `intake_sha256`,
`dossier_sha256`, and `resources`. Its sole resource `circuit_size` contains:

- `statement`: the complete unchanged accepted joint claim statement.
- `bound`: the exact canonical `bound.guarantee` AST, including theorem kind.
- `source_sha256`: the first accepted score-source identity, which binds the
  coefficient claim asserting the ENTIRE joint theorem and source snapshot.
- `evidence`: the same validated common proof-file hashes as the rich metric.

The table wire has exactly one resource, `circuit_size`; it does not inherit
the DCP comparator's open resource set. Statement length stays at 8,000
characters. Evidence includes `bound.json` and `proof.md`, with at most 32
entries and the pinned verifier's safe relative path limits: 512 characters,
12 components, 128 ASCII characters per component, and no `.`, `..` or `.git`
components. Files named `constructor`, `prototype` and `__proto__` remain
valid evidence basenames and are treated as inert own data properties; this
does not allow those names as metric/configuration keys. These consistency checks do not authenticate provenance or prove
the mathematics.

The second score record and joint claim are still verified by the existing
acceptance gate before either metric is emitted. The table namespace adds no
independent acceptance path or candidate-supplied comparison hint. The original
rich namespace retains both source records, full declaration/projection and
selected review identity for audit.

Read n as actual input bit length; all logarithms are base 2. A
`polynomial_log` bound represents a liminf guarantee with its exact rational
coefficient, not equality to actual circuit complexity. A lane endpoint is
omega(n^p): its ratio to n^p diverges without asserting a particular logarithmic
rate. Version 2's `eventual_lower_bound` means there exists N such that
C(n)>=L(n)-sum(terms(n)) for every integer n>=N for an explicit loss; for a
Big-O loss there exist K>0,N such that C(n)>=L(n)-K*n^a*(log2 n)^b*(log2 log2 n)^d
for every integer n>=N. A liminf assertion alone never certifies either loss.
`superpolynomial` means all fixed powers for the SAME family and all
sufficiently large lengths, not exponential growth or n^infinity. The numeric
graph remains the declared [coefficient, exponent] projection. Endpoints retain
null native scores while participating in the Complexity table.

Only current, access-authorized, valid results from one compatible track enter
the table. Unreviewed, revoked and stale results are filtered by the platform's
current-result policy; malformed or incompatible envelopes are omitted. These
are separate from mathematical comparison. Unsupported bounds compare as
unknown, never as equivalent. Equal full guarantees retain ties even if their
weaker graph projections differ.

## Comparison without a UI dependency

After authenticating and validating accepted metric records through the trusted
organizer channel, use:

```python
from expressions import compare

relation = compare(
    left["bound"]["guarantee"], right["bound"]["guarantee"],
    left_convention=left["convention"],
    right_convention=right["convention"],
)
```

This is a pure implication comparator. It does not authenticate records, load
candidate code, consult model hints, compare proof quality or establish new
mathematics. Callers must require the supported metric `version` and current
accepted status; unsupported versions, missing data, revoked acceptance or
incompatible conventions cannot enter a frontier. `normalize` validates and
copies the canonical AST; it does not round or repair it.

The comparator applies sufficient implication rules to real-valued guaranteed
envelopes. `stronger` means those rules establish the left-to-right direction
but not the reverse. It does not prove semantic inequivalence for integer circuit
sizes: 4n-1/4 and 4n-3/4 both imply the integer bound C(n)>=4n. There is no integer
rounding equivalence policy. `equivalent` means both directions follow under
the supported rules, even when families or proofs differ. Keep tied rows. `unknown` never means
equality, dominance, or permission to suppress a row. The native coordinatewise
Pareto frontier and this expression ordering are intentionally different.
Compare complete current accepted history, including null-score endpoints,
not just one page, native-frontier rows, or non-null scores.

The expression grammar and exact comparison rules are in [rules.md](rules.md).
All logarithms are base 2. Hypothetical n²/log2(n) is below a lane endpoint at 2
and above every fixed lane below 2. Its greatest qualifying exponent is not
attained. Changing only its weaker displayed projection creates no new
scientific progress. A bare lane theorem does not invent a log rate; a global
theorem does not invent an exponential rate. The p=1 negative-leading-log
case has no legal main projection and belongs in an unscored lemma.

## Trusted export and generation

The production metrics job consumes only the current workflow's
`llmv-package`, `llmv-intake` and `verifier-dossier` artifacts. It checks the
separately supplied intake identity, reconstructs all three or five phases at
the recorded finalization timestamp, and checks config, atomic acceptance,
source records, statement/evidence identity and algebraic projection agreement.
The benchmark selects `circuitbreaker-score`, not the shared raw
`benchmark-score`. Candidate uploads must never be treated as these artifacts;
hashes alone cannot prove reviewer authenticity.

For a local organizer pipeline, `finalize_score.py` delegates to the shared
finalizer and then exports metrics. It clears the final score on failure. The
ordinary `llm-verifier finalize` command alone intentionally produces only the
shared raw projection. To export an already finalized trusted dossier:

```sh
python export_metrics.py --intake intake.json --dossier dossier.json --out score.json
```

Regenerate hosted files with `python tools/generate_hosted.py`. The script runs
the pinned upstream generator and applies only the checked-in metric and test
job overlays plus the artifact/CLI selection in `benchmark.json`. It also checks
that the separate complete presentation still binds the generated objective
names/order/directions, coefficient chart, exponent grouping and table contract.
The template is a checked-in companion; the generator does not import a page or
insert presentation fields into the benchmark manifest. Run with `--check` for a
read-only drift check. Running upstream `github` alone removes
the overlays; CI detects that drift. No unsupported platform comparator setting
or native endpoint type is configured.

PR/push checker CI uses the same-run prepared source archive and verifier wheel
to run the custom declaration check, regression tests and generator drift check
in a job with no supplied credentials. Benchmark execution does not run test
fixtures or synthesize acceptance. Production synthetic-test examples are not
accepted results: tests build explicitly labelled unsigned records solely to
exercise exporter validation. No provider request or credential is used.


## Cross-language conformance corpus

`tests/fixtures/circuitbreaker-frontier-v1.json` contains 58 eligible synthetic
exporter outputs, 3,364 independently assigned pairwise expectations, selected
frontier/tie scenarios and labelled omitted cases. Every eligible envelope was
produced by the v1 exporter after a deterministic synthetic phase
chain at its original revision; there are no copied live results or relabelled old-runtime dossiers.
Expected strength tiers are explicitly arranged from weaker to stronger without
calling `expressions.compare`; tests compare the authoritative Python algebra
against those expectations and the platform uses the same corpus.

Only the table namespace is retained in the corpus to avoid duplicating rich
records. Separate regressions freeze the legacy statement/semantic bytes and validate
the current exporter's exact joint binding. Current review/rules changes necessarily
change fresh intake and dossier hashes; those do not rewrite historical fixtures. `mutated-exporter` cases deliberately
corrupt a copy for consumer validation; `consumer-state` cases retain valid
metric bytes while requiring the current-result policy to omit them. A valid
platform result awaiting manual promotion is distinct from an unreviewed claim.
Simulation is rejected by the exporter and has no output.

`python tools/generate_frontier_fixtures.py --check` verifies the frozen v1 file's
SHA-256, not reconstruction by the current exporter. In-place v1 regeneration
is disabled; historical reconstruction requires its original source revision.
Current exporter replay belongs to the v2 corpus, generated by
`python tools/generate_remainder_fixtures.py`; its `--check` mode reconstructs
the exact current exporter output and compares bytes. CI retains the 3,364 legacy
semantic comparisons as well as v2 schema/export tests. The current pin is
`c3196332cbd0ea75a8f478673e04a78ae03e581d`; its optional `structured_claims`
protocol remains disabled. Existing joint statements and file binding are
unchanged. Neither these tests nor this table integration requalifies the
baseline proof mathematically.
