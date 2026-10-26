# Accepted expression interface, version 1

`export_metrics.py` rebuilds a trusted shared-verifier dossier and appends
`metrics.circuitbreaker` to its unchanged native projection. This record follows
[the metric schema](schemas/accepted-metric-v1.schema.json); its `bound` follows
[the declaration schema](schemas/bound-v1.schema.json). `expressions.py` supplies
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

`stronger` means the left guaranteed lower bound strictly implies the right
within this supported language. `equivalent` means equal guarantee strength,
even when the families or proofs differ. Keep tied rows. `unknown` never means
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
job overlays plus the artifact/CLI selection in `benchmark.json`. Run with
`--check` for a read-only drift check. Running upstream `github` alone removes
the overlays; CI detects that drift. No unsupported platform comparator setting
or native endpoint type is configured.

PR/push checker CI uses the same-run prepared source archive and verifier wheel
to run the custom declaration check, regression tests and generator drift check
in a job with no supplied credentials. Benchmark execution does not run test
fixtures or synthesize acceptance. Production synthetic-test examples are not
accepted results: tests build explicitly labelled unsigned records solely to
exercise exporter validation. No provider request or credential is used.
