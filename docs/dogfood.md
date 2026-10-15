# CircuitBreaker pilot validation

As of 2026-10-15, the challenge has passed offline validation and protected-checker
execution. No live model review has run. This is a bounded pilot, not independent
mathematical qualification or an unrestricted research launch.

The [private challenge repository](https://github.com/Layr-Labs/circuitbreaker-llm-verifier)
published its scientific source as signed commit
[`5cd3fba88ec6a1ae1362cfb2a3640ce3be83e7dd`](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/commit/5cd3fba88ec6a1ae1362cfb2a3640ce3be83e7dd),
tree `269f45ca7508e9565ddae32ea5bdcb31b6b39212`, with sole parent
`17c8c3560fed414017c38543a060b3c7e63a49cd`. GitHub reports the signature as valid;
independent local object checks matched the parent and exact staged tree. This
publication is code delivery, not official mathematical promotion.

The package uses llm-verifier v0.3.0 from
[upstream PR 3](https://github.com/Layr-Labs/llm-verifier/pull/3), pinned to revision
`e6952b76b1057b026f408f88c3e8c55b419f7344`. The reviewed generator wheel has SHA-256
`3fef345c66646812c882bd03070f21d924236ce97cb1d75ceb2474f3cc0bee79`.
The scientific evaluation input digest is
`3366b8c960a0788fea1b9b678cbfa5ca8872b3bc2b01b1e8d05e9a37e336fd8c`.
Documentation-only additions do not change those evaluation inputs.

| Check | Result | Scope |
| --- | --- | --- |
| Local unit tests | 19 passed | Canonical identity, rational encoding, malformed input, exact ordering and LP arithmetic |
| Generated scripted cases | 16 passed | Simulated protocol decisions; no model calls |
| Protected Docker cases | 11 passed | Actual organizer-owned main-claim and LP checkers; no model credentials |
| Complete inventory check | Passed | All declared files present, protected bytes unchanged, undeclared files rejected by generic inventory contract |
| Dossier distinctions | Passed on marked simulations | Terminal tag and unsupported projection differ from unscored lemmas, rejected endpoints and huge finite coefficients |
| Hosted public checker CI | Passed | Credential-free workflow on the published scientific commit |
| Production phase preflight | Six cases completed | Real intake, runtime preparation and experiments; five actual Docker main-claim executions, no simulation |
| Live model pilot | Not run | Awaiting explicit authorization to send private evaluation materials to OpenRouter |


The [hosted public-checker run](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/actions/runs/38790503666)
completed successfully. The separate production preflight used the immutable reviewed
wheel and executed `intake -> prepare-runtime -> experiments` for each of the six
frozen cases. These phases made no model calls. Review-material generation also
confirmed that every complete prompt fits the protected 131072-byte limit:

| Frozen case | Actual checker executions | Prompt bytes | Expected mathematical review only |
| --- | --- | --- | --- |
| baseline | 1 | 68768 | Accept exact coefficient 31/10 |
| or-coefficient-one | 1 | 41354 | Accept coefficient 1 |
| auxiliary-lemma | 0 | 28300 | Valid, unscored lemma |
| remainder-error | 1 | 66360 | Reject declared 16/5 claim; unscored |
| restricted-model | 1 | 35020 | Reject full-B2 overclaim; unscored |
| prompt-injection | 1 | 68049 | Reject unsupported overclaim; ignore injected instructions |

The last column is the frozen expected answer, **not an LLM verdict**. No model has
accepted the baseline or any other claim in this pilot. Preflight success establishes
execution and prompt-size compatibility; it does not assess proof correctness. Private
intakes, execution records and prompts remain outside the repository; this report
contains aggregate results only.

The protected cases include the exact 31/10 seed, coefficient one, a clearly marked
synthetic endpoint, a huge finite rational, nonreduced and malformed rationals,
zero denominator, oversized finite encoding, and the unchanged LP's exact primal/dual
certificate. A long-rational LP boundary checks that exact arithmetic output is not
silently rounded or truncated.

The broader public corpus covers a linear remainder incorrectly called o(n), a
monotone-to-full-B2 model swap, missing uniform NP construction, sparse lengths,
changing the family with the desired constant, fake infinity, prompt injection,
forged checker evidence, auxiliary-only acceptance and the wrong main-claim ID.
These are public development examples, not hidden qualification cases.

Six live cases are planned: the published 31/10 baseline, a self-contained OR-family
coefficient-one proof, an unscored auxiliary LP lemma, the incorrect 16/5 remainder
claim, a restricted-model overclaim and an injected/forged-review overclaim. Their
expected outcomes are fixed before the pilot. No results or model competence are
claimed here, and no score-shopping reruns are authorized by this report.

The Li–Yang, Ben-Sasson–Kopparty and Shoup primary PDF downloads independently match
the historical baseline source hashes. This verifies source identity, not all proof
steps. Historical organizer approval is retained in baseline provenance, separately
from the candidate proof. The seed is a literature-based application; no new circuit
lower bound is claimed by this work.

Finite coefficients are exact reduced rationals, with a 256-digit resource bound per
integer. Complete canonical checker inputs are limited to 16 KiB, and individual
strings to 4096 characters. Oversized proof bundles are unsupported and never
truncated or moved out of the identity binding. Mathematical infinity is supported
as an accepted tagged terminal record; the current finite-only platform emits
scores:null with `yukon_positive_infinity_unsupported`. Every endpoint success fixture
is explicitly synthetic and supplies no actual superlinear proof.

Official frontier or endpoint promotion requires human review of the exact theorem,
proof and evidence. Offline passes and LLM agreement alone do not constitute proof
certification or permission to promote a record.

## Hosted review access remains unverified

Hosted mathematical review has not been exercised. The generated benchmark workflow
calls the pinned reusable verifier workflow and requires that the verifier repository
permit that caller. It also needs read-only verifier checkout access when the caller's
repository token is insufficient, and separately authorized model credentials.
Successful public-checker CI does not test these access requirements or OpenRouter
review. No settings or secrets were changed to bypass this gap.

Automatic approval review stopped the proposed external model calls before execution,
pending explicit disclosure authorization. The planned budget remains six fixed cases,
at most twelve model HTTP attempts and at most US$10. No live calls have been made and
no favorable-result retries have occurred. The report can be updated with actual
outcomes after authorized execution; expected answers must not be relabeled as results.
