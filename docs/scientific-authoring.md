# Scientific authoring contract
This challenge was built with generic llm-verifier init, adapting the requested
Yukon build-challenge workflow. Read generated AUTHORING.md, .agents builder guidance
and the context command before editing. Their explicit edit inventory is authoritative.
The generator owns the runtime workflow and manifests. Challenge authors supply the
approved definition, protected rules, checker implementations, submissions and public
cases; there is no challenge-specific runtime fork.

The scientific contract fixes one total NP family, full acyclic fan-in-two B2 circuits,
actual input bits, all sufficiently large lengths and unconditional lower bounds.
The finite mathematical score is exact p/q with no normalization or fixed threshold.
A proved superlinear bound for one family is the semantic endpoint. The platform's
finite-only transport limitation is documented in docs/transport.md and rules.md.

The protected main-claim checker reconstructs complete typed canonical metadata and
proof text. The referee must compare that payload to the submitted artifacts and
accept its exact claim ID; acceptance of a separate lemma never scores the main claim.
The optional lp-certificate checker supplies an exact primal/dual arithmetic check,
not a lower-bound theorem. Candidate code and forged output are never trusted checks.

Public conformance scripts test orchestration, including synthetic terminal behavior;
protected checker known-answer cases test arithmetic. Neither is hidden qualification
or proof verification. Keep credentials, private cases and independent review receipts
outside this export. Official promotion needs human approval of exact artifacts.

The seed adapts scientific data from the existing private CircuitBreaker repository
at its pinned revision recorded in docs/sources.json. Historical approval is provenance
only. The source PDFs are cited, not redistributed. No new mathematical proof is
claimed by synthetic cases.
