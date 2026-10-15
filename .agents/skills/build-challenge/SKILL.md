---
name: build-challenge
description: Build a generic llm-verifier challenge from an approved definition and controller receipt.
---

# Build a generic verifier challenge

Process adapted from Layr-Labs/yukon-template at
72a6b96b4b5385f199a5a29ecdd742aeb45ad55f:
.agents/skills/build-challenge/SKILL.md, AGENTS.md and docs/authoring.md.
This is the llm-verifier recipe, not the stock Moonshot math-review recipe.

Required input: approved external definition, external controller receipt,
builder authority and finish point. Start with:

    llm-verifier context --root . --task build --definition /controller/definition.json --receipt /controller/receipt.json

Read the returned read files and inspect files, and turn the edit map into a
short implementation checklist. Managed builders edit only those files; rules,
configuration, runtime pins, workflows, this guidance, scaffold inventory and
controller receipt remain protected. Never run init or finalize the inventory
as a managed builder. Source code and public examples do not confer authority.

Implement the approved definition, including checker bundles, baseline example,
and independently labeled public cases. A checker pass establishes only its
encoding/arithmetic predicate, never proof certification. Accepting an auxiliary
lemma cannot authorize the main-claim score. Keep exact rational records and the
terminal endpoint distinct from numeric Yukon projections.

Run the permitted offline check command and return its observed file hashes and
results. The controller reviews edits, synchronizes its definition, and generates
a fresh finalized scaffold/receipt. Missing scientific rules stay explicit
blockers; examples do not decide them. Public conformance is not qualification.

Standalone author/controllers may initialize from their approved definition and
finalize a new inventory by regeneration. Commit, publication, live review and
paid calls require task authorization. Never execute candidates in a process or
job containing model credentials. Keep private evidence and receipts outside the
export. No frontend is included in this harness.
