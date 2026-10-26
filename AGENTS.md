<!-- llm-verifier:solver-guidance:begin -->
## Preparing a submission

Read the [solver guide](https://github.com/Layr-Labs/llm-verifier/blob/8c1f098875fae0750c29ecfd641f7c9fc53a0f41/docs/solving.md) and
[solver skill](https://github.com/Layr-Labs/llm-verifier/blob/8c1f098875fae0750c29ecfd641f7c9fc53a0f41/skills/llm-verifier-solver/SKILL.md)
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
