# Source-only migration to coefficient/exponent objectives

This change prepares a **new registration**. Do not reimport or reorg the live
single-objective `circuitbreaker` registration (`main` track) using this manifest: the existing benchmark cannot
change its objectives through a normal reorg. No live registration, benchmark,
remote setting, UI, merge or deployment is part of this source change.

After independent exact-source review and authorization to register, the new
registration must use the frozen source containing:

- Objectives `gate_coefficient` and `gate_exponent`, both maximized, in that order.
- The new scoring/review versions, current fixed-family validity policy,
  atomic acceptance, and unchanged provider/resource configuration.
- The pinned shared verifier and regenerated hosted caller with the bounded
  metrics overlay. The runner must select `circuitbreaker-score`.
- The published (31/10,1) seed declaration and its unchanged proof evidence,
  reviewed under the current supplied-source policy. Passing structural tests
  or accepting an older baseline does not guarantee scientific acceptance.

Choose the new registration's name/slug and archival/cross-link treatment for
the existing benchmark at that later authorized step. Preserve historical
records and their policy versions. Do not rewrite old single-coordinate
results into accepted two-coordinate results without review of the exact joint
theorem and source. Results tied under full-expression comparison retain both
proofs; a weaker projection of an existing expression is not fresh progress.

The generic full-expression table is available; the coordinated platform
extension registers `circuitbreaker-v1` for these lower-bound guarantees. Once
that extension is available in the target environment, the reviewed
[Complexity tab fragment](examples/complexity-result-tab.json) can be selected
in page or track `presentation.resultTabs`. It reads the additive table envelope
specified in [metrics.md](metrics.md), including accepted null-score endpoints
and complete relevant current accepted history. The original rich metric stays
available for audit. The numeric graph remains a [coefficient, exponent]
projection, while full expressions determine table ordering and ties.

This source only prepares the producer and presentation fragment. Publishing
the source, releasing the registered extension, selecting the page fragment and
creating the two-objective registration are separate authorized actions. Do not
apply the fragment to historical one-objective data as a substitute for its
missing joint review, and do not change the live registration or page as part
of source preparation.

No model qualification was performed for this restructure. Before launching a
new registration, separately authorize and run the required reviewer/scientific
qualification on its exact source, then retain the selected evaluation evidence
and independent approval. A new evaluation must not be presented as a
mathematical discovery because a synthetic fixture passed.
