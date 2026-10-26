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

The future DCP-style full-expression table should consume the versioned metric
and trusted provenance documented in [metrics.md](metrics.md), including
accepted null-score endpoints and the entire relevant accepted history. Its
feature and comparator-registration interface have not landed; this source
adds no imagined registration key and makes no UI changes. Integration will
need a separately reviewed adapter for supported versions, current acceptance,
ties, unknown comparisons and endpoint rows. The numeric graph remains a
projection, so a table absence today must not be interpreted as no endpoint.

No model qualification was performed for this restructure. Before launching a
new registration, separately authorize and run the required reviewer/scientific
qualification on its exact source, then retain the selected evaluation evidence
and independent approval. A new evaluation must not be presented as a
mathematical discovery because a synthetic fixture passed.
