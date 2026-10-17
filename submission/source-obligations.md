# Sources, scope, and response to the earlier objections

The main theorem is in proof.md; its gate-elimination appendix is elimination.md.
Together they supply the mathematical arguments used in the score. Citations
below identify provenance and correspondence, rather than replace material
premises. All constants and hypotheses used by this instance are checked in F,
D, C and E of that proof package.

## Affine disperser

Eli Ben-Sasson and Swastik Kopparty, *Affine Dispersers from Subspace Polynomials*,
author manuscript dated March 12, 2010.
Source: https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
SHA-256: 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf

Sections 3, 6 and 7, particularly Theorem 6.1, Lemmas 7.1-7.2 and Theorem 7.3,
are the sources of the construction and proof strategy. Module D supplies the
characteristic-two supporting arguments: image subspace polynomials; cyclic
zero-block bounds; permanent coefficients with cyclic trace shifts; explicit
nonzero column choices; block echelon parameterization; preservation of maximal
binary weights; and the resource-assignment argument excluding other terms.
D6 proves its coefficient assertion for every positive e_i, including e_i=1;
it does not rely on the preprint Lemma 7.2's printed e_i>1 hypothesis.

The exact D8 hypotheses are prime k, 1<=t<r, k*t>=r*(r+1), and nonzero mu.
F supplies these with mu=1, an odd prime k, t=ceil(sqrt(r)), and the stated
exceptional-length tests. Zero padding is injective. The dimension threshold
actually used is B=k*r/t+k*t+r*(r+1), not an unverified numerical corollary.
The previous draft misnamed Ben-Sasson as Bourgain; the attribution is corrected.

## Circuit argument

Jiatu Li and Tianqi Yang, *3.1n-o(n) Circuit Lower Bounds for Explicit Functions*,
ECCC TR21-023, February 20, 2021.
Source: https://eccc.weizmann.ac.il/report/2021/023/download/
SHA-256: 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062

Theorem 1.1 is the published coefficient, and Theorem 4.1 with its Cases 0-8 is
the source of the elimination argument. Sections 2-3 supply the semicircuit,
restriction, packing, and normalization framework. C and E re-prove the
specialized argument at integer weights (5,48,9,1;64), including the nonempty
packing cases, all later branches, and the exceptional reconstruction cases.

C1-C2 address all sixteen tables and terminal negated projections. C3 counts
actual pins and shared consumers. C5 supplies a rank-preserving way to eliminate
semantic constants in the cyclic prefix, rather than assuming that deleting an
arbitrary constant equation preserves fairness. C4 proves the required initial
potential bound by path/cycle matching and a constant affine restriction;
no external potential-counting lemma remains. C7 explicitly converts the
measure to a gate bound with remainder (261/10)d+128/5.

## Eligibility and actual input length

The former submission used a fixed deterministic Shoup field algorithm.
This revision changes to the family G_N whose k-bit prefix describes the field.
F1-F2 prove existence and an exact polynomial-time validity test directly.
The main claim expressly names this new family. Invalid prefixes are assigned
zero. F3 proves the weak prime-growth estimate needed for every sufficiently
large length. F5 gives a fixed deterministic evaluator and both directions of
the empty-witness NP equivalence. F6 charges all k prefix bits and all n=N-k
data bits and proves the epsilon quantifiers with no subsequence or advice.

Thus the prior objections to citation-only specialized premises are addressed
with supporting arguments in the submitted text. This does not assert that the
previous review was an approval, or that local checks certify the theorem.
The optional canonical LP checker establishes only its exact arithmetic.
The finite diagnostics in verification.py cover only the predicates they print;
the universal mathematical coverage is supplied by the written proofs.
