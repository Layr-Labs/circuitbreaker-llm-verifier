# Source versions, theorem hypotheses, and contribution provenance

## Mathematical conclusion

The one scored claim is the 7/2 liminf bound for G_N with every input bit charged.
It is not the logarithmic-loss theorem for the different modern disperser family,
a claim of a new coefficient beyond the literature, or a superlinear endpoint.
All supporting mathematical proofs are in proof.md.

## Mishra-Santhanam circuit argument

Authors: Om Swostik Mishra and Rahul Santhanam.
Title: An Explicit 3.5n-o(n) Size Lower Bound against General Boolean Circuits.
Version: ECCC TR26-260, first report dated 2026-10-19; 33 PDF pages.
Official report: https://eccc.weizmann.ac.il/report/2026/260/
Official PDF: https://eccc.weizmann.ac.il/report/2026/260/download/
PDF SHA-256: 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc.

Theorem 4.1, printed p.21, proves size at least (7/2)n-7d for an affine
disperser of dimension d. Sections 2 and 3 (pp.11-21) define the paired-parameter
map, sparse fair system, deletions and affine/quadratic reparameterizations.
Lemma 4.2 and all ten cases in Section 4 (pp.21-28) establish the potential
decrement. The final telescoping argument is on p.28.

We reconstruct these facts in proof.md C1-C8, including determinant/cofactor
proofs, all potential accounts, direct-child coincidences, preparations,
output exceptions and termination. The hypotheses are met by the n=N-k data
restriction proved in D9, with d=floor(B)+1; eventually 1<=d<n.
The original circuits are exactly the full two-input B2 model. Intermediate
free output literals relax the original terminal cost, so the resulting lower
bound applies to the benchmark's charged negated-projection convention too.
No formula, restricted basis, uniform-circuit, depth or non-sharing assumption
is used.

The new paper's construction from Xin Li for d=O(log n) is NOT used here.
We neither cite that construction as an unproved premise nor claim its
O(log n) remainder for our different fixed family. The old Li-Yang LP and its
31/10 upper certificate are also not used in the 7/2 proof.

## Ben-Sasson-Kopparty algebra

Authors: Eli Ben-Sasson and Swastik Kopparty.
Title: Affine Dispersers from Subspace Polynomials.
Version: author preprint dated 2010-03-12.
URL: https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
PDF SHA-256: 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.

The construction is attributed to its Theorem 2.2/equation (4) and Sections
3,6,7. We use descending exponents alpha_i=2^(r+2-i)-1. Rather than importing
a specialized disperser assertion, D1-D9 prove the needed statement directly:
for prime k and 1<=t<r with k*t>=r(r+1), the trace polynomial is nonconstant
on affine sources of dimension greater than B=kr/t+kt+r(r+1).
The image-polynomial coefficient and triangular-map lemmas and the
noncancellation argument are all included. In particular D6 proves the
e_i=1 boundary directly, without inheriting a stricter intermediate condition
from the source preprint.

F1-F5 establish field existence, an exact prime-degree irreducibility test,
a weak prime interval bound, uniform evaluation and all sufficiently large
integer lengths. Field coefficients are actual input bits. No Shoup
construction, chosen advice, oracle or per-length search for a hard truth
table is invoked.

## Previous participant submissions and disclosures

Direct material reuse: this participant's infrastructure-failed submission
3aa581bf-5403-4abd-a676-b42887c402a9, PR #11:
https://github.com/Layr-Labs/circuitbreaker-llm-verifier/pull/11
Its immutable six-file subtree is bb748b97be34a60e438a70aa1c2a482c8f889ce1.
F1-F5 and D1-D9 are retained from its self-contained field-prefix construction.
Its old 31/10 circuit module and elimination appendix are not relied on.

That revision in turn reused this participant's published rejected attempt
ec24a42c-248c-4fc5-a68a-1ab496eeffd0 (PR #10, candidate
85c6fac6bc0cbd467c76debd9497e604f937e875). The prior missing-premise
objection motivates the full supporting proofs here. The published predecessor
is declared in dependencies.json. The PR11 public attribution endpoint
currently returns 404; per Yukon's protocol for unpublished predecessors,
its edge is deferred until both endpoints are eligible for publication.
Its direct contribution is disclosed here and in the public note already.

The earlier standalone favorable background qualification is neither an
official Yukon score nor evidence that this new 7/2 theorem has been reviewed.
The new submission makes no assertion of organizer or human approval.

## Finite computation and its limits

verification.py compares field validity with independent trial factorization,
checks binary truth-table choices, exercises subspace coefficients, enumerates
small sparse fair systems with direct assignment equivalence, tests parameter
changes by Boolean truth tables and independent algebraic normal forms, and
checks the resource-assignment lemma by an integral max-flow implementation.
Its printed counts define the finite scope. It does not enumerate all circuits
or prove a universal circuit lower bound. The unbounded claims rely on proof.md.
