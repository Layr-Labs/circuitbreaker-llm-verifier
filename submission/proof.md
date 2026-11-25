# A superlinear Boolean circuit lower bound for padded list distinctness

## Result and attribution

Let m_w=2w·2^w for integers w≥1. For an input x of N bits with N≥4, choose the largest w with m_w≤N. Parse the first m_w bits into 2^w consecutive records of 2w bits each. Accept exactly when these records are pairwise distinct. For N<4 always reject. Write f_N for this Boolean function and C_N for its minimum exact nonuniform B2 gate count.

We prove that C_N/(N log₂ log₂ N) tends to positive infinity over all integer lengths. In particular C_N≥N log₂ log₂ N for all sufficiently large N, the requested full guarantee has leading coefficient one, and the matching lane projection is divergence of C_N/N at exponent one. The input-hardness result is OpenAI's November 24, 2026 list-distinctness theorem. Our contribution here is its explicit specialization and transfer to this benchmark. We use neither a matrix-multiplication exponent assumption nor the arithmetic-circuit conclusion of the same paper.

## One polynomial-time language

An evaluator computes N, doubles the record count while tracking the next m_w, and compares every pair of records bit by bit. There are at most N records and each comparison reads at most N bits, so even a straightforward tape implementation takes polynomial time. This single deterministic algorithm handles all small lengths, malformed tails (all tail strings are ignored), repeated words and zeros. A verifier rejects every nonempty witness and returns the evaluator result on the empty witness. For x in the language that witness accepts; conversely any accepting witness is empty and the tested list has no repeated word. The witness bound is the zero polynomial. All N input bits are counted in the denominator.

## Single-list network theorem

Work over F₂. Write D_w for the indicator, as an element of F₂, of pairwise distinctness of 2^w binary records of length 2w. For some fixed real σ>0 and every sufficiently large integer w, any exact trivalent bit tensor network in the formal model described below that represents D_w and has v vertices satisfies

    m_w w^σ < v + m_w.                                         (1)

We use the field-specialized formal theorem for this normalization; a constant-factor change between vertex-plus-exterior size and trivalent volume would also suffice. The underlying published result is Theorem 2.2, Volume of list distinctness, in OpenAI, The matrix multiplication exponent is nine fourths, November 24, 2026, Sections 2.1–2.2 and The cost of testing many lists, revision 3344da59991ed94f9ba711217b29c1d760badfa6. The theorem allows t copies, and (1) only needs t=1. See network-proof.md for the exact premise, the supporting arguments, its formal proof path and the verification scope.

A bit tensor is a field-valued table on assignments of its incident binary indices. Each vertex here has at most three such indices. Joining wires means summing over their bit values the product of the table entries. Only the exterior bit values are supplied as inputs. All table entries in the circuit construction below are zero or one. Addition in F₂ is harmless because the number of consistent internal assignments at fixed exterior input will always be zero or one.

## A linear-size simulation of every benchmark circuit

Fix N≥4, its associated w, and any deterministic acyclic one-output B2 circuit computing f_N exactly with s gates. Fix the last N−m_w input bits to zero. The resulting network is required to expose precisely the m_w record bits, in their original order. Construct tables directly from the restricted truth tables and wires as follows.

* Keep only gates on an input-to-output dependency of the output, including any constant subexpressions on which the output formally depends; equivalently take all gate ancestors of the output. There are at most s. For each retained gate with truth table g, use a two-input, one-output table with entry one exactly when z=g(u,v). This works for every one of the sixteen truth tables, including constant tables, both projections and their negations. Retained gates use one vertex each, irrespective of fan-out. Add one unary test accepting exactly a one at the circuit's selected output. The test also works when that output selects a raw bit or constant directly.
* For each of the m_w free input bits, use a one-in, one-out identity table if it is read, or a one-in, zero-out table identically one if it is unused. This gives exactly one exterior wire for every named bit. A slot fed by a fixed padding bit or by a free constant source may instead receive its own one-output generator forced to the appropriate bit. There are at most 2s+1 such generators, since the retained gate inputs and the output test together have at most that many slots.
* If one producer must supply d≥1 slots, use a binary copy tree. Each copy has one input and two outputs and enforces equality of all three incident bits. Its d−1 internal vertices supply the d slots. Summing over producers costs at most 2s+1 copy vertices because the sum of all consumer slot counts is at most 2s+1. Two pins of the same gate reading the same producer count as two slots. Arbitrary sharing and arbitrarily large fan-out are covered. No extra gate is charged to the Boolean circuit for fixing a constant, distributing an input, choosing an output or eliminating a negation.

Thus the network has v≤5s+m_w+3 vertices. The network makes no independence assumption about shared values. For fixed free inputs, evaluate in a topological order to obtain the unique set of internal bit values consistent with the gate, input, constant and copy tables. Any satisfying contraction assignment agrees successively with these values. The final test admits that assignment exactly for output one. Consequently the contracted tensor is exactly D_w on every exterior bit string.

This network also has a literal construction in the formal Net syntax. First take the parallel product of the vertex tables, with all internal wire ends temporarily exterior. Use free input and output rewirings to order each internal wire's output end and input end at the same coordinate in a common bundle B. The entire product can now be typed from (lists w 1,B) to (unit,B), with all original named input positions in lists w 1. Trace B, which sums over all simultaneous internal wire assignments without changing the vertex count. The free rewires are generated by neighboring swaps, reassociation and unit insertion/removal, so arbitrary finite permutations and the balanced list bracketing are allowed. This construction also covers a direct output with no gates.

## All integer lengths and an explicit score

By (1), for every large enough w and the w selected at length N,

    m_w w^σ < 5s + 2m_w + 3.

As w tends to infinity we eventually have w^σ≥8 and m_w≥4, giving s≥m_w w^σ/10. Moreover m_{w+1}/m_w=2(w+1)/w≤4, so m_w≤N<m_{w+1}≤4m_w and w tends to infinity with N. For w≥1 we have log₂ N < w+2+log₂(w+1) ≤ 4w. Hence eventually

    C_N ≥ N (log₂ N)^σ / (40·4^σ).

All constants are fixed independently of N and of the circuit. For fixed σ>0, writing u=ln(log₂ N) gives (log₂ N)^σ/(log₂ log₂ N)=(ln 2)exp(σu)/u, which tends to infinity (already exp(σu)≥(σu)²/2). Dividing the last displayed lower bound by N log₂ log₂ N proves the asserted full guarantee and its lane projection. No numerical lower estimate on σ or practical cutoff length enters the score.
