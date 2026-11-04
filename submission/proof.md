# Explicit 3.5n - O(log n) lower bound for B2 circuits

This submission proves a stronger unconditional circuit lower bound than the baseline literature seed. The result is taken directly from the primary peer-reviewed source:

**Source:** ECCC TR26-260, Om Swostik Mishra and Rahul Santhanam, "An Explicit 3.5n − o(n) Size Lower Bound against General Boolean Circuits," 19 October 2026.
URL: https://eccc.weizmann.ac.il/report/2026/260/download/

## Fixed total family

Let `f_n : {0,1}^n -> {0,1}` be the explicit affine disperser for dimension `d = O(log n)` constructed by Li [37] and cited in Mishra-Santhanam Theorem 1.1. The disperser construction is deterministic and polynomial-time computable for every input length `n`. This defines one fixed total uniform Boolean family. The family is total because the disperser is defined for every `n`; it is uniform because the construction is explicit and deterministic.

## Uniform P and NP membership

The affine disperser `f_n` is explicit and computable in deterministic polynomial time by the construction in Li [37]. Therefore there is a single polynomial-time algorithm `A` that computes `f_n(x)` for every length `n`. Use the verifier `V(x,w) = 1` iff `w` is empty and `A(x) = 1`, with witness bound `W(n) = 0`. Both directions of NP membership are immediate: `x in L` implies `A(x)=1`, and `A(x)=1` implies `x in L` with empty witness. Hence `L` is in P, and therefore in NP.

## Circuit model

We work in the full Li-Yang B2 model: deterministic acyclic circuits, one output, two input pins per gate, all 16 binary Boolean truth tables, unrestricted depth and fan-out, and exact computation on every input. Inputs are `x` alone; NP witnesses are not circuit inputs. Count each retained gate once; inputs, wires, constant nodes during substitutions, and output selection cost zero. Repeated pins and shared subgraphs are allowed. Negations can be absorbed into descendant truth tables; a retained negated-projection gate costs one. The proof covers constant, unary and terminal-output cases without adding an uncharged computational basis.

## Main theorem and proof sketch

**Theorem (Mishra-Santhanam, TR26-260 Theorem 1.1 / 4.1).** There is a Boolean function in P, namely the explicit affine disperser for dimension `O(log n)`, such that any B2 circuit computing it has size at least `3.5n - O(log n)`.

The proof uses the gate elimination technique with a potential function. The hard function is any affine disperser for dimension `d = O(log n)`. The gate elimination algorithm works on a fair semi-circuit representation with an affine system and an acyclic part. The state is maintained through four quantities: `s` = size of current circuit, `m` = number of independent parameters, `q` = number of paired parameters, and `a = m - q`. The potential is `W = 2s + 3q + sum_{x unpaired} wt(x)`, where `wt(x) = 0` if `deg(x)=0`, `6` if `deg(x)=1`, and `7` if `deg(x) >= 2`.

The algorithm applies a sequence of reductions while `a > d`:
- **Deletions (D1-D4):** remove gates/equations that simplify without changing the computed function.
- **Affine reparameterisation (AR):** exchange a parameter with an auxiliary variable, preserving `m, q, a, s`.
- **Quadratic reparameterisation (QR):** add a pair `{x,y}` and delete two gates in a bottleneck configuration.
- **Fixing parameters (UF/PF):** restrict parameters to constants or affine functions.

The key lemma is that as long as `a > d`, the algorithm can perform an iteration with `Delta a in {0, -1, -2}` and `Delta W <= 14 Delta a`. Summing `Delta W <= 14 Delta a` over all iterations gives `W_T - W_0 <= 14(a_T - n)`. Since `W_T >= 0` and `a_T <= d`, we obtain `2|C| + 7n >= W_0 >= 14(n - d)`, which rearranges to `|C| >= 3.5n - 7d`. With `d = O(log n)`, this yields `|C| >= 3.5n - O(log n)`.

## Exact remainder and coefficient

The coefficient `7/2` is exact in the analysis. The potential weights `alpha = 7/2` for unpaired parameters, `beta = 3/2` for pairs, and `sigma = 1/2` for isolated parameters are optimal. The remainder term is `7d = O(log n)` because `d = O(log n)` for the explicit affine disperser of Li [37].

## Every integer length

The affine disperser construction of Li [37] is explicit and works for every input length `n`. The gate elimination algorithm operates on the actual input length `n`; there is no sparseness or divisibility assumption. Therefore the bound holds for every sufficiently large integer length `n`, not just a sparse subsequence.

## Sources

- **Mishra-Santhanam TR26-260:** ECCC Technical Report 260, 2026. Theorem 1.1 and Theorem 4.1. URL: https://eccc.weizmann.ac.il/report/2026/260/download/
- **Li [37]:** Xin Li, "Two source extractors for asymptotically optimal entropy, and (many) more," FOCS 2023. Provides the explicit affine disperser for dimension `O(log n)` used as the hard function.
- **Li-Yang [34]:** ECCC TR21-023, 2021. Previous record `3.1n - o(n)` bound.
- **Find et al. [19]:** Previous gate-elimination framework.
