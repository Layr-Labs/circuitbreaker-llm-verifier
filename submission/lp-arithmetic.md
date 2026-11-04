# Exact arithmetic: the LP optimum 31/10, the remainder constants, and the checker

Every number that enters the final inequality of this submission is computed here
with exact rational arithmetic (`fractions.Fraction`, no floating point). Nothing
in this file is a circuit theorem: it certifies arithmetic only. In particular
the LP optimality below says that `31/10` is the exact optimum of the
**unchanged** Li–Yang Theorem 4.1 linear program; it is not a barrier, not a
lower bound, and not a claim that no other technique can do better.

## 1. The LP (Li–Yang Theorem 4.1, unchanged)

With `u = δ − α_I`, Theorem 4.1 gives, for every packing `P` of a `B2` circuit
computing an affine disperser for dimension `d`,

    u ≤ min{ α_I/3,  2 − 2α_φ + α_Q,  4 − 4α_φ,  3 + α_φ,  5 − α_Q,
             (5 − 2α_φ + α_Q)/2 }.

The final coefficient is `c = δ − α_I − α_φ/2 = u − α_φ/2`. The optimisation
problem solved in the cited proof is

    maximise  u − α_φ/2
    subject to u ≤ α_I/3,  u ≤ 2 − 2α_φ + α_Q,  u ≤ 4 − 4α_φ,
               u ≤ 3 + α_φ,  u ≤ 5 − α_Q,  u ≤ (5 − 2α_φ + α_Q)/2,
               0 ≤ α_φ ≤ 1/2,  α_I, α_Q ≥ 0,  0 ≤ α_Q ≤ α_I,
    with δ = α_I + u ≥ 0.

## 2. The feasible point used in this submission

    α_I = 48/5,   α_Q = 9/5,   α_φ = 1/5,   δ = 64/5,   u = δ − α_I = 16/5.

The six right-hand sides, by exact rational arithmetic:

    48/5 / 3                  = 16/5
    2 − 2/5 + 9/5             = (10 − 2 + 9)/5      = 17/5
    4 − 4/5                   = 16/5
    3 + 1/5                   = 16/5
    5 − 9/5                   = 16/5
    (5 − 2/5 + 9/5)/2         = (32/5)/2            = 16/5

Five of the six terms equal `16/5`; the second term `2 − 2α_φ + α_Q` equals
`17/5`. Hence the minimum is `16/5 = u`, and the slacks `term − u` are
`(0, 1/5, 0, 0, 0, 0)`, all nonnegative, so the point is feasible with five
tight rows:

    δ = α_I + u = 48/5 + 16/5 = 64/5  as required. Admissibility:
`0 ≤ 9/5 ≤ 48/5`, `0 ≤ 1/5 ≤ 1/2`, `δ = 64/5 ≥ 0`. The resulting coefficient is

    c = u − α_φ/2 = 16/5 − 1/10 = 32/10 − 1/10 = 31/10.

## 3. Upper certificate: no admissible point beats 31/10

Take the three rows with right-hand sides `3 + α_φ`, `5 − α_Q` and
`(5 − 2α_φ + α_Q)/2`, i.e. as linear upper bounds on `u`,

    u ≤ 3 + α_φ,     u ≤ 5 − α_Q,     u ≤ (5 − 2α_φ + α_Q)/2.

Weight them by `7/10`, `1/10`, `1/5` (weights are nonnegative and sum to `1`), so
that the inequalities may be combined:

    (7/10)(3 + α_φ) = 21/10 + (7/10)α_φ
    (1/10)(5 − α_Q) =  5/10 − (1/10)α_Q
    (1/5)(5 − 2α_φ + α_Q)/2 = (5 − 2α_φ + α_Q)/10 = 5/10 − (1/5)α_φ + (1/10)α_Q

Sum: constant term `21/10 + 5/10 + 5/10 = 31/10`; `α_φ` coefficient
`7/10 − 1/5 = 7/10 − 2/10 = 1/2`; `α_Q` coefficient `−1/10 + 1/10 = 0`. Hence

    u ≤ 31/10 + α_φ/2      and therefore      c = u − α_φ/2 ≤ 31/10

for **every** admissible parameter tuple. Combined with §2, the exact optimum is
`31/10`, and it is attained at `(48/5, 9/5, 1/5, 64/5)`. This is the unscored
lemma B of `proof.md` §8: finite rational arithmetic, awarding no circuit
coefficient by itself.

## 4. The constants of the final inequality

Used in `gate-elimination.md` §3 and `proof.md` §5:

    coefficient:   64/5 − 48/5 − 1/10 = 16/5 − 1/10 = 31/10
    d-coefficient: 2δ + α_φ·(5/2) = 2·(64/5) + (1/5)(5/2) = 128/5 + 1/2 = 261/10
    constant:      2δ = 128/5

so

    |C| ≥ (31/10)·n − (261/10)·d(n) − 128/5.

With `d(n) = ⌊6n^{4/5}⌋ + 1`:

    R(n)/n = (261/10)·d(n)/n + (128/5)/n ≤ (261/10)·(6n^{−1/5} + 1/n) + (128/5)/n → 0,

so for every `ε > 0` there is `N_ε` with `R(n) ≤ εn` for all `n ≥ N_ε`, giving
`C_B2(f_n) ≥ (31/10 − ε)n` for every integer `n ≥ N_ε`. (Variant with the smaller
potential bound `d` displayed in Theorem 3.9: `128/5 + 1/5 = 129/5`; same
coefficient, as recorded in `gate-elimination.md` §5.)

## 5. Unscored lemma A, arithmetic part: the OR family

For `g_n = OR_n` with `g_0 = 0` (proof in `proof.md` §7): `C_B2(g_n) = n − 1` for
`n ≥ 2`, so `liminf C_B2(g_n)/n = 1`, the pair `(1,1)` is valid for that family,
and `(2,1)` is false because a binary OR tree has `n − 1` gates. This is a
different family from the main theorem's and supplies neither coordinate of the
main claim.

## 6. The supplied checker, run

`submission/submission.json` declares one `lp-certificate` check with the exact
rational inputs `alpha_I = 48/5`, `alpha_Q = 9/5`, `alpha_phi = 1/5`,
`delta = 64/5`. Running the supplied `checkers/lp_certificate.py` on them (no
model access, no network) prints

    {"passed":true,"coefficient":"31/10",
     "constraint_slacks":["0","1/5","0","0","0","0"],
     "dual_weights":["7/10","1/10","1/5"],"dual_bound":"31/10",
     "message":"Feasible LP point."}

(line breaks added for readability; the five zero slacks are the five tight rows
of §2, the `1/5` is the slack of the `17/5` row). `passed` confirms feasibility of the declared
point and of the protected dual certificate `(31/10, 1/2, 0)`; `coefficient`
confirms the exact value `31/10`. This certifies arithmetic only. It is not an
eligibility argument and not a circuit theorem, and the challenge score does not
depend on it.

## 7. What is deliberately not claimed here

* No assertion that `31/10` is optimal for anything other than the unchanged
  Li–Yang LP described in §1.
* No numeric `Infinity`, no endpoint, no `p > 1`, no logarithmic factor.
* No claim that passing the checker establishes the main theorem; the main
  theorem's premises are E1 (`disperser.md`) and E2 (`gate-elimination.md`), and
  their validity is a mathematical question for the reviewer, not an arithmetic
  one.
