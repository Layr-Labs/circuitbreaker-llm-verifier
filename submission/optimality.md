# An unscored lemma: the leading coefficient 4 is optimal for this potential family

This file contains a **lemma only**. It carries no score, no coefficient and
no exponent. It is stated because it is a relevant, proved partial result:
it says exactly how much room the proof method of `elimination.md` leaves,
and where the remaining room would have to come from.

## 1. The statement

> **Lemma (O).** Let the gate-elimination argument of `elimination.md` be
> analysed with any potential of the form
>
>     V = s + beta * q + alpha * I ,          alpha >= 0, beta >= 0 ,
>
> where `s` is the number of equations plus acyclic gates, `q` is the number
> of parameter pairs, and `I` is the number of used unpaired parameters, and
> let `theta` be the constant in the inequality `Delta V <= theta * Delta a`
> that the argument must achieve. Then the leading coefficient of the
> resulting lower bound on the number of gates is at most `4`.
>
> More precisely: if the row estimates of `elimination.md` Section 7 are the
> only constraints, then the coefficient is at most `4`, and `4` is attained
> (for example at `alpha = 6`, `beta = 4`, `theta = 10`, which is the point
> used in `elimination.md`).

The lower bound that the framework produces has the form

    |C| >= (theta - alpha) * n - theta * d ,                       (O.1)

because `V_0 <= |C| + alpha * n` initially, `V` stays nonnegative, and the
telescoping of `Delta V <= theta * Delta a` over the `n - d` units of `a`
gives `|C| + alpha n >= theta (n - d)`.

## 2. The proof

We use exactly the row estimates of `elimination.md` Section 7, in the form
`(Delta s, Delta I, Delta q, Delta a)`.

**The rows with `Delta a = 0`** must satisfy `Delta V <= 0`, since
`theta * Delta a = 0`:

| Row | `Delta s` | `Delta I` | `Delta q` | constraint |
|---|---:|---:|---:|---|
| Case 1, `deg(p) = 0` | 0 | 0 | -1 | `-beta <= 0`, automatic |
| Case 1, `deg(p) >= 2` | -2 | +1 | -1 | `-2 + alpha - beta <= 0` |
| Case 2 | -1 | 0 | -1 | `-1 - beta <= 0`, automatic |
| Case 3 | -2 | +1 | -1 | `-2 + alpha - beta <= 0` |
| (D1)-(D4) | -1 | 0 | 0 | `-1 <= 0`, automatic |

The two non-automatic constraints are therefore both

    alpha <= beta + 2 .                                           (O.2)

**The rows with `Delta a = -1`** must satisfy `Delta V <= -theta`, i.e.
`theta <= -Delta V`:

| Row | `Delta s` | `Delta I` | `Delta q` | bound on `theta` |
|---|---:|---:|---:|---|
| GP, `k = 1` child | -2 | -2 | +1 | `theta <= 2 - beta + 2 alpha` |
| GP, `k >= 2` children | `-(1+k)` | -2 | +1 | `theta <= 1 + k - beta + 2 alpha` |
| Case 5, `k = 1` paired input | -2 | -1 | -1 | `theta <= 2 + beta + alpha` |
| Case 5, `k = 2` paired inputs | -2 | 0 | -2 | `theta <= 2 + 2 beta` |

The `k >= 2` row is weaker than the `k = 1` row, since `1 + k >= 2`. So

    theta <= min( 2 - beta + 2 alpha ,  2 + beta + alpha ,  2 + 2 beta ) .   (O.3)

By (O.1) the coefficient is `theta - alpha`, and by (O.3) and (O.2),

    theta - alpha <= (2 - beta + 2 alpha) - alpha
                  =  2 - beta + alpha
                  <=  2 - beta + (beta + 2)
                  =  4 .                                              (O.4)

This proves the first half of (O).

**Attainment.** Take `alpha = 6`, `beta = 4`. Then (O.2) holds with slack
(`6 <= 6`), and (O.3) gives

    theta <= min( 2 - 4 + 12 ,  2 + 4 + 6 ,  2 + 8 ) = min(10, 12, 10) = 10 ,

so `theta = 10` is feasible, and the coefficient is `theta - alpha = 10 - 6 =
4`. These are exactly the weights used in `elimination.md`
(`W = s + 4q + 6I`, `Delta W <= 10 Delta a`), so the bound of (O) is attained
by the proof actually given.

For completeness, the whole family of points attaining `4` is: `alpha =
beta + 2` with `beta >= 4`. Indeed, substituting `alpha = beta + 2` into
(O.3) gives `theta <= min(6 + beta, 4 + 2 beta, 2 + 2 beta) = min(6 + beta,
2 + 2 beta)`, since `2 + 2 beta <= 4 + 2 beta`. The coefficient
`theta - alpha` is then `min(6 + beta, 2 + 2 beta) - (beta + 2) = min(4,
beta)`, which equals `4` exactly when `beta >= 4`. (The grid point `alpha =
10, beta = 8, theta = 14` found by exhaustive search is one such point; it
is equivalent to `alpha = 6, beta = 4` in the sense that both give coefficient
`4`, differing only in which row binds: at `(6,4)` the `GP` row and the
`Case 5, k = 2` row both bind, while at `(10,8)` only the `GP` row binds.)

## 3. Exact rational verification

For completeness, the constraint system above was also checked by exhaustive
enumeration over the grid `alpha, beta in {0, 1/8, 2/8, ..., 80/8}` with
exact rational arithmetic, using `Delta V = Delta s + beta * Delta q + alpha *
Delta I` and the requirement `Delta V <= theta * Delta a` row by row. The
maximum coefficient found was `4`, attained at `(alpha, beta, theta) =
(6, 4, 10)`, `(10, 8, 14)` and the other points of the family described
above. This finite check is a diagnostic for the algebra of Section 2; the
argument of Section 2 is the proof.

## 4. What this lemma does and does not say

* It says: *within this family of potentials and these row estimates*, no
  reweighting improves the leading coefficient above `4`. The same statement
  appears, for the earlier `3.5n` analysis, as Remark 4.3 of the cited
  source.
* It does **not** say that no circuit lower bound with coefficient above `4`
  exists. It is an optimisation bound for one proof technique, and a
  different technique, a different potential, or a different premise on the
  function could in principle do better.
* It does **not** say that the constant `10` multiplying `d` in the class
  theorem of `elimination.md` is optimal; only that the *leading* coefficient
  is.
* It is an unscored lemma. The scored result is the joint theorem in
  `bound.json`, with coefficient `4` and exponent `1`.
