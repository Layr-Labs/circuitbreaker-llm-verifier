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
> that the argument must achieve. Then the **leading gate coefficient**
> `theta - alpha` of the resulting lower bound is at most `4`.
>
> More precisely: if the row estimates of `elimination.md` Section 7 are the
> only constraints, then `theta - alpha <= 4`, and equality holds exactly on
> the family
>
>     alpha = beta + 2 ,   beta >= 4 ,   theta = beta + 6 ,
>
> so the coefficient `4` is attained, for example at `alpha = 6`, `beta = 4`,
> `theta = 10`, which is the point used in `elimination.md`.

**What is bounded, and what is not.** The bounded quantity is the
*difference* `theta - alpha`, which by (O.1) below is the leading coefficient
of the gate lower bound. The constant `theta` on its own is **not** bounded:
the constraints derived in Section 2 are compatible with arbitrarily large
`theta`, for instance with the whole family `alpha = beta + 2`, `beta >= 4`,
`theta = beta + 6`, for which `theta` grows without bound while
`theta - alpha = 4` stays fixed. A claim of the form "the largest admissible
`theta` is `4`" would be false and is not made here; the claim is about
`theta - alpha`.

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
| (D1)-(D4) | -1 | <= 0 | 0 | `-1 <= 0`, automatic |
| move affine gate | 0 | <= 0 | 0 | `0 <= 0`, automatic |
| (AR) | 0 | 0 | 0 | `0 <= 0`, automatic |

These are all ten rows of the table of `elimination.md` Section 7, with
`Delta a = 0`; the three rows with `Delta a = -1` are given below. The two
non-automatic constraints are therefore both

    alpha <= beta + 2 .                                           (O.2)

**The rows with `Delta a = -1`** must satisfy `Delta V <= -theta`, i.e.
`theta <= -Delta V`. These are exactly the three rows of `elimination.md`
Section 7 with `Delta a = -1`:

| Row | `Delta s` | `Delta I` | `Delta q` | `Delta a` | bound on `theta` |
|---|---:|---:|---:|---:|---|
| Case 4, both inputs unpaired | <= -2 | -2 | +1 | -1 | `theta <= 2 - beta + 2 alpha` |
| Case 5, `k = 1` paired input | <= -2 | -1 | -1 | -1 | `theta <= 2 + beta + alpha` |
| Case 5, `k = 2` paired inputs | <= -2 | 0 | -2 | -1 | `theta <= 2 + 2 beta` |

(In every row the tabulated value is the one that binds, because
`Delta s <= -2` is an upper bound on `Delta s` and `Delta V` is increasing
in `Delta s`; the weaker values `3 - beta + 2 alpha`, `3 + beta + alpha` and
`3 + 2 beta` follow from the slack and are never the minimum.)

Since `alpha >= 0`, the second of the three bounds is dominated by the first.
So the constraints on `theta` are

    theta <= min( 2 - beta + 2 alpha ,  2 + 2 beta ) .               (O.3)

**Bounding the leading coefficient.** By (O.1) the coefficient is
`theta - alpha`, and by (O.3),

    theta - alpha <= min( 2 - beta + 2 alpha ,  2 + 2 beta ) - alpha
                  = min( 2 - beta + alpha ,  2 + 2 beta - alpha ).   (O.4)

The first term of (O.4) is increasing in `alpha` and the second is
decreasing in `alpha`, so for a fixed `beta` the minimum of the two is
maximised where they cross:

    2 - beta + alpha = 2 + 2 beta - alpha
    <=>   alpha = 3 beta / 2 .                                    (O.5)

At that point the common value is `2 + beta / 2`, and (O.2) requires
`3 beta / 2 <= beta + 2`, i.e. `beta <= 4`. So for `beta <= 4` the largest
coefficient is `2 + beta / 2 <= 4`, with equality only at `beta = 4`.

For `beta > 4` the constraint (O.2) binds before (O.5) can be reached, so we
must take `alpha = beta + 2`; substituting into (O.4),

    theta - alpha <= min( 2 - beta + beta + 2 ,  2 + 2 beta - beta - 2 )
                  = min( 4 , beta )
                  = 4 ,                    for  beta >= 4 .       (O.6)

Combining the two ranges, the leading coefficient is at most `4` in all
cases, and it equals `4` exactly when `alpha = beta + 2` with `beta >= 4`.
This proves (O).

**Attainment.** Take `alpha = 6`, `beta = 4`. Then (O.2) holds with equality
(`6 = 4 + 2`), and (O.3) gives

    theta <= min( 2 - 4 + 12 ,  2 + 8 ) = min(10, 10) = 10 ,

so `theta = 10` is feasible, and the coefficient is `theta - alpha = 10 - 6 =
4`. These are exactly the weights used in `elimination.md`
(`W = s + 4q + 6I`, `Delta W <= 10 Delta a`), so the bound of (O) is attained
by the proof actually given.

The whole family of attaining points is `alpha = beta + 2`, `beta >= 4`,
`theta = beta + 6`: by (O.2) `alpha <= beta + 2` is met with equality, and
by (O.3) the three row bounds become `beta + 6`, `2 beta + 4` and
`2 + 2 beta`, whose minimum is `beta + 6` exactly when `beta >= 4`. Note
that `theta = beta + 6` is unbounded as `beta` grows, while the coefficient
`theta - alpha = 4` does not change; this is the point made in the note above
and is the reason the lemma is phrased in terms of `theta - alpha`.

## 3. Exact rational verification

For completeness, the constraint system above was also checked by exhaustive
enumeration over the grid `alpha, beta in {0, 1/8, 2/8, ..., 400/8}` with
exact rational arithmetic, using `Delta V = Delta s + beta * Delta q + alpha *
Delta I` and the requirement `Delta V <= theta * Delta a` row by row. The
maximum of `theta - alpha` found was exactly `4`, attained at
`(alpha, beta, theta) = (6, 4, 10)` and at every other point of the family
`alpha = beta + 2`, `beta >= 4`, `theta = beta + 6`; no point of the grid
gave a larger coefficient. This finite check is a diagnostic for the algebra
of Section 2; the argument of Section 2 is the proof.

## 4. What this lemma does and does not say

* It says: *within this family of potentials and these row estimates*, no
  reweighting improves the leading coefficient above `4`. The same statement,
  for the potential family and the row estimates of the earlier `7/2` analysis,
  appears as **Remark 4.3 of Mishra–Santhanam ECCC TR26-260 Revision 0**
  (19 October 2026), which reads: "Consider any potential of the form
  `V = s + βq + αI − σℓ`, where `I` is the number of used unpaired parameters,
  `ℓ` is the number of unpaired, isolated parameters, and the weights satisfy
  `α >= σ >= 0` and `β >= 0`. ... the weights chosen in our proof,
  `α = 7/2`, `β = 3/2`, and `σ = 1/2`, are optimal, that is, no other choice
  of weights can improve the coefficient `7/2`." The corresponding remark in
  **Revision 1** is Remark 4.3 as well ("Why this analysis gives lower bound
  coefficient 4"), for the two-term family `V = s + βq + αI` with `α = 6`,
  `β = 4`. This file reproduces that remark's argument from scratch, over the
  exact row estimates of `elimination.md`, and additionally computes the whole
  family of attaining points.
* It does **not** say that no circuit lower bound with coefficient above `4`
  exists. It is an optimisation bound for one proof technique, and a
  different technique, a different potential, or a different premise on the
  function could in principle do better.
* It does **not** say that the constant `10` multiplying `d` in the class
  theorem of `elimination.md` is optimal; only that the *leading* coefficient
  is. It also does not bound `theta` itself, which is unbounded above on this
  family of row estimates.
* It is an unscored lemma. The scored result is the joint theorem in
  `bound.json`, with coefficient `4` and exponent `1`.
