"""Exact feasibility and upper certificate for Li–Yang Theorem 4.1."""

from fractions import Fraction as F
import json
import re
import sys


def run(inputs):
    names = ("alpha_I", "alpha_Q", "alpha_phi", "delta")
    if type(inputs) is not dict or set(inputs) != set(names):
        raise ValueError("provide alpha_I, alpha_Q, alpha_phi and delta")
    point = {}
    for name in names:
        value = inputs[name]
        if (type(value) is not str or len(value) > 64
                or not re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value)):
            raise ValueError("parameters must be bounded canonical rational strings")
        point[name] = F(value)
        if str(point[name]) != value:
            raise ValueError("parameters must be reduced canonical rational strings")

    i, q, phi, delta = (point[name] for name in names)
    u = delta - i
    terms = (i / 3, 2 - 2 * phi + q, 4 - 4 * phi,
             3 + phi, 5 - q, (5 - 2 * phi + q) / 2)
    slacks = [term - u for term in terms]
    passed = (0 <= q <= i and 0 <= phi <= F(1, 2) and delta >= 0
              and all(slack >= 0 for slack in slacks))

    # The last three rows combine to u <= 31/10 + phi/2.
    # Canceling phi/2 bounds the circuit coefficient for this unchanged LP.
    weights = (F(7, 10), F(1, 10), F(1, 5))
    rows = ((F(3), F(1), F(0)),
            (F(5), F(0), F(-1)),
            (F(5, 2), F(-1), F(1, 2)))
    combined = tuple(sum(w * row[j] for w, row in zip(weights, rows))
                     for j in range(3))
    if sum(weights) != 1 or combined != (F(31, 10), F(1, 2), F(0)):
        raise ValueError("invalid protected upper certificate")
    return {
        "passed": passed,
        "coefficient": str(u - phi / 2),
        "constraint_slacks": [str(slack) for slack in slacks],
        "dual_weights": [str(weight) for weight in weights],
        "dual_bound": str(combined[0]),
        "message": "Feasible LP point." if passed else "Infeasible LP point.",
    }


if __name__ == "__main__":
    try:
        result = run(json.load(sys.stdin))
    except (ValueError, TypeError, ZeroDivisionError, RecursionError) as exc:
        result = {
            "passed": False,
            "coefficient": "",
            "constraint_slacks": [],
            "dual_weights": [],
            "dual_bound": "",
            "message": str(exc)[:1000],
        }
    print(json.dumps(result, separators=(",", ":"), allow_nan=False))
