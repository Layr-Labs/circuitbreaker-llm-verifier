"""Replay the reference enclosures, exact R=1 algebraic check, and cavity integrals."""
import json
import sys
from fractions import Fraction
from pathlib import Path

script = Path(__file__).resolve()
files = script.parent / "files" if script == Path("/snapshot/entrypoint.py") else script.parent
sys.path.insert(0, str(files))

import cavity as cavity_mod
from cavity import cavity
from reference import pi_enclosure, reference_enclosures


def verify_cavity_constants(pi_range):
    """Verify all seven integer scalar constants in cavity.py by exact rational arithmetic."""
    pi_hi = Fraction(pi_range.high)
    s = cavity_mod.SCALE
    assert cavity_mod.ROOT_LO ** 2 < 2 * s * s < cavity_mod.ROOT_HI ** 2
    assert 2 * cavity_mod.INV_R2_LO ** 2 < s * s
    assert 9 * cavity_mod.RHO_HI ** 2 > 8 * s * s
    assert cavity_mod.CG_NUM ** 2 * pi_hi < 24 * s * s
    assert 2 * cavity_mod.C2_NUM * pi_hi < 3 * s
    assert 8 * cavity_mod.C3_NUM ** 2 * pi_hi ** 3 < 108 * s * s


def verify_unit_layout_factor(lam, delta):
    """Verify by exact rational arithmetic that R(lam, delta) = 1 for all p in [0, 1/2]."""
    assert Fraction(0) < delta <= lam <= Fraction(1)
    r = (1 - lam) / (3 - lam)
    alpha_d = Fraction(3, 4) - lam * (1 + r) / 4
    alpha_s = Fraction(3, 4) + lam * (1 + r) / 4

    # Auxiliary inequalities guaranteeing nonnegative coefficients on Ab, Bb, zs, and Regime-1 zd:
    c_Ab = 1 - r / 4 - delta / (1 + delta)
    c_Bb = lam - r * (2 - lam) / 2 - (2 - r) * (1 - lam) / 4
    c_zs = lam - alpha_s / 2 - lam * delta * (1 - r / 4) / (1 + delta)
    c_zd_reg1 = (7 * lam - 3) / 4
    assert c_Ab >= 0 and c_Bb >= 0 and c_zs >= 0 and c_zd_reg1 >= 0

    # Regime-2 zd quadratic A*p^2 - B*p + C >= 0 for all real p via B^2 <= 4*A*C:
    A = alpha_d + lam * delta * r / (2 * (1 + delta))
    B = alpha_d - lam + lam * delta * (1 + r / 2) / (1 + delta)
    C = lam * delta / (1 + delta)
    disc = B * B - 4 * A * C
    assert A > 0 and (B <= 0 or disc <= 0)
    return {
        "r": str(r),
        "alpha_d": str(alpha_d),
        "alpha_s": str(alpha_s),
        "c_Ab": str(c_Ab),
        "c_Bb": str(c_Bb),
        "c_zs": str(c_zs),
        "c_zd_reg1": str(c_zd_reg1),
        "discriminant": str(disc),
    }


def main():
    if not __debug__:
        raise RuntimeError("Run the certificate with assertions enabled.")

    verify_cavity_constants(pi_enclosure())

    ref = reference_enclosures()
    U = ref["U"]
    W1 = ref["W1"]
    assert U == Fraction(1633094581, 12500000000)
    assert W1 == Fraction(4056504299, 25000000000)

    (c_pos1, c_neg1), den1, (c_pos2, c_neg2), den2 = cavity(56, 16)
    assert (c_pos1, c_neg1) == (590713922694, 370397261616)
    assert den1 == 321126400000000
    assert (c_pos2, c_neg2) == (6791760899379, 4741602879954)
    assert den2 == 8991539200000000

    s1 = Fraction(c_pos1 + c_neg1, den1)
    s2 = Fraction(c_pos2 + c_neg2, den2)

    lam = Fraction(27, 50)
    delta = Fraction(41, 20000)
    layout_checks = verify_unit_layout_factor(lam, delta)

    width = (1 - delta) * U + delta * W1 - (1 - lam - 2 * delta) * s1 - (1 - lam) * s2
    limiting = 1 + Fraction(1, 2) / width

    with (files / "bound.json").open() as stream:
        coefficient = json.load(stream)["guarantee"]["leading"]["coefficient"]
    declared = Fraction(int(coefficient["numerator"]), int(coefficient["denominator"]))
    assert declared == Fraction(19533, 4000) and declared < limiting

    print(
        json.dumps(
            {
                "U_upper": str(U),
                "W1_upper": str(W1),
                "grid": 56,
                "nt1": 16,
                "kind1_counts": [c_pos1, c_neg1],
                "kind1_denominator": den1,
                "kind2_counts": [c_pos2, c_neg2],
                "kind2_denominator": den2,
                "s1_lower": str(s1),
                "s2_lower": str(s2),
                "witness_total_lower": str(s1 + s2),
                "lambda": str(lam),
                "delta": str(delta),
                "layout_factor": "1",
                "layout_checks": layout_checks,
                "width_upper_limit": str(width),
                "coefficient_limit": str(limiting),
                "declared": str(declared),
            }
        )
    )


if __name__ == "__main__":
    main()
