"""Replay the reference path enclosure, scalar constant checks, and cavity integrals."""
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


def main():
    if not __debug__:
        raise RuntimeError("Run the certificate with assertions enabled.")

    verify_cavity_constants(pi_enclosure())

    ref = reference_enclosures()
    U = ref["U"]
    assert U == Fraction(163309458029, 1250000000000)

    (c_pos1, c_neg1), den1, (c_pos2, c_neg2), den2 = cavity(72, 18, 20)
    assert (c_pos1, c_neg1) == (1105971666807, 719008242156)
    assert den1 == 597196800000000
    assert (c_pos2, c_neg2) == (14563630181805, 10375498660161)
    assert den2 == 19110297600000000

    p0 = Fraction(c_pos1, den1)
    p1 = Fraction(c_neg1, den1)
    p2 = Fraction(c_pos2, den2)
    p3 = Fraction(c_neg2, den2)
    witness = p0 + p1 + p2 + p3
    assert witness == Fraction(13889747654797, 3185049600000000)

    penalty = Fraction(35, 71)
    debit = max(0, 1 - 2 * penalty)
    inflation = 1 + 9 * debit ** 2 / (16 * penalty * (3 - 2 * penalty))
    assert inflation == Fraction(80089, 80080)

    width = inflation * (U - (1 - penalty) * witness)
    assert width == Fraction(201921725829734301021559, 1571976806400000000000000)
    limiting = 1 + Fraction(1, 2) / width
    assert limiting == Fraction(987910129029734301021559, 201921725829734301021559)

    with (files / "bound.json").open() as stream:
        coefficient = json.load(stream)["guarantee"]["leading"]["coefficient"]
    declared = Fraction(int(coefficient["numerator"]), int(coefficient["denominator"]))
    assert declared == Fraction(489253, 100000) and declared < limiting

    print(
        json.dumps(
            {
                "reference_upper": str(U),
                "grid": 72,
                "nt1": 18,
                "nu_panels": 20,
                "kind1_counts": [c_pos1, c_neg1],
                "kind1_denominator": den1,
                "kind2_counts": [c_pos2, c_neg2],
                "kind2_denominator": den2,
                "event_fractions": [str(p0), str(p1), str(p2), str(p3)],
                "witness_lower": str(witness),
                "penalty": str(penalty),
                "layout_factor": str(inflation),
                "width_upper_limit": str(width),
                "coefficient_limit": str(limiting),
                "declared": str(declared),
            }
        )
    )


if __name__ == "__main__":
    main()
