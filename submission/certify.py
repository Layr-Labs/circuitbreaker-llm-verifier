"""Replay both parts of the numerical estimate and check the declared coefficient."""
import json
import sys
from fractions import Fraction
from pathlib import Path

script = Path(__file__).resolve()
files = script.parent / "files" if script == Path("/snapshot/entrypoint.py") else script.parent
sys.path.insert(0, str(files))
from cavity import cavity
from reference import reference_upper_bound


def main():
    reference = reference_upper_bound()
    counts, denominator = cavity(48)
    assert counts == [797485233, 439228269, 317972502, 199325451]
    assert denominator == 442368000000
    witness = Fraction(sum(counts), denominator)
    penalty = Fraction(17, 30)
    debit = max(0, 9 - 19 * penalty + 6 * penalty**2)
    inflation = 1 + debit**2 / (16 * penalty * (3 - penalty) * (9 - 7 * penalty + 2 * penalty**2))
    width = inflation * (reference - (1 - penalty) * witness)
    limiting = 1 + 1 / (2 * width)
    with (files / "bound.json").open() as stream:
        coefficient = json.load(stream)["guarantee"]["leading"]["coefficient"]
    declared = Fraction(int(coefficient["numerator"]), int(coefficient["denominator"]))
    assert declared == Fraction(487729, 100000) and declared < limiting
    print(json.dumps({"reference_upper": str(reference), "grid": 48,
                      "event_counts": counts, "common_denominator": denominator,
                      "witness_lower": str(witness), "penalty": str(penalty),
                      "layout_factor": str(inflation), "width_upper_limit": str(width),
                      "coefficient_limit": str(limiting), "declared": str(declared)}))


if __name__ == "__main__":
    main()
