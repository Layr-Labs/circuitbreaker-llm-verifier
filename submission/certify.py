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
    assert counts == [204011953695, 130843781187, 81407805552, 57106737168]
    assert denominator == 110592000000000
    witness = Fraction(sum(counts), denominator)
    penalty = Fraction(35, 71)
    debit = max(0, 1 - 2 * penalty)
    inflation = 1 + 9 * debit**2 / (16 * penalty * (3 - 2 * penalty))
    width = inflation * (reference - (1 - penalty) * witness)
    limiting = 1 + 1 / (2 * width)
    with (files / "bound.json").open() as stream:
        coefficient = json.load(stream)["guarantee"]["leading"]["coefficient"]
    declared = Fraction(int(coefficient["numerator"]), int(coefficient["denominator"]))
    assert declared == Fraction(48913, 10000) and declared < limiting
    print(json.dumps({"reference_upper": str(reference), "grid": 48,
                      "event_counts": counts, "common_denominator": denominator,
                      "witness_lower": str(witness), "penalty": str(penalty),
                      "layout_factor": str(inflation), "width_upper_limit": str(width),
                      "coefficient_limit": str(limiting), "declared": str(declared)}))


if __name__ == "__main__":
    main()
