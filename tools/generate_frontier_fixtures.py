"""Generate synthetic exporter parity cases; never a scientific results corpus.

Expected order is assigned by an explicit mathematical tier sequence, independently
of expressions.compare. Accepted envelopes come through the actual exporter and
full deterministic synthetic phase chain. Mutated/consumer-state cases are labelled.
"""

import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from export_metrics import projection
from expressions import CONVENTION, VERSION, claims_for, polynomial, ratio
from llm_verifier.errors import VerifierError
from tests.synthetic import synthetic_dossier
from tools.generate_hosted import VERIFIER_REVISION

FIXTURE = ROOT / "tests/fixtures/circuitbreaker-frontier-v1.json"


def finite(c, p):
    return {"kind": "finite", "coefficient": ratio(c), "exponent": ratio(p)}


def lane(p):
    return {"kind": "lane_endpoint", "exponent": ratio(p)}


def declarations():
    """Weakest to strongest, with an explicit shared tier for equal guarantees.

    Below n^p come negative leading log factors; then finite coefficients at
    n^p; then its bare omega endpoint; then positive leading log factors.
    All tiers at a smaller power precede those at a larger power. A global
    theorem is last. This independent oracle never calls the comparator.
    """
    entries = []
    tier = 0

    def add(identifier, guarantee, projected, *, same_tier=False, family=None):
        nonlocal tier
        if not same_tier:
            tier += 1
        entries.append((identifier, guarantee, projected, tier, family or identifier))

    def poly_group(p, q, r):
        for c in (F(1, 2), F(31, 10)):
            identifier = f"poly-{p}-{q}-{r}-{c}".replace("/", "_")
            if q > 0 or (q == 0 and r > 0):
                projected = lane(p)
            elif q == r == 0:
                projected = finite(c, p)
            else:
                projected = finite(1, (p + 1) / 2)
            add(identifier, polynomial(c, p, q, r), projected)

    for p in (F(1), F(3, 2), F(2), F(3)):
        if p > 1:
            poly_group(p, -1, -1)
            if p == 2:
                g = polynomial(1, 2, -1, 0)
                add("log-loss-lane", g, lane(F(3, 2)), family="shared-log-loss")
                add("log-loss-finite", g, finite(100, F(4, 3)), same_tier=True, family="shared-log-loss")
            poly_group(p, -1, 1)
            poly_group(p, 0, -1)
        if p < 3:
            poly_group(p, 0, 0)
            if p == 2:
                add("binary64-low", polynomial(2**54, p), finite(2**54, p))
                add("binary64-high", polynomial(2**54 + 1, p), finite(2**54 + 1, p))
            add(f"omega-{p}".replace("/", "_"), lane(p), lane(p))
        poly_group(p, 0, 1)
        poly_group(p, 1, -1)
        if p == 1:
            poly_group(p, 1, 0)
        poly_group(p, 1, 1)
    add("global", {"kind": "superpolynomial"}, {"kind": "global_endpoint"})
    assert len(entries) == 58
    return entries


def generate():
    contract = json.loads((ROOT / "examples/complexity-result-tab.json").read_text())["resultTabs"][0]
    cases = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "submission").mkdir()
        for filename in ("verifier.json", "rules.md"):
            shutil.copyfile(ROOT / filename, root / filename)
        shutil.copytree(ROOT / "checkers", root / "checkers")

        def emit(identifier, g, projected, family, *, valid=True, lemma=False, simulation=False):
            bound = {"version": VERSION, "convention": CONVENTION,
                     "family": f"SYNTHETIC transport family {family}; no circuit construction asserted.",
                     "theorem": "SYNTHETIC scripted evidence only; not a mathematical discovery or live review.",
                     "guarantee": g, "projection": projected}
            # These are valid filesystem basenames, not executable JS prototype
            # accessors. Include them in one actual accepted exporter case, then
            # explicitly remove them before subsequent snapshots.
            special = identifier == "poly-1-0-0-1_2"
            evidence = ["bound.json", "proof.md"]
            for name in ("constructor", "prototype", "__proto__"):
                path = root / "submission" / name
                if special:
                    path.write_text("SYNTHETIC inert evidence file.\n")
                    evidence.append(name)
                else:
                    path.unlink(missing_ok=True)
            manifest = {"version": 1, "title": "SYNTHETIC frontier parity case " + identifier,
                        "claims": ["SYNTHETIC unscored arithmetic lemma."] if lemma else claims_for(bound, sorted(evidence)),
                        "experiments": [], "checks": []}
            record = synthetic_dossier(root, bound, manifest, valid=valid, simulation=simulation)
            try:
                output = projection(record)
            except VerifierError as error:
                if not simulation:
                    raise
                assert str(error) == "simulation_cannot_publish_score"
                return None
            # Keep the exact native output and selected table namespace. Rich
            # metrics are checked independently, not duplicated in every row.
            if "metrics" in output:
                output["metrics"] = {"circuitbreaker_frontier": output["metrics"]["circuitbreaker_frontier"]}
            return output

        for identifier, g, projected, tier, family in declarations():
            cases.append({"id": identifier, "provenanceKind": "exporter", "consumerState": "accepted",
                          "expectedEligible": True, "strengthTier": tier,
                          "output": emit(identifier, g, projected, family)})
        for identifier, options in [("rejected", {"valid": False}), ("lemma", {"lemma": True}),
                                    ("simulation", {"simulation": True})]:
            cases.append({"id": identifier, "provenanceKind": "exporter", "consumerState": identifier,
                          "expectedEligible": False,
                          "output": emit(identifier, polynomial(1, 1), finite(1, 1), identifier, **options)})

    # Deliberate wire mutations test the consumer parser, never the producer's
    # willingness to export a forged dossier. Original exporter rows stay intact.
    original = cases[0]["output"]
    mutations = {
        "wrong-version": lambda m: m.update(version="unsupported-v2"),
        "wrong-comparator": lambda m: m.update(comparator="asymptotic-v1"),
        "wrong-scope": lambda m: m.update(scope="different-family-convention"),
        "wrong-model": lambda m: m.update(model="monotone-circuits"),
        "wrong-notation": lambda m: m.update(notation="big-o"),
        "bad-intake-hash": lambda m: m.update(intake_sha256="invalid"),
        "bad-dossier-hash": lambda m: m.update(dossier_sha256="A" * 64),
        "extra-envelope-field": lambda m: m.update(dominates_all=True),
        "bad-source-hash": lambda m: m["resources"]["circuit_size"].update(source_sha256="not-a-hash"),
        "empty-evidence": lambda m: m["resources"]["circuit_size"].update(evidence={}),
        "unsafe-evidence-path": lambda m: m["resources"]["circuit_size"].update(evidence={"../proof.md": "0" * 64}),
        "spaced-evidence-path": lambda m: m["resources"]["circuit_size"].update(evidence={"bound.json": "0" * 64, "proof.md": "1" * 64, "a b": "2" * 64}),
        "overlong-evidence-component": lambda m: m["resources"]["circuit_size"].update(evidence={"bound.json": "0" * 64, "proof.md": "1" * 64, "a" * 129: "2" * 64}),
        "overdeep-evidence-path": lambda m: m["resources"]["circuit_size"].update(evidence={"bound.json": "0" * 64, "proof.md": "1" * 64, "/".join(["a"] * 13): "2" * 64}),
        "git-evidence-path": lambda m: m["resources"]["circuit_size"].update(evidence={"bound.json": "0" * 64, "proof.md": "1" * 64, "a/.git/b": "2" * 64}),
        "missing-bound-evidence": lambda m: m["resources"]["circuit_size"].update(evidence={"proof.md": "1" * 64}),
        "missing-proof-evidence": lambda m: m["resources"]["circuit_size"].update(evidence={"bound.json": "0" * 64}),
        "overcount-evidence": lambda m: m["resources"]["circuit_size"].update(evidence={"bound.json": "0" * 64, "proof.md": "1" * 64, **{f"p{i}.md": "2" * 64 for i in range(31)}}),
        "long-statement": lambda m: m["resources"]["circuit_size"].update(statement="S" * 8001),
        "extra-resource": lambda m: m["resources"].update(other=copy.deepcopy(m["resources"]["circuit_size"])),
        "missing-circuit-size": lambda m: m.update(resources={}),
        "zero-coefficient": lambda m: m["resources"]["circuit_size"]["bound"].update(coefficient=ratio(0)),
        "unsupported-exponent": lambda m: m["resources"]["circuit_size"]["bound"].update(exponent=ratio(1, 2)),
        "nonreduced-coefficient": lambda m: m["resources"]["circuit_size"]["bound"].update(coefficient={"numerator": "2", "denominator": "2"}),
        "comparison-hint": lambda m: m["resources"]["circuit_size"]["bound"].update(comparison_hint="best"),
        "unsupported-kind": lambda m: m["resources"]["circuit_size"].update(bound={"kind": "exponential", "base": "2"}),
    }
    for identifier, mutate in mutations.items():
        output = copy.deepcopy(original)
        mutate(output["metrics"]["circuitbreaker_frontier"])
        cases.append({"id": identifier, "provenanceKind": "mutated-exporter", "consumerState": "accepted",
                      "expectedEligible": False, "output": output})
    for state in ("unreviewed", "revoked", "stale"):
        cases.append({"id": state, "provenanceKind": "consumer-state", "consumerState": state,
                      "expectedEligible": False, "output": copy.deepcopy(original)})

    eligible = [c for c in cases if c["expectedEligible"]]
    comparisons = [{"left": a["id"], "right": b["id"],
                    "expected": "stronger" if a["strengthTier"] > b["strengthTier"] else
                                "weaker" if a["strengthTier"] < b["strengthTier"] else "equivalent"}
                   for a in eligible for b in eligible]

    def scenario(identifier, ids):
        selected = [c for c in eligible if c["id"] in ids]
        best = max(c["strengthTier"] for c in selected)
        return {"id": identifier, "caseIds": ids,
                "expectedFrontierIds": [c["id"] for c in selected if c["strengthTier"] == best]}

    return {
        "contract": contract,
        "provenance": {"synthetic": True, "generator": "tools/generate_frontier_fixtures.py",
                       "verifierRevision": VERIFIER_REVISION,
                       "description": "Explicitly synthetic scripted phase chains, never live or reusable accepted results. Expected tiers are assigned independently of the Python comparator. Consumer-state cases require the current-result policy gate before metric parsing."},
        "cases": cases,
        "expectedEligibleIds": [c["id"] for c in eligible],
        "expectedFrontierIds": ["global"],
        "comparisons": comparisons,
        "frontierScenarios": [
            scenario("same-guarantee-weaker-projection-tie", ["log-loss-lane", "log-loss-finite"]),
            scenario("binary64-collision", ["binary64-low", "binary64-high"]),
            scenario("unattained-exponent", ["log-loss-lane", "omega-3_2", "omega-2"]),
            scenario("without-global", [c["id"] for c in eligible if c["id"] != "global"]),
        ],
    }


def encoded():
    return (json.dumps(generate(), ensure_ascii=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = encoded()
    if args.check:
        if FIXTURE.read_bytes() != content:
            raise SystemExit("Frontier fixture drift; regenerate with the pinned verifier")
        print("Frontier fixtures reproduce exactly")
    else:
        FIXTURE.write_bytes(content)
        print(f"Wrote {FIXTURE.relative_to(ROOT)} ({len(content)} bytes)")


if __name__ == "__main__":
    main()
