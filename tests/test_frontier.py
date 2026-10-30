"""Additive table wire, exact cross-language fixtures, and schema correspondence."""

import copy
from fractions import Fraction as F
import hashlib
import json
import re
from pathlib import Path
import shutil
import tempfile
import unittest

from export_metrics import projection, validate_contract
from expressions import CONVENTION, VERSION, claims_for, compare, polynomial, ratio
from llm_verifier.errors import VerifierError
from llm_verifier.files import relative
from llm_verifier.jsonutil import canonical
from tests.synthetic import synthetic_dossier
from tools.generate_hosted import VERIFIER_REVISION

ROOT = Path(__file__).resolve().parents[1]
PINNED_RICH_HASHES = {
    "finite": "9c30149d6222d2bfc65d742993e7c7af287e906867f39856386425ba415fd3c3",
    "lane": "44f0ac59180ad2194a27e30ab628fc5f0605265ddb290fa786e1a968d969fd5a",
    "global": "904098358937d647f902435ce5dd15a8e6ab108a7349d5aa3bbb8bc97446d974"
}


def sample_bound(kind):
    g, p = {
        "finite": (polynomial(F(31, 10), 1), {"kind": "finite", "coefficient": ratio(31, 10), "exponent": ratio(1)}),
        "lane": ({"kind": "lane_endpoint", "exponent": ratio(2)}, {"kind": "lane_endpoint", "exponent": ratio(2)}),
        "global": ({"kind": "superpolynomial"}, {"kind": "global_endpoint"}),
    }[kind]
    return {"version": VERSION, "convention": CONVENTION,
            "family": "SYNTHETIC table compatibility fixture; no claimed circuit construction.",
            "theorem": "SYNTHETIC scripted evidence; not a mathematical discovery.",
            "guarantee": g, "projection": p}


class FrontierTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        (self.root / "submission").mkdir()
        for filename in ("verifier.json", "rules.md"):
            shutil.copyfile(ROOT / filename, self.root / filename)
        shutil.copytree(ROOT / "checkers", self.root / "checkers")

    def dossier(self, kind="finite", *, valid=True, lemma=False, simulation=False):
        bound = sample_bound(kind)
        manifest = {"version": 1, "title": "SYNTHETIC table compatibility test",
                    "claims": ["SYNTHETIC unscored lemma."] if lemma else claims_for(bound),
                    "experiments": [], "checks": []}
        return synthetic_dossier(self.root, bound, manifest, valid=valid, simulation=simulation)

    def test_additive_wire_preserves_legacy_metric_and_exact_joint_binding(self):
        # Fresh verifier-revision fixture hashes, including new intake/review
        # provenance. The scientific/exporter implementation is unchanged.
        for kind in PINNED_RICH_HASHES:
            with self.subTest(kind=kind):
                dossier = self.dossier(kind)
                result = projection(dossier)
                rich, table = result["metrics"]["circuitbreaker"], result["metrics"]["circuitbreaker_frontier"]
                self.assertEqual(hashlib.sha256(canonical(rich)).hexdigest(), PINNED_RICH_HASHES[kind])
                self.assertEqual(set(table), {"version", "comparator", "scope", "model", "notation", "intake_sha256", "dossier_sha256", "resources"})
                self.assertEqual(set(table["resources"]), {"circuit_size"})
                resource = table["resources"]["circuit_size"]
                self.assertEqual(resource, {"statement": rich["statements"][0], "bound": rich["bound"]["guarantee"],
                                            "source_sha256": rich["score_records"][0]["source_sha256"], "evidence": rich["evidence"]})
                self.assertEqual(rich["statements"][0], rich["statements"][1])
                self.assertEqual(table["intake_sha256"], rich["intake_sha256"])
                self.assertEqual(table["dossier_sha256"], rich["dossier_sha256"])
                self.assertEqual(result["scores"], [3.1, 1.0] if kind == "finite" else None)

    def test_no_table_from_rejection_lemma_simulation_or_opt_in_protocol(self):
        for options in ({"valid": False}, {"lemma": True}):
            self.assertNotIn("metrics", projection(self.dossier(**options)))
        with self.assertRaisesRegex(VerifierError, "simulation"):
            projection(self.dossier(simulation=True))
        config = json.loads((ROOT / "verifier.json").read_text())
        self.assertNotIn("structured_claims", config)
        config["structured_claims"] = "atomic-v1"
        with self.assertRaisesRegex(VerifierError, "objective_contract"):
            validate_contract(config)

    def test_schema_and_page_fragment_match_exact_emitted_contract(self):
        schema = json.loads((ROOT / "schemas/frontier-v1.schema.json").read_text())
        tab = json.loads((ROOT / "examples/complexity-result-tab.json").read_text())["resultTabs"][0]
        table = projection(self.dossier())["metrics"]["circuitbreaker_frontier"]
        self.assertEqual(set(schema["required"]), set(table))
        self.assertEqual(set(schema["properties"]), set(table))
        self.assertFalse(schema["additionalProperties"])
        for field, config_key in [("version", "metricsVersion"), ("comparator", "comparator"), ("scope", "scope"), ("model", "model"), ("notation", "notation")]:
            self.assertEqual(schema["properties"][field]["const"], table[field])
            self.assertEqual(table[field], tab[config_key])
        self.assertEqual(tab["metricsPath"], "circuitbreaker_frontier")
        self.assertEqual(tab["columns"], [{"key": "circuit_size", "label": "Circuit size lower bound", "direction": "+"}])
        self.assertEqual(tab["defaultSort"], {"by": "metric", "key": "circuit_size", "direction": "desc"})
        resources = schema["properties"]["resources"]
        self.assertEqual(resources["required"], ["circuit_size"])
        resource = resources["properties"]["circuit_size"]
        self.assertEqual(set(resource["required"]), set(table["resources"]["circuit_size"]))
        self.assertFalse(resource["additionalProperties"])
        self.assertEqual(resource["properties"]["bound"], {"$ref": "bound-v1.schema.json#/properties/guarantee"})
        self.assertEqual(resource["properties"]["statement"]["maxLength"], 8000)

    def test_schema_evidence_paths_match_pinned_verifier(self):
        schema = json.loads((ROOT / "schemas/frontier-v1.schema.json").read_text())
        paths = schema["properties"]["resources"]["properties"]["circuit_size"]["properties"]["evidence"]["propertyNames"]
        examples = ["constructor", "prototype", "__proto__", "proof.md", "a" * 128, "a" * 129, "nested/proof.md", "proof with spaces.md",
                    "../proof.md", "a/./b", "a/.git/b", ".git", ".github/proof.md",
                    "proof.md\n", "/proof.md", "a/", "a\\b", "é.md",
                    "/".join(["a"] * 12), "/".join(["a"] * 13),
                    "/".join(["a" * 128] * 4)]
        for path in examples:
            matches_schema = len(path) <= paths["maxLength"] and re.search(paths["pattern"], path) is not None
            try:
                relative(path)
                matches_verifier = True
            except VerifierError:
                matches_verifier = False
            self.assertEqual(matches_schema, matches_verifier, repr(path))

    def test_shared_corpus_independent_order_and_frontier_membership(self):
        fixture = json.loads((ROOT / "tests/fixtures/circuitbreaker-frontier-v1.json").read_text())
        self.assertTrue(fixture["provenance"]["synthetic"])
        self.assertEqual(fixture["provenance"]["verifierRevision"], VERIFIER_REVISION)
        eligible = {c["id"]: c for c in fixture["cases"] if c["expectedEligible"]}
        self.assertEqual(len(eligible), 58)
        self.assertEqual(list(eligible), fixture["expectedEligibleIds"])
        self.assertEqual(len(fixture["comparisons"]), 3364)
        first = eligible["poly-1-0-0-1_2"]["output"]["metrics"]["circuitbreaker_frontier"]["resources"]["circuit_size"]
        self.assertEqual(set(first["evidence"]), {"bound.json", "proof.md", "constructor", "prototype", "__proto__"})
        second = eligible["poly-1-0-0-31_10"]["output"]["metrics"]["circuitbreaker_frontier"]["resources"]["circuit_size"]
        self.assertEqual(set(second["evidence"]), {"bound.json", "proof.md"})
        def guarantee(identifier):
            return eligible[identifier]["output"]["metrics"]["circuitbreaker_frontier"]["resources"]["circuit_size"]["bound"]
        for entry in fixture["comparisons"]:
            self.assertEqual(compare(guarantee(entry["left"]), guarantee(entry["right"])), entry["expected"], entry)
        for scenario in fixture["frontierScenarios"] + [{"caseIds": list(eligible), "expectedFrontierIds": fixture["expectedFrontierIds"]}]:
            frontier = [a for a in scenario["caseIds"] if not any(compare(guarantee(b), guarantee(a)) == "stronger" for b in scenario["caseIds"])]
            self.assertEqual(frontier, scenario["expectedFrontierIds"])
        self.assertEqual(eligible["binary64-low"]["output"]["scores"], eligible["binary64-high"]["output"]["scores"])
        self.assertNotEqual(guarantee("binary64-low"), guarantee("binary64-high"))
        self.assertIsNone(eligible["log-loss-lane"]["output"]["scores"])
        self.assertIsNotNone(eligible["log-loss-finite"]["output"]["scores"])
        for identifier in ["omega-1", "omega-3_2", "omega-2", "global"]:
            self.assertIsNone(eligible[identifier]["output"]["scores"])
        for case in fixture["cases"]:
            if case["id"] in {"rejected", "lemma"}:
                self.assertNotIn("metrics", case["output"])
            if case["id"] == "simulation":
                self.assertIsNone(case["output"])
            if case["provenanceKind"] != "exporter":
                self.assertFalse(case["expectedEligible"])

    def test_malformed_bounds_are_unknown_not_ties(self):
        fixture = json.loads((ROOT / "tests/fixtures/circuitbreaker-frontier-v1.json").read_text())
        cases = {c["id"]: c for c in fixture["cases"]}
        for identifier in ("zero-coefficient", "unsupported-exponent", "nonreduced-coefficient", "comparison-hint", "unsupported-kind"):
            bad = cases[identifier]["output"]["metrics"]["circuitbreaker_frontier"]["resources"]["circuit_size"]["bound"]
            self.assertEqual(compare(bad, bad), "unknown")
