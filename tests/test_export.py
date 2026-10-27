"""Synthetic sealed phase chains test transport only. No provider, credentials,
mathematical review, or claimed discovery. Production-shaped synthetic records
are built ONLY inside tests to exercise the trusted-channel consumer; hashes
are integrity checks, not proof of real reviewer origin.
"""

import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from llm_verifier import pipeline
from llm_verifier.artifacts import seal
from llm_verifier.errors import VerifierError
from llm_verifier.files import write_json
from llm_verifier.jsonutil import canonical, digest
from llm_verifier.runner import DockerRunner, output_record

from export_metrics import export, projection
from finalize_score import main as finalize_main
from expressions import (CONVENTION, VERSION, claims_for, compare, polynomial,
                         ratio)
from tests.synthetic import synthetic_dossier
from tests.test_expressions import GLOBAL, GLOBAL_PROJECTION, finite, lane

ROOT = Path(__file__).resolve().parents[1]


def reseal(record):
    return seal(record["stage"], record["bindings"], record["data"], record["parents"], record["simulation"])


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "submission").mkdir()
        shutil.copyfile(ROOT / "verifier.json", self.root / "verifier.json")
        shutil.copyfile(ROOT / "rules.md", self.root / "rules.md")
        shutil.copytree(ROOT / "checkers", self.root / "checkers")
        self.bound = {
            "version": VERSION, "convention": CONVENTION,
            "family": "Synthetic declaration for transport tests, not a claimed circuit construction.",
            "theorem": "SYNTHETIC: no mathematical acceptance is asserted by this test.",
            "guarantee": polynomial(1, 1), "projection": finite(),
        }
        self.reset()

    def reset(self):
        self.manifest = {"version": 1, "title": "Synthetic test only", "claims": claims_for(self.bound),
                         "experiments": [], "checks": []}

    def dossier(self, **options):
        return synthetic_dossier(self.root, self.bound, self.manifest, **options)

    def test_finite_full_record_and_source_binding(self):
        d = self.dossier()
        result = projection(d, expected_intake=d["data"]["phases"][0]["sha256"])
        self.assertEqual(result["scores"], [1.0, 1.0])
        metric = result["metrics"]["circuitbreaker"]
        self.assertEqual(metric["bound"], self.bound)
        self.assertEqual(metric["statements"][0], metric["statements"][1])
        self.assertNotEqual(metric["score_records"][0]["source_sha256"], metric["score_records"][1]["source_sha256"])
        self.assertEqual(metric["evidence"], {p: d["data"]["phases"][0]["data"]["files"][p]["sha256"] for p in ["bound.json", "proof.md"]})
        self.assertEqual(d, self.dossier())

    def test_endpoint_metrics_survive_null_native_scores(self):
        for g, p in [(lane(2), lane(2)), (GLOBAL, GLOBAL_PROJECTION)]:
            self.bound.update(guarantee=g, projection=p); self.reset()
            result = projection(self.dossier())
            self.assertIsNone(result["scores"])
            self.assertEqual(result["metrics"]["circuitbreaker"]["bound"]["guarantee"], g)
            self.assertEqual(result["metrics"]["circuitbreaker"]["score_records"][0]["record"], "positive_infinity")

    def test_same_full_expression_has_equal_strength_for_weaker_projections_and_families(self):
        self.bound.update(guarantee=polynomial(1, 2, -1), projection=lane(ratio_to_fraction(3, 2))); self.reset()
        a = projection(self.dossier())["metrics"]["circuitbreaker"]["bound"]
        self.bound.update(family="Another synthetic eligible-family declaration.", projection=finite(100, ratio_to_fraction(4, 3))); self.reset()
        b = projection(self.dossier())["metrics"]["circuitbreaker"]["bound"]
        self.assertEqual(compare(a["guarantee"], b["guarantee"]), "equivalent")
        self.assertNotEqual(a["family"], b["family"])

    def test_string_lemmas_and_zero_bound_never_export_main_metrics(self):
        for lemma in ["A relevant exact LP arithmetic lemma.", "C_B2(f_n)>=0; this buys no exponent."]:
            self.manifest["claims"] = [lemma]
            result = projection(self.dossier())
            self.assertIsNone(result["scores"])
            self.assertNotIn("metrics", result)  # Even with an unused bound.json.

    def test_invalid_main_claim_exports_nothing(self):
        result = projection(self.dossier(valid=False))
        self.assertFalse(result["valid"])
        self.assertIsNone(result["scores"])
        self.assertNotIn("metrics", result)

    def test_independent_auxiliary_scores_do_not_unlock_metrics(self):
        for statement in ["An LP lemma is true.", "For every p there exists a different family L_p."]:
            self.reset()
            self.manifest["claims"][1]["statement"] = statement
            with self.assertRaisesRegex(VerifierError, "joint statement"):
                projection(self.dossier())

    def test_global_family_swap_in_either_claim_is_rejected(self):
        self.bound.update(guarantee=GLOBAL, projection=GLOBAL_PROJECTION); self.reset()
        for index in [0, 1]:
            self.reset(); self.manifest["claims"][index]["statement"] = "For every p choose a new NP language L_p."
            with self.assertRaisesRegex(VerifierError, "joint statement"):
                projection(self.dossier())

    def test_joint_same_family_scientific_rejection(self):
        self.bound.update(guarantee=GLOBAL, projection=GLOBAL_PROJECTION,
                          theorem="For each exponent choose a different family; this does not establish a global theorem.")
        self.reset()
        # Only a reviewer can evaluate conflicting mathematical prose. A known
        # rejection is preserved; the exporter is not a proof checker.
        self.assertNotIn("metrics", projection(self.dossier(valid=False)))

    def test_mismatched_evidence_scores_and_projection(self):
        mutations = [lambda: self.manifest["claims"][1].update(evidence=["proof.md"]),
                     lambda: self.manifest["claims"][0].update(score=ratio(2)),
                     lambda: self.bound.update(projection=finite(2)),
                     lambda: self.bound["guarantee"].update(coefficient=ratio(0)),
                     lambda: self.bound["guarantee"].update(exponent=ratio(1, 2)),
                     lambda: self.bound.update(dominates_all=True)]
        for mutate in mutations:
            self.setUp()
            mutate()
            with self.assertRaises(VerifierError):
                projection(self.dossier())

    def test_simulation_cannot_export(self):
        with self.assertRaisesRegex(VerifierError, "simulation"):
            projection(self.dossier(simulation=True))

    def test_tampered_resealed_dossier_is_not_a_new_review(self):
        d = self.dossier()
        for field in ["scores", "score_records", "bound_file", "verifier_identity"]:
            changed = copy.deepcopy(d)
            if field == "scores": changed["data"]["yukon_projection"]["scores"] = [2, 1]
            if field == "score_records": changed["data"]["score_records"][0]["record"] = ratio(2)
            if field == "bound_file": changed["data"]["phases"][0]["data"]["files"]["bound.json"]["text"] = "{}"
            if field == "verifier_identity": changed["data"]["phases"][0]["data"]["implementation"]["content_sha256"] = "0" * 64
            with self.subTest(field=field), self.assertRaises(VerifierError):
                projection(reseal(changed))
        with self.assertRaisesRegex(VerifierError, "intake_mismatch"):
            projection(d, expected_intake="0" * 64)

    def test_changing_evidence_invalidates_source_identity(self):
        d = self.dossier()
        changed = copy.deepcopy(d)
        changed["data"]["phases"][0]["data"]["files"]["proof.md"]["text"] = "Different proof."
        with self.assertRaises(VerifierError): projection(reseal(changed))

    def test_output_is_cleared_on_failed_export(self):
        d = self.dossier(simulation=True)
        source, output = self.root / "dossier.json", self.root / "score.json"
        write_json(source, d); write_json(output, {"stale": True})
        with self.assertRaises(VerifierError): export(source, output)
        self.assertFalse(output.exists())

    def test_selected_final_review_controls_five_phase_export(self):
        base = self.dossier()
        intake, initial, first = copy.deepcopy(base["data"]["phases"])
        fixture = json.loads((ROOT / "checker-tests.json").read_text())[0]
        first["data"]["decision"] = {"kind": "followup", "judgment": None, "requests": [{
            "check": "lp-certificate", "inputs": fixture["inputs"], "reason": "Synthetic LP request.",
        }]}
        first = reseal(first)
        output = output_record("completed", 0, canonical({**fixture["expected"], "message": "Synthetic fixture"}), b"")
        with patch.object(DockerRunner, "run", return_value=output):
            followup = pipeline.execute(intake, first, initial)
        messages, schema, _ = pipeline.review_material(intake, initial, first, followup)
        for valid in (True, False):
            decision = copy.deepcopy(base["data"]["phases"][2]["data"]["decision"])
            decision["judgment"]["valid"] = valid
            provenance = {**first["data"]["provenance"], "prompt_sha256": digest(messages), "schema_sha256": digest(schema)}
            final = seal("final_review", intake["bindings"], {"decision": decision, "provenance": provenance},
                         [p["sha256"] for p in [intake, initial, first, followup]])
            record = pipeline.dossier(intake, initial, first, followup, final)
            result = projection(record)
            self.assertEqual(result["valid"], valid)
            self.assertEqual("metrics" in result, valid)
            if valid:
                self.assertEqual(result["metrics"]["circuitbreaker"]["selected_review_sha256"], final["sha256"])
            altered = copy.deepcopy(record)
            altered["data"]["selected_review_sha256"] = first["sha256"]
            with self.assertRaises(VerifierError): projection(reseal(altered))

    def test_adversarial_configurations_fail_closed(self):
        original = json.loads((ROOT / "verifier.json").read_text())
        mutations = [lambda c: c.update(id="other-challenge"),
                     lambda c: c.update(review_standard="other-review"),
                     lambda c: c.update(allow_unscored=False),
                     lambda c: c["objectives"][0].update(name="other_coefficient"),
                     lambda c: c["objectives"][1].update(score_source={"claim": 0}),
                     lambda c: c["objectives"][0].update(score_minimum=ratio(-1)),
                     lambda c: c["objectives"][1].update(direction="-")]
        for mutate in mutations:
            config = copy.deepcopy(original); mutate(config)
            write_json(self.root / "verifier.json", config)
            with self.assertRaisesRegex(VerifierError, "objective_contract"):
                projection(self.dossier())

    def test_native_float_tie_keeps_distinct_exact_guarantees(self):
        self.bound.update(guarantee=polynomial(2**54, 1), projection=finite(2**54, 1)); self.reset()
        a = projection(self.dossier())
        self.bound.update(guarantee=polynomial(2**54 + 1, 1), projection=finite(2**54 + 1, 1)); self.reset()
        b = projection(self.dossier())
        self.assertEqual(a["scores"], b["scores"])
        self.assertEqual(compare(b["metrics"]["circuitbreaker"]["bound"]["guarantee"],
                                 a["metrics"]["circuitbreaker"]["bound"]["guarantee"]), "stronger")

    def test_output_path_aliases_preserve_inputs(self):
        source = self.root / "dossier.json"
        write_json(source, self.dossier())
        original = source.read_bytes()
        (self.root / "nested").mkdir()
        alias = self.root / "nested/../dossier.json"
        with self.assertRaisesRegex(VerifierError, "distinct_input_output"):
            export(source, alias)
        self.assertEqual(source.read_bytes(), original)
        hardlink = self.root / "hardlink.json"
        os.link(source, hardlink)
        with self.assertRaisesRegex(VerifierError, "distinct_input_output"):
            export(source, hardlink)
        self.assertEqual(source.read_bytes(), original)
        self.assertEqual(hardlink.read_bytes(), original)
        for score, dossier in [(alias, self.root / "other.json"),
                               (hardlink, self.root / "other.json"),
                               (self.root / "new.json", self.root / "nested/../new.json")]:
            args = ["finalize_score.py", "--intake", str(source), "--initial", str(source),
                    "--initial-review", str(source), "--score", str(score), "--dossier", str(dossier)]
            with patch("sys.argv", args), self.assertRaisesRegex(VerifierError, "distinct_input_output"):
                finalize_main()
            self.assertEqual(source.read_bytes(), original)

    def test_local_benchmark_adapter_emits_metrics(self):
        d = self.dossier()
        args = []
        for flag, phase in zip(["intake", "initial", "initial-review"], d["data"]["phases"]):
            path = self.root / (flag + ".json"); write_json(path, phase)
            args += ["--" + flag, str(path)]
        out, dossier = self.root / "score.json", self.root / "final.json"
        subprocess.run([sys.executable, str(ROOT / "finalize_score.py"), *args,
                        "--score", str(out), "--dossier", str(dossier)], check=True)
        self.assertIn("circuitbreaker_frontier", json.loads(out.read_text())["metrics"])


def ratio_to_fraction(n, d):
    from fractions import Fraction
    return Fraction(n, d)
