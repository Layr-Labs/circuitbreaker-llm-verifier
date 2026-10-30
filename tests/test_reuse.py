"""Accepted-result transport with synthetic judges, never scientific acceptance."""

from fractions import Fraction
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from llm_verifier import pipeline
from llm_verifier.accepted import PlatformTransport
from llm_verifier.artifacts import intake_snapshot, verify_intake
from llm_verifier.config import load_config
from llm_verifier.errors import VerifierError
from llm_verifier.files import write_json
from llm_verifier.jsonutil import canonical, sha
from llm_verifier.prior import FORMAT
from llm_verifier.provider import HttpTransport
from llm_verifier.runner import DockerRunner, output_record

from export_metrics import export, projection
from expressions import CONVENTION, VERSION, claims_for, polynomial
from tests.synthetic import synthetic_dossier, synthetic_review
from tests.test_expressions import finite

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "11111111-1111-4111-8111-111111111111"
ON = {"version": 1, "submissionReferencesEnabled": True}
OFF = {**ON, "submissionReferencesEnabled": False}
CONTEXT = {"version": 1, "api": "https://yukon.invalid", "audience": "yukon:accepted-results:synthetic",
           "submission_id": "22222222-2222-4222-8222-222222222222", "benchmark_id": "consumer-benchmark",
           "challenge_id": None, "evaluation_id": "consumer-evaluation", "job_id": "consumer-job"}


class ReuseTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def package(self, name, coefficient, prior=False, *, admission=None):
        root = self.root / name
        (root / "submission").mkdir(parents=True)
        for filename in ("verifier.json", "rules.md", "benchmark.json"):
            shutil.copyfile(ROOT / filename, root / filename)
        shutil.copytree(ROOT / "checkers", root / "checkers")
        if admission is not None:
            benchmark = json.loads((root / "benchmark.json").read_text())
            benchmark["tracks"][0]["submissionAdmission"] = admission
            write_json(root / "benchmark.json", benchmark)
        bound = {"version": VERSION, "convention": CONVENTION,
                 "family": "SYNTHETIC reuse fixture; no circuit theorem.",
                 "theorem": "SYNTHETIC current declaration, not a scientific judgment.",
                 "guarantee": polynomial(coefficient, 1), "projection": finite(coefficient)}
        manifest = {"version": 1, "title": "SYNTHETIC reuse integration",
                    "claims": claims_for(bound), "experiments": [], "checks": []}
        if prior:
            manifest["prior_results"] = [{"submission_id": SOURCE_ID}]
        write_json(root / "submission/bound.json", bound)
        write_json(root / "submission/submission.json", manifest)
        (root / "submission/proof.md").write_text("SYNTHETIC reuse evidence; no circuit theorem.\n")
        return root, bound, manifest

    def evaluation(self, root, admission):
        return {"policies": json.loads((root / "verifier.json").read_text())["policies"],
                "fields": {"scores": True, "validity": True}, "admission": admission}

    def test_enabled_prior_intake_current_checks_and_selected_audit_export(self):
        source, bound, manifest = self.package("source", 9)
        source_record = synthetic_dossier(source, bound, manifest, valid=False, audit_valid=True)
        source_bundle = pipeline.export_results(source_record)
        self.assertEqual(projection(source_record)["scores"], [9, 1])
        self.assertEqual(source_record["data"]["phases"][2]["data"]["authority"]["judge"], "escalation")
        raw_bundle = json.dumps(source_bundle, indent=2) + "\n"
        raw_bytes = raw_bundle.encode()
        snapshot = {"version": 1, "snapshot_id": "synthetic-frozen-snapshot",
                    "consumer": {k: CONTEXT[k] for k in ("submission_id", "benchmark_id", "challenge_id", "evaluation_id", "job_id")},
                    "references": [{"submission_id": SOURCE_ID, "bundle_sha256": sha(raw_bytes)}],
                    "sources": [{"submission_id": SOURCE_ID, "benchmark_id": "source-benchmark", "challenge_id": None,
                                 "evaluation_id": "source-evaluation", "validity_result_id": "source-validity",
                                 "verifier_evaluation_id": source_bundle["source"]["evaluation_id"], "accepted": True,
                                 "submitted_at": "2020-01-01T00:00:00Z", "completed_at": "2026-10-29T00:00:00Z",
                                 "bundle": {"format": FORMAT, "sha256": sha(raw_bytes), "bytes": len(raw_bytes), "text": raw_bundle}}],
                    "ancestors": []}
        consumer, _, manifest = self.package("consumer", Fraction(7, 2), prior=True, admission=ON)
        # The current checker still runs; prior score 9 cannot replace current 3.5.
        manifest["checks"] = json.loads((ROOT / "submission/submission.json").read_text())["checks"]
        write_json(consumer / "submission/submission.json", manifest)
        original_manifest = (consumer / "submission/submission.json").read_bytes()
        evaluation = self.evaluation(consumer, ON)
        effective = load_config(consumer, "verifier.json", "submission", evaluation)
        with patch.dict("os.environ", {"ACTIONS_ID_TOKEN_REQUEST_URL": "https://actions.invalid/token",
                                      "ACTIONS_ID_TOKEN_REQUEST_TOKEN": "synthetic-never-network"}), \
                patch.object(PlatformTransport, "request", return_value=snapshot) as resolver:
            intake = intake_snapshot(effective, consumer / "submission", accepted_results_context=CONTEXT)
        resolver.assert_called_once()
        self.assertEqual(json.loads(resolver.call_args.args[1])["references"], [{"submission_id": SOURCE_ID}])
        self.assertEqual(intake["data"]["prior_snapshot"]["sources"][0]["bundle"]["text"], raw_bundle)
        self.assertEqual(source_bundle["source"]["manifest"], source_record["data"]["phases"][0]["data"]["files"]["submission.json"])
        before = canonical(intake)
        current = json.loads((consumer / "benchmark.json").read_text())
        current["tracks"][0]["submissionAdmission"] = OFF
        write_json(consumer / "benchmark.json", current)
        evaluation["admission"] = OFF
        fixture = json.loads((ROOT / "checker-tests.json").read_text())[0]
        output = output_record("completed", 0, canonical({**fixture["expected"], "message": "Synthetic current checker"}), b"")
        with patch.dict("os.environ", {"YUKON_EVALUATION": json.dumps(evaluation)}), \
                patch.object(PlatformTransport, "request", side_effect=AssertionError("No post-intake resolution")), \
                patch.object(DockerRunner, "run", return_value=output) as checker:
            verify_intake(intake)
            initial = pipeline.execute(intake)
            checker.assert_called_once()
            review = synthetic_review(intake, initial, valid=False, audit_valid=True, reference="prior:" + SOURCE_ID)
            record = pipeline.finalize(intake, initial, review, consumer / "native-score.json",
                                       consumer / "dossier.json", results_path=consumer / "accepted-results.json")
            with patch.object(HttpTransport, "request", side_effect=AssertionError("No inference during export")):
                export(consumer / "dossier.json", consumer / "score.json")
                accepted = pipeline.export_results(record)
        self.assertEqual(canonical(intake), before)
        self.assertEqual(intake["data"]["effective"]["evaluation"]["admission"], ON)
        self.assertEqual(json.loads((consumer / "score.json").read_text())["scores"], [3.5, 1])
        self.assertEqual(accepted["ancestors"], [{"submission_id": SOURCE_ID}])
        self.assertEqual(accepted, json.loads((consumer / "accepted-results.json").read_text()))
        self.assertEqual(accepted["source"]["manifest"]["text"].encode(), original_manifest)
        self.assertEqual(accepted["source"]["review_sha256"], review["sha256"])
        self.assertEqual(record["data"]["selected_review_sha256"], review["sha256"])
        self.assertEqual(review["data"]["authority"]["judge"], "escalation")
        self.assertEqual(review["data"]["decision"], review["data"]["escalation"]["review"]["decision"])
        self.assertEqual(json.loads((consumer / "accepted-results-receipt.json").read_text()),
                         {"version": 1, "snapshot_id": "synthetic-frozen-snapshot"})

    def test_enabled_manifest_cannot_upgrade_original_disabled_admission(self):
        root, _, _ = self.package("disabled-original", Fraction(7, 2), prior=True, admission=ON)
        self.assertEqual(json.loads((root / "benchmark.json").read_text())["tracks"][0]["submissionAdmission"], ON)
        for admission in (OFF, None, "missing"):
            original = self.evaluation(root, admission)
            if admission == "missing":
                original.pop("admission")
            effective = load_config(root, "verifier.json", "submission", original)
            with self.subTest(original_admission=admission), \
                    patch.dict("os.environ", {"YUKON_EVALUATION": json.dumps(self.evaluation(root, ON))}), \
                    patch.object(PlatformTransport, "request") as resolver, \
                    patch.object(HttpTransport, "request") as model, patch.object(DockerRunner, "run") as checker:
                with self.assertRaisesRegex(VerifierError, "^submission_references_disabled$"):
                    intake_snapshot(effective, root / "submission", accepted_results_context=CONTEXT)
                resolver.assert_not_called()
                model.assert_not_called()
                checker.assert_not_called()
        with self.assertRaisesRegex(VerifierError, "prior_platform_context_required"):
            intake_snapshot(load_config(root, "verifier.json", "submission", self.evaluation(root, ON)), root / "submission")
