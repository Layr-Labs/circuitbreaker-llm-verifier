"""Offline generic verifier/adapter integration with synthetic accepted sources."""

from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

try:
    from llm_verifier import derivation
except ImportError:
    derivation = None  # The separate pinned feature CI job requires this module.
from llm_verifier.errors import VerifierError
from llm_verifier.files import write_json
from llm_verifier.jsonutil import sha
from llm_verifier.simulation import ScriptedProvider

from score_derivation import SCORING_POLICY, projection
from tests.test_refinements import synthetic_origin

ROOT = Path(__file__).resolve().parents[1]


def uid(number):
    return f"00000000-0000-4000-8000-{number:012x}"


def archive(root, artifact_id, name, path):
    raw = path.read_bytes()
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as target:
        target.writestr(path.name, raw)
    content = output.getvalue()
    root.mkdir(exist_ok=True)
    (root / (artifact_id + ".zip")).write_bytes(content)
    return {"repository": "synthetic/circuitbreaker", "run_id": "123", "run_attempt": "1",
            "artifact_id": artifact_id, "artifact_name": name, "archive_sha256": sha(content),
            "path": path.name, "file_sha256": sha(raw)}


@unittest.skipIf(derivation is None, "requires the score-derivation verifier; exercised by its dedicated CI job")
class GenericAdapterIntegration(unittest.TestCase):
    def setUp(self):
        global build
        from tools.prepare_score_derivation_plan import build
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source, bundle, self.row, self.target = synthetic_origin(self.root)
        packet = self.row["evaluation"]
        self.row.update(id=uid(1), benchmarkId=uid(2))
        packet.update(id=uid(3), publishedAt="2026-01-01T00:00:00Z")
        for offset, name in ((10, "scores"), (20, "validity")):
            packet[name]["evidence"].update(resultId=uid(offset), evaluatorId=uid(offset + 1),
                                           policyVersionId=uid(offset + 2), jobId=uid(offset + 3))
        self.original_score = json.loads((bundle / "circuitbreaker-score/.yukon/score.json").read_text())
        self.old = json.loads((bundle / "llmv-intake/intake.json").read_text())
        self.archives = self.root / "archives"
        self.index = {"version": "circuitbreaker-derivation-artifacts-v1", "records": [{
            "submission_id": self.row["id"], "source_files_sha256": self.old["bindings"]["submission"],
            "artifacts": {
                "intake": archive(self.archives, "101", "llmv-intake", bundle / "llmv-intake/intake.json"),
                "accepted_results": archive(self.archives, "102", "verifier-accepted-results",
                                            bundle / "verifier-accepted-results/accepted-results.json"),
            },
        }]}
        self.config_raw = (ROOT / "derivation.json").read_bytes()
        self.plan, self.report = build({"submissions": [self.row]}, self.index, self.config_raw, self.archives)
        (self.source / "derivation.json").write_bytes(self.config_raw)
        write_json(self.source / "reorg/score-derivation-plan.json", self.plan)
        retained = {**deepcopy(self.plan["entries"][0]["origin"]),
                    "destination": {"job_id": uid(30), "evaluation_id": uid(31),
                                    "packet_id": uid(32), "source_ref": "synthetic-new-source"}}
        self.evaluation = {"policies": {"scoring": deepcopy(SCORING_POLICY),
                                       "validity": {"version": "fixed-family-full-b2-v1"}},
                           "fields": {"scores": True, "validity": False}, "retained_result": retained}
        self.context = {"kind": "github-actions", "candidate_revision": "d" * 40,
                        "verifier_revision": "e" * 40, "run_id": "456:1",
                        "repository": "synthetic/circuitbreaker", "score_path": ".yukon/score.json"}

    def intake(self):
        return derivation.prepare(self.source, "derivation.json", "submission", self.evaluation,
                                  self.context, self.archives)

    def decision(self, intake, status="derived"):
        statement = derivation.sources(intake)["claim:0"]
        return {"status": status, "payload": {"target_bound": self.target} if status == "derived" else None,
                "reason": "Synthetic transcription fixture, not a mathematical judgment.",
                "source_refs": [{"source": "claim:0", "sha256": sha(statement.encode()), "start_line": 1,
                                 "end_line": len(statement.split("\n")), "quote": statement}]}

    def simulate(self, status="derived"):
        intake = self.intake()
        with patch("socket.socket", side_effect=AssertionError("offline test attempted network")):
            result = derivation.review(intake, provider=ScriptedProvider(
                intake["data"]["config"]["provider"], self.decision(intake, status)))
        return intake, result

    def test_complete_sealed_derivation_preserves_prior_rich_claim_and_blocks_simulated_publication(self):
        intake, result = self.simulate()
        preview = projection(result, expected_intake=intake["sha256"], preview=True)
        self.assertTrue(preview["simulation"])
        self.assertFalse(preview["publishable"])
        self.assertNotIn("scores", preview)
        score = preview["score_preview"]
        self.assertEqual(score["metrics"]["circuitbreaker"], self.original_score["metrics"]["circuitbreaker"])
        self.assertEqual(score["scores"], self.original_score["scores"])
        self.assertEqual(score["derived_from"], self.evaluation["retained_result"])
        self.assertEqual(score["metrics"]["circuitbreaker_frontier"]["resources"]["circuit_size"]["bound"],
                         self.target["guarantee"])
        self.assertNotIn("valid", score)
        with self.assertRaisesRegex(VerifierError, "simulation"):
            projection(result, expected_intake=intake["sha256"])

    def test_unsupported_derivation_preserves_original_v1_without_inventing_a_remainder(self):
        intake, result = self.simulate("unsupported")
        score = projection(result, expected_intake=intake["sha256"], preview=True)["score_preview"]
        self.assertEqual(score["metrics"]["circuitbreaker_frontier"]["version"], "circuitbreaker-frontier-v1")
        self.assertEqual(score["metrics"]["circuitbreaker"], self.original_score["metrics"]["circuitbreaker"])

    def test_plan_covers_retained_rows_and_reports_missing_validity_without_model_calls(self):
        rejected = deepcopy(self.row)
        rejected.update(id=uid(40), status="rejected")
        rejected["evaluation"]["validity"]["valid"] = False
        missing = {"id": uid(41), "status": "failed", "evaluation": None}
        cancelled = {"id": uid(42), "status": "cancelled"}
        index = deepcopy(self.index)
        index["records"].append({"submission_id": rejected["id"],
                                 "source_files_sha256": self.old["bindings"]["submission"], "artifacts": None})
        plan, report = build({"submissions": [self.row, rejected, missing, cancelled]}, index,
                             self.config_raw, self.archives)
        self.assertEqual((report["score_only"], report["derive"], report["not_applicable"]), (2, 1, 1))
        self.assertEqual(report["ordinary_submission_ids"], [missing["id"]])
        self.assertEqual(report["excluded_submission_ids"], [cancelled["id"]])
        self.assertEqual(report["model_calls"], 0)
        self.assertFalse(report["published"])
        write_json(self.source / "reorg/score-derivation-plan.json", plan)
        self.evaluation["retained_result"] = {
            **deepcopy(plan["entries"][1]["origin"]),
            "destination": self.evaluation["retained_result"]["destination"],
        }
        intake = self.intake()
        with patch("llm_verifier.provider.OpenAI", side_effect=AssertionError("unexpected provider")):
            result = derivation.review(intake)
        score = projection(result, expected_intake=intake["sha256"])
        self.assertIsNone(score["scores"])
        self.assertEqual(set(score["metrics"]), {"circuitbreaker_derivation"})
        self.assertNotIn("valid", score)

    def test_missing_extra_and_duplicate_artifact_records_reject(self):
        for records in ([], self.index["records"] * 2,
                        [*self.index["records"], {**self.index["records"][0], "submission_id": uid(90)}]):
            with self.subTest(records=len(records)), self.assertRaises(ValueError):
                build({"submissions": [self.row]}, {**self.index, "records": records},
                      self.config_raw, self.archives)

    def test_changed_candidate_prior_selection_archive_and_citation_fail_closed(self):
        self.evaluation["retained_result"]["validity"]["result_id"] = uid(99)
        with self.assertRaises(VerifierError):
            self.intake()
        self.evaluation["retained_result"]["validity"]["result_id"] = uid(20)
        proof = self.source / "submission/proof.md"
        original = proof.read_bytes()
        proof.write_bytes(original + b"\nchanged")
        with self.assertRaises(VerifierError):
            self.intake()
        proof.write_bytes(original)
        intake = self.intake()
        decision = self.decision(intake)
        decision["source_refs"][0]["quote"] += "invented"
        with self.assertRaises(VerifierError):
            derivation.review(intake, provider=ScriptedProvider(intake["data"]["config"]["provider"], decision))
        archive_path = self.archives / "101.zip"
        archive_path.write_bytes(archive_path.read_bytes() + b"changed")
        with self.assertRaises(VerifierError):
            self.intake()

    def test_materialized_profile_preserves_default_and_routes_explicit_fields(self):
        from llm_verifier.config import load_config
        from tools.derivation_profile import materialize, read_profile, regenerate
        from tools.generate_hosted import check_generated
        original_config = json.loads((ROOT / "verifier.json").read_text())
        original_workflow = (ROOT / ".github/workflows/benchmark.yml").read_bytes()
        output = self.root / "profile"
        profile = materialize(ROOT, output, {"submissions": [self.row]}, self.index, self.archives, "e" * 40)
        self.assertEqual(read_profile(output), profile)
        self.assertEqual(profile["retained_submission_ids"], [self.row["id"]])
        config = json.loads((output / "verifier.json").read_text())
        for field in ("review_standard", "judging", "acceptance_policy", "limits", "provider"):
            self.assertEqual(config[field], original_config[field])
        self.assertEqual(config["policies"]["validity"], original_config["policies"]["validity"])
        baseline = deepcopy(self.evaluation)
        baseline.pop("retained_result")
        baseline["fields"]["validity"] = True
        effective = load_config(output, "verifier.json", "submission", baseline)
        self.assertEqual(effective["config"]["policies"]["scoring"], SCORING_POLICY)
        from export_metrics import validate_contract
        self.assertFalse(validate_contract(effective["config"]))
        workflow = (output / ".github/workflows/benchmark.yml").read_text()
        self.assertIn("if: needs.route.outputs.route == 'ordinary'", workflow)
        self.assertIn("if: needs.route.outputs.route == 'derive'", workflow)
        self.assertIn("/derive-score.yml@" + "e" * 40, workflow)
        self.assertIn("/verify.yml@" + "e" * 40, workflow)
        self.assertIn("python3 select_scoring_route.py", workflow)
        self.assertNotIn("@VERIFIER_REVISION@", workflow)
        self.assertEqual(json.loads((ROOT / "verifier.json").read_text()), original_config)
        self.assertEqual((ROOT / ".github/workflows/benchmark.yml").read_bytes(), original_workflow)
        from tools.check_derivation_profile import validate
        self.assertEqual(validate(output)["retained_records"], 1)
        report_path = output / "reorg/score-derivation-report.json"
        report_bytes = report_path.read_bytes()
        metadata_path = output / "score-derivation-profile.json"
        metadata_bytes = metadata_path.read_bytes()
        report = json.loads(report_bytes)
        report["plan_sha256"] = "0" * 64
        write_json(report_path, report)
        metadata = json.loads(metadata_bytes)
        metadata["files"]["reorg/score-derivation-report.json"] = sha(report_path.read_bytes())
        write_json(metadata_path, metadata)
        with self.assertRaisesRegex(ValueError, "report/plan"):
            validate(output)
        report_path.write_bytes(report_bytes)
        metadata_path.write_bytes(metadata_bytes)
        self.assertFalse((output / ".github/workflows/score-derivation-tests.yml").exists())
        profile_checks = (output / ".github/workflows/checker-tests.yml").read_text()
        self.assertIn("/test-checkers.yml@" + "e" * 40, profile_checks)
        self.assertIn("tools/check_derivation_profile.py", profile_checks)
        self.assertIn("-m unittest tests.test_score_derivation -v", profile_checks)
        self.assertNotIn("-m unittest -v", profile_checks)
        self.assertNotIn("generate_remainder_fixtures.py", profile_checks)
        bad_workflow = output / ".github/workflows/checker-tests.yml"
        bad_workflow.write_text(profile_checks + "# drift\n")
        with self.assertRaises(ValueError):
            validate(output)
        regenerate(output)
        self.assertEqual((output / ".github/workflows/benchmark.yml").read_text(), workflow)
        check_generated(output)
        (output / "reorg/score-derivation-plan.json").write_text('{"entries": []}')
        with self.assertRaises(ValueError):
            read_profile(output)

    def test_incomplete_profile_plan_is_rejected_before_any_source_exists(self):
        from tools.derivation_profile import materialize
        for snapshot, index in (({"submissions": []}, {**self.index, "records": []}),
                                ({"submissions": [self.row]}, {**self.index, "records": []})):
            output = self.root / "incomplete-profile"
            with self.assertRaises(ValueError):
                materialize(ROOT, output, snapshot, index, self.archives, "e" * 40)
            self.assertFalse(output.exists())

    def test_same_run_export_authenticates_transfer_and_rejects_cross_run_or_changed_provenance(self):
        from export_hosted_derivation import export, NAMES
        from tools.derivation_profile import materialize
        # A retained-invalid source has a deterministic, non-simulated result;
        # no fabricated acceptance or fake model response is needed to test delivery.
        self.row["status"] = "rejected"
        self.row["evaluation"]["validity"]["valid"] = False
        self.index["records"][0]["artifacts"] = None
        output = self.root / "invalid-profile"
        materialize(ROOT, output, {"submissions": [self.row]}, self.index, self.archives, "e" * 40)
        plan = json.loads((output / "reorg/score-derivation-plan.json").read_text())
        self.evaluation["retained_result"] = {**plan["entries"][0]["origin"],
                                              "destination": self.evaluation["retained_result"]["destination"]}
        # Use the exact original submitted bytes, as the real Yukon job does.
        import shutil
        shutil.rmtree(output / "submission")
        shutil.copytree(self.source / "submission", output / "submission")
        intake = derivation.prepare(output, "derivation.json", "submission", self.evaluation,
                                    self.context, self.archives)
        with patch("llm_verifier.provider.OpenAI", side_effect=AssertionError("unexpected provider")):
            result = derivation.review(intake)
        transfer = self.root / "current-run"
        transfer.mkdir()
        wheel = transfer / "synthetic.whl"
        wheel.write_bytes(b"synthetic transfer fixture, never installed")
        intake_path, result_path = transfer / "intake.json", transfer / "derivation.json"
        write_json(intake_path, intake)
        write_json(result_path, result)
        descriptors = [{"id": str(501 + i), "name": name, "sha256": "f" * 64,
                        "run_id": "456", "run_attempt": "1"} for i, name in enumerate(NAMES)]
        receipt_path = transfer / "receipt.json"
        write_json(receipt_path, {"version": 1, "run_id": "456", "delivery_attempt": "1", "artifacts": [
            {"artifact": item, "path": name, "files": {path.name: sha(path.read_bytes())}}
            for item, name, path in zip(descriptors, ("wheel", "intake", "derivation"), (wheel, intake_path, result_path))]})
        env = {"GITHUB_RUN_ID": "456", "GITHUB_RUN_ATTEMPT": "1", "GITHUB_SHA": "d" * 40,
               "GITHUB_REPOSITORY": "synthetic/circuitbreaker", "LLMV_REVIEW_RUN_ID": "456",
               "LLMV_REVIEW_RUN_ATTEMPT": "1", "LLMV_EXPECTED_INTAKE": intake["sha256"]}
        env.update({key: json.dumps(item) for key, item in zip(
            ("LLMV_PACKAGE_ARTIFACT", "LLMV_INTAKE_ARTIFACT", "LLMV_DERIVATION_ARTIFACT"), descriptors)})
        score, lineage = transfer / "score.json", transfer / "lineage.json"
        export(intake_path, result_path, wheel, score, lineage, receipt_path, env=env, root=output)
        self.assertIsNone(json.loads(score.read_text())["scores"])
        self.assertNotIn("valid", json.loads(score.read_text()))
        self.assertEqual(json.loads(lineage.read_text())["derived_from"], self.evaluation["retained_result"])
        for key, value in (("GITHUB_RUN_ID", "999"), ("GITHUB_SHA", "f" * 40),
                           ("GITHUB_REPOSITORY", "wrong/repository"), ("LLMV_EXPECTED_INTAKE", "0" * 64)):
            with self.subTest(key=key), self.assertRaises(VerifierError):
                export(intake_path, result_path, wheel, score, lineage, receipt_path,
                       env={**env, key: value}, root=output)
            self.assertFalse(score.exists())
            self.assertFalse(lineage.exists())
        wheel.write_bytes(b"changed wheel")
        with self.assertRaises(VerifierError):
            export(intake_path, result_path, wheel, score, lineage, receipt_path, env=env, root=output)
