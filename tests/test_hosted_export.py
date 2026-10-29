"""Deterministic delivery tests; synthetic records make no scientific claim."""

import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from llm_verifier.artifacts import intake_snapshot
from llm_verifier.errors import InfrastructureError, VerifierError
from llm_verifier.files import write_json
from llm_verifier.provider import HttpTransport
from llm_verifier.runner import DockerRunner

from export_hosted import ARTIFACT_NAMES, export_hosted, file_hash
from export_metrics import export
from expressions import CONVENTION, VERSION, claims_for, polynomial
from tests.synthetic import synthetic_dossier
from tests.test_expressions import finite
from tools.generate_hosted import VERIFIER_REVISION

ROOT = Path(__file__).resolve().parents[1]
HEAD = "a" * 40


class HostedExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "submission").mkdir()
        for name in ("verifier.json", "rules.md"):
            shutil.copyfile(ROOT / name, self.root / name)
        shutil.copytree(ROOT / "checkers", self.root / "checkers")
        self.intake = self.root / "intake.json"
        self.dossier = self.root / "dossier.json"
        self.wheel = self.root / "llm_verifier-test-py3-none-any.whl"
        self.wheel.write_bytes(b"Synthetic downloaded wheel identity; never installed.\n")
        self.score = self.root / "score.json"
        self.lineage = self.root / "recovery.json"
        self.receipt = self.root / "transfer.json"

    def fixture(self, *, valid=True, simulation=False, original_attempt=1, audit_valid=None, audit_failure=None):
        bound = {"version": VERSION, "convention": CONVENTION,
                 "family": "Synthetic delivery fixture, not a circuit theorem.",
                 "theorem": "SYNTHETIC transport test only.",
                 "guarantee": polynomial(1, 1), "projection": finite()}
        manifest = {"version": 1, "title": "Synthetic delivery fixture",
                    "claims": claims_for(bound), "experiments": [], "checks": []}
        context = {"kind": "github-actions", "candidate_revision": HEAD,
                   "verifier_revision": VERIFIER_REVISION, "run_id": f"123:{original_attempt}"}

        def hosted_intake(*args, **kwargs):
            return intake_snapshot(*args, context=context, **kwargs)

        with patch("tests.synthetic.intake_snapshot", side_effect=hosted_intake):
            record = synthetic_dossier(self.root, bound, manifest, valid=valid, simulation=simulation,
                                       audit_valid=audit_valid, audit_failure=audit_failure)
        write_json(self.dossier, record)
        write_json(self.intake, record["data"]["phases"][0])
        self.descriptors = [{"id": str(100 + i), "name": name, "sha256": str(i + 1) * 64,
                             "run_id": "123", "run_attempt": str(original_attempt)}
                            for i, name in enumerate(ARTIFACT_NAMES)]
        self.env = {
            "GITHUB_RUN_ID": "123", "GITHUB_RUN_ATTEMPT": "2", "GITHUB_SHA": HEAD,
            "LLMV_REVIEW_RUN_ID": "123", "LLMV_REVIEW_RUN_ATTEMPT": str(original_attempt),
            "LLMV_EXPECTED_INTAKE": record["data"]["phases"][0]["sha256"],
            **{key: json.dumps(value) for key, value in zip(
                ("LLMV_PACKAGE_ARTIFACT", "LLMV_INTAKE_ARTIFACT", "LLMV_DOSSIER_ARTIFACT"), self.descriptors)},
        }
        self.transfers = {"version": 1, "run_id": "123", "delivery_attempt": "2", "artifacts": [
            {"artifact": descriptor, "path": destination, "files": {path.name: file_hash(path)}}
            for descriptor, destination, path in zip(self.descriptors, ("wheel", "intake", "dossier"),
                                                    (self.wheel, self.intake, self.dossier))
        ]}
        write_json(self.receipt, self.transfers)
        return record

    def run_export(self):
        with patch.object(HttpTransport, "request", side_effect=AssertionError("No model calls")), \
                patch.object(DockerRunner, "run", side_effect=AssertionError("No Docker calls")):
            export_hosted(self.intake, self.dossier, self.wheel, self.score, self.lineage,
                          self.receipt, env=self.env)

    def test_retry_exports_identical_score_with_original_review_identity(self):
        for valid in (True, False):
            with self.subTest(valid=valid):
                record = self.fixture(valid=valid)
                originals = {p: p.read_bytes() for p in (self.intake, self.dossier, self.wheel, self.receipt)}
                expected = self.root / "normal-score.json"
                export(self.dossier, expected, intake_path=self.intake)
                self.run_export()
                self.assertEqual(self.score.read_bytes(), expected.read_bytes())
                self.assertEqual(originals, {p: p.read_bytes() for p in originals})
                lineage = json.loads(self.lineage.read_text())
                self.assertEqual(set(lineage), {"version", "runId", "headSha", "reviewAttempt",
                    "deliveryAttempt", "artifacts", "inputFiles", "identity", "scoreSha256"})
                self.assertEqual((lineage["runId"], lineage["reviewAttempt"], lineage["deliveryAttempt"]), (123, 1, 2))
                self.assertEqual(lineage["scoreSha256"], file_hash(self.score))
                self.assertEqual(lineage["identity"]["dossierSha256"], record["sha256"])
                self.assertEqual(lineage["identity"]["candidateSha"], HEAD)
                self.assertEqual(lineage["identity"]["verifierRevision"], VERIFIER_REVISION)
                self.assertEqual([x["runAttempt"] for x in lineage["artifacts"]], [1, 1, 1])
                self.assertEqual(lineage["inputFiles"], [
                    {"artifactName": name, "path": path.name, "sha256": file_hash(path)}
                    for name, path in zip(ARTIFACT_NAMES, (self.wheel, self.intake, self.dossier))])

    def test_record_only_disagreement_exports_primary_and_retains_both_decisions(self):
        record = self.fixture(valid=True, audit_valid=False)
        phase = record["data"]["phases"][2]["data"]
        self.assertTrue(phase["selection"]["selected"])
        self.assertEqual(phase["authority"]["judge"], "primary")
        self.assertTrue(phase["primary"]["review"]["decision"]["judgment"]["valid"])
        self.assertFalse(phase["escalation"]["review"]["decision"]["judgment"]["valid"])
        self.assertTrue(phase["comparison"]["available"])
        self.assertTrue(phase["comparison"]["valid_disagreement"])
        original = self.dossier.read_bytes()
        self.run_export()
        score = json.loads(self.score.read_text())
        self.assertTrue(score["valid"])
        self.assertEqual(score["scores"], [1, 1])
        self.assertEqual(self.dossier.read_bytes(), original)

    def test_record_only_failed_audit_exports_primary_without_completed_comparison(self):
        record = self.fixture(valid=True, audit_failure=InfrastructureError("provider_timeout"))
        phase = record["data"]["phases"][2]["data"]
        self.assertEqual(phase["authority"]["judge"], "primary")
        self.assertEqual(phase["escalation"]["status"], "failed")
        failure = phase["escalation"]["failure"]
        self.assertEqual(failure["state"]["error"], "provider_timeout")
        self.assertEqual(len(failure["calls"]), 1)
        self.assertEqual(phase["comparison"], {"available": False,
                         "valid_disagreement": None, "scores_disagreement": None})
        original = self.dossier.read_bytes()
        self.run_export()
        self.assertEqual(json.loads(self.score.read_text())["scores"], [1, 1])
        self.assertEqual(self.dossier.read_bytes(), original)

    def test_full_rerun_preserves_actual_attempt_without_claiming_recovery_eligibility(self):
        self.fixture(original_attempt=2)
        self.run_export()
        lineage = json.loads(self.lineage.read_text())
        self.assertEqual(lineage["reviewAttempt"], 2)
        self.assertEqual([x["runAttempt"] for x in lineage["artifacts"]], [2, 2, 2])

    def test_wrong_source_run_intake_or_attempt_fails_closed(self):
        self.fixture()
        original = dict(self.env)
        for key, value in (("GITHUB_SHA", "b" * 40), ("LLMV_REVIEW_RUN_ID", "124"),
                           ("LLMV_EXPECTED_INTAKE", "0" * 64), ("GITHUB_RUN_ATTEMPT", "0"),
                           ("LLMV_REVIEW_RUN_ATTEMPT", "3")):
            with self.subTest(key=key):
                self.env = {**original, key: value}
                self.score.write_text("stale")
                self.lineage.write_text("stale")
                with self.assertRaises(VerifierError):
                    self.run_export()
                self.assertFalse(self.score.exists())
                self.assertFalse(self.lineage.exists())

    def test_artifact_identity_or_extracted_bytes_cannot_be_substituted(self):
        self.fixture()
        mutations = [lambda r: r["artifacts"][0]["artifact"].update(id="999"),
                     lambda r: r["artifacts"][1]["files"].update({"intake.json": "0" * 64}),
                     lambda r: r["artifacts"][2].update(path="other"),
                     lambda r: r["artifacts"].pop(),
                     lambda r: r.update(delivery_attempt="1")]
        for mutate in mutations:
            changed = copy.deepcopy(self.transfers)
            mutate(changed)
            write_json(self.receipt, changed)
            with self.subTest(mutation=mutate), self.assertRaises(VerifierError):
                self.run_export()
            self.assertFalse(self.score.exists())
            self.assertFalse(self.lineage.exists())

    def test_descriptor_requires_exact_names_run_ids_and_safe_integer_ids(self):
        self.fixture()
        for change in ({"run_id": "124"}, {"name": "candidate-package"}, {"id": "9007199254740992"},
                       {"sha256": "broken"}, {"run_attempt": "3"}):
            self.env["LLMV_PACKAGE_ARTIFACT"] = json.dumps({**self.descriptors[0], **change})
            with self.subTest(change=change), self.assertRaises(VerifierError):
                self.run_export()

    def test_simulation_rejected_and_outputs_cleared(self):
        self.fixture(simulation=True)
        self.score.write_text("stale")
        self.lineage.write_text("stale")
        with self.assertRaisesRegex(VerifierError, "simulation"):
            self.run_export()
        self.assertFalse(self.score.exists())
        self.assertFalse(self.lineage.exists())

    def test_retry_after_partial_delivery_keeps_score_and_replaces_only_delivery_identity(self):
        self.fixture()
        self.run_export()
        # First upload succeeded, then the lineage upload failed. No new review
        # is needed: rerun the same export with only its delivery attempt changed.
        delivered_score = self.score.read_bytes()
        previous_lineage = json.loads(self.lineage.read_text())
        originals = {p: p.read_bytes() for p in (self.intake, self.dossier, self.wheel)}
        self.env["GITHUB_RUN_ATTEMPT"] = "3"
        self.transfers["delivery_attempt"] = "3"
        write_json(self.receipt, self.transfers)
        self.run_export()
        self.assertEqual(self.score.read_bytes(), delivered_score)
        current = json.loads(self.lineage.read_text())
        self.assertEqual(current, {**previous_lineage, "deliveryAttempt": 3})
        self.assertEqual(originals, {p: p.read_bytes() for p in originals})
        # Both derived uploads must replace any partial prior delivery; the
        # immutable original input artifacts are never uploaded by this job.
        overlay = (ROOT / "tools/metrics-job.yml").read_text()
        uploads = overlay.split("      - uses: actions/upload-artifact@")[1:]
        self.assertEqual(len(uploads), 2)
        for block, name in zip(uploads, ("circuitbreaker-score", "circuitbreaker-export-lineage")):
            self.assertIn("name: " + name + "\n", block)
            self.assertIn("overwrite: true", block)

    def test_output_aliases_preserve_trusted_inputs(self):
        self.fixture()
        before = self.intake.read_bytes()
        with self.assertRaisesRegex(VerifierError, "distinct"):
            export_hosted(self.intake, self.dossier, self.wheel, self.intake, self.lineage,
                          self.receipt, env=self.env)
        self.assertEqual(before, self.intake.read_bytes())
        with self.assertRaisesRegex(VerifierError, "distinct"):
            export_hosted(self.intake, self.dossier, self.wheel, self.score, self.receipt,
                          self.receipt, env=self.env)


if __name__ == "__main__":
    unittest.main()
