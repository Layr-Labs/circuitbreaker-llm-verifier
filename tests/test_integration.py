"""Challenge-specific integration guards; shared runtime tests stay upstream."""

import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from check_submission import check_submission
from tools.generate_hosted import BASELINE_RECOVERY, VERIFIER_REVISION, check_generated

ROOT = Path(__file__).resolve().parents[1]


class IntegrationTests(unittest.TestCase):
    def test_baseline_and_scientific_model_preserved(self):
        # Frozen original proof: this restructure claims no scientific repair.
        self.assertEqual(hashlib.sha256((ROOT / "submission/proof.md").read_bytes()).hexdigest(),
                         "bcc793e84a7a44f5bbe1ed7c0d9057ea909346282793978d58743b63708e2f7e")
        rules = (ROOT / "rules.md").read_text()
        family = rules.split("## Family, witnesses and circuit model\n\n")[1].split("## Finite score")[0]
        self.assertEqual(hashlib.sha256(family.encode()).hexdigest(), "43b1042cc167db496a8438b6cd6e6ed5f811e54680e6d94781691ec5829a9d44")
        self.assertIn("every fixed real p>=1", rules)
        self.assertIn("every sufficiently large integer length", rules)
        self.assertIn("base 2", rules)
        self.assertIn("No numeric Infinity", rules)
        self.assertTrue(check_submission(ROOT)["valid"])

    def test_hosted_artifact_is_metrics_output_not_shared_raw_score(self):
        track = json.loads((ROOT / "benchmark.json").read_text())["tracks"][0]
        self.assertEqual(track["objectives"], [{"name": "gate_coefficient", "direction": "+"}, {"name": "gate_exponent", "direction": "+"}])
        self.assertEqual(track["runner"]["scoreArtifactName"], "circuitbreaker-score")
        self.assertEqual(track["benchmarkCommand"][:2], ["python3", "finalize_score.py"])
        workflow = (ROOT / ".github/workflows/benchmark.yml").read_text()
        self.assertIn("needs: [verify]", workflow)
        self.assertIn("name: circuitbreaker-score", workflow)
        self.assertEqual(track["runner"]["baselineRecovery"], BASELINE_RECOVERY)
        for output in ("intake_artifact", "package_artifact", "dossier_artifact",
                       "helper_base64", "helper_sha256", "review_run_id", "review_run_attempt"):
            self.assertIn("needs.verify.outputs." + output, workflow)
        self.assertNotIn("actions/download-artifact@", workflow)
        self.assertEqual(workflow.count("overwrite: true"), 2)
        self.assertIn("name: circuitbreaker-export-lineage", workflow)
        self.assertIn("--transfer-receipt", workflow)
        self.assertIn("--intake", workflow)
        self.assertIn("export_hosted.py", workflow)
        self.assertIn("include-hidden-files: true", workflow)
        self.assertIn("yukon_accepted_results", workflow)
        self.assertIn("@" + VERIFIER_REVISION, workflow)
        self.assertNotIn("pull_request", workflow)
        checker = (ROOT / ".github/workflows/checker-tests.yml").read_text()
        self.assertIn("pull_request", checker)
        self.assertIn("llmv-checker-package", checker)
        self.assertIn("python\" -m unittest -v", checker)
        self.assertIn("tools/generate_hosted.py --check", checker)
        self.assertNotIn("OPENAI_API_KEY", checker)

    def test_regeneration_has_no_drift(self):
        check_generated(ROOT)

    def test_drift_check_detects_missing_metrics_overlay(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "source"
            shutil.copytree(ROOT, root, ignore=shutil.ignore_patterns(".git", "__pycache__"))
            p = root / "benchmark.json"
            m = json.loads(p.read_text()); m["tracks"][0]["runner"]["scoreArtifactName"] = "benchmark-score"
            p.write_text(json.dumps(m))
            with self.assertRaisesRegex(SystemExit, "drift"):
                check_generated(root)
