"""Synthetic prompt-size regressions; no model calls or scientific judgments."""

import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from llm_verifier import pipeline
from llm_verifier.artifacts import intake_snapshot
from llm_verifier.config import load_config
from llm_verifier.errors import InfrastructureError, VerifierError
from llm_verifier.jsonutil import canonical
from llm_verifier.provider import HttpTransport
from llm_verifier.runner import DockerRunner, output_record

from export_metrics import projection
from tests.synthetic import scripted_http, source_plan, synthetic_review

ROOT = Path(__file__).resolve().parents[1]


class PromptBudgetTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        for name in ("verifier.json", "rules.md"):
            shutil.copyfile(ROOT / name, self.root / name)
        for name in ("checkers", "submission"):
            shutil.copytree(ROOT / name, self.root / name)
        manifest_path = self.root / "submission/submission.json"
        manifest = json.loads(manifest_path.read_text())
        manifest.update(experiments=[], checks=[])
        manifest_path.write_text(json.dumps(manifest))
        # Seven files with enough source to exercise the duplicated numbered
        # planner context. This is deliberately not a claimed circuit proof.
        for i in range(4):
            (self.root / "submission" / f"support-{i}.md").write_text(
                ("Synthetic source line for prompt-budget regression. " + "x" * 70 + "\n") * 100
            )

    def intake(self):
        return intake_snapshot(load_config(self.root, "verifier.json", "submission"),
                               self.root / "submission", simulation=True)

    def test_seven_file_planner_exceeds_old_cap_without_losing_source(self):
        intake = self.intake()
        requests = []
        initial = pipeline.execute(intake)
        review = synthetic_review(intake, initial, valid=False, requests=requests)
        self.assertEqual(len(intake["data"]["files"]), 7)
        self.assertEqual(set(review["data"]["region_plan"]["sources"]),
                         set(intake["data"]["files"]))
        size = len(canonical(requests[0]["input"]))
        self.assertGreater(size, 131072)
        self.assertLessEqual(size, intake["data"]["effective"]["config"]["limits"]["max_prompt_bytes"])

        config_path = self.root / "verifier.json"
        config = json.loads(config_path.read_text())
        config["limits"]["max_prompt_bytes"] = 131072
        config_path.write_text(json.dumps(config))
        old_intake = self.intake()
        old_initial = pipeline.execute(old_intake)
        with patch.object(HttpTransport, "request") as http:
            with self.assertRaisesRegex(InfrastructureError, "full_evidence_context_budget"):
                pipeline.review(old_intake, old_initial)
            http.assert_not_called()

    def test_eight_region_aggregated_reports_fit_without_expanding_other_budgets(self):
        (self.root / "submission/support-4.md").write_text("Synthetic eighth source.\n")
        intake = self.intake()
        initial = pipeline.execute(intake)
        requests = []
        first = synthetic_review(intake, initial, request_followup=True, requests=requests)
        fixture = json.loads((ROOT / "checker-tests.json").read_text())[0]
        output = output_record("completed", 0, canonical(
            {**fixture["expected"], "message": "Synthetic fixture"}), b"")
        with patch.object(DockerRunner, "run", return_value=output):
            followup = pipeline.execute(intake, first, initial)

        plan = source_plan(intake)
        findings = [{"id": f"f{i}", "kind": "blocking", "reason": "R" * 4000,
                     "evidence_refs": ["source:proof.md"]} for i in range(16)]
        reports = [{"assessment": "A" * 8000, "findings": findings, "requests": []}]
        names = [f"proof/f{i}" for i in range(16)]
        for region in plan["regions"]:
            reports.append({"assessment": "A" * 8000, "findings": findings, "requests": [],
                            "coverage": {"region_id": region["id"], "start_line": region["start_line"],
                                         "end_line": region["end_line"], "complete": True}})
            names.extend(f"region/{region['id']}/f{i}" for i in range(16))
        responses = [{"finding_id": name, "reason": "D" * 1000,
                      "evidence_refs": ["source:proof.md"]} for name in names]
        reports.append({"responses": responses})
        reports.append({"decision": {"kind": "judgment", "requests": [], "judgment": {
            "valid": False, "scores": None, "reason": "Synthetic rejection; no scientific judgment.",
            "evidence_refs": ["source:proof.md"],
            "review_standard": intake["data"]["effective"]["config"]["review_standard"]}},
            "resolutions": [{**response, "outcome": "confirmed"} for response in responses]})
        with scripted_http(reports, requests):
            final = pipeline.review(intake, initial, first, followup)
        record = pipeline.dossier(intake, initial, first, followup, final)
        self.assertEqual(first["data"]["region_plan"], final["data"]["region_plan"])
        self.assertEqual(len(requests), 21)
        self.assertEqual(sum(phase["data"]["provenance"]["allotted_completion_tokens"]
                             for phase in (first, final)), 688128)
        largest = max(len(canonical(wire["input"])) for wire in requests)
        self.assertGreater(largest, 1048576)
        self.assertLessEqual(largest, 2097152)
        self.assertTrue(all(phase["simulation"] for phase in record["data"]["phases"]))
        with self.assertRaisesRegex(VerifierError, "simulation"):
            projection(record)
