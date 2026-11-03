"""Prevent chart grouping and table configuration from drifting from the producer."""

import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from tools.generate_hosted import validate_presentation

ROOT = Path(__file__).resolve().parents[1]


class PresentationTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        shutil.copytree(ROOT / "examples", self.root / "examples")
        self.track = json.loads((ROOT / "benchmark.json").read_text())["tracks"][0]
        self.path = self.root / "examples/circuitbreaker-presentation.json"
        self.page = json.loads(self.path.read_text())

    def validate(self):
        self.path.write_text(json.dumps(self.page))
        validate_presentation(self.root, self.track)

    def test_combined_template_matches_current_producer_without_live_binding(self):
        self.validate()
        self.assertEqual(self.page["participation"]["benchmarkRef"], "unbound/circuitbreaker/main")
        self.assertNotIn("benchmarkId", self.page)
        self.assertNotIn("published", self.page)
        self.assertNotIn("presentation", self.track)

    def test_rejects_wrong_objective_order_direction_and_mode(self):
        original = copy.deepcopy(self.track)
        for change in ("order", "direction", "mode"):
            with self.subTest(change=change):
                self.track = copy.deepcopy(original)
                if change == "order":
                    self.track["objectives"].reverse()
                elif change == "direction":
                    self.track["objectives"][1]["direction"] = "-"
                else:
                    self.track["mode"] = "scalar"
                with self.assertRaisesRegex(ValueError, "objective/display"):
                    self.validate()

    def test_rejects_wrong_grouping_or_charted_coordinate(self):
        original = copy.deepcopy(self.page)
        for field, value in (("scoreTracks", {"objective": "gate_coefficient", "sort": "desc", "default": "highest"}),
                             ("defaultView", {"mode": "record", "objective": "gate_exponent"})):
            with self.subTest(field=field):
                self.page = copy.deepcopy(original)
                self.page["results"][field] = value
                with self.assertRaisesRegex(ValueError, "objective/display"):
                    self.validate()

    def test_rejects_incompatible_metric_path_version_and_comparator(self):
        original = copy.deepcopy(self.page)
        for field, value in (("metricsPath", "circuitbreaker"), ("metricsVersion", "wrong-v1"),
                             ("comparator", "asymptotic-v1")):
            with self.subTest(field=field):
                self.page = copy.deepcopy(original)
                self.page["resultTabs"][0][field] = value
                with self.assertRaisesRegex(ValueError, "Complexity tab"):
                    self.validate()

    def test_rejects_shared_wrong_contract_in_both_presentation_files(self):
        original = copy.deepcopy(self.page)
        fragment_path = self.root / "examples/complexity-result-tab.json"
        for field, value in (("metricsVersion", "circuitbreaker-frontier-v1"),
                             ("comparator", "circuitbreaker-v1"),
                             ("scope", "different-scope"), ("model", "monotone"),
                             ("notation", "upper-bound"), ("columns", [])):
            with self.subTest(field=field):
                self.page = copy.deepcopy(original)
                self.page["resultTabs"][0][field] = value
                fragment_path.write_text(json.dumps({"resultTabs": self.page["resultTabs"]}))
                with self.assertRaisesRegex(ValueError, "exact v2 producer contract"):
                    self.validate()
