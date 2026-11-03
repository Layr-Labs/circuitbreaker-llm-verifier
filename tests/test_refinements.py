"""Synthetic historical transport only: no real proofs, providers or approval.

The frozen old-source/actual-wheel end-to-end replay is separately checked on
organizer artifacts. Here sealed synthetic chains exercise replay rejection and
output assembly without depending on Git, authentic submissions or network.
"""

import copy
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from llm_verifier import pipeline
from llm_verifier.artifacts import intake_snapshot
from llm_verifier.files import write_json

import export_metrics
import export_refinements as exporter
import expressions
from expressions import CONVENTION, VERSION, VERSION_V2, claims_for, polynomial, ratio
from refinement import replay_original as replay
from refinement.common import RefinementError, canonical, digest, read_json, sha256
from tests.synthetic import synthetic_dossier
from tests.test_reuse import CONTEXT

ROOT = Path(__file__).resolve().parents[1]


def synthetic_origin(root):
    source, bundle = root / "source", root / "bundle"
    (source / "submission").mkdir(parents=True)
    for name in ("verifier.json", "rules.md"):
        shutil.copyfile(ROOT / "tests/fixtures/legacy-v1" / name, source / name)
    shutil.copytree(ROOT / "checkers", source / "checkers")
    bound = {"version": VERSION, "convention": CONVENTION,
             "family": "SYNTHETIC language — transport only 😀.",
             "theorem": "SYNTHETIC: eventual 4n - 7√n. DEL: \u007f. No theorem or approval asserted.",
             "guarantee": polynomial(4, 1),
             "projection": {"kind": "finite", "coefficient": ratio(4), "exponent": ratio(1)}}
    manifest = {"version": 1, "title": "SYNTHETIC offline refinement", "claims": claims_for(bound),
                "experiments": [], "checks": []}
    context = {"kind": "github-actions", "candidate_revision": "a" * 40,
               "verifier_revision": replay.VERIFIER_REF, "run_id": "123:1"}
    snapshot = {"version": 1, "snapshot_id": "synthetic-no-network",
                "consumer": {k: CONTEXT[k] for k in
                             ("submission_id", "benchmark_id", "challenge_id", "evaluation_id", "job_id")},
                "references": [], "sources": [], "ancestors": []}

    def intake(*args, **kwargs):
        return intake_snapshot(*args, context=context, accepted_results_context=CONTEXT, **kwargs)

    with patch("tests.synthetic.intake_snapshot", side_effect=intake), \
         patch("llm_verifier.accepted.resolve", return_value=snapshot):
        dossier = synthetic_dossier(source, bound, manifest, valid=False, audit_valid=True)
    intake, initial, review = dossier["data"]["phases"]
    outputs = {"llmv-intake/intake.json": intake, "llmv-initial-experiments/initial.json": initial,
               "llmv-initial-review/review.json": review, "verifier-dossier/dossier.json": dossier,
               "benchmark-score/.yukon/score.json": dossier["data"]["yukon_projection"],
               "verifier-accepted-results/accepted-results.json": pipeline.export_results(dossier)}
    for name, value in outputs.items():
        write_json(bundle / name, value)
    export_metrics.export(bundle / "verifier-dossier/dossier.json", bundle / "circuitbreaker-score/.yukon/score.json")
    score = read_json(bundle / "circuitbreaker-score/.yukon/score.json", numbers=True)
    evidence = {"resultId": "synthetic-score-result", "jobId": CONTEXT["job_id"],
                "evaluatorId": CONTEXT["evaluation_id"], "policyVersionId": "synthetic-policy",
                "commitSha": context["candidate_revision"], "sourceRef": replay.ORIGINAL_REF,
                "evaluatedAt": "2026-01-01T00:00:00Z"}
    public = {"id": CONTEXT["submission_id"], "benchmarkId": CONTEXT["benchmark_id"], "status": "review",
              "officialMetrics": score["metrics"], "evaluation": {
                  "id": "synthetic-packet", "candidateKey": "digest:" + "b" * 64,
                  "evaluationId": CONTEXT["evaluation_id"], "validityVersionId": "synthetic-validity-policy",
                  "human": None, "scores": {"values": score["scores"], "evidence": evidence},
                  "validity": {"valid": True, "reason": dossier["data"]["judgment"]["reason"],
                               "evidence": {**evidence, "resultId": "synthetic-validity-result"}}}}
    write_json(root / "platform.json", {"submission": public})
    lineage = {"scoreSha256": sha256((bundle / "circuitbreaker-score/.yukon/score.json").read_bytes()),
               "identity": {"intakeSha256": intake["sha256"], "dossierSha256": dossier["sha256"],
                            "verifierSha256": intake["bindings"]["verifier"],
                            "configSha256": intake["bindings"]["config"],
                            "submissionSha256": intake["bindings"]["submission"],
                            "candidateSha": context["candidate_revision"], "verifierRevision": replay.VERIFIER_REF}}
    write_json(bundle / "circuitbreaker-export-lineage/recovery.json", lineage)
    target = copy.deepcopy(bound)
    target.update(version=VERSION_V2, guarantee={"kind": "eventual_lower_bound", "leading": bound["guarantee"],
                  "remainder": {"kind": "explicit", "terms": [{"coefficient": ratio(7), "exponent": ratio(1, 2),
                                 "log_exponent": ratio(0), "loglog_exponent": ratio(0)}]}})
    return source, bundle, public, target


def replay_synthetic(source, bundle, public, output):
    output.mkdir()
    # Current exporter's explicit legacy branch has the same v1 semantics.
    # Real frozen source loading/wheel identity is tested separately, not faked
    # as an authenticated historical result by these synthetic fixtures.
    with patch.object(replay, "load_original", return_value=(export_metrics, expressions)):
        return replay.replay(source, bundle, public, output)


class CanonicalTests(unittest.TestCase):
    def test_cross_language_goldens(self):
        for vector in read_json(ROOT / "tests/fixtures/owner-refinement-canonical.json"):
            with self.subTest(vector=vector["label"]):
                self.assertEqual(canonical(vector["input"]).decode(), vector["canonical"])
                self.assertEqual(digest(vector["input"]), vector["sha256"])

    def test_synthetic_golden_namespace_target_and_record_hashes(self):
        fixture = read_json(ROOT / "tests/fixtures/owner-refinement-synthetic.json", numbers=True)
        self.assertTrue(fixture["synthetic"])
        record = fixture["representationRefinements"]["records"][0]
        metrics = fixture["submission"]["officialMetrics"]
        self.assertEqual(record["origin"]["originalEnvelopeSha256"], digest(metrics["circuitbreaker_frontier"]))
        self.assertEqual(record["origin"]["originalRichSha256"], digest(metrics["circuitbreaker"]))
        self.assertEqual(record["recordSha256"], digest({k: v for k, v in record.items() if k != "recordSha256"}))
        exporter.validate_target(metrics["circuitbreaker"]["bound"], record["targetBound"])

    def test_reject_numbers_non_json_and_lone_surrogates(self):
        for value in (1, 1.5, float("nan"), {"x": -1}, "\ud800", {"\udfff": "x"}, ("x",), {1: "x"}):
            with self.subTest(value=repr(value)), self.assertRaises(RefinementError):
                canonical(value)
        self.assertEqual(canonical([True, False, None]), b"[true,false,null]")

    def test_reject_duplicate_json_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.json"
            path.write_text('{"a":"first","a":"last"}')
            with self.assertRaisesRegex(RefinementError, "duplicate"):
                read_json(path)


class RefinementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.fixture = Path(cls.temp.name) / "fixture"
        cls.source, cls.bundle, cls.public, cls.target = synthetic_origin(cls.fixture)
        cls.report = replay_synthetic(cls.source, cls.bundle, cls.public, cls.fixture / "replayed")

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.item = {"id": "synthetic-refinement", "bundle": str(self.bundle),
                     "downloadReceiptSha256": "d" * 64, "platform": str(self.fixture / "platform.json"),
                     "platformSha256": sha256((self.fixture / "platform.json").read_bytes()),
                     "originalSource": str(self.source), "originalPython": "synthetic-venv/bin/python",
                     "targetBound": copy.deepcopy(self.target)}

    def write_plan(self, records=None):
        path = self.root / "plan.json"
        path.write_bytes(canonical({"version": exporter.PLAN_VERSION, "records": records or [self.item]}) + b"\n")
        return path

    def test_original_full_claims_selected_override_and_bytes_retained(self):
        self.assertEqual(self.report["originalStatus"], {"submission": "review", "human": None,
                                                       "metric": "accepted_provisional"})
        self.assertEqual(len(self.report["byteIdenticalOutputs"]), 4)
        self.assertEqual(self.report["origin"]["originalRichSha256"],
                         digest(self.public["officialMetrics"]["circuitbreaker"]))
        self.assertNotEqual(self.public["evaluation"]["candidateKey"],
                            "digest:" + self.report["lineage"]["submissionSha256"])
        exporter.validate_target(self.report["originalBound"], self.target)

    def test_changed_target_scope_leading_projection_or_version_rejected(self):
        changes = (lambda t: t.update(family="another family"), lambda t: t.update(theorem="new claim"),
                   lambda t: t.update(convention="different-model"), lambda t: t.update(version=VERSION),
                   lambda t: t["projection"].update(coefficient=ratio(3)),
                   lambda t: t["guarantee"]["leading"].update(coefficient=ratio(5)),
                   lambda t: t["guarantee"]["remainder"]["terms"][0].update(exponent=ratio(2)))
        for change in changes:
            target = copy.deepcopy(self.target)
            change(target)
            with self.assertRaises(ValueError):
                exporter.validate_target(self.report["originalBound"], target)

    def test_assembly_binds_raw_manifest_prior_receipt_and_immutable_origin(self):
        path = self.write_plan()
        before = copy.deepcopy(self.report)
        with patch.object(exporter, "run_replay", return_value=self.report):
            overlay = exporter.export(path, self.root / "out")
        record = overlay["records"][0]
        self.assertEqual(self.report, before)
        self.assertEqual(record["recordSha256"], digest({k: v for k, v in record.items() if k != "recordSha256"}))
        self.assertEqual(record["audit"]["manifestSha256"], sha256(path.read_bytes()))
        receipt_raw = (self.root / "out/receipt.json").read_bytes()
        self.assertEqual(record["audit"]["receiptSha256"], sha256(receipt_raw))
        self.assertNotIn(record["recordSha256"].encode(), receipt_raw)
        self.assertNotIn("approved", record)
        self.assertEqual(set(record), {"id", "origin", "targetBound", "audit", "recordSha256"})
        with self.assertRaisesRegex(RefinementError, "must not exist"):
            exporter.export(path, self.root / "out")

    def test_alias_duplicate_and_failed_replay_create_no_output(self):
        path = self.write_plan()
        with self.assertRaisesRegex(RefinementError, "overlaps"):
            exporter.export(path, self.source / "mapping-output")
        for records, expected in (([self.item, self.item], "duplicate refinement id"),
                                  ([self.item, {**self.item, "id": "another"}], "duplicate origin")):
            path = self.write_plan(records)
            with patch.object(exporter, "run_replay", return_value=self.report), self.assertRaisesRegex(RefinementError, expected):
                exporter.export(path, self.root / "out")
            self.assertFalse((self.root / "out").exists())
        path = self.write_plan()
        with patch.object(exporter, "run_replay", side_effect=RefinementError("unmatched original")), self.assertRaisesRegex(RefinementError, "unmatched"):
            exporter.export(path, self.root / "out")
        self.assertFalse((self.root / "out").exists())

    def test_selected_packet_tampering_and_human_promotion_rejected(self):
        for field in ("sourceRef", "commitSha", "jobId", "evaluatorId"):
            public = copy.deepcopy(self.public)
            public["evaluation"]["validity"]["evidence"][field] = "changed"
            with self.subTest(field=field), self.assertRaisesRegex(RefinementError, "evidence identity"):
                replay_synthetic(self.source, self.bundle, public, self.root / field)
        for name, update in (("invalid", lambda p: p["evaluation"]["validity"].update(valid=False)),
                             ("human", lambda p: p["evaluation"].update(human={"valid": True})),
                             ("status", lambda p: p.update(status="accepted")),
                             ("metric", lambda p: p["officialMetrics"]["circuitbreaker"].update(status="accepted"))):
            public = copy.deepcopy(self.public); update(public)
            with self.subTest(name=name), self.assertRaises(RefinementError):
                replay_synthetic(self.source, self.bundle, public, self.root / name)

    def test_altered_score_review_and_submitted_file_rejected(self):
        for name in ("circuitbreaker-score/.yukon/score.json", "llmv-initial-review/review.json"):
            bundle = self.root / name.split("/")[0]
            shutil.copytree(self.bundle, bundle)
            path = bundle / name
            path.write_bytes(path.read_bytes() + b" ")
            if "review" in name:
                value = read_json(path, numbers=True); value["sha256"] = "0" * 64; write_json(path, value)
            with self.subTest(name=name), self.assertRaises(RefinementError):
                replay_synthetic(self.source, bundle, self.public, self.root / (bundle.name + "-out"))
        source = self.root / "altered-source"; shutil.copytree(self.source, source)
        (source / "submission/proof.md").write_text("Changed proof")
        with self.assertRaisesRegex(RefinementError, "submitted source"):
            replay_synthetic(source, self.bundle, self.public, self.root / "altered-out")

    def test_venv_symlink_and_cache_isolation(self):
        interpreter = self.root / "venv/bin/python"
        interpreter.parent.mkdir(parents=True)
        interpreter.symlink_to("/usr/bin/python3")
        self.item["originalPython"] = str(interpreter)
        path = self.write_plan()
        with patch.object(exporter, "run_replay", return_value=self.report) as runner:
            exporter.export(path, self.root / "out")
        self.assertEqual(runner.call_args.args[0]["originalPython"], str(interpreter))
        scratch = self.root / "scratch"; scratch.mkdir()
        output = scratch / "replay"; output.mkdir(); (output / "replay.json").write_bytes(canonical(self.report))
        with patch.object(exporter.subprocess, "run", return_value=SimpleNamespace(returncode=0, stderr="")) as process:
            exporter.run_replay(self.item, scratch)
        args, kwargs = process.call_args
        self.assertIn("pycache_prefix=" + str(scratch / "empty-bytecode-cache"), args[0])
        self.assertFalse((scratch / "empty-bytecode-cache").exists())
        self.assertIn("-B", args[0]); self.assertIn("-I", args[0])
        self.assertEqual(set(kwargs["env"]), {"PATH"})


class SourceAndArchiveTests(unittest.TestCase):
    def test_source_manifest_checks_exact_bytes_before_import(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "x.py").write_text("original source\n")
            manifest = root / "manifest.json"
            manifest.write_bytes(canonical({"sourceRef": replay.ORIGINAL_REF,
                                           "files": {"x.py": sha256((root / "x.py").read_bytes())}}))
            with patch.object(replay, "SOURCE_MANIFEST", manifest):
                replay.verify_source(root)
                (root / "x.py").write_text("changed source\n")
                with self.assertRaisesRegex(RefinementError, "source mismatch"):
                    replay.verify_source(root)

    def test_archive_pin_and_extracted_bytes_are_independently_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifacts, listed = [], []
            for index, name in enumerate(sorted(replay.ARTIFACT_NAMES)):
                (root / name).mkdir(); (root / name / "file.txt").write_bytes(b"synthetic")
                path = root / (name + ".zip")
                with zipfile.ZipFile(path, "w") as archive:
                    archive.writestr("file.txt", b"synthetic")
                raw = path.read_bytes()
                artifacts.append({"id": index, "name": name, "archive": path.name,
                                  "archiveSha256": sha256(raw), "archiveBytes": len(raw),
                                  "files": [{"path": "file.txt", "sha256": sha256(b"synthetic"), "bytes": 9}]})
                listed.append({"id": index, "name": name, "digest": "sha256:" + sha256(raw),
                               "size_in_bytes": len(raw), "workflow_run": {"id": 123, "head_sha": "a" * 40}})
            write_json(root / "run.json", {"id": 123, "head_sha": "a" * 40, "run_attempt": 1,
                       "repository": {"full_name": replay.REPOSITORY}, "status": "completed",
                       "conclusion": "success", "event": "workflow_dispatch"})
            write_json(root / "artifacts.json", {"artifacts": listed})
            receipt = {"runId": 123, "headSha": "a" * 40, "runAttempt": 1, "repository": replay.REPOSITORY,
                       "runMetadataSha256": sha256((root / "run.json").read_bytes()),
                       "artifactListingSha256": sha256((root / "artifacts.json").read_bytes()), "artifacts": artifacts}
            write_json(root / "download-receipt.json", receipt)
            pin = sha256((root / "download-receipt.json").read_bytes())
            replay.verify_download(root, pin)
            with self.assertRaisesRegex(RefinementError, "receipt hash"):
                replay.verify_download(root, "f" * 64)
            (root / artifacts[0]["name"] / "file.txt").write_bytes(b"changed")
            with self.assertRaisesRegex(RefinementError, "extracted artifact"):
                replay.verify_download(root, pin)
            with self.assertRaisesRegex(RefinementError, "unsafe"):
                replay.inside(root, "../escape")


if __name__ == "__main__":
    unittest.main()
