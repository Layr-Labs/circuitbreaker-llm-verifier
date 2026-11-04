"""Materialize an explicit scoring source; no dispatch, inference or publication."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
PROFILE = "score-derivation-profile.json"
VERSION = "circuitbreaker-score-derivation-profile-v1"
SCORING = {"version": "coefficient-exponent-expression-v3", "retained_result": "yukon-retained-result-v1"}
VALIDITY = {"version": "fixed-family-full-b2-v1"}


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_profile(root):
    """Small stdlib-only preflight runs before either paid workflow branch."""
    root = Path(root)
    value = json.loads((root / PROFILE).read_text())
    if set(value) != {"version", "verifier_revision", "files", "retained_submission_ids",
                      "ordinary_submission_ids", "excluded_submission_ids"} \
            or value["version"] != VERSION \
            or not isinstance(value["verifier_revision"], str) \
            or not re.fullmatch(r"[0-9a-f]{40}", value["verifier_revision"]):
        raise ValueError("invalid score-derivation source profile")
    required = {"derivation.json", "verifier.json", "reorg/score-derivation-plan.json",
                "reorg/score-derivation-report.json"}
    if set(value["files"]) != required or any(file_hash(root / name) != expected
                                              for name, expected in value["files"].items()):
        raise ValueError("source profile configuration/plan digest mismatch")
    plan = json.loads((root / "reorg/score-derivation-plan.json").read_text())
    if plan.get("config_sha256") != value["files"]["derivation.json"] \
            or not plan.get("entries") \
            or [entry["origin"]["submission_id"] for entry in plan["entries"]] != value["retained_submission_ids"]:
        raise ValueError("source profile must contain the complete prepared retained-result plan")
    config = json.loads((root / "verifier.json").read_text())
    if config["policies"] != {"scoring": SCORING, "validity": VALIDITY}:
        raise ValueError("source profile scoring/validity policy mismatch")
    return value


def render_workflow(root, revision):
    root = Path(root)
    jobs = (root / "tools/derivation-jobs.yml").read_text().replace("@VERIFIER_REVISION@", revision)
    header = """name: Benchmark
on:
  workflow_dispatch:
    inputs:
      yukon_accepted_results:
        description: Yukon automatic accepted-result context
        required: false
        type: string
        default: ''
      yukon_evaluation:
        type: string
        required: false
        default: ''
permissions:
  contents: read
  id-token: write
jobs:
"""
    return header + jobs + (root / "tools/metrics-job.yml").read_text()


def render_checks(root, revision):
    return (Path(root) / "tools/derivation-checks.yml").read_text().replace("@VERIFIER_REVISION@", revision)


def regenerate(root):
    root = Path(root)
    profile = read_profile(root)
    benchmark_path = root / "benchmark.json"
    benchmark = json.loads(benchmark_path.read_text())
    if len(benchmark["tracks"]) != 1:
        raise ValueError("CircuitBreaker derivation profile requires the existing single track")
    benchmark["tracks"][0]["paretoPolicies"] = {"scoring": SCORING, "validity": VALIDITY}
    benchmark_path.write_text(json.dumps(benchmark, indent=2, ensure_ascii=False) + "\n")
    (root / ".github/workflows/benchmark.yml").write_text(render_workflow(root, profile["verifier_revision"]))
    (root / ".github/workflows/checker-tests.yml").write_text(render_checks(root, profile["verifier_revision"]))
    # The library CI materializes a fresh profile and pins historical default
    # hashes. A materialized source instead validates its own configuration
    # and registered checkers, without recursively materializing itself.
    (root / ".github/workflows/score-derivation-tests.yml").unlink(missing_ok=True)


def materialize(root, output, snapshot, artifact_index, archives, verifier_revision):
    from llm_verifier.config import load_config
    from llm_verifier.files import write_json
    from tools.prepare_score_derivation_plan import build

    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists() or not output.parent.is_dir() or output.is_relative_to(root) or root.is_relative_to(output):
        raise ValueError("profile output must be a new, separate directory")
    if not re.fullmatch(r"[0-9a-f]{40}", verifier_revision):
        raise ValueError("an immutable reviewed verifier commit is required")
    config_raw = (root / "derivation.json").read_bytes()
    # Completeness and every available original artifact are checked before an
    # output source exists. Missing evidence cannot start an ordinary baseline.
    plan, report = build(snapshot, artifact_index, config_raw, archives)
    with tempfile.TemporaryDirectory(dir=output.parent) as directory:
        staging = Path(directory) / "source"
        shutil.copytree(root, staging, ignore=shutil.ignore_patterns(
            ".git", "__pycache__", ".venv", ".yukon", ".work", "*.pyc"))
        config_path = staging / "verifier.json"
        config = json.loads(config_path.read_text())
        if config["policies"]["validity"] != VALIDITY:
            raise ValueError("unexpected validity contract")
        prior_selector = deepcopy(config["policies"]["scoring"])
        policies = config["supported_policies"]["scoring"]
        matched = [item for item in policies if item["selector"] == prior_selector]
        if len(matched) != 1:
            raise ValueError("ordinary scoring instructions are ambiguous")
        config["policies"]["scoring"] = deepcopy(SCORING)
        if any(item["selector"] == SCORING for item in policies):
            raise ValueError("input source already contains a derivation profile")
        policies.append({**deepcopy(matched[0]), "selector": deepcopy(SCORING)})
        write_json(config_path, config)
        # Confirm the ordinary branch supports the selected scoring policy.
        load_config(staging, "verifier.json", "submission")
        write_json(staging / "reorg/score-derivation-plan.json", plan)
        write_json(staging / "reorg/score-derivation-report.json", report)
        names = ("derivation.json", "verifier.json", "reorg/score-derivation-plan.json",
                 "reorg/score-derivation-report.json")
        profile = {"version": VERSION, "verifier_revision": verifier_revision,
                   "files": {name: file_hash(staging / name) for name in names},
                   "retained_submission_ids": [entry["origin"]["submission_id"] for entry in plan["entries"]],
                   "ordinary_submission_ids": report["ordinary_submission_ids"],
                   "excluded_submission_ids": report["excluded_submission_ids"]}
        write_json(staging / PROFILE, profile)
        regenerate(staging)
        staging.rename(output)
    return profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "artifact-index", "archive-root", "verifier-revision", "out"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    materialize(ROOT, args.out, json.loads(Path(args.snapshot).read_text()),
                json.loads(Path(args.artifact_index).read_text()), args.archive_root, args.verifier_revision)
    print("Prepared isolated source profile; no reorg, inference or publication executed.")


if __name__ == "__main__":
    main()
