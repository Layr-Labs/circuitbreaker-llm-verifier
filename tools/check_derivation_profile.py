"""Validate a materialized source profile offline, without provider or dispatch."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from llm_verifier import derivation
from llm_verifier.config import load_config
from llm_verifier.jsonutil import digest, loads

from check_submission import check_submission
from export_metrics import validate_contract
from tools.derivation_profile import read_profile, render_checks, render_workflow
from tools.generate_hosted import BASELINE_RECOVERY, validate_presentation


def validate(root):
    root = Path(root)
    profile = read_profile(root)
    config = derivation.normalize_config(loads((root / "derivation.json").read_bytes()))
    if config["plan_path"] != "reorg/score-derivation-plan.json":
        raise ValueError("profile plan path differs from prepared plan")
    plan = loads((root / config["plan_path"]).read_bytes())
    if plan["version"] != derivation.PLAN_VERSION:
        raise ValueError("unsupported profile plan version")
    origins = set()
    for entry in plan["entries"]:
        derivation.validate_entry(entry)
        key = (entry["origin"]["benchmark_id"], entry["origin"]["submission_id"])
        if key in origins:
            raise ValueError("duplicate retained origin")
        origins.add(key)
    report = loads((root / "reorg/score-derivation-report.json").read_bytes())
    if report["plan_sha256"] != digest(plan) or report["score_only"] != len(origins) \
            or report["derive"] != sum(entry["action"] == "derive" for entry in plan["entries"]) \
            or report["not_applicable"] != sum(entry["action"] == "not_applicable" for entry in plan["entries"]) \
            or report["ordinary_submission_ids"] != profile["ordinary_submission_ids"] \
            or report["excluded_submission_ids"] != profile["excluded_submission_ids"]:
        raise ValueError("profile report/plan binding mismatch")
    check_submission(root)
    effective = load_config(root, "verifier.json", "submission")
    validate_contract(effective["config"])
    benchmark = json.loads((root / "benchmark.json").read_text())
    if len(benchmark["tracks"]) != 1:
        raise ValueError("unexpected profile track count")
    track = benchmark["tracks"][0]
    validate_presentation(root, track)
    if track["paretoPolicies"] != effective["config"]["policies"] \
            or track["runner"]["scoreArtifactName"] != "circuitbreaker-score" \
            or track["runner"]["baselineRecovery"] != BASELINE_RECOVERY \
            or track["benchmarkCommand"][:2] != ["python3", "finalize_score.py"]:
        raise ValueError("profile benchmark/execution contract mismatch")
    revision = profile["verifier_revision"]
    if (root / ".github/workflows/benchmark.yml").read_text() != render_workflow(root, revision) \
            or (root / ".github/workflows/checker-tests.yml").read_text() != render_checks(root, revision) \
            or (root / ".github/workflows/score-derivation-tests.yml").exists():
        raise ValueError("profile workflow generation differs from prepared source")
    return {"valid": True, "retained_records": len(origins), "verifier_revision": revision,
            "provider_calls": 0, "dispatched": False}


if __name__ == "__main__":
    print(json.dumps(validate(ROOT), sort_keys=True))
