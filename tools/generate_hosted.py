"""Regenerate upstream callers, then apply the reviewed CircuitBreaker overlay."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path
import subprocess
import sys

VERIFIER_REVISION = "8a7b8efc6fd7fe43bd9aae6461dc8a5a11d863ce"
ROOT = Path(__file__).resolve().parents[1]
BASELINE_RECOVERY = {
    "version": 1,
    "deliveryJob": "Export accepted CircuitBreaker expressions",
    "lineageArtifactName": "circuitbreaker-export-lineage",
    "lineagePath": "recovery.json",
    "artifactNames": ["llmv-package", "llmv-intake", "verifier-dossier"],
}


def validate_presentation(root, track):
    """Keep the separate page template bound to the generated scoring contract."""
    root = Path(root)
    page = json.loads((root / "examples/circuitbreaker-presentation.json").read_text())
    results = page["results"]
    objectives = [{key: objective[key] for key in ("name", "direction")}
                  for objective in results["objectives"]]
    if (results["kind"] != track["mode"] or results["direction"] != track["direction"]
            or objectives != track["objectives"]
            or [objective["name"] for objective in objectives] != ["gate_coefficient", "gate_exponent"]
            or results["scoreTracks"] != {"objective": "gate_exponent", "sort": "desc", "default": "highest"}
            or results["defaultView"] != {"mode": "record", "objective": "gate_coefficient"}):
        raise ValueError("Presentation objective/display contract differs from generated benchmark")
    tabs = page["resultTabs"]
    expected = json.loads((root / "examples/complexity-result-tab.json").read_text())["resultTabs"]
    # The full page and table-only fragment intentionally use different explanatory prose.
    def contract(items):
        return [{key: value for key, value in tab.items() if key != "description"}
                for tab in items]
    if contract(tabs) != contract(expected):
        raise ValueError("Presentation Complexity tab differs from the producer metric contract")


def generate(root=ROOT):
    root = Path(root)
    subprocess.run([
        sys.executable, "-m", "llm_verifier", "github", "--root", str(root),
        "--verifier-revision", VERIFIER_REVISION,
    ], check=True)
    path = root / "benchmark.json"
    benchmark = json.loads(path.read_text())
    track = benchmark["tracks"][0]
    validate_presentation(root, track)
    track["runner"]["scoreArtifactName"] = "circuitbreaker-score"
    track["runner"]["baselineRecovery"] = BASELINE_RECOVERY
    track["benchmarkCommand"][:2] = ["python3", "finalize_score.py"]
    path.write_text(json.dumps(benchmark, indent=2, ensure_ascii=False) + "\n")
    workflow = root / ".github/workflows/benchmark.yml"
    workflow.write_text(workflow.read_text() + (root / "tools/metrics-job.yml").read_text())
    checks = root / ".github/workflows/checker-tests.yml"
    checks.write_text(checks.read_text() + (root / "tools/tests-job.yml").read_text())


def check_generated(root=ROOT):
    root = Path(root)
    with tempfile.TemporaryDirectory() as tmp:
        copy = Path(tmp) / "source"
        shutil.copytree(root, copy, ignore=shutil.ignore_patterns(".git", "__pycache__", ".venv", ".yukon", ".work"))
        generate(copy)
        paths = ("benchmark.json", ".github/workflows/benchmark.yml",
                 ".github/workflows/checker-tests.yml", "README.md", "AGENTS.md",
                 "examples/circuitbreaker-presentation.json")
        changed = [p for p in paths if (root / p).read_bytes() != (copy / p).read_bytes()]
        if changed:
            raise SystemExit("Generated integration drift: " + ", ".join(changed))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    check_generated() if args.check else generate()
