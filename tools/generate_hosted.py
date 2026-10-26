"""Regenerate upstream callers, then apply the reviewed CircuitBreaker overlay."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path
import subprocess
import sys

VERIFIER_REVISION = "8c1f098875fae0750c29ecfd641f7c9fc53a0f41"
ROOT = Path(__file__).resolve().parents[1]


def generate(root=ROOT):
    root = Path(root)
    subprocess.run([
        sys.executable, "-m", "llm_verifier", "github", "--root", str(root),
        "--verifier-revision", VERIFIER_REVISION,
    ], check=True)
    path = root / "benchmark.json"
    benchmark = json.loads(path.read_text())
    track = benchmark["tracks"][0]
    track["runner"]["scoreArtifactName"] = "circuitbreaker-score"
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
                 ".github/workflows/checker-tests.yml", "README.md", "AGENTS.md")
        changed = [p for p in paths if (root / p).read_bytes() != (copy / p).read_bytes()]
        if changed:
            raise SystemExit("Generated integration drift: " + ", ".join(changed))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    check_generated() if args.check else generate()
