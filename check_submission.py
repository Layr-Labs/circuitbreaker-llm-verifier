"""Check the shared submission contract and CircuitBreaker joint declaration."""

from pathlib import Path

from llm_verifier.config import load_config
from llm_verifier.files import snapshot
from llm_verifier.scaffold import check
from llm_verifier.submission import preflight_manifest

from export_metrics import accepted_bound, validate_contract


def check_submission(root):
    root = Path(root)
    result = check(root)
    effective = load_config(root, "verifier.json", "submission")
    config = effective["config"]
    validate_contract(config)
    files = snapshot(root / "submission", config["limits"])
    manifest, _ = preflight_manifest(files, config)
    accepted_bound(manifest, files)
    return result


if __name__ == "__main__":
    print(check_submission(Path(__file__).resolve().parent))
