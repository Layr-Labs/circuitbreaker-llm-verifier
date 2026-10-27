"""Deterministic test-only phase chains. No provider requests or credentials.

These deliberately synthetic unsigned records exercise the trusted-channel
consumer. They never represent real reviewer acceptance or a discovered theorem.
"""

from unittest.mock import patch
import uuid

from llm_verifier import pipeline
from llm_verifier.artifacts import intake_snapshot, seal
from llm_verifier.config import load_config
from llm_verifier.files import write_json
from llm_verifier.jsonutil import digest
from llm_verifier.provider import provider_endpoint


def synthetic_dossier(root, bound, manifest, *, valid=True, simulation=False,
                      reference="source:proof.md"):
    write_json(root / "submission/submission.json", manifest)
    write_json(root / "submission/bound.json", bound)
    (root / "submission/proof.md").write_text("Synthetic transport fixture; no circuit theorem.\n")
    effective = load_config(root, "verifier.json", "submission")
    # Freeze identity and time to make complete fixture chains deterministic.
    with patch("time.time_ns", return_value=1700000000000000000), patch("uuid.uuid4", return_value=uuid.UUID(int=1)), patch("time.monotonic", return_value=1.0):
        intake = intake_snapshot(effective, root / "submission", simulation=simulation)
        initial = pipeline.execute(intake)  # No declared runs, no Docker calls.
        messages, schema, _ = pipeline.review_material(intake, initial)
        config = effective["config"]
        decision = {"kind": "judgment", "requests": [], "judgment": {
            "valid": valid, "scores": None, "reason": "SYNTHETIC fixture; not a scientific verdict.",
            "evidence_refs": [reference], "review_standard": config["review_standard"],
        }}
        provenance = {
            "provider": "openai", "endpoint": provider_endpoint(config["provider"]),
            "configured_model": config["provider"]["model"], "returned_model": config["provider"]["model"],
            "request_id": "synthetic-no-provider-request", "usage": {},
            "prompt_sha256": digest(messages), "schema_sha256": digest(schema),
            "provider_config_sha256": digest(config["provider"]), "attempts": 1, "elapsed_seconds": 0,
        }
        review = seal("initial_review", intake["bindings"], {"decision": decision, "provenance": provenance},
                      [intake["sha256"], initial["sha256"]], simulation)
        return pipeline.dossier(intake, initial, review)
