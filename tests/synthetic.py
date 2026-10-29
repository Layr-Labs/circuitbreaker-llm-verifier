"""Deterministic test-only phase chains. No provider requests or credentials.

These deliberately synthetic unsigned records exercise the trusted-channel
consumer. They never represent real reviewer acceptance or a discovered theorem.
The real committee pipeline runs with its HTTP boundary replaced by fixed data;
production-shaped records below are test fixtures, not proof of model origin.
"""

from contextlib import contextmanager
import json
from pathlib import Path
from unittest.mock import patch
import uuid

from llm_verifier import pipeline, regions
from llm_verifier.artifacts import intake_snapshot
from llm_verifier.config import load_config
from llm_verifier.files import write_json
from llm_verifier.jsonutil import canonical
from llm_verifier.provider import HttpResponse, HttpTransport
from llm_verifier.runner import DockerRunner, output_record

ROOT = Path(__file__).resolve().parents[1]


@contextmanager
def scripted_http(reports, requests=None):
    """Replace the entire HTTP boundary, including background worker dispatch."""
    pending = list(reports)

    def request(body, *args):
        wire = json.loads(body)
        if requests is not None:
            requests.append(wire)
        if not pending:
            raise AssertionError("Unexpected synthetic provider call")
        report = pending.pop(0)
        # HttpTransport returns the normalized chat envelope for background
        # responses too. No subprocess or network operation runs in this test.
        return HttpResponse(200, canonical({
            "object": "chat.completion", "id": "synthetic-no-provider-request",
            "model": wire["model"], "usage": {}, "choices": [{
                "finish_reason": "stop", "message": {
                    "role": "assistant", "content": json.dumps(report),
                },
            }],
        }))

    with patch.dict("os.environ", {"OPENAI_API_KEY": "synthetic-never-network"}), patch.object(HttpTransport, "request", side_effect=request):
        yield
    if pending:
        raise AssertionError("Unused synthetic provider reports")


def source_plan(intake):
    return {"regions": [
        {"id": f"r{i}", "path": path, "start_line": 1,
         "end_line": len(regions.lines(source)), "purpose": "Synthetic complete source coverage."}
        for i, (path, source) in enumerate(sorted(intake["data"]["files"].items()), 1)
    ]}


def synthetic_review(intake, initial, first=None, followup=None, *, valid=True,
                     reference="source:proof.md", request_followup=False,
                     blocking=False, requests=None):
    """Run the installed decomposed protocol against deterministic reports."""
    config = intake["data"]["effective"]["config"]
    plan = source_plan(intake)
    reports = [] if first else [plan]
    requested = []
    if request_followup:
        fixture = json.loads((ROOT / "checker-tests.json").read_text())[0]
        requested = [{"check": "lp-certificate", "inputs": fixture["inputs"],
                      "reason": "Synthetic LP request."}]
    reports.append({"assessment": "Synthetic whole-proof report.", "findings": [], "requests": requested})
    for index, region in enumerate(plan["regions"]):
        findings = [{"id": "gap", "kind": "blocking", "reason": "Synthetic concern.",
                     "evidence_refs": [reference]}] if blocking and index == 0 else []
        reports.append({"assessment": "Synthetic region report.", "findings": findings, "requests": [],
                        "coverage": {"region_id": region["id"], "start_line": region["start_line"],
                                     "end_line": region["end_line"], "complete": True}})
    if not request_followup:
        resolutions = []
        if blocking:
            response = {"finding_id": f"region/{plan['regions'][0]['id']}/gap",
                        "reason": "Synthetic response to original evidence.", "evidence_refs": [reference]}
            reports.append({"responses": [response]})
            resolutions = [{**response, "outcome": "dismissed" if valid else "confirmed"}]
        reports.append({"decision": {"kind": "judgment", "requests": [], "judgment": {
            "valid": valid, "scores": None, "reason": "SYNTHETIC fixture; not a scientific verdict.",
            "evidence_refs": [reference], "review_standard": config["review_standard"],
        }}, "resolutions": resolutions})
    with scripted_http(reports, requests):
        return pipeline.review(intake, initial, first, followup)


def synthetic_dossier(root, bound, manifest, *, valid=True, simulation=False,
                      reference="source:proof.md", followup=False, blocking=False,
                      requests=None):
    write_json(root / "submission/submission.json", manifest)
    write_json(root / "submission/bound.json", bound)
    (root / "submission/proof.md").write_text("Synthetic transport fixture; no circuit theorem.\n")
    effective = load_config(root, "verifier.json", "submission")
    # Freeze identity and time to make complete fixture chains deterministic.
    with patch("time.time_ns", return_value=1700000000000000000), patch("uuid.uuid4", return_value=uuid.UUID(int=1)), patch("time.monotonic", return_value=1.0):
        intake = intake_snapshot(effective, root / "submission", simulation=simulation)
        initial = pipeline.execute(intake)  # No declared runs, no Docker calls.
        first = synthetic_review(intake, initial, valid=valid, reference=reference,
                                 request_followup=followup, blocking=blocking, requests=requests)
        if not followup:
            return pipeline.dossier(intake, initial, first)
        fixture = json.loads((ROOT / "checker-tests.json").read_text())[0]
        output = output_record("completed", 0, canonical({**fixture["expected"], "message": "Synthetic fixture"}), b"")
        with patch.object(DockerRunner, "run", return_value=output):
            executed = pipeline.execute(intake, first, initial)
        final = synthetic_review(intake, initial, first, executed, valid=valid,
                                 reference=reference, blocking=blocking, requests=requests)
        return pipeline.dossier(intake, initial, first, executed, final)
