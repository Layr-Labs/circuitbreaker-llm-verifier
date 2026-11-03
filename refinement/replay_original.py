"""Isolated, model-free replay of the frozen original producer and its wheel.

This is a subprocess helper for export_refinements.py, not a delivery recovery
route. Input hashes must be supplied through the organizer's trusted channel.
"""

import argparse
import importlib.util
import io
import json
from pathlib import Path
import socket
import sys
from types import ModuleType
from unittest.mock import patch
import zipfile

# -I omits both CWD and the script directory. Load only this checked-out helper.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from refinement.common import (digest, read_json, require, sha256, validate_origin)

SOURCE_MANIFEST = ROOT / "refinement/original-source-v1.json"
ORIGINAL_REF = "2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851"
VERIFIER_REF = "c3196332cbd0ea75a8f478673e04a78ae03e581d"
REPOSITORY = "Layr-Labs/circuitbreaker-llm-verifier"
ARTIFACT_NAMES = {
    "llmv-package", "llmv-intake", "llmv-initial-experiments", "llmv-initial-review",
    "verifier-dossier", "circuitbreaker-score", "benchmark-score",
    "verifier-accepted-results", "circuitbreaker-export-lineage",
}


def inside(root, relative):
    path = Path(relative)
    require(not path.is_absolute() and ".." not in path.parts, "unsafe input path")
    result = (root / path).resolve()
    require(result.is_relative_to(root.resolve()), "input path escapes root")
    return result


def verify_download(bundle, expected_sha):
    receipt_path = bundle / "download-receipt.json"
    require(sha256(receipt_path.read_bytes()) == expected_sha, "download receipt hash mismatch")
    receipt = read_json(receipt_path, numbers=True)
    run_path, listing_path = bundle / "run.json", bundle / "artifacts.json"
    require(sha256(run_path.read_bytes()) == receipt["runMetadataSha256"], "run metadata hash mismatch")
    require(sha256(listing_path.read_bytes()) == receipt["artifactListingSha256"], "artifact listing hash mismatch")
    run, listing = read_json(run_path, numbers=True), read_json(listing_path, numbers=True)
    require((run["status"], run["conclusion"], run["event"]) ==
            ("completed", "success", "workflow_dispatch"), "original run was not successful dispatch")
    require(run["repository"]["full_name"] == receipt["repository"] == REPOSITORY,
            "original repository mismatch")
    require(run["id"] == receipt["runId"] and run["head_sha"] == receipt["headSha"] and
            run["run_attempt"] == receipt["runAttempt"] == 1, "original run identity mismatch")
    require({x["name"] for x in receipt["artifacts"]} == ARTIFACT_NAMES and
            len(receipt["artifacts"]) == len(ARTIFACT_NAMES), "complete original artifact set required")
    listed = {item["id"]: item for item in listing["artifacts"]}
    require(len(listed) == len(listing["artifacts"]), "duplicate artifact id")
    archive_hashes = {}
    for item in receipt["artifacts"]:
        meta = listed[item["id"]]
        raw = inside(bundle, item["archive"]).read_bytes()
        require(sha256(raw) == item["archiveSha256"] == meta["digest"].removeprefix("sha256:"),
                "archive digest mismatch")
        require(len(raw) == item["archiveBytes"] == meta["size_in_bytes"], "archive length mismatch")
        require(meta["name"] == item["name"] and meta["workflow_run"]["id"] == run["id"] and
                meta["workflow_run"]["head_sha"] == run["head_sha"], "artifact run identity mismatch")
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            members = [info.filename for info in archive.infolist() if not info.is_dir()]
            expected = [file["path"] for file in item["files"]]
            require(len(members) == len(set(members)) and len(expected) == len(set(expected)) and
                    set(members) == set(expected), "archive member mismatch")
            for file in item["files"]:
                data = archive.read(file["path"])
                require(len(data) == file["bytes"] and sha256(data) == file["sha256"],
                        "archive file hash mismatch")
                require(data == inside(bundle / item["name"], file["path"]).read_bytes(),
                        "extracted artifact differs from original archive")
        archive_hashes[item["name"]] = item["archiveSha256"]
    return receipt, archive_hashes


def verify_source(source):
    manifest = read_json(SOURCE_MANIFEST)
    require(manifest["sourceRef"] == ORIGINAL_REF, "unsupported original producer source")
    for name, expected in manifest["files"].items():
        require(sha256(inside(source, name).read_bytes()) == expected, "original source mismatch: " + name)
    return manifest


def verify_wheel(bundle):
    # The interpreter is itself a trusted organizer input. Before importing the
    # verifier, compare every installed package file with the authenticated wheel.
    spec = importlib.util.find_spec("llm_verifier")
    require(spec is not None and spec.origin is not None, "original verifier is not installed")
    installed = Path(spec.origin).parent
    wheels = list((bundle / "llmv-package").glob("*.whl"))
    require(len(wheels) == 1, "one original verifier wheel required")
    with zipfile.ZipFile(wheels[0]) as archive:
        members = [name for name in archive.namelist()
                   if name.startswith("llm_verifier/") and not name.endswith("/")]
        require(len(members) == len(set(members)), "duplicate wheel member")
        actual = {str(p.relative_to(installed)) for p in installed.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
        require(actual == {n.removeprefix("llm_verifier/") for n in members}, "installed wheel file set mismatch")
        for name in members:
            require(archive.read(name) == inside(installed, name.removeprefix("llm_verifier/")).read_bytes(),
                    "installed verifier differs from original wheel")
    return sha256(wheels[0].read_bytes())


def load_original(source):
    # Execute these verified source bytes, never stale bytecode or arbitrary
    # modules from the original checkout. No candidate file is imported.
    for name in ("expressions", "export_metrics"):
        module = ModuleType(name)
        module.__file__ = str(source / (name + ".py"))
        sys.modules[name] = module
        exec(compile(Path(module.__file__).read_bytes(), module.__file__, "exec"), module.__dict__)
    return sys.modules["export_metrics"], sys.modules["expressions"]


def blocked(*args, **kwargs):
    raise RuntimeError("network, model review and experiment execution forbidden during replay")


def replay(source, bundle, public, output):
    from llm_verifier import pipeline
    from llm_verifier.artifacts import implementation, verify_intake
    from llm_verifier.config import load_config
    from llm_verifier.files import write_json
    from llm_verifier.jsonutil import canonical as original_canonical, digest as original_digest
    from llm_verifier.provider import HttpTransport
    from llm_verifier.runner import DockerRunner

    exporter, expressions = load_original(source)
    intake = read_json(bundle / "llmv-intake/intake.json", numbers=True)
    dossier = read_json(bundle / "verifier-dossier/dossier.json", numbers=True)
    with patch.object(HttpTransport, "request", side_effect=blocked), \
         patch.object(DockerRunner, "run", side_effect=blocked), \
         patch.object(pipeline, "review", side_effect=blocked), \
         patch.object(pipeline, "execute", side_effect=blocked), \
         patch.object(socket.socket, "connect", side_effect=blocked), \
         patch.object(socket, "create_connection", side_effect=blocked):
        verify_intake(intake)
        effective = load_config(source, "verifier.json", "submission",
                                evaluation=intake["data"]["effective"]["evaluation"])
        require(original_canonical(effective) == original_canonical(intake["data"]["effective"]),
                "original effective config mismatch")
        phases = dossier["data"]["phases"]
        require(len(phases) == 3, "only the original three-phase archive contract is supported")
        require(original_canonical(phases[0]) == original_canonical(intake), "intake phase mismatch")
        for phase, path in ((phases[1], "llmv-initial-experiments/initial.json"),
                            (phases[2], "llmv-initial-review/review.json")):
            require(original_canonical(phase) == original_canonical(read_json(bundle / path, numbers=True)),
                    "standalone phase artifact mismatch")
        require(not dossier["simulation"] and dossier["data"]["judgment"]["valid"] is True,
                "original real accepted judgment required")
        rebuilt = pipeline.dossier(*phases, finalized_at_ns=dossier["data"]["finalized_at_ns"])
        require(original_canonical(rebuilt) == original_canonical(dossier), "original dossier mismatch")
        write_json(output / "dossier.json", rebuilt)
        exporter.export(bundle / "verifier-dossier/dossier.json", output / "score.json",
                        intake_path=bundle / "llmv-intake/intake.json")
        write_json(output / "benchmark-score.json", rebuilt["data"]["yukon_projection"])
        write_json(output / "accepted-results.json", pipeline.export_results(rebuilt))
        originals = {"dossier.json": "verifier-dossier/dossier.json",
                     "score.json": "circuitbreaker-score/.yukon/score.json",
                     "benchmark-score.json": "benchmark-score/.yukon/score.json",
                     "accepted-results.json": "verifier-accepted-results/accepted-results.json"}
        for name, path in originals.items():
            require((output / name).read_bytes() == (bundle / path).read_bytes(), "original byte replay mismatch: " + name)
        result = read_json(output / "score.json", numbers=True)
        packet = public["evaluation"]
        require(original_canonical(result["metrics"]) == original_canonical(public["officialMetrics"]),
                "selected platform metrics differ from original output")
        require(result["scores"] == packet["scores"]["values"] and packet["validity"]["valid"] is True,
                "selected original scores/validity mismatch")
        require(packet["validity"]["reason"] == dossier["data"]["judgment"]["reason"],
                "selected validity reason mismatch")
        require(public["status"] == "review" and packet["human"] is None,
                "bounded path requires unchanged pending-human original status")
        context = intake["data"]["context"]
        require(context["verifier_revision"] == VERIFIER_REF, "original verifier revision mismatch")
        for field in ("scores", "validity"):
            evidence = packet[field]["evidence"]
            require(evidence["jobId"] == context["accepted_results"]["job_id"] and
                    evidence["evaluatorId"] == context["accepted_results"]["evaluation_id"] and
                    evidence["commitSha"] == context["candidate_revision"] and
                    evidence["sourceRef"] == ORIGINAL_REF, "selected evidence identity mismatch")
        require(context["accepted_results"]["submission_id"] == public["id"] and
                context["accepted_results"]["benchmark_id"] == public["benchmarkId"], "origin identity mismatch")
        # Yukon candidateKey identifies the platform admission snapshot; the
        # verifier submission digest hashes its own files envelope. These are
        # distinct namespaces, each pinned by its authenticated input.
        manifest, files = intake["data"]["manifest"], intake["data"]["files"]
        bound = exporter.accepted_bound(manifest, files)
        require(bound is not None and bound["version"] == expressions.VERSION and
                bound["guarantee"]["kind"] == "polynomial_log" and bound["projection"]["kind"] == "finite",
                "finite original v1 polynomial declaration required")
        require(original_canonical(manifest["claims"][:2]) == original_canonical(
            expressions.claims_for(bound, manifest["claims"][0]["evidence"])), "original full atomic claims mismatch")
        for name, file in files.items():
            require(sha256(file["text"].encode()) == file["sha256"] and
                    inside(source / "submission", name).read_bytes() == file["text"].encode(),
                    "original submitted source file mismatch")
        rich, frontier = result["metrics"]["circuitbreaker"], result["metrics"]["circuitbreaker_frontier"]
        require(rich["status"] == "accepted_provisional" and rich["bound"] == bound and
                rich["selected_review_sha256"] == phases[2]["sha256"], "selected original review mismatch")
        for claim in manifest["claims"][:2]:
            require(json.loads(claim["statement"].split("Declaration: ", 1)[1]) == bound,
                    "entire declaration must be in both accepted statements")
        sources = [original_digest({"source": {"claim": i}, "files": original_digest(files)}) for i in (0, 1)]
        require([r["source_sha256"] for r in rich["score_records"]] == sources, "atomic claim source mismatch")
        code = implementation()
        lineage = read_json(bundle / "circuitbreaker-export-lineage/recovery.json", numbers=True)
        require(lineage["scoreSha256"] == sha256((output / "score.json").read_bytes()), "original score lineage mismatch")
        expected = {"intakeSha256": intake["sha256"], "dossierSha256": dossier["sha256"],
                    "verifierSha256": code["content_sha256"], "configSha256": intake["bindings"]["config"],
                    "submissionSha256": intake["bindings"]["submission"], "candidateSha": context["candidate_revision"],
                    "verifierRevision": VERIFIER_REF}
        require(all(lineage["identity"][k] == v for k, v in expected.items()), "original lineage identity mismatch")
        origin = {"benchmarkId": public["benchmarkId"], "submissionId": public["id"], "packetId": packet["id"],
                  "candidateKey": packet["candidateKey"], "evaluationId": packet["evaluationId"],
                  "validityVersionId": packet["validityVersionId"], "scoreEvidence": packet["scores"]["evidence"],
                  "validityEvidence": packet["validity"]["evidence"], "originalEnvelopeSha256": digest(frontier),
                  "originalRichSha256": digest(rich), "intakeSha256": intake["sha256"],
                  "dossierSha256": dossier["sha256"], "sourceSha256": sources[0]}
        validate_origin(origin)
        return {"origin": origin, "originalBound": bound, "lineage": expected,
                "selectedReviewSha256": phases[2]["sha256"], "fullDeclarationSha256": digest(bound),
                "fullClaimsSha256": original_digest(manifest["claims"]),
                "submittedFiles": {name: file["sha256"] for name, file in files.items()},
                "originalStatus": {"submission": public["status"], "human": None, "metric": rich["status"]},
                "byteIdenticalOutputs": {name: sha256((output / name).read_bytes()) for name in originals}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    request = read_json(args.request)
    source, bundle, output = (Path(request["originalSource"]), Path(request["bundle"]), Path(args.out))
    output.mkdir()  # Never overwrite an existing artifact directory.
    download, archives = verify_download(bundle, request["downloadReceiptSha256"])
    require(sha256(Path(request["platform"]).read_bytes()) == request["platformSha256"], "platform snapshot hash mismatch")
    public = read_json(request["platform"], numbers=True)["submission"]
    intake_context = read_json(bundle / "llmv-intake/intake.json", numbers=True)["data"]["context"]
    require(intake_context["run_id"] == f"{download['runId']}:{download['runAttempt']}",
            "intake run differs from authenticated artifact run")
    source_manifest = verify_source(source)
    wheel_sha = verify_wheel(bundle)
    report = replay(source, bundle, public, output)
    require(download["submissionId"] == report["origin"]["submissionId"] and
            download["headSha"] == report["origin"]["scoreEvidence"]["commitSha"], "download origin mismatch")
    # Recheck immutable source and downloaded inputs after reconstruction.
    verify_source(source)
    verify_download(bundle, request["downloadReceiptSha256"])
    require(sha256(Path(request["platform"]).read_bytes()) == request["platformSha256"], "platform snapshot changed")
    verify_wheel(bundle)
    report.update(originalProducer=source_manifest, originalWheelSha256=wheel_sha,
                  archiveSha256=archives, downloadReceiptSha256=request["downloadReceiptSha256"],
                  platformSha256=request["platformSha256"], networkCalls="0", modelCalls="0", experimentExecutions="0")
    from refinement.common import canonical
    (output / "replay.json").write_bytes(canonical(report) + b"\n")


if __name__ == "__main__":
    main()
