"""Build a protected full-history derivation plan from organizer-staged evidence.

Inputs must be captured by an authenticated organizer. No API/model calls are
made here. Runtime admission independently compares this plan with Yukon's issued
retained-result descriptor and verifies the original artifact/file digests.
"""

import argparse
from copy import deepcopy
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from llm_verifier import derivation, retained_result
from llm_verifier.derivation_transport import read_artifact
from llm_verifier.files import read_json, write_json
from llm_verifier.jsonutil import digest, loads, sha
from llm_verifier.prior import validate_bundle
from llm_verifier.artifacts import verify_envelope

from export_metrics import accepted_bound, validate_contract
from expressions import VERSION


def evidence(value):
    names = {"result_id": "resultId", "evaluation_id": "evaluatorId",
             "policy_version_id": "policyVersionId", "source_ref": "sourceRef",
             "job_id": "jobId", "commit_sha": "commitSha"}
    return {key: value[source] for key, source in names.items()}


def selected_origin(row):
    packet = row["evaluation"]
    result = {
        "version": "yukon-retained-result-v1", "benchmark_id": row["benchmarkId"],
        "submission_id": row["id"], "candidate_key": packet["candidateKey"],
        "packet_id": packet["id"], "score": evidence(packet["scores"]["evidence"]),
        "validity": {**evidence(packet["validity"]["evidence"]), "valid": packet["validity"]["valid"]},
    }
    retained_result.origin(result)
    return result


def build(snapshot, artifact_index, config_raw, archive_root):
    derivation.normalize_config(loads(config_raw))
    rows = snapshot["submissions"]
    if type(rows) is not list or len(rows) > 10000:
        raise ValueError("bounded submission list required")
    if set(artifact_index) != {"version", "records"} \
            or artifact_index["version"] != "circuitbreaker-derivation-artifacts-v1":
        raise ValueError("invalid organizer artifact index")
    indexed = {}
    for item in artifact_index["records"]:
        if set(item) != {"submission_id", "source_files_sha256", "artifacts"}:
            raise ValueError("invalid organizer artifact record")
        if item["submission_id"] in indexed:
            raise ValueError("duplicate organizer artifact record")
        indexed[item["submission_id"]] = item
    entries, seen, ordinary, excluded = [], set(), [], []
    for row in rows:
        if row["id"] in seen:
            raise ValueError("duplicate submission in snapshot")
        seen.add(row["id"])
        if row.get("deletedAt") or row["status"] == "cancelled":
            excluded.append(row["id"])
            continue
        packet = row.get("evaluation")
        if not packet or not packet.get("scores") or not packet.get("validity") \
                or not packet.get("publishedAt"):
            ordinary.append(row["id"])
            continue
        origin = selected_origin(row)
        if row["id"] not in indexed:
            raise ValueError("missing organizer artifact evidence for " + row["id"])
        record = indexed[row["id"]]
        entry = {"origin": origin, "action": "derive",
                 "reason": "Transcribe the retained full CircuitBreaker declaration under the scoring policy.",
                 "source_files_sha256": record["source_files_sha256"],
                 "artifacts": deepcopy(record["artifacts"])}
        if not origin["validity"]["valid"]:
            entry.update(action="not_applicable", artifacts=None,
                         reason="The selected validity result is false; preserve rejection and emit null scores.")
            if record["artifacts"] is not None:
                raise ValueError("invalid retained result must not be presented as accepted artifacts")
        else:
            derivation.validate_entry(entry)
            intake_raw = read_artifact(archive_root, entry["artifacts"]["intake"])
            accepted_raw = read_artifact(archive_root, entry["artifacts"]["accepted_results"])
            intake, accepted = loads(intake_raw), loads(accepted_raw)
            validate_bundle(accepted)
            source = verify_envelope(intake, "intake", parents=[])
            if intake["simulation"] or accepted["source"]["intake_sha256"] != intake["sha256"] \
                    or accepted["source"]["bindings"] != intake["bindings"]:
                raise ValueError("original accepted bundle/intake mismatch")
            if source["context"]["candidate_revision"] != origin["validity"]["commit_sha"] \
                    or source["context"]["run_id"] != (entry["artifacts"]["intake"]["run_id"] + ":" + entry["artifacts"]["intake"]["run_attempt"]):
                raise ValueError("original validity execution mismatch")
            if digest(source["files"]) != entry["source_files_sha256"]:
                raise ValueError("original candidate file digest mismatch")
            legacy = validate_contract(source["effective"]["config"])
            bound = accepted_bound(source["manifest"], source["files"])
            if bound is not None and legacy and bound["version"] != VERSION:
                raise ValueError("legacy review cannot cover a v2 declaration")
            if bound is None:
                entry.update(action="not_applicable",
                             reason="The retained full claims are unscored supporting lemmas; preserve null scores.")
        derivation.validate_entry(entry)
        entries.append(entry)
    if set(indexed) != {entry["origin"]["submission_id"] for entry in entries}:
        raise ValueError("artifact index has records outside the retained-result plan")
    if not entries:
        raise ValueError("no retained score-only records in snapshot")
    plan = {"version": derivation.PLAN_VERSION, "config_sha256": sha(config_raw), "entries": entries}
    report = {"version": "circuitbreaker-derivation-plan-report-v1", "plan_sha256": digest(plan),
              "score_only": len(entries), "derive": sum(e["action"] == "derive" for e in entries),
              "not_applicable": sum(e["action"] == "not_applicable" for e in entries),
              "ordinary_submission_ids": ordinary, "excluded_submission_ids": excluded,
              "baseline": "separate ordinary scores-and-validity execution",
              "model_calls": 0, "published": False}
    return plan, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "artifact-index", "config", "archive-root", "out"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    output = Path(args.out)
    if output.exists() or not output.parent.is_dir():
        raise ValueError("output must be a new directory beneath an existing parent")
    plan, report = build(read_json(args.snapshot), read_json(args.artifact_index),
                         Path(args.config).read_bytes(), args.archive_root)
    output.mkdir()
    write_json(output / "plan.json", plan)
    write_json(output / "report.json", report)
    print("Prepared protected plan; no workflow, review or publication executed.")


if __name__ == "__main__":
    main()
