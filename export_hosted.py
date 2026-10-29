"""Export from trusted workflow artifacts and record deterministic-delivery lineage.

The workflow supplies authenticated artifact descriptors. These hashes record
integrity, not authenticity; this module must not accept candidate-owned inputs.
"""

import argparse
import hashlib
import os
from pathlib import Path
import re

from llm_verifier.errors import VerifierError
from llm_verifier.files import clear_output, distinct_paths, read_json, write_json
from llm_verifier.jsonutil import loads

from export_metrics import export
from tools.generate_hosted import VERIFIER_REVISION

ARTIFACT_NAMES = ("llmv-package", "llmv-intake", "verifier-dossier")
MAX_SAFE_INTEGER = 2**53 - 1


def require(condition, field):
    if not condition:
        raise VerifierError("circuitbreaker_export_lineage_" + field)


def positive_integer(value, field):
    require(type(value) is str and re.fullmatch(r"[1-9][0-9]{0,15}", value), field)
    number = int(value)
    require(number <= MAX_SAFE_INTEGER, field)
    return number


def sha256(value, field):
    require(type(value) is str and re.fullmatch(r"[a-f0-9]{64}", value), field)
    return value


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environment(env):
    """Read only explicitly named, non-secret workflow identity fields."""
    try:
        result = {
            "run_id": positive_integer(env["GITHUB_RUN_ID"], "run_id"),
            "delivery_attempt": positive_integer(env["GITHUB_RUN_ATTEMPT"], "delivery_attempt"),
            "head_sha": env["GITHUB_SHA"],
            "review_run_id": positive_integer(env["LLMV_REVIEW_RUN_ID"], "review_run_id"),
            "review_attempt": positive_integer(env["LLMV_REVIEW_RUN_ATTEMPT"], "review_attempt"),
            "expected_intake": sha256(env["LLMV_EXPECTED_INTAKE"], "expected_intake"),
            "artifacts": [loads(env[key].encode()) for key in
                          ("LLMV_PACKAGE_ARTIFACT", "LLMV_INTAKE_ARTIFACT", "LLMV_DOSSIER_ARTIFACT")],
        }
    except KeyError as exc:
        raise VerifierError("circuitbreaker_export_lineage_missing_environment") from exc
    require(re.fullmatch(r"[a-f0-9]{40}", result["head_sha"]), "head_sha")
    require(result["review_run_id"] == result["run_id"], "review_run_id")
    require(result["review_attempt"] <= result["delivery_attempt"], "review_attempt")
    return result


def artifact_records(descriptors, run_id, delivery_attempt):
    require(type(descriptors) is list and len(descriptors) == len(ARTIFACT_NAMES), "artifacts")
    records = []
    for descriptor, name in zip(descriptors, ARTIFACT_NAMES):
        require(type(descriptor) is dict and set(descriptor) ==
                {"id", "name", "sha256", "run_id", "run_attempt"}, "artifact_fields")
        require(descriptor["name"] == name, "artifact_name")
        require(positive_integer(descriptor["run_id"], "artifact_run_id") == run_id, "artifact_run_id")
        attempt = positive_integer(descriptor["run_attempt"], "artifact_run_attempt")
        require(attempt <= delivery_attempt, "artifact_run_attempt")
        records.append({"id": positive_integer(descriptor["id"], "artifact_id"), "name": name,
                        "sha256": sha256(descriptor["sha256"], "artifact_sha256"), "runAttempt": attempt})
    require(len({record["id"] for record in records}) == len(records), "artifact_ids")
    return records


def export_hosted(intake_path, dossier_path, wheel_path, output, lineage_path, transfer_receipt, *, env=None):
    """Run the ordinary exporter, then bind its exact bytes to original inputs.

    Original and delivery attempts intentionally remain separate. The platform
    applies its stricter first-attempt recovery policy when importing a rerun.
    """
    inputs = [Path(wheel_path), Path(intake_path), Path(dossier_path)]
    distinct_paths([*inputs, transfer_receipt], [output, lineage_path])
    clear_output(output)
    clear_output(lineage_path)
    try:
        context = environment(os.environ if env is None else env)
        artifacts = artifact_records(context["artifacts"], context["run_id"], context["delivery_attempt"])
        require(inputs[0].suffix == ".whl" and inputs[1].name == "intake.json" and
                inputs[2].name == "dossier.json", "input_paths")
        input_files = [{"artifactName": name, "path": path.name, "sha256": file_hash(path)}
                       for name, path in zip(ARTIFACT_NAMES, inputs)]
        receipt = read_json(transfer_receipt)
        require(type(receipt) is dict and set(receipt) ==
                {"version", "run_id", "delivery_attempt", "artifacts"}, "transfer_receipt")
        require(type(receipt["version"]) is int and receipt["version"] == 1 and
                receipt["run_id"] == str(context["run_id"]) and
                receipt["delivery_attempt"] == str(context["delivery_attempt"]), "transfer_run")
        expected_transfers = [
            {"artifact": descriptor, "path": destination, "files": {item["path"]: item["sha256"]}}
            for descriptor, destination, item in zip(context["artifacts"], ("wheel", "intake", "dossier"), input_files)
        ]
        require(receipt["artifacts"] == expected_transfers, "transfer_artifacts")
        intake, dossier = read_json(intake_path), read_json(dossier_path)
        require(intake["sha256"] == context["expected_intake"], "intake_sha256")
        source = intake["data"]["context"]
        require(source.get("kind") == "github-actions", "source_kind")
        require(source.get("candidate_revision") == context["head_sha"], "source_revision")
        require(source.get("verifier_revision") == VERIFIER_REVISION, "verifier_revision")
        require(source.get("run_id") == f'{context["run_id"]}:{artifacts[1]["runAttempt"]}', "intake_run")
        require(artifacts[0]["runAttempt"] == artifacts[1]["runAttempt"] <=
                context["review_attempt"] <= artifacts[2]["runAttempt"], "artifact_order")
        # This is the existing exporter, including full dossier replay, intake
        # equality and simulation rejection. No provider or experiment is run.
        export(dossier_path, output, intake_path=intake_path)
        require(all(file_hash(path) == item["sha256"] for path, item in zip(inputs, input_files)), "input_changed")
        bindings = intake["bindings"]
        identity = {
            "evaluationId": bindings["evaluation_id"],
            "configSha256": bindings["config"],
            "candidateSha": source["candidate_revision"],
            "intakeSha256": intake["sha256"],
            "dossierSha256": dossier["sha256"],
            "submissionSha256": bindings["submission"],
            "verifierSha256": bindings["verifier"],
            "verifierRevision": source["verifier_revision"],
        }
        write_json(lineage_path, {
            "version": 1, "runId": context["run_id"], "headSha": context["head_sha"],
            "reviewAttempt": context["review_attempt"], "deliveryAttempt": context["delivery_attempt"],
            "artifacts": artifacts, "inputFiles": input_files, "identity": identity,
            "scoreSha256": file_hash(output),
        })
    except Exception:
        clear_output(output)
        clear_output(lineage_path)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("intake", "dossier", "wheel", "out", "lineage", "transfer-receipt"):
        parser.add_argument("--" + field, required=True)
    args = parser.parse_args()
    export_hosted(args.intake, args.dossier, args.wheel, args.out, args.lineage, args.transfer_receipt)


if __name__ == "__main__":
    main()
