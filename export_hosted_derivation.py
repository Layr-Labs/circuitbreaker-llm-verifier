"""Export only an exact current-run derivation from the trusted artifact channel."""

import argparse
import os
from pathlib import Path

from llm_verifier.files import clear_output, distinct_paths, read_json, write_json
from llm_verifier.jsonutil import loads

from export_hosted import file_hash, positive_integer, require, sha256
from score_derivation import projection
from tools.derivation_profile import read_profile

ROOT = Path(__file__).resolve().parent
NAMES = ("llmv-derivation-package", "llmv-derivation-intake", "verifier-score-derivation")


def export(intake_path, result_path, wheel_path, output, lineage_path, transfer_receipt, *, env=None, root=ROOT):
    inputs = [Path(wheel_path), Path(intake_path), Path(result_path)]
    distinct_paths([*inputs, transfer_receipt], [output, lineage_path])
    clear_output(output)
    clear_output(lineage_path)
    try:
        env = os.environ if env is None else env
        profile = read_profile(root)
        run_id = positive_integer(env["GITHUB_RUN_ID"], "run_id")
        delivery = positive_integer(env["GITHUB_RUN_ATTEMPT"], "delivery_attempt")
        review = positive_integer(env["LLMV_REVIEW_RUN_ATTEMPT"], "review_attempt")
        require(positive_integer(env["LLMV_REVIEW_RUN_ID"], "review_run_id") == run_id, "review_run_id")
        require(review <= delivery, "review_attempt")
        expected = sha256(env["LLMV_EXPECTED_INTAKE"], "expected_intake")
        descriptors = [loads(env[key].encode()) for key in
                       ("LLMV_PACKAGE_ARTIFACT", "LLMV_INTAKE_ARTIFACT", "LLMV_DERIVATION_ARTIFACT")]
        attempts, ids = [], []
        for value, name in zip(descriptors, NAMES):
            require(type(value) is dict and set(value) == {"id", "name", "sha256", "run_id", "run_attempt"}, "artifact_fields")
            require(value["name"] == name, "artifact_name")
            require(positive_integer(value["run_id"], "artifact_run_id") == run_id, "artifact_run_id")
            attempts.append(positive_integer(value["run_attempt"], "artifact_attempt"))
            ids.append(positive_integer(value["id"], "artifact_id"))
            sha256(value["sha256"], "artifact_digest")
        require(len(set(ids)) == 3 and attempts[0] == attempts[1] <= review <= attempts[2] <= delivery, "artifact_order")
        require(inputs[0].suffix == ".whl" and inputs[1].name == "intake.json" and inputs[2].name == "derivation.json", "input_paths")
        hashes = [file_hash(path) for path in inputs]
        receipt = read_json(transfer_receipt)
        expected_receipt = {
            "version": 1, "run_id": str(run_id), "delivery_attempt": str(delivery),
            "artifacts": [{"artifact": descriptor, "path": destination,
                           "files": {path.name: digest}}
                          for descriptor, destination, path, digest in
                          zip(descriptors, ("wheel", "intake", "derivation"), inputs, hashes)],
        }
        require(type(receipt.get("version")) is int and receipt == expected_receipt, "transfer_receipt")
        intake, result = read_json(intake_path), read_json(result_path)
        require(intake["sha256"] == expected, "intake_sha256")
        context = intake["data"]["context"]
        require(context.get("kind") == "github-actions", "source_kind")
        require(context.get("candidate_revision") == env["GITHUB_SHA"], "source_revision")
        require(context.get("repository") == env["GITHUB_REPOSITORY"], "source_repository")
        require(context.get("verifier_revision") == profile["verifier_revision"], "verifier_revision")
        require(context.get("run_id") == f"{run_id}:{attempts[1]}", "intake_run")
        require(intake["data"]["config_file_sha256"] == profile["files"]["derivation.json"], "profile_config")
        # The artifact channel supplies authenticity; verification below binds
        # the exact intake, decision, citations and non-simulation result.
        score = projection(result, expected_intake=expected)
        require(hashes == [file_hash(path) for path in inputs], "input_changed")
        write_json(output, score)
        write_json(lineage_path, {
            "version": "circuitbreaker-derivation-export-lineage-v1",
            "run_id": str(run_id), "delivery_attempt": str(delivery), "review_attempt": str(review),
            "candidate_revision": context["candidate_revision"], "verifier_revision": profile["verifier_revision"],
            "artifacts": descriptors, "input_sha256": dict(zip(NAMES, hashes)),
            "intake_sha256": intake["sha256"], "derivation_sha256": result["sha256"],
            "derived_from": score["derived_from"], "score_sha256": file_hash(output),
        })
    except Exception:
        clear_output(output)
        clear_output(lineage_path)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("intake", "derivation", "wheel", "out", "lineage", "transfer-receipt"):
        parser.add_argument("--" + field, required=True)
    args = parser.parse_args()
    export(args.intake, args.derivation, args.wheel, args.out, args.lineage, args.transfer_receipt)


if __name__ == "__main__":
    main()
