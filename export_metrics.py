"""Export trusted accepted CircuitBreaker guarantees, including null-score endpoints.

Input must come from the organizer's same-run artifact channel. Digest checks
provide integrity, not authentication; never import candidate-created dossiers.
"""

import argparse

from llm_verifier import pipeline
from llm_verifier.artifacts import verify_envelope
from llm_verifier.errors import VerifierError
from llm_verifier.files import clear_output, distinct_paths, read_json, write_json
from llm_verifier.jsonutil import canonical, digest, loads
from llm_verifier.scoring import score_value

from expressions import (CONVENTION, POLICIES, REVIEW, VERSION, ExpressionError,
                         claims_for, ratio, validate_bound)


def validate_contract(config):
    objectives = config["objectives"]
    if (config["id"] != "circuitbreaker"
            or config["review_standard"] != REVIEW
            or config["policies"] != POLICIES
            or config.get("acceptance_policy", "atomic-claims-v1") != "atomic-claims-v1"
            or config.get("manifest_preprocessor") is not None
            or "structured_claims" in config
            or not config["allow_unscored"]
            or [o["name"] for o in objectives] != ["gate_coefficient", "gate_exponent"]
            or [o["direction"] for o in objectives] != ["+", "+"]
            or [o.get("score_source") for o in objectives] != [{"claim": 0}, {"claim": 1}]
            or [o.get("score_minimum") for o in objectives] != [ratio(0), ratio(1)]
            or any("score_maximum" in o for o in objectives)):
        raise VerifierError("circuitbreaker_objective_contract")


def accepted_bound(manifest, files):
    claims = manifest["claims"]
    # Only string-only declarations are unscored lemmas. A malformed or partial
    # main claim must not be silently converted into an accepted lemma.
    if all(type(c) is str for c in claims):
        return None
    if (len(claims) < 2 or any(type(c) is not dict for c in claims[:2])
            or any(type(c) is not str for c in claims[2:])
            or "bound.json" not in files):
        raise ExpressionError("two joint main claims required")
    bound = loads(files["bound.json"]["text"].encode())
    validate_bound(bound)
    evidence = claims[0]["evidence"]
    if (evidence != sorted(set(evidence))
            or not {"bound.json", "proof.md"} <= set(evidence)
            or any(p not in files for p in evidence)
            or canonical(claims[:2]) != canonical(claims_for(bound, evidence))):
        raise ExpressionError("joint statement/score/evidence mismatch")
    return bound


def projection(record, *, expected_intake=None):
    data = verify_envelope(record, "dossier")
    phases = data["phases"]
    if type(phases) is not list or len(phases) not in (3, 5):
        raise VerifierError("circuitbreaker_dossier_phases")
    rebuilt = pipeline.dossier(*phases, finalized_at_ns=data["finalized_at_ns"])
    if canonical(rebuilt) != canonical(record):
        raise VerifierError("circuitbreaker_dossier_mismatch")
    if record["simulation"]:
        raise VerifierError("simulation_cannot_publish_score")
    intake = phases[0]
    if expected_intake is not None and intake["sha256"] != expected_intake:
        raise VerifierError("circuitbreaker_intake_mismatch")
    config = intake["data"]["effective"]["config"]
    validate_contract(config)
    result = dict(data["yukon_projection"])
    if not data["judgment"]["valid"]:
        return result
    manifest, files = intake["data"]["manifest"], intake["data"]["files"]
    try:
        bound = accepted_bound(manifest, files)
        if bound is None:
            if data["score_records"] or data["judgment"]["scores"] is not None:
                raise ExpressionError("lemma must not carry main scores")
            return result
        scores = validate_bound(bound)
    except ExpressionError as exc:
        raise VerifierError("circuitbreaker_bound_contract: " + str(exc)) from exc
    # Compare source identity, full exact values and projections against the
    # rebuilt gate. Acceptance of an unrelated pair of claims is insufficient.
    expected_records = [{
        "objective": config["objectives"][i]["name"],
        "source_sha256": digest({"source": {"claim": i}, "files": digest(files)}),
        "record": score,
        "projection": ("unsupported_positive_infinity" if score == "positive_infinity" else "finite"),
    } for i, score in enumerate(scores)]
    numeric = ([score_value(s) for s in scores]
               if bound["projection"]["kind"] == "finite" else None)
    if (canonical(data["score_records"]) != canonical(expected_records)
            or canonical(data["judgment"]["scores"]) != canonical(numeric)):
        raise VerifierError("circuitbreaker_projection_disagreement")
    result["metrics"] = {"circuitbreaker": {
        "version": VERSION,
        "convention": CONVENTION,
        "status": "accepted_provisional",
        "bound": bound,
        "intake_sha256": intake["sha256"],
        "dossier_sha256": record["sha256"],
        "selected_review_sha256": data["selected_review_sha256"],
        "submission_sha256": intake["bindings"]["submission"],
        "verifier_sha256": intake["bindings"]["verifier"],
        "statements": [c["statement"] for c in manifest["claims"][:2]],
        "score_records": expected_records,
        "evidence": {p: files[p]["sha256"] for p in manifest["claims"][0]["evidence"]},
    }}
    result["metrics"]["circuitbreaker_frontier"] = {
        "version": "circuitbreaker-frontier-v1",
        "comparator": "circuitbreaker-v1",
        "scope": CONVENTION,
        "model": "full-b2-v1",
        "notation": "lower-bound",
        "intake_sha256": intake["sha256"],
        "dossier_sha256": record["sha256"],
        "resources": {"circuit_size": {
            "statement": manifest["claims"][0]["statement"],
            "bound": bound["guarantee"],
            "source_sha256": expected_records[0]["source_sha256"],
            "evidence": result["metrics"]["circuitbreaker"]["evidence"],
        }},
    }
    return result


def export(dossier_path, output, *, intake_path=None):
    inputs = [dossier_path] + ([intake_path] if intake_path else [])
    distinct_paths(inputs, [output])
    clear_output(output)
    expected = None
    if intake_path:
        intake = read_json(intake_path)
        verify_envelope(intake, "intake")
        expected = intake["sha256"]
    write_json(output, projection(read_json(dossier_path), expected_intake=expected))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dossier", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--intake", help="separate trusted same-run intake artifact")
    args = parser.parse_args()
    export(args.dossier, args.out, intake_path=args.intake)


if __name__ == "__main__":
    main()
