"""CircuitBreaker representation checks for retained-validity score derivation.

This module validates the domain boundary of a *separately reviewed* score
derivation. It does not authenticate inputs or establish the mathematical truth
of a proposed transcription. The generic verifier owns those input bindings and
the derivation judgment. Original declarations and claim text remain immutable.
"""

from copy import deepcopy

from llm_verifier.errors import VerifierError
from llm_verifier.jsonutil import canonical, digest
from llm_verifier.scoring import score_value

from expressions import VERSION, VERSION_V2, ExpressionError, validate_bound


DERIVATION_VERSION = "circuitbreaker-score-derivation-v1"
SCORING_POLICY = {
    "version": "coefficient-exponent-expression-v3",
    "retained_result": "yukon-retained-result-v1",
}


def validate_transcription(original, target):
    """Accept the same declaration or a remainder transcription at the same lead.

    Only the guarantee's representation and its version may change. In
    particular, an accepted leading projection is not authority to alter the
    family, theorem, leading constant, or any quantifier in the source text.
    Semantic support for the remainder must be established by the score
    derivation reviewer from the authenticated original full claims.
    """
    try:
        original_scores = validate_bound(original)
        target_scores = validate_bound(target)
    except (ExpressionError, KeyError, TypeError) as exc:
        raise VerifierError("circuitbreaker_derivation_bound") from exc
    if canonical(original) == canonical(target):
        return deepcopy(target)
    if (
        original["version"] not in (VERSION, VERSION_V2)
        or original["guarantee"]["kind"] != "polynomial_log"
        or target["version"] != VERSION_V2
        or target["guarantee"]["kind"] != "eventual_lower_bound"
        or canonical(original_scores) != canonical(target_scores)
    ):
        raise VerifierError("circuitbreaker_derivation_representation")
    for name in ("convention", "family", "theorem", "projection"):
        if canonical(original[name]) != canonical(target[name]):
            raise VerifierError("circuitbreaker_derivation_changed_" + name)
    if canonical(original["guarantee"]) != canonical(target["guarantee"]["leading"]):
        raise VerifierError("circuitbreaker_derivation_changed_leading")
    return deepcopy(target)



def _project_verified(result, material):
    """Project already-verified derivation data; never call a provider here."""
    from export_metrics import FRONTIER_VERSIONS, accepted_bound, validate_contract
    from expressions import CONVENTION

    data = result["data"]
    intake = data["intake"]
    decision = data["decision"]
    evaluation = intake["data"]["evaluation"]
    retained = evaluation["retained_result"]
    if evaluation["policies"]["scoring"] != SCORING_POLICY:
        raise VerifierError("circuitbreaker_derivation_scoring_policy")
    metadata = {
        "version": DERIVATION_VERSION,
        "status": decision["status"],
        "intake_sha256": intake["sha256"],
        "dossier_sha256": result["sha256"],
        "retained_result": deepcopy(retained),
        "reason": decision["reason"],
        "source_refs": deepcopy(decision["source_refs"]),
    }
    output = {
        "scores": None,
        "derived_from": deepcopy(retained),
        "metrics": {"circuitbreaker_derivation": metadata},
    }
    if not retained["validity"]["valid"]:
        if decision["status"] != "not_applicable" or decision["payload"] is not None:
            raise VerifierError("circuitbreaker_derivation_rejected_source")
        return output
    original_intake = material["intake"]
    accepted = material["accepted_results"]
    manifest = material["manifest"]
    files = original_intake["data"]["files"]
    legacy = validate_contract(original_intake["data"]["effective"]["config"])
    try:
        original = accepted_bound(manifest, files)
    except (ExpressionError, KeyError, TypeError) as exc:
        raise VerifierError("circuitbreaker_derivation_original_claims") from exc
    metadata["original"] = {
        "intake_sha256": original_intake["sha256"],
        "dossier_sha256": accepted["source"]["dossier_sha256"],
        "selected_review_sha256": accepted["source"]["review_sha256"],
        "accepted_results_sha256": accepted["sha256"],
    }
    if original is None:
        if decision["status"] != "not_applicable" or decision["payload"] is not None:
            raise VerifierError("circuitbreaker_derivation_lemma")
        return output
    if legacy and original["version"] != VERSION:
        raise VerifierError("circuitbreaker_derivation_legacy_v2")
    if decision["status"] == "derived":
        payload = decision["payload"]
        if type(payload) is not dict or set(payload) != {"target_bound"}:
            raise VerifierError("circuitbreaker_derivation_payload")
        target = validate_transcription(original, payload["target_bound"])
        if not any(ref["source"] in ("claim:0", "claim:1") for ref in decision["source_refs"]):
            raise VerifierError("circuitbreaker_derivation_joint_claim_citation")
    elif decision["status"] == "unsupported" and decision["payload"] is None:
        # This is the original accepted declaration, not a guessed remainder.
        target = deepcopy(original)
    else:
        raise VerifierError("circuitbreaker_derivation_main_claim")
    scores = validate_bound(original)
    config = original_intake["data"]["effective"]["config"]
    original_records = [{
        "objective": config["objectives"][index]["name"],
        "source_sha256": digest({"source": {"claim": index}, "files": digest(files)}),
        "record": score,
        "projection": "unsupported_positive_infinity" if score == "positive_infinity" else "finite",
    } for index, score in enumerate(scores)]
    evidence = {path: files[path]["sha256"] for path in manifest["claims"][0]["evidence"]}
    # Keep the original rich metric's meaning. The new representation and
    # its scoring review are explicitly recorded in separate fields.
    output["metrics"]["circuitbreaker"] = {
        "version": original["version"],
        "convention": CONVENTION,
        "status": "accepted_provisional",
        "bound": deepcopy(original),
        "intake_sha256": original_intake["sha256"],
        "dossier_sha256": accepted["source"]["dossier_sha256"],
        "selected_review_sha256": accepted["source"]["review_sha256"],
        "submission_sha256": original_intake["bindings"]["submission"],
        "verifier_sha256": original_intake["bindings"]["verifier"],
        "statements": [claim["statement"] for claim in manifest["claims"][:2]],
        "score_records": original_records,
        "evidence": evidence,
    }
    version, comparator = FRONTIER_VERSIONS[target["version"]]
    metadata["target_bound"] = target
    output["scores"] = ([score_value(score) for score in scores]
                        if target["projection"]["kind"] == "finite" else None)
    output["metrics"]["circuitbreaker_frontier"] = {
        "version": version,
        "comparator": comparator,
        "scope": CONVENTION,
        "model": "full-b2-v1",
        "notation": "lower-bound",
        "intake_sha256": intake["sha256"],
        "dossier_sha256": result["sha256"],
        "resources": {"circuit_size": {
            "statement": manifest["claims"][0]["statement"],
            "bound": deepcopy(target["guarantee"]),
            "source_sha256": digest({
                "original_source_sha256": original_records[0]["source_sha256"],
                "derivation_sha256": result["sha256"],
            }),
            "evidence": evidence,
        }},
    }
    return output


def projection(result, *, expected_intake=None, preview=False):
    """Export a verified scoring result or an explicitly nonpublishing preview."""
    from llm_verifier.derivation import source_material, verify_derivation

    data = verify_derivation(result, expected_intake=expected_intake)
    if result["simulation"] and not preview:
        raise VerifierError("simulation_cannot_publish_score")
    material = source_material(data["intake"])
    output = _project_verified(result, material)
    if preview:
        # Intentionally no top-level scores: this object is not a Yukon result.
        return {"simulation": result["simulation"], "publishable": False,
                "score_preview": output}
    return output
