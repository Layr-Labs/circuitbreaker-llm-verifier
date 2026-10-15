"""Protected CircuitBreaker encoding/arithmetic; never a mathematical proof oracle."""
from __future__ import annotations

import hashlib
import json
import math
import re
from fractions import Fraction

MODEL = "li-yang-b2-v1"
MAX_DIGITS = 256
MAX_INPUT_BYTES = 350_000
CLAIM_KEYS = {
    "version", "kind", "family", "model", "length_parameter", "quantifier",
    "coefficient", "theorem", "proof", "sources", "assumptions",
}
FAMILY_KEYS = {"id", "definition", "np_verifier", "witness_bound", "np_equivalence"}


class ClaimError(ValueError):
    pass


def closed_object(value, keys, where):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise ClaimError(f"{where} must have exactly {sorted(keys)}")


def text(value, where, limit=120_000):
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ClaimError(f"{where} must be nonempty text of at most {limit} characters")
    if any(0xD800 <= ord(c) <= 0xDFFF for c in value):
        raise ClaimError(f"{where} must contain Unicode scalar values")
    return value


def rational(value):
    closed_object(value, {"numerator", "denominator"}, "coefficient")
    for name, digits in value.items():
        if not isinstance(digits, str) or not re.fullmatch(r"0|[1-9][0-9]*", digits):
            raise ClaimError(f"{name} must be an unsigned canonical decimal integer string")
        if len(digits) > MAX_DIGITS:
            raise ClaimError(f"{name} exceeds the {MAX_DIGITS}-digit transport resource limit; not infinity")
    p, q = int(value["numerator"]), int(value["denominator"])
    if q == 0 or math.gcd(p, q) != 1:
        raise ClaimError("coefficient must be reduced with positive denominator (zero is 0/1)")
    return Fraction(p, q)


def _unique_pairs(pairs):
    result = {}
    for key, val in pairs:
        if key in result:
            raise ClaimError(f"duplicate JSON key: {key}")
        result[key] = val
    return result


def parse_json(raw):
    if not isinstance(raw, str):
        raise ClaimError("claim_json must be a string containing the complete claim JSON")
    try:
        raw_bytes = raw.encode("utf-8")
    except UnicodeError as exc:
        raise ClaimError("claim JSON must be valid UTF-8") from exc
    if len(raw_bytes) > MAX_INPUT_BYTES:
        raise ClaimError("claim exceeds canonical input byte limit")
    def reject_constant(token):
        raise ClaimError(f"non-JSON number {token} is forbidden")
    try:
        return json.loads(raw, object_pairs_hook=_unique_pairs, parse_constant=reject_constant)
    except (json.JSONDecodeError, RecursionError) as exc:
        raise ClaimError(f"invalid claim JSON: {exc}") from exc


def validate_claim(claim):
    closed_object(claim, CLAIM_KEYS, "claim")
    if type(claim["version"]) is not int or claim["version"] != 1:
        raise ClaimError("version must be integer 1")
    if claim["kind"] not in ("finite", "superlinear", "lemma"):
        raise ClaimError("kind must be finite, superlinear, or lemma")
    closed_object(claim["family"], FAMILY_KEYS, "family")
    for key, value in claim["family"].items():
        text(value, "family." + key, 20_000)
    if claim["model"] != MODEL:
        raise ClaimError("main claim must name the protected full B2 model")
    if claim["length_parameter"] != "n=|x|":
        raise ClaimError("length_parameter must be n=|x|")
    if claim["quantifier"] != "all-sufficiently-large-integer-lengths":
        raise ClaimError("main claim must declare all sufficiently large integer lengths")
    for key in ("theorem", "proof"):
        text(claim[key], key)
    for key in ("sources", "assumptions"):
        if not isinstance(claim[key], list) or len(claim[key]) > 100:
            raise ClaimError(f"{key} must be a bounded list")
        for item in claim[key]:
            text(item, key + " item", 4000)
    if claim["kind"] == "finite":
        rational(claim["coefficient"])
    elif claim["coefficient"] is not None:
        raise ClaimError("only finite claims have a coefficient; never use a finite sentinel")
    if claim["kind"] != "lemma" and claim["assumptions"]:
        raise ClaimError("a scored main theorem must be unconditional; submit a conditional lemma unscored")
    return claim


def canonical_bytes(claim):
    validate_claim(claim)
    return json.dumps(claim, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def check(raw):
    """All input is untrusted; passed means encoding, NOT mathematical validity."""
    claim = parse_json(raw)
    encoded = canonical_bytes(claim)
    claim_id = "claim:" + hashlib.sha256(encoded).hexdigest()
    result = {
        "passed": True,
        "claim_id": claim_id,
        "proof_sha256": hashlib.sha256(claim["proof"].encode("utf-8")).hexdigest(),
        "kind": claim["kind"],
        "scope": "canonical encoding and exact rational only; theorem requires review",
    }
    if claim["kind"] == "finite":
        result["score_record"] = {"kind": "finite", "claim_id": claim_id, **claim["coefficient"]}
    elif claim["kind"] == "superlinear":
        result["score_record"] = {"kind": "terminal", "claim_id": claim_id, "terminal": "positive_infinity"}
    if claim["kind"] == "lemma":
        raise ClaimError("lemmas use a submission with checks:[] and receive no main-claim record")
    return result


def compare_finite(left, right):
    """Authoritative record ordering; never order by rounded binary64 display."""
    a, b = rational(left), rational(right)
    return (a > b) - (a < b)


def check_inputs(inputs):
    """Runtime stdin is this inputs object; proof is canonical typed chunks, not files."""
    closed_object(inputs, {"claim_json", "proof_chunks"}, "inputs")
    if len(json.dumps(inputs, ensure_ascii=True).encode("utf-8")) > 16 * 1024:
        raise ClaimError("canonical checker inputs exceed 16 KiB; never truncate")
    text(inputs["claim_json"], "claim_json", 4096)
    chunks = inputs["proof_chunks"]
    if not isinstance(chunks, list) or not 1 <= len(chunks) <= 8:
        raise ClaimError("proof_chunks must contain 1 to 8 text chunks")
    for chunk in chunks:
        text(chunk, "proof chunk", 4096)
    meta = parse_json(inputs["claim_json"])
    closed_object(meta, CLAIM_KEYS - {"proof"}, "canonical claim metadata")
    meta["proof"] = "".join(chunks)
    result = check(json.dumps(meta, ensure_ascii=False, allow_nan=False))
    return {"passed": True, "score_record": result["score_record"], "message": result["scope"]}


def failed(message):
    # A schema-valid placeholder under passed:false. It cannot unlock a score.
    return {"passed": False, "score_record": {"kind": "finite", "claim_id": "claim:" + "0" * 64,
            "numerator": "0", "denominator": "1"}, "message": str(message)[:1000]}


if __name__ == "__main__":
    import sys
    try:
        data = json.load(sys.stdin, object_pairs_hook=_unique_pairs)
        output = check_inputs(data)
    except (ClaimError, ValueError, TypeError, RecursionError, UnicodeError) as exc:
        output = failed(exc)
    print(json.dumps(output, ensure_ascii=True, separators=(",", ":"), allow_nan=False))
