"""Shared hashing contract for owner page mappings (not an approval mechanism)."""

import hashlib
import json
from pathlib import Path
import re


class RefinementError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise RefinementError(message)


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def checked_json(value, *, numbers=False, depth=0, count=None):
    count = [0] if count is None else count
    count[0] += 1
    require(depth <= 40 and count[0] <= 100_000, "JSON complexity limit")
    if type(value) is str:
        try:
            value.encode("utf-8")
        except UnicodeError as exc:
            raise RefinementError("lone surrogate in JSON string") from exc
    elif type(value) is dict:
        for key, item in value.items():
            require(type(key) is str, "JSON object key must be a string")
            checked_json(key, numbers=numbers, depth=depth + 1, count=count)
            checked_json(item, numbers=numbers, depth=depth + 1, count=count)
    elif type(value) is list:
        for item in value:
            checked_json(item, numbers=numbers, depth=depth + 1, count=count)
    elif value is None or type(value) is bool:
        pass
    elif numbers and type(value) in (int, float):
        # json.dumps below rejects nonfinite floats; arbitrary numbers are never
        # admitted to the shared mapping/namespace digest format.
        pass
    else:
        raise RefinementError("JSON numbers/non-JSON values forbidden in mapping hashes")


def canonical(value):
    checked_json(value)
    return json.dumps(value, sort_keys=True, ensure_ascii=True,
                      separators=(",", ":"), allow_nan=False).encode("ascii")


def digest(value):
    return sha256(canonical(value))


def read_json(path, *, numbers=False):
    raw = Path(path).read_bytes()
    require(len(raw) <= 16 * 1024 * 1024, "JSON byte limit")

    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result

    def constant(_):
        raise RefinementError("nonfinite JSON")

    try:
        value = json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)
        checked_json(value, numbers=numbers)
        json.dumps(value, allow_nan=False)
        return value
    except (ValueError, UnicodeError, RecursionError) as exc:
        raise RefinementError(f"invalid JSON in {path}: {exc}") from exc


def fields(value, names):
    require(type(value) is dict and set(value) == set(names), "unexpected object fields")


def identity(value):
    require(type(value) is str and bool(value.strip()) and len(value) <= 512,
            "nonempty identity string required (at most 512 code points)")
    checked_json(value)


def hash_string(value):
    require(type(value) is str and re.fullmatch(r"[a-f0-9]{64}", value), "SHA256 required")


def record_id(value):
    require(type(value) is str and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,159}", value),
            "bounded ASCII refinement id required")


EVIDENCE_FIELDS = ("resultId", "evaluatorId", "policyVersionId", "sourceRef",
                   "jobId", "commitSha", "evaluatedAt")
ORIGIN_IDENTITIES = ("benchmarkId", "submissionId", "packetId", "candidateKey",
                     "evaluationId", "validityVersionId")
ORIGIN_HASHES = ("originalEnvelopeSha256", "originalRichSha256", "intakeSha256",
                 "dossierSha256", "sourceSha256")


def validate_origin(origin):
    fields(origin, (*ORIGIN_IDENTITIES, *ORIGIN_HASHES, "scoreEvidence", "validityEvidence"))
    for key in ORIGIN_IDENTITIES:
        identity(origin[key])
    for key in ORIGIN_HASHES:
        hash_string(origin[key])
    for key in ("scoreEvidence", "validityEvidence"):
        fields(origin[key], EVIDENCE_FIELDS)
        for value in origin[key].values():
            identity(value)
