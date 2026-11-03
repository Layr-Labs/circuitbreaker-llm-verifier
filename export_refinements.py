"""Prepare a finite owner-page overlay from authentic original v1 artifacts.

Offline only: no workflow dispatch, delivery recovery, reorg, publication, or new
review. Organizer-pinned input hashes provide the trust root. Publishing the
prepared mapping in an owner-controlled page revision supplies mapping authority;
the output does not certify mathematical truth or human acceptance.
"""

import argparse
import copy
import os
from pathlib import Path
import subprocess
import tempfile

from expressions import VERSION, VERSION_V2, validate_bound
from refinement.common import (RefinementError, canonical, digest, fields, hash_string,
                               read_json, record_id, require, sha256, validate_origin)

ROOT = Path(__file__).resolve().parent
PLAN_VERSION = "circuitbreaker-refinement-plan-v1"
OVERLAY_VERSION = "circuitbreaker-owner-refinements-v1"
PLAN_FIELDS = ("id", "bundle", "downloadReceiptSha256", "platform", "platformSha256",
               "originalSource", "originalPython", "targetBound")
EXPORTER_FILES = ("export_refinements.py", "expressions.py", "refinement/__init__.py",
                  "refinement/common.py", "refinement/replay_original.py",
                  "refinement/original-source-v1.json")


def validate_target(original, target):
    """Check structural transcription limits; never infer a rate from liminf."""
    validate_bound(original)
    validate_bound(target)
    require(original["version"] == VERSION and original["guarantee"]["kind"] == "polynomial_log",
            "original v1 polynomial declaration required")
    require(target["version"] == VERSION_V2 and target["guarantee"]["kind"] == "eventual_lower_bound",
            "target must be a full v2 eventual declaration")
    require(original["projection"]["kind"] == target["projection"]["kind"] == "finite",
            "finite projection required")
    for name in ("convention", "family", "theorem", "projection"):
        require(canonical(original[name]) == canonical(target[name]), "target changes original " + name)
    require(canonical(original["guarantee"]) == canonical(target["guarantee"]["leading"]),
            "target changes original leading guarantee")
    canonical(target)  # No numbers or lone surrogates in the shared format.


def payload(item, report):
    record_id(item["id"])
    validate_origin(report["origin"])
    validate_target(report["originalBound"], item["targetBound"])
    return {"id": item["id"], "origin": copy.deepcopy(report["origin"]),
            "targetBound": copy.deepcopy(item["targetBound"])}


def validate_plan(plan):
    fields(plan, ("version", "records"))
    require(plan["version"] == PLAN_VERSION and type(plan["records"]) is list and
            1 <= len(plan["records"]) <= 32, "bounded refinement plan required")
    ids = set()
    for item in plan["records"]:
        fields(item, PLAN_FIELDS)
        record_id(item["id"])
        require(item["id"] not in ids, "duplicate refinement id")
        ids.add(item["id"])
        for name in ("downloadReceiptSha256", "platformSha256"):
            hash_string(item[name])
        for name in ("bundle", "platform", "originalSource", "originalPython"):
            require(type(item[name]) is str and bool(item[name]), "input path required")
        validate_bound(item["targetBound"])
    canonical(plan)


def run_replay(item, scratch):
    request = scratch / "request.json"
    request.write_bytes(canonical(item) + b"\n")
    output = scratch / "replay"
    command = [item["originalPython"], "-I", "-B", "-X",
               "pycache_prefix=" + str(scratch / "empty-bytecode-cache"),
               str(ROOT / "refinement/replay_original.py"),
               "--request", str(request), "--out", str(output)]
    # Do not forward provider credentials, PYTHONPATH, shell configuration or
    # application environment into the original verifier process.
    completed = subprocess.run(command, env={"PATH": os.defpath}, cwd=scratch,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, timeout=180)
    require(completed.returncode == 0, "original replay failed: " + completed.stderr[-4000:])
    return read_json(output / "replay.json")


def exporter_identity():
    return {name: sha256((ROOT / name).read_bytes()) for name in EXPORTER_FILES}


def export(plan_path, output):
    plan_path, output = Path(plan_path).resolve(), Path(output).absolute()
    require(not output.exists(), "output directory must not exist")
    require(output.parent.is_dir(), "output parent directory must exist")
    raw_plan = plan_path.read_bytes()
    plan = read_json(plan_path)
    validate_plan(plan)
    source_files = exporter_identity()
    reports, prepared, origins = [], [], set()
    with tempfile.TemporaryDirectory(prefix="circuitbreaker-refinement-") as temporary:
        for index, entry in enumerate(plan["records"]):
            item = copy.deepcopy(entry)
            for name in ("bundle", "platform", "originalSource", "originalPython"):
                path = plan_path.parent / item[name]
                # Resolving the venv python symlink selects the system Python
                # and silently loses the exact original installed wheel.
                item[name] = str(path.absolute() if name == "originalPython" else path.resolve())
            # All original artifacts are read-only; a mapping must not be
            # emitted inside either original source or artifact directory.
            for name in ("bundle", "originalSource"):
                require(not output.resolve().is_relative_to(Path(item[name])), "output overlaps original inputs")
            scratch = Path(temporary) / str(index)
            scratch.mkdir()
            report = run_replay(item, scratch)
            current = payload(item, report)
            identity = (current["origin"]["benchmarkId"], current["origin"]["submissionId"])
            require(identity not in origins, "duplicate origin submission")
            origins.add(identity)
            reports.append({"id": item["id"], "targetBound": current["targetBound"], "replay": report})
            prepared.append(current)
    require(plan_path.read_bytes() == raw_plan, "plan changed during export")
    require(exporter_identity() == source_files, "exporter source changed during export")
    manifest_sha = sha256(raw_plan)
    receipt = {"version": "circuitbreaker-refinement-receipt-v1", "manifestSha256": manifest_sha,
               "exporterFiles": source_files, "exporterSourceSha256": digest(source_files),
               "records": reports,
               "scope": "Prepared owner-page mappings; original review and human status are unchanged. "
                        "The organizer must check the explicit target against the authenticated full theorem. "
                        "No inference from liminf, new scientific acceptance, publication or live mutation."}
    receipt_bytes = canonical(receipt) + b"\n"
    audit = {"manifestSha256": manifest_sha, "receiptSha256": sha256(receipt_bytes)}
    records = []
    for current in prepared:
        current["audit"] = dict(audit)
        records.append({**current, "recordSha256": digest(current)})
    overlay = {"version": OVERLAY_VERSION, "records": records}
    # Complete verification before creating any final artifact. x mode and a
    # fresh directory protect original inputs and previous prepared mappings.
    output.mkdir()
    for name, raw in (("manifest.json", raw_plan), ("receipt.json", receipt_bytes),
                      ("representation-refinements.json", canonical(overlay) + b"\n")):
        with (output / name).open("xb") as handle:
            handle.write(raw)
    return overlay


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="finite organizer-prepared plan JSON")
    parser.add_argument("--out", required=True, help="new output directory")
    args = parser.parse_args()
    try:
        overlay = export(args.manifest, args.out)
    except (RefinementError, ValueError, OSError, KeyError, subprocess.TimeoutExpired) as exc:
        parser.exit(1, f"refinement export rejected: {exc}\n")
    print(f"Prepared {len(overlay['records'])} owner-page mapping(s) in {args.out}; nothing published.")


if __name__ == "__main__":
    main()
