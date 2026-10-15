"""Public protected-checker known-answer tests; no model calls or score output."""
import argparse
import json
import os
from pathlib import Path
import re
import selectors
import shutil
import subprocess
import tempfile
import time
import uuid


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def capture(command, request, env):
    # Bounded pipe capture; never buffer arbitrary checker output on the host.
    with tempfile.TemporaryFile() as stdin:
        stdin.write(request); stdin.seek(0)
        process = subprocess.Popen(command, stdin=stdin, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, env=env)
        selector = selectors.DefaultSelector()
        output = {"stdout": bytearray(), "stderr": bytearray()}
        for label in output:
            selector.register(getattr(process, label), selectors.EVENT_READ, label)
        deadline = time.monotonic() + 30
        try:
            while selector.get_map():
                if time.monotonic() >= deadline:
                    raise RuntimeError("Public checker timeout")
                for key, _ in selector.select(0.1):
                    chunk = os.read(key.fd, 16384)
                    if not chunk:
                        selector.unregister(key.fileobj)
                    else:
                        output[key.data].extend(chunk)
                        if len(output[key.data]) > 262144:
                            raise RuntimeError("Public checker output limit")
            if process.wait(timeout=max(0.01, deadline-time.monotonic())) != 0:
                raise RuntimeError("Public checker execution failed")
            return bytes(output["stdout"])
        finally:
            if process.poll() is None:
                process.kill()
            process.wait(timeout=5)
            selector.close()
            process.stdout.close(); process.stderr.close()


def run(root):
    if os.environ.get("OPENROUTER_API_KEY"):
        raise RuntimeError("Run public execution without model credentials")
    config = json.loads((root / "verifier.json").read_text())
    cases = json.loads((root / "public-checks.json").read_text())
    checks = {c["id"]: c for c in config.get("checkers", [])}
    if checks and not cases:
        raise RuntimeError("Registered checkers need public known-answer cases")
    image = config["runtime"]["image"]
    if not re.fullmatch(r"(?:docker\.io/library/)?python:[A-Za-z0-9._-]+@sha256:[a-f0-9]{64}", image):
        raise RuntimeError("Pinned Python image required")
    for case in cases:
        check = checks[case["check"]]
        with tempfile.TemporaryDirectory(prefix="llmv-public-") as tmp:
            snapshot = Path(tmp) / "files"
            snapshot.mkdir()
            for path in check["files"]:
                parts = Path(path).parts
                if Path(path).is_absolute() or any(p in (".", "..", ".git") for p in parts):
                    raise RuntimeError("Unsafe public bundle path")
                current = root
                for part in parts:
                    current = current / part
                    if current.is_symlink():
                        raise RuntimeError("Symlink in public bundle")
                if not current.is_file() or current.stat().st_nlink != 1:
                    raise RuntimeError("Regular public bundle files required")
                target = snapshot / path
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(current, target)
            os.chmod(tmp, 0o755)
            name = "llmv-public-" + uuid.uuid4().hex
            env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "HOME": tmp}
            command = ["docker", "run", "--name", name, "--pull=never", "--platform=linux/amd64",
                       "--network=none", "--read-only", "--user=65534:65534", "--cap-drop=ALL",
                       "--security-opt=no-new-privileges", "--memory=256m", "--cpus=1", "--pids-limit=32",
                       "--tmpfs=/tmp:rw,noexec,nosuid,size=16m", "--log-driver=none",
                       "--mount", "type=bind,src=" + str(snapshot) + ",dst=/snapshot/files,readonly",
                       "-i", image, "python", "-I", "/snapshot/files/" + check["source"]]
            try:
                stdout = capture(command, canonical(case["inputs"]).encode(), env)
                def unique(pairs):
                    result = {}
                    for key, value in pairs:
                        if key in result:
                            raise ValueError("Duplicate JSON key")
                        result[key] = value
                    return result
                output = json.loads(stdout, object_pairs_hook=unique)
                if canonical(output) != canonical(case["expected"]):
                    raise RuntimeError("Public known-answer mismatch: " + case["id"])
            finally:
                subprocess.run(["docker", "rm", "-f", name], stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, timeout=10, env=env, check=False)
                probe = subprocess.run(["docker", "container", "ls", "--all", "--filter",
                                        "name=^/" + name + "$", "--format", "{{.Names}}"],
                                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                       timeout=10, env=env, check=False)
                if probe.returncode != 0 or probe.stdout != b"":
                    raise RuntimeError("Public checker container cleanup unverified")
    print("Public checker cases passed:", len(cases), "(not mathematical qualification)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    run(parser.parse_args().root.absolute())
