"""Execute the actual metrics bootstrap without network or action credentials."""

import base64
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest

from export_hosted import ARTIFACT_NAMES
from tools.generate_hosted import BASELINE_RECOVERY

ROOT = Path(__file__).resolve().parents[1]


class ArtifactOverlayTests(unittest.TestCase):
    def bootstrap(self, workspace, encoded, expected):
        overlay = (ROOT / "tools/metrics-job.yml").read_text()
        code = textwrap.dedent(overlay.split("          python3 -I - <<'PY'\n", 1)[1].split("          PY\n", 1)[0])
        return subprocess.run([sys.executable, "-I", "-c", code], capture_output=True,
                              env={"PATH": os.environ["PATH"], "GITHUB_WORKSPACE": str(workspace),
                                   "LLMV_TRANSFER_HELPER": encoded, "LLMV_TRANSFER_SHA256": expected})

    def test_bootstrap_requires_exact_bounded_helper_bytes_and_fresh_directory(self):
        payload = b'raise RuntimeError("Bootstrap must not execute downloaded bytes")\n'
        encoded = base64.b64encode(payload).decode()
        digest = hashlib.sha256(payload).hexdigest()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.bootstrap(root, encoded, digest)
            self.assertEqual(result.returncode, 0, result.stderr)
            target = root / ".llmv-artifact-transfer/transfer.py"
            self.assertEqual(target.read_bytes(), payload)
            self.assertNotEqual(self.bootstrap(root, encoded, digest).returncode, 0)
            self.assertEqual(target.read_bytes(), payload)
        for text, sha in ((encoded, "0" * 64), ("not-base64!", digest), ("", digest),
                          ("A" * 87385, digest)):
            with self.subTest(encoded_size=len(text)), tempfile.TemporaryDirectory() as tmp:
                self.assertNotEqual(self.bootstrap(tmp, text, sha).returncode, 0)
                self.assertFalse((Path(tmp) / ".llmv-artifact-transfer").exists())

    def test_declared_recovery_matches_metrics_only_inputs_and_outputs(self):
        self.assertEqual(BASELINE_RECOVERY["artifactNames"], list(ARTIFACT_NAMES))
        overlay = (ROOT / "tools/metrics-job.yml").read_text()
        self.assertIn("name: " + BASELINE_RECOVERY["deliveryJob"], overlay)
        self.assertIn("name: " + BASELINE_RECOVERY["lineageArtifactName"], overlay)
        self.assertIn("/" + BASELINE_RECOVERY["lineagePath"], overlay)
        self.assertNotIn("actions/download-artifact@", overlay)
        self.assertNotIn("actions: read", overlay)
        self.assertNotIn("ACTIONS_RUNTIME_TOKEN", overlay)
        self.assertNotIn("OPENAI_API_KEY", overlay)
        for name in ("package", "intake", "dossier"):
            self.assertIn("needs.verify.outputs." + name + "_artifact", overlay)


if __name__ == "__main__":
    unittest.main()
