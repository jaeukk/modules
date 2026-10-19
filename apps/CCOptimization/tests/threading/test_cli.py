"""Exercise the multicomponent schedule keywords in compiled and direct-only builds."""
import os
import re
from pathlib import Path
import subprocess
import tempfile

binary = Path(__file__).resolve().parents[2] / "soft_core_stealthy.out"
with tempfile.TemporaryDirectory(prefix="cco-sequential-cli-") as temporary:
    root = Path(temporary)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="4", OMP_THREAD_LIMIT="4", OMP_DYNAMIC="false", OMP_WAIT_POLICY="PASSIVE")
    for name, commands in [("sequential", "finufft-sequential"),
                           ("reset", "finufft-sequential\nfinufft")]:
        prefix = root / name
        text = f"random\n2\n1\n15 1\n0.01\n0\n3\n5 0 wavenumber 0 0.1\n5 0 wavenumber 0 0.1\n5 0 wavenumber 0 0.1\n4\n{prefix}\n1\nground\n{commands}\nmaxsteps 5\nrun\n"
        result = subprocess.run([str(binary), "multi", "1", "17"], input=text, text=True,
                                cwd=root, env=env, capture_output=True, timeout=30)
        if "backend not available" in result.stderr:
            assert result.returncode != 0
        else:
            assert result.returncode == 0, result.stdout + result.stderr
            assert re.findall(r"Species schedule: (\w+)", result.stdout) == (
                ["sequential"] if name == "sequential" else ["sequential", "concurrent"])
            assert Path(str(prefix) + "_Success.ConfigPack").is_file(), result.stdout
print("PASS: finufft-sequential, reset to finufft, and missing-backend behavior.")
