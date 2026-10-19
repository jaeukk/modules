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
                           ("default", "finufft"),
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
                ["sequential"] if name == "sequential" else
                ["concurrent"] if name == "default" else ["sequential", "concurrent"])
            assert "L-BFGS history in use = 40" in result.stdout
            assert Path(str(prefix) + "_Success.ConfigPack").is_file(), result.stdout
    # A bounded, constrained run verifies that the default and explicit history 40
    # reach the same endpoint, while a request for the former history is honored.
    endpoints = {}
    for history in (None, 40, 100):
        prefix = root / f"history-{history}"
        command = "" if history is None else f"history {history}\n"
        text = (f"random\n2\n1\n30 1\n0.1\n100\n3\n"
                f"10 0 chi 0 0.3\n10 0 chi 0 0.3\n10 0 chi 0 0.3\n"
                f"4\n{prefix}\n1\nground\n{command}combined\n"
                "maxsteps 5\ntolerance 1e10\nexport binary\nrun\n")
        result = subprocess.run([str(binary), "multi", "1", "17"], input=text,
                                text=True, cwd=root, env=env, capture_output=True, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        expected = 40 if history is None else history
        assert f"L-BFGS history in use = {expected}" in result.stdout, result.stdout
        endpoints[history] = Path(str(prefix) + "_Success_0.mspb").read_bytes()
    assert endpoints[None] == endpoints[40]
    for invalid in ("0", "-1", "1.5", "40junk", "4294967296", "9999999999999999999999999999999999"):
        result = subprocess.run([str(binary), "multi", "1", "17"],
                                input=text.replace("history 100", f"history {invalid}"),
                                text=True, cwd=root, env=env, capture_output=True, timeout=30)
        assert result.returncode != 0 and "Expected: history N" in result.stdout, result.stdout + result.stderr
    # Exercise MINOP's zero-iteration and candidate-energy paths. The loose export
    # threshold above captures bounded endpoints; it is not physical acceptance.
    for steps in (0, 5):
        commands = text.replace(str(prefix), str(root / f"minop-{steps}"))
        commands = commands.replace("maxsteps 5", f"algorithm MINOP\nmaxsteps {steps}")
        result = subprocess.run([str(binary), "multi", "1", "17"], input=commands,
                                text=True, cwd=root, env=env, capture_output=True, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        initial = float(re.search(r"E_init=([^\s]+)", result.stdout).group(1))
        final = float(re.search(r"E_relax=([^\s]+)", result.stdout).group(1))
        assert final <= initial, result.stdout
        assert (final == initial) if steps == 0 else (final < initial), result.stdout
print("PASS: schedule selection, history defaults/overrides/rejections, MINOP candidate paths, and missing-backend behavior.")
