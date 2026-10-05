"""Small, unconstrained CLI wiring tests; not a ground-state physics benchmark."""
import pathlib
import subprocess
import tempfile
import struct

executable = pathlib.Path(__file__).resolve().parents[1] / "soft_core_stealthy.out"
with tempfile.TemporaryDirectory(prefix="cco-cli-") as directory:
    root = pathlib.Path(directory)
    for mode in ("txt", "binary", "off"):
        prefix = root / mode
        # Ten particles, two components; no constrained nonzero wavevectors.
        text = f"""random
2
1
10 1
0.01
0
2
5 0 wavenumber 0 0.1
5 0 wavenumber 0 0.1
1
{prefix}
1
ground
combined
maxsteps 1
export {mode}
run
"""
        result = subprocess.run([str(executable), "multi", "1", "17"], input=text,
                                text=True, capture_output=True, cwd=root, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        txt = pathlib.Path(str(prefix) + "_Success_0.txt")
        binary = pathlib.Path(str(prefix) + "_Success_0.mspb")
        assert txt.exists() == (mode == "txt"), result.stdout
        assert binary.exists() == (mode == "binary"), result.stdout
        if txt.exists():
            assert "# coordinates relative" in txt.read_text()
            rows = [line.split() for line in txt.read_text().splitlines() if not line.startswith("#")]
            assert len(rows) == 10
            assert [int(row[1]) for row in rows] == [0]*5 + [1]*5
        if binary.exists():
            data = binary.read_bytes()
            assert struct.unpack("<8sIIQQI", data[:36]) == (b"MSPB\r\n\x1a\n", 1, 2, 10, 2, 3)
            assert len(data) == 36 + 8*4 + 10*36
print("CLI export tests passed: txt, binary, off.")
