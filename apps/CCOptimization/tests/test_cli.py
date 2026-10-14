"""Small, unconstrained CLI wiring tests; not a ground-state physics benchmark."""
import pathlib
import subprocess
import tempfile
import struct
import math

executable = pathlib.Path(__file__).resolve().parents[1] / "soft_core_stealthy.out"
with tempfile.TemporaryDirectory(prefix="cco-cli-") as directory:
    root = pathlib.Path(directory)
    for mode in ("txt", "binary", "off", "legacy"):
        prefix = root / mode
        export_command = "" if mode == "legacy" else f"export {mode}"
        analysis_command = "" if mode == "legacy" else "pairstat 3 0.2 0.4 0.2 16 17"
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
{export_command}
{analysis_command}
run
"""
        result = subprocess.run([str(executable), "multi", "1", "17"], input=text,
                                text=True, capture_output=True, cwd=root, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        txt = pathlib.Path(str(prefix) + "_Success_0.txt")
        binary = pathlib.Path(str(prefix) + "_Success_0.mspb")
        assert txt.exists() == (mode == "txt"), result.stdout
        assert binary.exists() == (mode == "binary"), result.stdout
        report_path = pathlib.Path(str(prefix) + "_PairStat_0.txt")
        assert report_path.exists() == (mode != "legacy")
        if report_path.exists():
            rows = [line.split() for line in report_path.read_text().splitlines() if not line.startswith("#")]
            assert {row[0] for row in rows} == {"S", "NV", "VFV"}
            assert {int(row[1]) for row in rows} == {-1, 0, 1}
            assert all(math.isfinite(float(value)) for row in rows for value in row[1:])
        if txt.exists():
            assert "# coordinates relative" in txt.read_text()
            rows = [line.split() for line in txt.read_text().splitlines() if not line.startswith("#")]
            assert len(rows) == 10
            assert [int(row[1]) for row in rows] == [0]*5 + [1]*5
        if binary.exists():
            data = binary.read_bytes()
            assert struct.unpack("<8sIIQQI", data[:36]) == (b"MSPB\r\n\x1a\n", 1, 2, 10, 2, 3)
            assert len(data) == 36 + 8*4 + 10*36
            resumed = subprocess.run([str(executable), "multi", "1", "17"],
                input=text.replace(f"{prefix}\n1\nground", f"{prefix}\n2\nground"),
                text=True, capture_output=True, cwd=root, timeout=30)
            assert resumed.returncode == 0, resumed.stdout + resumed.stderr
            assert binary.read_bytes() == data
            assert pathlib.Path(str(prefix) + "_Success_1.mspb").exists()
            assert pathlib.Path(str(prefix) + "_PairStat_1.txt").exists()
    for commands in ("export invalid\nrun\n", "pairstat 3 0.2 0.4 0.2 1 17\nrun\n", ""):
        bad = subprocess.run([str(executable), "multi", "1", "17"],
            input=text.replace("run\n", commands), text=True, capture_output=True,
            cwd=root, timeout=30)
        assert bad.returncode != 0, bad.stdout + bad.stderr
    # Single-component GetCCO ground run (random packing re-binned for the contact search), twice: same bytes.
    packs = []
    for run in ("cco_a", "cco_b"):
        text = f"2\n0\n1.0\n0\n1\n0.505\n0.2\n1\n200\n1\nrandom\n{root / run}\nground\nrun\n"
        result = subprocess.run([str(executable), "1", "17"], input=text,
                                text=True, capture_output=True, cwd=root, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        packs.append(pathlib.Path(str(root / run) + ".ConfigPack").read_bytes())
    assert packs[0] and packs[0] == packs[1]
    # GetCCO `finufft`: runs and repeats bitwise in a USE_FINUFFT=1 build, fails with a message otherwise.
    packs = []
    for run in ("cco_f_a", "cco_f_b"):
        text = f"2\n0\n1.0\n0\n1\n0.505\n0.2\n1\n200\n1\nrandom\n{root / run}\nground\nfinufft\nrun\n"
        result = subprocess.run([str(executable), "1", "17"], input=text,
                                text=True, capture_output=True, cwd=root, timeout=60)
        if "requires compiling" in result.stdout:
            assert result.returncode != 0, result.stdout
            break
        assert result.returncode == 0 and "FINUFFT backend" in result.stdout, result.stdout + result.stderr
        packs.append(pathlib.Path(str(root / run) + ".ConfigPack").read_bytes())
    assert not packs or packs[0] == packs[1]
print("CLI tests passed: text/binary/off, legacy defaults, PairStat, resume, invalid/EOF input, GetCCO repeat, GetCCO finufft.")
