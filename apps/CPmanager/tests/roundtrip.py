#!/usr/bin/env python3
"""Exercise the actual CLI using independent text fixtures and numeric checks."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile

exe = Path(sys.argv[1]).resolve()


def run(*args, ok=True):
    result = subprocess.run([str(exe), *map(str, args)], capture_output=True, text=True, timeout=20)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result.stdout


def compare(expected, actual):
    a, b = expected.split(), actual.split()
    assert len(a) == len(b), (a, b)
    for x, y in zip(a, b):
        assert math.isclose(float(x), float(y), rel_tol=2e-14, abs_tol=2e-14), (x, y)


with tempfile.TemporaryDirectory(prefix='spherepack-') as tmp:
    root = Path(tmp)
    # Skew boxes, different dimensions/counts, and unequal radii with >10 digits.
    fixtures = [
        '2\n2 0 0 0\n0.5 1.5 0 0\n0.3212345678901234 0.45 0 0 0.1234567890123456\n1.2 0.9 0 0 0.2345678901234567\n',
        '3\n2 0 0 0\n0.25 1.5 0 0\n0.125 0.25 2 0\n0.4 0.5 0.6 0 0.3456789012345678\n',
        '2\n1 0 0 0\n0 1 0 0\n',  # Empty packing with a valid box.
    ]
    for i, fixture in enumerate(fixtures):
        (root / f'sphere_{i}.txt').write_text(fixture)
    pattern, pack = root / 'sphere_%d', root / 'packed'
    run(5, pattern, 0, 2, pack)
    run(5, pattern, 2, 1, pack)  # Reopen and append.
    assert '# configs = 3' in run(3, pack)
    assert pack.with_suffix('.SpherePack').stat().st_size > 0
    assert pack.with_suffix('.SpherePackIndex').stat().st_size > 0
    assert not pack.with_suffix('.ConfigPack').exists()
    run(4, pack, 1, 99)  # Clamp range, preserving original indices.
    assert not (root / 'packed__0000.txt').exists()
    run(4, pack, 0, 1)
    for i, fixture in enumerate(fixtures):
        compare(fixture, (root / f'packed__{i:04d}.txt').read_text())
    restored = root / 'restored'
    run(5, root / 'packed__%04d', 0, 3, restored)
    run(4, restored, 0, 3)
    for i, fixture in enumerate(fixtures):
        compare(fixture, (root / f'restored__{i:04d}.txt').read_text())
    # Existing ConfigPack modes remain usable with Cartesian point text.
    points = '2\n2 0 0 0\n0 2 0 0\n0.25 0.5 0 0\n1.25 1.5 0 0\n'
    (root / 'points_0.txt').write_text(points)
    run(2, root / 'points_%d', 0, 1, root / 'points')
    assert '# configs = 1' in run(0, root / 'points')
    run(1, root / 'points', 0, 1)
    compare(points, (root / 'points__0000.txt').read_text())
    for args in [(6,), ('bad',), (5,), (4, pack, -1, 1), (4, pack, 'bad', 1),
                 (4, pack, 2147483647, 1), (5, root / '%s', 0, 1, pack)]:
        run(*args, ok=False)
print('PASS: SpherePack CLI round trip, append, range clamping, empty packing, and ConfigPack regression')
