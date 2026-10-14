"""End-to-end archive CLI and optimizer integration tests, with independent MSPB bytes."""
import math
import os
from pathlib import Path
import struct
import subprocess
import tempfile

app = Path(__file__).resolve().parents[1]
cli = app / 'mspb_pack.out'
optimizer = app / 'soft_core_stealthy.out'


def run(args, cwd, text=None, ok=True):
    result = subprocess.run(list(map(str, args)), cwd=cwd, input=text, text=True,
                            capture_output=True, timeout=60,
                            env={**os.environ, 'OMP_NUM_THREADS': '1'})
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result.stdout


def payload(dim):
    data = struct.pack('<8sIIQQI', b'MSPB\r\n\x1a\n', 1, dim, 3, 4, 3)
    for i in range(dim):
        data += struct.pack('<' + 'd'*dim, *[
            4.0 if i == j else 0.25 if i == 1 and j == 0 else 0.0 for j in range(dim)])
    for i in range(3):
        data += struct.pack('<QI' + 'd'*(dim+1), 100+i*7, 2 if i == 1 else 0,
                            *[(i+j+1)/8 for j in range(dim)], (i+1)*0.02)
    return data


def frames(path):
    data = path.read_bytes()
    assert data[:16] == struct.pack('<8sII', b'MSPP\r\n\x1a\n', 1, 0)
    pos, result = 16, []
    while pos < len(data):
        size, = struct.unpack_from('<Q', data, pos)
        pos += 8
        result.append(data[pos:pos+size])
        assert len(result[-1]) == size
        pos += size
    assert pos == len(data)
    return result


def rows(path):
    return [line.split() for line in path.read_text().splitlines() if not line.startswith('#')]


with tempfile.TemporaryDirectory(prefix='mspb-cli-') as tmp:
    root = Path(tmp)
    for d in (2, 3):
        (root / f'd{d}.mspb').write_bytes(payload(d))
    run([cli, 'pack', 'ensemble', 'd2.mspb'], root)
    run([cli, 'pack', 'ensemble', 'd3.mspb'], root)
    assert frames(root / 'ensemble.MSPBPack') == [payload(2), payload(3)]
    assert '# snapshots = 2' in run([cli, 'inspect', 'ensemble'], root)
    run([cli, 'unpack', 'ensemble', 0, 2, 'restored'], root)
    assert (root / 'restored_0000.mspb').read_bytes() == payload(2)
    assert (root / 'restored_0001.mspb').read_bytes() == payload(3)
    run([cli, 'pairstat', 'ensemble', 0, 2, 'stats', 3, .2, .4, .2, 16, 17], root)
    for i, d in enumerate((2, 3)):
        report = (root / f'stats_{i:04d}.txt').read_text()
        assert '# component 1 particles 0' in report
        assert '# component 3 particles 0' in report
        data = rows(root / f'stats_{i:04d}.txt')
        assert {r[0] for r in data} == {'S', 'NV', 'VFV'}
        assert {int(r[1]) for r in data} == {-1, 0, 2}
        for row in data:
            assert all(math.isfinite(float(v)) for v in row[1:])
            if row[1] == '-1':
                sphere_volume = math.pi if d == 2 else 4*math.pi/3
                expected = sphere_volume * sum((j*.02)**d for j in (1, 2, 3)) / (4**d)
                assert math.isclose(float(row[3]), expected, rel_tol=1e-12)
    before = (root / 'ensemble.MSPBPack').read_bytes()
    (root / 'bad.mspb').write_bytes(payload(2)[:-1])
    run([cli, 'pack', 'ensemble', 'bad.mspb'], root, ok=False)
    assert (root / 'ensemble.MSPBPack').read_bytes() == before
    for args in [('inspect', 'missing'), ('unpack', 'ensemble', -1, 1, 'bad'),
                 ('unpack', 'ensemble', 0, 3, 'bad'),
                 ('unpack', 'ensemble', 0, 2, 'restored'),
                 ('pairstat', 'ensemble', 0, 1, 'bad', 3, .2, .4, .2, 1, 17)]:
        run([cli, *args], root, ok=False)

    def ground(prefix, successes=1, mode='pack'):
        return f'''random
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
{successes}
ground
combined
maxsteps 1
export {mode}
pairstat 3 0.2 0.4 0.2 16 17
run
'''
    run([optimizer, 'multi', 1, 17], root, ground('optimized'))
    archive = root / 'optimized_Success.MSPBPack'
    first = frames(archive)
    assert len(first) == 1
    assert not (root / 'optimized_Success_0.mspb').exists()
    run([optimizer, 'multi', 1, 17], root, ground('optimized', 2))
    assert len(frames(archive)) == 2 and frames(archive)[0] == first[0]
    run([cli, 'pairstat', 'optimized_Success', 0, 2, 'offline', 3, .2, .4, .2, 16, 17], root)
    for i in range(2):
        online, offline = rows(root / f'optimized_PairStat_{i}.txt'), rows(root / f'offline_{i:04d}.txt')
        assert len(online) == len(offline)
        for a, b in zip(online, offline):
            assert a[0] == b[0]
            assert all(math.isclose(float(x), float(y), rel_tol=1e-11, abs_tol=1e-12)
                       for x, y in zip(a[1:], b[1:]))
    # Enabling archive output after legacy successes must fail before another optimization.
    run([optimizer, 'multi', 1, 17], root, ground('legacy', mode='binary'))
    run([optimizer, 'multi', 1, 17], root, ground('legacy', 2), ok=False)
    assert not (root / 'legacy_Success.MSPBPack').exists()
    # Import the missing snapshot, then resume with matching zero-based archive indices.
    run([cli, 'pack', 'legacy_Success', 'legacy_Success_0.mspb'], root)
    run([optimizer, 'multi', 1, 17], root, ground('legacy', 2))
    assert len(frames(root / 'legacy_Success.MSPBPack')) == 2
print('PASS: archive CLI, byte restoration, offline PairStat, optimizer export/resume, and mismatch rejection')
