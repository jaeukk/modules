#!/usr/bin/env python3
"""Generate repository-scoped Doxygen docs without modifying source files."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'generated'
AREAS = ['core', 'pairstat', 'potential', 'models', 'nonspherical', 'phonons',
         'lp', 'misc', 'visualization', 'apps', 'examples']
EXTENSIONS = {'.c', '.h', '.cpp', '.hpp', '.cc', '.hh', '.cxx', '.py'}


def is_test_source(path):
    """Exclude test directories and standalone test/smoke translation units."""
    source = Path(path)
    stem = source.stem.lower()
    return (any(part.lower() in {'test', 'tests'} for part in source.parts[:-1])
            or stem in {'test', 'tests', 'smoke'}
            or stem.startswith(('test_', 'smoke_'))
            or stem.endswith(('_test', '_tests')))


def main():
    env = os.environ.copy()
    executable = env.get('DOXYGEN_BIN') or shutil.which('doxygen')
    isolated = Path('/tmp/modules-doxygen-tool')
    if not executable and (isolated / 'usr/bin/doxygen').is_file():
        executable = str(isolated / 'usr/bin/doxygen')
        libs = str(isolated / 'usr/lib/x86_64-linux-gnu')
        env['LD_LIBRARY_PATH'] = libs + (':' + env['LD_LIBRARY_PATH']
                                       if env.get('LD_LIBRARY_PATH') else '')
    if not executable:
        raise SystemExit('Doxygen is missing. Install it or set DOXYGEN_BIN; see README.md.')
    if not shutil.which('dot'):
        raise SystemExit('Graphviz dot is required for diagrams.')
    raw = subprocess.check_output(
        ['git', 'ls-files', '-co', '--exclude-standard', '-z', '--', *AREAS], cwd=ROOT)
    candidates = sorted({p for p in raw.decode().split('\0')
                    if p and Path(p).suffix in EXTENSIONS
                    and (ROOT / p).is_file() and not (ROOT / p).is_symlink()})
    paths = [p for p in candidates if not is_test_source(p)]
    excluded_tests = [p for p in candidates if is_test_source(p)]
    if not paths:
        raise SystemExit('No repository source files found.')
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    # Doxygen does not remove pages for inputs that disappeared. Clear only its
    # generated site/index directories so excluded tests cannot leave old pages.
    for directory in ('html', 'xml'):
        target = OUTPUT / directory
        if target.is_symlink():
            raise SystemExit(f'Refusing to clean a symlink: {target}')
        if target.exists():
            shutil.rmtree(target)
    (OUTPUT / 'validation.json').unlink(missing_ok=True)
    provenance = {
        'generated_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'repository': str(ROOT),
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip(),
        'git_status': subprocess.check_output(['git', 'status', '--short'], cwd=ROOT).decode(),
        'doxygen': subprocess.check_output([executable, '--version'], env=env).decode().strip(),
        'graphviz': subprocess.run(['dot', '-V'], capture_output=True, text=True).stderr.strip(),
        'source_count': len(paths), 'source_sha256': hashes,
        'excluded_test_sources': excluded_tests,
    }
    def quote(value):
        return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'
    config = (HERE / 'Doxyfile').read_text()
    config += '\nPROJECT_NUMBER = ' + quote(provenance['git_head'][:12] + ' working tree')
    config += '\nINPUT = ' + ' '.join(quote(p) for p in ['docs/doxygen/mainpage.md', *paths]) + '\n'
    (OUTPUT / 'effective.Doxyfile').write_text(config)
    start = time.monotonic()
    with (OUTPUT / 'build.log').open('w') as log:
        result = subprocess.run([executable, '-'], input=config, text=True,
                                cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    provenance['elapsed_seconds'] = time.monotonic() - start
    provenance['exit_code'] = result.returncode
    provenance['changed_during_build'] = [p for p, h in hashes.items()
        if not (ROOT / p).is_file() or hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != h]
    (OUTPUT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    if result.returncode or provenance['changed_during_build']:
        raise SystemExit('Build failed or source changed during generation; inspect generated logs.')
    if not (OUTPUT / 'html/index.html').is_file():
        raise SystemExit('Doxygen did not produce an HTML index.')
    print(f'Generated documentation for {len(paths)} source files: {OUTPUT / "html/index.html"}')


if __name__ == '__main__':
    main()
