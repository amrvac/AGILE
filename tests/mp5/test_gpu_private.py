#!/usr/bin/env python3
"""Check conditional FV buffers in generated GPU directives; no GPU required."""
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    fypp = shutil.which('fypp') or str(ROOT / '.venv/bin/fypp')
    checked = 0
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / 'fv.f90'
        for physics in ('hd', 'mhd'):
            for backend in ('OPENACC', 'OPENMP'):
                for source in (None, 'SOURCE_LOCAL', 'SOURCE_NONLOCAL', 'SOURCE_COMPACT'):
                    flags = [f'-DPHYS={physics!r}', "-DGEOM='Cartesian'", f'-D{backend}']
                    if source:
                        flags.append(f'-D{source}')
                    subprocess.run([fypp, *flags, '-I', str(ROOT / 'src'),
                                    str(ROOT / 'src/mod_finite_volume.fpp'), str(output)],
                                   check=True)
                    kernels = re.findall(
                        r'  subroutine (finite_volume_local_\w+)\(.*?'
                        r'  end subroutine \1', output.read_text(), re.S)
                    assert len(kernels) == 18, (physics, backend, source, kernels)
                    text = output.read_text()
                    for name in kernels:
                        body = text.split('  subroutine ' + name + '(', 1)[1].split(
                            '  end subroutine ' + name, 1)[0]
                        declared = bool(re.search(r'::\s*uprim\(', body))
                        private = any('uprim' in re.findall(r'\w+', clause)
                                      for clause in re.findall(r'private\(([^)]*)\)', body))
                        expected = '_muscl_muscl_' in name or source is not None
                        assert declared == private == expected, (physics, backend, source, name)
                        checked += 1
    print(f'PASS: {checked} GPU kernels have consistent uprim declarations/private clauses')


if __name__ == '__main__':
    main()
