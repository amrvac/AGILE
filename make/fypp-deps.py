from collections.abc import Iterable
from io import TextIOBase
from pathlib import Path
from dataclasses import dataclass
from textwrap import indent

import re
import argparse
import sys


INCLUDE_EX = re.compile(r"\s*#:include\s+(?P<quote>['\"])(?P<filename>[^'\"]+)(?P=quote)")


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="fypp-deps",
        description="""Scans a FYPP file (Fortran Preprocessor Language) for
        #:include statements that fortdepend doesn't register. Outputs Makefile
        compatible dependency statements that can be directly included by Make.

        The output (if any) is written to a file with the same basename as the
        target file, but the suffix changed to `.d`.
        """,
        epilog="""
        Conditional Includes: Fypp processes #:include statements in a first pass,
        resulting in a larger "virtual" code, in which other preprocessing statements
        are evaluated. There is no such thing as conditional inclusion.

        Any #:include statement is unconditional, even if the evaluation of the
        resulting expanded code may be conditional, so an #:include always
        leads to a dependency even when the included content goes unused.
        """)
    parser.add_argument("input")
    parser.add_argument("target")
    parser.add_argument("-I", dest="include_dirs", action="append", default=[],
                        type=Path, metavar="DIR",
                        help="""Directory searched for included files, after the
                        directory of the including file (the same order fypp
                        uses for its own -I). May be given more than once.""")

    return parser.parse_args()


@dataclass
class CircularDependencyError(Exception):
    visited: list[Path]
    file: Path

    def __str__(self) -> str:
        cycle = "\n".join(f"includes `{f.relative_to(self.visited[0].parent, walk_up=True)}`" for f in self.visited[1:] + [self.file])
        return f"Circular dependency detected in `{self.file}`: starting from `{self.visited[0]}`:\n" + indent(cycle, "  - ")


def resolve_include(filename: str, including: Path, include_dirs: list[Path]) -> Path:
    """Find an included file as fypp does: next to the including file first,
    then in each -I directory in turn."""
    for directory in [including.parent, *include_dirs]:
        candidate = directory / filename
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(f"`{filename}` included from `{including}`, searched in "
                            + ", ".join(f"`{d}`" for d in [including.parent, *include_dirs]))


def scan_deps(input_path: Path, include_dirs: list[Path],
              visited: list[Path] | None = None) -> set[Path]:
    visited = visited or list()
    if input_path in visited:
        raise CircularDependencyError(visited, input_path)

    visited.append(input_path)
    deps = set()

    with input_path.open("r") as fid:
        for line in fid:
            if m := INCLUDE_EX.match(line):
                d = resolve_include(m["filename"], input_path, include_dirs)
                deps.add(d)
                deps.update(scan_deps(d, include_dirs, visited))

    visited.pop()

    return deps


if __name__ == "__main__":
    args = parse_arguments()
    input_path = Path(args.input)

    try:
        deps = scan_deps(input_path, args.include_dirs)
    except (CircularDependencyError, FileNotFoundError) as e:
        print(f"Error scanning `{input_path}`:\n   ", e)
        sys.exit(-1)

    if deps:
        target_path = Path(args.target)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        with target_path.with_suffix(".d").open("w") as fid:
            print(args.target + ":", " ".join(map(str, deps)), file=fid)
