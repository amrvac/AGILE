#!/usr/bin/env python3
"""
Doc-build-only workaround for a FORD parser bug (reproduces on FORD 6.2.5
and 7.0.13): a same-line trailing `!<`/`!>` doc comment on a declaration
with a parenthesized attribute (e.g. `intent(in)`) or on a derived-type
component crashes FORD's parser with "Preceding documentation lines can
not be inline", silently dropping the whole file from the generated docs
(this is what happened to mod_physics.f90).

FORD's own documented convention places such comments on their own line
after the statement instead, which is semantically identical and parses
cleanly. This script rewrites every trailing `!<`/`!>` doc comment onto
its own following line. It only touches the preprocessed .f90 copies
under doc/ford/src/ -- never the real .fpp source.
"""
import sys
from pathlib import Path

# The docmark/predocmark set in project.md, plus FORD's default
# docmark_alt ('*') and predocmark_alt ('|') -- project.md doesn't enable
# those, but FORD still special-cases them, so an ordinary trailing
# comment that happens to start with '*' or '|' (e.g. `!* foo` in the
# vendored octree-mg code) can trip the same parser bug.
DOCMARKS = ("<", ">", "*", "|")


def split_trailing_doc(line):
    in_single = False
    in_double = False
    for i, c in enumerate(line):
        if in_single:
            if c == "'":
                in_single = False
        elif in_double:
            if c == '"':
                in_double = False
        elif c == "'":
            in_single = True
        elif c == '"':
            in_double = True
        elif c == "!":
            code, comment = line[:i], line[i:]
            if code.strip() and len(comment) >= 2 and comment[1] in DOCMARKS:
                indent = line[: len(line) - len(line.lstrip())]
                return code.rstrip() + "\n" + indent + comment
            return line
    return line


def process(path):
    text = path.read_text()
    lines = [split_trailing_doc(l) for l in text.splitlines()]
    new_text = "\n".join(lines) + ("\n" if text.endswith("\n") else "")
    if new_text != text:
        path.write_text(new_text)


def main():
    src_dir = Path(sys.argv[1])
    for f in src_dir.glob("*.f90"):
        process(f)


if __name__ == "__main__":
    main()
