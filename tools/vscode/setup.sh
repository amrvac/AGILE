#!/usr/bin/env bash
#
# Set up VS Code for this repository. Idempotent -- re-run it after rebuilding
# .venv, or after building the lint case for the first time.
#
#   tools/vscode/setup.sh [--case <dir>] [--no-extension]
#
# It does three things:
#   1. makes sure .venv has fypp and fortls (uv sync --extra ide)
#   2. installs the fypp syntax grammar as a local VS Code extension
#   3. writes .vscode/settings.json from settings.template.json
#
# Diagnostics are linted against one build: its config.mk supplies the fypp
# defines and its f90/ directory the .mod files. By default that is the most
# recently built case anywhere in the tree; --case names one explicitly, e.g.
#   tools/vscode/setup.sh --case tests/mhd/spherical_pole
#
# See tools/vscode/README.md.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
want_extension=1
lint_case=""               # empty: discover the most recently built case

while [ $# -gt 0 ]; do
  case "$1" in
    --no-extension) want_extension=0; shift ;;
    --case)         lint_case="${2:?--case needs a directory}"; shift 2 ;;
    -h|--help)      sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)              echo "error: unknown argument '$1'" >&2; exit 2 ;;
  esac
done

say() { printf '  %-13s %s\n' "$1" "$2"; }

# 1. Python environment ------------------------------------------------------
if [ ! -x "$repo/.venv/bin/fypp" ] || [ ! -x "$repo/.venv/bin/fortls" ]; then
  uv="$(command -v uv || true)"
  [ -z "$uv" ] && [ -x "$HOME/.local/bin/uv" ] && uv="$HOME/.local/bin/uv"
  if [ -z "$uv" ]; then
    echo "error: .venv is missing fypp and/or fortls, and uv was not found." >&2
    echo "       See the Installation section of README.md, then run:" >&2
    echo "         cd $repo && uv sync --extra ide" >&2
    exit 1
  fi
  (cd "$repo" && "$uv" sync --extra ide)
fi
say "python" ".venv (fypp, fortls)"

# 2. fypp grammar ------------------------------------------------------------
if [ "$want_extension" -eq 0 ]; then
  say "extension" "skipped (--no-extension)"
elif ! command -v code >/dev/null 2>&1; then
  say "extension" "skipped, 'code' is not on PATH"
else
  vsix="$(mktemp -d)/fypp-fortran.vsix"
  python3 - "$here/fypp-fortran" "$vsix" <<'PY'
import json, sys, zipfile
from pathlib import Path

src, out = Path(sys.argv[1]), Path(sys.argv[2])
pkg = json.loads((src / "package.json").read_text())

(out.parent).mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("extension.vsixmanifest", f"""<?xml version="1.0" encoding="utf-8"?>
<PackageManifest Version="2.0.0" xmlns="http://schemas.microsoft.com/developer/vsx-schema/2011">
  <Metadata>
    <Identity Language="en-US" Id="{pkg['name']}" Version="{pkg['version']}" Publisher="{pkg['publisher']}" />
    <DisplayName>{pkg['displayName']}</DisplayName>
    <Description xml:space="preserve">{pkg['description']}</Description>
    <Categories>Programming Languages</Categories>
    <Properties>
      <Property Id="Microsoft.VisualStudio.Code.Engine" Value="{pkg['engines']['vscode']}" />
      <Property Id="Microsoft.VisualStudio.Code.ExtensionKind" Value="ui,workspace" />
    </Properties>
  </Metadata>
  <Installation><InstallationTarget Id="Microsoft.VisualStudio.Code" /></Installation>
  <Dependencies />
  <Assets>
    <Asset Type="Microsoft.VisualStudio.Code.Manifest" Path="extension/package.json" Addressable="true" />
  </Assets>
</PackageManifest>
""")
    z.writestr("[Content_Types].xml", """<?xml version="1.0" encoding="utf-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension=".json" ContentType="application/json" />
  <Default Extension=".xml" ContentType="text/xml" />
</Types>
""")
    for f in ["package.json", "syntaxes/fypp.injection.tmLanguage.json"]:
        z.write(src / f, f"extension/{f}")
PY
  # Quiet on success, but show everything if the install fails.
  if ! out="$(code --install-extension "$vsix" --force 2>&1)"; then
    printf '%s\n' "$out" >&2
    exit 1
  fi
  say "extension" "agile.fypp-fortran (reload the window to load it)"
fi

# 3. Workspace settings ------------------------------------------------------
# gfortran writes a .mod per module even under -fsyntax-only, and dies if the
# -J directory is missing, so create it here rather than on the first lint.
mkdir -p "$repo/build/lint-mods"

summary="$(python3 - "$here/settings.template.json" "$repo/.vscode/settings.json" \
                    "$repo" "$lint_case" <<'PY'
import json, re, sys
from pathlib import Path

tpl, out, repo, want_case = sys.argv[1:5]
repo = Path(repo)

# A build directory is recognised by its generated config.mk, which records the
# exact fypp defines the build used, next to the f90/ directory holding its
# .mod files. Reading both from the same place is what keeps them consistent.
def builds(root):
    # "build/latest" is a symlink to the newest hashed build dir, so the glob
    # sees each build twice; keep one entry per real directory.
    seen = {}
    for cfg in root.glob("**/build/*/config.mk"):
        if (cfg.parent / "f90").is_dir():
            seen.setdefault(cfg.parent.resolve(), cfg)
    return list(seen.values())


def module_dir(build):
    """Prefer <case>/build/latest/f90: it keeps pointing at the newest build
    of that case, so a rebuild that changes the flag hash needs no re-run."""
    latest = build.parent / "latest"
    if latest.is_dir() and latest.resolve() == build.resolve():
        return str(latest / "f90")
    return str(build / "f90")

if want_case:
    case = Path(want_case)
    case = case if case.is_absolute() else repo / case
    if not case.is_dir():
        sys.exit(f"error: no such directory: {case}")
    found = sorted(builds(case), key=lambda c: c.stat().st_mtime, reverse=True)
    if not found:
        sys.exit(f"error: no build under {case}. Run: cd {case} && make arch=gnu")
else:
    found = sorted(builds(repo), key=lambda c: c.stat().st_mtime, reverse=True)

# config.mk emits one of three shapes per define (make/config_reader.py):
#   fypp_flags += -DNAME                NAME=""      -> fypp gets -DNAME
#   fypp_flags += -DNAME=\'value\'       NAME="'value'"
#   fypp_flags += -DNAME=12             NAME="12"
DEFINE = re.compile(r"^fypp_flags \+= -D(\w+)(?:=(.*))?$")

if found:
    cfg = found[0]
    defines = {}
    for line in cfg.read_text().splitlines():
        m = DEFINE.match(line.strip())
        if m:
            defines[m.group(1)] = (m.group(2) or "").replace("\\", "")
    linter, mod_dirs = "gfortran", [module_dir(cfg.parent)]
    case_name = cfg.parent.parent.parent.relative_to(repo)
else:
    # Nothing to lint against. Leave both empty rather than naming a build
    # that does not exist: the next run fills them in from a real one.
    defines, linter, mod_dirs, case_name = {}, "Disabled", [], ""

def block(value, indent="  "):
    return json.dumps(value, indent=2).replace("\n", "\n" + indent)

text = (Path(tpl).read_text()
        .replace("@AGILE_DIR@", str(repo))
        .replace("@LINTER@", linter)
        .replace("@DEFINITIONS@", block(defines))
        .replace("@MODULE_DIRS@", block(mod_dirs)))
Path(out).parent.mkdir(parents=True, exist_ok=True)
Path(out).write_text(text)

print(f"{linter}\t{case_name}\t" + " ".join(f"{k}={v}" for k, v in defines.items()))
PY
)"

linter="${summary%%$'\t'*}"
rest="${summary#*$'\t'}"
case_name="${rest%%$'\t'*}"
defines="${rest#*$'\t'}"

if [ "$linter" = "gfortran" ]; then
  say "diagnostics" "gfortran on fypp output"
  say "" "configuration $defines"
  say "" "from $case_name"
else
  say "diagnostics" "off -- no build found in the tree"
  printf '\n  To switch them on, build any case, e.g.\n'
  printf '    cd %s/tests/hd/KH3D && make arch=gnu\n' "$repo"
  printf '  then re-run this script.\n'
fi

printf '\n  Wrote %s. Reload the VS Code window to pick everything up.\n' "$repo/.vscode/settings.json"
