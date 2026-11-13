#!/usr/bin/env bash
#
# Set up VS Code for this repository. Idempotent -- re-run it after rebuilding
# .venv, or after building the lint case for the first time.
#
#   tools/vscode/setup.sh [--no-extension]
#
# It does three things:
#   1. makes sure .venv has fypp and fortls (uv sync --extra ide)
#   2. installs the fypp syntax grammar as a local VS Code extension
#   3. writes .vscode/settings.json from settings.template.json
#
# See tools/vscode/README.md.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
lint_case="tests/hd/KH3D"   # the build whose .mod files gfortran lints against

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
if [ "${1:-}" = "--no-extension" ]; then
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

mod_dir="$repo/$lint_case/build/latest/f90"
if [ -d "$mod_dir" ]; then linter="gfortran"; else linter="Disabled"; fi

python3 - "$here/settings.template.json" "$repo/.vscode/settings.json" \
         "$repo" "$linter" "$mod_dir" <<'PY'
import sys
from pathlib import Path

tpl, out, repo, linter, mod_dir = sys.argv[1:6]
text = (Path(tpl).read_text()
        .replace("@AGILE_DIR@", repo)
        .replace("@LINTER@", linter)
        .replace("@MODULE_DIR@", mod_dir))
Path(out).parent.mkdir(parents=True, exist_ok=True)
Path(out).write_text(text)
PY

if [ "$linter" = "gfortran" ]; then
  say "diagnostics" "gfortran on fypp output, modules from $lint_case"
else
  say "diagnostics" "off -- no build found under $lint_case"
  printf '\n  To switch them on:  cd %s/%s && make arch=gnu\n' "$repo" "$lint_case"
  printf '  then re-run this script.\n'
fi

printf '\n  Wrote %s. Reload the VS Code window to pick everything up.\n' "$repo/.vscode/settings.json"
