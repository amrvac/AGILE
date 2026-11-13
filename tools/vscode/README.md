# VS Code setup

AGILE's sources are `.fpp`: Fortran plus [fypp](https://fypp.readthedocs.io)
preprocessor syntax. Out of the box an editor gets both halves of that wrong —
the syntax highlighting does not know fypp, and the linter tries to compile a
file that no compiler can read. This folder fixes both.

```bash
tools/vscode/setup.sh
```

Then reload the VS Code window. The script is idempotent; re-run it after
rebuilding `.venv` or after building the lint case for the first time.

It needs the [Modern Fortran](https://marketplace.visualstudio.com/items?itemName=fortran-lang.linter-gfortran)
extension (`fortran-lang.linter-gfortran`), which `.vscode/extensions.json`
recommends — VS Code offers to install it when you first open the repository.

## What the script does

**1. Python environment.** Runs `uv sync --extra ide` if `.venv` lacks `fypp`
or `fortls`. The `ide` extra in `pyproject.toml` exists only to add `fortls`.

**2. Syntax highlighting.** Packages `fypp-fortran/` into a `.vsix` with
`python3` (nothing to compile) and installs it. It contributes one TextMate
*injection* grammar, layered onto the Fortran grammar, so Fortran keeps being
highlighted as Fortran and only fypp is added:

| construct | example |
|---|---|
| comment | `#! note` |
| directive | `#:if GEOM != 'Cartesian'` |
| eval, line form | `$:FACE_FLUX(d, 0, 'fx', proc)` |
| direct call | `@:addsource_local()` |
| eval, inline | `${GPU_ROUTINE_SEQ()}$` |
| directive, inline | `F#{if defined('FLUX_AD')}#, phi#{endif}#` |

Inside a directive, strings, numbers, builtins, macro names and the upper-case
build defines (`PHYS`, `GEOM`, `LOG_RADIUS`, ...) are highlighted separately,
and a directive continued with a trailing `&` stays highlighted across the
continuation. Every rule is bounded to one line, so an unbalanced quote in a
directive cannot bleed its colour down the file.

To remove it: `code --uninstall-extension agile.fypp-fortran`.

**3. `.vscode/settings.json`,** rendered from `settings.template.json`. It is
gitignored, because two of its paths have to be absolute. Edit the template and
re-run the script rather than editing it directly.

## Diagnostics

`fortls` provides go-to-definition, hover, peek and the outline. It cannot
parse fypp directives, so its own diagnostics are turned off and the real ones
come from gfortran, which the extension runs *after* fypp — fypp's line markers
map any error back to the `.fpp` line.

Two limits follow from that, both of them inherent:

- **One configuration is linted.** `PHYS` and `GEOM` are fypp defines, so a
  lint pass only sees one branch of the tree. The template pins
  `PHYS='hd'`, `GEOM='Cartesian'`. To check another, change
  `fortran.linter.fypp.definitions` *and* point `fortran.linter.includePaths`
  at a build of that same configuration — `.mod` files are per-configuration
  too.
- **It needs a build.** `use` statements resolve against the `.mod` files of
  `tests/hd/KH3D`; without them every file reports "Cannot open module file".
  The script detects this and leaves `fortran.linter.compiler` at `Disabled`
  until you have run `make arch=gnu` there.

## Troubleshooting

| symptom | cause |
|---|---|
| "It is highly recommended to use the fortls..." | `fortls` not found. `.venv` was rebuilt without the `ide` extra, or the path in settings is stale — re-run the script. |
| "Installing fypp through pip with --user option" | Same, for `fypp`. The extension does *not* expand `${workspaceFolder}` for `fortran.linter.fypp.path`, which is why the script writes it absolute. |
| Every `#:` line flagged `invalid preprocessing directive` | `fortran.linter.fypp.enabled` is off, so gfortran is reading raw `.fpp`. |
| fypp constructs unhighlighted | The grammar extension is not installed or the window was not reloaded. Check with `code --list-extensions \| grep fypp`. |
