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

Linting needs a build, for two things: the fypp defines and the `.mod` files
that `use` statements resolve against. The script takes both from one build's
generated `config.mk` and its `f90/` directory, so they always describe the
same configuration. By default it picks the most recently built case anywhere
in the tree. Call `setup.sh` with a `--case` to choose a specific build for the linter to adapt:

```bash
tools/vscode/setup.sh --case tests/mhd/spherical_pole
```

The chosen configuration is printed, e.g.

```
  diagnostics   gfortran on fypp output
                configuration PHYS='hd' SPECIALBOUNDARY= REFINE_USR= GEOM='spherical'
                from tests/hd/spherical_pole
```

Two consequences, both inherent to a compile-time-configured code:

- **One configuration at a time.** `PHYS`, `GEOM` and the feature flags are
  fypp defines, so a lint pass only ever preprocesses one branch. Code behind
  the others is not checked — not wrongly flagged, just unvisited. Build the
  case you are working on and re-run the script, and the branch you are
  editing is the one that gets checked, that case's own `mod_usr.fpp`
  included.
- **With no build, diagnostics stay off.** Without `.mod` files every file
  would report "Cannot open module file", so the script leaves
  `fortran.linter.compiler` at `Disabled` and says which command fixes it.

The module path points at `build/latest`, which follows the newest build of
that case, so re-running the script is only needed when switching cases, not
after an ordinary rebuild.

## Troubleshooting

| symptom | cause |
|---|---|
| "It is highly recommended to use the fortls..." | `fortls` not found. `.venv` was rebuilt without the `ide` extra, or the path in settings is stale — re-run the script. |
| "Installing fypp through pip with --user option" | Same, for `fypp`. The extension does *not* expand `${workspaceFolder}` for `fortran.linter.fypp.path`, which is why the script writes it absolute. |
| Every `#:` line flagged `invalid preprocessing directive` | `fortran.linter.fypp.enabled` is off, so gfortran is reading raw `.fpp`. |
| fypp constructs unhighlighted | The grammar extension is not installed or the window was not reloaded. Check with `code --list-extensions \| grep fypp`. |
