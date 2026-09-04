# CCP4i2 - Django Branch

CCP4i2 provides an environment for crystallographic computing. This branch (ccp4i2-django) implements a Qt-free architecture using Django for the backend and React/Electron for the frontend.

## Architecture

### Backend
- **Django REST API** (`server/`) - Database interactions via ORM, REST API with Django REST Framework
- **Core Modules** (`server/ccp4i2/core/`) - Business logic, data containers, task management

### Frontend
- **Electron/React App** (`client/`) - Next.js 15, React 19, Moorhen structure viewer

### Compounds App (extracted)
The compounds registry/assays app (`apps/compounds/`, backend + frontend) was
**extracted to the `newcastleuniversity/materia` repository** in the 2026-04-29
cut (commit `a444f56db`); it is no longer developed here. A few client-side
compounds stubs remain under `client/renderer/lib/compounds/` (auth/rdkit/theme
contexts), and the shared API contract it consumes lives in
`packages/ccp4i2-api/`. For the compounds app itself, see the Materia repo.

## Key Directories

The Python package lives at `server/ccp4i2/`; paths below are given in full,
because most of these directories were top-level in earlier layouts and moved
under the package when `core/` was relocated (commit `e92d6b7bf`).

| Directory | Purpose |
|-----------|---------|
| `server/` | Django project root (`manage.py`, settings) |
| `server/ccp4i2/core/` | Core Python modules and business logic |
| `server/ccp4i2/api/` | REST API views and serializers |
| `server/ccp4i2/db/` | Models, migrations, job execution, gleaner |
| `server/ccp4i2/lib/` | Shared utilities (containers, jobs, parameters, reporting) |
| `server/ccp4i2/cli/` | Command-line tools (i2run, utilities) |
| `server/ccp4i2/wrappers/` | Task wrappers for crystallographic programs |
| `server/ccp4i2/pipelines/` | Multi-step crystallographic pipelines |
| `server/ccp4i2/pimple/` | Matplotlib-based graph generation |
| `server/ccp4i2/smartie/` | Log parsing utilities |
| `server/ccp4i2/report/` | Report generation and parsing |
| `server/ccp4i2/tests/` | Test suite |
| `server/ccp4i2/demo_data/` | Small datasets used by tests and demos |
| `client/` | Electron/React frontend for desktop |
| `packages/ccp4i2-api/` | Shared API contract (npm `@ccp4/ccp4i2-api`, PyPI `ccp4i2-api`) |
| `docs/` | Developer and user documentation |
| `deploy/` | Deployment tooling (test VMs, infra scripts) |

There is no longer a top-level `apps/` directory: it was emptied by the Materia
cut and removed in commit `1f8673cca`.

## Windows Compatibility: No Unicode in print()

**Do not use emoji or non-ASCII characters in Python `print()` statements** in runtime code (anything outside `tests/`). On Windows, Python's console output defaults to cp1252 encoding, which cannot encode emoji (e.g. `✅`, `❌`, `⚠️`). A `UnicodeEncodeError` from `print()` inside a `try/except Exception` block will be silently caught and can cause jobs to fail — this was the root cause of aimless_pipe failing on Windows while working on mac/linux.

- `print()` — ASCII only (use `[OK]`, `ERROR:`, `WARNING:` etc.)
- `logger.debug()` / `logger.info()` — emoji OK (log files use UTF-8)
- Test files — emoji OK (controlled environment)

## Environment Detection

The system automatically detects which backend to use:

1. **Explicit**: Set `CCP4I2_BACKEND=django` or `CCP4I2_BACKEND=qt`
2. **Django context**: Presence of `DJANGO_SETTINGS_MODULE`
3. **Auto-detection**: Try importing PySide2; if unavailable, use Django mode

## CCP4 Environment Setup

**IMPORTANT**: CCP4i2 requires the CCP4 suite with `ccp4-python` to run. Before running any Python commands:

```bash
# Source the CCP4 setup script (adjust path as needed)
source ../ccp4-20260702/bin/ccp4.setup-sh

# This sets up:
# - ccp4-python interpreter with all dependencies (gemmi, clipper, etc.)
# - CCP4 environment variables
# - CCP4 binaries in PATH
```

All Python commands below should use `ccp4-python` instead of `python`.

### Inspecting PHIL-based Tool Parameters

When wrapping a PHIL-based tool (Phenix, PhaserTNG, DIALS, etc.), inspect its `master_phil` to understand parameter names, types, and multiplicity:

```bash
source ../ccp4-20260702/bin/ccp4.setup-sh

# Dump the full PHIL tree
ccp4-python -c "from my_tool import master_phil; master_phil.show()"

# For tools needing custom PHIL types (e.g., PhaserTNG):
ccp4-python -c "
from phasertng.programs import picard
from iotbx.cli_parser import CCTBXParser
parser = CCTBXParser(program_class=picard.Program, logger=None, parse_phil=False)
parser.master_phil.show()
"

# Check if a parameter is repeatable (.multiple = True):
ccp4-python -c "
from phasertng.programs import picard
from iotbx.cli_parser import CCTBXParser
parser = CCTBXParser(program_class=picard.Program, logger=None, parse_phil=False)
for obj in parser.master_phil.all_definitions():
    if 'xyzin' in obj.path:
        print(f'{obj.path}  type={obj.object.type}  multiple={obj.object.multiple}')
"

# Show a scope with full attributes (multiple, expert_level, type, constraints):
ccp4-python -c "
...
scope = parser.master_phil.get('picard')
scope.show(attributes_level=2)
"
```

Parameters with `.multiple = True` need `CList` + `PdbFileListShim` in the wrapper, not a single file type. See `server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md` for the full guide.

## Running

### Django Mode (i2run CLI)
`i2run` is a Django management command. Run it from `server/`:
```bash
cd server
export CCP4I2_BACKEND=django
export DJANGO_SETTINGS_MODULE=ccp4i2.config.settings
ccp4-python manage.py i2run <task> --project_name <proj> [--PARAM value ...]
```

**Never run pytest from a shell that has that export.** On 2026-09-04 the API
unit suite was run with `DJANGO_SETTINGS_MODULE=ccp4i2.config.settings` still
exported. pytest-django ranks the environment variable above the setting in
`pytest.ini`, so the suite ran on the production settings, and its transactional
teardown flushed the developer's live `~/.ccp4i2-django/db.sqlite3` — every
project and job — while the desktop app was running on it. Recovery was
`manage.py restore_projects --scan <projects dir>` from the per-project
`DATABASE.db.xml` snapshots. The suite now defends itself (`--ds` in
`pytest.ini` outranks the variable, and `ccp4i2/tests/conftest.py` refuses any
other settings module or a database under a real user home), but keep the
habit: set the variable per command, `env DJANGO_SETTINGS_MODULE=... ccp4-python
manage.py i2run ...`, rather than exporting it.

### Development Server
```bash
cd server
ccp4-python manage.py runserver
```

### Frontend
```bash
cd client
npm install
npm run dev
```

### Tests

All commands run from the `server/` directory using `ccp4-python -m pytest` (cross-platform, works on Windows/Linux/macOS).

#### Test Layout

Tests are organised by what they test and how fast they run:

```
ccp4i2/tests/
├── unit/                        # Fast (<30s total), no CCP4 binaries needed
│   ├── containers/              # CContainer, CData, CList, CDataFile, types
│   ├── mtz/                     # MTZ join/split/columns/conversions (needs gemmi)
│   ├── pdb/                     # PDB/mmCIF loading
│   ├── phil/                    # PHIL parameter system
│   ├── converters/              # Data format converters
│   ├── plugins/                 # Plugin infrastructure, def.xml, registry
│   ├── validation/              # Validity checks, error reporting
│   ├── serialization/           # JSON/XML encoding
│   └── lib/                     # Utilities, reports, sequences, uploads
├── parity/                      # Native ports vs the CCP4 binary they
│                                #   replaced (needs the CCP4 binaries)
├── async/                       # Async execution infrastructure
├── db/                          # Database, project import/export
├── api/
│   ├── unit/                    # REST endpoint tests (Django test client)
│   └── e2e/                     # End-to-end pipeline tests via REST API
└── i2run/                       # End-to-end task tests via CLI
```

#### Reproducing CI locally: the CCP4-free venv

CI runs `tests/unit/` and `tests/api/unit/` on **stock Python with no CCP4**, to
prove those suites need none. Running them under `ccp4-python` gives strictly
*more* coverage and so cannot prove it: a test that needs CCP4 passes locally
and fails in CI. That has happened.

```bash
cd server
python3 -m venv .venv                 # .gitignore already expects .venv/
.venv/bin/pip install -e ".[dev]"     # exactly what CI installs

.venv/bin/python -m pytest ccp4i2/tests/unit/ -q
.venv/bin/python -m pytest ccp4i2/tests/api/unit/ -q
```

Run the two suites as **separate** commands — sharing one pytest process makes
the report-fixture tests fail against the API suite's database (the CI workflow
explains why at length).

Anything needing CCP4 must skip rather than fail: guard with
`pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")` or a
`skipif`. Better still, test the CCP4-free part against something stdlib so it
runs everywhere.

#### Running Tests

```bash
cd server
source ../../ccp4-20260702/bin/ccp4.setup-sh

# ── Fast unit tests (no CCP4 binaries, good for CI on all platforms) ──
ccp4-python -m pytest ccp4i2/tests/unit/ -v

# By subsystem
ccp4-python -m pytest ccp4i2/tests/unit/mtz/ -v
ccp4-python -m pytest ccp4i2/tests/unit/containers/ -v
ccp4-python -m pytest ccp4i2/tests/unit/plugins/ -v

# ── Parity with the CCP4 binaries we replaced (needs those binaries) ──
# Run after touching any gemmi/numpy port of a CCP4 program.
ccp4-python -m pytest ccp4i2/tests/parity/ -v

# ── End-to-end task tests via CLI (slow, needs CCP4) ──
ccp4-python -m pytest ccp4i2/tests/i2run/ -v
ccp4-python -m pytest ccp4i2/tests/i2run/test_servalcat.py -v
ccp4-python -m pytest ccp4i2/tests/i2run/test_servalcat.py::test_servalcat_basic -v
ccp4-python -m pytest ccp4i2/tests/i2run/ -k "test_aimless"

# ── End-to-end pipeline tests via API (slow, needs CCP4) ──
ccp4-python -m pytest ccp4i2/tests/api/e2e/ -v

# ── API unit tests ──
ccp4-python -m pytest ccp4i2/tests/api/unit/ -v

# ── Everything ──
ccp4-python -m pytest ccp4i2/tests/ -v
```

#### Testing Philosophy

| Layer | What it tests | Speed | CCP4 needed? | When to run |
|-------|--------------|-------|-------------|-------------|
| `unit/` | Core data classes, gemmi utilities, MTZ/PDB operations, PHIL, converters, plugin infrastructure, validation, serialization | Fast (~5s) | No (just ccp4-python) | Every commit, CI on all platforms |
| `async/` | Async execution framework | Fast | No | With unit tests |
| `db/` | Database operations, project import/export | Medium | No | When touching DB code |
| `api/unit/` | REST endpoint behaviour | Fast | No | When touching API code |
| `parity/` | Native gemmi/numpy implementations against the CCP4 binary they replaced (freerflag, chltofom, matthews, clipper cell checks) | Fast | Yes — the binaries | After changing any gemmi-native port |
| `api/e2e/` | Full pipeline execution via REST API | Slow | Yes | Before release, after pipeline changes |
| `i2run/` | Full task execution via CLI | Slow | Yes | Before release, after wrapper changes |

**Key principles:**
- The suite runs on `ccp4i2.config.test_settings` and nothing else. `pytest.ini`
  forces it with `--ds`, and `ccp4i2/tests/conftest.py` aborts the session if
  the settings module is any other or the database path lies under
  `~/.ccp4i2x`, `~/.ccp4i2-django` or `~/.ccp4i2`. Do not weaken either: the
  API suite flushes whatever database it runs on (see the i2run section above)
- Unit tests must not depend on CCP4 binaries — only `ccp4-python` (for gemmi etc.)
- Unit tests must not depend on external data (test101, ProjectZips). Use `demo_data/` from the repo or `I2_TOP` for paths
- E2e tests (i2run, api/e2e) download test data via session-scoped fixtures from PDBe/RCSB/PDB-REDO
- All tests use `ccp4-python -m pytest` (not `run_test.sh`) for cross-platform consistency
- Failed test directories are preserved in `~/.cache/ccp4i2-tests/` for debugging.
  They are never pruned, so this grows: **2.1 GB across 133 directories** after a
  few weeks of runs. Delete the old ones when disk gets tight
- Downloaded fixtures are cached in `~/.cache/ccp4i2-tests/downloads/`, so a run
  does not depend on every external host being up. Files over 100 MB are *not*
  cached — the xia2 fixture pulls a 335 MB archive to use twenty frames of it.
  `CCP4I2_TEST_CACHE_MAX_MB=400` includes it, `=0` removes the limit, and
  `CCP4I2_TEST_REFETCH=1` ignores what is cached
- **The parity tier guards the slim-server ports.** Replacing a CCP4 binary with
  gemmi/numpy is only as good as the evidence it behaves the same, and that
  evidence goes stale quietly: the freerflag parity tests sat in the tree for two
  months while *no documented command collected them*, and in that time the
  completion path diverged from the binary unnoticed. Run this tier whenever you
  touch a port
- Each test gets an isolated SQLite database and project directory

### Compounds App Tests

The compounds app (and its test suite) was extracted to the
`newcastleuniversity/materia` repository in the 2026-04-29 cut — its tests no
longer live here. Run them from a Materia checkout, not this repo.

## Task Validation: `validity()` and `runTimeValidity()`

The server is the sole authority for job validation. The frontend displays what the server reports — there is no client-side validation logic.

### Two-tier architecture

| Method | When called | Speed requirement | Example checks |
|--------|------------|-------------------|----------------|
| `validity()` | Polled during parameter editing (`GET /validation/`) | Fast — no file I/O | Required files set, contentFlag match, cross-parameter logic |
| `runTimeValidity()` | Once at submission (`GET /run_time_validation/`) and by `process()` | Can be expensive | Monomer dictionary coverage (reads files with gemmi) |

The base `validity()` automatically checks `allowUndefined`, `mustExist`, `requiredContentFlag` on all container children, and filters out `outputData` errors (stale from cloned jobs).

### Overriding `validity()` for task-specific checks

```python
def validity(self):
    error = super(my_task, self).validity()

    # Non-blocking advisory (orange in UI, Confirm stays enabled)
    if not self.container.inputData.FREERFLAG.isSet():
        error.append(klass=self.TASKNAME, code=200,
            details='Free R flag is strongly recommended',
            name=f'{self.TASKNAME}.container.inputData.FREERFLAG',
            severity=CCP4ErrorHandling.SEVERITY_WARNING)

    # Blocking error (red in UI, Confirm disabled)
    if str(self.container.inputData.COMP_BY) == 'ASU':
        if not self.container.inputData.ASUFILE.isSet():
            error.append(klass=self.TASKNAME, code=201,
                details='ASU file required when COMP_BY is ASU',
                name=f'{self.TASKNAME}.container.inputData.ASUFILE',
                severity=CCP4ErrorHandling.SEVERITY_ERROR)
    return error
```

### Overriding `runTimeValidity()` for expensive pre-flight checks

```python
def runTimeValidity(self):
    error = super(my_task, self).runTimeValidity()
    if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
        return error  # skip expensive checks if already failing
    # ... expensive checks (e.g. checkMonomeCoverage) ...
    return error
```

### Severity mapping

| Backend constant | Value | Frontend | Run dialog |
|-----------------|-------|----------|------------|
| `SEVERITY_WARNING` | 2 | Orange indicator | Advisory, does not block |
| `SEVERITY_ERROR` | 4 | Red indicator | Blocks Confirm button |

### Error `name` field format

Must be `{taskName}.container.{section}.{fieldName}` (e.g. `servalcat_pipe.container.inputData.FREERFLAG`) to enable field-level error display in the frontend. Check the def.xml to confirm which section (`inputData` or `controlParameters`) the field is in.

## Core Data Model (CData)

Every task parameter is a `CData` object — Django-shaped declarative classes
(`content()` fields, `class Meta:`) with a per-instance qualifier override
layer. Before adding or changing a CData class, read
[docs/cdata.md](docs/cdata.md): it covers the syntax, the semantics
(set/unset/default states, assignment coercion, the two qualifier layers), the
class-construction machinery, and a how-to for new classes with the known
traps (`__setattr__` coercion, `CBoolean.__bool__`, the CString hash mismatch,
the def.xml CString fallback for unresolvable class names).

## Task Registry

Tasks are registered in a single static module, `server/ccp4i2/core/tasks.py`. It
holds a `TASKS` dict keyed by task name, each value a `Task` dataclass:

```python
@dataclass
class Task:
    title: str = None
    description: str = None
    shortTitle: str = None
    autogenerated: bool = False
    pluginPath: str = None      # "ccp4i2.wrappers.acorn.script.acorn:acorn"
    defXmlPath: str = None      # "wrappers/acorn/script/acorn.def.xml"
    reportPath: str = None      # optional "...:acorn_report"
    runningReport: bool = False
    watchedFile: str = None
```

Accessors (lazy, cached — they import the plugin module on first use) live in the
same module:

```python
from ccp4i2.core.tasks import (
    TASKS, get_plugin_class, get_report_class, get_task_title, locate_def_xml,
)

list(TASKS)                       # all registered task names
acorn_class = get_plugin_class('acorn')   # the CPluginScript subclass
def_xml = locate_def_xml('acorn')         # path to its .def.xml
```

### Adding a new task

1. Create `server/ccp4i2/wrappers/<Name>/script/<Name>.py` (a `CPluginScript`
   subclass) and `<Name>.def.xml`, plus empty `__init__.py` in the task dir and
   its `script/`.
2. Add one entry to the `TASKS` dict in `server/ccp4i2/core/tasks.py`
   (`pluginPath` + `defXmlPath`; `reportPath` only if you wrote a report class).
3. (Optional) Register a frontend interface in
   `client/renderer/components/task/task-interfaces/task-container.tsx` and add
   the task name to a category in
   `client/renderer/components/task/task-chooser.tsx` so it appears in the task
   chooser. Tasks with no registered interface fall back to `GenericInterface`,
   which renders the def.xml automatically.

There is **no registry-regeneration step** — adding the `TASKS` entry is the
whole registration. (The historical `plugin_registry.py` / `plugin_lookup.json`
scan and the `CCP4Modules.TASKMANAGER()` accessor were removed when task data was
consolidated into `tasks.py`.) Output files are persisted to the project by the
gleaner (`server/ccp4i2/db/async_db_handler.py: glean_job_files`) for any
`outputData` `CDataFile` that is set and exists on disk — there is no `saveToDb` gate.

### Making a task's data files exportable

To let the job panel's **Export MTZ** button offer a task's principal reflection
file, either declare a tracked `outputData` `CMtzDataFile` (preferred — the
generic fallback then offers it automatically) or implement the module-level
`exportJobFileMenu` / `exportJobFile` contract (for pipelines needing subjob
lookup or reconstruction). Never shell out to `cad` — use
`server/ccp4i2/lib/utils/jobs/export.py: combine_mtz_files` for gemmi-based
joins. Full guide: `server/ccp4i2/wrappers/EXPORT_TASK_GUIDE.md`.

## Files Preserved from Legacy

- `server/ccp4i2/wrappers/`, `server/ccp4i2/pipelines/` - Crystallographic logic
- `server/ccp4i2/pimple/`, `server/ccp4i2/smartie/` - Utility modules
- `client/assets/qticons/`, `client/assets/svgicons/` - Task icons (served from
  `client/renderer/public/` by the `copy-qticons` / `copy-svgicons` npm scripts)
- `server/ccp4i2/tipsOfTheDay/` - User tips
- `docs/` - Documentation

## Deployment (Docker / Azure)

> **The Docker build and Azure deployment tooling no longer lives in this repo.**
> It was extracted to the `newcastleuniversity/materia` repository as part of the
> compounds-app cut (commit `eb130b163`, "Post-cut tidy: Docker/ deleted"). The
> old `Docker/` directory — the layered CCP4 base/ARP-warp images, the
> compounds-overlay web image, and the `azure-uksouth/` container-app deploy
> scripts — all moved with it, because those images are the compounds/DDU-Database
> deployment rather than a CCP4i2-core concern. Look in the Materia repo for the
> current build-and-push / deploy-applications scripts and `.env.deployment`
> files.
>
> What remains in this repo under `deploy/` is CCP4i2's own infra tooling — the
> standing Windows/Linux test VMs and their bicep/cloud-init.

### TODO: lab-wide, Docker-based CCP4i2 deployment

There is currently **no** self-contained Docker/Compose deployment path in this
repo for standing up CCP4i2 (server + worker + web) for a lab or group. The
desktop Electron app (pip-into-`ccp4-python`) is the only first-class install
route here. Rebuilding a minimal, documented lab-wide Docker deployment —
independent of the compounds overlay, ideally driven by the same
`ccp4-python`-based server image — is a worthwhile future addition.
