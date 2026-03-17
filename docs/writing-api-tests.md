# Writing API Tests for CCP4i2

This guide covers how to write API-driven integration tests for CCP4i2 task
wrappers and pipelines.  The infrastructure lives under
`server/ccp4i2/tests/api/` and mirrors the GUI workflow: create a project,
create a job, upload files, set parameters, run the job, and validate outputs.

Companion guide: [docs/writing-i2run-tests.md](writing-i2run-tests.md)
Coverage tracker: [docs/test-coverage-tracker.md](test-coverage-tracker.md)

---

## How API tests differ from i2run tests

| Aspect | i2run tests | API tests |
|--------|-------------|-----------|
| **Entry point** | `i2run(args)` context manager | Django REST API via `APIClient` |
| **Parameter syntax** | CLI flags: `--NCYCLES 2` | Object paths: `set_param("controlParameters.NCYCLES", 2)` |
| **File input** | File paths on CLI | `upload_file()` / `upload_file_with_columns()` |
| **Job execution** | Synchronous, in-process | `run_job()` + `wait_for_completion()` polling |
| **Database** | Isolated SQLite per test (via conftest) | Same — isolated SQLite per test (via conftest) |
| **What it exercises** | Wrapper/pipeline logic only | Full stack: API routing, serialisation, parameter setting, file upload, job runner, output registration |
| **Class structure** | Bare functions | Classes inheriting `APITestBase` |
| **Test data** | `demoData()` helper | Pytest fixtures (`gamma_mtz`, `cif8xfm`, etc.) |

The key advantage of API tests is that they exercise the same code path the
frontend uses, catching issues in serialisation, permission checks, file
registration, and the REST layer that i2run tests cannot reach.

---

## Quick start

```python
# server/ccp4i2/tests/api/test_mytask_api.py
import pytest
from .base import APITestBase

pytestmark = pytest.mark.pipeline


@pytest.mark.usefixtures("file_based_db")
class TestMyTaskAPI(APITestBase):
    task_name = "mytask"
    timeout = 120

    def test_basic(self, gamma_mtz, gamma_model_pdb):
        self.create_project("test_mytask")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file("inputData.XYZIN", gamma_model_pdb)
        self.set_param("controlParameters.NCYCLES", 2)

        self.run_and_wait()

        self.assert_file_exists("XYZOUT.cif")
```

Run it:

```bash
cd server
./run_test.sh ccp4i2/tests/api/test_mytask_api.py -v
```

---

## How the harness works

Each API test inherits from `APITestBase` (defined in `base.py`) and follows
a fixed lifecycle:

1. **`conftest.py`** automatically creates an isolated SQLite database and
   project directory per test (`isolated_test_db` fixture, autouse)
2. **`setup_client`** (autouse) gives you `self.client` — a DRF `APIClient`
3. You call `self.create_project(name)` → POST to `/api/ccp4i2/projects/`
4. You call `self.create_job()` → POST to `/api/ccp4i2/projects/{id}/create_task/`
5. Upload files and set parameters via helper methods
6. `self.run_and_wait()` → POST to `/api/ccp4i2/jobs/{id}/run/` then polls
   the database until a terminal status is reached
7. On success, the project directory is cleaned up; on failure it is preserved

Permission checks are bypassed via the `bypass_api_permissions` autouse
fixture — all viewsets get `AllowAny` during tests.

---

## Test structure

### Class boilerplate

Every API test file follows this pattern:

```python
import pytest
from .base import APITestBase

# Mark all tests as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline


@pytest.mark.usefixtures("file_based_db")
class TestMyPipelineAPI(APITestBase):
    """API tests for my_pipeline."""

    task_name = "my_pipeline"   # Must match plugin registry name
    timeout = 300               # Max seconds to wait for job completion
    # poll_interval = 2         # Seconds between status checks (default 2)

    def test_case_one(self, gamma_mtz, gamma_model_pdb):
        ...

    def test_case_two(self, cif8xfm, mtz8xfm):
        ...
```

Key points:
- `@pytest.mark.usefixtures("file_based_db")` — required for subprocess job
  execution (provides the database path via environment variable)
- `task_name` — must exactly match the plugin name in the registry
- `timeout` — adjust per task; EP and MR pipelines need 300-600s
- One class per task, multiple test methods per class

---

## Test data sources

### Local demo data fixtures

The `conftest.py` provides session-scoped fixtures for common demo data files.
These are preferred for fast, reliable tests:

```python
def test_basic(self, gamma_mtz, gamma_model_pdb):
    ...

def test_anomalous(self, gamma_mtz, gamma_freer_mtz, gamma_asu_xml):
    ...

def test_mdm2(self, mdm2_unmerged_mtz, mdm2_model_cif):
    ...
```

Available fixtures (see `conftest.py` for full list):

| Fixture | File | Description |
|---------|------|-------------|
| `gamma_mtz` | `gamma/merged_intensities_Xe.mtz` | Anomalous intensities |
| `gamma_native_mtz` | `gamma/merged_intensities_native.mtz` | Native intensities |
| `gamma_unmerged_mtz` | `gamma/gamma_native.mtz` | Unmerged reflections |
| `gamma_model_pdb` | `gamma/gamma_model.pdb` | Search model |
| `gamma_freer_mtz` | `gamma/freeR.mtz` | Free R flags |
| `gamma_phases_mtz` | `gamma/initial_phases.mtz` | Initial phases |
| `gamma_asu_xml` | `gamma/gamma.asu.xml` | ASU contents |
| `gamma_heavy_atoms_pdb` | `gamma/heavy_atoms.pdb` | Heavy atom positions |
| `mdm2_unmerged_mtz` | `mdm2/mdm2_unmerged.mtz` | MDM2 unmerged data |
| `mdm2_model_cif` | `mdm2/4hg7.cif` | MDM2 search model |
| `rnase_model_pdb` | `rnase/rnase_model.pdb` | RNase model (for mismatch tests) |
| `demo_data_dir` | `demo_data/` | Root directory Path |

### Downloaded data fixtures

Session-scoped fixtures that download from PDB/PDB-REDO (cached per session):

| Fixture | Source | Description |
|---------|--------|-------------|
| `cif8xfm` | PDB-REDO | 8xfm refined mmCIF |
| `mtz8xfm` | PDB-REDO | 8xfm refined MTZ (FP, SIGFP, FREE) |
| `seq8xfm` | PDBe | 8xfm FASTA sequence |
| `cif7beq` | RCSB | 7beq mmCIF |
| `mtz7beq` | PDB-REDO | 7beq refined MTZ |

### Adding new fixtures

If your test needs data not covered above, add a new session-scoped fixture
to `conftest.py`:

```python
from ccp4i2.tests.i2run.urls import redo_mtz, pdbe_pdb
from ccp4i2.tests.i2run.utils import download

@pytest.fixture(scope="session")
def mtz1abc():
    """Download 1abc MTZ from PDB-REDO."""
    with download(redo_mtz("1abc")) as path:
        yield path
```

For one-off downloads inside a test method, use the `download` and `URLs`
helpers from `base.py`:

```python
from .base import download, URLs

def test_with_extra_data(self):
    with download(URLs.redo_cif("7ber")) as cif_path:
        self.upload_file("inputData.REFERENCE", cif_path)
        ...
```

---

## Setting parameters

### Simple parameters

```python
self.set_param("controlParameters.NCYCLES", 2)
self.set_param("controlParameters.USEANOMALOUS", True)
self.set_param("controlParameters.RESOLUTION_HIGH", 2.0)
self.set_param("inputData.EXPTYPE", "SIRAS")
```

The `set_param` method automatically prepends `{task_name}.container.` to the
object path, so `"controlParameters.NCYCLES"` becomes
`"shelx.container.controlParameters.NCYCLES"`.

### Uploading plain files

```python
self.upload_file("inputData.XYZIN", gamma_model_pdb)
self.upload_file("inputData.SEQIN", str(pir_path))
```

### Uploading MTZ files with column selection

When an MTZ contains many columns, use `upload_file_with_columns` to select
specific ones.  This triggers server-side MTZ splitting:

```python
self.upload_file_with_columns(
    "inputData.F_SIGF", gamma_mtz,
    column_labels="/*/*/[FP,SIGFP]"
)
self.upload_file_with_columns(
    "inputData.FREERFLAG", mtz8xfm,
    column_labels="/*/*/[FREE]"
)
```

The `column_labels` syntax is `/*/*/[COL1,COL2,...]`.

### List parameters (CList items)

For list-type inputs (e.g., unmerged file lists, ensemble lists), first add an
empty list item, then upload/set on the indexed item:

```python
# Add an item to the UNMERGEDFILES list
self.add_list_item("inputData.UNMERGEDFILES")

# Upload a file to the first item
self.upload_file("inputData.UNMERGEDFILES[0].file", unmerged_mtz)

# For multiple items, repeat:
self.add_list_item("inputData.UNMERGEDFILES")
self.upload_file("inputData.UNMERGEDFILES[1].file", second_mtz)
```

You can also pass structured data directly:

```python
self.set_param("inputData.ASU_CONTENT", [{
    "sequence": "MKTAYIAKQ...",
    "nCopies": 1,
    "name": "my_protein",
    "polymerType": "PROTEIN"
}])
```

### Nested/compound parameters

For parameters under sub-containers, use dot notation:

```python
self.set_param("prosmartProtein.TOGGLE", True)
self.set_param("prosmartProtein.MODE", "SELECTED")
self.upload_file("prosmartProtein.REFERENCE_MODELS[0]", cif_path)
self.set_param("prosmartProtein.CHAINLIST_1", "A")
```

### Overriding the object path prefix

In rare cases, a parameter may not live under `{task_name}.container`.  Use
the `prefix` argument:

```python
self.set_param("some.other.path", value, prefix="custom.prefix")
```

---

## Running and waiting

### Standard pattern

```python
self.run_and_wait()  # Runs job and asserts success
```

This calls `self.run_job()` then `self.wait_for_completion()`.  If the job
does not reach a terminal status within `self.timeout` seconds, a
`TimeoutError` is raised.

### Expecting failure

```python
# For tests where the job should fail (e.g., wrong model, bad data)
self.run_and_wait(expect_success=False)
```

### Custom timeout

```python
self.run_and_wait(timeout=600)  # Override for slow pipelines
```

---

## Output validation

### Check output files exist

```python
self.assert_file_exists("XYZOUT.cif")
self.assert_file_exists("program.xml")
path = self.assert_file_exists("HKLOUT.mtz")  # Returns Path object
```

### Validate with gemmi (local)

```python
self.validate_pdb("XYZOUT.pdb")       # Reads with gemmi.read_pdb
self.validate_mtz("FPHIOUT.mtz")       # Reads with gemmi.read_mtz_file
st = self.validate_structure("XYZOUT.cif")  # Returns gemmi.Structure
```

### Validate via API digest (preferred for new tests)

The digest endpoint validates files server-side and returns crystallographic
metadata — no local gemmi needed:

```python
# Basic: assert the file is valid and registered in the database
digest = self.assert_valid_mtz_output("HKLOUT")
digest = self.assert_valid_coords_output("XYZOUT")

# Detailed: check specific properties
self.validate_mtz_via_api("HKLOUT",
    expected_spacegroup="P 21 21 21",
    expected_resolution_high=1.81,
    required_columns=["FP", "SIGFP"]
)

self.validate_coords_via_api("XYZOUT",
    expected_spacegroup="P 21 21 21",
    expected_chains=["A", "B"],
    min_residues=200
)
```

### Parse program.xml

```python
xml = self.read_program_xml()

# Refinement cycles
cycles = xml.findall(".//RefmacInProgress/Cycle")
rworks = [float(c.find("r_factor").text) for c in cycles]
assert rworks[-1] < rworks[0]  # R-factor improved
assert rworks[-1] < 0.27

# Phaser LLG
llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
assert max(llgs) > 1000
```

### Read output files

```python
log = self.read_output_file("log.txt")
foms = [float(x) for x in re.findall(r"FOM is (0\.\d+)", log)]
assert max(foms) > 0.7
```

### Access job directory directly

```python
job_dir = self.get_job_directory()
assert (job_dir / "job_1" / "RESTRAINTS.txt").exists()
```

---

## Mapping i2run parameters to API parameters

The main translation task when converting an i2run test to an API test is
mapping CLI arguments to `set_param` / `upload_file` calls.

### i2run CLI → API

| i2run CLI | API equivalent |
|-----------|----------------|
| `args += ["--XYZIN", path]` | `self.upload_file("inputData.XYZIN", path)` |
| `args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]` | `self.upload_file_with_columns("inputData.F_SIGF", mtz, column_labels="/*/*/[FP,SIGFP]")` |
| `args += ["--NCYCLES", "2"]` | `self.set_param("controlParameters.NCYCLES", 2)` |
| `args += ["--USEANOMALOUS", "True"]` | `self.set_param("controlParameters.USEANOMALOUS", True)` |
| `args += ["--ADD_WATERS", "False"]` | `self.set_param("controlParameters.ADD_WATERS", False)` |
| `args += ["--SELECTIONS", "text=A/"]` | `self.set_param("inputData.XYZIN.selection.text", "A")` |

### Where to find parameter paths

1. **Look at the task's `.def.xml`** — this defines the container structure.
   The CLI keyword maps to a leaf node; the API path is the full dot-separated
   path from the container root.

2. **Look at an existing i2run test** for the same task — the `--KEYWORD`
   flags correspond to leaf names.  Wrap file paths with `upload_file()`,
   scalars with `set_param()`.

3. **Look at the existing API tests** — `test_utilities_api.py` and
   `test_refinement_api.py` cover many common patterns.

4. **Check the frontend task interface** — the React component shows the
   object paths used for `setParameter` calls.

### Common container prefixes

Most tasks follow this structure:

| Prefix | Contents |
|--------|----------|
| `inputData.*` | Input files and data-related parameters |
| `controlParameters.*` | Algorithm options, flags, numbers |
| `monitor.*` | Monitoring and validation toggles |

Some tasks use custom prefixes (e.g., `prosmartProtein.*`, `options.*`,
`modelParameters.*`).  Check the `.def.xml` if in doubt.

---

## Disabling slow validation steps

Many pipelines include optional validation (MolProbity, Iris, Ramachandran,
B-factor analysis).  Disable these in tests to speed things up:

```python
# Common refinement validation toggles
self.set_param("controlParameters.VALIDATE_MOLPROBITY", False)
self.set_param("controlParameters.VALIDATE_IRIS", False)
self.set_param("controlParameters.VALIDATE_BAVERAGE", False)
self.set_param("controlParameters.VALIDATE_RAMACHANDRAN", False)
self.set_param("controlParameters.RUN_ADP_ANALYSIS", False)
self.set_param("monitor.RUN_COORDADPDEV_ANALYSIS", False)
```

---

## Complete worked example

Here is `test_refinement_api.py::TestProsmartRefmacAPI` annotated:

```python
import pytest
from .base import APITestBase

pytestmark = pytest.mark.pipeline


@pytest.mark.usefixtures("file_based_db")
class TestProsmartRefmacAPI(APITestBase):
    task_name = "prosmart_refmac"       # 1. Plugin name
    timeout = 300                        # 2. Generous timeout

    def test_gamma_basic(self, gamma_mtz, gamma_model_pdb):
        # 3. Create project and job
        self.create_project("test_refmac_gamma")
        self.create_job()

        # 4. Upload files (anomalous intensities + model)
        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        # 5. Set parameters
        self.set_param("controlParameters.NCYCLES", 4)
        self.set_param("controlParameters.USEANOMALOUS", True)
        self.set_param("controlParameters.VALIDATE_MOLPROBITY", False)

        # 6. Run and wait for completion
        self.run_and_wait()

        # 7. Validate output files
        self.validate_pdb("XYZOUT.pdb")
        self.validate_mtz("ABCDOUT.mtz")
        self.validate_mtz("FPHIOUT.mtz")

        # 8. Check refinement improvement via program.xml
        xml = self.read_program_xml()
        cycles = xml.findall(".//RefmacInProgress/Cycle")
        rworks = [float(c.find("r_factor").text) for c in cycles]
        assert len(rworks) == 5         # Initial + 4 cycles
        assert rworks[-1] < rworks[0]   # R-factor improved
        assert rworks[-1] < 0.27        # Reasonable final R-factor
```

---

## Running tests

```bash
cd server

# Single test file
./run_test.sh ccp4i2/tests/api/test_refinement_api.py -v

# Single test method
./run_test.sh ccp4i2/tests/api/test_refinement_api.py::TestProsmartRefmacAPI::test_gamma_basic -v

# All API tests
./run_test.sh ccp4i2/tests/api/ -v

# Parallel (recommended)
./run_test.sh ccp4i2/tests/api/ -n 4

# With stdout visible
./run_test.sh ccp4i2/tests/api/ -v -s
```

---

## Debugging failed tests

Failed tests preserve their project directory under `~/.cache/ccp4i2-tests/`.
The directory name includes `api_` prefix for easy identification.

Key files to inspect:

| File | Contents |
|------|----------|
| `project.sqlite` | Test database (can query with `sqlite3`) |
| `CCP4_JOBS/job_1/diagnostic.xml` | Error reports and warnings |
| `CCP4_JOBS/job_1/program.xml` | Task results and statistics |
| `CCP4_JOBS/job_1/input_params.xml` | Parameters as received by the task |
| `CCP4_JOBS/job_1/log.txt` | Full program log output |
| `CCP4_JOBS/job_1/job_N/` | Sub-job directories (for pipelines) |

```bash
# List preserved test directories
ls -la ~/.cache/ccp4i2-tests/*api*

# Query the test database
sqlite3 ~/.cache/ccp4i2-tests/20260317_*/project.sqlite \
  "SELECT id, number, task_name, status FROM ccp4x_job;"

# Clean up
rm -rf ~/.cache/ccp4i2-tests/*
```

---

## Checklist for a new API test

1. **Identify the task name** as it appears in the plugin registry
2. **Find a corresponding i2run test** (if one exists) — it shows which
   parameters and data files are needed
3. **Choose test data** — use existing conftest fixtures where possible;
   add new session-scoped download fixtures if needed
4. **Create the test file** — `test_{category}_api.py` or add to an existing
   category file (refinement, utilities, mr_pipelines, etc.)
5. **Write the class** inheriting `APITestBase` with `task_name` and `timeout`
6. **Translate i2run args** to `upload_file` / `upload_file_with_columns` /
   `set_param` calls
7. **Disable slow validation** (MolProbity, Iris, etc.) unless you specifically
   need to test validation
8. **Assert on outputs** — prefer API digest validation for new tests;
   use gemmi or program.xml parsing for detailed checks
9. **Run locally**: `./run_test.sh ccp4i2/tests/api/test_mytask_api.py -v`
10. **Run full suite** to check for regressions:
    `./run_test.sh ccp4i2/tests/api/ -v`
11. **Update** `docs/test-coverage-tracker.md` with the new test

---

## APITestBase method reference

### Project & Job

| Method | Description |
|--------|-------------|
| `create_project(name)` | Create a project via API |
| `create_job(task_name=None)` | Create a job (uses `self.task_name` by default) |

### Parameters & Files

| Method | Description |
|--------|-------------|
| `set_param(path, value)` | Set a scalar/dict/list parameter |
| `upload_file(path, file_path)` | Upload a file to a parameter |
| `upload_file_with_columns(path, file_path, column_labels)` | Upload MTZ with column selection |
| `add_list_item(path, item_data=None)` | Add an empty item to a CList parameter |

### Execution

| Method | Description |
|--------|-------------|
| `run_job()` | Start job execution |
| `wait_for_completion(timeout=None)` | Poll until terminal status |
| `run_and_wait(timeout=None, expect_success=True)` | Run + wait + assert success |

### Validation (local)

| Method | Description |
|--------|-------------|
| `assert_file_exists(filename)` | Assert file exists, return Path |
| `read_output_file(filename)` | Read file as string |
| `read_program_xml(filename="program.xml")` | Parse XML, return ElementTree |
| `validate_pdb(filename)` | Read with gemmi, return Structure |
| `validate_mtz(filename)` | Read with gemmi, return Mtz |
| `validate_structure(filename)` | Read with gemmi (auto-detect format) |
| `get_job_directory()` | Return Path to job directory |

### Validation (API digest)

| Method | Description |
|--------|-------------|
| `assert_valid_mtz_output(param_name)` | Validate MTZ via digest endpoint |
| `assert_valid_coords_output(param_name)` | Validate coordinates via digest endpoint |
| `validate_mtz_via_api(param_name, ...)` | Detailed MTZ validation (spacegroup, resolution, columns) |
| `validate_coords_via_api(param_name, ...)` | Detailed coordinate validation (spacegroup, chains, residues) |
| `validate_cell(digest, ...)` | Check unit cell dimensions |
| `get_output_digest(param_name)` | Get raw digest dict for an output |

### Other

| Method | Description |
|--------|-------------|
| `get_validation()` | Get pre-run validation status |
| `check_validation(allow_warnings=True)` | Check validation passes |
| `get_job_status()` | Get current job status dict |
| `get_job_files()` | List registered output files |
| `get_job_report()` | Get report data |
