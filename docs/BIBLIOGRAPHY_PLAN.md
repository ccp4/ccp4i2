# Bibliography — design & wiring plan

Status: **IMPLEMENTED** (2026-07-11, `django` branch). Util + endpoints + modal
wired; reference coverage raised from 33 to 128 of 173 tasks.

## What shipped

- **`lib/utils/jobs/bibliography.py`** — parses the MedLine `references/` files
  directly (no report-layer dependency); `references_for_tasks`,
  `task_names_for_job` (subjob burrowing + optional progenitors),
  `task_names_for_project`, `bibliography_for_job/project`.
- **`TASK_CITES` map** — the key design refinement. Two classes of task:
  - **Subjob-composing pipelines** (prosmart_refmac, servalcat_pipe,
    aimless_pipe, SubstituteLigand, import_merged…) whose steps are real i2
    plugins: **not** listed — `task_names_for_job` burrows `Job.children` and
    unions their tasks, so constituents surface from what actually ran.
  - **Monolithic pipelines** that drive programs *internally* (crank2,
    modelcraft, mrbump_basic, mrparse, slicendice, arcimboldo): these DON'T
    create i2 subjobs, so they enumerate their toolchain here (grounded in each
    pipeline's own source). Thin variants (phaser_MR_*, crank2_*, xia2_*…) alias
    to their canonical program key, so one physical `.medline.txt` per program.
- **Reference files** — imported the 11 bundled `crank2/.../references/*.nbib`
  (buccaneer, shelxc/d/e, bp3, afro, crunch2, prasa, solomon, multicomb, crank2)
  into `references/` under canonical keys, and authored ~40 new
  `{program}.medline.txt` for citable CCP4/crystallography programs that lacked
  them (ctruncate, pointless, aimless, freerflag, fft, phaser, molrep, mosflm,
  dials, xia2, coot, shelx, prosmart, sheetbend, nautilus, acedrg, chainsaw,
  sculptor, mrbump, mrparse, servalcat, metalCoord, nucleofind, matthews,
  pisa, pairef, dimple, zanuda, acorn, clustalw, AUSPEX, …).
- **API** — `GET /api/ccp4i2/jobs/{id}/bibliography/?progenitors=0|1` and
  `GET /api/ccp4i2/projects/{id}/bibliography/`, returning
  `{success, data:{result:[{pmid,title,authors[],source,link,cited_by}]}}`.
- **Frontend** — `components/bibliography-dialog.tsx` modal; `tool-bar.tsx`
  Bibliography button wired to open it for the current job. NB reads
  `data.data?.result ?? data.result` — `get_endpoint` does NOT unwrap the
  `api_success` envelope (same gotcha the Export MTZ button hit). The modal
  offers **Copy to clipboard**, **Download as text**, and **Download as BibTeX**
  (BibTeX entries generated client-side from author/title/source/year/DOI).
- **Tests** — `tests/unit/lib/test_bibliography.py` (parsing, cites-map, dedup,
  every-key-resolves, every-file-parses) and
  `tests/api/unit/test_bibliography_api.py` (job subjob-burrow, monolithic
  cites-map, project union, 404).

Original design notes retained below.

---

The job-panel toolbar (`client/renderer/components/tool-bar.tsx`) has a
**Bibliography** button that is currently a no-op (`onClick: () => {}`). This doc
captures the truth about how bibliographic references already flow into reports,
and the plan to expose them on demand (modally) at job and project scope.

## The truth: references already power task reports

Every task report already renders a bibliography via the report metadata layer:

```python
# server/ccp4i2/report/core.py (add-standard-sections)
from ccp4i2.report.metadata import ReferenceGroup
referenceGroup = ReferenceGroup()
referenceGroup.loadFromMedLine(self.TASKNAME)   # <- the mechanism
```

`ReferenceGroup.loadFromMedLine(taskName)`
(`server/ccp4i2/report/metadata.py:214`) reads:

```
I2_TOP / "references" / f"{taskName}.medline.txt"
```

and parses **MedLine** format (`PMID- / TI - / SO - / AU - / URL -`) into
`Reference` objects (`articleTitle`, `authorList`, `source`, `articleLink`).

- **34 `*.medline.txt` files** live in `server/ccp4i2/references/` — this is the
  live, working source of truth used by every report today.
- There is a **second, partly-redundant** source: per-plugin
  `*.bibtex.txt` / `*.bibtxt.txt` files in the wrapper/pipeline `script/` dirs
  (22 files, BibTeX-with-PMIDs). These are **not** what `loadFromMedLine` reads.
  For this feature we standardise on the **MedLine `references/` dir** (the one
  reports actually consume); the bibtex files are a legacy parallel we can
  reconcile later, not depend on.
- The `CBibReference*` container classes in `core/CCP4Annotation.py` are **empty
  stubs** on this branch — ignore them; the report `metadata.ReferenceGroup` is
  the real machinery.

## Scope semantics (per the toolbar context)

The button lives on the job panel but should also work project-scoped:

| Context | Bibliography returned |
|---------|-----------------------|
| **Job selected** | Union of references for the job's `task_name` **+ all subjobs** (`Job.children`, recursively). This is the honest bibliography of what actually ran. |
| **Job + progenitors** *(optional richer mode)* | Also include tasks of `FileUse`-linked input progenitors (walk the input side of the `FileUse` graph). Gives "everything that contributed to this result." |
| **Project scope (no job)** | Union over the `task_name` of **all jobs in the project** → compiled project bibliography. |

Model facts that make this expressible (`server/ccp4i2/db/models.py`):
- `Job.parent` / `Job.children` (subjobs, numbered `1.1`) — for subjob burrowing.
- `FileUse(file, job, role)` with `Role` (input/output) — for progenitor walking.
  `JobViewSet.dependent_jobs` already walks the *forward* (consumer) direction;
  progenitors are the reverse (input-role) walk.

## Plan

1. **`lib/utils/jobs/bibliography.py`** (new) — pure, report-free helper:
   - `references_for_tasks(task_names: set[str]) -> list[dict]`: call
     `ReferenceGroup.loadFromMedLine` per task, collect `Reference` objects,
     dedupe (by PMID / title), return plain dicts
     (`title`, `authors[]`, `source`, `link`).
   - `task_names_for_job(job, include_subjobs=True, include_progenitors=False)`.
   - `task_names_for_project(project)`.
2. **API** on `JobViewSet` / `ProjectViewSet`:
   - `GET /api/jobs/{id}/bibliography/?progenitors=0|1` → job (+subjobs
     [+progenitors]) references.
   - `GET /api/projects/{id}/bibliography/` → compiled project references.
   - Return structured JSON (not rendered HTML) so the modal owns presentation.
3. **Frontend** — `tool-bar.tsx`:
   - Wire the Bibliography button to open a **modal** listing references
     (title, authors, source, clickable DOI/PubMed link).
   - When a job is selected → job endpoint; when project-scoped → project
     endpoint. Optionally a toggle for "include progenitors".
   - Reuse an existing MUI `Dialog` pattern (cf. the `I2RunDialog` /
     export-success dialog) for consistency.

## Notes / open questions

- ~~**Missing `references/{task}.medline.txt`**: `loadFromMedLine` appends error
  code 100 and returns empty — treat as "no refs for this task", don't fail the
  whole request.~~ **Superseded 2026-08-24.** Returning empty is still right —
  the request must not fail — but doing so *silently* is what hid a real bug for
  months. It now also logs a warning naming the task, unless the task is
  declared in `NON_CITABLE`. See "How a citation resolves" below.
- Reconcile the **bibtex vs medline** duplication eventually (single source);
  out of scope for the first wiring.
- Dedup key: prefer PMID where present; fall back to normalised title.
- This shares its shape with Export MTZ: *call plugin/report metadata via a
  server util, expose through an endpoint, render in a frontend modal.*

## How a citation resolves

*Reference section, added 2026-08-24 after three related bugs. If you are adding
a task, this is the part you need.*

A task's references come from MedLine-format files in
`server/ccp4i2/references/`. Resolution has three steps, and **both** consumers —
the Bibliography button and the per-report bibliography — now go through all
three, from the same maps in `server/ccp4i2/core/citations.py`:

1. **Task name → citation keys.** `TASK_CITES` maps a task to the keys it cites.
   Unmapped tasks cite themselves. The map is one-to-many: `ShelxCD` →
   `shelxc` + `shelxd`; every `crank2_*` step → the whole toolchain it can drive.
2. **Key → file.** `CCP4Utils.findReferenceFile(key)` looks for
   `references/{key}.medline.txt`, **exact name first, then case-insensitively.**
3. **File → references.** A record needs a title *or* a source line. `SO` alone is
   not required, because a program is a legitimate citation and has no journal.

A task with nothing to cite — i2 plumbing, format shims, glue — belongs in
`NON_CITABLE`. A task whose citation simply lives under another key belongs in
`TASK_CITES`. Putting the latter in `NON_CITABLE` silences a warning by throwing
away a real citation; that is what the two lists are for, and they are not
interchangeable.

### Adding a citation for a new task

- If the program already has a `references/*.medline.txt`, add your task to
  `TASK_CITES` pointing at that key. Do not copy the file.
- If it does not, add `references/{key}.medline.txt` (MedLine: `PMID-`, `TI  -`,
  `SO  -`, `AU  -`, `URL -`), naming it for the program rather than the task, so
  other tasks can cite it too.
- If it genuinely has nothing to cite, add it to `NON_CITABLE`.
- Do nothing and the task's report warns at load, naming itself. `pytest
  ccp4i2/tests/unit/lib/test_reference_missing_is_loud.py` fails if any report
  class would render an empty bibliography.

### Why both maps live in `core/citations.py`

`lib/utils/jobs/bibliography.py` imports `ccp4i2.db.models` at module level;
`report/metadata.py` is Django-free and must stay that way, since reports render
outside a configured Django. Importing the maps from `bibliography` — even lazily
inside a method — drags Django into report rendering. `core/citations.py` has no
imports at all, so both sides can share it. A test renders a bibliography in a
subprocess and inspects `sys.modules` to keep that honest; the import-time check
alone does not catch a lazy import.

### Three bugs this shape has already caused

All three were one file being reached two ways, and all three were **silent**:

| Symptom | Cause |
|---|---|
| AceDRG citation missing from every report **on Linux only** | Report asked for `Acedrg.medline.txt` (its `TASKNAME`), the file is `acedrg.medline.txt` (the citation key). macOS/Windows resolve case-insensitively; Linux does not. |
| 52 of 143 report classes rendered an **empty** references section | The Bibliography button expanded `TASK_CITES`; the report side did not. |
| BUSTER cited nothing despite having a file | `ReferenceGroup` required an `SO` line; BUSTER's entry is a program, not a paper. |

The lesson worth carrying: **never test this class of bug with `Path.exists()`**.
It consults the filesystem, so on macOS and Windows it answers *yes* to a name
that fails on Linux — passing on exactly the machines where the bug is invisible.
Compare against the real directory listing instead.

## Implementation note

Per agreed plan, **Bibliography is implemented in a separate conversation** to
keep context parsimonious; this doc is the handoff.
