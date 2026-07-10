# Bibliography — design & wiring plan

Status: **design, not yet implemented** (2026-07-10, `django` branch).

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

- **Missing `references/{task}.medline.txt`**: `loadFromMedLine` appends error
  code 100 and returns empty — treat as "no refs for this task", don't fail the
  whole request.
- Reconcile the **bibtex vs medline** duplication eventually (single source);
  out of scope for the first wiring.
- Dedup key: prefer PMID where present; fall back to normalised title.
- This shares its shape with Export MTZ: *call plugin/report metadata via a
  server util, expose through an endpoint, render in a frontend modal.*

## Implementation note

Per agreed plan, **Bibliography is implemented in a separate conversation** to
keep context parsimonious; this doc is the handoff.
