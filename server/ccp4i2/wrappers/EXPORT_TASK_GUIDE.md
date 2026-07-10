# Making a task's data files exportable (the "Export MTZ" button)

The job panel's **Export MTZ** button lets a user download a task's principal
reflection file. This guide is for task (wrapper/pipeline) authors who want their
task to offer a meaningful export. See `docs/EXPORT_MTZ_PLAN.md` for the wider
design.

## How the button decides what to offer

For the current job the frontend calls
`GET /api/jobs/{id}/export_job_file_menu/`. The server
(`lib/utils/jobs/export.py`) builds the menu in this order:

1. **Plugin contract** — if your plugin *module* defines `exportJobFileMenu`
   **and** `exportJobFile`, they are used.
2. **Generic fallback** — otherwise (or if the plugin offers nothing usable),
   the job's **tracked output MTZ `File`s** are offered directly.

Every declared mode is **validated**: the server calls `exportJobFile` for it and
only keeps the mode if it resolves to a file that exists on disk. So an empty or
unusable result simply means the button is hidden — never a broken download.

You therefore have two ways to make a task exportable. Prefer the first.

## Option A (recommended): declare a tracked output MTZ

If your task produces a principal reflection file, declare it as an `outputData`
`CMtzDataFile` so it is gleaned into the database as a `File`. The generic
fallback then offers it automatically — **no `exportJobFile` code needed** — and,
being a tracked output, it survives project export and is protected from purge.

1. **def.xml** — add to `<container id="outputData">`:

   ```xml
   <content id="COMPLETE_MTZ">
       <className>CMtzDataFile</className>
       <qualifiers>
           <saveToDb>True</saveToDb>
           <label>Complete reflection file</label>
       </qualifiers>
   </content>
   ```

2. **wrapper `processOutputFiles`** — point it at the unsplit file and set an
   annotation, *before* you split it into mini-MTZs:

   ```python
   self.container.outputData.COMPLETE_MTZ.setFullPath(
       os.path.join(self.getWorkDirectory(), "hklout.mtz"))
   self.container.outputData.COMPLETE_MTZ.annotation.set(
       "Complete unsplit reflection file")
   ```

   The gleaner persists any `outputData` `CDataFile` that is set and exists on
   disk (`db/async_db_handler.py: glean_job_files`).

3. **Keep it away from the purge scratch list.** `core/CPurgeProject.py`
   `SEARCHLIST` marks **`hklout.mtz` and `hklin.mtz` as category 1**
   ("scratch — always safe to delete", purged in *every* context). If your
   tracked file is literally `hklout.mtz`, it will be gleaned and then deleted
   by the next purge. Two ways to avoid that:

   - **Rename** the tracked copy to a non-scratch name (e.g. write/point at
     `complete.mtz`), which no default glob matches; **or**
   - **Protect it** by giving your plugin class a `PURGESEARCHLIST` entry with
     **category 0** (override/keep), which removes the matching default rule:

     ```python
     class my_task(CPluginScript):
         PURGESEARCHLIST = [["hklout.mtz", 0]]   # 0 = keep, overrides the default
     ```

   Category 0 only affects *your* task, so this does not change scratch handling
   for tasks that don't export their `hklout.mtz`.

## Option B: the imperative `exportJobFile` contract

Use this for **pipelines** where the exportable file lives in a subjob, or must
be **reconstructed** from several subjob outputs. Define two *module-level*
functions (not methods) in your plugin script:

```python
def exportJobFileMenu(jobId=None):
    # One row per exportable file: [mode, menu label, mime type].
    return [["complete_mtz", "MTZ file", "application/CCP4-mtz"]]

def exportJobFile(jobId=None, mode=None, fileInfo={}):
    # Return an absolute path to the file to serve (or None if unavailable).
    ...
```

`jobId` is the job's UUID. Inside these functions you may use the legacy DB
facade, which is backed by the Django ORM
(`PROJECTSMANAGER().db()` → `CCP4i2DjangoDbApi`):

- `getChildJobs(jobId, details=True)` → `[(number, uuid, task_name), ...]`
- `jobDirectory(jobId=uuid)` → the job's directory
- `getJobFilesInfo(jobId=uuid, jobParamName=...)`

### Locator pattern — return an existing unsplit file

Find the subjob that produced the file and return its path. Be robust: search
**all** matching subjobs rather than assuming a fixed position, and cover every
variant filename. Example (servalcat_pipe):

```python
def exportJobFile(jobId=None, mode=None, fileInfo={}):
    from ccp4i2.core import CCP4Modules
    db = CCP4Modules.PROJECTSMANAGER().db()
    if mode != "complete_mtz":
        return None
    servalcat_jobs = [cj for cj in db.getChildJobs(jobId=jobId, details=True)
                      if cj[2] == "servalcat"]
    for cj in reversed(servalcat_jobs):
        jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=cj[1], create=False)
        for name in ("refined.mtz", "refined_diffmap.mtz"):   # xtal / spa
            p = os.path.join(jobDir, name)
            if os.path.exists(p):
                return p
    return None
```

### Reconstruct pattern — combine subjob outputs

When there is no single unsplit file (import/scaling pipelines), combine the
relevant subjob outputs. **Do not shell out to `cad`** — the server is CCP4-free.
Use the gemmi helper:

```python
from ccp4i2.lib.utils.jobs.export import combine_mtz_files
# observed-data MTZ first so it wins on any column-label clash; FreeR after.
return str(combine_mtz_files([obsOut, freerflagOut], exportFile))
```

`combine_mtz_files` copies all data columns from each source into one
HKL-matched MTZ using gemmi (see `lib/utils/jobs/export.py`).

## Checklist

- [ ] If your task has a natural unsplit output → **Option A** (tracked output),
      and keep it off the category-1 purge list.
- [ ] If it's a pipeline needing subjob lookup / reconstruction → **Option B**.
- [ ] `exportJobFile` returns `None` (not raise) when the file is unavailable —
      the server validates and simply hides the button.
- [ ] Never call the `cad` binary; use `combine_mtz_files` for joins.
