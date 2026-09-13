# `import_merged`: architecture, as-is and target

**A design note for reworking the merged-reflection import task so its
crystallographic logic lives, once, in Python — where its owners can maintain
it.** Prepared for the Stuart McNicholas / Paul Bond / Phil Evans discussion.

## Why this task is special

`import_merged` is not a batch converter with a form in front of it. Its *real*
job is **intent-capture**: read a reflection file, diagnose what is in it, and
get the user to resolve the ambiguities — which mmCIF block, which columns,
intensities or amplitudes, mean or anomalous, where the cell and space group
come from, whether to cut resolution, and how to handle FreeR. Only once that
intent is resolved does the pipeline run.

Two consequences follow:

- **The `def.xml` is a contract between an *interactive resolver* and a *batch
  executor*.** The GUI resolves the file + the user's choices into concrete
  parameters (`MMCIF_SELECTED_BLOCK/_CONTENT/_COLUMNS`, `HKLIN_OBS_COLUMNS`,
  `HKLIN_OBS_CONTENT_FLAG`, cell, SG, resolution range, FreeR mode); the
  pipeline then consumes them non-interactively.
- **`i2run`/CLI cannot cover the bases.** It writes those resolved parameters
  *directly*, so it exercises the executor with the intent already decided and
  never touches the resolution step — which is the actual task. Testing and
  reasoning about this task means testing the resolver, not the CLI.

## As-is: one job, diagnosed four times

The Django port of the pipeline is faithful and *current* with the Qt CIF work
(`mmcifutils.py`/`mmcifconvert.py` are essentially the Qt code, kept up to
date). The problem is not the conversion; it is that **"diagnose a reflection
file → resolve the user's intent" is implemented in about four places, with
three engines and two vocabularies**:

| Where | What it does | Language |
|---|---|---|
| Pipeline (`import_merged.py`) | `bestcolumns()`/`columnthings()`; **`cmtzsplit` binary** (no-cut MTZ) *or* gemmi `ImportMTZ` (cut MTZ) *or* `ConvertCIF` (mmCIF) | Python |
| Digest endpoint (`digest.py`) | gemmi + `mmcifutils.CifBlockInfo`/`getColumnGroups`, serialised to JSON for the UI | Python |
| Standard import (`upload_param.py`) | `find_column_selections()` + `gemmi_split_mtz()` | Python |
| React interface (`import_merged.tsx`) | `groupColumnsByPattern()` / `COLUMN_PATTERNS` / `parseMtzColumns` | **TypeScript** |

Compounding factors:

- **Two content-flag vocabularies** — `CObsDataFile` (1 = I±, 2 = F±, 3 = Imean,
  4 = Fmean) vs `mmcifutils.CIFLabelSets` (1 = Imean, 2 = I±, 3 = Fmean,
  4 = F±). The 2026-08-24 "sfCIF defaulted to Fmean" fix only band-aided the
  *frontend* (`MMCIF_TYPE_PREFERENCE = [2,1,4,3]`); the underlying ordering is
  still divergent.
- **Format detection is filename-extension only**, and two signals can disagree:
  dispatch uses `self.fformat` (extension), validation uses the frontend-set
  `HKLIN_FORMAT`.
- **Two MTZ engines** selected only by whether a resolution cut was set:
  `cmtzsplit` binary vs gemmi `ImportMTZ` — different label mapping, different
  content-flag inference, for the same file.

### Intent-capture the Qt GUI does and the React UI drops

Because the crystallographic decisions were forked into TypeScript, several were
simply not carried over:

- **No unmerged-data warning + abort** — the whole reason "import *merged*" is a
  separate task. Qt blocks it with a modal; Django proceeds silently.
  (`getMerged()` is also a stub returning `True`.)
- **StarAniso handling entirely absent** — Qt detects it (mmCIF `_software.name`,
  or an `SA_flag` MTZ column) and forces `SKIP_FREER` / `COMPLETE=False` + swaps
  the advice. The digest never computes it; the UI never sees it. *Largest
  functional gap.*
- **`COMPLETE` never set** for an external FreeR; **resolution plumbing**
  (`CAN_CUT_RESOLUTION` / `MAXIMUM_RESOLUTION` / `RESOLUTION_RANGE_SET`) never
  wired; `MMCIF_SELECTED_DETAILS` never populated.
- **No conditional visibility** of FreeR options, so contradictory combinations
  are selectable.
- **Generic MTZ is driven by two competing mechanisms at once** — an in-page
  obs-group auto-selector *and* a modal that re-downloads and re-parses the file,
  both writing `HKLIN_OBS_*`.

## The governing principle: put the logic where its owners live

The crystallographic decision-making in this task — what counts as observations,
I vs F, mean vs anomalous, ASU/Friedel handling, FreeR validity, StarAniso,
merged vs unmerged, and the *recommended default* in each case — is the domain
of Phil Evans, Kathryn Cowtan, Paul Bond. Phil authored the foundational Python
selection logic (the mmCIF-data-selection work). That expertise can inspect,
correct and extend **Python**; it faces a real learning curve in **TypeScript**.

So the boundary is a maintainability contract, not just a technical one:

- **Crystallographic *decisions* → Python.** Owned by the crystallographers.
- **UX *presentation* → TypeScript.** Layout, dialog flow, rendering the options
  Python computed. Owned by the frontend developers.
- **Litmus test:** *could Phil change this rule without opening a `.tsx` file?*
  Anything that fails is domain logic in the wrong language and belongs
  server-side.

This is also the longer game: crystallographic knowledge kept legible *in Python*
is exactly what later makes it legible to an agent (the "expert knowledge"
substrate). The same boundary serves Phil now and an agent later.

## Target: one Python diagnosis, behind the round-trip

Reuse Python, do not re-express it in the browser. The desktop app always runs
the local Django server, so a round-trip is free; import is a deliberate,
one-file action, so the latency of a single digest call is imperceptible. And a
server round-trip lets the UI, the pipeline and the standard-import path call the
**same Python function** — shared *code*, the strongest form of one-source-of-
truth (stronger than gemmi-WASM, which shares only the *library* and keeps two
code paths).

1. **One `diagnose_reflection_file()` in Python** — returns a structured, fully
   *resolved* proposal: format; blocks/datasets with cell/SG/λ; column groups →
   *canonical* content type; FreeR presence + validity; StarAniso; merged vs
   unmerged; resolution range; **and the recommended default in each case**. Uses
   gemmi + `mmcifutils` (Phil's lineage), one canonical content-flag vocabulary.
2. **The digest endpoint serves it** — richer than today: not "here are the
   columns" but "here is the proposed resolution + every alternative + warnings."
3. **The React interface goes thin** — delete `groupColumnsByPattern`,
   `parseMtzColumns`, `COLUMN_PATTERNS` and the duplicate label tables. React
   renders the server's options and defaults, and writes the chosen `def.xml`
   params back. **No crystallographic logic in TypeScript.**
4. **The pipeline executes resolved intent** — remove `bestcolumns()` /
   `isIntensity()` run-time guessing; trust the params or call the *same*
   `diagnose_reflection_file()`. Guaranteed to agree with what the UI showed.
5. **`upload_param` calls the same function** — so "associate a merged MTZ
   directly" and "import it via import_merged" agree about what is in the file,
   differing only in downstream depth (import adds FreeR completion + aimless QC
   + resolution cut).
6. **One content-flag vocabulary; one gemmi *MTZ/mmCIF* engine** — collapse
   `cmtzsplit`-vs-`ImportMTZ` into a single gemmi split that always applies a
   possibly-null resolution range (the same cad→gemmi move already made for
   export). This removes the last gemmi-replaceable binary on the *MTZ* path; see
   the format-coverage section for SHELX/Scalepack, which is different.

### Format coverage: metadata-rich vs metadata-poor

The formats in scope split into two kinds, and the diagnosis must serve both.
**Critical constraint:** the diagnosis runs on the request/validation path, which
is the **CCP4-free slim boundary** — it may not shell out to any CCP4 binary. So
`f2mtz`, `combat` and `scalepack2mtz` (the classic Qt converters) are **not
available** and cannot be used, on the diagnosis path *or* anywhere the slim
server must run.

gemmi 0.7.5 (verified) natively reads **MTZ, mmCIF and XDS ASCII**, and has **no
SHELX or Scalepack reader**:

| Format | Metadata in the file | Reader (slim-safe) | Intent-capture is about… |
|---|---|---|---|
| **MTZ** | cell, SG, datasets, columns, types | gemmi | *disambiguating* what's in it |
| **mmCIF** (`.cif/.ent/.mmcif`) | cell, SG, λ, labels, FreeR status | gemmi | *disambiguating* what's in it |
| **XDS** (`XDS_ASCII.HKL`) | cell, SG, λ | gemmi (`XdsAscii`) | mostly *disambiguating* |
| **SHELX** (`.hkl`) | almost none — positional columns, no cell/SG, data-type undeclared | **new pure-Python/numpy reader** | *supplying* what's missing + *declaring* I-vs-F |
| **Scalepack** (`.sca`) | cell + SG in the merged-`.sca` header; data-type intensities | **new pure-Python/numpy reader** | reading the header + confirming |

Three design consequences:

- **SHELX/Scalepack need small pure-Python/numpy readers — not the CCP4
  binaries.** Both are simple fixed-format text; a reader parses the columns
  (and, for `.sca`, the cell/SG header) and builds a `gemmi.Mtz` via
  `Mtz.add_dataset`/`add_column`/`set_data` with the (parsed or user-supplied)
  cell/SG/data-type. This is the **same move already made in the slim server** for
  `freerflag`, `matthews`, `chltofom` — a gemmi/numpy port of a CCP4 binary — and
  it carries the same obligation: a **parity test** against `f2mtz`/`scalepack2mtz`
  output on sample files (they exist in a full CCP4 for test generation, just not
  in the slim runtime). This is the one genuinely *new* code in the plan, and the
  one part with real correctness risk.
- **The diagnosis must be completeness-aware, not just content-aware.**
  `diagnose_reflection_file()` returns, per format, both *what is embedded* and
  *what is absent and must be supplied* — for `.hkl`:
  `{cell: absent, spaceGroup: absent, dataType: undeclared}` (for `.sca`, cell/SG
  come from the header). The resolver then *requires* the user to provide
  `UNITCELL`/`SPACEGROUP`/wavelength and to *declare* the data-type before the job
  may run. The UI seed exists (the "Crystal Information" card renders for formats ∉
  {MTZ, mmCIF}), but it needs the **data-type declaration** and **hard
  validation** — you cannot proceed with a `.hkl` and no cell.
- **This is the reverse of the MTZ/mmCIF flow**, and it needs *more*
  crystallographic judgment (is a bare SHELX file intensities or amplitudes? is
  the supplied cell plausible?), which is exactly why it belongs in Python,
  behind the same one diagnosis surface — never in TypeScript.

### The round-trip budget is ~one call per file

The digest returns *all* blocks and *all* groups up front, so the interactive
choices — pick a block, pick a content type, toggle resolution — are **local
selections among already-delivered options**, not new round-trips. A small
"re-resolve" endpoint is only needed if some derived value genuinely requires
re-reading the file on a choice, which is rare.

### Why not client-side gemmi-WASM here

gemmi-WASM is bundled and *could* diagnose in the browser. But for this task:

| | Python round-trip | gemmi-WASM |
|---|---|---|
| Diagnosis implementations | **1** (Python) | 2 (Python batch + WASM UI) |
| Maintainable by the crystallographers | **yes** | no (logic in TS/WASM) |
| Agreement guarantee | shared *code* | shared *library* (can drift) |
| Latency on load | one digest call (imperceptible) | instant |
| Needs the local server | yes (always present in the desktop app) | no |

The only thing WASM buys — zero-latency, offline diagnosis — is low value for a
one-shot import, and it costs exactly what we are trying to avoid: crystallographic
logic outside Python. Reserve gemmi-WASM for continuously-interactive surfaces
(map/model manipulation), not this.

## Phased path

1. **Extract `diagnose_reflection_file()`** from the scattered Python
   (`digest.py`, `import_merged.py`, `upload_param.py`) into one module with one
   content-flag vocabulary; make `digest.py` a thin caller. *No behaviour change,
   but one authority.*
2. **Enrich the digest** to emit the resolved proposal + the missing signals
   (StarAniso, merged/unmerged, resolution range, per-dataset cell, defaults).
3. **Thin the React interface** — render the enriched digest; delete the TS
   diagnosis heuristics; restore the dropped intent-capture (unmerged abort,
   StarAniso, `COMPLETE`, resolution plumbing, conditional visibility).
4. **Make the pipeline an executor** — trust resolved params / call the shared
   function; collapse the two MTZ engines onto gemmi; retire `cmtzsplit` here.
5. **Point `upload_param` at the shared function** so both import doors agree.

## Sharp bugs to fix regardless

- `getMerged()` is a stub returning `True` (an unmerged MTZ is mis-reported).
- `importmtz()` swallows failure — `status` is overwritten to `SUCCEEDED`.
- `freerfile` is bound only inside `if not wdir.exists()` → `NameError` when
  `job_1/` already exists (both `columnthings()` and `convertmmcif()`).
- Dead fields: `HKLIN_FREER_COLUMN` and `outputData.HKLOUT` are never read on the
  MTZ path.
- `getFormat()` returns `CString` or `str` inconsistently (flagged in a code
  comment); guard uniformly.

---

*A design note, not a settled plan — a concrete shape to push against. Grounded
in the code as of this writing: the pipeline
(`server/ccp4i2/pipelines/import_merged/`), the digest
(`server/ccp4i2/lib/utils/files/digest.py`), the standard import
(`server/ccp4i2/lib/utils/files/upload_param.py`) and the React interface
(`client/renderer/components/task/task-interfaces/import_merged.tsx`).*
