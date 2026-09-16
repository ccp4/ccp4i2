# File subType integrity (issue #524)

A deep dive into how a data file's `subType` is (and isn't) set as it enters the
CCP4i2 database, why that matters, and how to close the gaps — both forward
(stop making badly-typed data) and retrospectively (repair what's there).

Motivating case: PR #523. servalcat's half-map inputs declared a permissive
`requiredSubType` that silently disabled subtype filtering, so autopopulation
offered a *trimmed map* (subType 1) for a *half-map* slot (subType 5), the same
file twice. Fixing the slot exposed the deeper truth — **a file whose subType is
wrong or unset is invisible to the very filtering that should have caught it**,
and half maps (and masks) are exactly the kinds that can't be recovered from the
file alone.

---

## 1. The subType model

`subType` is a small, static enum describing the *kind* of a data object. It is a
**different axis** from `contentFlag` (which describes the columns/representation
and *is* recoverable from file content — see §5).

Real-space maps (`CMapDataFile`, `CCP4XtalData.py:2384-2388`):

| subType | meaning | recoverable from content? |
|--------:|---------|---------------------------|
| 1 | normal density | — (the default/assumption) |
| 2 | difference (Fo–Fc) | weakly (mean≈0, symmetric ±) — no code does this |
| 3 | anomalous difference | provenance only |
| 4 | **mask** (mode-0 region map) | **yes, in principle** (CCP4 mode 0 + values in a tiny set) — *not implemented* |
| 5 | **half map** (one of an FSC pair) | **no — genuinely unknowable from one file** |

Map coefficients (`CMapCoeffsDataFile`) reuse 1/2/3. Observations
(`CObsDataFile`) and phases (`CPhsDataFile`) have their own provenance enums
(OBSERVED/DERIVED/REFERENCE; UNBIASED/BIASED) — orthogonal to F-vs-I, which is a
`contentFlag`.

The key property: **subType 5 (half map) is unknowable from the bytes** — a lone
half map is bit-compatible with a normal map; "half-ness" is a relational
property (one of a pair). So a *missing* half-map declaration is undetectable and
silently degrades to "normal". That makes it the highest-integrity-risk subtype.

---

## 2. How subType reaches the DB — the four routes

`File.sub_type` is an `IntegerField(null=True)` (`db/models.py:383`), so "unset"
is `NULL`, not 0. Both the gleaner and project import copy subType **verbatim,
with no content inference anywhere** (`async_db_handler.py:743`;
`import_i2xml.py:522-523`). Correctness is therefore decided entirely *upstream*
of persistence.

### Route 1a — dedicated `Import<Type>` tasks
A wrapper produces an output file; the gleaner persists whatever subType the
output CData object carries.

| Task | subType set? | Gap |
|------|--------------|-----|
| **ImportMap** | wrapper honours `MAP_SUBTYPE` (`ImportMap.py:28-30`) | **The React UI never renders `MAP_SUBTYPE`** (`ImportMap.tsx` shows only `MAPIN`). So every map imported through the desktop app is silently **subType 1** — half maps and masks mis-typed. Only i2run/CLI (`--MAP_SUBTYPE`) can set it. **Headline forward gap.** |
| **ImportMapCoeffs** | **inferred from MTZ column labels** (`infer_map_subtype`, `import_common.py:135-151`) | Correct by design — the model to emulate. |
| **ImportObs / ImportPhases** | class default 1 | Can't record non-default; 1 is the common case (low severity). |
| **ImportFreeR** | none → 0/unset | Benign — subType is meaningless for free-R. |
| **ImportCoordinate / ImportDictionary** | n/a | subType not applicable. |

### Route 1b — generic browse / fetch-from-URL (`upload_file_param`)
Every file parameter, when a user browses a file or fetches one (PDBe / RCSB /
AlphaFold / PDB-REDO), goes through `upload_param.py`. Neither browse nor fetch
carries any subtype intent (`utils.ts uploadFileParam`,
`fetch-file-for-param.tsx`) — correctness is entirely the **destination slot's
`requiredSubType`**:

```python
# upload_param.py:781-789  (added in PR #523)
ownSubType = int(param_object.subType)         # 0 if unset
subType = ownSubType if ownSubType > 0 \
          else _primary_required_subtype(param_object.get_qualifier("requiredSubType"))
```

- **Correct** when the file already carries a subType (a previous-job output), or
  the destination slot declares a meaningful `requiredSubType` (e.g. a half-map
  slot → 5, a mask slot → 4).
- **Falls back to 1** when the slot has **no `requiredSubType`** — so a mask or
  half map browsed/fetched into an un-annotated `CMapDataFile` slot is recorded
  as an ordinary map. #523 fixed servalcat's slots; every other bare map/obs/phase
  input slot still yields blanket 1.

### Route 2 — output gleaning
`glean_job_files` copies `CDataFile.subType` verbatim. That subType is correct
only if the wrapper set it — from the def.xml `<default><subType>` **or**
programmatically in the `.py`. A bare output that the wrapper never types is
gleaned as `NULL`, which downstream reads as "normal". See §3 for the audit.

### Route 3 — project import
`import_file` writes `sub_type` from the `filesubtype` XML attribute, or leaves
it `NULL` if absent (`import_i2xml.py:522-523`). **Legacy/Qt-era projects** whose
`DATABASE.db.xml` lacks `filesubtype` (or has 0) import verbatim with no
inference — the same blind spot as gleaning, on historical data we don't control.

### Route 4 — i2run
Same as Route 2 (runs the wrapper, gleans the output), plus it *can* pass
`--MAP_SUBTYPE` and other controls the UI omits — which is why CLI-imported maps
can be correctly typed while UI-imported ones are not.

---

## 3. Output-declaration audit (Route 2 exposure)

163 `outputData` entries across 5 file classes; **86 declare a subType, 77 bare**.
Excluding full-MTZ `CMtzDataFile` (24 — not a mini-MTZ, subType not meaningful),
the map/obs/mapcoeffs picture is **84/126 typed (~67%)**, but that is dominated
by well-curated map-coefficients masking almost-entirely-bare observations:

| class | declares subType | bare |
|-------|-----------------:|-----:|
| CMapCoeffs | 72 | 2 |
| **CObs** | 4 | **31** |
| CMapDataFile (real-space) | 8 | 9 |
| **CPhs** | 2 | **11** |

Many wrappers compensate by setting subType **programmatically** in `.py`
(servalcat, refmac, lorestr, phasertng, prosmart_refmac, cmapcoeff, splitMtz,
dm_multidomain masks, …), so "bare def.xml" overstates the runtime gap — but
**bare-and-not-backfilled** outputs default to NORMAL. The reference pattern to
copy is `refmac.def.xml` (FPHIOUT=1, DIFFPHIOUT=2, ANOMFPHIOUT=3) and
`molrep_map.def.xml` (the only fully-curated real-space set: trimmed=1, mask=4,
halfmap=5).

**Exposure ranking:** observations (31/35 bare) and phases (11/13 bare) are the
widest, but low-consequence (their subType is provenance, rarely filtered on).
**Real-space maps are the sharp end** — 9/17 bare, and that's where the
unknowable subtypes (mask, half-map) live. `servalcat`'s own `MAP_FO`/`MAP_FOFC`
and `ImportMap`'s `MAPOUT` are bare examples.

---

## 4. Inference taxonomy — what can we recover, and how

A ladder from most to least authoritative:

1. **Authoritative — declaration lookup (the retrospective lever).** subType is a
   static property of an *output parameter*, so for a File gleaned from a known
   task+param it is a constant in the def.xml. `File → Job.task_name →
   locate_def_xml → outputData <content id=job_param_name> → <subType>`. Verified
   feasible and low-risk (§6). Recovers *any* declared subtype, including
   half-map — the only route that can.
2. **Content-recoverable, already implemented.** Map-coefficient diff(2)/anom(3)
   from MTZ column labels (`infer_map_subtype`, `splitMtz` patterns); and, on a
   separate axis, `CObsDataFile` F-vs-I `contentFlag` via
   `CMiniMtzDataFile.miniMtzType()` reading columns with gemmi. These are solved.
3. **Content-recoverable in principle — masks — but we deliberately don't.** A
   mask is conventionally a CCP4 mode-0 map with values in {0,1}, so a heuristic
   *could* tag subType 4. We reject it: a content guess would second-guess the
   authoritative declaration, it scales badly (a full-DB scan touches every map),
   and it buys nothing once the declarations are curated (§6). **Inference must
   never override registered intent.**
4. **Unknowable — half maps (5).** Not recoverable from one file; must be
   declared at the producing/importing parameter. Only `molrep_map` does. Every
   other half-map-producing or half-map-accepting slot must declare it or the
   data is silently mis-typed with no way to detect it later.

**The governing principle.** subType comes from *declared intent*, never from a
guess about a file's bytes. The def.xml declaration (and, on inputs,
`requiredSubType`) is the single source of truth; every mechanism below either
reads that declaration or leaves the value alone — nothing overwrites an existing
subType, and nothing infers one from content.

---

## 5. A note on contentFlag (why it's *not* the problem)

F-vs-I for observations is a `contentFlag`, and it *is* content-recovered:
`miniMtzType()` reads the columns and matches each class's
`CONTENT_SIGNATURE_LIST` — "the answer comes from the content rather than from a
contentFlag someone may not have set." So the column/representation axis is
self-healing; it's the **provenance/kind `subType` axis** — especially the
non-inferable mask and half-map — that needs the work below.

---

## 6. Recommendations

The whole plan is **curate the declarations, then recover from them** — no
content inference at all.

**Curation scope.** Declare subType in the def.xml (or set it in the wrapper)
**only where the intended subType is not the recognisable default or sole subType
for that parameter.** A plain normal map (the fallback), or a param that can only
ever be one kind, stays bare — the default already types it correctly. The edits
land exactly where a bare declaration silently mis-types: the non-default
difference / anom / mask / half-map outputs (and the corresponding inputs). This
is a small, targeted set, and it turns the def.xml into the complete authoritative
source the backfill recovers from.

### Forward (stop generating badly-typed data)
- **F1 — Expose `MAP_SUBTYPE` in ImportMap's interface** (`ImportMap.tsx`). One
  control, closes the headline UI gap so a user importing a half map/mask can say
  so. *(Small, high value.)*
- **F2 — Declare subType on the non-default *output* slots** (difference, anom,
  mask, half-map) that are currently bare, per the scope above. Prioritise the
  unknowable ones (half-map, mask) since a missing declaration there is
  undetectable. Use `refmac`/`molrep_map` as the pattern. This is the curation
  that gives the backfill its authoritative data.
- **F3 — Annotate `requiredSubType` on the non-default *input* slots** that accept
  maps/obs/phases, generalising #523 beyond servalcat, so browse/fetch captures
  the right subtype.

*(No content detector. A mask/binary-map heuristic is deliberately excluded — it
would second-guess the declaration, scale badly across a large DB, and add
nothing once the declarations are curated.)*

### Retrospective (repair existing databases) — authoritative only
A **report-first management command** (`--dry-run` default) that, per File where
`sub_type` is NULL/0 (never touching a row that already has a subType):
- Resolve `(task_name, job_param_name)` to the historical task's def.xml subType
  and backfill it. Rules (all verified): strip the `[n]` suffix and read the CList
  `<subItem>` subType for list params; skip params that declare no subType; use
  the task that actually ran (do **not** follow `.successor`); skip tasks absent
  from `TASKS` (log them).
- Report every `(task, param, old→new)` and everything skipped, so a human can
  eyeball before `--apply`.

Because it reads only the (now-curated) declaration and only fills holes, it is
safe, reversible, and never overrides registered intent. Run the same pass **on
project import** so legacy zips are repaired as they land, not just existing DBs.

---

## 7. Priorities

| Item | Integrity risk closed | Effort | Notes |
|------|----------------------|--------|-------|
| F1 ImportMap UI | high (half-map/mask mis-typed on every UI import) | small | do first |
| F2 declare subType on non-default (half-map/mask/diff) outputs | high (unknowable, silent) | small–medium | the curation; audit the bare real-space outputs |
| Retrospective authoritative backfill | high (repairs history + legacy imports) | medium | report-first; recovers from the F2 curation |
| F3 `requiredSubType` on non-default input slots | medium (browse/fetch) | medium | generalises #523 |
| CObs/CPhs bare outputs | low (provenance, rarely filtered) | large | lowest priority — mostly leave bare |

The through-line: **half maps and masks are the real exposure** because they're
the subtypes you can't recover from the file — so the leverage is entirely in
getting the *declaration* right (F1/F2) and then recovering from it (the
backfill). No content is ever inspected; declared intent is the sole source of
truth.
