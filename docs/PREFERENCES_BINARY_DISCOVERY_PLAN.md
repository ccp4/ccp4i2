# Preferences & binary discovery — design & wiring plan

Status: **IMPLEMENTED** (2026-07-11, `django` branch, branch
`feat/preferences-binary-discovery`).

## What shipped

- **`config/program_discovery.py`** — resolution: explicit `{PROG}_EXECUTABLE`
  → `{SUITE}DIR`+exe → `exePaths` list → `PATH` → missing. Reports the source.
  Pure stdlib (slim-server safe; imported by both the server probe and the
  ccp4-python job env).
- **`CCP4PluginScript`** — exec site now resolves via `resolve_program` (was a
  bare `shutil.which`), with an actionable not-found log. Plus a
  `_checkProgramAvailable` **advisory warning** in `runTimeValidity` so a
  missing program surfaces at **job-config time** (in the run dialog), not just
  at execution. Pipelines with no single `TASKCOMMAND` (crank2) are skipped —
  and crank2's internal shelx already reads the same `SHELXDIR` via
  `PREFERENCES().SHELXDIR` → `user_preference`, so no change was needed there.
- **`config/preferences_migration.py`** — one-time, **desktop-only**,
  non-destructive, idempotent import of the legacy Qt GUI's program-location
  prefs (`~/.CCP4I2/configs/guipreferences.params.xml`: COOT_EXECUTABLE,
  SHELXDIR, …, and the EXEPATHLIST `CExePath` dirs → `exePaths`). Runs lazily on
  the first `program-preferences` GET. Guarded by `preferences.is_desktop()`
  (`CCP4I2_LOCAL_SESSION_TOKEN`).
- **API** — `GET /config/discover-programs/` (universal probe),
  `GET /config/program-preferences/` (values + `editable`),
  `PATCH /config/program-preferences/set/` (desktop-only; 409 in cloud). Plus
  `GET /tips/` (random tip HTML, image srcs rewritten) and
  `GET /tips/image/<name>` (traversal-guarded).
- **Frontend** — new `/ccp4i2/preferences` route (Edit→Preferences now points
  here, not the launch screen) with a **Program locations** panel: explicit
  fields + editable `exePaths` list + live found/not-found status chips. New
  Help → **Tip of the day** dialog.
- **Cloud posture** — discovery reads work in cloud via env vars
  (`user_preference` env>file>default); migration + file writes are desktop-only
  (409 write in cloud, `editable:false`). Tips ship in the wheel/image
  (`tipsOfTheDay/` is not MANIFEST-pruned) and serve through the authed API.
- **Tests** — `tests/unit/config/test_program_discovery.py` (11: resolution
  order, precedence, env-override, migration desktop/cloud/idempotent/
  non-destructive) + `tests/api/unit/test_config_and_tips_api.py` (6: discover,
  desktop roundtrip, cloud read-only 409, tips random/specific/traversal).

Original design notes retained below.

---

The Edit → Preferences menu item routes to the config page
(`/ccp4i2/config`, `components/config-content.tsx`), which currently exposes
only CCP4 dir, Projects dir, Python, backend, Theme, Developer mode. The legacy
`main` branch had a much richer preferences set — and, crucially, **program /
binary discovery controls that django lacks entirely**.

## The gap (verified)

**Django task-binary resolution is a bare PATH lookup with no override.**
`core/CCP4PluginScript.py:1697`:

```python
exe_path = shutil.which(self.TASKCOMMAND)   # PATH only; no preference override
```

If a program (coot, ccp4mg, shelx, dials, buster) is not on `PATH`, the task
fails and the user has **no supported way to point CCP4i2 at it**.

**`main` solved this** in `CPluginScript.getCommand` (core/CCP4PluginScript.py
~230-255), resolving per-program against preferences:

| Pref (`main`) | Type | Meaning |
|---------------|------|---------|
| `COOT_EXECUTABLE` | path | explicit coot executable |
| `CCP4MG_EXECUTABLE` | path | explicit ccp4mg executable |
| `SHELXDIR` | dir | dir holding shelx* binaries (`join(dir, exeName)`) |
| `DIALSDIR` / `BUSTERDIR` | dir | external-suite install dirs |
| `EXEPATHLIST` | `CExePathList` | general extra search dirs; `getExecutable(exeName)` |

Django `config/preferences.py` resolves only `ccp4Dir`, `ccp4i2Root`,
`database`, `projectsDir`, `version`, `userPreferences` — none of the above.
`COOT_EXECUTABLE` appears once, in a **docstring example only** (line 141), not
consumed anywhere.

## Plan

### 1. Backend — a resolution layer

- **`config/preferences.py`**: recognise new keys — `exePaths` (list, the
  django name for `EXEPATHLIST`), and optional explicit `COOT_EXECUTABLE`,
  `CCP4MG_EXECUTABLE`, `SHELXDIR`, `DIALSDIR`, `BUSTERDIR`. Same
  env > preferences.json > default precedence as every other pref.
- **`CCP4PluginScript`**: a `resolve_task_command(taskcommand)` helper used at
  line 1697 in place of the bare `shutil.which`. Resolution order:
  1. an explicit `{PROG}_EXECUTABLE` pref for known programs (coot, ccp4mg);
  2. a `{SUITE}DIR` pref joined with the exe name (shelx*, dials, buster);
  3. each dir in `exePaths` (the general `EXEPATHLIST` equivalent);
  4. `shutil.which(taskcommand)` (PATH) — unchanged default;
  5. fail with an **actionable** message naming the program and pointing at the
     Preferences → Program locations UI.
- Keep it slim-server-safe: pure path resolution, no CCP4 import; the jobs run
  in ccp4-python where this executes.

### 2. Backend — a discovery/probe endpoint

- `GET /api/config/discover_programs/` (or a config action): for a known list of
  programs, report `{name, resolved_path | null, source: pref|dir|PATH|missing}`.
  Powers the live "found at /…" indicator in the UI. Read-only; no execution.

### 3. Frontend — "Program locations" in the config page

- Extend `components/config-content.tsx` with a section:
  - editable `exePaths` list (add/remove dirs);
  - optional explicit fields for coot / ccp4mg / shelx / dials / buster;
  - per-program live status from the discover endpoint (green "found at …" /
    red "not found"), so users see immediately whether a path fixed discovery.
- Persist through the same config write path the page already uses.

### 4. Migration note

Legacy `main` stored these in `configs/guipreferences.params.xml`. Django uses
`~/.ccp4i2/preferences.json`. A one-time importer (read the old XML → seed the
JSON keys) is a nice-to-have, not required for the first cut.

## Out of scope / decided

- Pure external-link menu items (Tutorials, Quickstart, YouTube, Send error
  report, Task documentation/info) stay **dropped** — no bundled content, just
  outbound URLs; per-job Help already links task docs.
- The full `main` preferences set (toolbar-button visibility toggles, font/zoom,
  invalid-frame styling, etc.) is a **later** pass; this plan is scoped to the
  binary-discovery gap the owner flagged, which is the one with runtime impact.

## Companion quick win (separate, parallel): Tips of the Day

`server/ccp4i2/tipsOfTheDay/*.html` (35 files) already exist. Add a tiny
`GET /api/tips/random/` (or serve one at random) + a small Help → "Tip of the
day" dialog. Independent of the preferences work; low effort, real content.
