# Preferences Handling: Proposal for Harmonization

**Status:** Proposal
**Date:** 2025-03-11
**Authors:** Martin Noble, Claude

## Problem

Legacy CCP4i2 managed 60+ preferences via a hierarchical XML system:
site defaults (`<install>/local_setup/`) merged with user overrides (`~/.CCP4I2/configs/`).
Many of these preferences are *functional* — wrapper code queries them at runtime
for paths to external programs (Coot, SHELX, BUSTER, DIALS), PDB-REDO credentials,
diagnostic flags, and more.

The current Electron/Django branch has only 6 preferences (CCP4Dir, projects dir,
zoom, theme, devMode, projectRoot), all stored in Electron's own config.json.
The Django backend provides a `_PreferencesStub` that returns `None` for every
attribute — meaning wrapper code that asks for `PREFERENCES().SHELXDIR` silently
gets nothing.

This gap blocks two things:
1. **Migration** — users cannot bring their settings from legacy to new.
2. **Functionality** — features that depend on configured paths or tokens are broken.

## Preferences That Matter

Not all 60+ legacy preferences need to survive. Many were Qt GUI concerns
(toolbar button visibility, font sizes, frame styles) that the new React UI
handles differently. The preferences that *do* matter fall into two groups:

**Functional (consumed by Python wrappers/pipelines):**
- `SHELXDIR`, `DIALSDIR`, `BUSTERDIR` — external program directories
- `COOT_EXECUTABLE`, `CCP4MG_EXECUTABLE` — external program paths
- `EXEPATHLIST` — additional binary search paths
- `PDB_REDO_TOKEN_ID`, `PDB_REDO_TOKEN_SECRET` — service credentials
- `RETAIN_DIAGNOSTIC_FILES` — runtime behavior flag

**UI (consumed by the frontend):**
- `theme`, `zoomLevel`, `devMode` — already handled by electron-store
- Any new UI preferences the React app needs

## Options

### Option A: Database-backed preferences (Django model + CLI + API)

Store functional preferences in the Django database. Provide access via:
- **REST API** (`GET/PATCH /api/preferences/`) for the frontend
- **CLI** (`i2 preferences list|get|set|reset|migrate`) for terminal use
- **Python accessor** (replace `_PreferencesStub` with a real object)

Resolution order: environment variable > user DB value > site defaults file > None.

Electron-only UI preferences (zoom, theme) stay in electron-store where they belong.

**Pros:**
- Single source of truth accessible from Python, CLI, API, and GUI
- Works in Docker/cloud deployment (no filesystem assumptions)
- `i2 preferences set SHELXDIR /opt/shelx` is simple and scriptable
- Site admins can seed defaults via fixture or `local_setup/site_preferences.json`
- Natural extension of the existing `i2` CLI and Django REST patterns

**Cons:**
- Preference values live in SQLite/PostgreSQL, not in a directly editable text file
- Users cannot simply open a file in a text editor to inspect or change preferences
- Adds a database migration

### Option B: File-based preferences (XML or JSON on disk)

Keep preferences in a structured file on disk, similar to legacy. Either:
- **B1: Retain the legacy XML format** (`guipreferences.params.xml`)
- **B2: Use JSON** (simpler parsing, aligns with Electron conventions)

The file would live at `~/.CCP4I2/configs/preferences.{xml,json}`.
Site defaults would remain in `<install>/local_setup/`.

Replace `_PreferencesStub` with a file-reading object that loads and merges
site + user files on startup.

**Pros:**
- Human-readable and directly editable with any text editor
- Familiar to legacy users — same location, similar format (if XML)
- No database dependency
- Simple to inspect, back up, copy between machines

**Cons:**
- No API access without additional plumbing (would need a read/write endpoint
  that proxies the file, with locking/conflict concerns)
- CLI management still desirable — would need a script that parses and rewrites the file
- Docker/cloud deployment needs a volume mount or alternative mechanism for the config file
- Two sources of truth if Electron preferences remain separate

### Option C: Hybrid — file as source of truth, cached in DB

Preferences stored on disk (XML or JSON), but loaded into the database on startup
and exposed via API. Edits through the API or CLI write back to the file.

**Pros:**
- Human-readable file *and* API/CLI access
- File remains the canonical, portable artifact

**Cons:**
- Complexity of keeping file and DB in sync
- Race conditions if file is edited externally while app is running
- More code to maintain for marginal benefit

## CLI Management (applies to all options)

Regardless of storage backend, CLI access to preferences is valuable.
The `i2` command already handles `projects`, `jobs`, `files`, etc. Adding
`preferences` as a resource is natural:

```bash
i2 preferences list                          # show all with current values
i2 preferences list --category paths         # filter by category
i2 preferences get SHELXDIR                  # show one preference
i2 preferences set SHELXDIR /opt/shelx       # set a value
i2 preferences reset SHELXDIR                # revert to default
i2 preferences migrate                       # import from legacy XML
```

For Option A this calls Django management commands. For Option B/C this
reads/writes the config file directly.

## Migration Path (applies to all options)

A one-time import from legacy files:

1. Detect `~/.CCP4I2/configs/guipreferences.params.xml`
2. Parse the XML (reuse existing `CContainer.loadDataFromXml`)
3. Map legacy keys to new keys (simple dictionary)
4. Write to the chosen store (DB or file)

This could run automatically on first launch or explicitly via
`i2 preferences migrate`.

## Site-Level Defaults (applies to all options)

Institutional deployments need a way for admins to set defaults.
A file in `<install>/local_setup/` (JSON or XML) serves as the base layer,
overridden by user-level settings. This preserves the existing admin workflow.

## Recommendation

**Option A** (database-backed) is the cleanest fit for the architecture we're
building — it aligns with how projects, jobs, and files are already managed,
works across desktop and cloud deployments, and the CLI provides the
"human interface" that compensates for preferences not being in a plain text file.

However, the readability concern is legitimate: being able to `cat` a config file
and see what's set, or copy it to another machine, is a real workflow. If this
is a strong requirement, **Option B2** (JSON file) is the pragmatic choice —
it's human-readable, simple to parse, and avoids the XML baggage.

**Option C** adds complexity without proportionate benefit and is not recommended.

## Questions for Discussion

1. How important is direct text-file editability of preferences in practice?
   (vs. using `i2 preferences set` or the GUI)
2. Do we need to support the legacy XML format specifically, or is JSON acceptable?
3. Should migration from legacy be automatic on first launch, or explicit?
4. Are there other legacy preferences (beyond the functional ones listed) that
   need to survive?
