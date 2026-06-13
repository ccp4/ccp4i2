# Preferences Handling: Design

**Status:** Accepted — supersedes the 2025-03-11 proposal (in git history).
**Date:** 2026-06-13
**Authors:** Martin Noble, Claude

This document replaces the earlier *Proposal for Harmonization*. That proposal
debated database-backed (A) vs file (B) vs hybrid (C) storage and leaned toward A.
We now have three concrete deployment perspectives in hand — **cloud/container**,
**desktop Electron app**, and the **`i2`/`i2run` CLI** — and together they collapse
that debate. The good ideas from the proposal are carried forward below.

## What the fuller overview settles

1. **Bootstrap settings cannot be database-backed.** Where the database is, where
   projects live, and where CCP4 is are needed *before* you can open the database —
   chicken-and-egg. So they must come from a **file** (or env). This rules out pure
   Option A for the settings that matter most.
2. **Cloud doesn't need the file.** Containers set everything via environment
   variables (including secrets from Key Vault), so as long as **env wins**, the file
   is inert in cloud. The proposal's main objection to a file ("needs a volume mount
   in cloud") dissolves. Verified against Materia: it sets `DATABASE_URL` /
   `CCP4I2_PROJECTS_DIR` in `entrypoint.sh` and ships no `preferences.json`.
3. **i2run and the GUI must agree, and both boot the same `settings.py`.** A file read
   during settings load makes the CLI and the control plane coherent for free.

**Decision: file-based JSON (the proposal's Option B2), with environment-variable
override and an optional site-defaults layer.** This is B2 with its two cons —
"no cloud story" and "no API/CLI" — neutralised by env-precedence and by adding
accessors over the same file.

## The store

| Layer | Location | Holds | Audience |
|-------|----------|-------|----------|
| **User preferences** | `~/.ccp4i2/preferences.json` (home overridable via `CCP4I2_HOME`) | bootstrap + functional preferences | GUI + CLI + control plane + worker |
| **Site defaults** (optional) | `CCP4I2_SITE_PREFERENCES` env, else `$CCP4I2_ROOT/local_setup/site_preferences.json` | institutional defaults | admins |
| **UI chrome** | Electron `electron-store` (`userData`) | window geometry, zoom, theme | GUI only |

`~/.ccp4i2` is the modern, lowercase home that already holds the SQLite database and
the project store — one directory for all per-user state. JSON so the Electron/JS
side reads and writes it as easily as Python. It is the descendant of classic (Qt)
CCP4i2's `~/.CCP4I2` XML preferences.

**Why the file (not the DB) for the canonical store:** the bootstrap subset is
file-resident by necessity (point 1 above); keeping the functional preferences in the
*same* file avoids two sources of truth, and is human-readable, portable, and
diff-able. The DB/API is a *projection* for the GUI (see Consumers), not the source.

## Resolution precedence (every setting)

```
environment variable  >  user preferences.json  >  site defaults  >  built-in default
```

- **Cloud:** env vars are set for everything → files never consulted → strict no-op.
- **Desktop:** the user file is the persistence layer; site defaults seed institutions.

## Three preference classes

1. **Bootstrap / environment** — `ccp4Dir`, `projectsDir`, `database`, `ccp4i2Root`.
   Consumed by `config/settings.py`. File-resident by necessity. **Implemented.**
2. **Functional** — consumed by wrappers/pipelines via a `PREFERENCES()` accessor:
   `SHELXDIR`, `DIALSDIR`, `BUSTERDIR`, `COOT_EXECUTABLE`, `CCP4MG_EXECUTABLE`,
   `EXEPATHLIST`, `PDB_REDO_TOKEN_ID`, `PDB_REDO_TOKEN_SECRET`,
   `RETAIN_DIAGNOSTIC_FILES`, … In cloud these are supplied as **env vars / secrets**;
   on desktop they live in the file.
3. **UI chrome** — `theme`, `zoomLevel`, `devMode` — stay in `electron-store`; the CLI
   has no interest in them.

## Schema (`preferences.json`)

```json
{
  "version": 1,
  "ccp4Dir":     "/path/to/ccp4",
  "projectsDir": "/path/to/projects",
  "database":    "sqlite:////abs/db.sqlite3",   // or postgres://…  (env DATABASE_URL)
  "ccp4i2Root":  "/path/to/.../ccp4i2",
  "userPreferences": {                          // functional class (class 2)
    "SHELXDIR": "/opt/shelx",
    "PDB_REDO_TOKEN_ID": "…",
    "RETAIN_DIAGNOSTIC_FILES": false
  }
}
```

Top-level keys are the bootstrap class (mapped to env vars in `settings.py`);
`userPreferences` is the functional bag the `PREFERENCES()` accessor exposes.

## Consumers — one store, several readers

| # | Consumer | Mechanism | Status |
|---|----------|-----------|--------|
| 1 | `config/settings.py` (bootstrap) | `preferences.resolve(key, env, default)` | **Done** |
| 2 | `PREFERENCES()` accessor (wrappers) | real object over the same file+precedence, replacing the `None`-returning `CPreferencesStub` | Planned |
| 3 | `i2 preferences get/set/list/reset/migrate` | reads/writes the file | Planned |
| 4 | `GET/PATCH /api/preferences/` (GUI settings panel) | projects the file over REST | Planned |
| 5 | Electron main (CCP4-dir / projects pickers) | writes bootstrap keys to the file | Planned (desktop step 3) |

Single store; the API and CLI are interfaces onto it, not separate truths.

## CLI (carried from the proposal)

```bash
i2 preferences list                      # all keys + effective values (+ source layer)
i2 preferences get SHELXDIR
i2 preferences set SHELXDIR /opt/shelx   # writes user preferences.json
i2 preferences reset SHELXDIR            # drop user value (fall back to site/default)
i2 preferences migrate                   # import from legacy ~/.CCP4I2 XML
```

Because storage is the file, these are thin read/modify/write operations (atomic
write already provided by `config/preferences.save_preferences`).

## Migration from legacy (carried from the proposal)

One-time import, on first launch or via `i2 preferences migrate`:
detect `~/.CCP4I2/configs/guipreferences.params.xml` → parse → map legacy keys to the
schema above → write `~/.ccp4i2/preferences.json`. Note the case change
(`.CCP4I2` → `.ccp4i2`); on case-insensitive filesystems they coincide.

## Implementation status

- **Done:** `ccp4i2/config/preferences.py` (load/save/atomic/`resolve`, `CCP4I2_HOME`
  override) + `config/settings.py` resolving bootstrap settings as
  `env > preferences.json > default`; unit + end-to-end tests passing on slim CPython
  and ccp4-python; cloud path verified unaffected.
- **Next:** site-defaults layer in `resolve()`; the `PREFERENCES()` accessor over the
  functional class (closes the wrapper gap); Electron write-side; `i2 preferences`
  CLI; `/api/preferences/`; legacy migration.

## Cloud invariance (new, and load-bearing)

The env-first precedence guarantees containers are unaffected: Materia/cloud sets
config (and secrets) via env, ships no `preferences.json`, and the file layer returns
empty. The file is exclusively the desktop persistence layer; functional secrets like
`PDB_REDO_TOKEN_*` come from env/Key Vault in cloud and from the user file on desktop.
