# Project store, database location, and per-project relocation

**Status:** design note / decision record · **Audience:** desktop + backend devs
(written to share with Paul re: the `django-edit-project` branch)

## TL;DR

- `CCP4I2_PROJECTS_DIR` should **stay**. It is the "preferred location for new
  projects", and it is already a *user preference* (via `preferences.json`), with
  the environment variable kept as an override. Removing it (as the
  `django-edit-project` branch began to) is unnecessary and has been superseded by
  the `preferences.json` layer that landed afterwards.
- **Database location and projects location are already decoupled in the backend**
  — separate preference keys, separate defaults. The only coupling is a
  *desktop-launcher policy* that derives the db path from the chosen projects dir.
  We should relax that so the db defaults to the user home (the legacy layout),
  independent of `projectsDir`.
- **The bootstrap trio (`ccp4Dir`, `projectsDir`, `database`) is frozen at server
  start.** None can be changed "hot". Their home is the pre-launch config screen;
  any in-app change must be an explicit **"save & restart backend"** action.
- **Per-project relocation is the one genuinely useful idea to lift** from
  `django-edit-project`: a `move_directory` action that moves *one* project's
  directory and rewrites its DB `directory` column. This is hot-safe because it
  never touches the bootstrap settings or the db's own location.

---

## Background: how config resolves today

There is a shared preferences layer, `server/ccp4i2/config/preferences.py`, that is
the single source of truth across the GUI, the `i2`/`i2run` CLI, and the Django
control plane. It lives at `~/.ccp4i2-django/preferences.json` and resolves every
setting with the precedence:

```
environment variable  >  preferences.json  >  built-in default
```

`settings.py` uses it for the relevant paths:

| Setting | Env override | `preferences.json` key | Default |
|---|---|---|---|
| CCP4 install | `CCP4Dir` (store) | `ccp4Dir` | discovered |
| Projects root | `CCP4I2_PROJECTS_DIR` | `projectsDir` | `~/.ccp4i2-django/CCP4X_PROJECTS` |
| Database | `DATABASE_URL` / `CCP4I2_DB_FILE` | `database` | `~/.ccp4i2-django/db.sqlite3` |

Two important consequences:

1. **`CCP4I2_PROJECTS_DIR` is already "a user preference, not just an env var."**
   You set the preferred projects location by writing `projectsDir` into
   `preferences.json` (the Electron directory picker does exactly this). The env
   var is the override; cloud deployments drive everything through env vars and are
   unaffected by the file. So the goal the `django-edit-project` branch was reaching
   for (make it a preference) is already met — and met more cleanly than by
   deleting the setting.

2. **DB location and projects location default to *separate* places.** With no
   overrides, the db is `~/.ccp4i2-django/db.sqlite3` and projects live under
   `~/.ccp4i2-django/CCP4X_PROJECTS/` — i.e. the backend already matches the
   **legacy CCP4i2 layout** (db in the user's config home, projects elsewhere).

## The coupling, and where it comes from

Despite the backend defaults above, on the **desktop** the db ends up *inside*
`projectsDir`. That coupling is injected at two Electron points, not in the backend:

1. The projects-dir picker derives the db path from the chosen folder:
   `client/main/ccp4i2-ipc.ts` —
   `updatePreferences({ projectsDir, database: sqliteUrl(join(projectsDir, "db.sqlite3")) })`.
2. `startDjangoServer` sets `CCP4I2_DB_FILE = join(CCP4I2_PROJECTS_DIR, "db.sqlite3")`
   as an env var, which then *overrides* the settings default.

So **uncoupling is a desktop-side change only** — no backend work is required. Stop
deriving `database`/`CCP4I2_DB_FILE` from `projectsDir`; let the db fall to the
`USER_DIR` default (or an independently-chosen `database` preference).

### Tradeoff (why it's a decision, not just a bug)

| | Coupled (current desktop) | Uncoupled (legacy default) |
|---|---|---|
| Point at a folder | Get its db **and** projects as one portable/syncable unit | db stays in user home; folder is just project storage |
| One central index | ✗ each store root has its own db | ✓ one stable db indexes projects anywhere |
| Move the store root | db travels with it | db unaffected |

**Recommendation:** default to the **uncoupled / legacy** layout (db in
`~/.ccp4i2-django`, independent of `projectsDir`), while *keeping* the `database`
preference so anyone who wants a portable self-contained folder can still place the
db inside a store root explicitly. The preference schema already separates the keys,
so this is a small, contained change.

## The hard constraint: bootstrap settings are not "hot"

`settings.py` executes **once**, when Django boots. `CCP4I2_PROJECTS_DIR`, the
`DATABASES` dict, `CCP4Dir` etc. are module-level values bound at that moment.
Reading `settings.CCP4I2_PROJECTS_DIR` per request returns the value bound at
startup — editing `preferences.json` while the server runs changes the *file*, not
the *live value*. Therefore:

- **Database** — cannot be repointed hot (it is an open connection *and* a
  startup-bound setting).
- **projectsDir** — also not hot: a running server keeps the value it booted with
  until relaunch.

This means the natural (and honest) home for all three bootstrap settings is the
**pre-launch config / "get ready" screen** — which is already where `projectsDir`
and `ccp4Dir` live. A selector for these in the *running* app would mislead
("I changed it but nothing happened").

### If we ever offer in-app changes: "save & restart backend"

Electron owns the uvicorn lifecycle (`start-uvicorn`, and it holds the process
handle), so the only correct mechanism is save-and-restart:

1. Stop the uvicorn/Django process.
2. (For an actual *move* of `db.sqlite3`) relocate the file on disk while nothing
   holds it open.
3. Write the new `database` / `projectsDir` preference.
4. Restart uvicorn — it re-reads settings and binds the new location.

Presented to the user as **"Apply & restart backend"**, never a live toggle.

## Per-project relocation — the piece worth lifting

Moving a **single project's** directory is a different, *hot-safe* operation: it
moves that project's folder on disk and rewrites the project row's `directory`
column. It never touches the bootstrap settings or the db's own location, so it is
safe while the server runs. This is the one genuinely useful idea in
`django-edit-project`:

- **Server:** a `move_directory` action on `ProjectViewSet` (`shutil.move` +
  update the `directory` field). Add guards it currently lacks: require an absolute
  target, refuse if the project has running jobs, handle
  cross-device/permission errors, and — when no target is given — default to the
  preferred `CCP4I2_PROJECTS_DIR / slug` (the "reset to preferred location" case,
  making move symmetric with create).
- **Frontend:** the dedicated edit-project page + directory-field + hook from the
  branch can be ported as an additive capability.

## What to leave behind from `django-edit-project`

- The removal of `CCP4I2_PROJECTS_DIR` from `settings.py` / the Electron store /
  `startDjangoServer` — superseded by the `preferences.json` layer; would break the
  ~20 call sites that still use the setting.
- The deletion of the tags UI (`tags-of-project.tsx`, still live on `django`) and
  its replacement with a non-functional `NewTagDialog` stub.

## Suggested sequence

1. **Decide + document** the db/projectsDir decoupling (this note).
2. **Uncouple** the db default from `projectsDir` (Electron-side change only).
3. **Land `move_directory`** + edit-project page as an additive feature on a fresh
   branch off `django`, leaving tags and `CCP4I2_PROJECTS_DIR` untouched.

## Pointers

- `server/ccp4i2/config/preferences.py` — shared preference resolver + schema
- `server/ccp4i2/config/settings.py` — `CCP4I2_PROJECTS_DIR`, `DATABASES`, `USER_DIR`
- `server/ccp4i2/api/serializers.py` — `ProjectSerializer` default-directory fill-in
- `client/main/ccp4i2-ipc.ts` / `client/main/ccp4i2-preferences.ts` — desktop
  preference read/write + directory picker + server launch
- `client/renderer/components/config-content.tsx` — the pre-launch config screen
  (where the projects-dir picker lives)
