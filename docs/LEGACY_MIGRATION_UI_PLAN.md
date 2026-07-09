# Legacy CCP4i2 → Django Migration: UI + Hardening Plan (a2)

**Goal:** Iron out wrinkles in, and provide a UI for, migrating a legacy (Qt-era)
CCP4i2 installation into the new Django-backed CCP4i2.

**Primary audience (this round):** a **desktop user migrating their own
`~/.CCP4I2` database** into their new install. The server is their own machine;
legacy project directories are local and already present. A cloud/admin variant
(upload + remap + real admin gating) is a later extension, not this round.

**Decisions taken:** desktop self-migration first; **add file-copying** to the
importer so a migration can produce a clean, self-contained new install; nested
projects are copied out to `CCP4I2_PROJECTS_DIR` and repointed (mixed-mode).

**Guiding principle — this whole thing is UI-driven.** The wizard is the front
door and the sole intended entry point for a normal user. Every decision the
migration makes — source, copy vs schlurp, the per-project plan, and how each
structural issue (§2.5) is resolved — must be **presented and confirmed in the
UI**, never assumed from a default or hidden behind a flag. Consequences that
shape the design below:
- The backend does **no irreversible action without a preceding validate the UI
  has shown the user.** Validate → dry-run → import is a hard sequence in the UI,
  not an optional nicety.
- The REST layer must return everything the UI needs to *render* a decision (the
  full report, the per-project `plan`, per-issue `{type, severity, detail,
  resolution}`), not just a pass/fail. No decision is UI-invisible.
- The `manage.py import_sqlite` / `validate_sqlite` commands stay as a
  power-user / scripting / test path, but they are **not** the primary UX and must
  not grow behaviour the UI can't drive. Anything the CLI can do, the UI can do.

---

## 1. What already exists (backend is done)

Records-only migration is fully implemented and exposed over REST — there is just
no frontend.

| Layer | File | What it does |
|-------|------|--------------|
| Core validator | [`import_sqlite.py`](../server/ccp4i2/db/import_sqlite.py) `SQLiteValidator` | Read-only. Checks project dirs, `CCP4_JOBS/job_N` dirs, `CCP4_IMPORTED_FILES`, import-source files exist on disk; referential integrity; data quality. Returns a report with a top-level `summary`. |
| Core importer | `import_sqlite.py` `SQLiteImporter` | 7-phase transactional import (lookups → projects → groups → jobs → files → history → comments). Supports `dry_run`, `remap_dirs`, `continue_on_error`, de-dup by dir/uuid. Returns `{dry_run, stats, errors, source}`. |
| CLI | `manage.py import_sqlite` / `validate_sqlite` | Command wrappers over the two classes. |
| REST (admin-only) | [`admin_views.py`](../server/ccp4i2/api/admin_views.py) | `admin/validate-sqlite/`, `admin/import-sqlite/`, `admin/import-status/` (+ `admin/import-legacy/` for dumpdata JSON). Each accepts an **uploaded** `sqlite_db` **or** a server-side `db_path`, plus `dry_run`, `remap_from`, `remap_to`. |

**REST response shapes the UI consumes**

- `validate-sqlite/` → the full validator report:
  `{source, counts{Projects,Jobs,Files,…}, projects{total,dir_exists,dir_missing[],dir_empty[]}, jobs{total,dir_exists,dir_missing_count,dir_missing[]}, files{total,exists,missing_count,missing[]}, imported_files{total,source_exists}, integrity{ok,issues[]}, data_quality{ok,issues[]}, summary{ok, projects_on_disk:"n/total", jobs_on_disk, files_on_disk, import_sources_on_disk, integrity_issues, data_quality_issues}}`
- `import-sqlite/` → `{dry_run, stats{projects,jobs,files,…}, errors[], source}` (HTTP 400 if `errors`).
- `import-status/` → `{projects{total,groups,…}, jobs{…}, files{…}}`.

---

## 2. The wrinkles to iron out

1. **Records-only — nothing copies files.** `SQLiteImporter._import_projects`
   sets `Project.directory` to the (remapped) *legacy* directory verbatim
   ([import_sqlite.py:736](../server/ccp4i2/db/import_sqlite.py#L736)). No
   `shutil`/copy anywhere. Result: the new install reads jobs/files straight out
   of the old `~/.CCP4I2` tree. Fine short-term, but the old and new installs
   then **share the same directories** — the user can't safely delete the old
   install, and any legacy dir that's missing yields dangling DB rows with no
   error at import time.
2. **Validation is optional and easy to skip.** The validator detects (1), but
   nothing makes the user look before importing. The UI must surface the disk
   report *before* the real import.
3. **Concrete bug:** [`import_i2xml.py`](../server/ccp4i2/db/management/commands/import_i2xml.py)
   reads `options["program_xml"]` but the arg is `project_xml` — crashes if run.
   (Not on the sqlite path, but it's in the migration family; trivial fix.)
4. **Admin gating unenforced in the frontend**, and the renderer's `useAuth()` is
   a mock-admin stub. For desktop self-migration this is effectively moot — flag
   it, don't block on it.

---

## 2.5 Structural pre-flight checks (edge cases before schlurp/copy)

Terminology used here:
- **schlurp** = in-place migration — the new install adopts the legacy directory
  as-is (`copy_files=false`).
- **copy** = duplicate the legacy tree into `CCP4I2_PROJECTS_DIR` and repoint
  (`copy_files=true`).

Before *either* operation, the validator must decide whether each project's
directory is in a shape that's safe to schlurp/copy. Today it only checks
existence, integrity and data quality — it does **no structural reasoning about
the directory tree itself**. The key facts that make this necessary:

- `Project.directory` is `unique=True` in the model, so *exact-duplicate* dirs are
  rejected at DB level — but **nesting is not**.
- A project's real content lives under `<directory>/CCP4_JOBS` and
  `<directory>/CCP4_IMPORTED_FILES`. So if one project's directory sits inside
  another's tree, the outer project physically *contains* the inner one.

New `SQLiteValidator` section `structure` (each returns a list of offenders so the
UI can show exactly which projects are affected):

1. **Nested projects (the headline case).** For every ordered pair (A, B) of
   legacy project directories, flag when `B` is inside `A`
   (`Path(B).resolve()` is relative to `Path(A).resolve()`, using
   `is_relative_to` / `commonpath`, after remap).
   - Report `{outer, inner, relationship: "B under A/CCP4_JOBS/…"}`.
   - **Resolution (decided): copy the nested project out.** Any project involved
     in a nesting relationship is **always copied** to a fresh dir under
     `CCP4I2_PROJECTS_DIR` and repointed — *regardless of the overall
     schlurp/copy choice*. There is already a project root that is the default
     home for newly generated projects, so an un-nested clone belongs there. This
     turns the dangerous case into the safe one:
     - **copy mode:** already copying; just ensure the outer project's copy
       **excludes** the inner subtree (so no duplication) and the inner project
       gets its own dest. Both dests end up self-contained.
     - **schlurp mode:** the run is otherwise in-place, but nested projects are
       promoted to copy-out (**mixed mode**). Only the entangled projects move;
       everything else is still adopted in place. The UI must state plainly which
       projects will be copied and where.
   - Which end gets copied out? Default: copy the **inner** project(s) out and
     leave the outer in place (fewer/smaller trees to move, and the outer is
     usually the "real" project). If the *outer* is the smaller/less-important one
     that's a rare inversion — not worth a toggle in a2; document the rule.

2. **Directory collisions at the destination (copy mode only).** Two legacy
   projects whose **basename** is identical map to the same
   `CCP4I2_PROJECTS_DIR/<name>` even though their source dirs differ. Also flag
   when a computed destination already exists on disk (pre-existing project or a
   prior partial migration). Report intended dest per project + collisions.

3. **Shared / overlapping directories not caught by `unique`.** Symlinks or
   `..`-relative paths that resolve to the same real path as another project
   (canonicalise with `Path.resolve()` before comparing). `unique=True` compares
   the raw string, so `~/x` and `/home/u/x` and a symlink all slip past it.

4. **Directory is not actually a CCP4i2 project root.** Legacy dir exists but has
   no `CCP4_JOBS` subdir, or `CCP4_JOBS` is empty while the DB lists jobs →
   mismatch between DB and disk. Cheap sanity check that catches a wrong/moved
   remap target.

5. **Legacy directory outside expected roots / non-writable dest.** For copy:
   `dest_root` not writable, or insufficient free space vs. estimated tree size
   (sum sizes during the walk we already do for nesting). For schlurp: legacy dir
   under a path that won't survive (e.g. a temp mount) — advisory only.

6. **`I1ProjectDirectory` overlap.** Each project also carries a legacy-i1
   directory; if present it can nest/collide too. Lower priority — check for
   containment against project dirs, warn only.

7. **Case-insensitive filesystem name clashes (copy mode).** On macOS (default)
   and Windows, dest names `Foo` and `foo` collide even though they're distinct on
   the source Linux tree. Compare intended dest basenames **case-folded** (and
   NFC-normalised, see §8) when detecting §2 collisions — a naive
   case-sensitive check misses this and the second `copytree` silently merges into
   the first. Also flag two *source* project dirs that differ only by case when
   the migration host is case-insensitive.

8. **Unicode-normalisation clashes.** Two names equal under NFC/NFD (common when a
   tree has been round-tripped through macOS, which stores NFD) collide at the
   dest. Normalise to NFC before dest-name comparison; report as a §2 collision
   variant.

9. **Shared / symlinked `CCP4_IMPORTED_FILES` (or `CCP4_JOBS`).** Two projects may
   point their import/job subdirs at a common location via symlink even when the
   project roots differ (so §3's root-level canonicalisation misses it). Walk each
   project's `CCP4_IMPORTED_FILES` and `CCP4_JOBS`, resolve symlinks, and flag when
   subtrees of *different* projects resolve into a shared real directory. In copy
   mode this means files get duplicated or, worse, a `copytree` follows a symlink
   into another project. Default `copytree(symlinks=True)` (preserve links) plus a
   report of cross-project shared targets; blocking only if a shared target sits
   under another migrating project's dest.

10. **Read-only / permission-denied legacy trees (copy mode).** A legacy dir the
    server process can't read (or individual unreadable files within) makes
    `copytree` fail partway. Pre-flight an access check (walk + `os.access(...,
    R_OK)` on a sample, or catch during the size walk) and report unreadable
    paths per project so the copy doesn't die mid-import. For schlurp this is only
    a problem if the running install also can't read them at runtime — advisory.

Severity model:
- **Blocking in copy mode** (until acknowledged / auto-resolved): dest collisions
  (2), case/NFC clashes (7, 8), shared targets landing under a migrating dest (9),
  unreadable source trees (10).
- **Nesting (1) is auto-resolved** (copy-out to the project root) rather than
  blocking — but the UI must surface the resulting copy plan and let the user
  confirm before Step 4.
- Everything else — (3)–(6), non-blocking (9) — are **warnings**.
- In schlurp mode, nesting still triggers the mixed-mode copy-out (1); the other
  copy-only checks (2, 7, 8, 10) don't apply to projects staying in place.

All feed the same `structure.issues[]` the UI renders in Step 2, each carrying
`{type, severity, projects[], detail, resolution}` so the UI can render both the
problem and what the migration will do about it (e.g. "will be copied to
`CCP4X_PROJECTS/Foo_2`").

**Where this runs:** add `_validate_structure(cur)` to `SQLiteValidator`,
included in `run()`'s report and folded into `summary` (`structure_issues` count,
and `summary.ok` must go false on any *blocking* structural issue). The validator
computes, and returns, a **per-project migration plan**:
`plan[legacy_project_id] = {mode: "in_place" | "copy", dest, reason}`. The importer
consumes the same plan so validate and import cannot disagree — the importer must
not re-derive it independently. A project's plan `mode` is `copy` if the run is in
copy mode **or** the project is nesting-entangled (§2.5.1 mixed-mode), else
`in_place`. `dest` is a de-duplicated, case/NFC-safe name under `dest_root`.

---

## 3. Backend work

### 3a. Add `copy_files` to `SQLiteImporter` (the file-copying decision)

Give the importer an optional copy-in mode so a desktop migration can produce a
**self-contained** new install rather than sharing the legacy tree.

- New constructor args: `copy_files: bool = False`, `dest_root: Path | None`
  (defaults to `settings.CCP4I2_PROJECTS_DIR`, the same root native projects use —
  [settings.py:241](../server/ccp4i2/config/settings.py#L241)).
- **Plan-driven, not flag-driven.** The importer consumes the validator's
  per-project `plan` (§2.5 "Where this runs"): for each project it does what the
  plan says (`in_place` → adopt legacy dir; `copy` → `shutil.copytree` to
  `plan.dest`, repoint `Project.directory`). This is what lets a single run mix
  in-place and copy (nesting mixed-mode) without special-casing in the copy loop.
- `copy_files=true` simply means "every project's plan defaults to `copy`";
  nesting can flip an otherwise-`in_place` project to `copy` (§2.5.1).
- `copytree` uses `symlinks=True` (preserve links, §2.5.9) and an `ignore` that
  **excludes any other project's dir nested under this one** (§2.5.1) so nested
  content is never duplicated. Skip if source missing — record a skip stat.
- Copy happens **inside** the existing `transaction.atomic()` phase ordering but
  is a filesystem side-effect: on `dry_run`, copy nothing and just report intended
  destinations + source-missing warnings. On real run, copy before/at project
  creation so `directory` is correct.
- Idempotency: if `dest` already exists, skip-copy + reuse (mirrors the existing
  dir-based de-dup). Record `projects_copied` / `projects_copy_skipped` /
  `copy_errors` in `stats`.
- Surface `copy_files`, `remap`, and dest root in the returned `summary`/`stats`
  so the UI can show "N projects copied to <root>".

Plumb the flag through `admin_views.import_sqlite` (`copy_files` form field) and
the `manage.py import_sqlite` command (`--copy-files`, `--dest-root`).

**Guardrail for (1):** when **not** copying and a project's legacy dir is missing
(and no working remap), the importer should count it and the API should return it
prominently so the UI can warn rather than silently create dangling rows.

**Blocking-issue gate.** If the plan contains *blocking* structural issues that
are **not** auto-resolved (dest collisions, case/NFC clashes, shared targets,
unreadable trees — §2.5 severity model), `import-sqlite/` returns 400 with the
offenders and does nothing, unless the request carries an explicit
`allow_structural_issues=true` (set by the UI only after the user acknowledges the
Step-2 warning). Nesting is **not** in this gate — it's auto-resolved via the plan
(copied out to `dest_root`), so it never blocks; it only needs the user to confirm
the copy plan.

**Stats to record:** `projects_copied`, `projects_in_place`, `projects_nested`
(copied-out due to nesting), `copy_excluded_subtrees`, `projects_copy_skipped`,
`copy_errors`. Idempotency: if a computed `dest` already exists, skip-copy + reuse
(mirrors the existing dir-based de-dup). Surface the plan summary in the response
so the UI can show "N in place, M copied to `<root>` (K un-nested)".

### 3b. Fix `import_i2xml` command

`options["program_xml"]` → `options["project_xml"]`. One line.

### 3c. (Deferred) De-dupe importers

`import_legacy_ccp4i2` (dumpdata JSON) reimplements the same 7 phases as
`SQLiteImporter`. Out of scope for a2 (refactor risk, touches tests) — noted so it
isn't forgotten.

---

## 4. Frontend work

### Route & scaffold
New admin route, following house style
([import-project/page.tsx](../client/renderer/app/ccp4i2/(authed)/import-project/page.tsx)):

- `client/renderer/app/ccp4i2/(authed)/admin/migrate-legacy/page.tsx` — thin
  `"use client"` shell: `CCP4i2TopBar title="Migrate legacy CCP4i2"` + `Paper` +
  `<MigrateLegacyContent/>`.
- `client/renderer/components/migrate-legacy-content.tsx` — the wizard logic.

### The wizard (MUI `Stepper`, introduced fresh — no house pattern exists)

**Step 1 — Source & mode**
- Default: server-side **`db_path`**, prefilled `~/.CCP4I2/db/database.sqlite`
  (the desktop-self case). Secondary option: upload a `.sqlite` (for a DB that
  lives elsewhere).
- **Copy mode toggle** (the key UX for the file-copying decision):
  - **Copy projects into new install** (default, recommended) → `copy_files=true`,
    dest = `CCP4I2_PROJECTS_DIR`. "Your old install stays untouched; you can
    delete it afterward."
  - **Migrate in place** → `copy_files=false`. "Faster, no extra disk, but the
    new install reads from your old `~/.CCP4I2` directories — don't delete them."
- Optional `remap_from`/`remap_to` (advanced; for a DB whose recorded paths no
  longer match, e.g. moved home dir).

**Step 2 — Validate** → POST `admin/validate-sqlite/`.
- Render `summary`: `projects_on_disk`, `jobs_on_disk`, `files_on_disk`,
  `import_sources_on_disk` as `n/total` chips; integrity + data-quality issue
  counts.
- If anything is missing on disk → amber/red banner listing `projects.dir_missing`
  / `jobs.dir_missing` (this is wrinkle (1)/(2) made visible). In *copy* mode a
  missing source dir means that project can't be copied — say so here.
- **Structural section (§2.5)** — render `structure.issues` grouped by type:
  - **Nested projects** table: outer → inner, with the relationship path. This is
    the headline edge case; in copy mode it's blocking. Offer an explicit
    "I understand, proceed anyway" checkbox that sets `allow_structural_issues`
    for Steps 3–4; leave it off by default.
  - **Destination collisions** (copy mode): projects mapping to the same
    `CCP4I2_PROJECTS_DIR/<name>` or an existing dest.
  - Shared/canonicalised-duplicate dirs, non-project-root dirs, i1 overlaps —
    as warnings.
  - When `summary.ok` is false on a blocking issue, disable **Next → Dry run**
    until acknowledged (mirrors the run-dialog SEVERITY_ERROR blocking pattern).
- Reuse the status-table pattern from
  [import-project-content.tsx:119](../client/renderer/components/import-project-content.tsx#L119).

**Step 3 — Dry run** → POST `admin/import-sqlite/` with `dry_run=true` (+
`copy_files`, remap). Show per-model `stats` (rolled back) and any `errors`. In
copy mode, show intended destination root and projects-to-copy count.

**Step 4 — Import** → same endpoint, `dry_run=false`.
- On success, GET `admin/import-status/` and show resulting DB counts as
  confirmation; `usePopcorn().setMessage("Migration complete","success")`.
- On `errors` (HTTP 400), `usePopcorn().setError(...)` and keep the wizard on this
  step so the user can retry.

### API calls
`useApi().post("/admin/import-sqlite/", formData)` etc.
([api.ts](../client/renderer/api.ts); multipart via `FormData`, same as
`import-project-content.tsx`). Form fields:
`db_path`/`sqlite_db`, `remap_from`, `remap_to`, `dry_run`, `copy_files`,
`dest_root` (optional override), and `allow_structural_issues` (only sent after
the user ticks the Step-2 acknowledgement).

### Nav & gating
- Add an entry button to
  [projects-toolbar.tsx](../client/renderer/components/projects-toolbar.tsx)
  (`router.push("/ccp4i2/admin/migrate-legacy")`), styled like New/Import.
- Wrap the button and gate the page on `useAuth().canAdminister`
  ([auth-context.tsx](../client/renderer/lib/compounds/auth-context.tsx)). Note
  the renderer stub returns mock-admin, so this is a no-op on desktop and only
  bites in Docker — acceptable for the desktop-first target; verify gating in the
  Docker build before relying on it.

---

## 5. Tests

- **Backend (fast, no CCP4):** unit-test `SQLiteImporter(copy_files=True)` against
  a synthesised legacy sqlite + a temp legacy project tree; assert files land in
  `dest_root/<name>`, `Project.directory` is repointed, `stats` reports copies,
  and dry-run copies nothing. Test the missing-source-dir warning path. Extend
  `server/ccp4i2/tests/db/` (this is DB-layer, not CCP4-binary).
- **Structural checks (§2.5), fast, no CCP4** — synthesise a legacy sqlite whose
  Projects table exercises each case: (a) B nested under A's `CCP4_JOBS`; (b) same
  basename → dest collision; (c) symlink/`..` resolving to another project's dir;
  (d) dir with no `CCP4_JOBS`; (e) case-only difference `Foo`/`foo`; (f) NFC/NFD
  name pair; (g) two projects whose `CCP4_IMPORTED_FILES` symlink to a shared real
  dir; (h) an unreadable source dir. Assert `_validate_structure` flags each with
  the right `type`/`severity`, `summary.ok` reflects blocking vs warning, and the
  returned `plan` has the expected `mode`/`dest`/`reason` per project.
- **Plan/mixed-mode behaviour** — assert nesting produces a `copy` plan entry even
  when `copy_files=false` (only the entangled project moves; siblings stay
  `in_place`); assert copy-with-nesting excludes the inner subtree (no duplicated
  files, both dests self-contained); assert `allow_structural_issues=false` makes
  `import-sqlite/` refuse on a *blocking* issue but proceed on nesting alone;
  assert dest-name de-dup is case/NFC-safe on a case-insensitive tmpdir.
- **Backend regression:** the `import_i2xml` arg-name fix — smoke test the command
  parses.
- **API:** `admin/import-sqlite/` with `copy_files` + `db_path` (unit, Django test
  client, admin user) → asserts response `stats`/`summary` shape the UI relies on.
- **Frontend:** component test of `migrate-legacy-content` stepping through
  validate → dry-run → import with a mocked `useApi` (mirroring existing component
  tests, if any; otherwise manual `/run` verification of the wizard).

---

## 6. Sequencing

The backend exists to serve the wizard, so build each backend capability with its
UI consumer defined, and verify through the UI — not the CLI.

1. **Define the wire contract first** — pin the exact JSON the wizard needs from
   `validate-sqlite/` and `import-sqlite/`: the full report, the per-project
   `plan`, and per-issue `{type, severity, projects[], detail, resolution}`. This
   is the UI/backend handshake; everything else follows from it.
   **→ Drafted in [`LEGACY_MIGRATION_WIRE_CONTRACT.md`](./LEGACY_MIGRATION_WIRE_CONTRACT.md)**
   (marks every field current vs new; additive-only so the UI can be built first).
2. Backend to satisfy that contract: `_validate_structure` + per-project plan
   (§2.5); `SQLiteImporter` plan-consumption + `copy_files`/mixed-mode (§3a); the
   `import_i2xml` fix (§3b). Land with tests.
3. Frontend: route + wizard (validate → dry-run → import) rendering the plan and
   issues, with the acknowledgement gate. This is the deliverable.
4. Nav button + admin gating (§4).
5. **Verify through the UI**, not the command line: drive the wizard end-to-end
   via `/run` against a real `~/.CCP4I2` (or a fixture) — including a nested-project
   fixture so the mixed-mode plan is exercised in the actual UI.

**Explicitly out of scope this round:** cloud/admin upload flow hardening,
importer de-dup (3c), zip import/export polling.
