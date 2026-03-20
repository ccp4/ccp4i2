# CCP4i2-Django Release Checklist

Tracking document for countdown to initial release. Each section captures
outstanding work, open questions, and acceptance criteria.

---

## 1. Core Functionality Gaps

Must-have features that are missing or incomplete before day-0 release.

- [ ] **Task coverage** — all 63 GUI-visible tasks run end-to-end via i2run
      and REST API (currently 45 at full coverage; see
      [test-coverage-tracker.md](test-coverage-tracker.md) for the live matrix)
- [ ] **Project management** — create / open / archive / delete projects
- [ ] **Data import** — drag-and-drop / file-browser import of MTZ, PDB, mmCIF,
      sequence files
- [ ] **Job history & provenance** — browse previous jobs, re-run with modified
      parameters, trace data lineage
- [ ] **Report viewing** — HTML reports render correctly in Electron and web
      (Moorhen, graphs, tables)
- [ ] **Export** — export project or selected results (PDB/mmCIF + MTZ + report)
- [ ] **User preferences** — persistent per-user settings (default space group,
      column labels, etc.)
- [ ] **Undo / job cancellation** — cancel a running job; clean up partial
      outputs
- [ ] **Notifications** — desktop/web notification on job completion or failure
- [ ] **Help integration** — context-sensitive help links from task dialogs to
      CCP4 docs

> **Action:** audit the above against the classic CCP4i2 feature list and flag
> anything else users would expect on day 0.

---

## 2. Edge-Case & Deployment Scenarios

Environments that must work at launch, beyond the happy-path desktop case.

- [ ] **SSH / remote X** — Electron app launched over SSH with X-forwarding;
      verify rendering, file-dialog behaviour, and latency
- [ ] **Headless CLI (i2run)** — runs without display (`DISPLAY` unset); all
      tasks callable from scripts
- [ ] **Web / cloud mode** — Docker-deployed server + web client accessible from
      browser (Chrome, Firefox, Safari); auth via Azure AD
- [ ] **i2remote over SSH tunnel** — `i2remote` CLI targeting a server via
      `ssh -L` port forward; token refresh, upload/download of large files
- [ ] **Shared filesystems (NFS / GPFS)** — concurrent project access from
      multiple machines; file-locking, symlink handling
- [ ] **Mixed CCP4 versions** — user has both CCP4-8 and CCP4-10 installed;
      verify PATH/LD_LIBRARY_PATH isolation
- [ ] **Non-Latin paths** — project directory containing spaces, Unicode
      characters, or very long path names
- [ ] **Low-bandwidth / high-latency** — web mode over a slow connection;
      sensible loading states, no silent timeouts
- [ ] **Windows** — if Windows support is in scope for day 0, test Electron
      build + CCP4 environment detection on Windows

> **Action:** decide which of the above are day-0 vs. post-release and mark
> accordingly.

---

## 3. Data Migration from CCP4i2-Classic

Strategy for bringing forward existing user projects and results.

- [ ] **Migration script** — command-line tool that reads a classic CCP4i2
      SQLite database + project directory and imports into the Django database
  - Maps classic task UUIDs → new task records
  - Copies or symlinks data files
  - Preserves job history and parent–child relationships
- [ ] **Selective import** — option to import only specific projects or jobs
- [ ] **Validation report** — after migration, produce a summary of what was
      imported, what was skipped, and any warnings
- [ ] **Round-trip safety** — migration is non-destructive to the source; users
      can keep running classic CCP4i2 in parallel
- [ ] **Documentation** — step-by-step migration guide for users
- [ ] **Testing** — run migration against at least 3 real-world project
      databases of varying size (small demo, medium lab, large synchrotron
      session)

> **Open question:** do we migrate just metadata, or also re-register data
> files? What about jobs that reference programs no longer available?

---

## 4. Testing Rollout

Plan for expanding automated and manual test coverage, and recruiting external
testers.

### 4a. Automated Testing

- [ ] **CI pipeline** — GitHub Actions (or equivalent) running the full test
      suite on every push to `django` branch
- [ ] **i2run tests** — target 100% coverage of GUI-visible tasks
      (currently 49/63)
- [ ] **API tests** — target 100% coverage (currently 54/63)
- [ ] **mmCIF compatibility** — all wrappers tested with mmCIF input where
      applicable; eliminate all `!` entries in tracker
- [ ] **Integration tests** — multi-step pipeline tests (e.g. import → process
      → refine → report)
- [ ] **Frontend tests** — basic smoke tests for React components (task
      dialogs render, form submission works)
- [ ] **Performance benchmarks** — baseline timings for common workflows;
      detect regressions

### 4b. Manual / Exploratory Testing

- [ ] **Internal dogfooding** — team uses django branch for day-to-day work
      for at least 2 weeks before external release
- [ ] **Friendly testers (alpha)** — recruit 3–5 crystallographers who are
      comfortable reporting bugs; provide them a feedback channel
      (GitHub Issues? shared spreadsheet?)
- [ ] **Wider beta** — open to a broader group (e.g. CCP4 workshop
      participants); collect structured feedback via a form
- [ ] **Hostile testers** — if any known power-users or sceptics are willing,
      invite them early; their edge-case usage patterns are invaluable

### 4c. Test Infrastructure

- [ ] **Test data repository** — curated set of small MTZ/PDB/mmCIF files
      covering common space groups, anomalous data, twinned data, etc.
- [ ] **Reproducible environment** — Docker image or script that sets up a
      clean CCP4 + ccp4i2-django for testing

> **Action:** set a target date for alpha tester onboarding.

---

## 5. Documentation & Communications Strategy

Drive understanding among developers and uptake among users.

### 5a. Developer-Facing

- [ ] **Architecture overview** — update/expand CLAUDE.md into a standalone
      developer guide covering backend, frontend, data flow
- [ ] **Contributing guide** — how to set up a dev environment, run tests,
      submit PRs
- [ ] **Wrapper writing guide** — expand existing
      [pipeline_best_practices.md](pipeline_best_practices.md) and
      `PHIL_TASK_GUIDE.md`
- [ ] **API reference** — auto-generated or hand-written docs for the REST API
      endpoints
- [ ] **Plugin system docs** — how to register a new task, where files go,
      how the registry works

### 5b. User-Facing

- [ ] **"What's new" summary** — one-page document explaining what
      ccp4i2-django is, why it exists, and what's different from classic
- [ ] **Quick-start guide** — install CCP4-10, launch ccp4i2, process your
      first dataset (with screenshots)
- [ ] **Migration guide** — see section 3 above
- [ ] **FAQ / troubleshooting** — common issues (CCP4 not found, port
      conflicts, browser compatibility)
- [ ] **Video walkthrough** — short (5–10 min) screencast of a typical workflow

### 5c. Outreach

- [ ] **CCP4 newsletter article** — announce the new interface
- [ ] **Workshop materials** — hands-on tutorial for CCP4 study weekends /
      workshops
- [ ] **Mailing list announcement** — ccp4bb post at appropriate moment

> **Action:** identify who owns each documentation deliverable.

---

## 6. CCP4-10 Build Integration

Plan for incorporating ccp4i2-django into official CCP4 distribution builds.

### 6a. Build System

- [ ] **Build script** — ccp4i2-django installs cleanly via CCP4 `configure` /
      build system (autotools? CMake? pip install?)
- [ ] **Node.js dependency** — CCP4 builds currently don't ship Node.js;
      decide whether to bundle it, use a pre-built Electron binary, or ship
      only the web mode in CCP4-10
- [ ] **Electron packaging** — produce platform-specific Electron app bundles
      (macOS .app, Linux AppImage/deb, Windows .exe) as part of the build
- [ ] **Entry point** — `ccp4i2` command launches the new interface; classic
      interface available as `ccp4i2-classic` (or vice versa during transition)

### 6b. Python Dependency Conflicts

Known or anticipated clashes with CCP4's bundled Python environment:

| Package | CCP4-bundled version | ccp4i2-django needs | Risk |
|---------|---------------------|---------------------|------|
| Django | not bundled | 4.2+ LTS | New dependency — check against CCP4 Python version (3.9+) |
| djangorestframework | not bundled | 3.14+ | New dependency |
| numpy | varies | compatible | Low — usually fine |
| gemmi | CCP4-bundled | compatible | Low — we use whatever CCP4 provides |
| requests | may be present | 2.28+ | Low — check version |
| msal | not bundled | optional (auth) | Only needed for cloud/Azure mode |
| websockets | not bundled | if used | Check — may conflict with other CCP4 tools |

- [ ] **Audit** — run `pip list` in ccp4-python and diff against our
      `requirements.txt`; flag version conflicts
- [ ] **Virtual environment strategy** — decide whether ccp4i2-django gets its
      own venv inside CCP4 or installs into the shared environment
- [ ] **Pinning** — pin our dependencies to versions tested against the target
      CCP4-10 Python (currently 3.11/3.12)

### 6c. Coexistence with Classic CCP4i2

- [ ] **Side-by-side installation** — both classic (Qt) and django versions
      installable; no file conflicts
- [ ] **Database isolation** — django uses its own DB; does not touch classic
      SQLite files
- [ ] **Shared data directory** — both versions can read the same project data
      files (read-only from classic's perspective)

> **Decision:** CCP4-10 will ship *only* the Django/React ccp4i2. Classic
> (Qt) CCP4i2 will not be included.

---

## Milestone Targets

| Milestone | Target Date | Key Criteria |
|-----------|-------------|--------------|
| Internal alpha | TBD | All core tasks run; team dogfooding |
| External alpha | TBD | Migration script works; 5 friendly testers |
| Beta | TBD | CI green; docs draft complete; wider testing |
| RC1 | TBD | All blockers resolved; migration tested on real data |
| Release (CCP4-10.x) | TBD | Integrated into build; announcement published |

---

## Post-Release

- [ ] **Disentangle crank2** — extract crank2 from the ccp4i2 codebase into
      its own standalone package / repository

*This is a living document. Update it as items are completed or priorities
change.*
