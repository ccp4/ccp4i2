# CCP4i2 Documentation Map

The single entry point to CCP4i2's developer documentation. Docs live in a few
places by convention:

- **`docs/`** (this folder) — task authoring, CLI, deployment/service, project notes
- **`mddocs/`** — setup, testing, architecture, REST API, pipeline patterns
- **`client/`** — frontend (React/Electron) developer docs
- **`server/ccp4i2/…`** — some deep-dive guides live next to the code they describe

> New to the project? Start with the root [README](../README.md) for setup, then
> come here.

---

## I want to…

### …just try it on my machine
- **As a user (no build):** **[Give it a try](give-it-a-try.md)** — download the
  desktop app, point it at CCP4, let it install the backend.
- **As a developer (from source):**
  [Development Setup](../mddocs/setup/DEVELOPMENT_SETUP.md) — clone, install
  editable into `ccp4-python`, run the dev server + client.

### …add a new task (wrapper)
**Start here → [Authoring a Task](authoring-a-task.md)** — the end-to-end path.
Supporting references:
- [def.xml Reference](def-xml-reference.md) — declare the data model
- [PHIL Task Guide](../server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md) — the PHIL
  alternative for Phenix/PhaserTNG/DIALS-style tools
- [Task Registry: `core/tasks.py`](../server/ccp4i2/core/tasks.py) — the one-line
  registration
- [Task Interface Implementation Guide](../client/renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md)
  — build a bespoke UI (optional; tasks auto-render without one)

### …set up a dev environment / run tests
- [Development Setup](../mddocs/setup/DEVELOPMENT_SETUP.md)
- [Testing Guide](../mddocs/setup/TESTING.md)
- [Writing API tests](writing-api-tests.md) · [Writing i2run tests](writing-i2run-tests.md)

### …work on a pipeline
- [Pipeline Best Practices](pipeline_best_practices.md)
- [Error Handling Remediation](error-handling-remediation.md) — the tracked
  backlog of core-machinery and per-pipeline error-handling defects, with a
  `scripts/scan_error_handling.py` burn-down. **Read before touching a pipeline's
  failure paths**; several sections of Pipeline Best Practices are pending
  amendment against it.
- [Error Handling Patterns](../mddocs/pipeline/ERROR_HANDLING_PATTERNS.md)
- [Validity Patterns](../mddocs/pipeline/VALIDITY_PATTERNS.md)

### …use the REST API or CLI
- [API Overview](../mddocs/api/API_OVERVIEW.md)
- [i2run](i2run.md) · [i2remote](i2remote.md) · [CLI Reference](../mddocs/cli/CLI.md)
- [Value State (design note)](value-state-design.md) — what NOT_SET / DEFAULT /
  EXPLICITLY_SET claim, what they do (a CString is EXPLICITLY_SET from birth,
  and `isSet()` reports non-emptiness for strings but state for ints), and the
  fourth state that is missing: DERIVED, for what the system worked out rather
  than what a user chose. The dependency for leaves ceasing to be
  HierarchicalObjects.
- [CData Simplification (design note)](cdata-simplification.md) — what the
  magic in CData buys, measured, and what a plain-attribute or dataclass
  foundation would have to replace. Distinguishes *composition* (22 classes,
  fixed fields, a record) from *containment* (`CContainer` alone, 70 shapes,
  built from def.xml). Addresses the two standing objections: two-way
  navigation, and construction from def.xml.
- [Container Construction Defects (brief)](container-construction-defects.md) —
  two defects that make jobs run and be wrong: the def.xml parser duplicates
  nested containers onto the root (2,304 ghost parameters), and `children()`
  iterates a weakref set so sibling order — and therefore which parameter a bare
  `--FLAG` means — is a coin flip. Evidence, measured blast radius, proposed
  fixes. Read alongside
  [error-handling-remediation.md](error-handling-remediation.md).
- [i2run Parameter Syntax (proposal)](i2run-parameter-syntax.md) — a single
  addressing rule to replace the current three, a declared default field so
  `--XYZIN beta.pdb` works at every depth, ambiguity and undeclared fields as
  errors, and `--describe`/`--echo`. Assumes the container fixes above. Includes
  a measured argument for why PHIL's *behaviours* are worth adopting and its
  object model is not.

### …work on the frontend
- [Frontend README](../client/README.md)
- [Frontend Development Guide](../client/FRONTEND_DEVELOPMENT.md)
- [Task Interface Implementation Guide](../client/renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md)
- Moorhen 3D viewer: `client/renderer/MOORHEN_SCENES_SCHEMA_V1_DESIGN.md`

### …organise projects (tags, groups, campaigns)
- [Organising Projects](organising-projects.md) — `ProjectTag` vs `ProjectGroup`:
  which to reach for, why hierarchy lives on tags, and where legacy Qt-i2
  "folders" map to

### …understand the architecture
- Project [`CLAUDE.md`](../CLAUDE.md) — architecture, task registry, validation,
  testing, Docker/Azure (the densest single reference)
- [Architecture notes](../mddocs/architecture/)
- [CCP4i Classic Mode Thinking](../CCP4I_CLASSIC_MODE_THINKING.md) — strategy for
  absorbing the legacy CCP4i interface

### …handle configuration and secrets
- [Preferences Handling](preferences-proposal.md) — where non-secret config lives
  (`~/.ccp4i2-django/preferences.json`, env-first precedence)
- [Handling Secrets](CREDENTIALS_DESIGN.md) — API tokens, ssh passwords and key
  passphrases: the credential store, its REST surface, and why a credential is
  never a task parameter
- [Web-API Task Guide](../server/ccp4i2/wrappers/WEB_API_TASK_GUIDE.md) — the
  practical how-to for a task that calls an authenticated service
- [Remote Job Execution Plan](REMOTE_JOB_EXECUTION_PLAN.md) — ssh/qsub/SLURM
  dispatch, which reuses the same credential machinery

### …cut a release
- [Releasing CCP4i2](RELEASING.md) — one `v*` tag → `ccp4i2` on PyPI + desktop
  installers on a GitHub Release (OIDC, no tokens).
- [GitHub Actions](github-actions.md) — what the four workflows run and when,
  what is and is not gated, and the gotchas that bite when changing them.

### …deploy
- Docker / Azure sections in [`CLAUDE.md`](../CLAUDE.md)
- [Service Contract](CCP4I2_SERVICE_CONTRACT.md)

---

## Legacy / historical (do not follow as current)

These describe removed or superseded mechanisms and are kept only for reference:

- [Plugin Registry](../mddocs/PLUGIN_REGISTRY_README.md) — removed
  `plugin_registry.py` / `plugin_lookup.json` scan → now `core/tasks.py`
- [Stub Modules](../mddocs/STUBS_README.md) — for the removed `plugin_lookup.py`
- [Qt Task GUI Guide](../mddocs/qt_task_gui_guide.md) — Qt UIs → now React /
  `GenericInterface`
- [`mddocs/migration/`](../mddocs/migration/) — in-flight refactor logs
