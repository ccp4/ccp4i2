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
- [Error Handling Patterns](../mddocs/pipeline/ERROR_HANDLING_PATTERNS.md)
- [Validity Patterns](../mddocs/pipeline/VALIDITY_PATTERNS.md)

### …use the REST API or CLI
- [API Overview](../mddocs/api/API_OVERVIEW.md)
- [i2run](i2run.md) · [i2remote](i2remote.md) · [CLI Reference](../mddocs/cli/CLI.md)

### …work on the frontend
- [Frontend README](../client/README.md)
- [Frontend Development Guide](../client/FRONTEND_DEVELOPMENT.md)
- [Task Interface Implementation Guide](../client/renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md)
- Moorhen 3D viewer: `client/renderer/MOORHEN_SCENES_SCHEMA_V1_DESIGN.md`

### …understand the architecture
- Project [`CLAUDE.md`](../CLAUDE.md) — architecture, task registry, validation,
  testing, Docker/Azure (the densest single reference)
- [Architecture notes](../mddocs/architecture/)
- [CCP4i Classic Mode Thinking](../CCP4I_CLASSIC_MODE_THINKING.md) — strategy for
  absorbing the legacy CCP4i interface

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
