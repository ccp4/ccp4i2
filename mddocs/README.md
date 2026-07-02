# CCP4i2 Documentation

Reference documentation for the ccp4i2-django project. This folder is for live reference material only; in-flight planning notes, progress reports, and one-time refactor logs are not kept here (they live in commits and PR descriptions).

> **Looking for the full index?** See the top-level
> **[Documentation Map](../docs/README.md)**, which links every doc tree
> (`docs/`, `mddocs/`, `client/`, and code-adjacent guides) in one place.

## Getting Started

- **[Development Setup](setup/DEVELOPMENT_SETUP.md)** — Environment setup for development and testing
- **[Testing Guide](setup/TESTING.md)** — Running and writing tests

## Authoring a Task

- **[Authoring a Task](../docs/authoring-a-task.md)** — end-to-end guide (def.xml → wrapper → registration → UI)
- **[def.xml Reference](../docs/def-xml-reference.md)** — declaring a task's data model
- **[PHIL Task Guide](../server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md)** — the PHIL alternative for tools with a `master_phil`
- Registration: the `TASKS` dict in [`server/ccp4i2/core/tasks.py`](../server/ccp4i2/core/tasks.py)

## Architecture

- **[Quick Reference](QUICK_REFERENCE.md)** — Async execution infrastructure reference
- [architecture/](architecture/) — System design documentation

### Legacy / historical (not current — kept for reference)

- **[Plugin Registry](PLUGIN_REGISTRY_README.md)** ⚠️ — describes the removed `plugin_registry.py` / `plugin_lookup.json` scan; registration is now the `TASKS` dict in `core/tasks.py`
- **[Stub Modules for Plugin Discovery](STUBS_README.md)** ⚠️ — stubs for the removed `plugin_lookup.py`
- **[Stub/Implementation Inheritance Pattern](STUB_IMPLEMENTATION_INHERITANCE_PATTERN.md)** — full-fat classes inheriting from generated stub + parent (MRO concern)
- **[Qt-Based Task GUI Guide](qt_task_gui_guide.md)** ⚠️ — the legacy Qt task GUI system; UIs are now React / `GenericInterface`
- [migration/](migration/) — in-flight refactor logs (historical)

## REST API

- **[API Overview](api/API_OVERVIEW.md)** — REST API endpoints and usage
- [api/](api/) — Detailed API documentation

## Pipeline Development

- **[Error Handling Patterns](pipeline/ERROR_HANDLING_PATTERNS.md)** — CErrorReport usage, try/except, ERROR_CODES
- **[Validity Patterns](pipeline/VALIDITY_PATTERNS.md)** — Content-aware validation via `validity()` overrides

## Frontend Development

- **[Frontend README](../client/README.md)** — Quick start and overview
- **[Frontend Development Guide](../client/FRONTEND_DEVELOPMENT.md)** — Comprehensive developer documentation

## CLI Documentation

- **[CLI Reference](cli/CLI.md)** — Full command-line interface documentation
- **[CLI Quick Reference](cli/CLI_QUICK_REFERENCE.md)** — One-page CLI cheat sheet
- **[i2run Guide](cli/I2RUN_GUIDE.md)** — Task runner with complex data object construction
