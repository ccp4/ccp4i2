# ccp4i2

**ccp4i2** is the Python backend of [CCP4i2](https://www.ccp4.ac.uk/) — the
graphical interface and scripting environment for the CCP4 macromolecular
crystallography suite. This is the Django/Qt-free generation of CCP4i2: a REST
API, a task/pipeline framework, and the data-management core that drive the
CCP4i2 desktop (Electron/React) and web deployments.

## CCP4-free control plane

The control plane — the REST API, task configuration, parameter validation, and
the container/digest machinery — imports and runs on a **stock Python with no
CCP4 installation**. Only actual **job execution** (running Refmac, Phaser,
Servalcat, …) needs the full CCP4 suite, and that always happens out-of-process
on a CCP4-bearing worker.

In practice that means you can `pip install ccp4i2` into a plain virtual
environment to configure and validate tasks, serve the API, and drive the
frontend; the crystallographic binaries and their toolkits (clipper, iotbx,
mmtbx, ccp4mg, phaser, dials, xia2, …) are supplied by CCP4 at execution time and
are deliberately **not** pip dependencies.

## Installation

```bash
pip install ccp4i2
```

Optional extras:

```bash
pip install ccp4i2[science]   # rdkit, scipy, pandas — fuller local execution
pip install ccp4i2[dev]       # test/lint toolchain (pytest, mypy, black, …)
```

`ccp4i2` depends on the shared
[`ccp4i2-api`](https://pypi.org/project/ccp4i2-api/) library (authentication +
DRF integration), which is installed automatically.

> **Note** — `demo_data` (the ~400 MB of example datasets used by the test/demo
> workflows) is **not** bundled in the package; it is distributed out-of-band.

## Running

```bash
export DJANGO_SETTINGS_MODULE=ccp4i2.config.settings

# i2run — run a task from the command line
i2run <task> --project_name <proj> [--PARAM value ...]

# development server
python -m uvicorn ccp4i2.config.asgi:application
```

## Links

- CCP4: https://www.ccp4.ac.uk/
- Source: https://github.com/ccp4/ccp4i2
- Documentation: https://www.ccp4.ac.uk/html/INDEX.html

## Licence

LGPL-3.0. See the repository for full licensing details.
