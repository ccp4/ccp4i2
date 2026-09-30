# CCP4i2 user documentation

The task help, as Sphinx/RST. Imported from
[gitlab.com/ccp4i2/rstdocs](https://gitlab.com/ccp4i2/rstdocs), which stays the
help for the Qt interface; this copy is the help for this one.

## Building

```bash
python3 -m venv .venv-docs
.venv-docs/bin/pip install -r docs/user/requirements.txt
.venv-docs/bin/sphinx-build -b html docs/user/source /tmp/userdocs
```

## Pictures from the running app

A page whose figures come from the app has a `shots.json` next to its
`index.rst`, naming what to capture, and a scenario script in `tools/` that
builds the project with real jobs the pictures are taken from. Nothing is
added to `demo_data`; the scenario fetches its data from the PDB and PDB-REDO.

```bash
# 1. The scenario, into a scratch home (never a live one: the script refuses).
cd server
env CCP4I2_HOME=/tmp/docs-home DJANGO_SETTINGS_MODULE=ccp4i2.config.settings \
    ccp4-python manage.py migrate
env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_parrot.py

# 2. The app against that home, on ports of its own.
env CCP4I2_HOME=/tmp/docs-home DJANGO_SETTINGS_MODULE=ccp4i2.config.settings \
    ccp4-python manage.py runserver 127.0.0.1:3421
cd ../client/renderer && env BUILD_TARGET=web \
    NEXT_PUBLIC_API_BASE_URL=http://127.0.0.1:3421 npx next dev -p 3420

# 3. The pictures (Node 22+, Google Chrome; CHROME=path if not the macOS default).
cd docs/user/source/tasks/parrot && node ../../../tools/capture.mjs shots.json
```

`capture.mjs` opens each job page, turns developer mode off, crops to one
section and numbers the fields the text refers to. It finds fields by the label
the user sees, so when a label changes the capture fails and names it: that is
the moment to reread the page's text too.

Pages converted so far: `tasks/parrot`. The rest still show the Qt interface.
