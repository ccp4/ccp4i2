# Give it a try (desktop app — no build)

The fastest way to try the new CCP4i2 on your own machine. No cloning, no
building — download the desktop app, point it at your CCP4 installation, and let
it set up the backend for you.

> **Two ways to "give it a go" — pick your path:**
> - **As a user (this page):** download the packaged desktop app. Nothing to
>   build; the app installs the backend for you.
> - **As a developer (run from source):** clone the repo and install editable
>   into `ccp4-python`, then run the dev server + client — see
>   [Development Setup](../mddocs/setup/DEVELOPMENT_SETUP.md).

> **You do not run `pip` yourself.** The CCP4i2 backend must live inside CCP4's
> own Python (`ccp4-python`), not your system Python. The app installs it into
> `ccp4-python` for you (see step 4). A bare `pip install ccp4i2` against system
> Python will **not** work.

---

## 1. Install CCP4 (provides `ccp4-python`)

Download and install a **CCP4 10** build from
<https://ccp4serv6.rc-harwell.ac.uk/10/downloads/>. This gives you the CCP4
suite and, with it, `ccp4-python` — the interpreter the backend runs in.

You do **not** need to source any setup script for the desktop app; it locates
your CCP4 installation itself (step 3).

## 2. Download the desktop app (from CI)

Builds are produced by the **Electron Multiplatform Build** workflow on the
`django` branch. On GitHub (you'll need to be signed in):

**ccp4/ccp4i2 → Actions → "Electron Multiplatform Build" → latest run on
`django` → Artifacts**, then download the one for your OS:

| OS | Artifact | Contains |
|---|---|---|
| macOS | `macOS-dmg` | a `.dmg` |
| Windows | `windows-installers` | a `.exe` installer |
| Linux | `linux-appimages` | an `.AppImage` |

GitHub delivers artifacts as a `.zip` — unzip it to get the installer/app.

### macOS: clear the quarantine flag

The build isn't notarised, so macOS Gatekeeper will quarantine it. After
mounting the `.dmg` and copying **ccp4i2-django.app** to `/Applications` (or
wherever you like), clear the quarantine attribute:

```bash
xattr -cr "/Applications/ccp4i2-django.app"
```

Then launch it (first launch: right-click → **Open** if prompted).

## 3. Launch and point it at CCP4

Start the app. It auto-detects CCP4 in the usual locations; if it can't find
yours, point it at your CCP4 installation directory when asked. It uses that
installation's `ccp4-python` to run the Django backend — no Python is bundled in
the app.

## 4. Let the app install the backend

On startup the app checks whether a new-enough `ccp4i2` is present in
`ccp4-python`. If it's **missing or too old**, you'll see an **Install** button.
Click it — the app installs/upgrades `ccp4i2` (and a curated dependency lock)
**into `ccp4-python`** for you. Wait for it to finish.

> This is why you don't run pip yourself: the app targets the right interpreter
> and pins a compatible version automatically.

## 5. Play

Create a project, import some data (or use the demo data), pick a task from the
task chooser, and run it. Try the MTZ previewer and the Moorhen 3D viewer on the
outputs.

---

## Versions — do I need to match them?

**No, not by hand.** The desktop app carries a **minimum backend version floor**
and installs/upgrades `ccp4i2` in `ccp4-python` to satisfy it (step 4). You don't
pick versions.

Two things worth knowing:

- **Always let the app do the install.** Installing `ccp4i2` yourself risks
  putting it in the wrong Python or at an incompatible version.
- **The published PyPI package can lag the `django` branch.** The app installs
  `ccp4i2` from PyPI, while the desktop build comes from `django` CI. They're
  normally aligned, but if you hit an unexpected mismatch, that's the likely
  cause — grab the latest desktop build and re-run the in-app Install.

## Troubleshooting

- **"Installed ccp4i2 … is older than the required …"** — click **Install** to
  upgrade; that's the version floor doing its job.
- **App won't find CCP4** — point it explicitly at your CCP4 install directory.
- **macOS "app is damaged / can't be opened"** — you missed the `xattr -cr`
  step above (or ran it on the wrong path).
- **Install step fails** — confirm the CCP4 install is complete and its
  `ccp4-python` runs (`<ccp4>/bin/ccp4-python --version`).

## Prefer to run from source? (developer path)

This page is the **user** path. To try CCP4i2 **as a developer** — clone the
repo, install editable into `ccp4-python`, run the dev server and Electron client
— follow [Development Setup](../mddocs/setup/DEVELOPMENT_SETUP.md), then
[Authoring a Task](authoring-a-task.md) when you're ready to add one.
