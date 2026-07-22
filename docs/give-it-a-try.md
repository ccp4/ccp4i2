# Give it a try (desktop app — no build)

The fastest way to try the new CCP4i2 on your own machine. No cloning, no
building — download the desktop app, point it at your CCP4 installation, and let
it set up the backend for you.

> 🐣 **Early-bird alpha.** The current build is **`3.1.0a1`** — a pre-release for
> early testers. Expect rough edges, and please report what breaks. The desktop
> app and its backend are **version-locked**: an alpha app installs and runs
> **only** its exact matching `ccp4i2` from PyPI, so it can't get crossed with
> any other version. (Any earlier `3.0.x` you may have installed is a separate,
> superseded line — this alpha won't use it.)

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

## 2. Download the desktop app

Go to **[Releases](https://github.com/ccp4/ccp4i2/releases)** and pick the latest
build marked **Pre-release** (the alpha, e.g. `v3.1.0a1`), then download the
installer for your OS — no GitHub login needed:

| OS | Asset |
|---|---|
| macOS | `.dmg` |
| Windows | `.exe` installer |
| Linux | `.deb` (Ubuntu/Debian — **recommended**) or `.AppImage` (portable) |

> During the alpha, releases are **pre-releases**, so they won't appear under
> "Latest release" — pick the top entry tagged *Pre-release* on the Releases
> page.

<details>
<summary>Want the bleeding-edge build instead of a release?</summary>

Every push to `django` also produces installers as CI artifacts (these require a
GitHub login and expire): **ccp4/ccp4i2 → Actions → "Electron Multiplatform
Build" → latest run → Artifacts** (`macOS-dmg`, `windows-installers`,
`linux-appimages`). GitHub delivers artifacts as a `.zip` — unzip to get the
installer.
</details>

### macOS: clear the quarantine flag (unsigned builds only)

If the build you downloaded is **unsigned/unnotarised**, macOS Gatekeeper will
quarantine it. After mounting the `.dmg` and copying **ccp4i2-django.app** to
`/Applications` (or wherever you like), clear the quarantine attribute:

```bash
xattr -cr "/Applications/ccp4i2-django.app"
```

Then launch it (first launch: right-click → **Open** if prompted). If the app
opens normally without any warning, it was signed and notarised and you can skip
this step.

### Linux: install the `.deb` (recommended) or run the AppImage

**Ubuntu / Debian — install the `.deb` (recommended):**

```bash
sudo apt install ./ccp4i2-django_*.deb
```

This installs to `/opt/ccp4i2-django`, adds an application-menu entry and a
`ccp4i2-django` command, and sets up the Chromium sandbox the same way Chrome
and VS Code do — so it runs on **modern locked-down kernels (Ubuntu 24.04+)**
with a full sandbox and **no workarounds**. This is the best option on a
desktop; prefer it if you can `sudo`.

**Other distros (Fedora, RHEL/Rocky, openSUSE, Arch, …), or no root — the
portable `.AppImage`:** downloaded AppImages aren't executable by default (your
browser doesn't set the flag), so mark it once:

```bash
chmod +x ccp4i2-django-*.AppImage
```

(or right-click → **Properties** → **Permissions** → **Allow executing file as
program**), then double-click or run `./ccp4i2-django-*.AppImage`. On these
distros the kernel permits the unprivileged-user-namespace sandbox, so the app
runs **fully sandboxed with no flags and no system-wide change**.

> **On Ubuntu 24.04+ / Debian 13, use the `.deb`, not the AppImage.** These
> kernels restrict the namespace sandbox (see below), and the AppImage cannot
> supply the privileged fallback, so launching it aborts with
> `The SUID sandbox helper binary … is not configured correctly`. The AppImage
> deliberately does **not** force `--no-sandbox` (that would silently strip the
> Chromium sandbox). To run the AppImage here you must re-enable unprivileged
> user namespaces yourself:
>
> ```bash
> sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
> ```
>
> after which it runs **fully sandboxed**. The `.deb` sandboxes correctly with no
> workaround at all, which is why it is the recommended artifact here.

<details>
<summary><strong>Why the <code>.deb</code>? (and why the AppImage struggles on hardened kernels)</strong></summary>

Chromium — which Electron, and therefore this app, is built on — refuses to run
its rendering processes unless it can put them in a **sandbox**. It can build
that sandbox one of two ways, and it needs *one* of them:

1. a small helper (`chrome-sandbox`) that is **allowed to run as root** (the
   "SUID sandbox"), or
2. **unprivileged user namespaces** — a kernel feature letting ordinary programs
   create isolated environments without being root.

Modern "locked-down" distros (Ubuntu 24.04+) **switch off option 2** by default
for security (`kernel.apparmor_restrict_unprivileged_userns=1`). That leaves the
app needing option 1 — a root-privileged helper — and here the packaging format
is decisive:

| | AppImage | `.deb` |
|---|---|---|
| Runs an install step **as root**? | No — it's a self-mounting single file | Yes — `dpkg`/`apt` runs its `postinst` as root |
| Can make `chrome-sandbox` setuid-root? | No — it unpacks to a `nosuid` FUSE mount in `/tmp`, and nothing runs as root to `chown`/`chmod` it | **Yes** — `postinst` does `chown root:root` + `chmod 4755` |
| SUID sandbox available? | No | **Yes** (the helper is privileged, so it needs no unprivileged userns) |
| Namespace sandbox available on 24.04+? | No (kernel restriction) | Not needed |
| Result on a locked-down kernel | Sandbox can't start, so it **aborts** — unless you re-enable unprivileged user namespaces (`sysctl … apparmor_restrict_unprivileged_userns=0`), after which it runs fully sandboxed | **Full sandbox, no flags, no user action** |

The load-bearing sentence: **the sandbox needs a privilege it can obtain only two
ways; the locked-down kernel removes one of them, and only a package with a root
install step can supply the other.** An AppImage has no root install step; a
`.deb` does.

The proof it works: **Google Chrome and VS Code run fine on these same kernels**
— because they ship as `.deb`s whose installer does this exact setuid step. Our
`.deb` does the same via electron-builder's standard Linux `postinst` (it
setuids `chrome-sandbox` only on kernels that need it, installs an AppArmor
profile, and links `ccp4i2-django` onto `PATH`). So the AppImage isn't broken so
much as structurally unable to sandbox on a hardened desktop; the `.deb` is the
right tool, and the AppImage remains a portable fallback for other distros or
where you can't `sudo`.

</details>

## 3. Launch and point it at CCP4

Start the app. It auto-detects CCP4 in the usual locations; if it can't find
yours, point it at your CCP4 installation directory when asked. It uses that
installation's `ccp4-python` to run the Django backend — no Python is bundled in
the app.

## 4. Let the app install the backend

On startup the app checks whether the **exact** `ccp4i2` it needs is present in
`ccp4-python`. If it's **missing or a different version**, you'll see an
**Install** button. Click it — the app installs the matching `ccp4i2` (and a
curated dependency lock) **into `ccp4-python`** for you. Wait for it to finish.

> This is why you don't run pip yourself: the app targets the right interpreter
> and pins a compatible version automatically.

## 5. Play

Create a project, import some data (or use the demo data), pick a task from the
task chooser, and run it. Try the MTZ previewer and the Moorhen 3D viewer on the
outputs.

---

## Versions — do I need to match them?

**No, not by hand — and during the alpha they're locked together.** Each alpha
desktop build is pinned to one **exact** `ccp4i2` version (e.g. `3.1.0a1`), built
from the same release tag, and its Install button installs *precisely* that
version into `ccp4-python` (step 4). So the app and backend can't drift apart,
and the app won't run against a different `ccp4i2`.

- **Always let the app do the install.** Installing `ccp4i2` yourself risks
  putting it in the wrong Python or at the wrong version. (And a plain
  `pip install ccp4i2` won't even fetch an alpha — pre-releases are hidden from
  normal installs by design.)
- **To move to a newer alpha, grab the newer desktop build** from Releases and
  click Install again — it switches `ccp4-python` to the new matching backend.

## Troubleshooting

- **"Installed ccp4i2 … does not match the version this app requires …"** —
  click **Install**; the app is version-locked to one exact backend and is
  setting it up.
- **App won't find CCP4** — point it explicitly at your CCP4 install directory.
- **macOS "app is damaged / can't be opened"** — you missed the `xattr -cr`
  step above (or ran it on the wrong path).
- **Linux: double-click does nothing / only offers "Open With"** — the AppImage
  isn't executable yet; `chmod +x` it (see the Linux step above).
- **Linux: the AppImage crashes or shows a blank white window on a locked-down
  kernel** (e.g. Ubuntu 24.04+, "The SUID sandbox helper binary … is not
  configured correctly") — **install the `.deb` instead** (`sudo apt install
  ./ccp4i2-django_*.deb`). The `.deb` gives Chromium a proper sandbox the way
  Chrome/VS Code do; the AppImage can't, and is only a portable fallback. Do
  **not** weaken the kernel with `sudo sysctl … unprivileged_userns=0`.
- **Install step fails** — confirm the CCP4 install is complete and its
  `ccp4-python` runs (`<ccp4>/bin/ccp4-python --version`).

## Prefer to run from source? (developer path)

This page is the **user** path. To try CCP4i2 **as a developer** — clone the
repo, install editable into `ccp4-python`, run the dev server and Electron client
— follow [Development Setup](../mddocs/setup/DEVELOPMENT_SETUP.md), then
[Authoring a Task](authoring-a-task.md) when you're ready to add one.
