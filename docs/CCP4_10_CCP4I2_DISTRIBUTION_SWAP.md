# CCP4i2 in CCP4 10: legacy → pip-installable swap

**Reverse-dependency manifest: what is safe to remove, and what must stay.**
Prepared for the distribution conversation with Charles Ballard.

## TL;DR

- Replace the legacy Qt CCP4i2 (2.6.2) with `pip install ccp4i2` (3.0.0, the
  Django/Qt-free generation).
- **`pip install ccp4i2` is a light, low-risk addition, not a rewrite.** The
  heavy scientific dependencies (`gemmi`, `numpy`, `lxml`, `biopython`) are
  already in `ccp4-python` and are CCP4-native, so they are satisfied with no
  change. The only things added are a handful of **small, pure-Python web
  packages** (`django`, `djangorestframework`, `django-cors-headers`,
  `django-filter`, `uvicorn`, `requests-cache`, `xmlschema`, plus `ccp4i2-api`) —
  no native builds, no version conflict with the scientific stack. Confirm with
  `pip install --dry-run` on a pristine build (see §6).
- **Removable** = the legacy CCP4i2 Qt *application* and the Python Qt *bindings*
  (`PySide2`/`shiboken2`, plus `PyQt4`/`sip` if present) that nothing else uses.
- **Must stay** = the Qt *C++ runtime* (CCP4MG compiles against it), GTK (Coot),
  and all crystallographic toolkits (the new backend lazy-loads them at
  execution).
- So *"remove Qt and its dependencies"* is the one phrase to avoid: it would
  wrongly threaten CCP4MG's Qt C++ runtime, and would not touch Coot (GTK) at
  all. The removable set is narrower and more precise than any toolkit label.
- Net effect: primarily a **cleanliness / maintainability** win. It reclaims the
  legacy i2 tree plus the Python Qt bindings (a real but not dominant chunk); the
  GB-scale crystallographic toolkits stay regardless.

## 1. What is being swapped

| | Legacy | New |
|---|---|---|
| Version | **2.6.2** (`ccp4/ccp4i2:main`) | **3.0.0** (`ccp4/ccp4i2:django`) |
| Architecture | Qt (PySide2) desktop app | Django REST backend + Electron/React frontend |
| GUI toolkit (Python) | PySide2 / shiboken2 | none (frontend is Electron/web) |
| Install | bundled tree in the distribution | `pip install ccp4i2` into `ccp4-python` |

## 2. Safe to remove — legacy-i2-exclusive

**2a. The legacy CCP4i2 Qt application.** On `main`, the GUI application surface
is the `qtcore/` and `qtgui/` trees plus the **235 `.py` modules that import
PySide/PySide2** (and 13 importing PyQt4, 4 importing sip). This is the desktop
application that the Electron/React frontend replaces.

**2b. The Python Qt bindings: `PySide2`, `shiboken2`** (and `PyQt4`/`sip` if
still present). These exist in the distribution to serve legacy i2's Python GUI.

- The **new backend does not import them**: it has **0** module-level Qt imports.
  (It ships 4 files that `import PySide2`, but these are *payload* scripts that
  are copied into a job's work directory and executed *inside CCP4MG's own
  interpreter* — they are never imported by the server, and they travel with the
  `ccp4i2` package, not the interpreter.)
- **Already reflected in the candidate:** the CCP4 10 candidate `ccp4-python`
  tested here has **no `PySide2`, `shiboken2`, `PyQt5` or `sip` installed**. So
  dropping the Python Qt bindings is partly already done — removing legacy i2 is
  consistent with the direction the build has already taken.

## 3. Must stay — shared / load-bearing for the new backend

The new backend runs the crystallographic work out-of-process and still drives
several GUI programs at execution time, so the following are **not** removable:

- **Qt C++ runtime libraries.** CCP4MG is a compiled C++/Qt program, and the new
  i2 still launches it (`ccp4mg_edit_model`). It needs the Qt **C++** libraries,
  *not* the Python bindings — so "remove Qt" must be read as "remove the Python
  Qt *bindings*", never the Qt C++ runtime.
- **GTK.** Coot is a gtkmm program; the new i2 still drives it
  (`coot1`, `coot_rebuild`, the `coot_*` wrappers). GTK is orthogonal to the Qt
  question entirely.
- **Crystallographic toolkits.** `clipper`, `iotbx`, `mmtbx`, `ccp4srs`,
  `phaser`, `dxtbx`, `dials`, `xia2`, `ccp4mg` (the program), `ccp4srs`, etc. The
  new backend **lazy-loads** these only when a job executes. They are undeclared
  by the pip package *by design*, precisely because CCP4 supplies them.

## 4. The pip side: a light, low-risk addition (not a no-op)

> **Important caveat on how this was measured.** The `ccp4-python` inspected here
> already has `ccp4i2` **editable-installed** (`pip install -e .`), so a naive
> "everything is present" reading would be circular — the editable install had
> already pulled `ccp4i2`'s own dependencies in. The table below corrects for
> that using pip's `Required-by` graph to separate genuinely-CCP4-native packages
> from ones the `ccp4i2` install introduced.

**Genuinely CCP4-native — present regardless of ccp4i2, satisfied with no change:**

| dependency | in ccp4-python | required-by (evidence it is CCP4-native) |
|---|---|---|
| gemmi | 0.7.5 | servalcat, modelcraft, nucleofind, ridl, adding_stats_to_mmcif |
| numpy | 1.26.4 | biopython, matplotlib, scipy, scikit-learn, modelcraft, … (≈20 pkgs) |
| lxml | 6.1.1 | widely used; native (libxml2) |
| biopython | 1.87 | adding_stats_to_mmcif, conkit |
| psutil, pytz, future, requests, openpyxl, python-docx | present | general-purpose, broadly used |

**Introduced by the `ccp4i2` / `ccp4i2-api` install — i.e. *would be added* on a
pristine CCP4 distribution** (all small, pure-Python, no native build, no
conflict with the scientific stack):

| package | required-by (shows it rides on ccp4i2) |
|---|---|
| django | ccp4i2-api, djangorestframework, django-cors-headers, django-filter |
| djangorestframework | ccp4i2, ccp4i2-api |
| django-cors-headers | ccp4i2 |
| django-filter | ccp4i2 |
| uvicorn | ccp4i2 |
| requests-cache | ccp4i2 |
| xmlschema | ccp4i2 |
| asgiref | Django (transitive) |
| **ccp4i2-api** | ccp4i2 |

So on a pristine CCP4 `ccp4-python`, `pip install ccp4i2` adds the ~50 MB package
plus that Django web stack (a handful of small pure-Python wheels). It does **not**
touch or rebuild the scientific stack, and nothing pins a version that would force
an upgrade. The definitive check is `pip install --dry-run` on a clean build (§6).

`demo_data` (~400 MB of example datasets) is **not** in the package; it is
distributed out-of-band.

## 5. Electron frontend

Distribute the CI-built, signed apps via download links (keeps the CCP4
distribution lean). Three things to get right:

1. Genuine signing **and notarization** on all three platforms (macOS Gatekeeper
   in particular).
2. A **stable host** for the links (GitHub Releases, not ephemeral CI artifacts).
3. **Version alignment** — the app build must match the `ccp4i2` backend API it
   talks to, so the 3.0.0 backend ships with its corresponding app.

## 6. Verification to run on a *full* CCP4 distribution

The candidate build tested here is a partial Mac unpack (no CCP4MG/Qt present to
inventory), so the toolkit/runtime claims in §3 are architectural, not measured.
Before removing anything, confirm on a complete build:

```bash
# (a) Confirm the Python Qt bindings are used ONLY by legacy i2 (so removable):
grep -rIlE "import (PySide2|shiboken2|PyQt|sip)" "$CCP4" --include='*.py' \
  | grep -v "/ccp4i2/"        # expect: only legacy-i2 paths, nothing else

# (b) Confirm the Qt C++ runtime + CCP4MG + Coot are present and untouched:
find "$CCP4" \( -iname 'libQt*' -o -iname 'QtCore*' -o -iname '*ccp4mg*' \
            -o -iname 'coot*' -o -iname 'libgtk*' \) | head

# (c) See exactly what pip would add on a PRISTINE build (run where ccp4i2 is
#     NOT already editable-installed, or it will report nothing to do):
ccp4-python -m pip install --dry-run ccp4i2
#     expect: ccp4i2 + the small pure-Python Django stack (django, DRF,
#     cors-headers, filter, uvicorn, requests-cache, xmlschema, asgiref,
#     ccp4i2-api); the scientific stack (gemmi/numpy/lxml/biopython) already
#     satisfied. No native builds, no forced upgrades.
```

If (a) returns only legacy-i2 paths, `PySide2`/`shiboken2` are safe to drop with
legacy i2. If anything outside i2 appears, keep them.

## 7. Open question for the release

Is CCP4 10 a **clean cutover** (legacy i2 removed entirely) or a **coexistence**
period (both 2.6.2 and 3.0.0 available for a release)? "Remove all Qt CCP4i2"
implies the former, which is defensible for a major release — but it is better
made an explicit decision than an implied one, given any user workflows still
pinned to the old GUI.
