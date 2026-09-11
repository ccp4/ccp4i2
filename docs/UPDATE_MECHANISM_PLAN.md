# Post-distribution update mechanism: planning document

**How CCP4i2 stays up to date once it ships inside CCP4 — the Python backend
and the Electron desktop app, which update on different terms.**

Prepared for the core-team distribution conversation. Sibling to
[CCP4i2 in CCP4 10: legacy → pip-installable swap](CCP4_10_CCP4I2_DISTRIBUTION_SWAP.md),
which covers the *initial* swap; this covers what happens *after* it ships.

## TL;DR

- CCP4i2 is two artifacts that must stay compatible: the **Python backend**
  (`ccp4i2` in `ccp4-python`'s site-packages) and the **Electron desktop app**
  (a per-platform, checksummed installer). They differ in size, signing needs,
  and cadence, so they should update on **independent channels**.
- The first decision is not "which downloader" but **coupling policy**. Today
  the app pins the backend *exactly* and refuses any other version — an
  alpha-phase discipline. Updates cannot move independently until this is
  relaxed to a **compatible range** under a stable API contract. The code
  already anticipates this ([`ccp4i2-server-version.ts`](../client/main/ccp4i2-server-version.ts)).
- **Python side → CCP4 UM.** Once `ccp4i2` lives in CCP4 site-packages, its
  updates *are* CCP4 updates; the CCP4 Update Manager's
  background-download-and-surface model is the right owner. Conditional on the
  upstream metadata/dependency fixes (see §2) so UM's update needs none of the
  install-time repair the desktop button does today.
- **Electron side → its own path.** No auto-update exists today. Start with a
  lightweight **version-gate prompt** (no new infrastructure), and treat
  `electron-updater` from a published feed as the real target — **gated on code
  signing** (macOS notarization + Windows sign, [#238](https://github.com/ccp4/ccp4i2/issues/238)).
- Do **not** fold the Electron app into UM. One channel for both buys tidiness
  at a real cost (UM would have to understand per-platform Electron packaging
  and swap a large signed bundle).

## 1. Current state (what actually exists today)

| Concern | Today |
|---|---|
| Backend version binding | **Exact pin.** The packaged app installs and requires `ccp4i2==<exact>`; it refuses any other version. `CCP4I2_REQUIRED_SERVER_VERSION`, stamped at build time from the release tag. |
| Backend install route | The desktop **Install** button runs `pip` into `ccp4-python` with a metadata-repair preflight and a resolver-free two-step install (see §2 and the [distribution-swap doc](CCP4_10_CCP4I2_DISTRIBUTION_SWAP.md)). |
| App auto-update | **None.** No `build.publish` feed, no `electron-updater` dependency, no `autoUpdater` code. Updating the app means downloading a new installer by hand. |
| Release production | One `v*` tag → `release.yml` builds the PyPI wheel + mac/win/linux installers on a GitHub Release. See [RELEASING.md](RELEASING.md). |
| App vs Python version | Already **decoupled** as numbers — the Electron `package.json` version (`0.0.1`) is deliberately independent of the Python package version. |

The exact pin is intentional and documented: *"Exact pinning trades
auto-pick-up-newer for cannot-silently-mix-versions, which is what we want while
the software is not yet world-ready. Relax to a floor / compatible-range once a
stable line is declared."* That sentence is the hinge this whole plan turns on.

## 2. Why the Python update is not just "pip install ccp4i2"

The desktop install today performs repairs because a hand-rolled `ccp4-python`
breaks stock pip in two ways (full write-up:
[`client/main/ccp4i2-ipc.ts`](../client/main/ccp4i2-ipc.ts), `REPAIR_METADATA_PY`
and `runInstall`):

1. **Corrupt distribution metadata.** Several distributions already in CCP4's
   site-packages ship `*.dist-info`/`*.egg-info` metadata with no parseable
   `Version:` field (observed: `gyp_next`, `meson`, `SCons`,
   `typing_extensions` — varies build-to-build). pip enumerates the *entire*
   installed set and parses every version; a single missing one raises
   `TypeError: … got 'NoneType'` and aborts pip outright. The installer stamps a
   synthetic `Metadata-Version`/`Name`/`Version` header into any such file first
   (idempotent, build-agnostic).
2. **A stale dependency stack.** CCP4's Python is pinned behind what the Django
   backend needs (notably Django/asgiref), so the installer never invokes pip's
   resolver: `--no-deps` for the wheel, then `--no-deps --ignore-installed -r
   requirements-runtime.txt` for a curated lock that deliberately *excludes*
   CCP4's ABI-native packages (numpy/gemmi/lxml), so nothing compiled is
   overwritten. Success is judged by an **import probe**, not pip's exit code
   (pip can exit non-zero on the benign post-install metadata crash after the
   package is already in place).

**Implication for the integrated build.** Both problems are upstream in the
CCP4 build, not in CCP4i2. If UM is to own the Python update cleanly, the two
cleanest wins are to **fix the missing-`Version` metadata at the CCP4 build
source** and to **align the Django/asgiref pins** in the shipped environment.
Then a UM update of `ccp4i2` needs no repair machinery at all, and the desktop
Install button's repair path can be **relegated to a fallback** (for
non-integrated CCP4s and self-repair) rather than the primary route.

## 3. The coupling decision (do this first)

Everything downstream depends on this one choice.

- **Option A — stay locked (exact pin).** Every update event ships *both* app
  and Python together. Correctness is trivial (versions can never mismatch), but
  every backend fix drags a full per-platform app rebuild and redownload. Fine
  for alpha; too heavy for a shipped product.
- **Option B — decouple to a compatible range.** The backend moves under a
  declared-stable API contract while the app stays put; the app checks a
  *floor/range* (`>=`) instead of `==`. This is what makes "UM updates Python
  quietly, app updates rarely" possible.

**Recommendation: B.** The Electron app is a thin shell (UI + a process that
spawns the `ccp4-python` Django server); the crystallography all lives in the
Python wheel, which will churn far more than the shell. Binding the shell's
cadence to the backend's is the expensive choice. The enabling discipline is
already in place: the shared **`@ccp4/ccp4i2-api` contract** is where
range-compatibility is promised and tested. Switching the pin is a small code
change (`==` → `>=` via the existing `compareVersions`); the *work* is the
contract discipline that makes it safe.

## 4. Python side — CCP4 UM

Once `ccp4i2` is in CCP4 site-packages, its updates are CCP4 updates and UM is
the natural owner: background download, surfacing, and the update cadence CCP4
users already expect. Preconditions:

1. The §2 upstream fixes land, so UM's update path needs no pip repair.
2. The desktop **Install/repair path becomes a fallback**, not the primary
   route. Post-integration UM owns the Python update; the button remains a
   "repair / reinstall if broken" affordance for environments UM does not manage.

Open question for the core team: **does UM update individual Python packages, or
whole CCP4 components?** If the latter, `ccp4i2` rides the CCP4 release train;
if the former, `ccp4i2` can have its own cadence within a CCP4 line. Either is
compatible with Option B, but it changes how often the backend can move.

## 5. Electron side — options

The app is a per-platform packaged artifact tied to a checksum; there is no
auto-update today. Options, cheapest first:

1. **Version-gate + prompt (no new infrastructure).** The app already knows the
   backend version it needs; invert it — when a newer supported app line exists,
   show "a new CCP4i2x is available" with a download link. Manual, honest,
   zero-maintenance. Recommended **first step**.
2. **`electron-updater` against a feed** (GitHub Releases or a CCP4-hosted URL).
   The "proper" in-app updater: publishes `latest*.yml` carrying **sha512 +
   size + blockmap** per artifact (the checksum/integrity tie), with delta
   downloads. `release.yml` already produces the installers; this adds the
   update metadata and a `build.publish` target. **Blocker: signing.** Without
   macOS notarization and the Windows code-sign ([#238](https://github.com/ccp4/ccp4i2/issues/238)),
   the updater works but users hit Gatekeeper/SmartScreen. So **#238 is a
   prerequisite** for a nag-free app auto-update.
3. **UM manages the app bundle too.** One CCP4 update channel for everything.
   Conceptually tidy, but UM would have to understand per-platform Electron
   packaging and swap a large signed bundle in place — more work than letting
   `electron-updater` do what it is built for. **Not recommended.**

## 6. Recommended plan (phased, each step independently useful)

1. **Declare a stable line and relax the pin** (`==` → `>=`/range) under the
   `@ccp4/ccp4i2-api` contract. *Unblocks independent updates.*
2. **Fix the upstream metadata + dependency pins** in the CCP4 build so the
   Python update needs no repair. *Unblocks clean UM ownership.*
3. **Hand the Python update to CCP4 UM**; relegate the desktop Install/repair
   path to a fallback. *Retains the update model CCP4 users expect.*
4. **Ship the app version-gate prompt** — small, no infrastructure. *Users stop
   silently running stale apps.*
5. **Sign the installers** ([#238](https://github.com/ccp4/ccp4i2/issues/238),
   plus [macOS signing setup](macos-signing-setup.md)). *Prerequisite for (6).*
6. **Add `electron-updater` from a published feed.** *Real in-app app updates.*

Steps 1–4 are low-cost and deliver value on their own; 5–6 are the larger
investment and can wait until the app's cadence and the signing story justify
them.

## 7. Open questions for the core team

- **Coupling:** endorse decoupling to a compatible range (Option B), or keep
  lockstep for now?
- **UM granularity:** does UM update individual Python packages or whole CCP4
  components — i.e. can `ccp4i2` have its own cadence within a CCP4 line?
- **Upstream fixes:** who owns fixing the corrupt-metadata / stale-pin issues in
  the CCP4 build (§2)? They are the precondition for clean UM ownership.
- **App channel:** version-gate-and-prompt now, `electron-updater` later — or is
  there appetite to fold the app into UM despite the cost?
- **Signing:** is macOS notarization + Windows signing ([#238](https://github.com/ccp4/ccp4i2/issues/238))
  resourced? It gates the only nag-free app-update path.

---

*Straw man for discussion — concrete enough to push against, not a settled
design. Grounded in the current code: the exact-pin contract
([`ccp4i2-server-version.ts`](../client/main/ccp4i2-server-version.ts)) and the
install/repair path ([`ccp4i2-ipc.ts`](../client/main/ccp4i2-ipc.ts)).*
