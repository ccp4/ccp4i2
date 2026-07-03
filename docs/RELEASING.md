# Releasing CCP4i2

One version tag produces a complete release: the **`ccp4i2` backend wheel on
PyPI** and the **desktop installers attached to a GitHub Release**. Driven by
[`.github/workflows/release.yml`](../.github/workflows/release.yml).

> This is the app/backend release. The shared `ccp4i2-api` library has its own,
> separate flow — see [`packages/ccp4i2-api/RELEASING.md`](../packages/ccp4i2-api/RELEASING.md).

## Single source of truth

The version lives in [`server/ccp4i2/__init__.py`](../server/ccp4i2/__init__.py)
(`MAJOR`/`MINOR`/`PATCH` + `PRERELEASE`). Everything else derives from it:

- the PyPI wheel version (setuptools `attr: ccp4i2.__version__`);
- the **exact** backend version the Electron app pins to (stamped into the build
  from the tag);
- the GitHub Release name.

The workflow **asserts the tag equals `ccp4i2.__version__`** and fails fast on a
mismatch, so the desktop app and PyPI can never drift apart at release time.

## Cutting a release

```bash
# 1. Bump the version
#    edit server/ccp4i2/__init__.py -> MAJOR/MINOR/PATCH + PRERELEASE (+ date).
#    During the alpha, bump PRERELEASE: "a1" -> "a2" -> ... (or "" for a final).
git add server/ccp4i2/__init__.py
git commit -m "release: v3.1.0a1"

# 2. Tag and push (tag must be v<the same version>, e.g. v3.1.0a1)
git tag v3.1.0a1
git push origin django          # or your release branch
git push origin v3.1.0a1        # <-- this fires the Release workflow
```

The workflow then:

1. **verify-version** — tag `v3.1.0a1` must match `ccp4i2.__version__` (3.1.0a1);
   a pre-release version (`a`/`b`/`rc`) is flagged so the Release is marked
   **Pre-release**.
2. **publish-pypi** — builds `server/` sdist+wheel and publishes `ccp4i2` to
   PyPI via OIDC.
3. **build-desktop** — builds macOS `.dmg`, Windows `.exe`, Linux `.AppImage`
   with `CCP4I2_SERVER_VERSION_FLOOR=3.1.0a1` stamped in (the env var name is
   historical; it now carries the **exact** required version).
4. **release** — creates the GitHub Release and attaches all three installers as
   public, permanent, no-login assets.

A tag like `v3.0.2-rc1` is published as a **pre-release** (GitHub marks it as
such automatically; PyPI treats `rcN` as a pre-release version).

## One-time setup (before the first `v*` tag)

PyPI Trusted Publishers are **per project**. The existing `ccp4i2-api` publisher
does **not** cover `ccp4i2`. Add one for the main package:

- <https://pypi.org/manage/project/ccp4i2/settings/publishing/> →
  - **Owner:** `ccp4`
  - **Repository:** `ccp4i2`
  - **Workflow:** `release.yml`
  - **Environment:** `pypi` (reuses the existing GitHub environment — no new
    secret, no token)

No `PYPI_TOKEN` is required or used — publishing is OIDC trusted publishing, the
same as `ccp4i2-api`.

## macOS signing & notarisation (optional but recommended)

The release **degrades gracefully**: with no Apple secrets configured, the macOS
build ships **unsigned** (users clear quarantine with `xattr -cr`, exactly as
today). Add the secrets below and the Release workflow **signs and notarises**
the `.dmg` automatically — no `xattr` needed by users, no "app is damaged"
scare.

You need **two** Apple credentials — the API key notarises but does **not**
sign:

1. **Developer ID Application certificate** (`.p12`) — for code signing.
2. **App Store Connect API key** (`.p8` + Key ID + Issuer ID) — for notarising.

Add these as repository (or org) **secrets**:

| Secret | What it is |
|---|---|
| `MAC_CSC_LINK` | base64 of the Developer ID Application `.p12` (`base64 -i cert.p12 \| pbcopy`) |
| `MAC_CSC_KEY_PASSWORD` | the `.p12` export password |
| `APPLE_API_KEY_P8` | the **contents** of the App Store Connect `AuthKey_XXXX.p8` |
| `APPLE_API_KEY_ID` | the API Key ID |
| `APPLE_API_ISSUER` | the API Issuer ID |

The `build-desktop` job detects these: when `MAC_CSC_LINK` **and**
`APPLE_API_KEY_P8` are present it signs (`hardenedRuntime` + entitlements from
[`client/assets/entitlements.mac.plist`](../client/assets/entitlements.mac.plist))
and runs `electron-builder … -c.mac.notarize=true`; otherwise it builds unsigned
with `CSC_IDENTITY_AUTO_DISCOVERY=false`.

> **Validate with a pre-release tag first**, e.g. `v3.0.2-rc1` — signing +
> notarisation can only be truly confirmed by building on CI and opening the
> resulting `.dmg` on a clean Mac (no `xattr`). Once it's green, drop the `xattr`
> note from `docs/give-it-a-try.md`.

## Notes & gotchas

- **A PyPI version can be published only once.** If a release fails after the
  wheel is uploaded, bump the pre-release counter (`a1`→`a2`) and cut a new tag —
  you cannot re-upload the same version.
- **Pre-release versioning + exact pin (alpha discipline).** The version carries
  a PEP 440 pre-release segment (`PRERELEASE = "a1"` in `server/ccp4i2/__init__.py`).
  Pre-releases are invisible to a plain `pip install ccp4i2`, and the desktop app
  pins the backend **exactly** (`ccp4i2==<version>`, see
  [`client/main/ccp4i2-server-version.ts`](../client/main/ccp4i2-server-version.ts)),
  so an alpha app and backend are strictly bound and can't mix with the escaped
  `3.0.x` line. A pre-release tag (`vX.Y.Za1`) is auto-marked **Pre-release** on
  GitHub. To declare a stable line, set `PRERELEASE = ""`, tag `vX.Y.Z`, and (if
  desired) relax the app back to a floor/compatible-range check.
- **Optionally yank the superseded finals.** `3.0.0`/`3.0.1` can be *yanked* on
  PyPI (PEP 592) so normal resolution won't pick them, without deleting them
  (exact pins still work). Alpha apps are unaffected (they pin exact anyway).
- **Backend build needs no CCP4.** `server/ccp4i2/__init__.py` imports only the
  standard library, so the wheel builds on a bare runner.
