# Releasing CCP4i2

One version tag produces a complete release: the **`ccp4i2` backend wheel on
PyPI** and the **desktop installers attached to a GitHub Release**. Driven by
[`.github/workflows/release.yml`](../.github/workflows/release.yml).

> This is the app/backend release. The shared `ccp4i2-api` library has its own,
> separate flow — see [`packages/ccp4i2-api/RELEASING.md`](../packages/ccp4i2-api/RELEASING.md).

## Single source of truth

The version lives in [`server/ccp4i2/__init__.py`](../server/ccp4i2/__init__.py)
(`MAJOR`/`MINOR`/`PATCH`). Everything else derives from it:

- the PyPI wheel version (setuptools `attr: ccp4i2.__version__`);
- the Electron backend-version **floor** (stamped into the build from the tag);
- the GitHub Release name.

The workflow **asserts the tag equals `ccp4i2.__version__`** and fails fast on a
mismatch, so the desktop app and PyPI can never drift apart at release time.

## Cutting a release

```bash
# 1. Bump the version
#    edit server/ccp4i2/__init__.py  ->  MAJOR/MINOR/PATCH (and __version_date__)
git add server/ccp4i2/__init__.py
git commit -m "release: v3.0.2"

# 2. Tag and push (tag must be v<the same version>)
git tag v3.0.2
git push origin django          # or your release branch
git push origin v3.0.2          # <-- this fires the Release workflow
```

The workflow then:

1. **verify-version** — tag `v3.0.2` must match `ccp4i2.__version__` (3.0.2).
2. **publish-pypi** — builds `server/` sdist+wheel and publishes `ccp4i2` to
   PyPI via OIDC.
3. **build-desktop** — builds macOS `.dmg`, Windows `.exe`, Linux `.AppImage`
   with `CCP4I2_SERVER_VERSION_FLOOR=3.0.2` stamped in.
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

## Notes & gotchas

- **A PyPI version can be published only once.** If a release fails after the
  wheel is uploaded, bump the patch and cut a new tag — you cannot re-upload the
  same version.
- **The version floor is a lower bound, not an exact pin.** The released app
  requires `ccp4i2 >= <release version>`; it will accept a newer backend.
- **macOS is not notarised.** Users still run `xattr -cr` on first install (the
  Release notes say so). Adding Apple notarisation (Developer ID cert + secrets
  in the build) is a future improvement that would remove that step.
- **Backend build needs no CCP4.** `server/ccp4i2/__init__.py` imports only the
  standard library, so the wheel builds on a bare runner.
