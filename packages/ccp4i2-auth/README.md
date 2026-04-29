# `@ccp4/ccp4i2-auth` / `ccp4i2-auth`

Shared auth + api-fetch contract for CCP4i2 and consumers (CCP4i2 Compounds,
third-party integrators). One package, two language artifacts: TypeScript on
the client side (browser, Electron) and Python on the server side (Django
middleware, DRF authentication). Both halves agree on the canonical bearer-
token format and 401 response shape.

## Status

**Draft v0 — published to npm + PyPI from the in-tree workspace.** The
package will move to its own repo (`ccp4/ccp4i2-auth`) once the GitHub-side
ask lands; until then the source of truth is
[`packages/ccp4i2-auth/`](https://github.com/ccp4/ccp4i2/tree/django-sliced/packages/ccp4i2-auth)
inside the `ccp4/ccp4i2` monorepo on the `django-sliced` branch.

Versioning follows semver from `0.x.y` onwards. The v0 contract is
documented in
[`apps/compounds/docs/CCP4I2_SERVICE_CONTRACT.md`](https://github.com/ccp4/ccp4i2/blob/django-sliced/apps/compounds/docs/CCP4I2_SERVICE_CONTRACT.md);
field stability promises take effect from this version.

## Layout

| Path | Purpose |
|---|---|
| `src/` | TypeScript source. Built to `lib/` via `npm run build`. |
| `lib/` | Built TypeScript output. Generated; gitignored. |
| `dist/` | Python distribution output (`python -m build`). Generated; gitignored. Kept distinct from `lib/` so `twine upload dist/*` doesn't accidentally pick up TypeScript artefacts. |
| `ccp4i2_auth/` | Python source. Installed editable via `pip install -e .`. |
| `tests/js/` | TypeScript tests (vitest, when added). |
| `tests/python/` | Python tests (pytest, when added). |

## Consumer wiring

During the monorepo phase, consumers reference this package via local paths.

**TypeScript** (`client/package.json`, `apps/compounds/frontend/package.json`):

```json
"dependencies": {
  "@ccp4/ccp4i2-auth": "file:../../packages/ccp4i2-auth"
}
```

(Path depth varies by consumer location — see each consumer's `package.json`.)

**Python** (`Docker/server/Dockerfile`, local dev setup):

```bash
pip install -e packages/ccp4i2-auth/
```

At repo split, consumers swap `file:` to a published version (`^0.1.0` for
TS, `>=0.1` for Python). The package's exports, build artifacts, and
contract are unchanged — only the resolution mechanism switches.

## Development

```bash
# TypeScript
cd packages/ccp4i2-auth
npm install
npm run build       # produces lib/
npm run watch       # rebuilds on change

# Python
cd packages/ccp4i2-auth
ccp4-python -m pip install -e .
ccp4-python -c "import ccp4i2_auth; print(ccp4i2_auth.__version__)"
```
