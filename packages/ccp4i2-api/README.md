# `@ccp4/ccp4i2-api` / `ccp4i2-api`

Shared API contract for CCP4i2 and consumers (CCP4i2 Compounds, third-party
integrators) — auth handshake, api-fetch helpers, and request/response
types. One package, two language artifacts: TypeScript on the client side
(browser, Electron) and Python on the server side (Django middleware, DRF
authentication). Both halves agree on the canonical bearer-token format,
401 response shape, and the typed payloads carried over the authenticated
channel.

## Status

**Draft v0 — published to npm + PyPI from the in-tree workspace.** Source
of truth lives at
[`packages/ccp4i2-api/`](https://github.com/ccp4/ccp4i2/tree/django-sliced/packages/ccp4i2-api)
inside the `ccp4/ccp4i2` monorepo on the `django-sliced` branch. The npm
scope `@ccp4` and the PyPI name `ccp4i2-api` are CCP4-owned (claimed
April 2026); a dedicated `ccp4/ccp4i2-api` GitHub repo may follow later
but is not required while the monorepo hosts the source.

Versions `0.1.0`–`0.3.0` were published under the previous name
`@ccp4/ccp4i2-auth` / `ccp4i2-auth`; the package was renamed at `0.3.0`
because its scope had grown beyond auth to cover the broader API contract.
The old name is unpublished/yanked; consumers should depend on
`@ccp4/ccp4i2-api` and `ccp4i2-api` from `0.3.0` onward.

Versioning follows semver from `0.x.y` onwards. The v0 contract is
documented in
[`docs/CCP4I2_SERVICE_CONTRACT.md`](https://github.com/ccp4/ccp4i2/blob/django-sliced/docs/CCP4I2_SERVICE_CONTRACT.md);
field stability promises take effect from this version.

## Layout

| Path | Purpose |
|---|---|
| `src/` | TypeScript source. Built to `lib/` via `npm run build`. |
| `lib/` | Built TypeScript output. Generated; gitignored. |
| `dist/` | Python distribution output (`python -m build`). Generated; gitignored. Kept distinct from `lib/` so `twine upload dist/*` doesn't accidentally pick up TypeScript artefacts. |
| `ccp4i2_api/` | Python source. Installed editable via `pip install -e .`. |
| `tests/js/` | TypeScript tests (vitest, when added). |
| `tests/python/` | Python tests (pytest, when added). |

## Consumer wiring

In-monorepo consumers can reference this package by local path for fast
iteration; out-of-monorepo consumers pull the published versions.

**TypeScript** — in-monorepo (`client/package.json`):

```json
"dependencies": {
  "@ccp4/ccp4i2-api": "file:../packages/ccp4i2-api"
}
```

(Path depth varies by consumer location.) Out-of-monorepo consumers use the
published range, e.g. `"@ccp4/ccp4i2-api": "^0.3.0"`.

**Python** — in-monorepo (`Docker/server/Dockerfile`, local dev setup):

```bash
pip install -e packages/ccp4i2-api/
```

Out-of-monorepo consumers `pip install ccp4i2-api>=0.3`.

## Development

```bash
# TypeScript
cd packages/ccp4i2-api
npm install
npm run build       # produces lib/
npm run watch       # rebuilds on change

# Python
cd packages/ccp4i2-api
ccp4-python -m pip install -e .
ccp4-python -c "import ccp4i2_api; print(ccp4i2_api.__version__)"
```
