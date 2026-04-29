# `@ccp4/ccp4i2-auth` / `ccp4i2-auth`

Shared auth + api-fetch contract for CCP4i2 and consumers (CCP4i2 Compounds,
third-party integrators). One package, two language artifacts: TypeScript on
the client side (browser, Electron) and Python on the server side (Django
middleware, DRF authentication). Both halves agree on the canonical bearer-
token format and 401 response shape.

## Status

**Pre-publication scaffolding.** The package lives as a workspace inside the
`ccp4/ccp4i2` repo during the monorepo phase. It will be extracted to its
own repo (`ccp4/ccp4i2-auth`) and published to npm + PyPI as part of the
CCP4i2 / Compounds split. See
[`apps/compounds/docs/CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md`](../../apps/compounds/docs/CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md)
for the architectural rationale.

## Layout

| Path | Purpose |
|---|---|
| `src/` | TypeScript source. Built to `dist/` via `npm run build`. |
| `dist/` | Built TypeScript output. Generated; gitignored. |
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
npm run build       # produces dist/
npm run watch       # rebuilds on change

# Python
cd packages/ccp4i2-auth
ccp4-python -m pip install -e .
ccp4-python -c "import ccp4i2_auth; print(ccp4i2_auth.__version__)"
```
