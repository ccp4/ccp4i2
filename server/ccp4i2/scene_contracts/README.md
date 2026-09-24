# Moorhen scene contracts — server-readable mirror

The LLM-facing contracts for the Moorhen scene format. Two are **generated
artifacts mirrored from the frontend**, where the Zod source of truth lives
(`client/renderer/lib/scene/`); the third is a hand-written authoring guide (see
the table below). They live here so a slim, frontend-free Django deployment —
and any external adopter — can read them from the installed `ccp4i2` package.
`server/ccp4i2/**/*` ships as package data, so:

```python
from importlib.resources import files
text = (files("ccp4i2") / "scene_contracts" / "moorhen-scene.system-prompt.v1.md").read_text()
```

| File | Role | Generated? |
|------|------|-----------|
| `moorhen-scene.system-prompt.v1.md` | the static LLM **system** message for scene authoring (Materia's `nlp_scene` endpoint) | yes — from Zod |
| `moorhen-scene.structured.v1.json` | strict **OpenAI Structured Outputs** profile (production `json_schema` constraint) | yes — from Zod |
| `moorhen-scene-authoring.v1.md` | the **authoring guide**: pitfalls, worked examples, and the empirical findings behind them | **no — hand-written** |

## Using these with any model

The three files are provider-neutral and ship as `ccp4i2` package data, so an
adopter can drive whichever model they prefer:

```python
from importlib.resources import files
contracts = files("ccp4i2") / "scene_contracts"
system = (
    (contracts / "moorhen-scene.system-prompt.v1.md").read_text()
    + "\n\n"
    + (contracts / "moorhen-scene-authoring.v1.md").read_text()
)
# Constrain the response with moorhen-scene.structured.v1.json wherever the
# provider supports JSON-Schema-constrained output.
```

The two generated files say what is **well-formed**; the authoring guide says
what is **correct in practice**. A scene can pass the schema and still render
nothing — an unparenthesised `//*/LIG`, a comma list of residue numbers, a
domain range that lies outside the modelled residues. Those traps are only in
the authoring guide, so ship it alongside the grammar rather than instead of it.

## Do not edit the generated files by hand

`moorhen-scene-authoring.v1.md` **is** hand-maintained — it is prose, not
derived from the schema, so edit it directly and keep it honest against the
validator (its scene examples are checked by the frontend test suite). The
other two are byte-equality-checked against the canonical frontend copies by
`client/renderer/__tests__/scene-schema.test.ts`. To regenerate after an
intentional schema/prompt change, run (from `client/`):

```bash
UPDATE_SCHEMA=1 npx vitest run renderer/__tests__/scene-schema.test.ts
```

That rewrites both the `client/renderer/lib/scene/` originals and these mirrors
from the one Zod source, so they can never drift. See
`client/renderer/MOORHEN_SCENES_SCHEMA_V1_DESIGN.md` §12.
