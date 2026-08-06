# Moorhen scene contracts — server-readable mirror

These files are **generated artifacts mirrored from the frontend**, where the Zod
source of truth lives (`client/renderer/lib/scene/`). They are duplicated here so
a slim, frontend-free Django deployment can read them from the installed `ccp4i2`
package — `server/ccp4i2/**/*` ships as package data, so:

```python
from importlib.resources import files
text = (files("ccp4i2") / "scene_contracts" / "moorhen-scene.system-prompt.v1.md").read_text()
```

| File | Role |
|------|------|
| `moorhen-scene.system-prompt.v1.md` | the static LLM **system** message for scene authoring (Materia's `nlp_scene` endpoint) |
| `moorhen-scene.structured.v1.json` | strict **OpenAI Structured Outputs** profile (production `json_schema` constraint) |

## Do not edit by hand

They are byte-equality-checked against the canonical frontend copies by
`client/renderer/__tests__/scene-schema.test.ts`. To regenerate after an
intentional schema/prompt change, run (from `client/`):

```bash
UPDATE_SCHEMA=1 npx vitest run renderer/__tests__/scene-schema.test.ts
```

That rewrites both the `client/renderer/lib/scene/` originals and these mirrors
from the one Zod source, so they can never drift. See
`client/renderer/MOORHEN_SCENES_SCHEMA_V1_DESIGN.md` §12.
