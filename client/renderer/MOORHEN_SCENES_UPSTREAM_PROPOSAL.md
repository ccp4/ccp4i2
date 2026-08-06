# Moorhen scenes: a portable, validated description of a view

## Summary

A "scene" is a small YAML document that says how to render one or more
structures in Moorhen: which files to load, what to draw on them, how to colour
it, and where to put the camera. It is readable, version-controllable, and
because residue ranges are resolved against whatever structure is actually
loaded, it can be applied to different structures of the same protein, not only
the one it was written against. ccp4i2 has been using it in production. This
note asks whether Moorhen would want it upstream.

## Two uses, one mechanism

The format pays its way in two ways. The second is the reason for writing this.

### 1. Authoring by a program, including a language model

A scene is declarative. It records intent ("chain A as a ribbon coloured by
domain; the ligand as sticks; a translucent surface; the camera here"), not a
sequence of API calls. That makes it a practical thing for a program to produce,
and a language model is just a program that can take the intent in English.

This generally does not work with the alternatives. A model cannot drive a GUI,
and asking it to emit coot/Moorhen API calls invites invented method names and,
worse, one bad call leaving the viewer in a broken state. A scene avoids both,
because it is data rather than code, and because there is a validator and a
forgiving resolver between the author's output and the renderer:

- `parseScene` rejects malformed input and reports where the error is. That is a
  usable signal for a generate, check, repair loop.
- The resolver's `onMissingResidues: clamp-and-log` policy clamps a slightly
  wrong residue range to the residues that are present and logs it, instead of
  failing. The author does not have to be exactly right.

So the scene is a checked intermediate form between a fallible author and the
viewer. The deterministic resolver does the rendering; the author only has to
produce valid YAML.

We already depend on this without any model involved. ccp4i2 builds a
fragment-campaign summary scene server-side, in Python, and hands it to the same
resolver. A scene written by a model is the same pattern with a different author.

#### What works now, and what it needs

Given the grammar (one Markdown file) and a short description of a structure's
contents, its chains, ligand codes and residue ranges, which gemmi produces in
milliseconds, a model drafts a working scene today. The honest limitation is
that it has to be told those contents; it cannot guess that the ligand sits on
chain H at residue 1. The loop is therefore: read the structure, draft the
scene, validate, repair on error. The clamp-and-log resolver absorbs the small
inaccuracies that remain.

#### A demonstration

Prompt = the grammar file, plus `{chains, ligands, ranges}` from gemmi, plus
"show the binding site with the ligand as sticks and a translucent protein
surface." The model returns YAML, `parseScene` validates it, `applyScene`
renders it. Then feed it a residue range that runs off the end of the chain and
show the log line that records the clamp, in place of a crash. That is the whole
argument in one run.

### 2. Sharing and reproducibility

The same properties, declarative and validated and portable, make scenes useful
between people. A scene is a recipe for a view that can be committed next to a
paper, read in review, and re-applied to a refined or related structure. It is
not a frozen session. Moorhen's `backupSession` already covers the frozen case
well; a scene is the editable, re-applicable complement to it.

## How it would sit in Moorhen

The resolver takes the host application's fetchers as callbacks
(`SceneFileFetcher`, `SceneDictionaryFetcher`, `SceneMapFetcher`). The lifter and
resolver otherwise touch only Moorhen's own molecule, map and representation
APIs. There is no CCP4 dependency in the logic. The ccp4i2-specific parts are
narrow: a couple of file-reference kinds (`fileId`+`projectId`, `job`+`param`),
the fetcher implementations, and the routine that recovers a file id from a
loader URL.

A core version would keep the universal reference kinds (`pdb`, `url`, `bundle`,
inline `cifText`) and turn the application-specific ones (`fileId`+`projectId`,
`job`+`param`, the origin-relative `relativeUrl`) into injected extension points.
That same seam is where an application plugs in its data backend and, if it wants,
a program that authors scenes. ccp4i2 would then consume the core and register its
own reference kinds, which is a fair proof that the extension API holds.

The format is now formalised as a versioned JSON Schema contract generated from a
single Zod source — see `MOORHEN_SCENES_FOR_MOORHEN.md` for the data model, the
mapping onto Moorhen's own API (`glRefSlice`, `m2tParameters`,
`RepresentationStyles`), and the published `moorhen-scene.core.v1.json`.

## What we would bring

- The format and its grammar document.
- The lifter (Moorhen state to scene), the resolver (scene to Moorhen), the
  validator, and the `.scene.zip` bundle for self-contained sharing.
- The Scenes panel: capture, edit, open, apply, save.
- The test suite.

## Open questions

- Whether Moorhen would rather grow `backupSession` than carry a second format.
  They are different things, but that is a fair design discussion to have first.
- Who maintains it upstream. That is usually the real decision, not the code.
- Naming, and a few remaining ccp4i2-isms (proxy URLs, the id-from-URL routine)
  to clean up in the carve-out.
