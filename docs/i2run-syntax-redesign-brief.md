# Brief: a consistent parameter syntax for i2run

**Status:** brief only. No design, no implementation. Written 2026-08-26 out of
the `mergeMtz` investigation, which cost half a day largely because nobody could
say what a valid parameter key was.

## The task

Propose a parameter syntax for `i2run` that is consistent, unambiguous,
discoverable and expressible against CCP4i2's `CData` object model. Say what it
would break, and what a migration would look like. Reference:
[docs/i2run.md](i2run.md) is the current documentation;
`server/ccp4i2/cli/i2run/i2run_components.py` is the parser.

## Why now

Three defects this week traced back to the syntax rather than to any task:

- `mergeMtz` silently merged nothing for years. The test passed `filename=`
  where the field is `fileName`; `CData` accepts dynamic attributes, so the
  parser created a junk attribute and the real field stayed unset. Fixed
  2026-08-26 by refusing undeclared fields — which is a patch, not a design.
- The documentation lists four keys as though that were the set. The real rule
  is "any field the target declares", which appears nowhere.
- There is no way to ask what a parameter accepts. Since the fix, the *error
  message* is the discovery mechanism.

## What exists today

Three addressing mechanisms, acquired separately:

| Form | Meaning | Example |
|---|---|---|
| repetition | append an item to a top-level list | `--ENSEMBLES … --ENSEMBLES …` |
| `[n]` | index into a list *inside* a parameter | `pdbItemList[1]/structure=…` |
| `/` | descend into a sub-object | `imageFile/baseName=…` |

Plus special keys handled before the general path: `fullPath=`, `file=`,
`columnLabels=`, `annotation=`, `seqFile=` (which converts a sequence file into
ASU XML).

Verified reachable: `--ENSEMBLES use=True pdbItemList[0]/structure=A.pdb
pdbItemList[1]/structure=B.pdb` sets the second PDB of that ensemble correctly.
So the syntax is *expressive*; the problems are elsewhere.

## Known problems, with evidence

1. **Asymmetric addressing.** Top-level lists are append-only (position implied
   by repetition order); sub-lists are random-access by `[n]`. "Set one field of
   the third ensemble" requires constructing the first two.
2. **Two spellings, two semantics.** `--XYZIN fullPath=X` goes through file
   handling; nested `pdbItemList[0]/structure=X` assigns the whole path to
   `baseName` directly, bypassing the project-relative splitting that
   `setFullPath()` exists to do.
3. **A heuristic where a grammar should be.** Whether tokens are `key=value`
   pairs or plain values is decided by
   `sum(has_equals) >= len(has_equals) * 0.5`, which exists because SMILES
   strings contain `=`. A value's interpretation depends on its neighbours.
4. **No unset.** Nothing expresses "like that job, but without the FREERFLAG".
5. **No introspection.** `--help` lists top-level parameters only.

## Constraints a proposal must respect

- **`CData` is dynamic.** Attributes can be created by assignment, so "does this
  field exist" must be asked explicitly — `[c.objectName() for c in obj.children()]`
  is the available answer, and it returns `[]` for some containers (e.g.
  `CAssembly`), which any check must tolerate.
- **Lists are typed by `makeItem()`**, and items may be containers with their own
  lists.
- **The def.xml is the source of truth** for what a task's parameters are, and
  `core/tasks.py` can enumerate all 173 tasks cheaply — so a `--describe` is
  buildable from what already exists.
- **The GUI and the REST API populate the same containers** by other routes.
  Whatever is proposed should not need `CData` itself to change.
- **Existing scripts and the whole i2run test suite use the current syntax.** A
  proposal that cannot be adopted incrementally is not adoptable.

## What a good answer contains

1. A written grammar — one page, unambiguous, no heuristics.
2. A single addressing rule covering top-level lists, nested lists and
   sub-objects, or an explicit argument for why two are needed.
3. A discovery command (`i2run <task> --describe <PARAM>` or similar) generated
   from the def.xml.
4. A statement of what breaks, and a migration path — including whether the old
   forms stay accepted, and for how long.
5. Whether the special keys survive as sugar, or become ordinary fields.

## Non-goals

Changing `CData`, changing the def.xml format, or redesigning the REST API. The
brief is about the command line only.

## Where to look

- `server/ccp4i2/cli/i2run/i2run_components.py` — `PluginPopulator`, the parser
- `server/ccp4i2/cli/i2run/CCP4i2RunnerBase.py` — argument construction
- `docs/i2run.md` — current reference, corrected 2026-08-26
- `server/ccp4i2/tests/i2run/` — ~160 tests, the de-facto specification of the
  syntax in use
- `manage.py i2run <task> --i2run_configure …` configures a job without running
  it, and prints the populated container as XML. That is the cheap way to test
  whether a form expresses what you meant.
