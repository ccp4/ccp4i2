# Proposal: a consistent parameter syntax for i2run

**Status:** proposal. Written 2026-08-26 against the brief of the same date;
revised 2026-08-27 to assume the container fixes.
Companion to [i2run.md](i2run.md), which documents the syntax as it is.

**Assumes the two container-construction fixes in
[container-construction-defects.md](container-construction-defects.md) have
landed** — the def.xml hoisting fix and deterministic `children()`. Nothing here
depends on *how* they were done, only that container trees are now free of ghost
duplicates and that sibling order is stable and follows declaration order. Where
a figure would differ without them, it is marked.

Every number below was measured against this tree; the method is in
[Appendix A](#appendix-a-how-the-numbers-were-obtained).

---

## Summary

Keep i2run's shape — `--PARAM` flags carrying tokens — and replace the three
ad-hoc addressing mechanisms with one. The changes, in order of value:

1. **Ambiguity becomes an error.** 112 parameters across 7 tasks still collide,
   and i2run silently picks one by path length. None of them has a defensible
   automatic answer.
2. **Undeclared fields become errors that name the declared ones.**
3. **One separator (`.`), one index rule (`[n]`, `[+]`), at every depth.**
4. **A declared default field per type, followed transitively** — so
   `--XYZIN beta.pdb` is a *rule*, not a special case, and it reaches into list
   items and sub-objects where today it does not.
5. **`None` unsets.** **`--` ends assignment parsing.** The `>= 0.5` heuristic goes.
6. **`--describe` and `--echo`**, generated from the container, printed in the
   input syntax.

On adopting PHIL wholesale: **no — adopt its behaviours, not its object model.**
The reason is short and structural, and it is exactly the parsimony requirement:
PHIL forbids a scope from carrying a value, so `XYZIN=beta.pdb` and
`XYZIN.annotation=native` cannot coexist in PHIL at all. Detail in
[§7](#7-on-adopting-phil-wholesale).

---

## 1. What is actually wrong

The brief lists five problems. Measuring them changed the priority order.

### 1.1 Ambiguous flags resolve silently

`ArgumentBuilder` registers each parameter under its minimum unique suffix and
adds a bare-name alias for whichever colliding parameter has the shortest path.
Nothing warns. With the container fixes in place, 112 parameters across 7 tasks
still collide, and every one of them is a *real* collision:

| Task | Params | The collision |
|---|---:|---|
| `servalcat_pipe` | 53 | `coordination.COORD16` vs `metalCoordWrapper.coordination.COORD16`, and `inputData.XYZIN` vs `metalCoordWrapper.inputData.XYZIN` |
| `prosmart_refmac` | 43 | `ADVANCED`, `OCCUPANCY` and friends in `prosmartProtein`, `prosmartNucleicAcid` *and* `libg` |
| `molrep_pipe`, `dr_mr_modelbuild_pipeline`, `adding_stats_to_mmcif_i2` | 4 each | the same name in `inputData` and `outputData` |
| `shelxeMR`, `xia2_ssx_reduce` | 2 each | ditto |

None of these has a defensible automatic answer. `prosmart_refmac --ADVANCED`
genuinely could mean any of three parameters; the user has to say which. What is
wrong is not that they collide but that i2run picks one and says nothing:

```
$ i2run prosmart_refmac --ADVANCED True     # silently -> whichever path is shortest
```

There is one honest resolution and PHIL already ships it: refuse, and name the
candidates (§7, R1).

The `inputData`/`outputData` pairs deserve a note. Once declaration order is
preserved they resolve to `inputData`, because 158 of 163 def.xml trees declare
`inputData` first and none declares it the other way. That is the right answer
almost always — but it is right by convention, not by rule, and a task that
declared its containers the other way round would silently invert. An error
costs five call sites and removes the question.

> **Historical note.** Before the container fixes this section read very
> differently: 2,412 parameters collided, 96% of them ghosts manufactured by the
> def.xml parser, and the tie-break between equally-deep candidates was
> non-deterministic run to run. Both defects, their evidence and their measured
> blast radius are written up in
> [container-construction-defects.md](container-construction-defects.md).

### 1.2 No default field below the top level

`--XYZIN beta.pdb` works: `_handle_single_value` sees a `CDataFile` and calls
`setFullPath`. But this is a branch, not a rule, and it exists only for a
directly-addressed `CDataFile`. Consequences:

- **`mergeMtz`.** `MINIMTZINLIST` items are `CMergeMiniMtz` (fields `fileName`,
  `columnTag`, `columnNames`) — a composite, not a file. A bare path therefore
  sets nothing:

  ```
  len=1 item=CMergeMiniMtz fileName.isSet=False
  ```

  The test wrote `filename=` to work around it, which created a junk attribute,
  which is the half-day the brief refers to.

- **Nested files bypass `setFullPath`.** `pdbItemList[0]/structure=/data/proj/beta.pdb`
  goes through the *list* branch of `_handle_single_value`, which ends in a plain
  `setattr`. Verified:

  ```
  baseName: '/data/proj/beta.pdb'   relPath: ''
  ```

  The whole path lands in `baseName`. The top-level branch would have split it.

### 1.3 Asymmetric addressing, for twelve parameters

Top-level lists append by flag repetition; sub-lists index with `[n]`. Across all
171 tasks that instantiate, there are **284 CLists and exactly 12 nested ones**,
in two shapes: `ENSEMBLES[n]/pdbItemList` (11 tasks) and
`splitMtz:COLUMNGROUPLIST[n]/columnList`. Maximum nesting depth is 8.

That is the whole cost of unifying the rule. Twelve parameters is not enough
structure to justify two mechanisms.

### 1.4 The heuristic

`sum(has_equals) >= len(has_equals) * 0.5` decides whether a group's tokens are
assignments. A value's interpretation depends on its neighbours. It exists to
protect SMILES strings, and there is already a better rule in the same file —
`is_fundamental_type` — which the design below promotes to *the* rule.

### 1.5 No unset, no introspection

Confirmed as stated in the brief. `--help` lists flags only.

### 1.6 One correction to the originating brief

The originating brief warns that `children()` returns `[]` for some containers and cites
`CAssembly`. `CAssembly` does not exist in `core/`, and walking all 171 tasks
(19,155 nodes) produced **no** node where `children()` was empty but fields were
declared. So the concern does not reproduce.

`dataOrder()` is still the right primitive, for a different reason: it unions
`children()` with `_data_order` and with the `@cdata_class attributes` metadata
walked over the MRO, so it reports fields that are **declared but not yet
instantiated**. "Is this field declared?" is exactly that question.

```
CMergeMiniMtz.dataOrder() -> ['fileName', 'columnTag', 'columnNames']
'filename' declared? False
```

---

## 2. The grammar

```ebnf
invocation  ::= "i2run" task_name { global_opt } { group }

group       ::= "--" address { setting } [ "--" { literal } ]

setting     ::= scalar                        (* bare: assigns the default field *)
              | address "=" scalar            (* named: assigns that field *)

address     ::= segment { "." segment }
segment     ::= NAME [ "[" index "]" ]
index       ::= INTEGER | "+"                 (* "+" appends a new item *)

scalar      ::= WORD | QUOTED | "None"
literal     ::= WORD | QUOTED                 (* never parsed as a setting *)
```

Seven rules give it meaning.

**R1 — Resolution.** A `--` flag's address is any unambiguous *suffix* of a
container path, `.`-separated. `--ADVANCED` is ambiguous in `prosmart_refmac`,
but `--libg.ADVANCED` and `--container.libg.ADVANCED` both name one parameter. If a suffix
matches more than one parameter it is an **error** listing every match. (PHIL's
rule and PHIL's diagnostic; see §7.)

**R2 — Descent.** `.` descends one level, everywhere: in flag addresses, in
setting addresses, in default-field declarations. `/` is retired.

**R3 — Indexing.** `[n]` selects item *n* of a list; `[+]` appends one. Both work
at every depth, top-level lists included:

```bash
--ENSEMBLES[1].pdbItemList[0].structure=blip.pdb
```

**R4 — Implied append.** A group whose address is a list and whose settings are
*named* applies them to a freshly appended item — `--PARAM x=1` means
`--PARAM[+] x=1`. This is exactly today's "one flag invocation, one list item",
so repeating `--PARAM` keeps working unchanged.

**R5 — The default field.** A `CData` class may declare `CLI_DEFAULT_FIELD`, an
*address* (not just a name). A bare scalar assigns it. The chain is followed
transitively until it reaches something that accepts a scalar:

| Class | `CLI_DEFAULT_FIELD` | Reach |
|---|---|---|
| `CDataFile` (and every subclass) | `fullPath` | 1,441 file params, ~145 file-typed list slots |
| `CMergeMiniMtz` | `fileName` | → then `CDataFile.fullPath` |
| `CPdbEnsembleItem` | `structure` | → then `fullPath` |
| `CEnsemble` | `pdbItemList[+].structure` | → then `fullPath` |
| `CString`/`CInt`/`CFloat`/`CBoolean` | *(itself)* | 7,298 params — 82% of all |

So all of these become legal and mean the obvious thing:

```bash
--XYZIN beta.pdb                       # XYZIN.fullPath
--MINIMTZINLIST[+] a.mtz               # .fileName.fullPath   <- fixes mergeMtz
--ENSEMBLES[+].pdbItemList[+] beta.pdb # .structure.fullPath
--ENSEMBLES beta.pdb blip.pdb          # two ensembles, one PDB each
```

There are 36 distinct list-item classes and 82 distinct parameter classes, but
after `CDataFile` and the four scalars the residue needing a bespoke entry is
about a dozen: `CEnsemble`, `CAtomRefmacSelectionGroups`, `CAltSpaceGroup`,
`CImportUnmerged`, `CAtomRefmacSelection`, `CAtomRefmacSelectionOccupancy`,
`CRunBatchRange`, `CXia2ImageSelection`, `CAsuContentSeq`, `CTLSRange`,
`CSMILESString`, `CDataReflFile`. A class with no declaration simply rejects
bare scalars, with a message listing its fields.

**R6 — Arity.** A list target takes any number of bare scalars, each appending
one item. A non-list target takes at most one; a second is an error.

**R7 — Literals, and the end of the heuristic.** A token is an assignment **iff**
the target type admits assignments. Fundamental types never do — so
`--SMILESIN "CN1CCC(=O)CC4"` and `--SMILESIN C=C` are literal, always,
independent of neighbouring tokens. Where a composite's default field is a
string that may contain `=`, the token `--` inside a group ends assignment
parsing:

```bash
--ANNOTATION -- key=value looks like an assignment but is not
```

`None` unsets (`unSet()`), at any depth: `--FREERFLAG None`,
`--XYZIN.annotation None`. A literal `"None"` string goes after `--`.

---

## 3. Errors

The two new errors are the point of the exercise; both replace a silent wrong answer.

```
$ i2run prosmart_refmac --ADVANCED True
error: '--ADVANCED' is ambiguous in prosmart_refmac. It matches:
    --prosmartProtein.ADVANCED
    --prosmartNucleicAcid.ADVANCED
    --libg.ADVANCED
  Give enough of the path to disambiguate.
```

```
$ i2run mergeMtz --MINIMTZINLIST filename=a.mtz
error: CMergeMiniMtz has no field 'filename'.
  Declared fields: fileName, columnTag, columnNames
  Did you mean 'fileName'?
  This object's default field is fileName.fullPath, so you can write:
    --MINIMTZINLIST a.mtz
```

---

## 4. Discovery

**`i2run <task> --describe [ADDRESS]`**, generated by walking the instantiated
container — `dataOrder()` for fields, `get_merged_metadata("qualifiers")` for
types, defaults, enumerators and tooltips. Both already exist; nothing new is
needed from `CData` or the def.xml.

```
$ i2run mergeMtz --describe MINIMTZINLIST
MINIMTZINLIST : CMergeMiniMtzList   (list, append with --MINIMTZINLIST or [+])
  item: CMergeMiniMtz              default field: fileName.fullPath
    fileName    CMiniMtzDataFile   default field: fullPath
      fullPath      path           (virtual: setFullPath)
      columnLabels  str            (virtual: splits the MTZ with gemmi)
      annotation    str
      ...
    columnTag   CString
    columnNames CString
```

Output is in the input syntax, so it doubles as a template. `--describe` with no
address prints the tree, eliding expert-level parameters unless `--all`.

**`i2run <task> ... --echo`** replaces `--i2run_configure`'s XML dump with the
canonical command line that reproduces the configured job — every value that
differs from its default, in the new grammar. That is what makes §6 mechanical.

---

## 5. What breaks

Measured against the 181 tests in `server/ccp4i2/tests/i2run/`:

| Form | Uses | Fate |
|---|---|---|
| bare path after a file flag | pervasive | unchanged — now R5 rather than a branch |
| `fullPath=` | 67 | unchanged |
| `columnLabels=` | 61 | unchanged (declared virtual field, see §5.1) |
| `seqFile=` | 6 | unchanged (declared virtual field) |
| repeated `--PARAM` for a list | pervasive | unchanged (R4) |
| `a/b=` | ~10 | `/` → `.`; old form warns, then errors |
| `columnList[0]/x=` | 3 sites | `columnList[0].x=` |
| `pdbItemList/structure=` (no index) | 2 | must write `[+]` or `[0]` |
| `file=` | 7 | **removed** |
| `filename=` (mergeMtz) | 2 | already broken; now a diagnosed error |
| ambiguous bare names | 112 params / 7 tasks | error where previously silent |

Two forms genuinely break — `file=` and bare `/` — plus the index-less
"apply to the last item" convenience. All three are mechanical rewrites and
`--echo` prints the replacement.

`file=` deserves its removal: it is an undeclared alias honoured by "some file
types", which is precisely the shape of the `filename=` bug. One spelling,
declared.

### 5.1 The special keys survive, but stop being special

`fullPath`, `columnLabels`, `seqFile` are not fields — they are actions
(`setFullPath`; a gemmi MTZ split that rewrites `fullPath` and `annotation`; a
sequence→ASU-XML conversion). Today they are `if` branches in the populator,
invisible to `--help` and indistinguishable from typos.

Make them **declared virtual fields**: entries in the owning class's field table
carrying a setter. Then there is one lookup rule (is this declared?), one error
path, and they appear in `--describe`. `annotation` is already a real field and
should simply drop off the "special keys" list in the docs.

---

## 6. Migration

Four phases, each shippable alone. Old forms stay accepted through phase 2 —
about one alpha cycle at the current cadence.

**Phase 0 — diagnostics only, no syntax change.** Ambiguity → error; undeclared
field → error naming the declared ones. No valid command line changes meaning;
some currently-silent-and-wrong ones start failing, which is the point. Exactly
one of the 374 `(task, flag)` pairs the i2run suite exercises uses an ambiguous
flag — `servalcat_pipe --XYZIN` — so expect a single site to need the long form.
Run `server/run_i2run_baseline.sh` before and after and diff `results.xml`.

**Phase 1 — `--describe` and `--echo`.** Read-only. This is what makes phases 2
and 3 checkable rather than hopeful.

**Phase 2 — new grammar accepted alongside old.** `.`, `[n]`/`[+]` everywhere,
`None`, `--`, the default-field chain. `/`, `file=` and index-less list
navigation warn, naming their replacement. Migrate `tests/i2run/` in one commit:
the suite is the de-facto specification, so it should be the first consumer of
the new one.

**Phase 3 — one release later.** The deprecated forms error. The `>= 0.5`
heuristic is deleted.

There is no installed base beyond this repo's test suite, so the window can be
short. What there *is* is a documentation debt: [i2run.md](i2run.md) currently
lists four keys "as though that were the set", and Phase 1 is what makes the
real answer printable instead of prose.

---

## 7. On adopting PHIL wholesale

Worth taking seriously — ten CCP4i2 tasks already have PHIL underneath, and
`utils/phil_to_cdata.py` already converts one way. I tested the other direction.

**What PHIL gets right, verified against `ccp4-20260702`:**

- Unique-suffix resolution, with **ambiguity as a hard error** naming the
  candidates — the fix from §3, already specified and battle-tested:
  ```
  Sorry: Ambiguous parameter definition: cycles = 20
  Best matches:  refine.cycles   scale.cycles
  ```
- `None` unsets.
- `.show(attributes_level=2)` is `--describe`, in the input syntax.
- `fetch_diff` is `--echo`: what you set, minus the defaults, ready to re-run.
- Repeated `{ }` blocks build `.multiple` scopes.

**And the packaging objection is weaker than it looks.** `libtbx.phil` is pure
Python: it imports on stock CPython 3.13 in the CCP4-free venv with `six` as its
only extra dependency, pulling in 13 modules totalling 281 KB and loading no
compiled extensions. It would not breach the slim-server boundary.

**Two things stop it anyway.**

*First, and decisively: PHIL cannot express the parsimony we want.* PHIL enforces
a strict scope/definition dichotomy — a scope has children and no value; a
definition has a value and no children. Both directions are refused:

```
xyzin { baseName=... }  +  xyzin=beta.pdb
  -> RuntimeError: Incompatible parameter objects: scope "xyzin" vs. definition "xyzin"

xyzin=None (.type=path)  +  xyzin.annotation="native"
  -> RuntimeError: Incompatible parameter objects: definition "xyzin" vs. scope "xyzin"
```

A CCP4i2 file parameter is both at once: it has a natural scalar (the path) *and*
sub-fields (`columnLabels`, `annotation`, `selection`). `--XYZIN beta.pdb`
alongside `--XYZIN.annotation native` is the headline requirement, and it is the
one thing PHIL's grammar structurally forbids. R5 is our own invention because it
has to be.

*Second, the CLI layer is not free even if we take it.*
`argument_interpreter.process_arg()` mis-splits multi-definition blocks for
`.multiple` scopes — on the ENSEMBLES shape specifically. Two blocks became four
scopes:

```
ai.process('ensemble { pdb=beta.pdb \n identity=0.8 }')   # via argument_interpreter
  -> 4 ensembles:  (beta.pdb, 0.9)  (None, 0.8)  (blip.pdb, 0.9)  (None, 0.7)

libtbx.phil.parse('ensemble { pdb=beta.pdb \n identity=0.8 }')   # direct
  -> 2 ensembles:  (beta.pdb, 0.8)  (blip.pdb, 0.7)     <- correct
```

So a PHIL-based i2run would have to route block arguments through `parse()` and
reserve `argument_interpreter` for bare `suffix=value` — i.e. write the
command-line layer regardless.

Behind both sits a third: PHIL's model is scalars and scopes, and CCP4i2
parameters are objects with behaviour. `setFullPath`'s project-relative
splitting, the gemmi MTZ split, the ASU conversion — none has a PHIL home, and
all would end up in a post-processing layer. Our own `Phil2CData` already shows
the seam: it maps `.multiple` to a *qualifier*, not to a `CList`.

**Recommendation.** Take PHIL's behaviours — suffix resolution, ambiguity as an
error with candidates, `None`, describe-in-input-syntax, echo-the-diff — and
implement them against `CData`. Where a task *is* PHIL-backed, the surface then
reads like the PHIL the underlying tool takes, which is the compatibility that
actually matters. Reject the object model, for the reason above.

---

## Appendix A: how the numbers were obtained

All figures come from instantiating every task in `core/tasks.py` in
`server/.venv` (stock CPython 3.13, no CCP4) and walking the containers. 171 of
173 tasks instantiate; `KeywordExtractor` introspects 168.

| Figure | Value |
|---|---|
| Tasks instantiated / introspected | 171 / 168 |
| Leaf parameters seen by argparse | 8,926 |
| ...with a colliding simple name | **112, in 7 tasks** (2,412 in 16 before the container fixes) |
| ...that are `CString`/`CInt`/`CFloat`/`CBoolean` | 7,298 (82%) |
| Nodes reached descending into files and list items | 19,155 |
| `CDataFile` instances | 1,441 |
| `CList` parameters | 284 |
| ...nested inside another parameter | **12**, in 2 shapes |
| Maximum nesting depth | 8 |
| Distinct top-level parameter classes | 82 |
| Distinct `CList` item classes | 36 |
| `(task, flag)` pairs exercised by `tests/i2run/` | 374 |
| ...using a flag that is ambiguous for that task | 1 (`servalcat_pipe --XYZIN`) |
| Nodes where `children()` was empty but fields declared | 0 |

PHIL behaviour was tested with `ccp4-20260702/bin/ccp4-python`; the pure-Python
claim by appending CCP4's `site-packages` to `sys.path` in `server/.venv` and
importing `libtbx.phil` and `libtbx.phil.command_line`.

Test-suite counts are token censuses over `server/ccp4i2/tests/i2run/*.py`
(101 files, 181 tests).
