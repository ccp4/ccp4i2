# CData: what the magic buys, and what it would take to stop paying for it

**Status:** design note, 2026-08-28. Written to be argued with. Every number in
it was measured on this tree; the scripts are named so they can be re-run.

**Related:** [container-construction-defects.md](container-construction-defects.md)
recorded the defects that prompted this and holds the conformance harness
design. This note is the destination that work is heading towards.

---

## The proposition

Paul's position, and it is the right instinct: a CData object should carry as
close to no magic as possible. Children and values might collapse into ordinary
Python attributes, or dataclasses.

Two objections have kept it from happening:

1. **Two-way navigation.** A leaf walks *up* the tree to find out where it
   lives. Plain attributes cannot do that.
2. **Construction from def.xml.** The shapes are not known until a file is
   read, and dataclasses are declared in source.

Both are real. Both turn out to fall on one side of a distinction the current
model does not make.

---

## What a leaf costs today

A `CInt` leaf carries **16 instance attributes**:

    _children              _children_by_name      _child_storage
    _lock                  _event_handlers        _properties
    _sigmgr                _qualifiers            _default_values
    _value_states          _parent_ref            _name
    _state                 _hierarchy_initialized _skip_validation
    _value

That is **~352 bytes against 28 for a plain int**, and **86%** of a task tree's
objects are leaves --- 335 of 391 in `phaser_pipeline`. Across the registry:
about 18,000 leaves, each with five dicts and an `RLock`, none of which they
use, because a leaf never has children.

The machinery is paid overwhelmingly by the objects that do not need it.

---

## The distinction the model does not make

Measured over every task in the registry (`/tmp` script, reproduced in the
appendix): 23 classes own children, 289 instances; 23 classes are always
leaves, 1,936 instances. `CContainer` appears in **both** lists --- 189 with
children, 168 empty --- so even "is a container" is not a stable property of a
class.

The revealing measurement is what those 23 owner classes look like:

| | |
|---|---:|
| classes whose children are **always the same set** | **22** |
| classes whose children **vary between instances** | **1** (`CContainer`, 70 shapes) |

`CObsDataFile` is always exactly those 7 fields. `CPdbDataFile` always 8.
`CI2XmlHeader` always 13. Only `CContainer` varies, because its contents come
from def.xml at runtime.

So one word --- "children" --- is doing two jobs:

**Composition.** A `CObsDataFile` *is made of* `baseName`, `relPath`,
`project`, `annotation`, `dbFileId`, `contentFlag`, `subType`. The *shape* is
fixed and declared in Python source, identical for every instance. This is a
**record** --- though see the note on qualifiers below: what each field *means*
can still be adjusted by the def.xml that uses it.

**Containment.** A task's `inputData` *contains* `XYZIN`. Dynamic, arbitrary
depth, addressable by path, defined by data outside the program. This is a
**tree**.

Conflating them is why the model needs machinery general enough for the tree
on every object that only needs the record.

---

## The fixed field list is itself the weirdness

A record class's fields look like attributes and are not. `CDataFile` declares
them as entries in a decorator argument:

```python
@cdata_class(
    attributes={
        "baseName":    attribute(AttributeType.CUSTOM, custom_class="CFilePath"),
        "relPath":     attribute(AttributeType.CUSTOM, custom_class="CFilePath"),
        "dbFileId":    attribute(AttributeType.CUSTOM, custom_class="CUUID"),
        "annotation":  attribute(AttributeType.CUSTOM, custom_class="CString"),
        "subType":     attribute(AttributeType.CUSTOM, custom_class="CInt"),
        "contentFlag": attribute(AttributeType.CUSTOM, custom_class="CInt"),
    },
    qualifiers={...})
class CDataFile(CData):
    """Attributes are automatically created from embedded metadata..."""
```

Nothing in the class body names a field. The docstring lists them in prose,
which is the tell: the source cannot say what the object is made of, so a
human wrote it down twice.

And the declaration is not even in the same file as the class. `CObsDataFile`
is `class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile)` --- data model in a
stub, methods in the derived class:

| | |
|---|---:|
| classes decorated `@cdata_class` | 229 |
| `*Stub` classes | 210 |
| classes inheriting a stub | 253 |
| lines of auto-generated stubs | 21,472 |

`CObsDataFileStub` is at line 5,742 of an auto-generated
`core/cdata_stubs/CCP4XtalData.py` marked **DO NOT EDIT**, and carries four
separate kinds of metadata --- `attributes`, `qualifiers` (a list of names),
`qualifiers_definition`, `content_qualifiers`.

So answering "what is a `CObsDataFile` made of" means finding the derived class
in one file, its stub in another, and merging decorator metadata up an MRO that
runs through `CMiniMtzDataFileStub`, `CDataFileStub`, `CDataStub`.

**This is why `dataOrder()` merges four sources and `value_dict_for_object`
tries five strategies.** When a field list is assembled from decorator metadata
across two parallel class hierarchies, there is no single place to ask --- so
every consumer invents its own way of asking, and they agree by luck.

A dataclass removes exactly this. Fields are visible in the body of the class
that has them; `dataclasses.fields(cls)` is the one introspection point;
`field(metadata=...)` carries what the `attribute(...)` entries and the
qualifier dicts carry now. The 21,472 lines of generated stubs exist to
work around Python not being told about the fields --- which is the problem
dataclasses were added to the language to solve.

---

## Why the model went into metadata, and what became of that goal

The decorator was not an accident or an operational necessity. The intent, in
migrating away from QObject, was to **capture each class's data model somewhere
that is not Python code**, so that it could later be harvested to inform
frontend items. Discipline rather than mechanism.

The ambition is right, and it deserves recording because the rest of this note
reads as criticism otherwise. Two things happened to it.

**The harvest never took place.** Nothing in the client reads
`get_merged_metadata`, the `attributes` dicts, or the generated stubs. What the
frontend actually consumes is the **live container JSON** --- `{_class,
_value, _qualifiers}` at every node --- fetched over REST and rendered by
`GenericInterface`. So the goal was met by a different route: the *instance*
turned out to be a perfectly good harvest of the *model*, and it arrives
already, per task, at runtime.

**And a better mechanism for the original goal exists in this repository.**
The Moorhen scene format declares its schema in Zod and *generates* JSON Schema
from it (`client/renderer/lib/scene/`). Author in a language the type checker
and the IDE understand; emit the language-neutral artefact for whoever needs
it. That is the same ambition without putting the model beyond Python's reach
--- and `dataclasses.fields()` plus `field(metadata=...)` is directly
harvestable in exactly that way, at build time, to JSON Schema or anything
else.

So the cost of keeping the model out of Python --- 21,472 generated lines, no
single introspection point, consumers each reconstructing the field list --- was
paid, and the benefit it was paid for has been taken by another route. That is
the strongest reason to think the model can come back into the class body
without losing anything that was wanted.

---

## Objection 1: two-way navigation

Real, and specific. `_find_plugin_parent()` is how a **file value answers
questions about where it lives**:

- am I in database-aware mode (`cdata_file.py:203`)
- what is my database handler (`:285`)
- what is my project root, or my job's working directory, for resolving
  `relPath` (`:715`)

That is why `container.inputData.XYZIN = "/some/path"` can create a File record
in the right project: the leaf walks up and asks. Upward navigation appears
about 350 times across its spellings.

**But the walk is already known not to work.** `_find_plugin_parent` opens with
an escape hatch:

```python
if hasattr(self, '_temp_plugin_ref'):
    return self._temp_plugin_ref
```

and the gleaner bolts that attribute on before use and strips it off afterwards
(`db/async_db_handler.py:676, 697`). A hierarchy walk with a manual override
for the cases it cannot reach is already passing context explicitly --- by
monkey-patching, rather than by design.

**What the leaf needs is not its parent. It is its context**: project root,
working directory, whether a database is attached, which handler. One object.
The parent chain is a *lookup mechanism* for it, not the requirement.

Hand the context down instead of walking up and:

- no depth limit, no weak reference to a parent, no `_temp_plugin_ref`
- **children need not be discoverable from below at all** --- which is the
  property that blocks plain attributes
- the requirement lands on **file values**, one level, not on all 18,000 leaves.
  `baseName` and `annotation` need no context whatsoever.

---

## Objection 2: construction from def.xml

The resolution is to stop trying. **Nothing generates a class from def.xml.**

def.xml describes *containment*, which is data; a container built from it stays
a dynamic mapping. Dataclasses are for the 22 *composition* classes, which are
declared in Python source and never vary.

A sketch, built from the real `freerflag.def.xml`:

```python
@dataclass
class Context:
    """What a file value needs to resolve itself. Handed down, never walked up."""
    project_root: str = ''
    work_directory: str = ''
    db_handler: Any = None

@dataclass
class DataFile:
    """The 7 fields every CDataFile has --- fixed for every instance."""
    baseName: str = ''
    relPath: str = ''
    project: str = ''
    annotation: str = ''
    dbFileId: Optional[str] = None
    contentFlag: int = 0
    subType: int = 0
    _context: Optional[Context] = field(default=None, repr=False, compare=False)

    def full_path(self):
        root = self._context.project_root if self._context else ''
        return os.path.join(root, self.relPath, self.baseName)

class Container(dict):
    """Dynamic, ordered, addressable. The shape is data; pretending otherwise
    buys nothing."""
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)
```

Parsing the actual def.xml and walking it gives:

    sections      : ['inputData', 'outputData', 'controlParameters']
    dotted access : DataFile
    addressing    : ['container.inputData.F_SIGF', 'container.inputData.FREERFLAG',
                     'container.outputData.FREEROUT', 'container.controlParameters.CUTRESOLUTION']
    full_path     : /projects/demo/CCP4_IMPORTED_FILES/gamma.mtz

Dotted access, path addressing, and a file value resolving its own full path,
with no generated classes and no upward navigation. The context was handed to
the value when it was made.

---

## What the sketch does not yet carry

Stated plainly, because these are the work:

**Value state.** The dependency for step 2, and the one that has to be
*decided* rather than moved. It has its own note ---
[value-state-design.md](value-state-design.md) --- because measuring it turned
up that the three-state model is not what the system runs on: a `CString` is
`EXPLICITLY_SET` from birth, and `isSet()` reports non-emptiness for strings
but state for ints, so the predicate 719 call sites depend on means different
things by type. A fourth state is also missing --- `DERIVED`, for what the
system worked out rather than what a user chose --- which is why `contentFlag`,
`cell` and output paths are recorded as though somebody typed them.

**Qualifiers.** Currently a per-instance dict, and the obvious simplification
--- put them on the class --- does not work, because **def.xml overrides them
per field**. Measured across the registry: 18 of 44 classes carry more than one
qualifier set, `CString` alone in **109 distinct** variants; and the same class
with the same field name still differs, `CFilePath.baseName` in 17 variants,
`CObsDataFile.F_SIGF` in 3. A `baseName` under one task's file input is not
qualified the same way as under another's.

They are nonetheless not *instance* data. Two instantiations of the same task
agree at every path, in **all 170** tasks checked --- so qualifiers are a
property of **(task, path)**, and belong in a per-task schema computed once
when the def.xml is read, layering its overrides over the class defaults.

That is a larger structure than a class attribute, and it is the honest cost of
the override facility. It is still ~21,000 schema entries held once per task
rather than once per object, and it removes a dict from every one of the 1,936
leaf instances.

**Assignment coercion.** `container.inputData.XYZIN = "/some/path"` setting a
file from a string --- 28 uses, against 3,711 dotted reads. Making that work
means something intercepts assignment.

An earlier draft said "the container intercepts assignment --- magic, but in one
place rather than eighteen thousand". That is wrong twice, and one concrete path
shows both. Setting the annotation of the third model in the second assembly of
a `phaser_MR` input:

```python
container.inputData.ENSEMBLES[1].pdbItemList[2].structure.annotation = "..."
```

**It is not one place, because the interception happens on a `CPdbDataFile`.**
`structure` is not a `CContainer`; neither is the `CEnsemble` above it nor the
`CPdbEnsembleItem` between them. Every type that has `CData` children has to
intercept, so the magic lives on `CContainer`, on `CList`, and on each of the
composed classes --- three of them on this single path. That is still a real
reduction, from every object to the 14% that are parents, but it is ~22 classes,
not one.

**And leaves cannot shed it either, because the leaf owns its own state.**
`_value_states` is written inside `__setattr__` *on the object being assigned
to*, so a leaf must intercept to record that it was set:

    frac.value = 0.3          # no parent involved anywhere
    frac.getValueState('value')   ->  EXPLICITLY_SET

That is not a rare back door: **55 sites** assign `leaf.value` directly, against
228 that go through a parent. So "state on the container" cannot simply be
declared --- either those 55 sites change, or the leaf reaches up to its parent
to record the change, which is precisely the two-way navigation the sketch was
trying to remove.

The honest position: assignment coercion is reducible from eighteen thousand
objects to about twenty classes, and value-state tracking is **not** reducible
without first deciding what `leaf.value = x` should mean when the leaf does not
know its parent. That question is upstream of the dataclass sketch, not a detail
inside it.

**CList.** Items are addressed by index and named `[n]`; see Defect C in the
companion note.

**Cached file content, which is not a parameter at all.** A `CDataFile` caches
a digest of what is *inside* the file --- cell, spaceGroup, wavelengths --- on
read, and it is then serialised as though it were data the user supplied.
Measured in one real project: **27 of 42** `params.xml` files carry it, and
**7** `input_params.xml` files do, including one configured by hand in the GUI
(8,099 bytes for a form somebody filled in).

It is derived, it can go stale --- the file may change while `params.xml` keeps
the old cell --- and `setEtree` restores it, so a reloaded or cloned job can
carry a cached cell that no longer matches its file without anything reading
the file to find out.

This is the signals leak one level up: **machinery surfacing as data because
nothing distinguishes them.** `value_dict_for_object` needed a name-based skip
list to keep signals out; `getEtree` has no way to know that `content` is a
cache rather than a value. Where a record's fields are declared, a cache simply
is not one of them and the question does not arise.

---

## Why this is attemptable now, and was not before

Everything built this week pins the **facade**, not the implementation:

- 21,014 container paths, and the ordering of every one
- **9,142 i2run addressing rows** --- which candidate a bare `--FLAG` resolves
  to, derived purely from container structure
- the GUI-rendered JSON shape for all 171 tasks
- the validity report for every task
- that serialisation reports exactly the children the tree has
- five tiers of conformance: access agreement, leaf operations, mutation,
  persistence, identity and lifetime
- the same invariants over the PHIL construction route

All of it is stated in terms of what a caller sees. **An implementation that
satisfies those numbers is a correct implementation, by construction** --- and
`scripts/snapshot_containers.py --diff` answers in 25 seconds, without running
a job.

That is the difference between this being a rewrite and being a refactor.

---

## Conformance is accountability, not frozenness

An earlier draft of this note demanded that nothing derived from a container
should change. That is wrong, and it would forbid every fix --- including ones
already shipped. The rule actually being followed is:

> **Every difference must be predicted before it is observed.** An unpredicted
> difference is a defect. A predicted one is the deliverable.

Defect A's acceptance test was never "no change". It was *0 paths gained, 1,358
removed, in these ten tasks* --- a shape stated in advance, so that the diff
could confirm the prediction and rule out anything riding along with it.
Defect B's was *13 unstable tasks become 0, total paths unchanged*.

The demand therefore gets **stronger** when a change is a genuine fix, because
that is exactly when something unintended can hide behind the intended
difference.

`snapshot_containers.py --diff` now reflects this:

- **gains fail the run outright** --- nothing should invent a parameter
- **losses and derived changes are reported**, for a human to check against
  what they predicted
- **`--strict` fails on any difference**, which is what a pure refactor asserts

Verified on a real case: defect A's removal passes by default and fails under
`--strict`.

### Match the instrument to the change

The container snapshot is not universal, and assuming it is would be the
subtler mistake. Cached file content --- a file's cell and spacegroup ---
**never appears in it**, because the snapshot builds fresh containers that have
read nothing. That fix belongs to the persistence tier and to a before/after
over real `params.xml` files.

Its expected shape, when someone takes it:

| | |
|---|---|
| container snapshot | **no difference** --- content is not present at construction |
| `params.xml` / `input_params.xml` | content elements **removed**: 27 of 42 and 7 of the files in one project |
| reading an old params file | **must still work** --- tier 4 holds this already |
| a cloned job | **recomputes** the cell rather than restoring a stale one; desirable, and visible |

---

## Suggested order

1. **`_children` key and contents.** Container-side only. Ordered, keyed on the
   `weakref.ref` rather than `id()` --- CPython reuses ids after collection, and
   tier 5 exists to catch that.
2. **Leaves stop being `HierarchicalObject`s.** 86% of objects, none of the
   hierarchy semantics. Value plus shared qualifiers; state on the container.
3. **Composition classes become dataclasses.** The 22 with fixed field lists.
4. **Context handed down, `_find_plugin_parent` and `_temp_plugin_ref` retired.**
5. **`contents()` as one declared answer**, so `value_dict_for_object`'s five
   strategies, `CCP4i2JsonEncoder`'s walk and `find_all_files` collapse onto one
   call instead of agreeing by coincidence.

Steps 1 and 2 are worthwhile on their own and do not commit to the rest.

### Progress

| Step | State | |
|---|---|---|
| 1. `_children` key and contents | **landed** | #305 |
| 2. Leaves stop being `HierarchicalObject`s | not started | groundwork done, see below |
| 3. Composition classes become dataclasses | not started | scope revised, see below |
| 4. Context handed down | not started | |
| 5. `contents()` as one declared answer | not started | |

Landed alongside, and not foreseen when this note was written:

| | | |
|---|---|---|
| Signals belong to the thing with subscribers | #304 | |
| An empty `CString` is not a user choice | #306 | 812 leaves, 14%, wrongly `EXPLICITLY_SET` at construction |
| `<sendWhen>` declares transmission policy | #307 | retires phaser's `requiredDefaultList` |
| Qualifier template merged once per class | `qualifier-cache` | ~6% off construction |

Step 2 acquired a prerequisite that took most of the effort: *"state on the
container"* presupposes knowing what the state model **is**, and it turned out
to be wrong at construction for 14% of leaves. That is
[value-state-design.md](value-state-design.md), and #306 is the fix. Nothing can
read `ValueState` to decide behaviour until it means what it says.

### Two findings that change the plan

**Leaves cannot become dataclasses, and the note should stop implying they
might.** Wrappers call `param.isSet()`, `str(param)`, arithmetic dunders,
everywhere --- `phaser_MR.setKeyword` leans on several in one expression. That
API is what the codebase promises plugin authors. The win in step 2 is therefore
*lighter objects*, not the absence of objects: **six of the sixteen attributes
--- `_children`, `_children_by_name`, `_child_storage`, `_event_handlers`,
`_properties`, `_sigmgr` --- are empty on all 1,454 leaves measured across 60
tasks, without exception**, so they can be removed one at a time rather than in
one replatform. Step 3's 22 composition classes remain genuine dataclass
candidates.

**Qualifiers are shareable, and are not shared.** 16,989 leaves across the
registry hold 16,989 separate `_qualifiers` dicts with only **1,214 distinct
contents** between them --- one identical `CString` set duplicated 2,452 times,
`CInt` 1,218, `CBoolean` 1,130 --- for 6.3 MB. They are a property of (class +
def.xml declaration), identical on every request, so they can be interned behind
copy-on-write. The audit for that is done: `set_qualifier` is the single
mutation entry point, nothing mutates the dict `qualifiers()` returns, and
`MakeLink` proves copy-on-write is *needed* rather than precautionary, since it
sets `allowUndefined` on ~20 parameters conditionally inside a method.

### Step 2 was aimed at the wrong cost

The note argues from a leaf's **16 instance attributes** and 352 bytes against 28
for a plain int. Both true, and both nearly irrelevant, because the profiling
was never done. Measured:

| | |
|---|---|
| the five eager dicts plus the `RLock` | 384 B, **232 ns** per object |
| per API request (~100 leaves) | **0.023 ms**, against 40--56 ms construction |
| across the whole registry | 3.9 ms, 6.2 MB |

Removing them buys about **0.05%** of construction, in exchange for touching
`_children_by_name`, which is read on every attribute access, and `_lock`. The
memory number is worse than it looks too: the registry is never resident, since
containers are built per request and discarded, so the live saving is ~38 kB.

Leaves are cheap because they are few and short-lived. The note reasoned from
the size of one object and never multiplied by how many exist at once.

### Where the time actually goes

Profiling one `servalcat_pipe` construction:

| tottime | calls | |
|---|---|---|
| **44.5 ms** | **146,790** | `cdata.py:107 __getattribute__` |
| 11.2 ms | 160,707 | `str.startswith` --- the `name.startswith('_')` inside it |
| 8.3 ms | 6,533 | `load_xml.py:534 _build_id_path_map` |

`CData.__getattribute__` runs on **every attribute access of every CData
object**. Against default lookup it costs 192 ns versus 15 ns --- **12.8x** ---
which at that call count is about **26 ms of a 56 ms construction**. Half the
cost of building a container is one method, and it is paid on every parameter
edit because containers do not persist.

It exists for one reason: children are stored in the hierarchy rather than in
`__dict__`, so ordinary attribute lookup cannot find them and every access has
to be intercepted to check. That was a deliberate choice --- one source of truth
for the tree --- and the interception is what it costs.

**So the real step 2 is to delete the override, by making `__dict__` the place
children live.** Default lookup then finds them with no Python-level code at
all, and the structure stays single because `__dict__` *is* the structure ---
the same argument that made `_children` an ordered dict in #305, applied one
level out. `_child_storage` already holds strong references by name, so a
per-name mapping exists; this would make it the one Python already consults.

Not attempted here. Before it is, three things need answering: what
`__setattr__` (already overridden, 5.4 ms) must do to keep the hierarchy in
step; whether anything depends on a child being absent from `__dict__`, notably
the serialisers that walk `vars()`; and what happens to a child whose name is
not a valid identifier. The prize is large enough to justify asking properly:
roughly half of every container build, on the hot path of the API.

### Removing `__getattribute__`: attempted twice, not yet viable

**First attempt** --- delete the override, store children in `__dict__` --- failed
32 unit tests. Two naming faults were in the way, and finding them was the
useful part.

The `fileContent` backing store was `self.content`, not underscore-prefixed, so
`__setattr__` registered it as a *child* named `content` while the property had
already registered the same object as a child named `fileContent`. One object
under two keys, with `__dict__['content']` left a stale `None` that only the
override concealed. That also **disabled the property**: the override consults
`_children_by_name` before the descriptor protocol, so the getter ran once and
never again --- its reload-if-unloaded branch is dead code beyond first access.
Fixed by renaming to `_content`, which lands on its own.

Second, `__setattr__` keys `_children_by_name` by the *attribute* name while
`_add_child` keys it by the child's `_name`, and they differ whenever a caller
writes `container.input_file = CDataFile(name="XYZIN")`.

**Second attempt**, with both fixed. Everything cheap said yes:

| | |
|---|---|
| unit tests | **1946 passed** |
| api/unit | 79 passed |
| container snapshot | byte-identical, `--strict`, **including `params.xml`** |
| API round-trip | identical across 173 tasks |
| construction | `prosmart_refmac` 40.6 → **23.0 ms**, `aimless_pipe` 41.6 → **22.7**, `servalcat_pipe` 53.2 → **30.7** --- about **43%** |
| unit suite wall clock | 28.7 s → 19.3 s |

**i2run said no: 7 of 182 failed**, four of them phaser, plus
`provide_sequence`, `servalcat` neutron and `substitute_ligand` with SMILES.
The symptom is a validity error --- *"Ensemble item requires either Identity or
RMS to be set"* --- and the failing job's `params.xml` holds the ensemble item's
`<structure>` but no `<identity_to_target>`, so a value the pipeline sets with
`pdbItem.identity_to_target.set(0.9)` never reaches serialisation.

Bisected: restoring **only** the override, with the `__dict__` storage and both
renames still in place, makes it pass. So the cause is the override itself, not
the supporting changes.

**It is not an i2run artefact.** The same scenario driven through the REST API
--- `tests/api/e2e/test_mr_pipelines_api.py::TestPhaserSimpleAPI::test_gamma_basic`
--- passes on the working tree in 9 s and fails on the refactor in 3 s, with the
identical error. Two independent drivers agreeing is what makes the bug real
rather than a harness fault, and the API test is much the faster reproduction
for whoever picks this up.

Not isolated further, and these are ruled out --- each measured identical on
both trees:

- `ENSEMBLES` is empty at construction, and after loading the failing job's
  `input_params.xml`; the container-level validity report is the same three
  errors either way
- `identity_to_target.set(0.9)` on a bare `CPdbEnsembleItem` gives
  `EXPLICITLY_SET`, `isSet(allowDefault=False)` True, and serialises
- the pipeline's own sequence run standalone --- `makeItem`, the
  `remove`/`append` loop, `structure.set`, `identity_to_target.set`, save ---
  keeps the value

What remains unexamined is the subjob path: `phaser_simple` populates the
ensemble in `process()` and hands it to `phaser_MR`, and the failing job's
`params.xml` carries the item's `<structure>` but not its
`<identity_to_target>`. Validation runs *before* `process()`, so the error is
raised against a container the pipeline has not filled in yet --- which is why a
probe placed at the `identity_to_target.set(0.9)` line never fires.

**What this establishes.** The prize is real and large --- about 43% of
container construction, on the hot path of every parameter edit --- and every
cheap instrument available, including two added this week, called the change
clean. Only running actual jobs found the problem. Worth remembering when
judging how far the fast tiers can be trusted for a change of this shape.

### `_children` cannot collapse onto `__dict__`

`_child_storage` went, because it was `__dict__` keyed identically --- 19,529
entries, all present, all the same object. `_children` looked like the same
kind of duplication: derived from `__dict__` by taking non-underscore entries
holding a `HierarchicalObject`, it reproduced `children()` **exactly, in
membership and in order, for all 19,700 objects** across the registry.

It still cannot go, and the conformance tier said so within a minute:

    test_a_child_whose_name_collided_is_still_returned
    test_a_collided_child_is_not_returned_twice
    test_two_children_with_the_same_name_...
    test_children_sharing_a_name_are_not_all_kept_alive

`__dict__` is keyed by name; `_children` is keyed by the weakref, which hashes
by referent *identity*. Two children sharing a name --- which every `CList`
item did before Defect C gave them distinct ones --- are two entries in
`_children` and one in `__dict__`. The second silently replaces the first.

So the 19,700 agreement was a fact about today's data, not an invariant, and
the distinction between those two things is the whole reason that tier exists.
A container can hold two children with one name; a `__dict__` cannot represent
that.

What would make it collapsible is not more measurement but a decided rule:
that a child's name is unique among its siblings, enforced at registration
rather than hoped for. Defect C moved `CList` towards that by numbering items,
and `_declared_content` gives containers somewhere to state their model, but
`HierarchicalObject` still permits collisions and four tests pin that it does.

### A caution about the evidence in this note

Every "byte-identical across 171 trees" recorded here and in the PRs above was
first measured with a `snapshot_containers.py` that could import the **wrong
tree**: it inserted `os.getcwd()` and trusted `import ccp4i2` to follow, and the
dev venv carries an editable install pinned to the main checkout. Run from a
worktree root it compared the main tree with itself and reported success.

All affected claims were re-measured from a clean worktree with the resolved
import verified, and all held. The tool now imports the tree it lives in and
refuses to start otherwise. More to the point, it has a **positive control** it
never had: renaming `freerflag`'s `FRAC` is shown to be caught. An instrument
nobody has watched fail is not evidence, and this note leans on that instrument
throughout.

---

## Appendix: reproducing the measurements

    cd server
    ccp4-python ../scripts/snapshot_containers.py before.json
    ccp4-python ../scripts/snapshot_containers.py --diff before.json after.json

The class-shape and leaf-cost figures came from short scripts walking
`get_plugin_class(task)(...).container` for every task in `TASKS`, recording
`type(obj).__name__` against `tuple(dataOrder())`. They are quick to rewrite
and were deliberately not committed as tests: they measure the present shape,
which is expected to change.

## Direction of travel: retire the stub / implementation split

Recorded 2026-08-31, as a preference rather than a scheduled step.

Every CData type is currently two classes: a generated `CFooStub` in
`core/cdata_stubs/` carrying the data model, and a hand-written `CFoo`
carrying the methods. That split was right when it was made --- the stubs were
machine-produced from the Qt-era definitions, and a "do not edit" half kept
regeneration from clobbering hand-written code.

The regeneration it protects is over:

| | |
|---|---|
| stub classes | 209 |
| with an implementation subclass | 208 |
| whose implementation is **empty** (docstring only) | **153 --- 74%** |
| whose implementation has 3+ members | 20 |
| things that invoke `migration/CData/stub_generator.py` | **none** |
| commits hand-editing the stubs since | 18, most recently 2026-06-23 |

The generator sits under `migration/`, nothing calls it, and the stubs have
been hand-edited eighteen times --- adding a KPI class, adding mask fields.
The banner is already false, so the contract it names is not being kept.

### It is not only a navigation cost

The two halves interleave in the MRO, and the implementation can drop out of
it entirely:

    CObsDataFile -> CObsDataFileStub -> CMiniMtzDataFile
                 -> CMiniMtzDataFileStub -> CMtzDataFileStub

`CMtzDataFile` is not in that chain. So `isinstance(obs, CMtzDataFile)` is
`False` for a file that plainly is one, and every method defined on
`CMtzDataFile` is unreachable from `CObsDataFile`. Twelve implementation
classes are bypassed this way, across 30 subclass relationships --- among them
`CSpaceGroup`, `CXmlDataFile` and `CI2XmlDataFile`.

Two modules already work around it in comments rather than fix it:

    # Import stub class for isinstance checks - subclasses like CObsDataFile
    # inherit from stubs (CMtzDataFileStub) not implementations (CMtzDataFile)

That is this note's recurring shape: a mechanism that conceals a defect rather
than causing one. Anywhere outside `lib/utils/files/` that tests
`isinstance(x, CMtzDataFile)` is quietly wrong, and nothing reports it.

### Preferred end state

One class per type. `CFoo` declares its fields as annotations and carries its
own methods.

1. **Collapse the 153 empty implementations** --- mechanical, no behaviour to
   preserve, removes 153 classes and one level of indirection.
2. **Merge the 55 that have members** --- stub fields become annotations on
   the implementation. This is what fixes the MRO bypass, and it is the half
   that needs the snapshot under `--strict`.
3. **Retire the "do not edit" banner regardless of when 1 and 2 happen.** It
   currently misleads about a contract nobody keeps, which is its own small
   version of the problem.

The annotations step (`5ad7c00b0`) is what makes this cheap: a stub's whole
content is now annotations plus a decorator, and the annotations are the
declaration. Collapsing two classes into one is a merge, not a re-derivation.

Blast radius is small but not zero: `lib/utils/files/digest.py` and
`upload_param.py` import stub names deliberately, and both do so *because* of
the bypass --- they get simpler, not harder, when it goes.

### `contents_order` is not the free step it looked like

Measured 2026-08-31, before attempting it. The expectation was that a
dataclass preserves declaration order, so `contents_order=` would become a
duplicate of the annotation order and could simply be deleted. It is not, for
a reason worth recording.

`dataOrder()` (`cdata.py:798`) resolves in three stages: an explicit
`CONTENT_ORDER`, then the MRO attribute walk, then `children()`. The second
stage **sorts alphabetically**. So declaration order is not preserved anywhere
today, and `contents_order` is the only thing holding a meaningful order:

| | |
|---|---|
| stub classes carrying `contents_order` | 44 |
| **load-bearing** --- order changes if removed | **33** |
| redundant --- order unchanged if removed | 11 |

Delete it from `CCellStub` and a unit cell renders `a, alpha, b, beta, c,
gamma`. Delete it from `CSpaceGroupCellStub` and the cell precedes the space
group. These are not cosmetic.

The clean fix is to make stage two use declaration order --- which is exactly
what a dataclass gives, so it is the same direction of travel. But it is a
behavioural change, not a refactor:

| | |
|---|---|
| stub classes whose order would change | **88 of 125** |
| ...that have a `contents_order` today | 16 |
| ...that have none, and so render alphabetically today | **72** |

Much of that change is an improvement. `CCootHistoryDataFile` currently renders
`annotation, baseName, contentFlag, dbFileId, project, relPath` --- alphabetical,
and worse than the `project, baseName, relPath, ...` its declaration states.
But 88 classes shifting is a visible change to every task interface built on
them, and it needs saying out loud rather than arriving inside a refactor.

**Recommendation.** Take it as a deliberate, separately-justified step:
switch stage two to declaration order, take the snapshot diff *without*
`--strict`, and read the changed orders as the deliverable. Then the 33
load-bearing `contents_order` entries can be checked against it and the
redundant ones deleted. Do not attempt it as a silent equivalence --- it is
not one.

#### Resolved: the alphabetical order was manufactured by a property

The measurement above was right about the symptom and wrong about the cause,
which turned a feared 88-class rewrite into one line.

`CONTENTS_ORDER` is a **property**, and its last resort was
`sorted(self.CONTENTS)`. `dataOrder()` stage 1 asks `hasattr(self,
'CONTENTS_ORDER')` --- always true for a property --- and takes the value if
non-empty. So stage 1 always won, the alphabetical fallback was imposed on
every CData that did not set `contents_order`, and **stage 3's MRO walk, which
already computed declaration order correctly, was unreachable**.

`children()` is `__dict__` insertion order, which is the order
`apply_metadata_to_instance` created the fields, which is the order they are
declared. So the fix is `sorted(self.CONTENTS)` -> `list(self.CONTENTS)`.

What it changes, from the snapshot over 171 tasks:

| surface | difference |
|---|---|
| paths | **0** --- nothing invented or lost |
| i2run addressing | **0 tasks** |
| validity report | **0 tasks** |
| `params.xml` written | **0 tasks** --- saved jobs and reruns unaffected |
| GUI-rendered shape | **164 tasks** --- the deliverable |

The scale is explained by `CDataFile`, whose declared order is `project,
baseName, relPath, dbFileId, annotation, subType, contentFlag` and which
rendered as `annotation, baseName, contentFlag, dbFileId, project, relPath,
subType` --- in every task, for every file parameter.

#### What `contents_order` is actually for

It orders an unordered bag, and it never contracted to span every field: named
fields come first in the given order, the rest follow. Only the meaning of "the
rest" changed here. With declaration order supplying it:

| | before | after |
|---|---|---|
| redundant --- restates declaration order | 11 | **28** |
| load-bearing | 33 | **16** |

The 16 survivors all do the one thing declaration order cannot express:
hoisting a subclass's own field ahead of the inherited ones, as
`CAsuDataFile` does putting `selection` before `project, baseName, ...`. MRO
order puts a parent's fields first, so no arrangement of declarations reaches
it. That is a real residual job, and the reason the argument stays.

Deleting the 28 redundant ones is a follow-on, and now provably a no-op.
