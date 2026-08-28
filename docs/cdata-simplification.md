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
means the container intercepts assignment. That is magic, but in one place
rather than eighteen thousand.

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
