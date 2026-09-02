# CData: the CCP4i2 data model

Every parameter of every task is a `CData` object. A job's inputs and outputs,
the files it reads, the numbers a user types into the interface, the KPIs it
reports — all of them are nodes in a tree of `CData`, and that tree is what
gets rendered as a task interface, serialised into `params.xml`, and validated
before a job runs.

This document covers what a CData class looks like, what its parts mean, and
how it is built. For how to write a *task*, see
[authoring-a-task.md](authoring-a-task.md); for the XML that describes a task's
parameters, see [def-xml-reference.md](def-xml-reference.md).

## The one-sentence version

**CData is a Django-shaped declarative model plus a per-instance override
layer.**

That is worth unpacking, because the first half will feel familiar and the
second half is where CData genuinely differs from anything you have used
before — and where most surprises live.

## Declaring a class

```python
class CCell(CData):
    """A unit cell"""

    class Meta:
        qualifiers = {"toolTip": 'Cell lengths and angles',
                      "helpFile": 'crystal_data#cell'}

    a = content(CCellLength, toolTip='Cell length a in A', guiLabel='a')
    b = content(CCellLength, toolTip='Cell length b in A', guiLabel='b')
    c = content(CCellLength, toolTip='Cell length c in A', guiLabel='c')
    alpha = content(CCellAngle, toolTip='Cell angle alpha in degrees', guiLabel='alpha')
    beta = content(CCellAngle, toolTip='Cell angle beta in degrees', guiLabel='beta')
    gamma = content(CCellAngle, toolTip='Cell angle gamma in degrees', guiLabel='gamma')
```

Three things are being said:

- **the docstring** describes the class, and is read by the API reference
- **`class Meta:`** carries options about the class as a whole
- **`content(...)`** declares one field: its type, and the qualifiers that
  belong to *that field*

If you know Django models, this is the same shape: fields as class attributes,
options in an inner `Meta`. It is deliberately that shape, because this is a
Django codebase and the convention needs no explanation.

### `content()`

```python
name = content(TypeOrItsName, **qualifiers)
```

The type may be the class or its name as a string. A string is resolved when
the field is built, which is how a field can name a class defined further down
the same module, or in a module that would import circularly. Both forms are
common; use the class where it is already in scope.

The qualifiers given here are the **declared** qualifiers for that field —
its label, its tooltip, its bounds, its permitted values. They are inherited by
subclasses and can be overridden by them.

### `class Meta:`

Options about the class rather than any one field:

| key | meaning |
|---|---|
| `qualifiers` | qualifiers for the class itself, e.g. its own `toolTip` |
| `error_codes` | codes this class adds, merged with its ancestors' |
| `contents_order` | display order, where declaration order will not do |
| `content_qualifiers` | qualifiers for a field the class *inherits* |

`contents_order` is needed rarely and for one specific reason: a class's own
fields come *after* the ones it inherits, because that is what MRO order gives.
Where a subclass needs its own field to come first, only `contents_order` can
say so. `CAsuDataFile` puts `selection` before `CDataFile`'s `project` and
`baseName` that way. Eighteen classes need it; the other two hundred do not.

`content_qualifiers` exists for the same asymmetry: to qualify a field the
class inherits, there is no local `content()` call to attach the qualifier to.
`CEnsemblePdbDataFile` gives the inherited `subType` its enumerators
(`unknown/model/homolog/fragment/heavy atoms`) that way.

## Semantics

### A field holds an object, not a value

```python
cell = CCell()
type(cell.a)          # CCellLength, not float
float(cell.a)         # the number
```

Every field is itself a `CData` — with a name, a parent, its own qualifiers and
its own validity. That is what lets the interface address `cell.a` as a widget,
`params.xml` record it by name, and validation report against it individually.

### Assignment coerces; it does not replace

```python
cell.a = 55.0         # cell.a is still a CCellLength, now holding 55.0
```

This is the plugin API, used at hundreds of call sites, and it is why the
declared type describes what a field *holds* rather than what you may assign to
it. A type checker reading `a = content(CCellLength, ...)` would object to
`cell.a = 55.0`; the runtime does the right thing.

The corollary is a real trap:

```python
container.RESOLUTION = other.RESOLUTION      # assigns the VALUE, correctly
```

but historically some code assigned a `CData` where a number was meant and the
object was silently adopted as a child, leaving the value unset. `set()`
unwraps a CData argument to its value for exactly this reason.

### Two layers of qualifiers

This is the part with no Django or dataclass equivalent, and the one to
understand before changing anything.

| layer | where | lifetime | mutable |
|---|---|---|---|
| **declaration** | `content(...)` and `Meta.qualifiers` | one per class, merged down the MRO and cached | no |
| **instance** | `self._qualifiers`, seeded from the declaration | one per object | **yes** |

A task can legitimately change a qualifier on one object at runtime:

```python
self.container.inputData.ENSEMBLES.setQualifiers({'listMinLength': 0})
```

`phaser_simple` does this before emptying a list that is declared non-empty.
`MakeLink` sets `allowUndefined` on about twenty parameters conditionally,
inside a method. `xia2_xds` rewrites `enumerators` and `default` at run time.

That is why a class's declared qualifiers are shared and immutable while an
instance's are private and mutable, and why `test_qualifiers_are_per_instance`
guards the isolation: containers are rebuilt on every API request, so a leak
between instances would carry one job's state into the next parameter edit by
anyone.

It is also why `dataclasses.field(metadata=...)` cannot host these. That
mapping is a read-only `mappingproxy`, shared by every instance — the right
model for a database column, whose constraints belong to the schema, and the
wrong one for a task parameter, whose constraints depend on what the user has
chosen so far.

### Set, unset, and default

A field distinguishes three states, and the distinction is load-bearing for
serialisation: only what a user actually chose should be written to
`params.xml`.

| state | meaning |
|---|---|
| `NOT_SET` | nothing has been assigned |
| `DEFAULT` | the declared default was applied |
| `EXPLICITLY_SET` | assigned, by a user or by a plugin |

`isSet()` takes `allowUndefined`, `allowDefault` and `allSet` to ask the
question different ways. Note that `CBoolean(False)` is falsy, so
`if container.FLAG:` skips a deliberately-false value — test
`if container.FLAG is not None:` or `isSet()`.

## Implementation

Three hooks build a class, in this order, all at class-creation time:

1. **`Content.__set_name__`** fires for every `content()` as the class body is
   executed, recording the declaration on the class. It has to be a hook rather
   than a later read of the class dictionary, because a subsequent statement in
   the same body can replace the attribute — `CResolutionRange` declares fields
   and then defines properties of the same name.
2. **`CData.__init_subclass__`** fires next, reads `Meta` if there is one, and
   builds the class's metadata: merged attributes, merged error codes, the
   qualifier template.
3. **`apply_metadata_to_instance`** runs per object, walking the MRO to
   assemble the full field set and creating each child with its parent, its
   name and its qualifiers.

Everything is merged along the MRO with the most-derived class winning, and a
class's own hand-written values winning over both. Error codes merge rather
than shadow, so a class declares only what it adds.

## What it deliberately is not

CData is **not** a dataclass, and `is_dataclass()` on a CData class is `False`
on purpose. A dataclass would make four promises CData cannot keep:

| a dataclass would | CData |
|---|---|
| compare field by field | compares by identity; `isSet()` compares a value against its default |
| disable `__hash__` when `eq=True` | CData objects are hashable |
| support `asdict()` | raises: a tree of live parented objects holding a lock |
| mean "assign a `CCellLength` here" | `cell.a = 55.0` coerces |

`content()` is a CCP4i2 declaration that reads like a dataclass field and
promises only what it delivers. That is a deliberate trade: unfamiliar rather
than misleadingly familiar.

## Reference

The fields of any class can be listed in the Sphinx documentation with

```rst
.. cdata-fields:: ccp4i2.core.CCP4XtalData.CCell
```

which reads the live class at build time — including which ancestor declares
each field — so it cannot drift from the code.
