# Value state: what it claims, what it does, and what it should be

**Status:** design note, 2026-08-28. A decision is needed before leaves stop
being `HierarchicalObject`s, because value state is the thing that cannot
simply move --- it has to be *chosen*.

**Related:** [cdata-simplification.md](cdata-simplification.md) has the wider
argument; this is the dependency it names.

---

## What it claims to be

```python
class ValueState(Enum):
    NOT_SET        = auto()   # never explicitly set
    DEFAULT        = auto()   # using a default from qualifiers
    EXPLICITLY_SET = auto()   # explicitly assigned
```

Three states, and a good three: they distinguish *a user chose this* from *this
is what the def.xml says* from *nobody has said*. That distinction is
load-bearing. `phaser_MR.setKeywords` skips keywords sitting at their default,
so a user who types the default value and a user who leaves the field alone are
meant to produce different command lines.

## What it actually does

Measured, on this tree:

```
fresh CString      state=EXPLICITLY_SET   isSet=False   value=''
set to ''          state=EXPLICITLY_SET   isSet=False
set to 'x'         state=EXPLICITLY_SET   isSet=True
fresh CInt         state=NOT_SET          isSet=False
CInt set to 0      state=EXPLICITLY_SET   isSet=True
fresh baseName     state=EXPLICITLY_SET   isSet=False   value=''
```

Three faults, and they compound.

**A `CString` is `EXPLICITLY_SET` from birth.** Before anything touches it. So
any logic that asks "did somebody choose this?" is answered wrongly for every
string in the tree --- and on a fresh `CDataFile`, `baseName`, `annotation` and
`relPath` are all in that state.

**`isSet()` does not report the state.** For a string it reports
*non-emptiness*; for an int it reports the state, so `0` counts as set. The
predicate that **719 call sites** depend on therefore means different things
depending on the type of the thing being asked.

**Emptiness is conflated with absence.** A user who deliberately clears a field
and a user who never saw it are indistinguishable.

The three-state model is, in practice, not what the system runs on. What it
runs on is `isSet()`, which is a type-dependent mixture of state and emptiness.

## And a fourth state is missing

`input_params.xml` for a job configured by hand contains, with
`exclude_unset=True` already applied:

    UNMERGEDFILES...file.contentFlag = 0        never chosen by the user
    UNMERGEDFILES...file.subType     = 0        never chosen by the user
    UNMERGEDFILES...crystalName      = New      never typed
    UNMERGEDFILES...cell.a           = 34.15    read out of the MTZ

None of those is a user decision. `contentFlag` is *determined from the file*;
`cell` is a digest *cached on read*; output paths are assigned by
`checkOutputData`. They are all facts the system worked out, and the model has
nowhere to put them, so they are recorded as though the user had typed them.

That is the missing distinction: **who set this**, not merely *whether* it was
set.

---

## What the states need to be

Four, and the fourth is the new one:

| state | meaning | persisted? | recomputed? |
|---|---|---|---|
| `UNSET` | nobody has given it a value | no | --- |
| `DEFAULT` | the value the schema specifies | see *provenance is not policy* | on load |
| `DERIVED` | the system worked it out from evidence | **as a cache, marked as such** | when the evidence changes |
| `USER` | chosen through the GUI, the API or i2run | **always** | never |

`USER` and `DEFAULT` are today's `EXPLICITLY_SET` and `DEFAULT`. `UNSET`
becomes honest for strings. `DERIVED` is what does not exist, and its absence
is why `contentFlag`, `cell` and output paths are indistinguishable from things
somebody typed.

The immediate payoffs:

- **`params.xml` records decisions, not deductions.** A reloaded or cloned job
  recomputes `contentFlag` and `cell` from the file rather than trusting a
  cached value that may no longer match it.
- **The stale-digest bug goes away by construction**, rather than being fixed.
- **`isSet()` can mean one thing.** See below.
- **An empty string a user typed is a value**, distinct from a field never
  touched.

## Provenance is not policy

A correction to an earlier draft of this note, which had `DEFAULT` implying
"need not be transmitted". Those are two different questions, and conflating
them is already causing trouble.

`phaser_MR.setKeyword` sends a keyword only when:

```python
if (not parameterObject.isDefault() and not str(parameterObject) == 'Auto') \
        or parameterName in self.requiredDefaultList:
```

Some parameters **must be given to phaser even at their default value**, or the
script does not work. `phaser_EP_LLG` names `PART_VARI`, `PART_DEVI` and
`LLGM`; `phaser_EP_AUTO` names the first two. The `requiredDefaultList` is the
patch for a distinction the model does not make.

So:

- **state** says where a value came from --- schema, system, or user
- **transmission** says whether the consuming program needs to be told
- **persistence** says whether the job's record needs it to be reproducible

Phaser needs `PART_VARI` on its command line whether or not anybody chose it.
And if a program requires a value, then a job's `params.xml` needs it too, or
the job cannot be reproduced from its own record --- so "never persist
defaults" is wrong for the same reason.

The existing mechanism is sound as far as it goes: `phaser_MR` is the base
class, `phaser_EP_AUTO` extends it and `phaser_EP_LLG` extends that, so the
empty list in `phaser_MR` is a base-class default the EP wrappers override.
phaser_MR genuinely requires none.

What it is not is *visible*. The knowledge that phaser must be told `PART_VARI`
even at its default lives in a Python class attribute inside one wrapper
family, where nothing else can see it: not the GUI deciding what to show as
meaningfully set, not `params.xml` deciding what a reproducible record
contains, not a reader wondering why a keyword appears in one job and not
another. A **qualifier in the def.xml, beside the default it concerns** ---
`alwaysSend`, or its inverse --- puts it where the parameter is defined and
where every consumer already looks.

This also means `isDefault()` cannot simply be deleted in favour of provenance.
Its nine callers are asking a policy question, and they need an answer whatever
the storage becomes.

## Replacing `isSet()`

719 call sites cannot be examined individually with confidence, so the
replacement must be mechanical:

    has_value()      any state other than UNSET --- "is there something to use"
    was_chosen()     USER only --- "did a person decide this"

Today's `isSet()` is *approximately* `has_value()` for ints and *approximately*
"non-empty" for strings. The mechanical move is:

1. make `isSet()` mean `has_value()` exactly, for every type
2. measure what changes --- the conformance snapshot reports validity and i2run
   addressing across 171 tasks, and `tests/i2run` exercises the rest
3. examine only the sites the measurement implicates

Step 2 is the point. The risk in this change is not mechanical, it is that some
caller depends on the emptiness reading; and the way to find those callers is
to make the change and see what moves, not to read 719 sites.

## Where the state lives

On the **container**, keyed by path --- not on the value.

That is what unburdens the leaves: `_value_states` and `_default_values` are
two of the sixteen attributes every `CInt` carries, and 86% of the objects in a
task tree are leaves. A single mapping per container replaces ~18,000 dicts.

It also matches how state is *used*: `getEtree(excludeUnset=True)` walks a
container asking each child; validity walks a container; the GUI renders a
container. Nothing asks a lone leaf about its state without having walked from
a container to reach it.

## What this note does not decide

**Whether `DERIVED` is persisted at all.** Recording it as a marked cache makes
reload fast and lets a stale value be detected; omitting it entirely makes
staleness impossible but means every load re-reads files. The measurement to
inform that: how long a digest takes against how often a job is reloaded.

**What happens to existing `params.xml` files.** 27 of 42 in one project carry
cached content, and 7 `input_params.xml` do. They must keep loading --- tier 4
holds that --- but a value that arrives with no provenance has to be assigned
one. Treating unmarked values as `USER` is safe and wrong; treating them as
`DERIVED` is correct and risks discarding a real choice. This needs deciding
explicitly rather than by default.

**Whether `DEFAULT` needs storing at all.** The schema knows the value, so
storing it looks redundant --- but see *provenance is not policy*: a job whose
program requires a defaulted parameter needs that parameter in its own record
to be reproducible. The decision is probably per-parameter rather than global,
which is an argument for the `alwaysSend` qualifier carrying it.

---

## Acceptance

Per [cdata-simplification.md](cdata-simplification.md): every difference
predicted before it is observed.

For the state model itself the predicted diff is **not empty**, and that is the
point --- `input_params.xml` should stop carrying `contentFlag`, `subType` and
`cell`. Expected shape:

| | |
|---|---|
| container snapshot | no difference --- state is not visible at construction |
| `validity` in the snapshot | no difference; `allowUndefined` checks read `has_value()` |
| `i2run addressing` | no difference; addressing is structural |
| `params.xml` | derived values removed or marked; user choices unchanged |
| old params files | still load (tier 4) |
| i2run | unchanged; a failure here means a keyword changed, which is the risk |

The last row is the one to watch. If a keyword reaching a program changes, it
will show as an i2run failure and not as a unit-test failure.
