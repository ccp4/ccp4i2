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

## What plugin authors have already invented

Worth looking, because where implementers have hand-rolled a mechanism they
have told us what the model failed to provide. There are two, and they solve
**opposite halves of the same missing distinction**.

**phaser: send it even though it is default.** `requiredDefaultList`, a class
attribute naming parameters the program must be given regardless ---
`PART_VARI`, `PART_DEVI`, `LLGM`. Inherited properly through
`phaser_MR` → `phaser_EP_AUTO` → `phaser_EP_LLG`.

**aimless: only send it if the user opted in.** Four booleans in the def.xml,
each labelled "override default …":

    OUTLIER_OVERRIDE        override default outlier rejection
    ANALYSIS_OVERRIDE       override default outlier rejection
    SDCORRECTION_OVERRIDE   override default SDcorrection
    INTENSITIES_OVERRIDE

```python
if not par.OUTLIER_OVERRIDE:
    return
if not par.OUTLIER_EMAX.isDefault():
    self.appendCommandScript("REJECT EMAX %f" % par.OUTLIER_EMAX)
```

That checkbox exists so a user can say *the values below are mine, not
defaults* --- a hand-rolled provenance signal, presented to the user as though
it were science. With reliable provenance it is unnecessary semantically,
though it may well be worth keeping as a **UI affordance**: "let me adjust
outlier rejection" is a reasonable thing to offer, and removing the checkbox
would change the interface. The point is that it should stop carrying meaning
the model ought to supply.

Two other things that block should be read for:

- `isSet()` and `isDefault()` are used within the same function for similar
  decisions --- `OUTLIER_COMBINE.isSet()` beside `OUTLIER_SDMAX.isDefault()`.
  Given that `isSet()` means non-emptiness for strings and state for numbers,
  adjacent lines are asking materially different questions.
- The whole file is only reachable behind an override flag, so the
  `isDefault()` calls inside it are a second filter on an already-filtered
  block.

### `isSet` already takes `allowDefault`

The vocabulary is not missing. `isSet()` has taken an `allowDefault` argument
all along --- *"if False, consider values that equal the default as not set"*
--- and it is load-bearing where it matters most: `PhilPluginScript.py:312`
uses it to decide what reaches a PHIL program at all, `i2run.py:110` to decide
what to echo, `CCP4Data.py` to decide what a copy carries, and `cdata.py:1100`
to decide what persists.

So the gap is not the concept. It is that the choice is made **at each call
site**, by whoever wrote that line, rather than declared once beside the
parameter --- which is why two wrappers reached for two different mechanisms
instead of the argument that was already there.

**And why aimless could not simply use it: the polarity of its default is
inverted between the base class and the leaves.**

| | `allowDefault` defaults to |
|---|---|
| `CData` (base) | `False` --- a defaulted value is *not* set |
| `CInt`, `CFloat`, `CString`, `CBoolean` | `True` --- a defaulted value *is* set |
| `CDataFile`, `CProgramColumnGroup`, `CPerformanceIndicator` | `False` |
| `CPdbEnsembleItem` | `True` |

The four types a control parameter actually is are the four that flip it. So a
bare `isSet()` on `OUTLIER_EMAX` answers "does it have a value", including the
default --- useless as a provenance test --- while the identical call on a
`CDataFile` excludes defaults and *is* one. `isDefault()` is what remains once
`isSet()` cannot be trusted to mean the same thing twice.

### The container does not persist --- so only the file carries provenance

In the client/server setting a container is not a living object. Every
parameter set over the API constructs the plugin afresh
(`get_plugin_with_context`), applies one value, rewrites the whole
`input_params.xml`, and discards it. `ValueState` therefore never survives a
request on its own account. Whatever provenance exists is whatever the file
encodes.

Measured, the file encodes it **soundly, by absence**:

| | |
|---|---|
| untouched container saved | 0 of 106 aimless control parameters written; 0 of 101 servalcat; 0 of 9 freerflag |
| set a value, save, reload | returns `EXPLICITLY_SET` |
| set a value *equal to the default*, save, reload | returns `EXPLICITLY_SET` --- the deliberate choice survives |
| repeat the round trip six times | stable; no ratchet |

So the rule on read is **present ⟹ explicit, absent ⟹ default**, and it is
correct in both directions. This is worth stating plainly because it retires
the concern that prompted the exercise: `params.xml` does *not* accumulate
values nobody set, and provenance does not decay toward "everything is
explicit" under repeated editing.

### Where the wrongly-set parameters actually come from

They are real, but they are not a persistence failure --- they are wrong the
moment the container is built, before anything touches it.

`CString` marks itself `EXPLICITLY_SET` at construction. Alone among the four
scalar leaves:

    CString    value=''      EXPLICITLY_SET
    CInt       value=0       NOT_SET
    CFloat     value=0.0     NOT_SET
    CBoolean   value=False   NOT_SET

Across all 171 constructible tasks that is **812 of 5,803 leaf parameters ---
14%, in 65 tasks --- claiming a user chose them before a user has seen them**
(790 `CString`, the rest its subclasses: `CFilePath`, `CSpaceGroup`,
`CSMILESString`, `CCrystalName`, `CDatasetName`). Aimless's
`CHOOSE_LAUEGROUP`, `CHOOSE_SPACEGROUP` and `EXCLUDED_BATCHES` are three of
them.

Nothing breaks today, for a reason worth noticing: **`CString.isSet()` does not
trust the state.** It overrides the base to return False for an empty string
regardless, so every consumer that asks `isSet()` gets the right answer from
the wrong reasoning, and the file stays clean because the write filter asks
`isSet()` too.

There are in fact **two write filters, and they disagree.** The live path for
`controlParameters` is `CContainer.getEtree(excludeUnset=True)`, which asks
`isSet()` (`cdata.py:845`) and so drops the empty strings. Beside it sits
`ParamsXmlHandler._is_explicitly_set`, which asks the *state* --- and returns
True for all 812. It governs a different branch today, so it never runs on
them; confirmed by calling it directly on aimless's `CHOOSE_LAUEGROUP`, which
answers True while the parameter is correctly absent from the file.

That is the load-bearing warning for this design. A model that reads
`ValueState` directly --- which is precisely what provenance-driven
transmission would do --- inherits all 812 as *"the user chose this"*. **The
state machine must be fixed before anything is built on it**, and the fix is
small: construct `CString` as `NOT_SET`, like its three siblings, and let the
emptiness override retire with the workaround it exists to be.

### The comparison branch is dormant, not absent

`isSet(allowDefault=False)` contains a second test beyond the state check: it
compares the value against the `default` qualifier and returns False on a match
even when the state is `EXPLICITLY_SET` (`cdata.py:519`). That would demote a
user who deliberately chose the default value.

It never fires. The def.xml `<default>` is applied to `value` and recorded as
`ValueState.DEFAULT`; it is not kept as a qualifier and not kept in
`_default_values`, both of which come back empty. So `get_qualifier('default')`
is None, the branch short-circuits, and `isSet(allowDefault=False)` is in
practice a *pure state test* --- which is the correct behaviour.

It is dormant rather than harmless. Anyone who later stores the default where
that lookup can find it --- a reasonable-looking tidy-up --- switches on
value-comparison demotion simultaneously at every `allowDefault=False` call
site, including the one deciding what a PHIL program is told. Delete the
branch rather than leave it armed.

### What that suggests for the def.xml

Move the decision from the call site to the declaration, beside the default it
governs:

| | |
|---|---|
| `always` | send it whatever its provenance --- phaser's `PART_VARI` |
| `ifChosen` | send it only when a user chose it --- aimless's outlier block |
| `ifSet` | send it whenever it has a value, today's leaf-type behaviour |

`always` replaces `requiredDefaultList`. `ifChosen` replaces the aimless flags
as *semantics* --- and, by consulting state rather than comparing values, stops
demoting the user who chose the default on purpose. Uniform declaration also
removes the reason the base/leaf polarity split exists.

### The GUI flags stay

The four aimless override booleans are not only backend gates: they drive
conditional visibility in the interface, through `useBoolToggle` and
`{outlierOverride.value && (…)}` in `ImportantOptionsTab` and
`AdditionalOptionsTab`. Revealing a block of controls is a real job, and a good
one --- "let me adjust outlier rejection" is a reasonable thing to offer.

So this is one flag doing two jobs, and only the second should go. Note that
the code already knows they are different: the flag reveals the block, and the
`isDefault()` test inside then filters what within it the user actually
touched. Ticking the box and changing nothing correctly sends nothing. The
proposal keeps the reveal and lets provenance do the filtering it is already
being asked to do by hand.

### One signature to fix on the way past

`CPdbEnsembleItem.isSet(attributeName=None, allowDefault=True)` does not accept
`field_name`, `allowUndefined` or `allSet` --- the three keywords the base
class passes in its own recursion over container children (`cdata.py:476`).
It does not fire today, because that loop short-circuits on an earlier set
child before reaching one. A landmine rather than a bug, but it should be
brought into line with the others whenever this area is touched.

## Decisions taken

**`DERIVED` is not persisted.** It has not been historically, and continuing
not to is the conservative choice: nothing can go stale if nothing is stored.
Recording it later --- as a marked cache, so staleness is detectable --- would
be an enhancement, and is only worth doing once we can say concretely how a
consumer would use it. Whatever is done must be backward compatible, since a
reader of today's files must keep working.

**Legacy `params.xml`: presence means explicit, absence means default.**
Not the "explicit unless it equals the default" rule considered earlier --- the
measurements above retire that. Absence already carries "default" unambiguously,
so comparing values adds nothing and costs something: it would demote the user
who deliberately chose the default, introducing exactly the loss the dormant
branch above is one tidy-up away from causing. The simpler rule is also the one
the code already implements, so legacy files need no special reader.

**`DEFAULT` is stored sometimes**, per parameter, governed by the vocabulary
above rather than by a global rule.

## What this note still does not decide

**Whether `ifChosen` should apply per parameter or per block.** Aimless's
override flags gate whole sections; expressing that as a qualifier on each
member is repetitive, and on the container is a different kind of statement.

**What to do about `isSet()`/`isDefault()` being mixed within a block.**
Mechanical replacement will preserve today's behaviour including its
inconsistencies; deciding which was *meant* needs someone who knows aimless.

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
