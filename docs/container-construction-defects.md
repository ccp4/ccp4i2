# Brief: two container-construction defects behind silent wrong answers

**Status:** FIXED, 2026-08-27. Both defects are corrected and the fixes are
measured below under *What landed*. The investigation that follows is kept as
written, because the evidence is the argument.

A third defect, C, was found while fixing these and is **not** yet fixed: a
CList item's name is not the handle it looks like. See *Defect C* at the end.

**Related:** [error-handling-remediation.md](error-handling-remediation.md) —
this is the same failure class that tracker exists for: the job runs, reports
success, and did something other than what was asked.
[i2run-parameter-syntax.md](i2run-parameter-syntax.md) is a separate proposal
that *assumes* these two fixes have landed.

---

## The task

Evaluate, and if you agree, fix two defects in how CCP4i2 builds a task's
container. Neither is an i2run defect, though i2run is where they surface most
sharply. Both produce answers that are wrong without being noisy.

- **A.** The def.xml parser attaches every nested container to the root a second
  time, creating a parallel ghost tree that nothing reads (2,304 parameters,
  15 tasks).
- **B.** `children()` iterates a `set` of weak references hashed on object
  identity, so sibling order is arbitrary — per process *and* per object. Where
  two parameters share a name, which one a bare `--FLAG` means is decided by a
  coin flip (108 parameters).

Everything below was measured on this tree, CCP4-free, in `server/.venv`.
Reproduction snippets are inline. Please verify rather than take on trust — I
got the first version of defect A's fix wrong (see *The trap*).

---

## Defect A — the def.xml parser duplicates nested containers onto the root

### Mechanism

`core/task_manager/def_xml_handler.py:133`, in `parse_def_xml`:

```python
# Parse any additional containers at root level (skip those already processed)
all_containers = root.findall(".//container[@id]")   # descendant-or-self: EVERY container
processed = set()                                    # ...but only direct body children go in
if body is not None:
    for c in body.findall("./container[@id]"):
        processed.add(c.get("id"))
    ...
for container in all_containers:
    if container_id not in processed:
        self._parse_container(container, root_container)   # re-parsed onto ROOT
```

`.//container[@id]` matches at any depth; `processed` holds only the *direct*
children of `<ccp4i2_body>`. Every container nested one level deeper is therefore
parsed a second time and attached to the root as an **independent object tree**.

### Evidence

`phaser_MR.def.xml` declares `resolution` only inside `keywords`, and no def.xml
anywhere declares a top-level `resolution`. Yet:

```python
p = get_plugin_class('phaser_pipeline')(parent=None)
p.container.resolution.RESO_HIGH is p.container.keywords.resolution.RESO_HIGH
# -> False
```

Scale, after expanding `<file>` inheritance with `load_nested_xml`:

| | |
|---|---:|
| Tasks affected | 15 |
| 9 phaser tasks — all 23 `keywords` sections ghosted | 252 params each |
| `servalcat_pipe` — `coordination` ghosted | 36 params |
| `cmapcoeff`, `comit`, `fft`, `validate_protein`, `cpatterson` | re-set names root already had (harmless) |
| Ghost parameters in total | **2,304** |

### Why it matters

`phaser_MR.setKeywords()` iterates `self.container.keywords.dataOrder()`, then
each section's `dataOrder()`, calling `inputObject.set<NAME>()`. It walks the
`keywords` subtree only — the root-level copy is read by nothing.

But i2run's bare-flag alias resolves a collision to the *shortest* path, and the
ghost is one segment shallower (3 against 4). So every bare phaser keyword flag
sets the copy phaser never reads:

```
--PEAK_ROTA_CUTO  ->  container.peaks_rotation.PEAK_ROTA_CUTO    (ghost)
--MACM_PROT       ->  container.macmr.MACM_PROT                  (ghost)
```

`--RESO_HIGH` is a special case worth knowing: the resolution parameters *were*
deliberately promoted. `phaser_MR_AUTO.def.xml` (and FRF/FTF/PAK/RNP) declare
`inputData.RESOLUTION_HIGH` / `RESOLUTION_LOW`, `phaser_MR_AUTO.py:87-92` reads
them into `setRESO`/`setHIRES`, and that is exactly why `RESO_HIGH`/`RESO_LOW`
sit in `phaser_MR.excludedKeywords`. `--RESOLUTION_HIGH` is unambiguous and
works today; `--RESO_HIGH` is a ghost of a switched-off keyword.

### The trap

**Do not simply delete the loop.** Five def.xml files — `cmapcoeff`, `comit`,
`cpatterson`, `fft`, `validate_protein` — open with a `<ccp4i2_body>` carrying
**no `id` attribute**, so `root.find(".//ccp4i2_body[@id]")` returns `None`,
`_parse_body` never runs, and the loop is the *only* thing building their
container. Deleting it strips `inputData`, `outputData` and `controlParameters`
from those five tasks entirely. I proposed the deletion before measuring it.

### Proposed fix

```python
body = root.find(".//ccp4i2_body[@id]")
if body is None:
    body = root.find(".//ccp4i2_body")          # (a) accept an id-less body

visited = set()
if body is not None:
    for c in body.findall(".//container[@id]"):
        visited.add(id(c))                      # (b) identity, not id string
    self._parse_body(body, root_container)

for c in root.findall(".//container[@id]"):
    if id(c) not in visited:                    # only genuinely unvisited containers
        self._parse_container(c, root_container)
```

Measured by building every task's full container tree both ways and diffing:

| | |
|---|---:|
| Task trees **byte-identical** after the fix | **153** |
| Task trees that change | 10 |
| Paths **gained** anywhere | 0 |
| each of 9 phaser tasks | −149 ghost paths |
| `servalcat_pipe` | −17 ghost paths |

### Consumers, checked

- **Python** — no module reads a root-level ghost section.
- **Frontend** — no binding references a ghost path. The phaser interfaces bind
  `itemName="keywords"` (the live subtree, by full path) plus unique names like
  `RESOLUTION_HIGH`.
- **Persisted jobs** — the ghosts *are* serialised into `params.xml` (29
  top-level elements for phaser). Re-loading an old params file into a fixed
  container is safe: `setEtree(..., ignore_missing=True)` drops the unknown
  elements and the real `keywords.*` values come back intact. The only loss is
  keyword values written *only* to the ghost, which never reached phaser.

### Open question for you

Is normalising the five id-less `<ccp4i2_body>` elements — adding the `id`
attribute to those def.xml files — better than teaching the parser to tolerate
them? That would let part (a) go away. It edits data rather than code, and I do
not know whether anything else keys off the body id.

---

## Defect B — sibling order is non-deterministic

### Mechanism

`core/base_object/hierarchy_system.py:62`:

```python
self._children: Set[weakref.ReferenceType] = set()
```

A `weakref.ref` hashes by its referent's hash, which for these objects is the
default identity hash. So `children()` iterates in allocation-address order:
arbitrary per process, and arbitrary per object *within* a process.

Two sites then propagate it:

- `hierarchy_system.py:247` — `children()` itself.
- `core/CCP4PluginScript.py:496` — `loadContentsFromXml` re-attaches the parsed
  def.xml tree with `for child in parsed_container.children()`, which is where
  document order is destroyed.

### Evidence

`i2run`'s `_compute_minimum_paths` breaks a name collision with
`min(indices, key=lambda idx: len(paths[idx]))`. When the candidates have the
same depth, that falls through to extraction order — i.e. to `children()`.
Running the same flag twelve times:

```
molrep_pipe              --F_SIGF            {'outputData': 7, 'inputData': 5}
xia2_ssx_reduce          --DIALS_INTEGRATED  {'outputData': 9, 'inputData': 3}
adding_stats_to_mmcif_i2 --FPHIOUT           {'inputData': 8, 'outputData': 4}
shelxeMR                 --PERFORMANCE       {'inputData': 11, 'outputData': 1}
```

`i2run molrep_pipe --F_SIGF native.mtz` populates an **output** slot on roughly
half its invocations. Same command, same tree, different meaning run to run.

The instability also reaches the serialised container: for the 29-section phaser
tasks, top-level `dataOrder()` differs between two instantiations of the same
task, so section order in the JSON the GUI consumes is arbitrary.

### This is not self-correcting

It is tempting to assume a mis-landed output assignment is harmless because the
machinery recomputes output paths. It does not. `checkOutputData` guards every
assignment with

```python
if not basename_is_set:
    ...  # setOutputPath / setFullPath into the job directory
```

so an explicitly set output `baseName` is left alone. And `CPluginScript.validity()`
separately filters `outputData` errors on the grounds that "outputs are set
during execution" (`CCP4PluginScript.py:1043-1051`). So the coin flip sticks, is
not recomputed, and is not reported.

### Proposed fix

The order is already recorded. `_add_child` maintains
`_children_by_name: Dict[str, weakref.ReferenceType]` alongside the set — a
dict, therefore insertion-ordered. Coverage, measured over all tasks with the
plugins held alive:

| | |
|---|---:|
| children visited | 10,141 |
| recoverable from `_children_by_name` in insertion order | **10,120 (99.79%)** |
| not recoverable | 21 |

All 21 misses are `CList` items, every one named `'[0]'` — they collide on the
name so only the last survives. That does not matter: a `CList` orders its items
from its own `_items`, not from `children()`.

```python
def children(self):
    with self._lock:
        self._cleanup_dead_children()
        out, seen = [], set()
        for _nm, ref in self._children_by_name.items():   # insertion order
            c = ref()
            if c is not None:
                out.append(c); seen.add(id(c))
        for ref in self._children:                        # unnamed / name-collided
            c = ref()
            if c is not None and id(c) not in seen:
                out.append(c)
        return out
```

Fixing `children()` also fixes `CCP4PluginScript.py:496`, because the parsed
container's name cache was filled in document order.

### Measured with the change applied

| | |
|---|---|
| `phaser_pipeline` top-level order, 5 fresh instances | identical |
| …and equal to def.xml document order | yes — `inputData, guiParameters, keywords, outputData` |
| Tasks whose `dataOrder()` still differs between two instantiations | **0** |
| `--F_SIGF`, `--DIALS_INTEGRATED`, `--FPHIOUT`, `--PERFORMANCE` × 12 runs | all stable, all `inputData` |
| `tests/unit/` | 1,696 passed, 35 skipped — identical to baseline |
| `tests/api/unit/` | 79 passed, 92 skipped — identical to baseline |

(`tests/db/` has a pre-existing collection error, unchanged either way.)

The last row of the tie-break table is the interesting one. **Declaration order
already encodes the priority you would otherwise have to invent.** Of the 163
def.xml trees, 158 declare `inputData` before `outputData` and **none** declares
it the other way (5 declare neither). So preserving order delivers
"inputData wins" without anyone writing that rule down, and covers
`controlParameters` and `keywords` at the same time.

### What it costs

Order becomes a contract where today there is none. def.xml declaration order
would start to govern GUI section order, `params.xml` element order and any
future introspection output, so re-ordering a def.xml becomes a visible change.
That looks like a fair price — and arguably the order authors already believe
they are choosing — but it is a decision, not a free win.

### Open question for you

Should the ordering guarantee be *documented* (a line in the def.xml reference
saying declaration order is preserved and observable), or left as an
implementation detail that merely happens to be stable? Documenting it is what
makes it safe to rely on for the tie-break rule.

---

---

## What landed

Both fixes are in, measured over the 171 loadable tasks, before and after.

**Defect B.** `children()` takes its order from `_children_by_name` (a dict, so
insertion-ordered) and its **membership** from `_children`. The two are not
equivalent, which cost a first attempt --- see *The membership trap* below.

| | |
|---|---:|
| tasks with unstable top-level order across three instantiations | 13 -> **0** |
| total paths | 22,372 -> 22,372 (order only) |
| `--F_SIGF`, `--DIALS_INTEGRATED`, `--FPHIOUT`, `--PERFORMANCE`, 12 runs each | all `inputData`, was 7/5, 9/3, 8/4, 11/1 |

"inputData wins" came free, and for a stronger reason than the 158-of-163
count: the base class pre-creates `inputData, outputData, controlParameters,
guiAdmin` in that order *before* the def.xml is parsed. `fft` declares
`inputData, controlParameters, outputData` and still builds `inputData,
outputData, controlParameters`. The tie-break is structural, not a fortunate
convention.

**Defect A.** The parser tracks element **identity** for every descendant of
the body, not id strings of its direct children, and accepts an id-less
`<ccp4i2_body>`.

| | |
|---|---:|
| trees byte-identical | 161 |
| trees that change | 10 (9 phaser tasks -149 each, servalcat_pipe -17) |
| **paths gained anywhere** | **0** |
| paths removed | 1,358 |

The five id-less-body tasks keep every path they had. Tolerating an id-less
body in the parser was preferred to normalising those five def.xml files: the
parser change covers a sixth file nobody has written yet.

**Both open questions closed.** (a) Tolerate in the parser, as above. (b) The
ordering guarantee is documented in `def-xml-reference.md`, because a tie-break
nobody has written down is not one anybody can rely on.

**Testing.** 25 unit tests; with the fixes reverted, 6 of 11 ordering tests and
7 of 14 ghost tests fail, so they detect rather than decorate. Unit suites
unchanged: 1863 under ccp4-python, 1813 CCP4-free, 79 API. i2run
**169 passed, 0 failed, 13 skipped** --- identical to the `post-stack`
comparator.

### The membership trap

The first version of the B fix read `_children_by_name` alone, and the i2run
suite caught it: five phaser tests failed with `'CFloat' object has no
attribute 'BOXS'`.

`_remove_child` always drops the set entry, but clears the cache entry only
when `_children_by_name[child._name]` still points at that child --- and a
child renamed since it was registered no longer matches, so the entry is
orphaned. Reading the cache for membership resurrected a removed child: a
CFloat keyword came back holding a child named after itself, `dataOrder()`
reported it, and `phaser_MR.setKeywords` --- which treats any child with a
non-empty `dataOrder()` as a sub-container --- descended into a leaf and asked
the float for a `BOXS` of its own.

So: the cache decides order, the set decides membership. Each is authoritative
for exactly one of them. The orphaned entry is the same stale-key root cause as
Defect C, which is what makes C worth doing rather than merely tidy.

## Why the test suite never caught either of these

Worth understanding before you change anything, because it tells you what a
regression test has to do.

**1. The ghosts are stable-wrong, not flaky.** The `min(len(path))` tie-break
only falls through to `children()` order when candidates have equal depth:

| Tie-break | Params | Behaviour |
|---|---:|---|
| Unequal depth — `min(len)` decides | **288** | deterministic |
| Equal depth — falls through to set order | **108** | non-deterministic |

All 252 phaser ghosts are unequal depth (3 against 4), as are servalcat's 36. So
defect A resolves the same way every run — always to the ghost. A stable wrong
answer is invisible to a test that does not assert on that parameter, and no
i2run test passes a phaser keyword-section parameter at all.

**2. The unstable 108 are barely exercised.** Harvesting every `--FLAG` the
suite passes and checking it against that task's ambiguous set: of **374
(task, flag) pairs, exactly one** uses an ambiguous flag —
`servalcat_pipe --XYZIN` — and that one is unequal depth, so it resolves stably
and to the right object. Two of the five equal-depth tasks are not in the suite
at all; two more collide on names the tests never pass.

**3. Where it did bite, it was worked around in the test.** The fifth is
`molrep_pipe`, whose collisions are exactly `F_SIGF` and `FREERFLAG`:

```python
args += ["--inputData.F_SIGF",    demoData("gamma", "merged_intensities_Xe.mtz")]
args += ["--inputData.FREERFLAG", demoData("gamma", "freeR.mtz")]
args += ["--XYZIN",               demoData("gamma", "gamma_model.pdb")]
```

Long form for the two that collide, short form for the one that does not.
`git log -L11,12:...test_molrep.py` puts it in the commit that created the suite
("Pytest testing to replace test101", #62) — so the flake was met and
disambiguated on day one, in the test rather than in the parser.

---

## What a good answer contains

1. **Independent verification of both defects**, ideally not by re-running my
   snippets. Everything here reproduces CCP4-free in `server/.venv`.
2. **A judgement on whether these are worth fixing now**, given the error strand's
   other priorities. My view: B is small, self-contained and a correctness fix; A
   is larger but its blast radius is measured and narrow.
3. **Regression tests that would have caught them.** Neither is covered today,
   and both are cheaply testable without CCP4:
   - A: assert no task's container has a top-level container that its def.xml
     declares only as a nested one — or simply pin the container tree of
     `phaser_pipeline` and `cmapcoeff`.
   - B: assert `dataOrder()` is identical across two instantiations of every
     task, and that it matches def.xml document order for a sample. This is a
     one-screen test over `TASKS` that would have failed from day one.
4. **Answers to the two open questions** (id-less `<ccp4i2_body>`; whether to
   document the ordering contract).
5. **An i2run baseline either side.** `server/run_i2run_baseline.sh` writes JUnit
   XML to `server/.test-baselines/<label>/`; diff `results.xml`. Expect no
   change from B, and none from A either — but A touches container construction
   for every task, so the baseline is the evidence, not the argument.

## Non-goals

Changing `CData`, changing the def.xml *format*, or any i2run syntax change.
Syntax is the separate proposal in
[i2run-parameter-syntax.md](i2run-parameter-syntax.md), which assumes this work
has landed and does not depend on any of it being done a particular way.

## Where to look

| | |
|---|---|
| `core/task_manager/def_xml_handler.py:133` | defect A, the hoisting loop |
| `core/CCP4PluginScript.py:496` | re-attach by `children()` — propagates B |
| `core/base_object/hierarchy_system.py:62,247` | defect B, the weakref set |
| `core/base_object/cdata.py:664` | `dataOrder()` — preferred order, then `children()` |
| `cli/i2run/i2run_components.py:164` | the `min(len(path))` tie-break |
| `core/CCP4PluginScript.py:1403` | `checkOutputData`'s `if not basename_is_set` guard |
| `pipelines/phaser_pipeline/wrappers/phaser_MR/script/phaser_MR.py:87` | `setKeywords` walks `container.keywords` |
| `server/.venv` | stock CPython, no CCP4 — everything above reproduces here |

---

## Defect C — a CList item's name is not the handle it looks like

For appending to docs/container-construction-defects.md once the A/B i2run run
finishes and the containers worktree can be rebased onto django (which now
carries the brief, via #295).

### What is already right

`CList.append`, `.insert`, `.remove` and `.pop` all renumber `_items`, so
`item._name` genuinely tracks the index, and `object_path()` already yields
`container.inputData.DICT_LIST[2]`. The coherence between index and name is
maintained; nobody has to build it.

### What is wrong

Measured on the patched worktree, servalcat_pipe's DICT_LIST:

    after 3 appends   _name=['[0]','[1]','[2]']  children=['[0]','[1]','[2]']  cache=['temp_item_0','temp_item_1','temp_item_2']
    after pop(0)      _name=['[0]','[1]']        children=['[0]','[0]','[1]']  cache=['temp_item_0','temp_item_1','temp_item_2']
    find_child("[0]") -> None

**C1. The name cache is keyed on the name the child had at registration.**
`append` calls `set_parent()` *before* the rename, so `_children_by_name` holds
`temp_item_0`, never re-keyed when `_name` changes. `find_child('[0]')` returns
None: the meaningful name cannot be looked up. It also means defect B's fix
orders list children by *registration* order rather than by index — correct
today, but by luck rather than construction.

**C2. `pop`/`remove` renumber the survivors and never detach the departing
item.** After `pop(0)` the list holds two items and three children, two of them
named `[0]` — a stale child colliding with a live one. `_child_storage` keeps a
strong reference under the old key, so it is a leak as well as a lie.

**C3. `__setitem__` is a third case, and the worst of the three.**

    value.set_parent(self)
    value.name = f"{self.name}[{index}]"   # ".name", not "._name" -- and "DICT_LIST[2]", not "[2]"
    self._items[index] = value             # the displaced occupant stays a child

An item assigned by index gets a name in a different *format* from one that was
appended, written to a different attribute, and the item it replaced is never
detached. Same leak as C2, plus a naming inconsistency that would make any
path-based addressing of list items ambiguous.

### Proposed fix (~30 lines, no new cost)

- `HierarchicalObject.rename()` — moves the parent's `_children_by_name` and
  `_child_storage` entry, O(1). CList's four renumber loops call it instead of
  assigning `_name`; those loops are already O(n), so complexity is unchanged.
- `pop`/`remove`/`clear` call `set_parent(None)` on the departing item, O(1),
  which closes the leak and the collision together.

- `__setitem__` uses `rename()` and the `[n]` format, and detaches what it
  displaces.

Then `[2]` means what it says: `find_child('[2]')` works, `children()` is in
index order by construction, and no child exists that is not in the list.

### Agreed

Its own before/after, after A/B is flattened into django. Not folded into the
A/B measurement, since it changes list identity semantics.

### Why this is the groundwork for reorderable list widgets

A move is a renumber of a contiguous span: the same objects stay children
throughout, nothing detaches. So `CList.move(from, to)` is a handful of lines
*once rename() exists*, and close to impossible without it --- the cache would
keep its registration-order keys, so `children()` would report the original
order while `_items` reported the new one, and the GUI and the backend would
quietly disagree.

Beyond C, reordering needs: `CList.move()`; an opt-in def.xml qualifier
alongside listMinLength (46 uses) and listMaxLength (3), meaning "order is
load-bearing in this task" --- DICTOUT_LIST[0] is *the* dictionary in
SubstituteLigand, whereas the dictionaries handed to servalcat are a set; and a
move endpoint, which is not a value edit. That last is the only real cost, and
it is client-side: a move changes what every path below it refers to, so the
frontend's path-keyed lookup for that subtree must be rebuilt rather than
patched.

Caveat to test rather than assume: paths are positional, so anything holding a
path to a list item goes stale across a move --- a validation error's `name`
field (`...DICT_LIST[1].fullPath`) being the obvious one. Errors are recomputed
each poll, so it should be benign.

---

## Notes towards a holistic CData reevaluation

Gathered 2026-08-27 while the A/B baseline ran. All measured on this tree.

### The access mechanisms plugins actually use

| mechanism | uses | files |
|---|---:|---:|
| dotted `container.a.b` | 3711 | 146 |
| `list[n]` | 2016 | 200 |
| `int()`/`float()`/`str()` on an attribute | 526 | 87 |
| `hasattr(obj, name)` | 328 | 60 |
| `[-1]` | 362 | 94 |
| `getattr(obj, name)` | 207 | 46 |
| `==` on an attribute | 102 | 28 |
| assignment `container.a.b = v` | 28 | 14 |
| `setattr(obj, name, v)` | 45 | 16 |
| **`dataOrder()`** | **13** | **8** |
| `children()` / `find_child()` | 0 | 0 |

The last two rows matter most. `children()`/`find_child()` are core-only, but
the core uses them on everyone's behalf. `dataOrder()` has 13 uses and is
load-bearing: it is how a plugin asks *is this a container or a leaf*, and
answering wrongly is what broke five phaser tests today. Frequency is a poor
guide to importance here.

### Protocol surface --- measured, not inferred

Nothing implements `__format__`, so **every** format spec raises:

    f"{x}"          works (via __str__)
    f"{x:.2f}"      TypeError   (CFloat)
    f"{x:03d}"      TypeError   (CInt)
    f"{x:>6}"       TypeError   (CString)

The 20 f-string uses in plugins today are all spec-less, so it is latent --- and
the first author to write `f"resolution {res:.2f}"` inside a try/except gets a
silent failure of exactly the kind the error-handling strand exists to kill.

Asymmetries that need a decision:

- `int(CInt)`, `float(CFloat)` work; `int(CBoolean)` raises, though `int(True)` is 1
- `abs(CFloat)` works; `abs(CInt)` raises; `-CInt` raises; `divmod` raises
- `CInt + 1` returns a **CInt**, so `x += 1` rebinds to a new *parentless* CInt.
  Zero uses in plugins today --- untested territory, not established behaviour
- `CString == "abc"` is True but `hash` differs, so `{"abc": 1}[CString]` raises
- `os.fspath` fails on CString and CDataFile: neither can go to `open()` or
  `Path()` without `str()`. That is what the 526 explicit conversions are paying for
- `CList` has no `__delitem__` (`del lst[0]` raises, `lst.pop(0)` works) and no
  `__eq__` (two empty CLists compare unequal)
- `CContainer` has `__getitem__` and `__len__` but no `__iter__`, so iteration
  falls back to the old integer-index protocol

### The signal overlay has no subscribers

| signal | emitted | connected |
|---|---:|---:|
| `child_added` | 1 | **0** |
| `child_removed` | 1 | **0** |
| `dataChanged` | 1 | **0** |

`connectSignal` is used for three names in the whole tree --- `finished` (16),
`jobId` (1), `progressUpdated` (1) --- all plugin lifecycle, none data change.
The frontend polls: useSWR in 18 files, refreshInterval in 14, no EventSource,
no WebSocket of ours (the hits are Next.js hot-reload artefacts), no
StreamingHttpResponse or event-stream endpoint. Confirmed with the user, who
built it.

Yet `HierarchicalObject.__init__` creates a SignalManager and four signals
(`destroyed`, `parent_changed`, `child_added`, `child_removed`) for **every**
object.

### What moving signals to CPluginScript would buy

Measured in CPU time (wall clock was contaminated by a concurrent i2run run),
three runs each, building all 171 loadable task trees:

| | CPU |
|---|---|
| current | **3.89 s** (3.89 / 3.85 / 3.94) |
| signals only where used | **3.44 s** (3.42 / 3.42 / 3.49) |
| | **-11.6%** |

23,093 nodes, 135 per task tree, so ~115,000 signal objects per full-registry
build with no subscribers. Only 0.09 s of the saving is the constructor itself
(3.9 us for a manager plus four signals); the rest is allocation, GC pressure
and the teardown path. For calibration, 23,093 RLocks cost 0.003 s and three
empty collections 0.002 s --- the signal overlay is the only per-object cost
worth removing.

Stubbing the four signals broke every task build with exactly one error:
`CPluginScript.__init__` line 323, `self._signal_manager.create_signal("finished", dict)`.
That is the only consumer, which is the whole argument.

Caveat: this measured *not creating them eagerly*, not the real refactor. A
genuine move should do slightly better (no slot at all).

### The conformance harness

Full cost is affordable: 171 tasks x 3 instantiations plus a walk of 21,014
paths is **18.5 s**. Precedent for the shape is
`test_task_registry_conformance.py` --- 823 assertions in seconds.

For every task, for every path: assert the mechanisms **agree with each other**,
which is stronger than any single-mechanism check and needs no golden file.

- `container.a.b` == `getattr(container.a, 'b')` == `container.a['b']`, and
  `hasattr` agrees with all three
- every name in `dataOrder()` is reachable by `getattr` --- this alone catches
  the BOXS bug
- a leaf reports no children; a container's `children()` matches its
  `dataOrder()` set
- round-trip: read, write back, read again, unchanged; `isSet()` transitions
- CLists after append / insert / delete / renumber: `_items[i]._name == f"[{i}]"`,
  `find_child("[i]")` finds it, `list[i]` and `children()[i]` agree, and no
  child exists that is not in `_items`

Three additions from today's evidence:

1. **Test before *and* after mutation.** Every bug today was a *stale* state ---
   a cache key that was right when written. A tree only ever built and read
   looks perfect.
2. **Hold the plugin, and also test not holding it.** A container outliving its
   plugin silently empties; make that an assertion rather than a trap.
3. **Assert agreement, not expected values.** Expected values need a golden file
   that rots.
