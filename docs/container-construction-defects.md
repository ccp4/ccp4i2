# Brief: two container-construction defects behind silent wrong answers

**Status:** brief. Investigation is done and the evidence is below; the fixes are
*not* written. Handed over 2026-08-27 for evaluation on the error-handling /
i2run strand.

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
