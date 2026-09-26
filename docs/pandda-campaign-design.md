# PanDDA over a campaign: an orchestrator task and a per-dataset receipt

**Status:** design note, 2026-09-23, revised 2026-09-24 against `django` at
`ca0ec99ff` (#392) after review. Nothing here is implemented; work proceeds on
the `pandda-campaign` branch (§14.7). Sections marked **OPEN** are decisions to take before
coding, not details to settle during it.

**Scope:** running PanDDA2 across a fragment campaign from inside CCP4i2, and
landing its results in the per-crystal projects as first-class data.

**Stance: desktop/laptop-local execution is the target.** CCP4i2 code exists to
serve a user running on their own machine. But PanDDA is the one task where
local will routinely *not fit* (§6), and a large-compute mechanism already
exists elsewhere (§5). So this design does two things at once: implement the
local run, and place the seams so that the same staged inputs and the same
results-landing work unchanged when the run happens somewhere bigger. It does
not implement the big-compute path, and does not make local execution
subordinate to it. Done right, the receipts it lands are what the campaign
views review — and a separate review application stops being necessary
(§5.6).

---

## 0. Provenance

| Source | What it gives this design |
|---|---|
| `lib/pandda_export.py` | Working knowledge of what PanDDA2 needs on disk, and the FreeR relabelling that makes it accept our MTZs |
| `pandda-inspect-api` ("reinspect") | A relational model of PanDDA output built against real runs — its *findings* are load-bearing here even though its architecture is not adopted |
| `materia`: `CCP4I2_PANDDA_INVOCATION_CONTRACT.md`, `PANDDA2_ON_AZURE.md`, `pandda-batch.bicep` | A specified, runner-agnostic invocation contract; measured memory behaviour; a working Azure Batch compute plane. §4 and §6 are adopted from it rather than re-derived |
| The CCP4 bundle itself (`ccp4-20260702`) | Ships PanDDA2 in its own micromamba env — the zero-install local execution story, and an older upstream snapshot whose gaps §4.6 enumerates |
| `django` merges #598–#605, reviewed 2026-09-24 | #603 landed the `{SUCCEEDED, UNSATISFACTORY}` glean gate §7.2 relies on; #600 gave `interactive.py` a direct test harness §5.5 builds on; #602 left `CDmDomain` and the registry untouched; #392 left the `core.ts` anchors §11 cites where they were. Nothing changed the design; two things confirmed it |
| Materia review, 2026-09-23 | A verification pass against this repo that found the composed-type registration requirement (§7.1), the glean-gate collision (§7.2), and the scheduling of §10.4. Its findings are incorporated, with two remedies changed — see those sections |
| `pipelines/MakeProjectsAndDoLigandPipeline` | A previous attempt at a job that creates jobs in other projects. It no longer runs; its four failure modes are the specification for §8 (Appendix A) |

---

## 1. The shape of the problem: computation and ownership have different units

One PanDDA run spans N crystals. A CCP4i2 job belongs to exactly one project.
So no single job can honestly own the per-dataset results, and a design that
tries will either put N datasets' outputs in one project or invent a parallel
store outside the job model.

That gives the architecture:

```
  fan-in               run                    fan-out
  (job construction)   (orchestrator task)    (management command)
  campaign-aware       campaign-blind         campaign-aware
  DB code              pure container         DB code
  local, cheap         LOCAL or ELSEWHERE     local, cheap
        │                     │                      │
        ▼                     ▼                      ▼
   DATASETS list ──▶ contract-conformant ──▶ one receipt job in each
   prefilled          input tree, run,         member project
                      output tree
```

**Two tasks, three steps.** The middle step is a task because it is a
computation with inputs and outputs. The outer two are not tasks, because they
are operations on the database — and keeping them out of the task layer is what
lets the task stay database-free (§3).

**The middle step is also the only one that might not run here.** Fan-in and
fan-out are cheap, local and database-shaped; they belong to CCP4i2 on the
user's machine regardless of where the compute happened. That asymmetry is the
whole of §5.

---

## 2. What exists today, and what we keep

| Thing | State | Disposition |
|---|---|---|
| `export_pandda` → ZIP + `Projects.csv` | Works; user runs PanDDA2 elsewhere. Already the input source for the Materia/Azure path | Becomes staging inside the orchestrator job (§3.2), keeping the same output shape |
| `prepare_mtz_for_pandda()` | Relabels FreeR to a PanDDA-accepted name | **Keep verbatim.** PanDDA2 accepts only `FREE`, `FreeR_flag`, `R-free-flags` |
| `collect_pandda_datasets()` | Finds each member's latest finished dimple + acedrg | Becomes the fan-in default, not the only source (§3.1) |
| `ProjectGroup` + `one_parent_per_group` | Campaign with exactly one reference project | The orchestrator job's home |
| `CampaignSite` / `SiteEvaluation` | Sites with origins in the parent frame; per-(project, site) human verdicts | The run-independent half of the event model (§9) |
| Campaign Moorhen page, site scenes, verdict chips | Works | The inspection UI. This design feeds it; it does not get a new one |

---

## 3. The orchestrator task is a pure function of a declared dataset list

### 3.1 It takes a list, not a campaign

`inputData` declares `DATASETS`, a `CList` of a composed type holding, per
dataset: `DTAG` (`CString`), `XYZIN` (`CPdbDataFile`), `HKLIN`
(`CObsDataFile`), `DICT` (`CDictDataFile`, optional).

The task therefore never asks the database anything. Consequences, all wanted:

- **`validity()` stays pure container logic** — "is the list non-empty, are the
  files typed and present" — so the standing rule that a task can run in a
  non-database-backed setting survives, and the meta-task exception we thought
  we needed is not needed.
- **It is runnable under `i2run`** against a hand-built list of three datasets,
  which is the only realistic way to test it.
- The list, in `params.xml`, is the **record of exactly what was submitted** —
  which matters for a run whose inputs may be edited underneath it over hours.
- **Batching is free** (§6.3): splitting a campaign into runs over disjoint
  subsets needs no new mechanism, because the subset *is* the parameter.

**Campaign-awareness lives in fan-in**: a campaign-aware endpoint creates the
job with `DATASETS` prefilled from `collect_pandda_datasets()`, exactly as the
campaign client already creates jobs via `projects/{id}/create_task/`. The user
sees and can edit the list before running. A dataset the default rule misses is
added by hand — a property of prefilling a real parameter rather than computing
the list at run time.

### 3.2 Staging produces the contract's input tree, exactly

PanDDA2's `--data_dirs` must name a directory containing *only* directories,
one per crystal, each holding a model and reflections matched by regex. **The
directory name is the dtag**, so a staging tree is not avoidable — it is how
dtags are assigned — and it is built inside the orchestrator's job directory.

We emit the shape the invocation contract already specifies, unchanged:

```
<staging>/
├── datasets/
│   ├── xtal-0000/{final.pdb, final.mtz, dict.cif?, compound/dict.cif?}
│   └── xtal-0001/…
└── Projects.csv          # header "Dataset, Project", one row per xtal-NNNN
```

(`compound/dict.cif` is a hardlink of `dict.cif`, for the CCP4-bundled
reader — §4.6 gap 6. The contract's tree is unchanged; this adds to it.)

**Names must be clean `xtal-NNNN`.** PanDDA infers a crystal number from the
last run of integers in the directory name, so a project-derived name with a
date in it parses as ~20 million and the range filter silently drops every
dataset. This is a known failure mode with a catalogue entry
(`dataset_range_zeroed`), not a hypothetical.

`Projects.csv` is kept **as well as** the richer manifest of §3.7 — it costs
nothing, it is what the existing Reinspect ingest reads to recover crystal
identity, and dropping it would break a working consumer for no gain.

Staging uses one helper, `link_or_copy()`: hardlink when source and destination
share a filesystem, copy otherwise. Not symlinks:

- symlink creation on Windows needs privilege or developer mode (§12);
- project export, move and restore would each have to decide whether to follow
  them;
- PanDDA2's own `-pandda-input.pdb`/`.mtz` are already symlinks into a sibling
  `data/`, so we would be linking to links.

To be fair to symlinks: they demonstrably work for the *PanDDA* read path —
the BAZ2B staged tree on `/Volumes/LocalStore` is built from them and runs
fine. The objections above are all CCP4i2-side, not PanDDA-side, and a hardlink
gives the same disk saving on the same filesystem without any of them.

Pure linking is impossible anyway: `prepare_mtz_for_pandda()` writes a *new*
MTZ whenever FreeR needs relabelling, so some entries are real files by
construction.

**Disk cost, measured rather than feared:** a staged dataset is 1.8 MB
(96 KB model + 1.6 MB reflections + 12 KB dictionary), so 120 datasets is
~216 MB — and only the relabelled MTZs arecopies; the rest hardlink to zero.
The tree lives in the orchestrator's job directory and is removed when that job
is deleted, which is the ordinary CCP4i2 lifecycle and needs no reaping policy
of its own.

The same helper is used again at the far end (§8.3). Bulk bytes move in exactly
two places, through one function.

### 3.3 Outputs

- The **run manifest**: a plain JSON `CDataFile` mapping dtag → project uuid →
  source job uuid → staged file digests, plus the run's provenance (§4.4).
  Keyed on **uuid, never name**: `Projects.csv` keys on project name, and names
  are renameable.
- The **output tree** location, so fan-out and Resume can find it.

No new `CDataFile` *subclass* for the manifest. The test for whether a subclass
is justified: **name the task that would take it as a typed input.** Fan-out is
not a task, so nothing would autopopulate it, and subType matching would be
decoration. §10.3 applies the same test to the event set.

---

## 4. The invocation contract is the interface, not the PanDDA CLI

Materia's `CCP4I2_PANDDA_INVOCATION_CONTRACT.md` already specifies, in
runner-agnostic terms, how PanDDA2 is invoked and what it produces. **We adopt
it rather than inventing a second invocation**, for a reason that is not
politeness: two of its arguments are defensive against failure modes that are
invisible until they have cost you a six-hour run.

### 4.1 The argv, and why

```bash
pandda2.analyse \
  --data_dirs        <staging>/datasets \
  --out_dir          <output> \
  --local_cpus       N \
  --pdb_regex        final.pdb \
  --mtz_regex        final.mtz \
  --ligand_cif_regex dict.cif \
  --ligand_pdb_regex ligand.pdb \
  --dataset_range    0-999999999
```

| Argument | Why it is what it is |
|---|---|
| `--ligand_pdb_regex ligand.pdb` | **Defensive.** PanDDA's default matches `final.pdb` and takes the *protein model* as a ligand, producing nonsense. We never emit `ligand.pdb`, so this literal never matches — which is the point |
| `--dataset_range 0-999999999` | **Defensive.** Guards the integer-parsing behaviour described in §3.2 against any tree not staged by us |
| `--local_cpus N` | Real, but it changes *concurrency only*. It does **not** lower peak RAM — see §6 |
| the three other regexes | Match what §3.2 stages |

### 4.2 Environment: `RAY_TMPDIR` is not optional

PanDDA runs Ray, which spills to `RAY_TMPDIR`. The image default is `/tmp/ray`,
which the contract explicitly marks *do not use in production*: filling the OS
root is a known SIGKILL trigger with its own catalogue entry
(`ray_scratch_full`). It needs **tens of GB on a fast local disk** — the Batch
pool picks an SKU with ~600 GB of local NVMe specifically for this.

On a laptop this is the setting most likely to be wrong, so: a declared
`SCRATCH_DIR` control parameter, defaulting into the job directory, with a
`runTimeValidity()` check that it exists and has headroom. A run that dies at
hour four because `/tmp` filled is the worst failure this task can have, and
it is entirely preventable at submit time.

### 4.3 The output tree shape is contract, not convention

PanDDA writes `pandda2_out/` under `--out_dir`, containing `analyses/`,
`processed_datasets/<dtag>/` and `pandda_log.txt`. **The `pandda2_out/` name is
part of the contract; do not rename it.** Fan-out reads this shape, and so does
Reinspect's ingest — which is what lets both consume the same tree (§5.3).

### 4.4 Record what produced the tree

The contract pins runs to an explicit image tag so failure classification stays
version-aligned. The same discipline applies locally: the job records the
PanDDA version / fork ref / image tag and the contract version in its params.
Without it, a failure catalogue entry cannot be trusted to apply, and a
two-year-old receipt cannot be interpreted.

There is a second reason, stronger than hygiene: **built poses are
path-dependent.** PanDDA's pose search is stochastic — which is why it has a
seed parameter at all — so two runs over identical data can build the same
event differently, plausibly equally well but not identically. A receipt whose
provenance does not record the build path and seed cannot be compared with
another: any difference between them is unattributable between the data and the
search. So the recorded provenance is what makes two receipts commensurable,
not merely traceable.

This is also the second argument for §7.3's rule that the **apo model is the
model of record and poses are candidates**. Treating a pose as the answer is
treating one draw as the distribution, and it is why nothing in this design
merges a pose automatically.

### 4.5 Progress and failure are contract symbols, not log scraping

- **Progress**: PanDDA emits `PANDDA_PROGRESS: dataset <j>/<N>`, matched by
  `^PANDDA_PROGRESS: dataset (\d+)/(\d+)$`. This drives the running report
  (`runningReport=True` + `watchedFile` in the `TASKS` entry). Absence of the
  line means *progress unknown*, never failure — older builds lack the patch.
- **Failure**: classify by **stderr regex against the catalogue**, never by
  exit-code value; the exit code carries only the success/failure split. The
  seed catalogue (`oom`, `free_r_label`, `ligand_block`,
  `dataset_range_zeroed`, `ray_scratch_full`, `ccp4_missing`) maps onto
  `ERROR_CODES` with recovery prompts. Adopt it; do not re-derive it.

### 4.6 Which PanDDA — there are three, and they disagree

| Where | What it is | Provenance |
|---|---|---|
| `$CCP4/bin/pandda2.analyse` | PanDDA**2** in a micromamba env (`$CCP4/share/mamba/envs/pandda2`, Python 3.9), CNN checkpoints bundled. A shim: `micromamba run -r $CCP4/share/mamba -n pandda2 pandda2.analyse` | An upstream `ConorFWild` snapshot **predating 2026-05-07** — it has `cnn/resnet.py`, which upstream renamed to `resnets.py` on that date. Version string is `0.0.1` and `direct_url.json` records only a Jenkins build path, so **there is no commit recorded and no version to test** |
| `$CCP4/bin/pandda.analyse`, `pandda.inspect` | PanDDA**1** (`panddas-1.0.0`, under ccp4-python 3.11) | A different program with a different output tree. Not a target here |
| `ConorFWild/pandda_2_gemmi:master` | **Current upstream.** HEAD `b18389e6`, 2026-06-17 | Carries every fix below: PRs **#96** (free-R by column type), **#97** (ligand comp block by content), **#98** (`PANDDA_PROGRESS`), **#99** (ligand dict bonds via gemmi `ChemComp`) are all **merged** |
| `martinemnoble1/pandda_2_gemmi:crowther-local-autobuild` | A personal branch, not upstreamed | Upstream plus the unvalidated memory/fit experiment of §6.5 — the only thing on it that is not now in upstream |

**The CCP4 build is what users will have, so it is the build this design
targets.** Zero install, no conda, no container, no ~200 MB checkpoint
download.

The important thing about the gaps below is that **they are not a fork-versus-
upstream question.** Every one of the four fixes is merged into upstream
`master`; the CCP4 build simply predates them. So this is a staleness problem
with a one-line answer (§4.8), and in the meantime four of the five gaps are
things CCP4i2 can handle on its own side.

The five gaps, one of which bites *precisely because the inputs come from
CCP4i2*:

1. **Ligand-dictionary bond orders — silent when it bites.** The CCP4 build
   reads bond order only from `_chem_comp_bond.type`
   (`autobuild/inbuilt.py:139`). A dictionary that instead names that column
   `_chem_comp_bond.value_order` (with `pdbx_aromatic_flag`) makes the bond
   loop return empty, the `zip()` over bond columns run zero times, and the
   molecule get built **with no bonds** — a bag of disconnected atoms that
   RDKit's conformer generation stacks on top of each other. Collapsed poses,
   no error. Upstream fixed it in PR #99 by reading through gemmi `ChemComp`
   (`dataset/small.py`), which understands both spellings.

   **Which producers emit which spelling is version-dependent, and must be
   checked against the acedrg in the CCP4 we target.** Sampled evidence from
   this machine: acedrg/CCP4-monomer-style dictionaries (including the BAZ2B
   campaign dictionaries) carry `type`; wwPDB/PDBx component files carry
   `value_order`. So the CCP4 build is *not* broken on every CCP4i2 dictionary
   today — but it is one acedrg output-format change away from being so, and it
   is already broken on any PDBx-derived dictionary a user supplies. That is
   what §4.7 insures against, cheaply, rather than depending on which way the
   version lottery lands.
2. **FreeR detection** by label only, raising `No RFree Flag found!` on a
   mismatch; the fork detects by MTZ column *type*. We are protected either
   way, because §3.2's staging already relabels — worth noting that **our own
   staging is what compensates**, which is the pattern §4.7 generalises.
3. **`--dataset_range` defaults to `0-99999`** (the fork made it `None`). So
   the contract's explicit `--dataset_range 0-999999999` is **required**
   against the CCP4 build, not merely defensive: any dtag whose trailing
   integer exceeds 99999 is silently dropped.
4. **No `PANDDA_PROGRESS`**, so no running-report progress. Per §4.5 absence
   means *unknown*, so this degrades rather than breaks.
5. **No bounded-memory autobuild path**, so §6's memory story is the
   unmitigated one. This is the **only one of the five we cannot neutralise**
   (§6.5).
6. **Ligand files are read only from a `compound/` subdirectory** of each
   dataset (`fs/pandda_input.py`, `ligand_dir_regex`); a `dict.cif` beside
   the model, which is what the contract specifies and what upstream reads
   first, is invisible to it: *"Filtered because no ligand data!"* for every
   dataset. Found on the first real run (2026-09-24). Staging now puts the
   dictionary in both places, the second a hardlink, and the argv names
   `--ligand_dir_regex compound` explicitly. One ligand either way.
7. **The dictionary block is opened by name, `comp_LIG`** (fallback
   `comp_XXX`; `dataset/small.py`). acedrg names the block after the
   ligand — `comp_MZ0` — so every CCP4i2 dictionary fails here; upstream
   finds the block by content (PR #97, the `ligand_block` catalogue entry).
   `prepare_dict_for_pandda()` now appends a `comp_LIG` copy of the
   restraint block after the original, so the name-based reader finds it
   and the content-based reader still takes the true residue first. Nothing
   is renamed.

Tally: gaps 1, 6 and 7 are removed by §4.7 and staging, gap 2 is already
removed by staging, gap 3 is removed by the contract's argv, gap 4 degrades
gracefully. **So CCP4i2 needs no
PanDDA fork and no particular PanDDA version** — which is the property that
makes this task shippable. A CCP4 rebuild from current upstream (§4.8) removes
gaps 1–4 at source and costs us nothing either way.

**Resolution order** — discover, do not hard-code, using the existing
[binary-discovery mechanism](PREFERENCES_BINARY_DISCOVERY_PLAN.md):

1. an explicit preference or control parameter;
2. `pandda2.analyse` on `PATH` (a user's own env or the fork);
3. `$CCP4/bin/pandda2.analyse`.

And record the resolved path in params (§4.4). Since both candidates report
version `0.0.1` and carry no commit, **the resolved path plus a capability
probe is the only usable provenance there is.** Probe for behaviour, never for
a version number.

### 4.7 Fix what we can at our end: `prepare_dict_for_pandda()`

`prepare_mtz_for_pandda()` already exists and is exactly the right precedent:
when PanDDA cannot read something we produced, **normalise it at staging rather
than requiring a different PanDDA.** We control the staging code; we do not
control which build the user has.

So the dictionary problem in §4.6.1 is solvable on our side, and the right shape
is a **normaliser, not a converter**: read the dictionary with gemmi `ChemComp`
and emit a `dict.cif` that carries **both** spellings consistently, whichever
the producer used. The old reader maps `_chem_comp_bond.type` through
`{single, double, triple, SINGLE, DOUBLE, TRIPLE, aromatic, deloc}`, so a
`value_order` dictionary gains an equivalent `type` column
(`pdbx_aromatic_flag == 'Y'` → `aromatic`); a `type`-only dictionary is left
readable by the new reader, which already handles it.

Written that way it is correct for every combination of producer and PanDDA
version, including the ones that do not exist yet — which is the property worth
having, given §4.6.1 shows the spelling is a version lottery. And it lives next
to the FreeR relabelling it mirrors.

That turns a version incompatibility into a staging concern, which is where it
belongs: the task then works with whatever PanDDA the user has, and gets better
rather than *working* when they upgrade.

### 4.8 Two asks on the CCP4 bundle, neither blocking

**Ask 1: rebuild the bundled `pandda2` from current upstream `master`.** The
bundle's snapshot predates 2026-05-07; upstream HEAD is 2026-06-17 and contains
all four fixes. Nothing about the packaging changes — same micromamba env, same
bundled checkpoints — so this is a version bump, and it removes gaps 1–4 of
§4.6 including the silent ligand-bond one. Worth raising with whoever owns
`devtools/co/pandda2` in the CCP4 build.

**Ask 2 (longer term): the Python 3.9 pin is one dependency deep.** It reads as
a PanDDA constraint and is repeated as one, but it is not:

- `pyproject.toml` declares `requires-python = ">=3.9"` — a floor, not a pin.
- `requirements.txt` is entirely commented out, so `setup.py` contributes no
  constraint at all.
- `environment.yml` pins `python=3.9` with the comment *"need to pin to avoid
  issues with builds"* — a 2022 build-convenience decision.
- The real constraint is **`torch==1.13.1`** (December 2022), whose last wheels
  are cp310. That is also what forces `numpy<2.0`, since it predates the NumPy
  2.0 ABI break.

Everything *else* in the CCP4 env has already moved on — `pytorch_lightning
2.2.5`, `ray 2.51.2`, `rdkit 2025.9.2`, `gemmi 0.7.5`, `numpy 1.26.4`. Torch is
the lone straggler, and PanDDA's torch surface is inference-only and tiny:
`torch.nn`, `torch.from_numpy`, `torch.ones`, `torch.flatten`, and a single
`torch.load`.

So relaxing the pin is a bounded job whose one real hazard is that single
`torch.load`: `torch` ≥2.6 flipped its default to `weights_only=True`, which a
Lightning checkpoint will not survive without an explicit flag or Lightning's
own loader. The test is correspondingly concrete — build a Python 3.11 env with
current torch, load both checkpoints, and compare per-event scores against the
3.9 baseline on a small public dataset.

**Why CCP4i2 should care:** a `torch==1.13.1` environment cannot be installed
into `ccp4-python` (3.11), and that is the entire reason PanDDA2 needs a second
interpreter inside CCP4 at all. Relaxing the pin is what would let it live in
`ccp4-python` — one interpreter, no micromamba shim, and PanDDA becomes an
ordinary CCP4 dependency.

**Both asks are CCP4-bundle work, not CCP4i2 work.** They are recorded here
because this is where the evidence was gathered, and they are out of scope for
every work package in §14. Note also that ask 2 is wider than one `torch.load`
flag: `numpy<2.0` travels with the torch pin, so relaxing it is cross-cutting
through the whole dependency set rather than a single call site.

---

## 5. One contract, two deployments — and CCP4i2 owns the ends, not the middle

### 5.1 What already exists, and where CCP4i2 already sits in it

Materia's large-compute path is built and decided (2026-06-06): a Materia
frontend button calls **CCP4i2's `export_pandda`** to stage the tree, POSTs a
run to Reinspect, whose `PANDDA_JOB_RUNNER` factory (`LocalProcessRunner` |
`AzureBatchRunner`) submits to an Azure Batch pool that scales 0→1 on a
128–256 GiB node, and whose ingest reads `pandda2_out/` in place. Reinspect
owns that run lifecycle: the Batch SDK, the run-status model, sizing policy,
retry semantics and the merge-vs-replace ingest rules.

Note the first step: **CCP4i2 is already the stager.** This design does not
insert CCP4i2 into that flow; it makes the part CCP4i2 already does a proper
job rather than an endpoint side effect.

### 5.2 The division that follows

| Step | Who | Where it runs |
|---|---|---|
| Fan-in + staging | CCP4i2 | Always local |
| The run | a runner | Local (implemented here) **or** elsewhere (exists already) |
| Fan-out into member projects | CCP4i2 | Always local |
| Triage and review | The campaign Moorhen views, consuming receipts (v2, §5.6); Reinspect where it is already deployed, until then | Wherever the projects are |

CCP4i2 owns the ends because they are database work about *its* projects.
CCP4i2 does not own the middle in every deployment, and should not try to: it
would mean reimplementing the Batch SDK, the sizing catalogue and the failure
taxonomy that already exist and are already maintained.

> **Reversed 2026-09-26 (decision 18).** The argument above was conditional on
> Reinspect existing to own the middle. Materia intends to deprecate Reinspect
> in favour of a CCP4i2-based launch mechanism plus CCP4i2's own campaign
> views, so the premise goes with it. CCP4i2 now owns the *contract* for the
> middle — a run-target registry with two narrow interfaces — while each
> deployment registers its implementation. What CCP4i2 still does not own is
> any platform: no Batch SDK, no Service Bus, no sizing catalogue beyond the
> `sizing_hint` it already computes.

### 5.3 The two consumers do not conflict — and one of them is transitional

Reinspect's ingest and CCP4i2's fan-out read the same `pandda2_out/` tree and
build different things — a cross-dataset triage model, and per-crystal receipts
in per-crystal projects. Neither writes into the tree. A deployment that has
both gets both; a desktop has only the second.

But this is a transitional state, not a division of labour. §5.6 says why.

### 5.4 What "preempting" costs us here: two requirements, no code

The big-compute path stays available without building any of it, provided:

1. **§3.2's staged tree is contract-conformant.** It is — that is why §3.2
   emits the contract's exact shape including `Projects.csv`, rather than a
   shape of our own choosing.
2. **Fan-out accepts any conformant output tree**, whoever produced it (§8.4).

Both are free. Neither requires an Azure dependency, a runner abstraction with
one implementation, or any code in CCP4i2 that talks to Batch.

### 5.5 The escape hatch: a run that happened elsewhere

Given §6, the common case for a real campaign on a laptop is *it will not fit*.
The graceful answer is an **external-run mode**: the orchestrator job stages the
tree, reports that the run must happen elsewhere and where the tree is, and
finishes in a state that can later be handed the resulting `pandda2_out/` —
after which fan-out proceeds normally.

The job then still exists in CCP4i2 with its `DATASETS` list, its manifest, its
provenance and its receipts. The precedent is the recorded Moorhen task: a job
whose work happens outside the plugin process and is recorded on return.

**DECIDED: v1.** Two arguments, the second decisive.

First, §6 spends its longest section establishing that local will routinely not
fit — a 50-dataset campaign OOMing 128 GiB is not a marginal case. A v1 that
implements only local execution therefore serves the case the document itself
argues is the minority.

Second, it is cheaper than this section first assumed, because **the lifecycle
already exists**. `Task.interactive` ([tasks.py:38](../server/ccp4i2/core/tasks.py))
already means *Run does not dispatch; the job is finished later by an external
event*, and `lib/utils/jobs/interactive.py` provides `open_session`,
`session_state`, `finish_session` and `cancel_session` against existing schema
— no migration — with a direct test harness in
`tests/db/test_interactive_session.py` (#600) that drives every disposition
without a thread or a window.

Be precise about what is reused: **the session lifecycle, not the drop
mechanism.** `drop_file` accepts only `kind in ("model", "dictionary")` and
writes the `output<N>.<ext>` contract that GUI tasks harvest; a `pandda2_out/`
tree is neither. External-run mode finishes a session by recording *where the
tree is* — a path the receipt reader and fan-out then consume exactly as they
would a locally produced one (§8.4). So the new code is the finish-with-a-path
step and its validation, on top of a lifecycle that already has tests.

### 5.6 Done right, the receipts make a separate review application unnecessary

Reinspect exists because PanDDA's output tree was not first-class data
anywhere: it had to be ingested into a relational model, served through a
`DataStore` seam, reviewed in a bespoke Moorhen client, and its poses refined
through a hand-built job service that shells out to refmac. Every one of those
is something CCP4i2 already has natively, once the output tree becomes
receipts:

| Reinspect built | CCP4i2 has |
|---|---|
| Ingest into `Dataset`/`Event`/`Artifact` rows | Fan-out into typed receipts, gleaned to `File` rows (§7, §8) |
| `Finding` vs `Observation`, moved-peak guard | `SiteEvaluation` vs receipt events, §9 |
| Artifact serving via `source_root` | Ordinary job-file serving |
| Merge pose onto apo | "Place the ligand here", generalised to *accept this pose* |
| Per-event refinement via `jobservice` | `SubstituteLigand` / `servalcat_pipe`, with dictionaries, provenance and export already correct |
| Decision persistence | `SiteEvaluation`, backed by the snapshot once §14.0 lands |

So the strategic shape is: **the campaign Moorhen views expand to consume the
receipts**, and Reinspect's review role is absorbed. Its *findings* — the
`events.yaml` shape, the BDC token, absolute contour units, the apo start
model, unstable event ordinals — are already lifted into this document. Its
*code* does not need to survive.

Two honest limits on that claim:

- **The cross-dataset triage queue is new.** The campaign views today are
  per-dataset. "Every event across 200 datasets, ranked, step through them
  fast" is the view Reinspect is actually *for*, and it is exactly the query
  §10.2's `CampaignEvent` projection exists to answer. That deferral therefore
  flips: **`CampaignEvent` is v2-essential, not optional** (§10.2, §14.3).
- **Batch submission is the one Reinspect function that does not fold into
  "views consume receipts".** Where a deployment runs PanDDA on Azure Batch,
  *something* submits; today that is Reinspect's `AzureBatchRunner`. It is a
  small piece — a couple of hundred lines of SDK — and the natural home for it
  is a runner behind the orchestrator (§5.4's seam), but that is Materia's
  deployment to decide and is not in this document's scope.

Sequencing matters for how this is said to Materia: v1 lands receipts; v2
lands the triage view over `CampaignEvent` and *accept-this-pose*. At that
point Reinspect's review UX is redundant for any deployment running CCP4i2's
campaign views, and the earlier question of "who owns the human verdict" is
answered in the best way — there is one store.

---

## 6. Memory is the binding constraint, and on the shipped build nothing moves it

This is the fact that most shapes the local story, and it is measured, not
estimated.

### 6.1 The numbers

Peak RAM ≈ (comparators in a shell) × (map size) × (parallel workers), plus a
load-everything baseline. Map size scales with unit-cell volume: ~22–25 MB per
map for a long-axis cell (185–197 Å) against ~4 MB for a small bromodomain —
a 5–6× difference. Comparator count scales with dataset count.

Observed on ~120 large-cell datasets: at `--local_cpus 6`, swap ran 14→27 GB
and the OS killed it at ~3.5 min; at `--local_cpus 2` it still climbed
11→34 GB and was killed at ~5 min. **A single shell of that job intrinsically
needs 30–40 GB+.**

A production data point is sharper still, and worse: a **50**-dataset
CDK4/CyclinD1 campaign (58 × 64 × 186 Å) with autobuild is on record **OOMing a
128 GiB node**, and completing in 2 h 29 min on 256 GiB. So "≥128 GB for ~120
datasets" understates it — 50 datasets on a long-axis cell can exceed 128 GB.
Dataset count is not the variable to reason from on its own; cell volume moves
it more.

A small cell with ≤60 datasets fits in 16–32 GB, and that remains the laptop
case.

### 6.2 Every *CLI* flag that looks like a lever is not one

Verified against the code: `--memory_availability` and `--low_memory` are
stored and **never read**; `--grid_spacing` is not used in the crystallographic
path; `--sample_rate` is overridden by a hardcoded `resolution/0.4999`;
`--max_shell_datasets` is a *minimum floor* on comparators despite its name.
`--local_cpus` only changes how fast the wall is hit.

Do not offer these as tuning parameters in the task interface. Offering a
control that does nothing is worse than offering none.

### 6.3 Fewer datasets per run is a lever, and §3.1 already gives it

Batching a campaign into runs over disjoint subsets needs no new mechanism: the
subset is the `DATASETS` parameter, so it is N orchestrator jobs, each with its
own receipts, and §8.4's "latest finished receipt" rule already handles members
appearing in more than one.

**But batching is a compromise, not a scaling strategy**, and the task
interface should say so: fewer datasets per run means fewer comparators, which
degrades the ground-state model PanDDA's statistics rest on. It is a way to
make a run fit, not a way to make it better. PanDDA is single-node by design;
there is no distributed mode to grow into.

### 6.4 `runTimeValidity()` earns its keep here

The two-tier split lands perfectly: `validity()` stays cheap and pure (§3.1),
while `runTimeValidity()` — which may read files — estimates peak memory from
dataset count and cell volumes (gemmi, from the staged MTZs) and warns when the
machine will not fit, *before* six hours are spent finding out.

That estimate is the **same `sizing_hint: {datasets, cell_volume_class}`** the
run contract already defines. Compute it once: locally it is a warning, and if
the run is delegated it is the payload field. One calculation, two uses,
defined by someone else's already-working contract.

Severity: `SEVERITY_WARNING`, not error. The user may know something we do not,
and a blocked Confirm on an estimate is the wrong trade.

### 6.5 An experimental bounded-memory path exists; this design does not rest on it

For the record, because it will otherwise be rediscovered: a personal branch
(`martinemnoble1/pandda_2_gemmi:crowther-local-autobuild`) carries an
`PANDDA_LOCAL_AUTOBUILD` env switch that cuts local density boxes around each
event instead of unmasking the full cell, which would make autobuild-stage
memory independent of unit-cell size. A sibling `PANDDA_CROWTHER_FIT` replaces
the pose search outright.

**Neither is upstreamed and neither is validated end-to-end** — they are the
only part of that branch not now in upstream `master` (§4.6). The local-cut is
unit-tested against the full-cell cut and the fit work reproduces a brute-force
search on one public dataset, but there is no completed pipeline comparison. So:

- **This design does not depend on them, does not enable them, and does not
  expose bespoke controls for them.** §6.3's batching is the memory lever a
  shipped task can honestly offer.
- If a generic "extra environment for the PanDDA process" escape hatch is
  provided at all, it is an unlabelled pass-through, and the documentation
  records one trap: these switches are **presence-checked, not value-checked**,
  so `PANDDA_LOCAL_AUTOBUILD=0` is **ON**. Anything mapping a `CBoolean` onto
  them inverts its own meaning for `False`; the variable must be *absent* when
  off. (The `CBoolean` pitfall in mirror image.)
- If the work is ever validated and upstreamed, §6 gets better and §6.6's third
  row changes. Until then, treating it as available would be promising the user
  something we cannot support.

### 6.6 So what does "slow but I want to run it on my laptop" actually mean?

| Shape | Verdict |
|---|---|
| Small cell, ≤60 datasets, 16–32 GB | Runs as-is, on the CCP4 build. This is the supported laptop case |
| Large cell, modest dataset count | Autobuild-stage memory scales with the cell and there is no supported lever (§6.5). Batch it down, or run it elsewhere |
| Large cell, ~120 datasets | Comparator loading alone wants 30–40 GB+. Batch it (§6.3), accepting the statistical cost, or run it elsewhere (§5.5) |

The design's job is to make each of those *legible before the run starts*, not
to promise the third. Concretely: `runTimeValidity()` estimates and warns
(§6.4), `--local_cpus` defaults low rather than to core count, `SCRATCH_DIR` is
checked for headroom (§4.2), and the resolved PanDDA's capabilities are reported
so the size of the job is something the user is told rather than something they
discover at hour four.

Slow is an acceptable answer. Being killed at hour four is not, and almost every
instance of it is predictable at submit time.

---

## 7. The receipt task carries no logic, but it does carry a declared shape

One job per dataset per run, in the member project — a *receipt* for that
dataset's share of one PanDDA run.

### 7.1 Declared, not a bucket

`outputData` declares:

| Field | Class | Notes |
|---|---|---|
| `XYZIN_APO` | `CPdbDataFile` | The `-pandda-input.pdb`. **The model of record** — §7.3 |
| `ZMAP` | `CMapDataFile` | |
| `EVENTS` | `CList` of a composed type | One item per event |
| `PANDDA_MODEL` | `CPdbDataFile`, optional | The merged `pandda-model.pdb`, as a machine opinion only |

Each `EVENTS` item holds `EVENT_IDX` (`CInt`), `BDC`, `RSCC`, `BUILD_SCORE`,
`HIT_PROBABILITY`, `OPTIMAL_CONTOUR` (floats), `CENTROID`, a site reference,
and — the point — `EVENT_MAP` (`CMapDataFile`) and `POSE` (`CPdbDataFile`) **as
members of the composed type**.

#### The composed type must be a registered core class

**This is the one place the design costs more than it looks.** A `<subItem>`'s
`<className>` is resolved against a registry built by importing a **hard-coded
list of `ccp4i2.core.*` modules**
([def_xml_handler.py:54-69](../server/ccp4i2/core/task_manager/def_xml_handler.py#L54-L69)).
There is no plugin-side class discovery and no inline nested content in a
def.xml. So `CPanddaEvent` cannot live in the plugin directory: it needs a
`ccp4i2/core/CPanddaEvent.py` plus one line in `implementation_modules`.

> **Superseded 2026-09-26 (decision 18, item 15).** A task now declares the
> modules of its own CData classes in `Task.dataTypes`, and one registry
> (`core/cdata_registry.py`) serves both the def.xml handler and the file
> digest. `CPanddaEvent` lives in `wrappers/pandda_events/script/
> pandda_events_types.py` and `CPanddaDataset` in `wrappers/pandda_campaign/
> script/pandda_campaign_types.py`; `core/` carries nothing PanDDA-specific.
> The v1 rows below record what was true when v1 landed.

**And the failure mode is silent.** An unresolvable class name does not raise —
it logs `Warning: Unknown class '<name>', using CString as fallback`
([def_xml_handler.py:403](../server/ccp4i2/core/task_manager/def_xml_handler.py#L403))
and the task parses cleanly with a string where the composed type should be.
Getting the registration wrong produces something that looks like it works.

**There is a precedent, and it is cheap.** `CDmDomain`
([core/CDmDomain.py](../server/ccp4i2/core/CDmDomain.py), 72 lines) is exactly
this shape — a task-specific composed `CData` used as a `CList` subItem in
`dm_multidomain.def.xml:105`, with a docstring that says in as many words
"Resolvable by the def.xml class-name lookup via `ccp4i2.core.CDmDomain`". So
the pattern is established and the cost is one small module plus one line.

**What has no precedent is narrower: a composed subItem whose fields are
`CDataFile`s.** `CDmDomain` holds only strings. Nothing in the tree exercises
nested files in a list end to end, which is why §14.2 opens with a spike.

The parts of the path that *do* check out: `find_all_files()`
([ccontainer.py:746](../server/ccp4i2/core/base_object/ccontainer.py#L746))
recurses through `children()`; `glean_job_files` derives `param_name` from
`objectPath()` and already handles list indices; and `job_param_name` is a
`CharField(255)`, so `EVENTS[3].EVENT_MAP` fits with nothing downstream parsing
it.

Why declared rather than a bucket: outputs being *typed* is what drives glean,
autopopulation, mimeType-driven UI, Moorhen loading and Export MTZ. An
undeclared directory of dropped files is filesystem-shaped thinking with extra
steps.

### 7.2 Its entire logic budget — and the status it returns when short

Check that what was delivered matches what was declared, and **say so loudly
when it does not**. Over hundreds of datasets the silent partial is the failure
mode that costs weeks.

But "loudly" must not mean `FAILED`, because glean runs only for
`GLEANING_PLUGIN_STATUSES`
([async_db_handler.py:35](../server/ccp4i2/db/async_db_handler.py#L35)) — a
failed receipt would publish *nothing*, losing four good events to make a point
about a missing fifth. That would contradict at dataset level the argument §8.4
makes at run level.

The mechanism already exists and needs no invention — it landed the day
before this section was written (#603, 2026-09-23, "a job that ran to the end
publishes what it made"). That frozenset is `{SUCCEEDED, UNSATISFACTORY}`, and
the comment above it describes this case exactly — *"UNSATISFACTORY is 'I got there, but look at the report', and its
partial outputs are often the whole point"*. So **a short receipt returns
`UNSATISFACTORY`**: everything that arrived is gleaned and usable, and the job
is unmistakably not clean in the UI, in its report, and in a KPI carrying the
expected-versus-delivered counts. `FAILED` is reserved for a receipt that can
make no sense of the tree at all.

### 7.3 Two facts that must not be rediscovered

Both were paid for once already, in reinspect:

- **The apo `-pandda-input.pdb` is the start model.** Every event is a
  *candidate* pose merged onto apo. The merged `pandda-model.pdb` is the
  machine's opinion and must be unmistakably separate, or something picking
  "the structure" picks the wrong one.
- **`Optimal Contour` is in absolute map units, not σ.** PanDDA computes it as
  a threshold on raw BDC-corrected sample values. Moorhen's stored
  `contourLevel` is also absolute; the scene format's is **rmsd-relative**
  ([core.ts:316](../client/renderer/lib/scene/core.ts#L316)). §11 says where the
  conversion belongs.

### 7.4 Event index is a field, never a position

Carry `EVENT_IDX` explicitly; never rely on list position, given the known
nondeterminism in `children()` ordering. And from reinspect: `(dtag,
event_idx)` is an **unstable per-run ordinal** — across two runs, A's event 1
can be B's event 2. It identifies an event *within* a run and nowhere else.
Cross-run identity is §9's job.

---

## 8. Fan-out is a separate, re-runnable step

Not the tail of the orchestrator's `process()`. Appendix A is the evidence.
(A separate *task*, since 2026-09-24 — decision 1 as amended; the sections
below describe the behaviour, which the task and the command share.)

### 8.1 It calls the same entry points the UI calls

`projects/{id}/create_task/` then `jobs/{id}/run/` — never a parallel
implementation of job creation. The previous attempt hand-rolled ~120 lines of
framework internals and rotted when the framework moved (A.2).

### 8.2 Four properties, each answering a way the previous attempt failed

| Property | Shape |
|---|---|
| **Idempotent** | Keyed on `(orchestrator job uuid, dtag)`. Re-running skips what landed |
| **Previewable** | A dry run reporting what it *would* create |
| **Reported** | A result table per dtag: created / skipped / failed, with the reason |
| **Retryable** | Failures are addressable individually, because they are recorded |

### 8.3 Delivery copies bytes into the member project

Event maps and poses are `link_or_copy()`d from the output tree into the
receipt job's directory, not referenced in place. This costs disk — order a few
MB per event — and buys the property §10 depends on: **a project directory is
self-contained**, so export, move and restore work without reaching into
another project. Given the doctrine in §10, referencing in place is not an
option; hardlinking makes it nearly free on one filesystem.

### 8.4 It takes a tree and a manifest, not a runner

Fan-out's input is a conformant `pandda2_out/` plus the manifest that maps
dtags back to projects. **It does not care who produced the tree** — local run,
HPC, Azure Batch, or a colleague's USB stick. This is requirement 2 of §5.4,
and it is the single property that keeps the big-compute path open at zero
cost.

It also means fan-out can be pointed at a **partial** tree. PanDDA can write a
hundred processed datasets and then OOM; glean is gated on SUCCEEDED, so the
orchestrator job publishes nothing, which is right. But the hundred datasets
are real and cost six hours. Because fan-out is a separate operator-initiated
step, it can be run against a failed run's tree deliberately, with the receipts
recording that they came from an incomplete run. That is a direct dividend of
the separation, not an afterthought.

### 8.5 A rerun is a new receipt, and that is the feature

A second run produces a second orchestrator job and a second receipt per
dataset. Both stay visible. Nothing collapses or updates in place.

This is reinspect's multi-run problem solved by the job model rather than by a
reconciliation policy — *provided* nothing keys off "the" PanDDA receipt in a
project. Every consumer asks for the **latest finished** one, exactly as
`_latest_finished_job()` already does for dimple.

---

## 9. Sites land on `CampaignSite`, and verdicts survive reruns

Reinspect had to invent a run-independent `Finding` distinct from a run-scoped
observation, because re-ingesting clobbered human decisions. CCP4i2 already has
both halves:

- `CampaignSite` — a place, with identity and an origin in the parent's frame;
- `SiteEvaluation` — one person's verdict on one dataset at one site, where
  *absence of a row* means nobody looked.

So: PanDDA sites merge into `CampaignSite`; events reference the site they fall
in; human verdicts sit on `(project, site)` and **a rerun cannot touch them.**

Two rules, both from reinspect's scars:

- **Match sites by centroid proximity in the parent frame, never by name.**
  Names were what orphaned references before migration 0024.
- **Derive centroids from member-event coordinates, not from the
  `pandda_analyse_sites.csv` column**, which is frequently `(0,0,0)`.

Worth borrowing too: reinspect's **moved-peak guard** — before honouring any
match that would re-attribute a human call, cross-check the peak coordinate, so
a renumber cannot silently move a verdict onto different density.

### 9.1 Marrying the two models: what reconciliation needs, and what §14.0 must not preclude

§14.0 is about to serialise `CampaignSite` and `SiteEvaluation`. Before it
does, the question is whether the longer-term goal — reconciling PanDDA's
sites and scores with our own — wants those tables shaped differently. The
answer is that the *shapes* are right, and two *invariants* and one
*serialiser property* are what §14.0 must carry so v2 never has to reopen it.

**The entities line up, once machine and human are kept apart.**

| Concept | PanDDA / Reinspect | CCP4i2 | Nature |
|---|---|---|---|
| A place in the campaign | PanDDA site (per-run ordinal, unreliable centroid); Reinspect `Site` | `CampaignSite` — origin in the parent frame, stable identity | run-independent |
| A machine observation at (dataset, place, run) | PanDDA event; Reinspect `Observation` | the receipt's `EVENTS` item (§7.1); optionally `CampaignEvent` (§10.2) | run-scoped |
| A machine opinion | `hit_in_site_probability`, `interesting`, build/RSCC scores | KPIs on the receipt; `detect_ligands()` on refined coordinates | run-scoped, never a verdict |
| A human decision at (dataset, place) | Reinspect `Finding.decision` | `SiteEvaluation` — absence means nobody looked | run-independent |

Reinspect learned to keep the second and fourth rows distinct because
re-ingest clobbered decisions. Our model already separates them *structurally*:
machine observations live in receipts (one per run, per job) and human verdicts
live in `SiteEvaluation` (per dataset, per site). Nothing needs to change for
that to hold — but it must be *said*, because it is the one thing that makes
the marriage safe:

> **Invariant 1: nothing automated ever writes a `SiteEvaluation` row.**
> PanDDA proposes sites and records observations. Only a person records a
> verdict. `interesting = True` is an opinion, not a decision.

**Sites: PanDDA proposes, `CampaignSite` disposes.** A PanDDA site either
matches an existing `CampaignSite` by centroid proximity in the parent frame
or becomes a new one. Once matched or created, it *is* a `CampaignSite`, with
the same identity as one a human placed by hand. That gives a second
invariant, which is what makes reruns safe:

> **Invariant 2: a rerun may propose sites and observations; it may not delete
> or move a `CampaignSite` that any `SiteEvaluation` refers to.** Unmatched
> machine-proposed sites with no verdicts may be retired; sites with verdicts
> are the humans' now.

**Two things reconciliation will want that are not there today, and are
deliberately not §14.0's job:**

- *Provenance on a site.* Once some sites are machine-proposed, "was this
  placed by a person or by run X" matters — the retirement rule above depends
  on it. A nullable pointer to the proposing job, added in v2.
- *Evidence on a verdict.* The moved-peak guard needs to know what the human
  was looking at: which map, at which coordinate. Today `SiteEvaluation`
  records neither. A nullable evidence reference (a file uuid and an xyz),
  added in v2 — and it is broader than PanDDA, since a verdict made in the
  campaign viewer today is made against dimple's map, not an event map.

Both are *additive nullable columns*. The recovery format is designed to key
on uuids and tolerate model-shape change ([PROJECT_RECOVERY.md](PROJECT_RECOVERY.md)),
so adding them later is a small migration and nothing more — **provided the
§14.0 serialiser passes fields through generically rather than enumerating
them.** That is the one property §14.0 must have for v2's sake:

> **§14.0 serialises campaign rows by field, not by hand-written column list**,
> so that a nullable column added in v2 rides through the snapshot without
> touching the serialiser.

**One hazard to record now because it is easy to design around and expensive
to discover:** PanDDA chooses its own reference dataset and aligns everything
to *that* frame, which is not our parent project's frame unless we make it so.
Site matching by centroid therefore needs either the per-dataset transform
into the parent frame (`lib/superposition.py` exists; the scene format's
`matrix` block records exactly this with provenance, §11) or control over
PanDDA's reference choice. v2's matcher must do one or the other; v1 records
centroids in the dataset frame and says so.

**Net for §14.0:** uuid on `CampaignSite`, the ownership split of §10.4, a
generic field-level serialiser, and the two invariants written into the
docstrings. Nothing about reconciliation changes its shape; everything about
reconciliation depends on it landing first.

---

## 10. Everything gleanable is also backed by a persistent artefact

The standing doctrine, from [PROJECT_RECOVERY.md](PROJECT_RECOVERY.md):
*anything user-authored that exists only in the database* must have an on-disk
artefact, because the database is reconstructible from the project directories
and the user's judgement is not.

### 10.1 The event↔file mapping needs no manifest, because it is structural

The composed-CData design (§7.1) satisfies the doctrine for free: the container
is serialised to `params.xml` in the job directory, and the recovery path is
built on exactly that
([restore_project.py:16](../server/ccp4i2/db/restore_project.py#L16)).

A separate manifest file would be a *second artefact asserting the same
association*, and recovery would have to decide which wins when they disagree.
The association being structural — the `CDataFile` sits inside the event record
— is what stops that question arising.

### 10.2 `CampaignEvent` is a projection, not a source of truth — **v2, decided**

Campaign-wide questions ("every dataset with an unbuilt event above 0.5 at site
3") cannot be answered from `params.xml`, and KPIs are per-job, not per-event.
A small `CampaignEvent` table keyed `(receipt job, site, event_idx)` answers
them.

Under the doctrine it is explicitly a **cache**: every field is derivable from
finished receipts' `params.xml`. So it ships with a rebuild command, which
doubles as the post-restore re-glean path.

**Decided: v2, and essential there.** It was going to be deferred on the
grounds that nothing queried it. §5.6 names the thing that does: the
cross-dataset triage queue that lets the campaign views replace a separate
review application. Still not v1 — v1 has no schema footprint — but it is the
first thing v2 builds, not the last.

### 10.3 A typed event-set `CDataFile` — **OPEN, default no**

Apply the §3.3 test: name the task that would take an event set as a typed
input. If a "refine this event" or campaign-triage *task* is coming, a subclass
earns its place, because that is what subType matching is for. If not, it is
decoration.

### 10.4 The gap this design must not build on top of

`ProjectGroup`, `CampaignSite` and `SiteEvaluation` appear **nowhere** in
`project_snapshot.py`, `export_project.py`, `restore_project.py` or
`import_i2xml.py`, and the receivers in `signals.py` watch only `Job`, `File`,
`Project` and `ProjectTag`.

So today: lose `db.sqlite3`, and **every hit/empty/unclear verdict in every
campaign is gone**, with all the project directories intact. `SiteEvaluation`
is the definition of what the rule protects.

It was missed for a structural reason: the snapshot is per project, and a
campaign spans projects. There is now an unambiguous answer, because
`one_parent_per_group` guarantees a single home.

**Recommendation, split along ownership:**

- Site definitions and group membership → the **parent** project's snapshot.
- `SiteEvaluation` rows → **each member's own** snapshot, because a verdict is
  about that dataset. Exporting or moving one member then carries its verdicts,
  and re-importing restores them.

**Prerequisite:** `CampaignSite` has no `uuid` — it is an integer pk with
`unique_together = [group, name]`. The recovery format keys on uuids by design
(primary keys are meaningless after a rebuild) and name is the mutable thing
whose instability motivated migration 0024. **A site needs a uuid before it can
be snapshotted or referenced from an event record.**

This is **§14.0**: its own PR, ahead of this work. It is a live defect in
existing campaigns, not a work package of this feature, and it should land
whether or not PanDDA is ever built.

---

## 11. A scene recipe per task

The receipt's outputs are only useful if opening one shows the right thing:
this event map, at the right contour, centred here, with this dataset in the
reference frame.

Two findings mean the format is ready:

- The superposition block already carries `method: "matrix"` with a row-major
  `mat`, a `vec`, **and provenance recording what the matrix was derived from**
  ([core.ts:186-207](../client/renderer/lib/scene/core.ts#L186-L207)). PanDDA's
  per-dataset alignment is expressible today. No schema change.
- `contourLevel` is **rmsd-relative** while `Optimal Contour` is **absolute**
  (§7.3). The division belongs in the recipe — once, next to the task that
  knows where the number came from — not in each consumer.

**Landed differently, 2026-09-24: scenes are outputs, not callables.** The
receipt writes one scene per event (`EVENTS[i].SCENE`) and an overview
(`SCENE`) as typed `CMoorhenSceneDataFile` outputs, gleaned like any file,
referencing its own outputs by job number and parameter; the Moorhen job
view applies a job's `SCENE` output when it has one and infers only for
jobs that wrote none (Martin's rule: honour, don't infer). The contour unit
question resolved itself on reading the resolver: the scene's `contourLevel`
is passed straight into Moorhen's store, which is **absolute**; the
grammar's "rmsd-relative" annotation was wrong and is corrected. So no
conversion belongs anywhere: the Z-map reads at 3.0 and an event map at the
`DISPLAY_CONTOUR` the receipt records (1.5 times the map's non-zero spread,
capped by the optimal contour, which on a poorly characterised run can sit
above the map's peak). The overview opens on the highest-scoring event with
the other event maps hidden; the v2 triage view walks events across
datasets.

**Original shape considered: a callable registered in the `Task` dataclass
alongside `reportPath`.** A recipe is a producer, not a document: contour, centre and
transform are functions of *this* job's event. It is the deterministic sibling
of `buildAuthoringPrompt`, and should share its vocabulary. This generalises
past PanDDA — any task with a defensible "here is what you should be looking
at" gets one. **OPEN:** the exact registration shape.

---

## 12. Windows

PanDDA does not run there, and the orchestrator task is not `ccp4_free` and
will not be offered where it cannot run.

It must not *add* barriers to it ever running there, which costs nothing if
decided now: `link_or_copy()` rather than symlinks (§3.2), `pathlib`
throughout, no shell-quoting assumptions in the invocation, ASCII-only
`print()` per the house rule. Fan-in and fan-out are pure database and file
work and should be Windows-clean regardless — a user may well curate a campaign
on Windows and run the analysis elsewhere, which §5.4 makes a supported shape
rather than an accident.

---

## 13. Decisions

**Decided:**

1. ~~Two tasks, three steps; fan-in and fan-out are not tasks.~~ **Amended
   2026-09-24: three tasks.** Fan-out is the `pandda_fanout` task, taking the
   manifest as a typed input (`fromPreviousJob` from the orchestrator), with
   the management command kept as a second way in for scripting; both call
   `lib/utils/jobs/pandda_fanout.py`. Fan-in is two plugin methods over the
   generic `object_method` endpoint (§14.2 status). What the original
   decision guarded against was Appendix A's shape — fan-out *inside the
   orchestrator's run* — and a separate task keeps every §8.2 property while
   adding the two the command lacked: a record of each fan-out (which tree,
   which manifest, what it created, when) and a place in the job list.
   Conditions, from the challenge that preceded the change: the task is
   `ccp4_free` and runs the receipts one at a time to completion inside
   itself, so its report is a true record and a laptop is not asked to start
   two hundred processes at once; the report is explicitly a *creation*
   record — receipts are top-level jobs of their own projects and carry their
   own state, and the fan-out job is not their parent; the manifest stays a
   plain `CDataFile` until whether a `FileType` row needs a data migration is
   settled; and on Azure "created" can only mean queued.
2. The orchestrator takes a declared `DATASETS` list and never touches the
   database; campaign-awareness lives in job construction.
3. ~~Local execution is the implemented path. The big-compute path is kept open
   by two zero-cost requirements (§5.4), not by building a runner abstraction.~~
   **Amended 2026-09-26 (decision 18):** the two §5.4 requirements stand and
   are still what makes a tree produced elsewhere usable; the "no runner
   abstraction" half is reversed. It was right with one implementation. There
   are now two real mechanisms (local subprocess, Azure Service Bus) plus a
   null one (`stage_only`), and the deployment that runs PanDDA already
   chooses between two Batch pools, so the abstraction has work to do on day
   one.
4. The Materia invocation contract is adopted verbatim — argv including both
   defensive arguments, `RAY_TMPDIR`, output-tree shape, progress symbol,
   failure catalogue.
5. Staging emits the contract's input tree, `Projects.csv` included, with clean
   `xtal-NNNN` names.
6. The receipt is typed, not a bucket; event files are members of a composed
   type inside a `CList`.
7. Fan-out is separate, idempotent, previewable, reported, retryable, and takes
   a tree rather than a runner — so it works on a partial tree and on a tree
   produced elsewhere.
8. No manifest file for the event↔file mapping; it is structural in
   `params.xml`.
9. Staging and delivery by hardlink-or-copy; never symlinks.
10. PanDDA sites merge into `CampaignSite` by centroid, with a moved-peak
    guard; `SiteEvaluation` survives reruns untouched.
11. Receipt files are copied into the member project, for self-containment.
12. **CCP4i2 depends on no PanDDA fork and pins no PanDDA version.** It targets
    the build CCP4 ships; staging compensates for what that build mis-reads
    (§4.7); the one gap that cannot be compensated is reported rather than
    worked around (§6.5); and a bundle bump to current upstream (§4.8) improves
    matters without being a precondition.
13. Memory-tuning controls are **not** exposed: the CLI flags do nothing (§6.2)
    and the env-var path is unvalidated and un-upstreamed (§6.5). The offered
    lever is dataset count.
14. The PanDDA executable is **discovered, not hard-coded** (§4.6), the resolved
    path is recorded, and capability is established by probing behaviour — there
    is no version number to test.
15. Where a build cannot read what we produce, **staging normalises it**
    (§4.7): `prepare_dict_for_pandda()` joins `prepare_mtz_for_pandda()`.
16. `MakeProjectsAndDoLigandPipeline` is deleted (Appendix A).

17. **External-run mode ships in v1** (§5.5), over the existing
    `Task.interactive` lifecycle.

18. **Run-target dispatch (v2), settled with Materia 2026-09-26** in the shared
    brief "Run-target dispatch: a CCP4i2 v2 briefing". The governing rule, for
    this and everything after it: **CCP4i2 builds generic capabilities with
    hooks; each environment registers what it has.** Nothing in `ccp4/ccp4i2`
    branches on a platform, and the desktop developer setup and the built
    desktop apps must always work with no configuration at all.
    - **One registry, two interfaces.** `CCP4I2_RUN_TARGETS` in settings maps a
      name to a dotted class path; CCP4i2 ships only `local`. A target
      implements axis A, `run_job(job)` (the whole CCP4i2 job moves — what
      `context_run` does today), or axis B, `submit/poll/cancel/logs` (one
      heavy program moves, the job stays in `RUNNING_REMOTELY`), or both, and
      declares which. `context_run` resolves axis A through the registry
      (`CCP4I2_JOB_TARGET`, default `local`); the orchestrator resolves axis B.
    - **The Azure code leaves CCP4i2.** `run_job_azure` and the Service Bus
      import move to Materia's `azure_extensions` as its registered axis-A
      target, in the same step as the registry. No deprecated built-in for a
      transition release: a clean seam for backend attachment. `EXECUTION_MODE`
      and `SERVICE_BUS_*` stop being read by CCP4i2.
    - **`DISPATCH_TARGET` is a separate parameter**, not a generalised
      `RUN_MODE`: the def.xml's enumerators are a static list, so a target
      registered only in a deployment's settings could never be one of them.
      `RUN_MODE` becomes `local | stage_only | dispatch`; `DISPATCH_TARGET` is
      validated in `validity()` against the registry and defaults to the
      deployment's single non-local target when there is exactly one.
      `stage_only` stays an instruction to stop after staging, not a target.
    - **The dispatch record is typed `outputData`** beside `PROVENANCE_*` and
      `CONTRACT_VERSION`: target name plus a target-tagged handle. Not
      `process_id` (an integer that means a local pid and feeds
      `postProcessCheck()`), and not the Qt-era `ServerJob` table (ten fields,
      referenced only by `import_sqlite.py`), which is retired when the typed
      record lands.
    - **Targets never classify failure.** `logs(handle)` returns a path under
      the job directory, present by the time `poll` reports terminal; the §4.5
      catalogue classifies from it, unchanged for a local subprocess, a Batch
      task or an HPC job alike.
    - **Reconcile is idempotent and user-driven**: a plugin method over the
      generic `object_method` endpoint plus a management command wrapping the
      same function; `pandda_fanout` calls it first; it supersedes
      `cleanup_stale_jobs` for dispatched jobs. Not `JobInteractiveSession`
      (that models a user window, and `open_session` refuses a job in
      `RUNNING_REMOTELY`); what is shared is `finish_session`'s out-of-process
      completion, lifted into one "complete a job from outside its process"
      function both paths call.
    - **`cancel` is in the interface, not the UI**, mapping to `INTERRUPTED`.
    - **The UI offers dispatch only when the API reports a registered target**,
      the same pattern as the scenes Generate tier.
    - Separately, and on its own PR: plugins register their own CData classes,
      so `CPanddaEvent` and `CPanddaDataset` leave `core/` (the hard-coded list
      in `def_xml_handler._build_class_registry()` is the same closed-enum
      defect in a different place).
18. A short receipt returns **`UNSATISFACTORY`, never `FAILED`** (§7.2), so
    what arrived is still published.
19. The composed event type is a **registered core class** (§7.1), and a spike
    proves the shape before anything is built on it (§14.1).
20. The campaign-persistence defect (§10.4) is **its own PR, ahead of this
    work** (§14.0) — it is a live defect, not a work package of this feature.
21. **The campaign Moorhen views are the review application** (§5.6).
    `CampaignEvent` is v2-essential because the triage queue needs it;
    Reinspect's review role is absorbed once v2 lands, and only Batch
    submission remains a separate question, for Materia.

**Open:**

1. Registration shape for scene recipes in `TASKS` (§11).
2. Typed event-set `CDataFile` — default no, revisit when a consumer task is
   named (§10.3).
3. Campaign snapshot split — recommendation in §10.4 needs sign-off before
   §14.0 is written.
6. **What we tell a user whose job will not fit.** §4.7 removes the dictionary
   problem and §6.4 gives the estimate, so the remaining question is what the
   warning *says*. "Reduce the dataset count" is supportable; "install this
   branch" is not, while §6.5 stands. Worth noting for whoever writes the text
   that CCP4's `pandda2` env ships `python3.9` **and `pip`**, so a determined
   user can `micromamba run -r $CCP4/share/mamba -n pandda2 pip install …` into
   it and reuse CCP4's interpreter, dependencies and bundled checkpoints with no
   conda or container — but that is a thing a user may choose, not a thing this
   task should recommend.

---

## 14. Delivery

Three things ship in order: a defect fix that is not part of this feature, a
spike that retires the design's riskiest assumption, then v1 and v2.

### 14.0 First, and independently: the campaign-persistence defect

§10.4 is not a PanDDA work package. `ProjectGroup`, `CampaignSite` and
`SiteEvaluation` appear in **none** of `project_snapshot.py`,
`export_project.py`, `restore_project.py`, `import_i2xml.py` or `signals.py`,
and `CampaignSite` has no `uuid`. That is a **live data-loss defect in existing
campaigns today**, entirely independent of anything here: lose `db.sqlite3` and
every verdict is gone with the project directories intact.

Filing it as "v2 work package 0" buried it behind a feature it does not depend
on. It is **its own PR, now, ahead of this sequence** — and it happens to
de-risk v2 as a side effect.

Scope: `CampaignSite.uuid` + migration, campaign state into the snapshot per
§10.4's ownership split, and a destroy-and-restore round-trip test that does
not exist today for campaigns.

### 14.1 Then a spike, before the delivery plan can be trusted

**Build a throwaway task with a `CList` of a composed type containing
`CDataFile`s, run it through glean, and look at the `File` rows and the UI.**

Nothing in the codebase exercises that shape end to end (§7.1). `CDmDomain`
establishes the composed-subItem pattern but holds only strings; the three
mechanisms that must cooperate — registry resolution, `find_all_files()`
recursion, `objectPath()`-derived `param_name` — each check out in isolation
and have never been checked together. And the failure mode is a silent
downgrade to `CString`, not an error.

It is a day's work and it is the single highest-variance assumption in this
document. Everything in §14.2 items 3–5 depends on it; if it does not hold,
the receipt's shape changes and so does the delivery plan.

**Outcome (2026-09-24): it holds.** A throwaway task with `outputData.EVENTS`
a `CList` of a core-registered composed type holding a `CMapDataFile` and a
`CPdbDataFile`, run under `i2run`, gave one `File` row per nested file with
`job_param_name` exactly `EVENTS[0].EVENT_MAP`, `EVENTS[0].POSE`,
`EVENTS[1].EVENT_MAP`, `EVENTS[1].POSE` and the right mime types; `params.xml`
carried each nested file with its `dbFileId`; `get_job_plugin` (what the
`/container/` endpoint uses) read it back as typed items with existing files;
and the JSON encoder emitted the items as `CSpikeEventPair` with typed nested
files. The parser gave the list a real subItem class, asserted directly. Three
things learned on the way, none blocking:

- **A misregistered *subItem* is quieter than line 403 suggests.** That
  warning fires for a top-level `<content>`; for a `<subItem>` the parser
  simply sets no `subItem` qualifier and `makeItem()` returns `CString` with
  no message at all. The receipt's test must assert the item class.
- **`checkOutputData()` does not name nested files.** It defaults paths only
  for top-level `CDataFile`s and for lists *of* `CDataFile`s; a file inside a
  composed item gets no default. The receipt sets every nested path itself
  (it would anyway, since the names come from the PanDDA tree).
- `lib/utils/containers/get_container.get_job_container` (bare `CContainer`
  + `loadContentsFromXml` of a def.xml) returns an empty container and is not
  what the UI uses; the read path to test against is `get_job_plugin`.

The client already dispatches composed list items by `_class` to a
per-class element (`CDmDomain` has `cdmdomain.tsx`), so `CPanddaEvent` will
want one such element for the receipt's page; `GenericInterface` will not
fall over without it, but will not show the files well either.

### 14.2 v1 — the runnable path

v1 has **no database-schema footprint**: no model, no migration, no signal
receiver. That is bought by deferring site matching — v1's events carry
centroids and scores but **no `CampaignSite` reference**, so nothing touches
the campaign tables. Adding the reference later is additive: existing receipts
keep their events and the v2 matcher reads centroids that are already there.

| # | Work | Where it lives |
|---|---|---|
| **1** | `prepare_dict_for_pandda()` (§4.7) | orchestrator plugin dir |
| **2** | Staging: `link_or_copy()`, contract-conformant tree, `Projects.csv`, uuid-keyed manifest (§3.2–3.3) | orchestrator plugin dir; **imports** `prepare_mtz_for_pandda` and `collect_pandda_datasets` from `lib/pandda_export.py` rather than moving them |
| **3** | `CPanddaEvent` — the composed event type (§7.1) | **`ccp4i2/core/CPanddaEvent.py`** + one line in `implementation_modules` |
| **4** | Orchestrator task `pandda_campaign` — def.xml, executable discovery (§4.6), invocation (§4.1–4.2), progress parsing, failure classification, `runTimeValidity()` sizing warning (§6.4), provenance (§4.4) | orchestrator plugin dir |
| **5** | Receipt task `pandda_events` — def.xml, tree reader, verify-and-report-short (§7.2) | receipt plugin dir |
| **6** | `pandda_fanout` management command (§8) — idempotent, previewable, reported, retryable; takes a tree and a manifest | `db/management/commands/` |
| **7** | External-run mode (§5.5) over `Task.interactive` + `lib/utils/jobs/interactive.py` | orchestrator plugin dir + `TASKS` flag |

**Footprint outside the two plugin directories** — larger than an earlier draft
of this document claimed, and the correction is item 3:

| File | Change |
|---|---|
| **`ccp4i2/core/CPanddaEvent.py`** | **new core module** — the composed type cannot live plugin-side (§7.1) |
| **`core/task_manager/def_xml_handler.py`** | **one line** in `implementation_modules` |
| `core/tasks.py` | two `TASKS` entries, ~20 lines |
| `db/management/commands/pandda_fanout.py` | one new file |
| `client/.../task-chooser.tsx` | one category entry per task, 2 lines |
| `tests/unit/…`, `tests/i2run/…` | new files only |

Still no client interface (`GenericInterface` renders both def.xmls), still no
migration, and `lib/pandda_export.py` and the `export_pandda` endpoint are
still **not modified** — which defers the coordination question with Materia,
whose tooling scripts that endpoint.

**Status, 2026-09-24.** Items 1–6 are on `pandda-campaign` (commits
964607ac4, 083110d07/0c0db80fb, c7f3b9444/0d39d9c87, 243e05620, plus the
client elements). §14.0 is PR #608 and §14.4 is PR #610, both against
`django`. Two things differ from the table above:

- **Item 7 is done as a run mode, not an interactive session.**
  `RUN_MODE=stage_only` stages the tree, records provenance and finishes
  `SUCCEEDED`; the tree is shipped; fan-out later takes the resulting
  `pandda2_out/` with that job's manifest, whose provenance carries the job's
  uuid, so the receipts trace to it. That satisfies every requirement §5.5
  states — the job exists with its `DATASETS`, manifest, provenance and
  receipts — without a session. The `Task.interactive` route was not taken
  because *Run opens a session instead of dispatching*, so the staging the
  user needs before running elsewhere would not happen until the session
  finished; and a plugin that dispatches and then waits for the tree would
  hold a worker for days. If a "hand the tree back to this job" affordance
  is wanted in the UI, it is a fan-out button on the stage-only job, not a
  session.
- **Fan-out is a task** (`pandda_fanout`, decision 1 as amended) as well as
  a command; receipts run one at a time to completion inside it
  (`RUN_RECEIPTS`), a preview is `DRY_RUN`, and a tree produced elsewhere is
  `PANDDA_OUT_DIR`. On Azure that queues them.
- **Fan-in is v1 after all, as plugin methods.** Filling `DATASETS` from
  the campaign turned out to be the difference between a screenshot and a
  runnable job, so it landed with v1: `campaignCandidates()` and
  `fillDatasetsFromCampaign()` on the plugin, reached through the generic
  `object_method` endpoint (no PanDDA-specific route; Martin's call,
  2026-09-24), with the database reading in `lib/utils/jobs/pandda_fanin.py`.
  `CPanddaDataset` carries `PROJECT_UUID` and `SOURCE_JOB_UUID` explicitly,
  because the files are imported into the orchestrator's project when the
  job runs and would otherwise stamp the parent's uuid into the manifest.
  v2 item 12 is therefore only the campaign-page affordance that creates
  the job and calls the fill.

- **The orchestrator's report recapitulates PanDDA's own analysis.** The
  bundled PanDDA 2 writes no HTML summary (`analyses/html_summaries/` stays
  empty), so `pandda_run_summary.py` recovers what Reinspect recovered from
  PanDDA 1: the events table joined with each dataset's `events.yaml` (score,
  build score, RSCC), the sites table, and each `processed_dataset.yaml`
  (processing resolution, comparator count, selected model, candidate-event
  counts before and after the size and score filters). It goes into
  program.xml as `<analysis>` and the report draws events per site, binned
  histograms of event fraction, hit probability, resolution and R-free, and
  events, sites and per-dataset tables. Two facts worth keeping:
  `processed_dataset.yaml` carries numpy scalars as `!!python/object/apply`
  tags, which `yaml.safe_load` refuses, so it is read with a loader that
  yields `None` for them; and the report viewer draws `barchart` but not
  `histogram`, so the binning is done in Python.

Also found on the way: the gleaner stored only `float` KPIs, so every `CInt`
KPI in the tree was silently dropped (fixed in 0c0db80fb); imported input
files nested inside a composed list item are registered under their
`objectName` (`XYZIN`), not their path (`DATASETS[0].XYZIN`), unlike outputs
— harmless, noted.

### 14.3 v2 — the campaign-aware path

| # | Work | Gated on |
|---|---|---|
| 8 | Site matching: PanDDA sites → `CampaignSite` by centroid with the moved-peak guard (§9); events gain a site reference | §14.0 |
| 9 | `CampaignEvent` projection + rebuild command (§10.2) | §14.0, item 8 |
| 10 | **Campaign triage view**: every event across the campaign, ranked, stepped through in the existing campaign Moorhen page; *accept this pose* as the generalisation of place-ligand (§5.6) | items 9, 11 |
| 11 | Scene recipes (§11) — the triage view's per-event scene is the first consumer | — |
| 12 | Fan-in as a campaign-aware endpoint; `pandda_fanout` gains an API wrapper (the command remains the implementation) | — |

**Run-target dispatch (decision 18), in this order, each its own PR.** Every
one is judged first against the desktop invariant: an unpacked Electron dev
tree and a packaged app see no behaviour change, nothing new is imported at
request time, and the CCP4-free unit tier keeps passing on stock Python.

| # | Work | Gated on |
|---|---|---|
| 13 | This note: record decision 18 and the §5.2 / decision 3 reversals | — |
| 14 | `ccp4i2.lib.dispatch`: registry, the two interfaces, `LocalTarget`; `context_run` resolves axis A through it; `run_job_azure` and the Service Bus import removed (Materia lands its registered target in the same step) | — |
| 15 | CData classes registered by their plugin; `CPanddaEvent` / `CPanddaDataset` move out of `core/` | — |
| 16 | Axis B: `DISPATCH_TARGET`, the typed dispatch record, the reconcile (plugin method + command), the shared out-of-process completion, the API reporting registered targets, the UI affordance | items 14, 15 |

### 14.4 Unrelated, whenever convenient

**Delete `MakeProjectsAndDoLigandPipeline`** (Appendix A) — small, its own PR.

### 14.5 Suggested PR shape

Given `django` merges serially with up-to-date-required:

0. **§14.0** — the persistence defect. Ahead of everything, on its own merits.
1. **v1 items 1–2** — dictionary normaliser and staging. Pure gemmi and file
   operations, no PanDDA needed to test any of it, and standalone value: it is
   `export_pandda` done properly, and is the piece Materia could adopt.
2. **v1 items 3–7** — the composed type, both tasks, fan-out, external-run.

The spike (§14.1) sits between 1 and 2 and lands nothing.

### 14.6 What must be right in v1, because it is expensive later

1. The declared `DATASETS` list rather than a campaign lookup (§3.1) —
   retrofitting means rewriting the task, and it is what makes `i2run` testing
   possible at all.
2. The typed receipt with event files inside a composed `CList` (§7.1) —
   retrofitting a bucket means migrating receipts that already exist. **This is
   the item the spike exists to de-risk.**
3. Provenance in params (§4.4) — cannot be added retroactively, and without it
   two receipts are not comparable.
4. Fan-out takes `(tree, manifest)` and is idempotent (§8.4) — this is what
   keeps the big-compute path open and partial trees recoverable.
5. `prepare_dict_for_pandda()` (§4.7) — the only *correctness* item on the
   list rather than an architecture one.
6. `--dataset_range 0-999999999` and `--ligand_pdb_regex ligand.pdb` (§4.1) —
   two literals, each silent and catastrophic if omitted.

### 14.7 Branch discipline: a long-lived branch beside an advancing `django`

Work happens on `pandda-campaign`, in the worktree
`~/Developer/ccp4i2-with-pandda`, cut from `django` at `ca0ec99ff` (#392).
`django` will advance underneath it for weeks. Four rules keep the rebases
trivial:

1. **New files in new places.** Both plugin directories, `core/CPanddaEvent.py`,
   the fan-out command, the tests — all new paths, which never conflict.
2. **The shared-file touches live in one tiny commit, kept at the tip.** Only
   three existing files change: `core/tasks.py` (two entries),
   `core/task_manager/def_xml_handler.py` (one line), `task-chooser.tsx` (two
   lines). A conflict there is re-resolved in seconds if it is isolated, and a
   nuisance if it is entangled with real work. No reformatting of anything.
3. **Rebase, never merge, and often.** `git rebase origin/django` after each
   merge to `django` that touches anything nearby (#598/#603 would have been
   two such). Linear history keeps `git bisect` and the eventual squash honest.
4. **No migration ever rides the long branch.** §14.0 carries a migration, so
   it goes on its own short branch off `django`, merges, and `pandda-campaign`
   rebases onto it. Same for the deletion (§14.4). Migration numbering on a
   long-lived branch beside an advancing main is the one conflict that is
   *not* trivial to re-resolve.

The spike (§14.1) is throwaway and lands nothing; it can live on a scratch
branch or in this worktree's uncommitted state.

---

## 15. Tests, and the data we already have

### 15.1 Real fixtures exist, on `/Volumes/LocalStore/pandda`

This matters more than usual, because the expensive part of testing this task
is getting real PanDDA inputs and outputs.

| What | Where | Gives us |
|---|---|---|
| A **contract-conformant staged tree**, 201 datasets: `Projects.csv` + `datasets/xtal-NNNN/{final.pdb, final.mtz, dict.cif, ligand.pdb}` | `BAZ2B/` | The exact shape §3.2 must produce. Staging can be tested by *reproducing* it from CCP4i2 projects and diffing |
| A **real, complete output tree** over all 201 datasets — `analyses/pandda_analyse_events.csv`, `processed_datasets/`, 62 event maps | `BAZ2B/pandda2_subset60/` | The receipt reader's fixture (§7.1), and the fan-out fixture (§8.4), with no run required |
| A **partial output tree** (no events table, no event maps) | `CDK4CyclinD1/pandda2_out/` | The incomplete-run case for fan-out and receipts (§8.4) |
| A **large-cell campaign** (58 × 64 × 186 Å) | `CDK4CyclinD1/` | The sizing case §6.1 quotes |
| **A/B/C path comparison** + `ABC_report.md` | `pathx_runs/` | Evidence on the §6.5 experiment, including its own honest account of what is confounded |

Two things the fixtures do **not** give us:

- **A `value_order` dictionary.** BAZ2B's `dict.cif` files are `type`-style, so
  they do not exercise §4.7's problem case at all. v1 item 1's test needs a
  PDBx-derived or newer-acedrg dictionary constructed for the purpose — which
  is cheap, but must be done deliberately or the test passes vacuously.
- **A small run that finishes quickly.** There is no such thing. PanDDA
  characterises the ground state from comparator datasets and refuses to run
  below `--min_characterisation_datasets`, whose default is **25**
  (`pandda_gemmi/constants`); a 3-dataset subset cannot run, and a 25-dataset
  one is not representative. The end-to-end run test for the orchestrator
  (v1 item 4) therefore **needs the volume**: it runs over BAZ2B on
  `/Volumes/LocalStore/pandda`, is skipped when that is not mounted, and is
  an afternoon rather than a unit test. (Corrected 2026-09-24; an earlier
  draft claimed a 3-dataset run was the smallest honest check.)

The three-dataset subset cut into `demo_data/pandda_baz2b_mini/` (v1 item 2,
5.2 MB, kept deliberately) is a **staging** fixture: real dimple outputs to
stage and diff against the real tree. It is not, and cannot be, a PanDDA
input. A complete run over all 201 BAZ2B datasets already exists on the
volume at `BAZ2B/pandda2_subset60/` (the name notwithstanding: its
`input.yaml` lists 201 datasets, and it has the events table and 62 event
maps), which makes it the fixture for the receipt over a real run and for
fan-out (§8.4) with no run required. `CDK4CyclinD1/pandda2_out/`, named
above, is an *incomplete* run — no events table, no event maps — and is
useful only as the partial-tree case.

### 15.2 The assertions

- **The spike (§14.1), before anything else:** a `CList` of a composed type
  holding `CDataFile`s yields one `File` row per nested file after glean, with
  `param_name` of the form `EVENTS[3].EVENT_MAP`, and the class resolves to
  itself rather than silently to `CString`. **Assert the resolved class**, not
  just the behaviour — the fallback at
  [def_xml_handler.py:403](../server/ccp4i2/core/task_manager/def_xml_handler.py#L403)
  makes a misregistration look like success.
- **v1 item 1:** a `value_order` dictionary comes out carrying a `type` column
  whose tokens the old reader's map accepts; a `type`-only dictionary is passed
  through unchanged; both round-trip through gemmi `ChemComp` with identical
  bond orders.
- **v1 item 2:** the staged tree matches the contract — `datasets/xtal-NNNN/`
  with clean integer names, `Projects.csv` present and parseable, a relabelled
  FreeR column where the source needed one. Diffable against `BAZ2B/`.
- **v1 item 4:** `i2run` against three hand-specified datasets, no campaign, no
  database beyond the job itself. **If this needs a campaign, §3.1 has been
  violated** — this is the architectural test, not a functional one. Plus: the
  emitted argv contains both defensive arguments, and the child environment
  contains no PanDDA switch we did not mean to set (§6.5's presence-check
  semantics make a stray `"0"` an *enable*, so this is asserted on the
  environment we build, not on our intent).
- **v1 item 5:** a receipt short of a declared event map returns
  `UNSATISFACTORY` and **still gleans what arrived** (§7.2) — the assertion is
  the `File` rows, not the status alone. And glean finds every nested `CDataFile` — assert `File` rows for
  event maps *inside* the `CList`, since that is the load-bearing assumption of
  §7.1. Read against `CDK4CyclinD1/pandda2_out/`, including its
  zero-event datasets.
- **v1 item 6:** fan-out twice over the same manifest creates nothing the second
  time; fan-out over a **partial** tree creates receipts for what is there and
  reports the rest as absent; fan-out against a tree with no local run at all
  works, which is §5.4's second requirement as an executable assertion.
- **§14.0:** destroy-and-restore round trip; verdicts survive. The test that
  proves the doctrine holds, and it does not exist today for campaigns.
- **v2 item 11:** a recipe's emitted contour for a known absolute `Optimal
  Contour` and map rmsd equals the expected σ value. One assertion, pinning a
  unit conversion that has already cost debugging time once.

---

## Appendix A: `MakeProjectsAndDoLigandPipeline`, and why it is deleted

The previous attempt at a job that creates jobs in other projects. It is
registered in `TASKS` and has a bespoke client interface, but **it cannot
run**: the first statement of `startProcess` is

```python
from ccp4i2.core.CCP4Modules import JOBCONTROLLER, PROJECTSMANAGER
```

and `CCP4Modules` now exports only `PREFERENCES()`. Beneath that: `CException`
is never imported, so all three `except CException` clauses are a `NameError`;
the work-dir handler prints an undefined `workDirectory`; `Rfactor` is filled
from `rfreeNodeText`; and `updateProject(projectId, 'parentprojectid', …)`
targets a Qt-era column that `ProjectGroup` replaced.
[error-handling-remediation.md:955](error-handling-remediation.md#L955) already
records it as the highest defect density in the tree. It is not in the task
chooser, so it is unreachable from the UI.

**Its four failure modes are this design's specification:**

**A.1 It fanned out inside `startProcess`.** Every error path is a bare
`continue`: a failed dataset is skipped and the only trace is a `Warnings`
string in the running XML. Nothing records which datasets succeeded, so nothing
can retry the rest. → §8.2.

**A.2 It hand-rolled job creation** — `createJob`, `getJobInfo`, `mkdir`,
plugin instantiation, `setDbData`, `saveParams`, copy PARAMS to JOB_INPUT, then
`JOBCONTROLLER().runTask`. Every one a framework internal; the framework moved
and the copy did not. → §8.1.

**A.3 Concurrency was in-process `waitForFinished` polling with a cap of 10**,
pinning the parent job alive for the whole fan-out. → §8, fan-out has its own
lifecycle.

**A.4 Its summary scraped the children's XML by XPath** for resolution, R-free
and R-factor. One of the three was copy-pasted wrong and nobody noticed,
because ten bare excepts turned every failure into `"N/D"`. → §10.2: numbers
the campaign view needs are gleaned as KPIs and read from the database, not
re-extracted from subjob XML by the parent.

**Deletion touches:** the pipeline directory; the `TASKS` entry
([tasks.py:110](../server/ccp4i2/core/tasks.py#L110));
[citations.py:34](../server/ccp4i2/core/citations.py#L34); the report-port test
at `tests/unit/lib/test_report_lxml_ports.py:77`; the 89-line client interface
and its registry line in `task-container.tsx:233`; and the two icons. Query for
jobs with that `task_name` first — the house rule keeps dead tasks registered
so old jobs still open — though since it raises `ImportError` on entry, any
such job failed.

---

## Appendix B: external references

Not vendored here, because they are maintained elsewhere and copies drift:

| Document | Repo | What it settles |
|---|---|---|
| `CCP4I2_PANDDA_INVOCATION_CONTRACT.md` | `materia` | argv, env, input/output tree shape, progress symbol, failure catalogue, image pinning |
| `PANDDA2_ON_AZURE.md` | `materia` | Memory measurements (§6), dead tuning flags, single-node constraint, SKU sizing, the 2026-06-06 ownership decisions |
| `pandda-batch.bicep` | `materia` | The Batch pool: `Standard_E16ds_v4` (16 vCPU / 128 GiB / ~600 GB local NVMe for `RAY_TMPDIR`), scale-to-zero |
| `docs/RUN_LIFECYCLE.md`, `docs/MULTI_RUN_DATA_MODEL.md` | `pandda-inspect-api` | Run lifecycle, retry semantics, merge-vs-replace ingest, cross-run event identity |
| `ConorFWild/pandda_2_gemmi` `master` | upstream | The version to target. PRs #96–#99 merged; HEAD 2026-06-17. Compare against `$CCP4/share/mamba/envs/pandda2` to see what a bundle is missing |
| `CROWTHER_DEPLOYMENT_NOTE.md` | `~/Developer/` (not in a repo) | The two env-var switches of §6.5 and their presence-checked semantics. **Describes experimental, un-upstreamed, not-end-to-end-validated work**; cited here so the switches are recognisable, not as something to build on |
| `REMOTE_JOB_EXECUTION_PLAN.md` | this repo | The ssh/qsub shepherd, if §5.5 grows a third local-ish runner |
