# Accepting PanDDA poses across sites, without touching the rest of CCP4i2

**Status:** planning document, for refinement. Nothing here is built.
**Companion to:** [`pandda-campaign-design.md`](pandda-campaign-design.md), whose
§9.2 describes the run-site axis this builds on.

## 1. The problem

Site-ordered work is the point of the run-site view: you look at one site across
every dataset that has an event there, judge the poses together, and move on to
the next site. That ordering is what makes the view worth having, and it is also
what breaks the model-building path.

One dataset commonly has events at several sites. In campaign 68 (`CDK4CyclinD1Fragments`,
PanDDA job 16) there are 94 events over 36 datasets: site 1 has 15 events in 12
datasets, site 2 has 22 in 14, site 4 has 22 in 18. Datasets recur across sites by
construction.

So the sequence is: accept a pose for dataset X at site 1, work through twenty
other datasets, then arrive at site 2 where X appears again. The second
acceptance must land on **the model that already has the first ligand in it**,
not on the dataset's original DIMPLE output.

## 2. What breaks today

Three failure modes, which a fix should be judged against separately:

1. **Staleness.** The site-2 session resolved X's coordinates before the site-1
   edit existed, and holds a molecule read from the older file.
2. **Latency.** The model of record is
   `_latest_finished_job(project, REFINE_TASK_NAMES)` — and that tuple is
   `servalcat_pipe`, `prosmart_refmac`, `refmac`, `i2Refmac`. **`coordinate_selector`
   is not in it**, so today's "Push to CCP4i2" does not advance the model of
   record at all; only a finished refinement does. Site-ordered work crosses
   datasets in seconds and refinements take minutes, so this is not an edge case,
   it is the normal case.
3. **Forking.** Two acceptances derived from the same parent produce two siblings
   with one ligand each, and nothing merges them.

Chaining file to file fixes (1) and can detect (3). Nothing makes (2) go away,
because it is inherent in "the model of record is the output of a job".

## 3. Proposal: acceptance is the durable fact, the model is derived

Do not add a ligand to the file just augmented. Record the **acceptance**, and
make the built model a function of the accepted set:

```
built(dataset) = base + (accepted − contained(base))
```

* `base` — the dataset's current model of record (DIMPLE early on, the most
  recently refined ligand-bearing model later).
* `accepted` — poses a person has accepted, each stored with its coordinates.
* `contained(base)` — those already present in `base`.

**Set difference, not append.** Append is correct only for the first refinement;
after that the base already holds the earlier ligands and appending doubles them.

Properties that matter for this workflow: order-independent, idempotent, and with
no parent to race against. Visit sites in any order, revisit, undo — the model is
always the same function of the same facts, and it is correct *immediately*,
which is what removes failure mode 2.

This follows the grain of `SiteEvaluation`, which already makes the human verdict
the durable thing that survives a rerun. An accepted pose is that verdict's
physical counterpart.

## 4. Containment, from provenance that already exists

`FileUse(file, job, role, job_param_name)` already records the coordinate file
each job consumed. So containment needs no geometric matching:

> A pose is in the base **iff** a file it was composed into is, or is an ancestor
> of, the base — walked through `FileUse`.

Containment is therefore inherited along a refinement chain for free: refine,
refine again, and the lineage still passes through the file that first held the
ligand.

Geometric matching is explicitly rejected as the primary mechanism. Refinement
*moves* the ligand, and two poses at adjacent subsites can fall inside any
tolerance worth picking.

**Stamp at composition, not at refinement.** When the composed model is written we
know exactly which poses went into it and what file resulted, so the record is
written there. Nothing needs to hook refinement completion — which is what keeps
generic tasks out of this entirely.

**Degradation.** A refinement launched from the job list, outside any campaign
flow, still writes `FileUse`, so lineage still reaches a known artefact and
containment is still exact. It fails only for a coordinate file imported from
outside with no lineage at all. There the honest behaviour is to **refuse to
compose and say the base's contents are unknown** — a silent double-add produces a
model that looks plausible and refines to nonsense.

## 5. Footprint

This is the section for reviewers who do not consider PanDDA core. The intent is
that this work is **additive and campaign-scoped**, and that a reviewer can
satisfy themselves of that quickly.

**Added**

| | |
|---|---|
| `AcceptedPose` | new table: dataset, site, event, coordinates, provenance. Campaign-scoped, beside `CampaignSite` / `SiteEvaluation`. |
| pose ↔ composed-file rows | new table recording which poses went into which composed coordinate file. |
| composition service | new module under `lib/`, plus endpoints under the existing campaign viewset. |

**Not touched**

* **No column is added to `File`, `Job` or `FileUse`.** Every new foreign key lives
  on a new table; existing models gain at most a reverse accessor, which costs
  nothing and changes no query.
* **No migration alters an existing table.** New tables only — which also means no
  risk to DDU's production data, the failure mode that bit migration 0024.
* **No generic task is modified.** `servalcat_pipe`, `refmac`, `prosmart_refmac`
  and friends are unchanged and unaware. Composition happens *before* refinement
  and produces an ordinary coordinate file; the refinement task cannot tell it
  from a file a user picked by hand.
* **No change to `REFINE_TASK_NAMES`** or to how the model of record is chosen.

The composed model enters the normal CCP4i2 data flow as a job output, so
provenance, the job list and the file browser all work without special cases.

## 6. Decisions needed

1. **Carrier for the composed model.** Reuse the existing `coordinate_selector`
   task (XYZIN → XYZOUT; already what "Push to CCP4i2" creates), or add a
   campaign-scoped `campaign_compose` wrapper? Reuse adds no task to the task
   list but leaves the job's purpose legible only from its title; a new wrapper is
   self-describing and still purely additive, at the cost of one more entry in a
   task list colleagues already find long. **Leaning:** new wrapper, for honest
   provenance — but this is the decision most worth challenging.

2. **Where `AcceptedPose` sits relative to `SiteEvaluation`.** A `hit` verdict with
   no placed ligand is meaningful — judged a hit, not yet modelled — so they are
   not the same row. Separate table referencing the evaluation, or a nullable
   relation? **Leaning:** separate, referencing.

3. **Ligand dictionaries.** Refinement needs restraints for *every* ligand in the
   model. Composing N poses from N datasets means merging N dictionaries, or
   composing one. This is campaign-side work and touches no generic task, but it
   is real and is the most likely source of refinement failures. Needs its own
   section before implementation.

4. **Chain and residue numbering.** N fragments each called `LIG` or `DRG` in their
   own numbering need a deterministic assignment, or refinement sees duplicates.

5. **Clash policy.** Two accepted poses at adjacent subsites can overlap, and
   alternate conformations of one site are a legitimate case wanting altlocs.
   Refuse, flag, or model as altloc — but never silently interleave.

## 7. Traps

**Frame.** In the site view every pose is superposed into the exemplar's frame for
display. A pose nudged there and accepted is being read in *transformed*
coordinates, and must be stored in the dataset's own frame — the inverse fit
applied before writing. This is the same class of error as the site-origin
negation, and it will look entirely correct in the view that produced it.

**Historical coordinates.** Once a pose is incorporated and refined, the refined
position supersedes the accepted one, which must never be re-applied over it.
Keying containment on identity rather than geometry gives this for free.

**Asymmetric undo.** Before incorporation, un-accepting is dropping a row and
recomposing — free. After, the ligand is in the refined model and removing it is a
deletion against the base producing a new artefact. Acceptance is *provisional*
until refinement and *committed* after, and the interface should show which poses
are still provisional so a person can see what a refinement is about to make
permanent.

## 8. Not in scope

Adoption of a run site as a `CampaignSite`; the side panel; any change to how
refinement is run or to what counts as the model of record.
