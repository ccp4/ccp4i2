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

## 3. What was wrong with the first proposal

The first version of this document proposed storing accepted poses as fragments
and deriving the model as `base + (accepted - contained(base))`, composed server
side. Martin's objection killed it, and it is recorded here because the reasoning
generalises.

Composing server side means "Add ligand here" stops being a Coot operation and
becomes a round trip, after which a composed PDB has to **replace the coordinate
set in Coot**. That destroys the loop the work actually consists of: place the
ligand, jiggle it, adapt the neighbouring side chains, look again. Nobody models
a ligand in one shot, and a design that is only correct if they do is not a
design.

Worse, it makes the database's claim falsifiable by the next mouse click. If the
row says "pose accepted at site 2" and the person then deletes the ligand
coot-matically, the database asserts a provenance that the model does not
support, and nothing detects it.

The mistake was making the **modelling act** trigger the **commitment**. They are
different things happening at different times, and only the second belongs to the
server.

## 4. Proposal: edit locally, commit explicitly, derive the record

**The modelling act stays entirely in Coot.** "Add ligand here" does what it does
today: places the pose in the local coordinate set, instantly, and the person
jiggles it and adapts the protein around it. No round trip, no replacement of the
coordinate set, no change to that loop at all.

**The commitment is the push**, which already exists. At that moment the client
holds the authoritative molecule -- every adjustment included -- and sends it
whole. One round trip, at a moment the person chose.

**The record is derived from what was pushed, not asserted alongside it.** On
commit, the server reconciles the campaign's per-site record against the molecule
it just received: a site with a ligand in it is modelled, a site without one is
not. Delete the ligand in Coot and push, and the site-2 record withdraws itself,
because it was never an independent claim.

That is the inversion the objection forces, and it is worth stating plainly:

> Do not record the acceptance and derive the model. **Record the model and
> derive the acceptance.**

Geometry is the right tool *here*, where it was the wrong tool for deciding
whether a fragment was already incorporated, because the question genuinely is
spatial: is there a ligand at this site in this model.

**Origin is still worth keeping**, and geometry cannot supply it: "this ligand
came from event 3 of site 2 of run 16" is real provenance that a spatial check
cannot reconstruct. So the origin is recorded when a person takes a pose into
their model, and treated as a *hint about where the ligand came from*, reconciled
against the pushed molecule at every commit. It never defines the model.

## 5. The model of record, and why the latency goes away

Failure mode 2 (§2) does not need composition to fix. It needs the push to count.

For campaign purposes the head of a dataset is **the most recent campaign
coordinate artefact**: the latest push if there is one after the latest
refinement, else the latest refinement. A push is synchronous, so the moment you
finish at site 1, site 2 will load a model containing that ligand. No job has to
finish first.

This is a campaign-side resolver -- a different question from "what has been
refined", asked by campaign views only. `REFINE_TASK_NAMES` and every generic
consumer of it are untouched.

Staleness and forking (failure modes 1 and 3) are then handled by the ordinary
mechanism: the session records which artefact it loaded, the push declares that
parent, and the server rejects a push whose parent is no longer the head. On
conflict the person is told what changed, rather than silently producing a
sibling. There is no containment set, no lineage walk, and no set difference,
because there is no fragment store to reconcile -- the model is the model.

## 6. Footprint

This is the section for reviewers who do not consider PanDDA core, and the
objection in §3 made it smaller: with no server-side composer there is no
fragment store, no containment set and no lineage walk.

**Added**

| | |
|---|---|
| per-site modelling record | new table: dataset, site, event of origin, and the artefact it was last reconciled against. Campaign-scoped, beside `CampaignSite` / `SiteEvaluation`. |
| campaign head resolver | new function: the latest campaign coordinate artefact for a dataset. Read-only over existing tables. |
| reconciliation on push | new code on the existing push path, plus a parent check. |

**Not touched**

* **No column is added to `File`, `Job` or `FileUse`.** Every new foreign key lives
  on a new table; existing models gain at most a reverse accessor, which costs
  nothing and changes no query.
* **No migration alters an existing table.** New tables only -- so no risk to
  DDU's production data, the failure mode that bit migration 0024.
* **No generic task is modified, and none is added.** The push already creates a
  `coordinate_selector` job; this changes what is recorded about that job, not
  what the job is.
* **No change to `REFINE_TASK_NAMES`**, or to how any non-campaign view resolves a
  model. The campaign head is a separate question, asked by campaign views only.
* **No change to the Coot editing loop.** "Add ligand here" is untouched, and that
  is now a design requirement rather than an accident.

## 7. Decisions needed

1. **Does a push become the campaign head immediately, or on an explicit
   "commit"?** Immediately is simpler and matches what a person expects after
   pressing the button. The cost is that an exploratory push -- someone trying
   something and thinking better of it -- moves the head. **Leaning:** immediate,
   with an undo that pushes the previous artefact back, rather than a second
   ceremony before every save.

2. **What the reconciliation does when it disagrees with the origin record.** A
   ligand present at a site with no recorded origin is fine and needs no comment.
   A recorded origin with no ligand present means the person removed it:
   withdraw silently, or tell them what was withdrawn? **Leaning:** tell them, once,
   in the panel -- silent withdrawal of a decision is the behaviour Reinspect
   existed to stop.

3. **Tolerance for "a ligand is at this site".** A site is a point and a ligand is
   a cloud of atoms; the test needs a radius, and adjacent subsites make it
   matter. Probably the same environment radius the site scene already uses for
   pocket residues, so a person sees the same neighbourhood the check uses.

4. **Dictionaries, now a small question.** Composition never crosses datasets, and
   a dataset is one crystal soaked with one compound, so its ligands are copies of
   the same code and use the dictionary that project already holds. A co-frag soak
   is bounded at two codes (`DRG` + `LIG`), which the existing ingest models. The
   only open case is a dataset with no dictionary at all: refuse the push, or
   accept it and let refinement fail with a legible reason.

5. **Chain and residue numbering.** Several copies of one ligand code is ordinary
   crystallography, but they need deterministic assignment or refinement sees
   duplicates. Note this now happens in **Coot**, client side, as part of placing
   the pose -- not in a server-side composer.

## 8. Traps

**Frame -- avoidable by construction, if the site view stays read-mostly.** In
the site view every pose is superposed into the exemplar's frame for display. A
pose nudged *there* and taken into a model would be read in transformed
coordinates and would need the inverse fit applied before storage -- the same
class of error as the site-origin negation, and it would look entirely correct in
the view that produced it.

The way not to have this problem is to not put a modelling commit in the
superposed view at all. The site views judge and dispatch; modelling happens in
the dataset view, in the dataset's own frame. Then no coordinate ever crosses a
frame boundary on the write path, and the inverse fit is needed nowhere.

This is the strongest argument for the view model in
[`moorhen-view-model.md`](moorhen-view-model.md): it is not tidiness, it removes
a whole class of silent error. If a modelling commit is ever added to a
superposed view, the inverse fit becomes mandatory and this trap returns.

**A push is a whole model, so it can quietly lose work.** Pushing from a session
that loaded a stale head replaces newer content with older. The parent check (§5)
is what makes that visible; without it the failure is silent and looks like
someone else's edit vanishing.

**Two people at one dataset.** Site-ordered work makes collisions likelier than
dataset-ordered work did, because two people can be at the same site in the same
campaign at once. The parent check catches it; the message needs to name who.

## 9. Not in scope

Adoption of a run site as a `CampaignSite`; the side panel; any change to how
refinement is run or to what counts as the model of record.
