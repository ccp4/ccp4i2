# Campaign sites, evaluations, and the merged SubstituteLigand route

**Branch:** `campaign-site-verdicts-ui` — PR #564, on top of #563 (merged)
**Written:** 2026-09-20, against `django` at a70. Updated 2026-09-21.

This is a handoff note. It says what is finished, what is half-done, what was
learned the hard way, and what to do next.

---

## The shape of the work

Two threads, which became entangled because the second blocked testing the first.

1. **Campaign sites and evaluations** — a data model for recording *what was
   found where* in a fragment campaign, replacing the practice of tagging a
   project with a site's name.
2. **The merged SubstituteLigand route** — making a fragment campaign's jobs
   actually run, which turned out to need a free-R reconciliation step.

---

## Thread 1: sites and evaluations

### Landed in PR #563 (`campaign-sites-evaluations`)

* **`CampaignSite`** — sites were a JSON list on `ProjectGroup`, replaced
  wholesale on every write. That gave them no identity: a site was found by
  its position and its `name`, so renaming one orphaned every reference, and
  two people editing a campaign silently lost each other's work. Sites are now
  rows with stable ids; endpoints are `GET`/`POST` on the collection and
  `PATCH`/`DELETE` on `sites/<id>/`.

* **`SiteEvaluation`** — one verdict (`hit` / `empty` / `unclear`) for one
  dataset at one site. Per site, because a dataset is routinely a hit at one
  and empty at another. **No confidence field**, deliberately: the three
  verdicts carry the uncertainty.

* **The absence of a row is the "not yet evaluated" state.** This is the whole
  reason for a table rather than a tag: `empty` asserts somebody looked and
  found nothing, where no row means nobody looked. An untagged project was
  both at once.

* **Migration 0024** moves sites out of JSON (reversible; skips a site with no
  usable origin rather than inventing one, disambiguates duplicate names).
  **Migration 0025** converts site tags to `hit` evaluations — tagging a
  dataset at a site asserted a hit, so `hit` is the faithful reading — marks
  them `evaluator="migrated"`, and consumes the tag. The reverse restores the
  tags and deletes only rows carrying that marker.

* **`make_demo_campaign`** — builds a throwaway campaign from a real BAZ2B
  fragment series (5E9I, 5DYU, 5E9K, 5E9L, 5E9M, 5E9Y) via PDBe. Six datasets,
  same space group, cells within ~2 Å, **a different fragment in each**.

* **`delete_campaign <name>`** — the inverse: removes the group, its sites and
  evaluations, the parent and member projects, and their directories. Asks for
  the name to be typed back (`--yes` skips that, `--dry-run` only reports,
  `--keep-files` leaves the directories). A project that is also in another
  group is kept, and it refuses while jobs are queued or running.

### On this branch, not yet in a PR

* **Evaluation endpoints** (`67a3d8753`). `PUT`/`DELETE` on
  `sites/<id>/evaluation/<project>/`. The findings ride in `member_projects`
  rather than a separate endpoint, so the overview stays one request rather
  than N. Only `hit` and `unclear` are listed — a rich campaign has 30–40
  sites and most are empty, so sending every verdict would put a 40-entry list
  on every row to render two chips. `sites_evaluated` / `sites_total` give the
  completeness count. 14 tests.

### Landed on this branch: the UI

Built to the design below, which is left in place because it records *why*
each choice was made rather than only what was built.

* A **Sites column** in the campaign overview, showing chips **only for `hit`
  and `unclear`**. Empty and unevaluated render nothing.
* **`hit`** green filled; **`unclear`** amber outlined — deliberately
  downplayed.
* Hits sort before unclears so they are never pushed out of a capped list.
* A **`12/40` count shown only while incomplete**; denominator is the
  campaign's current site count, so it is comparable down the column.
* A chip click navigates to
  `/ccp4i2/moorhen-page/campaign/<id>?job=<jobId>&site=<siteId>`. The campaign
  page already reads `view`, `job` and `summary` and already loads sites, so
  the `site` param is a small addition.
* **A project with no job cannot have been evaluated**, so it shows no chips
  and the click case never arises.
* The **verdict control replaces the site-tag button** in the Moorhen campaign
  panel. Tag and verdict should not coexist: two records of one fact, with
  only one of them maintained.

Three things came out of building it that the design did not anticipate:

* **The client was still on the old sites API.** PR #563 moved sites to rows
  with ids and per-site endpoints, but the client went on PUTting the whole
  list back, so adding, renaming, moving and deleting a site were all broken
  on this branch. Writes now go one site at a time, addressed by id.
* **The overview payload cannot drive the verdict control.** It carries only
  hits and unclears by design, and a control that *records* verdicts has to
  distinguish "looked, found nothing" from "not looked at yet" or it writes
  the wrong one back. Added `GET evaluations/<project>/`: every verdict for
  one dataset, empties included. The overview stays lean; the per-dataset
  view can afford the whole picture.
* **Which job a chip opens** was unspecified. It opens the dataset's latest
  *finished* top-level job, falling back to its latest job of any status so a
  running or failed dataset still opens on something.

### Also requested, not started

A **"place ligand here"** button in the Moorhen campaign page. Today, spotting
a ligand at a site means opening Moorhen's generic "Get monomer" dialog,
selecting the protein as the dictionary source, typing the three-letter code
and clicking OK. A direct Coot API call would do it in one action.

---

## Thread 2: the merged route

### The problem

Campaign jobs failed part-way through, inside `i2Dimple`, with
`MtzMergeError: Incompatible unit cells`. Not a def.xml `sameCrystal`
assertion — a runtime check inside `merge_mtz_files`.

On the merged route (`OBSAS=MERGED`) aimless never runs, so the free-R set is
whatever the user supplies. In a campaign that is the reference crystal's: its
cell differs by a percent or more, and it need not cover this dataset.

### The fix (`6b1f2b509`)

**Reconcile once, up front**, rather than teaching every downstream merge to
ignore the difference. `freerflag` in `COMPLETE` mode does both halves:

* joins by **reflection index**, so existing flags keep the reflections they
  were assigned to;
* stamps the **data's cell**;
* **fills in flags** for reflections the input lacks.

Re-stamping a cell by hand does only the first half — which is why an earlier
`recell_free_r` helper was written and then deleted. `freerflag` is
gemmi-native, so this costs no subprocess, and the result is harvested to
`FREERFLAG_OUT`: a real project file the next job can use.

**It runs whenever a free set is supplied**, not only on cell mismatch.
Gamma's demo free-R set shares the data's cell and still leaves 264 observed
reflections unflagged; `COMPLETE` is a no-op when the set already fits.

### The bug that cost the most time

```python
self.freerToUse = someOtherCDataFile   # does NOT rebind
```

Once the attribute holds a `CDataFile`, assignment is intercepted by CData's
`__setattr__` and **coerced into the existing object**, which keeps its own
`baseName` and `relPath`. So the reconciled set was created correctly,
harvested correctly, and handed downstream as the original input — while the
code plainly said otherwise, and the annotation even came across.

What finally showed it: printing `id()` of both objects after the assignment.
They differed.

Reassignments now go through `_useFile()`, which uses `object.__setattr__` and
documents the trap. **The first assignment was always fine** (the attribute
starts as `None`, nothing to coerce into) — it is every reassignment after
that which silently does nothing.

**This pattern appears at about a dozen sites in `SubstituteLigand.py` alone,
and very likely elsewhere in the codebase.** Worth a sweep.

### Tests

`OBSAS=MERGED` had **no end-to-end coverage at all** — every existing
SubstituteLigand i2run test supplied `UNMERGEDFILES`. That is why so much
accumulated in it unseen. Two new tests:

* `test_substitute_ligand_merged_data` — the route runs.
* `test_merged_free_r_other_crystal` — builds a free-R set stamped ~2 Å out,
  asserts it really does fail the strict check, then that the job refines
  anyway and publishes a `FREERFLAG_OUT` carrying the data's cell.

All 5 tests in the file pass (~5 min). Unit suite 2570.

---

## Things learned that are not obvious from the code

* **`SMILES` vs `SMILESIN`.** `inputData` has both. `SMILESIN` is what the GUI
  binds and what `processInputFiles` reads; `SMILES` is legacy with a
  hardcoded default. Writing to the wrong one leaves the Run dialog reporting
  that no SMILES was given over a job that looks configured.

* **`set_parameter` returns a `Result`, it does not raise.** Unchecked calls
  leave the job holding its default. Control parameters also need a
  `container.` prefix that the file-import path does not.

* **Diagnostics propagate correctly.** An earlier claim here that they did not
  was wrong — the parent's `diagnostic.xml` does carry the child's cause. A
  regex matching only `<description>` misses the `<details>` where it lives.
  Minor nits remain: the parent's own summary is listed before the child's
  cause, so the least informative line is shown first, and
  "i2Dimple pipeline failed" is recorded twice.

* **PDB-REDO was unreachable** while this was written — its root 404s,
  including for the `8xfm` entry the i2run suite fetches. Those e2e tests will
  be failing for reasons unrelated to any change here. `make_demo_campaign`
  uses PDBe instead.

* **The merged route writes no `F_SIGF_OUT`** — the observations pass through
  unchanged, so there is nothing to re-export. `FREERFLAG_OUT` *is* written,
  because the free set is reconciled. Anything that goes looking for a
  member's reflections must therefore accept the imported `F_SIGF_IN`
  (`directory=2`) when there is no output: the campaign Moorhen page's
  "Run servalcat refinement" asked for `directory=1` only, and reported "no
  reflection data" for every demo member until that was relaxed. It looked
  like a consequence of the demo parent being empty, and was not.

* **`make_demo_campaign` populates the parent** with a `coordinate_selector`
  job (XYZOUT, protein only) and a `freerflag` job (FREEROUT), as the campaign
  page's import dialogs do. `parent_files` recognises reference data by
  exactly those param names and types, so an empty parent shows the campaign
  as unconfigured. Member jobs take their model and free set from those files.

* **`merge_mtz_files` builds a complete reflection list.** It takes the unique
  set for the resolution range and unions in every index every input holds, so
  it does not assume either file is complete. Reflections absent from a source
  get NaN (MNF) in that source's columns — verified: 44 such in a reconciled
  set, none of them observed.

---

## The `__setattr__` sweep (done)

The mechanism is worse than described above, and narrower.

**Worse:** the reassignment does not merely leave the old `baseName` in place
— the new file's path never arrives at all. `CData.__setattr__` dispatches a
CData-into-CData assignment to `_smart_assign_from_cdata`, which for a
CDataFile copies nothing that matters. The attribute is completely unchanged.

**Narrower:** it needs the *containing* object to be an initialised CData.
That is always true of a CPluginScript in a real run, but it means the trap
cannot be reproduced on a plain object, or on a plugin built with `__new__` —
the guard `not hasattr(self, "_hierarchy_initialized")` sends those straight
to the default assignment. A first attempt to reproduce it this way appears to
show no bug at all.

There is a **second silent no-op in the same code path**: assigning a CData
whose values are all unset is skipped entirely, so a stage that produced
nothing leaves the previous stage's file in place rather than clearing it.

An AST sweep of `wrappers/` and `pipelines/` for attributes assigned a
container-derived value more than once found 11 candidates. Six are false
positives — the two assignments are arms of one `if`/`else`, so only one ever
runs (`PrepareDeposit.coordsToUse`, `dr_mr_modelbuild_pipeline.dictToUse`,
and that pipeline's `coordinatesForCoot`/`mapToUse`). The rest were real:

* **`SubstituteLigand`** — Phase 3 (Dimple *or* Phaser RNP) binds
  `finalCoordinates`, `mapToUse` and `coordinatesForCoot`; Phase 4
  (Servalcat, which runs *always*) and Phase 5 (Coot ligand fitting) re-point
  them, and were doing nothing. What that actually cost, traced through to the
  consumers: **Coot fitted its ligand into Phase 3's coordinates and Phase 3's
  map**, not into the refined coordinates and the better map Servalcat had
  just produced — `_runCootLigandFitting` reads `coordinatesForCoot` and
  `mapToUse` directly. The *published* XYZOUT is not affected: Coot writes
  straight to `outputData.XYZOUT`'s path and Servalcat's is harvested into it,
  neither going through these attributes. And the Phase 5 re-point of
  `finalCoordinates` has no consumer at all — a trap lying in wait for whoever
  reads it next, rather than a live fault.
* **`phaser_rnp_pipeline`** — the most consequential. `F_SIGF_TOUSE` and
  `FREERFLAG_TOUSE` are set from the inputs, then re-pointed at pointless's
  reindexed output when the cells disagree; `runPhaser` and `runRefmac` both
  read them afterwards. So when pointless decided a reindex was needed, phaser
  and refmac ran on the original, un-reindexed data anyway — which is the one
  thing that branch exists to prevent.
* **`prosmart_refmac`** — `currentCoordinates` is re-pointed at Coot's output
  after water fitting and handed to the post-Coot Refmac, which was therefore
  refining the pre-Coot coordinates.

`_useFile` is now a method on `CPluginScript` rather than a local helper in
SubstituteLigand, so every plugin has the remedy, and the three files above
use it. Seven tests in `tests/unit/plugins/test_working_file_attributes.py`
pin the trap itself, the unset-source no-op, and that `_useFile` defeats both
— including one that should be deleted loudly if CData's semantics ever
change so that plain assignment rebinds.

## Suggested next steps

1. ~~Land PR #563~~ — merged 2026-09-21 as `e42bb3471`.
2. ~~Split this branch~~ — not done, and no longer worth doing. The split was
   proposed so the cell/reconcile work would not wait on the UI; with the UI
   finished they land together as #564, and separating them now would mean
   rewriting pushed history to save nothing.
3. **Decide on the permissive-merge commits.** Recommendation: **keep them.**
   They fix the failure at a different layer from the reconcile, which only
   runs inside SubstituteLigand — `i2Dimple` and `pointless_reindexToMatch`
   are reusable, and `phaser_pipeline`, `phaser_pipeline_phil` and
   `phaser_rnp_pipeline_phil` all reach the same join without passing through
   the reconcile. Dimple's own argument stands on its own terms too: fitting a
   model that does not quite match its data is what the program is *for*, so a
   1 Å cell test is the wrong default for it. `STRICT_CELL_MATCH` keeps the
   guard reachable. The cost is that data genuinely paired from the wrong
   crystal no longer fails early — worth a look at whether the report says
   loudly enough that the cells differed.
4. ~~Build the Sites column~~ — done, above.
5. ~~Sweep for the `__setattr__` coercion pattern~~ — done, above.
6. **Still not started: the "place ligand here" button** in the Moorhen
   campaign page (a direct Coot API call in place of the generic "Get
   monomer" dialog).
7. **`member-project-row.tsx` is dead code** — exported from the campaigns
   index but rendered nowhere; the virtualized table has its own row. It did
   not get the Sites column. Delete it, or wire it up.
