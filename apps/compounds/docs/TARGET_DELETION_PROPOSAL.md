# Admin-only Target Deletion

**Status**: Proposal
**Scope**: Adds an admin-gated way to delete a Target from the registry. Because `Compound.target` is `on_delete=PROTECT`, a plain delete currently raises `ProtectedError` whenever a Target has any compounds. This proposal introduces **two** delete modes, both admin-only:

1. **Reassign mode** (default / recommended): admin picks a replacement Target, the server bulk-reassigns affected compounds in a single transaction, then deletes the original.
2. **Cascade mode** (destructive, opt-in behind extra confirmation): affected compounds — and everything that cascades from them (batches, QC files, compound documents, etc.) — are deleted along with the target.

Reassign is the primary path because it's non-destructive; cascade exists for genuinely bad imports (a whole campaign's worth of junk rows landed under a garbage target name) where reassigning makes no sense.

## Motivation

Targets accumulate over time (typos, experimental renames, merged campaigns). Admins currently have no way to remove a Target without dropping into the DB. The ELN bulk-import flow also makes it easy to create an unwanted Target *and* fill it with compounds by supplying an unrecognised name in the spreadsheet — right now the only mitigation is to manually revert via the DB. Reassign handles the "wrong-target" case; cascade handles the "wholesale-bad-import" case.

## Related FK landscape

`Target` is referenced by several models. The only one that blocks deletion is `Compound`; everything else already has sensible `on_delete` behaviour and does **not** need special handling in this proposal.

| Model | FK | `on_delete` | Behaviour today |
|---|---|---|---|
| `Compound.target` | PROTECT | PROTECT | **Blocks delete — the only blocker, and the thing this proposal reassigns.** |
| `CompoundTemplate.target` | CASCADE | CASCADE | Templates belong to their target — cascade is correct. |
| `Assay.target` (line 383) | SET_NULL | SET_NULL | Historical assays become target-less — acceptable. |
| `Protocol.target` (line 498) | SET_NULL | SET_NULL | Same. |
| `Hypothesis.target` (line 737) | CASCADE | CASCADE | Hypotheses belong to their target. |
| `Target.parent` (self-FK) | SET_NULL | SET_NULL | Children become root-level — acceptable. |

So reassign mode only needs to reason about `Compound`; everything else the ORM already handles correctly via existing `on_delete` rules.

**Cascade mode implications**. Deleting a Compound is itself cascading: `Batch` (CASCADE from `Compound`), `BatchQCFile` (CASCADE from `Batch`), `CompoundDocument` (CASCADE from `Compound`), `MolecularProperties` (OneToOne from `Compound`). Any assay/data-series relations that FK to Compound or Batch with CASCADE also go. The preview endpoint must surface these second-order counts so the admin sees the full blast radius before confirming. Files on blob storage (QC uploads, compound documents) are removed by the existing `post_delete` signal handlers on `BatchQCFile` / `CompoundDocument`.

## Backend

### New endpoints

All three are admin-only (`is_platform_admin`). Two distinct delete endpoints rather than one mode-switched endpoint — the URL itself expresses intent, which makes audit logs grep-able and makes it harder to click the wrong one by accident.

1. **`GET /api/compounds/targets/<id>/deletion_preview/`** — returns counts for the modal, covering both direct FKs and cascade-mode second-order impact:
   ```json
   {
     "target": {"id": "...", "name": "NCL_ARd"},
     "compound_count": 47,
     "template_count": 3,
     "hypothesis_count": 1,
     "assay_count": 12,
     "protocol_count": 4,
     "child_target_count": 0,
     "cascade_impact": {
       "batches": 89,
       "batch_qc_files": 145,
       "compound_documents": 23,
       "assay_data_series": 312
     },
     "reassignable_targets": [
       {"id": "...", "name": "NCL_CDK2A"},
       ...
     ],
     "affected_sample": {
       "compounds": ["NCL-00026042", "NCL-00026043", "..."]
     }
   }
   ```
   `cascade_impact` is what the admin sees in cascade mode; `reassignable_targets` excludes the one being deleted. `affected_sample.compounds` lists the first 5 `formatted_id`s so the admin can eyeball what they're about to touch. Paginated if we ever have many targets; currently under 50 everywhere so no pagination in v1.

2. **`POST /api/compounds/targets/<id>/delete_with_reassign/`** — body:
   ```json
   {"replacement_target_id": "<uuid>"}
   ```
   Wrapped in `transaction.atomic`. Steps:
   1. Load both Target objects with `select_for_update()`.
   2. Validate: replacement exists, is not the same as the one being deleted, is not a descendant of it (to avoid breaking any parent-chain invariants).
   3. `Compound.objects.filter(target=target).update(target=replacement)` — bulk, no per-row ORM overhead.
   4. `target.delete()` — at this point `CompoundTemplate`, `Hypothesis`, etc. cascade/null as per existing FK rules.
   5. Emit a `reversion` revision so the reassignment can be inspected historically.
   Returns `{"compounds_reassigned": N, "target_deleted": true}`.

3. **`POST /api/compounds/targets/<id>/delete_with_cascade/`** — body:
   ```json
   {"confirmed_compound_count": 47}
   ```
   Wrapped in `transaction.atomic`. Steps:
   1. Load target with `select_for_update()`.
   2. Count current compounds; if it doesn't match `confirmed_compound_count`, return `409 {"error": "compound count changed — reload the preview"}`. This guards against another user adding compounds between the admin opening the modal and clicking confirm.
   3. Delete compounds in a loop (or `.delete()` on the queryset — ORM handles cascade to `Batch` → `BatchQCFile`, `CompoundDocument`, etc.). The `post_delete` signals on file-bearing models clean up blob storage.
   4. `target.delete()` — remaining FKs (templates, hypotheses, child targets) clean up as per existing rules.
   5. Emit a `reversion` revision summarising the cascade for audit.
   Returns `{"compounds_deleted": N, "batches_deleted": M, "files_removed": K, "target_deleted": true}`.

### Permission

Reuse the existing `IsPlatformAdmin` permission class (or equivalent — confirm name when implementing). Do **not** gate any of these on `canContribute`; deletion is strictly admin. Both delete endpoints share the same permission class; cascade's extra safety is in the payload (`confirmed_compound_count`) and the frontend confirmation UX, not in a weaker permission.

### Validation rules

For reassign:
- `400` if `replacement_target_id` missing, same as target being deleted, or unknown.
- `400` if replacement is a descendant of the target (parent FK chain would be broken).
- `404` if target id doesn't exist.
- `409` (or bespoke error) if another admin is mid-reassignment (`select_for_update()` should serialize, so this is mostly academic).

For cascade:
- `400` if `confirmed_compound_count` missing or not an integer.
- `409` if `confirmed_compound_count` doesn't match the current count (user is operating on stale data).
- `404` if target id doesn't exist.

## Frontend

### Entry point

On `/registry/targets` list page, add a `Delete` icon button at the end of each row, visible **only when `canAdminister` is true**. Tooltip when disabled: "Admin operating level required".

### Modal (new component, `TargetDeletionDialog.tsx`)

On open, calls `deletion_preview`. The modal has a **mode toggle** at the top — a two-option `ToggleButtonGroup` with "Reassign" (default, selected) and "Delete everything" (destructive-styled, red). Picking a mode changes the body below.

**Shared header**:
- `Delete target "NCL_ARd"?`
- Impact summary (always visible so admin sees both framings):
  - `47 compounds`, `89 batches`, `145 QC files`, `23 compound documents` — raw counts from `deletion_preview`.
  - `3 compound templates will be deleted (belong to this target)` *(only if >0)*
  - `12 assays + 4 protocols will have their target cleared` *(only if >0)*
  - `1 hypothesis will be deleted` *(only if >0)*
- Preview list: first 5 affected compound `formatted_id`s + "… and 42 more" if truncated.

**Reassign mode body**:
- **Replacement target selector** (required `<Autocomplete>`): the list from `reassignable_targets`.
- **Typed confirmation**: "Type `NCL_ARd` to confirm" — input must match exactly before the Delete button enables. Protects against slip-of-the-click.
- **Action buttons**: Cancel, Reassign & Delete (primary).

**Cascade mode body** (destructive-styled — red banner, amber icons):
- **Prominent warning banner**: "This will permanently delete 47 compounds, 89 batches, 145 QC files, and 23 documents. This cannot be undone."
- **Two confirmations required** (both must be satisfied before the Delete button enables):
  1. **Typed target name**: "Type `NCL_ARd` to confirm" — same as reassign.
  2. **Typed compound count**: "Type `47` to confirm the number of compounds to be deleted". Using the number from the preview makes stale-data bugs obvious — if another admin added a compound in between, the user sees the updated count and types the old one, which won't match.
  3. (Optional) A checkbox: "I have verified this target's compounds cannot be reassigned."
- **Action buttons**: Cancel, **Delete Everything** (destructive red). Button label is deliberately dire.

On Reassign submit:
1. POST to `delete_with_reassign/` with `{replacement_target_id}`.
2. On 200: toast `"Reassigned N compounds to <target> and deleted NCL_ARd"`, close modal, refresh target list (`mutate('targets/')`, plus any compound/template caches that show target name).
3. On error: surface the server message in an Alert at the bottom of the modal without closing it.

On Cascade submit:
1. POST to `delete_with_cascade/` with `{confirmed_compound_count}`.
2. On 200: toast `"Deleted N compounds, M batches, and target NCL_ARd"`.
3. On `409` (count mismatch): reload the preview, highlight the new count in red, force the admin to re-type the corrected value. Do not let them submit with stale data.
4. Other errors: surface in-modal Alert.

### Navigating away

If an admin is currently on the detail page of the target being deleted, the page should redirect to the target list after a successful delete (handled naturally if the admin deletes from the list view; detail-page delete button is a nice-to-have for v2).

## What's intentionally **out of scope for v1**

- **Soft delete / archive**. If admins want the target hidden but not gone, we can add an `is_archived` flag later. Both reassign and cascade delete remain permanent.
- **Delete-with-orphaning** (NULL the compound target instead of reassigning or deleting). Deferred — it breaks the existing assumption in the UI and aggregation queries that `compound.target_name` is always populated. If the admin wants something "in between" today, they can create a catch-all Target (e.g. "Unassigned") and reassign to it.
- **Bulk target deletion**. One at a time; bulk is rarely needed and hugely more dangerous.
- **Non-admin contributor deletion of empty targets**. Keep it admin-only even when `compound_count == 0`, to reduce the audit surface.
- **Per-compound cherry-pick in cascade mode** (e.g. "delete these 30 compounds, reassign the other 17"). If an admin needs that, they should first filter the compound list, bulk-reassign the keepers to a different target, then cascade-delete what remains.

## Open questions

1. Should the preview endpoint also list the first N compound `formatted_id`s in each bucket, so the admin can sanity-check? Probably yes — cheap to compute, makes the modal more confidence-inspiring. Add `affected_sample: {compounds: [...], templates: [...]}` with first 5 each.
2. Audit trail: `reversion` should pick this up automatically for models that register with it. Worth a quick check during implementation that `Compound.target` changes create a revision.
3. Should the replacement Target dropdown be filtered to only targets the admin is a member of (group-based ownership)? Not relevant on kawamura (everyone's in one group), but on shared instances it would matter. Flag for later.
