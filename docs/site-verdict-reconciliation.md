# Site verdicts from more than one reviewer

A design note for letting several people score the same campaign site in the
same dataset, finding the cells where they disagree, and recording how each
disagreement was settled.

**Baseline:** commit `a531c8852` (branch `django`), surveyed 2026-10-13.
**Status:** design note, nothing built. The five decisions were taken on
2026-10-13; see [Decisions](#decisions).

Companion documents:

- [docs/pandda-campaign-design.md](pandda-campaign-design.md), §9 and §9.3:
  `CampaignSite`, `SiteEvaluation`, the two invariants, and the dataset ×
  site matrix this note changes.
- [docs/PROJECT_RECOVERY.md](PROJECT_RECOVERY.md): how verdicts are written
  into each member's `DATABASE.db.xml` snapshot and restored from it.

## The problem

Rescoring is part of how a campaign is meant to be read: a second person goes
through the sites a first person scored, and the cells where they differ are
the ones worth discussing. These are also the best examples to show a new
reviewer. The current model cannot represent any of that.

1. **One verdict per cell.** `SiteEvaluation` is unique on
   `(project, site)`, and the write endpoint is an `update_or_create` on that
   pair. A rescorer's verdict *replaces* the first one. The original call, and
   the fact that there was a disagreement at all, are gone.
2. **No reviewer is recorded.** The model has an `evaluator` field, but
   nothing sets it. The endpoint takes `evaluator` from the request body, and
   the client (`campaign-page-client.tsx`, `handleSetVerdict`) never sends
   one. Every verdict recorded so far has `evaluator = ""`.

The obvious alternative would be a fourth verdict value, "disputed". That
would not do it. It mixes a verdict with where the cell is in the review,
still throws away the two opinions that need reconciling, and depends on
somebody remembering to set it.

## The proposal in one paragraph

Keep **one verdict per reviewer** (`SiteEvaluation` unique on
`(project, site, evaluator)`). The server, not the client, records who the
reviewer is. **Derive** a cell's state from its verdicts: a cell where
reviewers differ is *to reconcile*, without anyone flagging it. Settling a
disagreement writes a **resolution**, a separate human record that becomes
the cell's verdict. Everything that today reads "the" verdict reads the
cell's **effective verdict** instead.

## The model

### `SiteEvaluation`: one row per reviewer

| | Today | Proposed |
|---|---|---|
| Uniqueness | `(project, site)` | `(project, site, evaluator)` |
| `evaluator` | free text from the request body, never sent | the authenticated user, set by the server |
| `verdict`, `note`, `evaluated_at` | unchanged | unchanged |

The docstring's "one person's verdict on one dataset at one site" becomes
true rather than aspirational.

**Who the reviewer is.** Both deployments already authenticate a user on every
request:

- **Desktop:** `LocalSessionAuthMiddleware` resolves the OS user who launched
  the app, with an email computed by Electron (`CCP4I2_LOCAL_USER_EMAIL`).
- **Hosted:** `AzureADAuthMiddleware` resolves the Azure AD user from the
  token's claims.

So the write endpoint sets `evaluator = request.user.email` and ignores any
`evaluator` in the body. A reviewer can then neither be impersonated nor
forgotten, and the client has nothing to get wrong. The column stays a string
(an email) rather than a foreign key to `User`: a verdict travels in a
project snapshot to other installations, where that user row does not exist.

### `SiteResolution`: how a disagreement was settled

A new model, unique on `(project, site)`:

| Field | Meaning |
|---|---|
| `project`, `site` | the cell |
| `verdict` | `hit` / `empty` / `unclear`, the settled call |
| `resolved_by` | the authenticated user, set by the server |
| `resolved_at` | when |
| `note` | why, which is the point of an exemplar |

A separate model, not a flag on `SiteEvaluation`, because a resolution is
not one reviewer's opinion. It must not count as one more vote, and
withdrawing a reviewer's verdict must not withdraw it. **Invariant 1 of the
campaign design (§9.1) extends to it: nothing automated writes a
resolution.**

### A cell's state, derived

| Verdicts at the cell | Resolution | State | Effective verdict |
|---|---|---|---|
| none | none | **not evaluated** | none |
| one or more, all the same | none | **agreed** | that verdict |
| two or more, not all the same | none | **to reconcile** | none (see below) |
| any | present | **resolved** | the resolution's |

"Verdicts" here means attributed ones whenever the cell has any. An
unattributed verdict (one recorded before reviewers were tracked) counts only
at a cell with no attributed verdict. See
[Existing unattributed verdicts](#existing-unattributed-verdicts).

Note:

- **One reviewer is "agreed".** A cell scored once is in the same state as
  today, so a campaign nobody has rescored looks exactly as it does now.
- **A resolution stands when verdicts change later.** But if any verdict is
  newer than the resolution and differs from it, the cell carries a flag,
  *changed since resolved*. The resolution is not withdrawn automatically:
  that would be the machine overriding a person.
- **A "to reconcile" cell has no effective verdict.** Consumers that need one
  (below) treat it as `unclear`, which is honest: the reviewers have not
  settled it.

Any two different verdicts are a difference, so a hit beside an unclear is
"to reconcile" too (decision 1).

## What changes, consumer by consumer

Everything that reads "the" verdict today, and what it reads afterwards.

| Where | Today | Afterwards |
|---|---|---|
| `PUT sites/{site}/evaluation/{project}` | upserts the cell's verdict | upserts **the caller's** verdict; `evaluator` comes from `request.user` |
| `DELETE` on the same path | deletes the cell's verdict | withdraws **the caller's** verdict only |
| `GET evaluations/{project}` | one verdict per site | per site: every reviewer's verdict, the resolution, the derived state, and which verdict is the caller's |
| new: `PUT` / `DELETE sites/{site}/resolution/{project}` | — | settle a cell, or reopen it |
| `member_projects` → `site_cells` (`lib/campaign_events.py`) | `verdict` | `verdict` = the **effective** verdict, so an older client keeps working; adds `state`, `verdicts` (`[{evaluator, verdict}]`) and `mine` |
| `member_projects` → `site_evaluations` (the chips fallback) | hit/unclear verdicts | hit/unclear **effective** verdicts; a "to reconcile" cell counts as unclear |
| `sites_evaluated` ("Eval." column) | cells with a verdict | cells **the viewer** has scored, with everyone's count in the tooltip once there is more than one reviewer (decision 4) |
| Site and summary scenes (`lib/campaign_scene.py`): membership is the verdict | hit verdicts; `include=unclear` adds unclears | the same rule on the **effective** verdict, so "to reconcile" cells are drawn only with `include=unclear` |
| Site delete warning (`evaluation_count`) | counts verdicts | counts verdicts and resolutions |
| `delete_campaign` | deletes evaluations with sites | also deletes resolutions (cascade) |
| Snapshot export (`db/export_project.py`) | `siteevaluationTable`, by field | the same, so `evaluator` rides through; adds `siteresolutionTable` in the evaluated project's snapshot |
| Snapshot import (`db/import_i2xml.py`) | `update_or_create(project, site)` | `update_or_create(project, site, evaluator)`; resolutions on `(project, site)` |

The snapshot serialiser already writes campaign rows by field rather than by
a column list (the property §9.1 asked §14.0 for), so `evaluator` needs no
serialiser change. Only the import's lookup key moves.

**Mixed versions.** An older installation importing a newer snapshot keys on
`(project, site)`, so several reviewers' rows collapse to the last one in the
file. That loses information but not the project, and it only happens when a
snapshot travels backwards between versions. It should be stated in the
release note, not engineered around.

## The interface

### The campaign Moorhen panel (where verdicts are recorded)

- The three verdict buttons record **my** verdict, as now. Clicking the
  active one withdraws it, as now.
- Beneath them, the other reviewers' verdicts at this cell, as small chips:
  `JS · empty`, with the full email and the time in the tooltip.
- **Blind scoring is a switch** (decision 2). With it on, other reviewers'
  verdicts at a cell stay hidden until you have recorded your own there, so a
  rescore is not anchored by the first score. With it off, they show at once,
  which is what a discussion needs. The switch is per viewer and remembered in
  the browser, not a campaign setting, and it defaults to on. The matrix
  follows it too: with blinding on, a cell you have not scored shows only
  whether *anyone* has scored it, not what they said.
- In a "to reconcile" cell, a **Resolve** control: the three verdicts again,
  plus a note field. Anyone who can score can resolve (decision 3), and the
  resolution records who did. A resolved cell shows the resolution, who
  resolved it and the note, and a **Reopen** action.

### The overview matrix

The glyph rule from #811 extends: every state can be read without colour.

| State | Box | Mark |
|---|---|---|
| agreed / resolved | the verdict's colour, as now | ✓ / ✕ / ? |
| to reconcile | a colour outside red/amber/green (MUI `secondary`) | ≠ |
| resolved, changed since | the resolution's colour | the verdict's mark, with a small corner dot |

Fill keeps meaning "PanDDA found an event here". The tooltip lists each
reviewer's verdict, and the resolution if there is one.

### Finding the exemplars

- A **To reconcile** filter on the overview: only rows with at least one cell
  to reconcile. The page header counts them: "7 cells to reconcile".
- In the Moorhen campaign page, **Next to reconcile**, which steps to the next
  such cell across datasets, so working through them takes one click each.

## Existing unattributed verdicts

Every verdict recorded before PR 1 has `evaluator = ""`. Their history is
lost and is not reconstructed (decision 5). Nobody can say who gave them, and
code must not guess.

They are **kept, and outranked**:

- At a cell where **only** unattributed verdicts exist, the cell behaves as
  today: the unattributed verdict is the effective verdict, so the matrix of
  a campaign scored before this change looks the same afterwards.
- The first **attributed** verdict at a cell is treated as though it were the
  first verdict there. From then on the cell's state is derived from
  attributed verdicts only, so an earlier unattributed verdict can never put a
  cell "to reconcile", including against the same person's new verdict.
- The unattributed row is not deleted or rewritten. It stays in the database
  and the snapshot, and the panel and the cell tooltip show it as *earlier,
  unattributed: hit*, for the record.

No `attribute_verdicts` command and no data migration are needed.

## Agents

The MCP facade has no verdict tool today, and invariant 1 forbids one that
writes a verdict. This note adds a reason to keep it that way. An agent drives
the API with the desktop user's session token, so a verdict it wrote would be
stamped with the person's email and be indistinguishable from theirs. If an
agent ever proposes verdicts, it does so somewhere else (a receipt's KPIs, as
PanDDA's `interesting` is), never in `SiteEvaluation`.

## Delivery

Three PRs, each useful alone.

1. **Attribution (small, first).** The write endpoint sets `evaluator` from
   `request.user`. No schema change. From the moment it ships, verdicts are
   attributed, so fewer unattributed ones pile up while the rest is built.
2. **Model and API.** The migration (uniqueness to
   `(project, site, evaluator)`, the `SiteResolution` model), the derived
   state (with unattributed verdicts outranked), the endpoint changes,
   snapshot export and import, and scene membership on the effective verdict.
3. **Interface.** The panel's other reviewers' verdicts and the blind-scoring
   switch, Resolve and Reopen, the matrix's ≠ state, the per-viewer "Eval."
   count, the To reconcile filter and Next to reconcile.

### Tests

- **Migration, seeded** (per the project's rule that a migration passing on an
  empty database proves nothing): existing rows survive with their verdicts,
  and a second reviewer's row can then be added at the same cell.
- **API:** two authenticated users scoring the same cell give two rows and
  "to reconcile"; one withdrawing leaves the other's; a resolution sets the
  effective verdict and survives a later change, which raises *changed since
  resolved*; a body-supplied `evaluator` is ignored; an unattributed verdict
  is the effective verdict alone, but never makes a cell "to reconcile" once
  an attributed verdict exists.
- **Snapshot round trip** with two reviewers and a resolution.
- **Scenes:** a "to reconcile" cell is drawn only with `include=unclear`.
- **Client:** the cell-state rules in `lib/site-matrix.ts` (pure, as today),
  and the ≠ glyph rendered.

## Decisions

Taken 2026-10-13.

1. **Any difference is a disagreement.** Hit beside unclear is "to
   reconcile", like hit beside empty. The overview can still sort
   hit-against-empty first.
2. **Blind scoring is a per-viewer switch,** on by default. See
   [The interface](#the-interface).
3. **Anyone may resolve,** and the resolution records who acted as resolver
   (`resolved_by`).
4. **"Eval." counts the viewer's own scoring,** with everyone's count in the
   tooltip.
5. **History before attribution is accepted as lost.** Unattributed verdicts
   are kept, and an incoming attributed verdict is treated as though it were
   the first. See
   [Existing unattributed verdicts](#existing-unattributed-verdicts).
