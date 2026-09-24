# "Place the ligand here" in a fragment campaign

**Status:** implemented 2026-09-21 (client-only, uncommitted at time of
writing). Written 2026-09-21 against `django` at a71. The design is left as
written because it records *why* each choice was made; see *As built* at the
foot for the two places the implementation departed from it.
**Scope:** one button in the campaign Moorhen page. Independent of
[campaign-site-scene-design.md](campaign-site-scene-design.md) — neither
needs the other, and they touch disjoint files.

---

## The problem, counted in interactions

A campaign dataset is open in the viewer. There is density at a site, the
fragment is obviously there, and the user wants it modelled. Today:

1. Open the Moorhen menu, find **Get monomer**.
2. Leave **Source** on *Default*.
3. Pick the source molecule from a dropdown — and it must be the member's own
   coordinates, because that is the molecule the dictionary was attached to.
   Picking anything else silently fails with "missing dictionary".
4. Type the three-letter code.
5. **OK**. The monomer arrives *as a new molecule*, floating at the view centre.
6. **Merge molecules**: choose the new molecule, choose the member's
   coordinates as the destination.
7. Delete the now-empty leftover molecule, or live with it cluttering the list.

Seven interactions, of which three (2, 3, 4) restate something the page
already knows and one (7) is cleaning up after the tool.

**What the page already knows**, at the moment the button would be pressed:

| Fact | Where it lives today |
|---|---|
| Which dataset is selected | `selectedMemberProjectId` in `campaign-moorhen-wrapper.tsx` |
| Which job's coordinates are on screen | `fileSource` / `selectedJobId` |
| Which molecule those coordinates are | `molecules` selector; the molecule's `uniqueId` is its loader URL |
| Which dictionaries are attached to it | fetched in `fetchJobFiles` via `fetchJobDictionaryFiles(jobId)` |
| Which comp_ids those dictionaries define | `extractDictCompIds()` in `lib/moorhen-dictionaries.ts` |
| Where "here" is | the Redux `sceneSettings.origin` |

Every input to the seven steps is already in hand. The button is a wiring
job, not a new capability.

---

## The finding that shrinks this task

**Moorhen already implements the whole sequence as one method.**

`MoorhenMolecule.addLigandOfType(resType, fromMolNo)`
(`baby-gru/src/utils/MoorhenMolecule.ts`, ~line 2032) does:

* `get_monomer_and_position_at(TLC, fromMolNo, -x, -y, -z)` at the current
  view origin;
* if that returns `-1`, `loadMissingMonomer(TLC, fromMolNo)` and retries —
  so a code the attached dictionary does not cover is fetched from the
  monomer library rather than failing;
* wraps the result in a `MoorhenMolecule`, `mergeMolecules([it], doHide=true)`
  — which also calls `transferLigandDicts`, so the merged-in chemistry
  travels with the atoms — and then `delete()`s the leftover.

That is steps 5, 6 and 7 in one call. It is present in the installed build:
`client/node_modules/moorhen/types/src/utils/MoorhenMolecule.d.ts:509`, and
`moorhen.Molecule` is a type alias for `MoorhenMolecule`
(`types/moorhen.d.ts:40`), so molecules taken straight from the Redux
`molecules.moleculeList` selector have it, typed.

**So the remaining design question is not "how do we place a monomer". It is
"which code, into which molecule" — and both of those are questions about
CCP4i2's data model, not Moorhen's.**

---

## Which code: the dictionary, not the coordinates

The obvious server-side answer is wrong, and it is worth writing down why.

`lib/campaign_scene.detect_ligands()` already identifies a dataset's ligand —
but it does it by reading the *refined coordinates* for a fragment-like
residue name. At the moment this button is pressed, **the ligand is not in the
coordinates**. That is the entire reason the button is being pressed. So the
campaign-scene hit detector cannot answer this question, and must not be
reused for it.

The ligand a dataset is *looking for* is recorded in two places:

1. **`SubstituteLigand`'s `LIGAND_CODE` control parameter** (default `DRG`;
   `SubstituteLigand.def.xml`, "Code for the fitted ligand"). This is the
   statement of intent, and it is right even when nothing was ever fitted.
2. **The comp_ids in the job's restraint dictionary** — the acedrg output that
   the refinement took as input. `jobs/{id}/dictionaries/` already returns
   these files, the wrapper already downloads their text in `fetchJobFiles`,
   and `extractDictCompIds()` already parses `data_comp_<X>` out of them.

**Recommendation: use (2), and do not add an endpoint for (1).**

The reasons:

* (2) is already fetched. The wrapper has the dictionary *texts* in a local
  `dictionaries` variable inside `fetchJobFiles` and currently throws them
  away after attaching them. Keeping the comp_ids in state costs one
  `useState`.
* (2) works for datasets that never ran SubstituteLigand — a project whose
  ligand came from a bare acedrg job plus a refinement has no `LIGAND_CODE`
  to read.
* (2) is the code that is *actually in the dictionary Coot holds*, which is
  what `get_monomer_and_position_at` needs. A `LIGAND_CODE` that disagrees
  with the dictionary is a bug we would be propagating rather than catching.
* Reading (1) means a `container.xml` parse per member project, which is the
  expensive path — see `docs/container-construction-defects.md` on
  construction cost.

The failure mode of (2) is that a refmac `LIBOUT` may declare standard
monomers alongside the fragment. Filter with the same rule
`campaign_scene._is_fragment_code()` uses — `gemmi.find_tabulated_residue()`
plus the `COMMON_NON_LIGANDS` set — but **on the client**, where there is no
gemmi. Port the cheap half: exclude the twenty standard amino acids, the
standard nucleotides, water, and the additive list. That list is small,
static, and already written down in `server/ccp4i2/lib/campaign_scene.py`;
the design decision is to **duplicate it into
`client/renderer/lib/ligand-codes.ts` rather than fetch it**, because it is a
constant, and a round trip to learn that `GOL` is glycerol is absurd.

### When there is more than one candidate

Usually there is exactly one. When there are several, **show them**: the
button becomes a split button, or a small menu of codes. Do not guess, and do
not silently take the first — two fragments in one dictionary is exactly the
case where guessing produces a wrong model that looks right.

### When there are none

Disable the button with a tooltip saying why ("no ligand dictionary for this
job"). Do not fall back to prompting for a code — that is the Get monomer
dialog, which already exists and is one menu away. The button's whole claim
is that it knows the answer; a version that asks is worse than not having it.

---

## Which molecule: the one the dictionary is attached to

`addLigandOfType(code, fromMolNo)` takes `fromMolNo` as *where to look for the
dictionary*, and `this` as *what to merge into*. In our case they are the same
molecule, and that is not a coincidence — it is forced by the dictionary rule
recorded in `lib/moorhen-dictionaries.ts` and the `moorhen-scenes` skill: a
dictionary is attached to **one** molecule via `addDict`, never to Coot's
global store, precisely so that two datasets with a ligand called `DRG` do not
inherit each other's chemistry.

So: `target.addLigandOfType(code, target.molNo)`.

**Identifying `target`.** The member's coordinates are loaded by
`fetchMolecule`, which sets `newMolecule.uniqueId = url` where the url is
`/api/proxy/ccp4i2/files/{id}/download/`. Track that file id alongside
`ligandDictFileId` when `asCurrentMember` is true, and find the molecule by
matching `uniqueId`. Do **not** use "the active molecule" or
`molecules[molecules.length - 1]`: the campaign page can hold extra jobs
brought in through the project browser (`importJobFiles`), and merging a
fragment into somebody else's model is silent and not obviously wrong on
screen.

If the tracked molecule is gone (the user removed it), disable the button.

---

## Place, or place and fit?

`get_monomer_and_position_at` puts the monomer at the view centre in an
arbitrary orientation. It is *placed*, not *fitted*. The user then rotates and
real-space-refines it by hand.

There is an obvious next affordance — call Coot's ligand fitting into the
difference map after placing — and it should be a **separate, later decision**,
for two reasons:

* Fitting needs a map, a choice of which map, and a decision about what to do
  when the fit is poor. That is a dialog's worth of choices, which is the
  thing this button exists to avoid.
* `SubstituteLigand` *already* runs Coot ligand fitting as its Phase 5. If
  automated fitting is what is wanted, the job route is better: it is
  recorded, re-runnable and produces a project file. The button's niche is the
  case where the automatic fit did *not* happen or did not work, and a human
  is doing it by hand.

**Recommendation: place only, in v1.** Name the button accordingly — "Add
ligand here", not "Fit ligand here" — so it does not promise a fit.

---

## What happens to the result

The merged ligand lives in the browser. Nothing is saved.

`PushToCCP4i2Panel` (already in the campaign control panel) is the existing
route back: it creates a `coordinate_selector` job in a chosen project,
uploads the serialised molecule to `XYZIN`, and runs it. That is the right
mechanism and should not be duplicated.

The design decision is **whether the button nudges the user towards it**.
Recommendation: yes, but only as a snackbar line after a successful add —
"Ligand DRG added. Push to CCP4i2 to keep it." Do not auto-push: a placed,
unfitted ligand in an arbitrary orientation is not something anyone wants
silently written into their project.

---

## Undo

`mergeMolecules` changes the molecule in Coot. Moorhen has its own history
(`MoorhenHistory`), and `addLigandOfType` goes through `cootCommand` with
`changesMolecules: [this.molNo]` inside `mergeMolecules`, so it should land in
the session history and be undoable by Moorhen's own undo.

**This is the one claim in this document that has not been checked against a
running viewer.** Verify it before shipping; if undo does not cover it, say so
in the tooltip rather than building a bespoke undo.

---

## Implementation

All client-side. No server change, no new endpoint, no migration.

| Step | File | What |
|---|---|---|
| 1 | `client/renderer/lib/ligand-codes.ts` (new) | `NON_LIGAND_CODES` (ported from `campaign_scene.COMMON_NON_LIGANDS` + standard residues) and `candidateLigandCodes(dictTexts: string[]): string[]`, built on the existing `extractDictCompIds`. |
| 2 | `campaign-moorhen-wrapper.tsx` `fetchJobFiles` | When `asCurrentMember`, keep what is currently discarded: `setLigandCodes(candidateLigandCodes(dictionaries.map(d => d.text)))` and `setMemberCoordFileId(coordFile.id)`. |
| 3 | `campaign-moorhen-wrapper.tsx` | `handleAddLigand(code)`: find the molecule whose `uniqueId` ends in `/files/${memberCoordFileId}/download/`; `await mol.addLigandOfType(code, mol.molNo)`; `await mol.redraw()`; `dispatch(setRequestDrawScene(true))`; snackbar via `usePopcorn`. |
| 4 | `campaign-control-panel.tsx` | The button, next to the existing `Ligand2DView` — that panel already shows *which* ligand is sought, so "add it" belongs beside it. Props: `ligandCodes: string[]`, `onAddLigand: (code) => Promise<void>`, disabled with a reason when empty. |

Roughly 120 lines, most of it the code list in step 1.

### Tests

* `client/renderer/__tests__/ligand-codes.test.ts` — `candidateLigandCodes`
  against real dictionary text: an acedrg `DRG` dict yields `["DRG"]`; a
  refmac `LIBOUT` carrying `ALA`/`HOH`/`GOL` alongside `LIG` yields `["LIG"]`;
  a two-fragment dict yields both, sorted.
* A component test that the button is disabled with no codes, renders one
  click-target per code with two, and calls `addLigandOfType(code, molNo)` on
  the molecule matching the tracked file id — with a second, decoy molecule in
  the list to pin the "not the active molecule" rule.

No i2run or API test: nothing server-side changes.

---

## What this deliberately does not do

* **No dictionary import.** If the dataset has no dictionary, this button is
  not the place to acquire one; that is an acedrg job.
* **No code entry.** See above.
* **No fitting.** See above.
* **No new API.** Everything it needs is already fetched by the page.

## Open questions

1. **Undo** — unverified (above).
2. **Where the monomer lands when the user has not centred on the site.**
   `get_monomer_and_position_at` uses the view origin, so an uncentred view
   puts the fragment somewhere arbitrary. Worth considering whether the button
   should be enabled only when a site is selected, or whether the snackbar
   should just say where it went. Leaning towards: no constraint, because the
   user is looking at density when they press it, which *is* the centring.
3. **Chain and residue number of the merged ligand.** `merge_molecules` picks;
   we do not control it. If campaigns need a predictable chain (for the
   site-scene CIDs of the sibling design), that has to be checked and possibly
   forced afterwards.

---

## As built

Implemented as described, with two departures, both improvements:

* **The lookup lives in `lib/ligand-codes.ts`, not inline in the wrapper.**
  The design put "find the molecule by `uniqueId`, call `addLigandOfType`,
  redraw" inside `handleAddLigand`. Moving it into `moleculeForFile()` and
  `placeLigand()` makes the "not the active molecule" rule unit-testable
  against a decoy, which the design's own Tests section asks for and which a
  wrapper-embedded version cannot deliver.
* **The button is not inside the `Ligand2DView` block.** That block is gated
  on `ligandDictFileId`, so a button placed inside it would *vanish* in the
  no-dictionary case rather than being disabled with a reason — the opposite
  of what the design asks for. It now renders whenever a member dataset is
  selected, immediately below the 2D view.

Files: `lib/ligand-codes.ts`, `components/moorhen/add-ligand-button.tsx`
(both new), `campaign-moorhen-wrapper.tsx`, `campaign-control-panel.tsx`,
plus `__tests__/ligand-codes.test.ts` and `__tests__/add-ligand-button.test.tsx`
(12 tests). `tsc --noEmit` clean; full client suite 501 passing.

**A property worth knowing: it fails safe.** `memberCoordFileId` and
`ligandCodes` are not cleared when a dataset with no loadable job is selected,
but that cannot misfire — the action is gated on *finding the molecule*, and
`cleanupLoadedContent()` has already removed it, so the button simply
disables. Any future change that gates on the ids instead of the molecule
would lose this.

**Known limitation, by design:** a job brought in through the project browser
(`importJobFiles`, `asCurrentMember: false`) does not enable the button, even
though its coordinates are on screen. That is the wrong-molecule rule working:
adding the member's fragment to a different job's model is exactly what it
exists to prevent.

**Still unverified:** whether `addLigandOfType` lands in Moorhen's undo
history. It is a runtime property; neither `tsc` nor vitest can show it. Check
it in a running viewer before relying on it.
