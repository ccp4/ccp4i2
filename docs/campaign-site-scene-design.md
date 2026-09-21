# A site view for a fragment campaign, composed as a Moorhen scene

**Status:** design, not implemented. Written 2026-09-21 against `django` at a71.
**Scope:** one new server-side scene builder and one endpoint. Independent of
[campaign-place-ligand-design.md](campaign-place-ligand-design.md).

---

## What is wanted

Open one binding site of a campaign and see **what has been found there**: the
protein (or DNA/RNA) drawn once as a ribbon from an exemplar, and the ligand
of *every* dataset judged a hit at that site drawn on top of it as sticks,
each with its own restraint dictionary so that two fragments sharing a
residue code keep their own chemistry.

The campaign already has a whole-campaign version of this
(`?summary=1`). A site view is the same picture, narrowed to one place.

---

## Compose it as a scene

This is the right shape, and it is already the grain of the code. The
campaign summary is **not** a bespoke render path: `build_summary_scene()`
returns a plain scene dict, `campaign-page-client.tsx` hands it to
`CampaignMoorhenWrapper` as `summaryScene`, the wrapper runs it through
`serialiseScene()` into YAML (`campaign-moorhen-wrapper.tsx:1136`), and seeds
it into the **Scenes panel** as `initialSceneYaml` with
`autoApplyInitialScene`. The panel's own parse/apply path — the single one
shared with hand-edited scenes — does the rest.

So a site view built as a scene gets, for free:

* **the same rendering pathway**, already exercised and tested;
* **an editable artefact**: the generated YAML lands in the Scenes tab, so a
  user who wants the ribbon in a different colour, one hit hidden, or a map
  contoured, edits it in place rather than asking for a feature;
* **something saveable and shareable** — `.scene.yaml` / `.scene.zip` — which
  is the natural form for "the figure of site 3 in this campaign";
* **validation**: the scene validates against the published contract, so a
  builder that emits nonsense is caught at the boundary rather than producing
  an empty viewer.

Before writing the builder, read the grammar document
(`client/renderer/types/moorhen-scene.md`) and the `moorhen-scenes` skill; do
not reconstruct the grammar from this document, which paraphrases it.

---

## What exists, and the three gaps

`server/ccp4i2/lib/campaign_scene.py` already does most of the hard part.
What it does **not** do:

**Gap 1 — it is campaign-wide, with no notion of a site.** Every hit in the
campaign goes into one scene. At a real campaign's 30–40 sites that is an
unreadable pile.

**Gap 2 — `_first_site_view()` was dead code. FIXED 2026-09-21.** It read
`getattr(group, "sites", None)` and treated the result as a list of dicts. But
migration 0024 moved sites out of the `ProjectGroup` JSON field into the
`CampaignSite` table, related name **`site_set`**. `ProjectGroup` has no
`sites` attribute any more (see the comment at `db/models.py:79`), so the
`getattr` always returned `None` and **the summary scene carried no camera
from #563 until this fix**.

It now reads `group.site_set.first()` — a real attribute rather than a
defaulted `getattr`, so the next such move breaks loudly instead of quietly
returning "no sites". Two tests in `tests/api/unit/test_summary_scene_api.py`
pin it; the first was checked to fail (`KeyError: 'view'`) against the old
body before the fix went in.

**Gap 3 — hits are detected, never *located*.** `detect_ligands()` returns
residue *codes*, and the scene draws `//*/(CODE)` — every residue of that name
in that dataset, wherever it is. So the builder cannot say which copy of `DRG`
is at *this* site, and a dataset that is a hit at two sites contributes both
ligands to both site views.

**This design deliberately does not close gap 3** — see *Draw every copy* below
for why that is the right call and not a shortcut. Gap 2 is already fixed;
gap 1 is the work.

---

## Which datasets belong to a site

Two sources of truth now exist, and they answer different questions.

**`SiteEvaluation` — the human verdict.** `hit` / `empty` / `unclear` for one
dataset at one site, with the absence of a row meaning nobody has looked
(`db/models.py:171`). This is *exactly* the question the site view asks, and
the table was built to answer it.

**`detect_ligands()` — the automatic check.** Reads the refined coordinates
for a fragment-like residue. It knows nothing about sites.

**Recommendation: the verdict decides membership, and nothing else does.**
Geometry has no say in which datasets appear — it only chooses the pocket
residues drawn for context (see *The site environment*).

* A dataset is in the site scene iff it has a `hit` evaluation at that site.
* `unclear` verdicts are included **only on request** — a `include=unclear`
  query parameter — and drawn in a muted colour. A site view is a claim about
  what is there; mixing the confident and the doubtful without marking them
  is how an overlay figure becomes misleading.
* `empty` and no-row are excluded. Excluding them is not the same judgement:
  `empty` means somebody looked, and the count of those is worth reporting in
  `stats` even though nothing is drawn.

Deliberately **not** recommended: falling back to `detect_ligands()` when a
site has no verdicts. It looks helpful and it silently reintroduces gap 3 —
a scene that shows every ligand in the campaign under the name of one site.
When a site has no hits, return an honest empty scene (exemplar ribbon, no
sticks) and say so in `stats`. "Nobody has evaluated this site yet" is
information; a wrong picture is not.

---

## Draw every copy, and let the camera do the filtering

The tempting move is to close gap 3: find the ligand nearest the site origin
and draw only that residue, by chain and sequence number. Resist it, at least
in v1.

**Draw every residue of the hit's code — `//*/(CODE)`, exactly what the
summary scene already emits — and let the camera decide what is visible.**

The reasoning:

* **Distant copies do not intrude.** A site view is centred on the site origin
  at a zoom that frames a pocket, perhaps 30–40 Å across. A second copy 60 Å
  away is off-screen. `view.slab` clips it in depth as well — note that slab
  takes a *selection*, not a point, so it hangs off the pocket residues the
  builder computes below rather than off the origin directly.
* **Proximal copies are the interesting case, not the failure case.** A
  fragment bound in an adjacent subsite is precisely what a site view should
  show: it speaks to compatible or cooperative binding, and it is the
  observation that fragment linking and growing start from. A radius filter
  would throw that away, and throw it away *silently* — the user would never
  learn that the neighbouring density existed.
* **It embeds no judgement that can be wrong.** A 10 Å cutoff is a guess about
  what "at this site" means, applied invisibly, and a wrong one produces a
  picture that looks authoritative and is missing a ligand. "Everything this
  dataset has, framed on the site" is a statement with no hidden claim in it.
* **It is already implemented.** `detect_ligands()` and the `//*/(CODE)`
  selection exist and are tested. No new gemmi function, no radius constant in
  the selection path, no boundary cases.

### Keep the distance measurement anyway — as a diagnostic

Dropping the radius from the *selection* does not mean not computing it. The
builder reads each hit's coordinates with gemmi regardless (that is what
`detect_ligands` does), so the centroid of each fragment-like residue and its
distance to the site origin are essentially free.

Report the nearest such distance per dataset in `stats`. It answers a question
the picture cannot: a dataset with a `hit` verdict whose nearest fragment is
40 Å from the site origin is either a verdict recorded against density nobody
has modelled yet, or two structures in different frames (see below). That is
worth surfacing in the panel; it is not worth using to hide a ligand.

### The upgrade path, if clutter turns out to be real

If a real campaign shows that all-copies is too noisy — many sites, many
multi-copy datasets — the fix is a `detect_ligand_residues()` returning
`[(chain, seqid, code, centroid)]`, with `detect_ligands()` reimplemented on
top of it so the existing hit detection and its tests are untouched, and a
per-residue `//<chain>/<seqid>` selection. Build it then, on evidence of
clutter, not now on the anticipation of it.

## The site environment: residues around the origin

The complement to drawing the ligands is drawing what they bind *to*. A ribbon
alone does not show which residues line the pocket, and that is most of what
makes a site view readable.

**Draw the exemplar's residues within `ENVIRONMENT_RADIUS` of the site origin
as sticks, over the ribbon.** The exemplar is the right structure to take this
from: it is the campaign's reference frame, it is drawn once, and taking the
pocket from each hit in turn would overlay six near-identical copies of the
same side chains.

Mechanically: `gemmi.NeighborSearch` over the exemplar for atoms within the
radius of `(origin_x, origin_y, origin_z)`, collapse to whole residues, and
emit an explicit CID list — `//A/45||//A/47||//B/112` — as a second
representation on the exemplar element. Explicit, because a Coot CID cannot
express a sphere; this is exactly the kind of computation a server-side scene
builder exists to do, and it leaves a scene that says in the YAML which
residues it decided the pocket was, where a user can edit the list.

**This is where the arbitrary radius belongs**, and it is worth noticing what
the previous section did to it. A radius that decides *which ligand is a hit
at this site* is a correctness judgement, applied invisibly, and damaging when
wrong. A radius that decides *which side chains to draw as sticks* is cosmetic
— wrong only in the sense of showing a residue too many or too few, and
visibly so. Same constant, much safer job. Start at **8 Å**.

### `residue_environment` is not this

Moorhen has `environment`, `ligand_environment` and `residue_environment`
representation styles, and the scene grammar lists them. They are not a
substitute here. `getEnvironmentBuffers` passes its CID through `cidToSpec` to
a **single** residue and calls `make_exportable_environment_bond_box` against
`parentMolecule` — so it draws contacts and H-bonds around one residue, within
one molecule. Applied to a hit's ligand it would find that *hit dataset's* own
pocket residues, not the exemplar's.

That is a genuinely useful second affordance — "show me the contacts for this
hit" — and a good candidate for a later per-hit toggle. It is not the
site-wide pocket, and reaching for it expecting one would produce an overlay
of six datasets' side chains.

## Frames: superpose locally, on the site

**This section previously recommended not superposing by default. That was
wrong, and a picture settled it.** A screenshot of the campaign summary
(2026-09-21) shows the members' helices drawn as a fan of offset ribbons,
displaced from one another by on the order of an Ångström, with the fragments
spread correspondingly. Binding events that are probably equivalent are made
to look different, and the difference is the protein frame, not the chemistry.

The earlier reasoning was that campaign members are molecular-replaced or
rigid-body-fitted from a common reference and are therefore already in frame.
That holds for datasets processed *through* the campaign. It does not hold for
datasets imported from the PDB, each deposited with its own origin choice
within the space group — which is what a demo campaign is made of, and what a
real one becomes as soon as anyone brings in a published structure. And the
screenshot suggests it is not reliable even when it should be.

So: **superpose by default.** The cost argument was real but small — an LSQ fit
is cheap next to loading the structure it applies to — and it was being
weighed against a picture that is actively misleading, which is the wrong
trade.

### Local, not global

The important word is *local*. A global SSM onto the exemplar removes gross
frame differences, but it distributes the residual over the whole molecule,
and none of that residual is guaranteed to land anywhere but the pocket. For a
view whose entire purpose is comparing ligands at one site, the frame that
must coincide is **the site's**.

The scene format already expresses this. `superpose` takes an `lsq` method
with explicit residue ranges, and the resolver implements it
(`moorhen-scene-resolver.ts:1084`, `mov.lsqkbSuperpose(...)` → Coot's
`add_lsq_superpose_match` + `lsq_superpose`).

### What to fit on: CAs, in a sphere that grows until it has enough

Three quantities, deliberately distinct. An earlier draft of this document
claimed the pocket residues drawn as sticks could double as the fit target —
"one computation, two uses". That was too neat. Seeing and fitting want
different atoms at different distances:

| | Criterion | Radius | Job |
|---|---|---|---|
| Environment sticks | **any atom** within radius | `ENVIRONMENT_RADIUS`, 8 Å | what the user *sees* lining the pocket |
| LSQ fit target | **CA only** within radius | `FIT_RADIUS`, from 15 Å, grown | what defines the local *frame* |

**CA only, not any atom.** The fit is on main chain, so a residue earns its
place by where its backbone is. A long side chain — Arg, Lys, Glu — can reach
several Ångströms into a pocket from a backbone that is nowhere near it;
including it on the strength of that reach adds a point that does not belong
to the local frame. Conversely a residue whose CA is close but whose side
chain points away is exactly what should anchor the fit.

**From 15 Å, and grow if need be.** 15 Å around a site in a folded domain
typically encloses several dozen residues, far more than a stable fit needs.
But a shallow surface site, a small protein, or a site near a domain edge can
be sparse, so the radius grows in steps until the CA count is sufficient.

**Cap the growth.** Past roughly 25–30 Å it is no longer a *local*
superposition in any meaningful sense — it is a global fit wearing a sphere.
At that point stop, fall back to global `ssm`, and record it in `stats`.
Silently growing to 40 Å would produce the very thing this section exists to
avoid, while reporting success.

### Tolerating disorder: fit on the intersection, and re-count after

Residues present in the exemplar are routinely missing from a member — a
disordered loop at the pocket rim is ordinary, not exceptional. An earlier
draft said to check the exemplar's pocket residues exist in the mover and skip
the dataset's superposition where they do not. That is too brittle: one absent
residue would forfeit the whole fit.

The rule instead:

1. Take the exemplar's CAs within the current `FIT_RADIUS` of the site origin.
2. **Intersect** with the residues that are present *and have a CA* in the
   moving structure.
3. **Apply the count gate to the intersection, not to the exemplar's list.** A
   pocket that is ample in the exemplar can fall below threshold in a dataset
   with a disordered rim, and that dataset's fit is the one that would be bad.
   If it is short, grow the radius and repeat; if the cap is reached, fall
   back to `ssm`.
4. Build the match ranges **from the intersection**, collapsed into contiguous
   runs per chain.

Step 4 matters more than it looks. The format's `matches` are *ranges*, and
whether Coot pairs by residue number and skips residues absent from one side,
or does something less forgiving, is **not verified here** — it lives in
`lsq_superpose`'s handling of `add_lsq_superpose_match`. Building ranges from
the intersection means the question never arises: every range spans only
residues that exist on both sides. That is worth a little extra fragmentation
(`1858-1861`, `1863-1867` rather than `1858-1867`) to avoid depending on
behaviour nobody has checked.

**Each dataset therefore gets its own match list**, because each has its own
intersection with the exemplar. This is not one shared block of ranges
repeated per `move:` — a sketch that shows identical ranges for every dataset
is showing the easy case, not the general one.

`matchType: main` (the format's default) is right here. Side chains move
between datasets — sometimes *because* of the fragment — so fitting on them
would let the ligand's own effect pull the frame around.

### Compute the fit ourselves, and hand Coot a matrix

Everything above describes *which atoms* to fit on. It does not require that
Coot do the fitting, and there is a good case for doing it server-side instead.
Both halves have been checked against the real APIs:

* **gemmi can do the fit.** `gemmi.superpose_positions(pos1, pos2, weight=[])`
  returns a `SupResult` with `.transform` (`.mat`, `.vec`), `.rmsd` and
  `.count`. Verified on gemmi 0.7.5: fitting four points against the same four
  translated by +5 Å in x gives `rmsd 0.0` and `vec (-5, 0, 0)` — so the
  transform maps **pos2 onto pos1**. Call it as
  `superpose_positions(reference_CAs, moving_CAs)` and the result moves the
  member onto the exemplar. Worth stating explicitly: that argument order is
  an easy silent sign error.
* **Coot can apply an arbitrary transform.**
  `apply_transformation_to_atom_selection(imol, cid, n_atoms, m00…m22, c0 c1 c2,
  t0 t1 t2)` — `molecules-container.hh:2747`, exposed in the wasm bindings.
  Note the **rotation centre** `c0 c1 c2` is separate from the translation.
  gemmi's transform is `x' = mat·x + vec` about the origin, so feed
  `c = (0, 0, 0)` and `t = vec`. Passing a centroid as the centre and `vec` as
  the translation would apply the shift twice.

**Recommendation: compute the transform in Python, carry it in the scene, and
have the resolver apply it.** The reasons are about what can be tested and
what can be read:

* **It is testable without a browser.** Radius growth, the intersection, the
  CA count gate, the fallback — all of it becomes a pure function over two
  structures and a point, with unit tests in `tests/unit/lib/`. Routed through
  Coot, none of that logic can be tested at all; it can only be looked at.
* **It removes the range question entirely.** No `matches`, no contiguous
  runs, no dependence on how `add_lsq_superpose_match` treats a residue absent
  from one side.
* **It makes the fit inspectable.** The scene carries the matrix *and* what it
  was derived from — how many CAs, at what radius, to what RMSD. A reader can
  judge the superposition instead of trusting it, and a lifted scene
  reproduces it exactly rather than re-deriving something slightly different.
* **It allows outlier rejection**, which is what actually makes a local fit
  robust. Fit, drop CAs beyond ~2× RMSD (a rim loop that genuinely moved),
  refit. gemmi's `weight` argument supports the softer version. Coot's LSQ
  will not do this for us, and without it one shifted loop drags the whole
  frame.

The shape, as a third `method` alongside `ssm` and `lsq`:

```yaml
superpose:
  - method: matrix
    move: x0104
    mat: [1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0]   # row-major
    vec: [-0.31, 0.12, 0.05]
    fitted:            # provenance: what the matrix was derived from
      onto: reference
      atoms: 47        # CAs in the intersection
      radius: 15       # Å, after any growth
      rmsd: 0.21
```

**The honest cost is a schema change.** The scene format is generated and
CI-gated, so this means editing the Zod source and regenerating the
contracts — read `MOORHEN_SCENES_SCHEMA_V1_DESIGN.md` first, as the
`moorhen-scenes` skill directs, and do not hand-edit the generated files.

**If that cost is unwelcome, `lsq` with ranges is a legitimate first cut** —
it needs no schema change and the resolver already implements it. Take it
knowing what is given up: the selection logic stays untestable, the range
semantics stay unverified, and there is no outlier rejection. It is the
cheaper half of the same idea, not a different one.

### Constants, and what breaks them

* `FIT_RADIUS_START = 15 Å`, grown in 5 Å steps to `FIT_RADIUS_MAX = 30 Å`.
* `MIN_FIT_CAS`: enough for a conditioned rigid-body fit with margin. Three
  non-collinear points determine the transform in principle; that is not a
  usable floor. Start at **12**, and note it is a floor, not a target — the
  15 Å sphere will normally supply several times that.
* **Numbering must correspond** between exemplar and member. True for a
  campaign on one construct; false as soon as a differently-numbered PDB entry
  arrives. The intersection makes this fail *safely* rather than wrongly — a
  mismatched numbering simply yields a small intersection, which trips the
  count gate and falls back to `ssm`. That is the right failure: a silent fit
  on coincidentally-numbered residues would be far worse than a global one.
* **The superposition is site-specific.** Two site scenes of the same campaign
  align the same datasets differently, because each aligns on its own pocket.
  That is correct, and worth stating plainly: a site scene is a view *from* a
  site, not a general overlay that happens to be centred on one.

### The campaign summary has the same disease

The screenshot is of the whole-campaign summary, not a site view, so this is
not only a site-scene concern. The summary has no single pocket to fit on, so
the local argument does not apply — but a **global `ssm` onto the exemplar**
is still strictly better than the staggered fan it currently draws, and is a
small change to `build_summary_scene`. Worth doing independently of the site
view.

### Keep the diagnostic

The nearest-fragment distance is still worth reporting, and superposing makes
it *more* meaningful rather than less: measured after alignment, a hit whose
fragment is still far from the site origin is a real anomaly rather than a
frame artefact.

---

## The exemplar

**The parent project**, via `ProjectGroup.parent_project` — a campaign has at
most one, enforced by a constraint, and it *is* the reference frame the site
origin is expressed in (`db/models.py:85-97`). Drawing anything else as the
ribbon while positioning by the parent's frame would be quietly inconsistent.

`build_summary_scene` already uses it, through `_parent_coord_file`.

When there is no parent, fall back to the hit dataset with the lowest project
id (stable, arbitrary, and stated as arbitrary in `stats`) and record that the
frame is that dataset's. Do not silently produce a sticks-only scene: a
ribbon-less overlay of six fragments in space is not a site view.

---

## Sketch of the output

Illustrative only — the grammar document is authoritative.

```yaml
scene: BAZ2B campaign — site "Acetyl-lysine pocket"
version: 1
authoredIn: { projectName: BAZ2B fragments }
files:
  - { name: reference, kind: coordinates, fileId: 4012, projectId: <uuid> }
  - { name: x0104, kind: coordinates, fileId: 5120, projectId: <uuid> }
  - { name: x0104_dict, kind: dictionary, fileId: 5118, projectId: <uuid> }
  - { name: x0212, kind: coordinates, fileId: 5301, projectId: <uuid> }
  - { name: x0212_dict, kind: dictionary, fileId: 5299, projectId: <uuid> }
superpose:
  # Fitted server-side on the CAs near the site origin, present in BOTH
  # structures. Each entry is its own fit: x0212 is missing a disordered rim
  # loop, so fewer CAs and a slightly worse RMSD -- visible here rather than
  # buried.
  - method: matrix
    move: x0104
    mat: [0.9999, -0.0121, 0.0043,  0.0121, 0.9999, -0.0018,  -0.0043, 0.0018, 1.0000]
    vec: [-0.31, 0.12, 0.05]
    fitted: { onto: reference, atoms: 47, radius: 15, rmsd: 0.21 }
  - method: matrix
    move: x0212
    mat: [0.9998, 0.0184, -0.0072,  -0.0184, 0.9998, 0.0031,  0.0072, -0.0031, 1.0000]
    vec: [0.44, -0.19, -0.08]
    fitted: { onto: reference, atoms: 31, radius: 15, rmsd: 0.34 }
elements:
  - file: reference
    representations:
      # the whole exemplar, once
      - { style: CRs, selection: /*/*/*/*, colour: "#b0bec5" }
      # the pocket: residues within 8 A of the site origin, enumerated by the
      # builder because a CID cannot express a sphere
      - style: CBs
        selection: //A/1863||//A/1867||//A/1870||//A/1886||//A/1893
        colour: "#90a4ae"
  - file: x0104
    dictionaries: [x0104_dict]
    representations:
      # every copy of this dataset's fragment; the camera frames the site
      - { style: CBs, selection: //*/(DRG), colour: "#1f77b4" }
  - file: x0212
    dictionaries: [x0212_dict]
    representations:
      - { style: CBs, selection: //*/(DRG), colour: "#d62728" }
view:
  origin: [12.4, -3.1, 28.7]
  quat: [0.0, 0.0, 0.0, 1.0]
  zoom: 0.35
  slab: { file: reference, selection: "//A/1863||//A/1893", pad: 8 }
resolver: { onMissingResidues: clamp-and-log }
```

Notes on the choices visible there:

* **Each hit is fitted on the pocket before anything is drawn** — without it
  the datasets fan out by an Ångström or so and equivalent binding events read
  as different ones. See *Frames*.
* **Each hit contributes its ligand's every copy**, by code, not one located
  residue — see *Draw every copy*. The `view` is what makes that legible:
  `origin` and `zoom` frame the site in x/y, and `slab` clips it in depth.
* **The exemplar carries two representations**: the ribbon for context and the
  pocket residues as sticks. Only the exemplar's protein is drawn; the hits
  contribute ligands alone, or six copies of the same side chains would pile up.
* **Each hit gets its own colour** from a fixed categorical palette, so the
  hits are tellable apart against the grey ribbon. `PARENT_RIBBON_COLOUR`
  (`#b0bec5`) already exists and was chosen for this. Reuse it, and give the
  pocket sticks a near neighbour of it so the context reads as context.
* **`dictionaries` is per element**, never `globalDictionaries` — the rule
  recorded in `lib/moorhen-dictionaries.ts` and the `moorhen-scenes` skill. A
  global entry is inherited by every molecule that lacks its own, which is
  precisely the collision a site view of six `DRG`s would hit.
* **`view` comes from the site row**: `origin` from `origin_x/y/z`, `quat` and
  `zoom` from the saved camera when present. This is what `_first_site_view`
  was meant to do and no longer does. `slab` is the builder's own, derived
  from the pocket selection it just computed — `view.slab` sets depth only and
  `view.origin` sets the camera, so both are needed (see the grammar).
* **No maps.** See below.

---

## Maps

**Recommendation: none in v1.** N datasets' MTZs is a lot of bytes for a
picture whose point is the *positions* of the fragments, and a pile of
overlapping difference densities is unreadable.

The natural follow-up, once the view exists, is a map for **one** hit at a
time — which the Scenes tab can already express by hand (`maps:` plus
`columns:`), so the feature can be discovered before it is built. Leave the
question open rather than guessing at it now.

---

## API

```
GET /api/projectgroups/{id}/sites/{site_id}/scene/[?include=unclear][&superpose=none]
  -> {"scene": <MoorhenScene>, "stats": {...}}
```

Same envelope as `summary_scene` (`ProjectGroupViewSet.py:613`), and the same
DRF routing style as the existing site routes
(`url_path=r"sites/(?P<site_id>[0-9]+)/..."`), so `site_evaluation` next door
is the template for the 404 handling.

Superposition is on by default, so the parameter turns it **off**
(`superpose=none`) rather than on — for the rare case of wanting to see the
frames as deposited, and for diagnosing a fit that went wrong.

`stats` should carry: `hits_drawn`, `hits_claimed` (verdicts found),
`empty_verdicts`, `unclear_verdicts`, `parent_present`, whether the exemplar
was the parent or a fallback, and per dataset how it was
superposed — the CA count, the radius it took to reach it, the RMSD, and
whether it fell back to `ssm` because the count gate could not be met inside
the radius cap. A fit is only as trustworthy as those numbers, so they belong
in the payload rather than the server log. Also
`skipped: [{project, reason, nearest}]` —
where `nearest` is the distance from the site origin to the closest
fragment-like residue, because that number is what tells a user whether the
verdict is premature or the frames disagree.

### Client wiring

The campaign page already reads a `site` URL parameter
(`campaign-page-client.tsx`, `siteParam`) and already passes it down as
`initialSiteId`, and the verdict chips in the overview already navigate to
`...?job=<jobId>&site=<siteId>`. The change is narrow: when `summary=1` **and**
`site=<id>` are both present, fetch the site scene instead of the campaign
scene and hand it to the same `summaryScene` prop. Everything downstream —
serialise, seed the Scenes panel, auto-apply — is unchanged.

Add a way to reach it: a "View this site" action on each row of the site list
in `campaign-control-panel.tsx`, next to the existing go-to-site button.
Going *to* a site (move the camera in the current dataset) and *viewing* a
site (load every hit) are different acts and should be different controls.

---

## Refactor, don't fork

`build_summary_scene` and `build_site_scene` draw a hit dataset identically —
same job resolution, same coordinate file, same dictionary, same
`//*/(CODE)` sticks. They differ only in which datasets they include, and in
what the exemplar and the camera do. Extract the shared middle — resolving a
member's refinement job, its coordinate file, its dictionary file or inlined
CIF text, and the `files[]`/`elements[]` pair that results — into one helper,
and let each builder decide membership, the exemplar's representations and
the view.

Two near-identical 120-line builders in one module is how the campaign-scene
hit rule ends up fixed in one and not the other.

---

## Performance

The builder reads every candidate dataset's coordinates with gemmi, in a
synchronous request. `build_summary_scene` already does this for the whole
campaign, so a site view is strictly cheaper — but "cheaper than something
that is already slow" is not a defence.

Note it, measure it on the 40-dataset demo campaign
(`manage.py make_demo_campaign`), and if it is bad, cache per
`(job_id, mtime)` rather than making the endpoint async. The module docstring
already declares itself "deliberately Django-light ... callable from a
synchronous ViewSet today and a background worker tomorrow"; honour that and
keep the builder free of request objects.

---

## Tests

* `tests/unit/lib/test_campaign_scene.py` (extend) — the pocket selection, as
  a pure function over an exemplar path and an origin: residues within the
  radius are returned as CIDs, one entry per residue however many of its atoms
  are in range, sorted and deduplicated, and an origin in empty solvent
  returns an empty list (which must degrade to *no* pocket representation, not
  an empty `selection:` that draws everything).
* The nearest-fragment diagnostic, likewise pure: a structure and an origin in
  gives a distance out, and a structure with no fragment-like residue gives
  `None` rather than `inf`.
* **The fit itself**, as a pure function over two structures and a point —
  which is the whole reason for computing it server-side. Take a structure,
  apply a known rotation and translation, and assert the recovered transform
  inverts it to within floating-point noise and reports ~0 RMSD. A sign error
  in the `superpose_positions` argument order passes every vaguer test and
  fails this one.
* **Disorder tolerance**: delete a rim loop from the moving copy and assert
  the fit still succeeds, on the intersection, with a correspondingly lower
  `atoms` count — not that it is skipped.
* **The count gate and radius growth**: a site with too few CAs at 15 Å grows
  the radius; one that cannot reach `MIN_FIT_CAS` by the cap falls back to
  `ssm` and says so in `stats`, rather than emitting an ill-conditioned fit.
* **Numbering mismatch fails safely**: renumber the moving copy and assert the
  intersection collapses, the gate trips, and the result is an `ssm` fallback
  — never a confident fit on coincidentally-numbered residues.
* **Outlier rejection**: displace a few CAs far from their partners and assert
  they are dropped and the RMSD of the retained set is small, rather than the
  whole frame being dragged.
* `tests/api/unit/test_site_scene_api.py` — mirroring
  `test_summary_scene_api.py`: a site with two hits yields two ligand
  elements and one ribbon; `empty` and unevaluated members are absent;
  `include=unclear` adds them and omitting it does not; a site with no hits
  returns an exemplar-only scene (ribbon and pocket, no sticks) and
  `hits_drawn: 0`; a site of another campaign 404s. Per the memory note, **do not** `@django_db` these — `tests/api/`
  has its own DB fixtures.
* The `_first_site_view` tests are already in
  `tests/api/unit/test_summary_scene_api.py` (gap 2, fixed). A site scene's
  own camera should get the same treatment: the site addressed in the URL is
  the one in `view`, not the campaign's first.

---

## Open questions

1. **`ENVIRONMENT_RADIUS = 8 Å` is no longer purely cosmetic**, now that the
   same residue list drives the LSQ fit. Too small gives an ill-conditioned
   superposition; too large drags in residues that move for reasons unrelated
   to the site. It is still a far safer knob than a ligand-membership radius
   would have been, but it now wants checking against a real campaign on both
   counts, not just on how the sticks look.
2. **Whether all-copies really is quiet enough.** The argument above is that
   distant copies fall outside the frame. It has not been checked on a
   campaign with genuinely multi-site hits at a zoom a user would choose. If
   it turns out noisy, the per-residue upgrade path is written down; take it
   then.
3. **A ligand bridging two adjacent sites** appears in both site views. That
   is almost certainly right — it *is* at both — but it means the sum of the
   site views double-counts, so any "hits per site" figure must come from the
   verdicts, never from what the scenes drew.
4. **Whether the site view should offer the verdict control** for the datasets
   it shows, so that looking at a site and recording what is there is one act.
   Attractive, and out of scope for the first cut.
