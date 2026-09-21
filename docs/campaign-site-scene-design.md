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

## Frames: when to superpose

Everything above assumes the members' coordinates and the site origin share a
frame. In a fragment campaign that is usually true — the members are
molecular-replaced or rigid-body-fitted from the same reference — and the
existing summary scene already assumes it, drawing every hit on the parent
ribbon with no transformation.

But it is an assumption, and it fails on origin/indexing ambiguity in
polar and high-symmetry space groups.

The scene format has the remedy built in: a **`superpose`** block, `ssm` or
`lsq`, implemented in the resolver (`lib/moorhen-scene-resolver.ts:1063`,
`mov.SSMSuperpose(movChain, refMolNo, refChain, true)`).

**Recommendation: do not superpose by default; detect and report instead.**

* Superposing every dataset onto the exemplar costs an SSM run per dataset in
  the browser and moves atoms that, in the common case, were already right.
* The builder can *detect* the failure for free, through the diagnostic
  described above: if a dataset has a `hit` verdict and fragment-like residues
  but its **nearest one is far** from the site origin — tens of Ångströms,
  well beyond any pocket — the frames are a live suspect. Report the distance
  in `stats` and let the panel say so.
* Offer `superpose=1` as a query parameter, which emits one `ssm` entry per
  member onto the exemplar. The chain to superpose on is the exemplar's first
  polymer chain; that is a guess, which is another reason not to make it the
  default.

If real campaigns turn out to need it routinely, flip the default then, on
evidence.

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
GET /api/projectgroups/{id}/sites/{site_id}/scene/[?include=unclear][&superpose=1]
  -> {"scene": <MoorhenScene>, "stats": {...}}
```

Same envelope as `summary_scene` (`ProjectGroupViewSet.py:613`), and the same
DRF routing style as the existing site routes
(`url_path=r"sites/(?P<site_id>[0-9]+)/..."`), so `site_evaluation` next door
is the template for the 404 handling.

`stats` should carry: `hits_drawn`, `hits_claimed` (verdicts found),
`empty_verdicts`, `unclear_verdicts`, `parent_present`, whether the exemplar
was the parent or a fallback, and `skipped: [{project, reason, nearest}]` —
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

1. **`ENVIRONMENT_RADIUS = 8 Å`** is a judgement with no data behind it —
   but a cosmetic one now, which is the point. Check it against a real
   campaign; if the pocket looks sparse or bloated, change the number and
   nothing else moves.
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
