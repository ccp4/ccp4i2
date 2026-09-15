# `molrep_map`: fast cryo-EM map → model placement, prepared for refinement

**A design note and port brief.** Place a model into a cryo-EM map, resolving the
map's hand ambiguity, and emit a **refinement-ready package** — a trimmed,
valid-cell map, a mask, and the placed model, all in one frame — that drops
straight into both *interactive* (Coot/Moorhen, real space) and *(pseudo-)
reciprocal-space* (servalcat/refmac) refinement. This is the brief for porting
the Qt-era `molrep_map` task into the Django branch; it is complete enough to
implement from without further supervision, and it records the crystallographic
decisions that are easy to get subtly wrong.

## 0. Provenance

Ported from **`martinemnoble1/ccp4i2`, branch `chapi`,
`wrappers/molrep_map/`** (a `CPluginScript` that used `chapi.flip_hand` to invert
the map and ran `molrep` on both hands). The Django port keeps the two-hand logic
and the map-preprocessing intent but:

- does the **map surgery in gemmi**, not Coot's `chapi` (see §4 — `chapi`/
  `coot_headless_api` is not in the slim server, and the surgery this task needs
  is gemmi header/cell/grid work, not something `chapi` helps with);
- makes the **MR engine pluggable** — `molrep` by default, `phaser` as an
  alternative (§6);
- emits a **refinement-ready `(map, mask, model, Δ)`** package rather than just a
  placed model (§3, §7).

Fetch the four source files (`molrep_map.py`, `molrep_map.def.xml`,
`molrep_map_gui.py`, `molrep_map_report.py`) from that branch when implementing —
they are the reference for the molrep command construction and the report shape.

## 1. Purpose and the speed-first stance

The task exists to **accelerate map → model evaluation** — a triage step, not a
final refinement. Speed is a first-class requirement, and it drives two defaults:

- **Aggressive down-sampling** of the *placement* map: coarsen the grid as far as
  still gives MR enough signal to find the orientation. cryo-EM placement needs
  only the low-frequency envelope; a 2–4× down-sample (or a resolution cap of
  ~6–8 Å for the search) is typically ample and cuts the FFT cost by an order of
  magnitude.
- **A wall-clock time limit** on the MR engine, exposed as a control parameter,
  defaulting low. A run that has not placed the model within the budget returns a
  clear "no placement in the time allowed" rather than grinding.

Both are *controls with fast defaults*, not hard-coded — a user who wants a more
thorough search can loosen them.

Down-sampling is **only ever applied to the disposable placement map**. The
deliverable map is full-resolution (§2).

## 2. The deliverable is one shared-frame `(maps, mask, model)` package, two consumers

The output is not "a placed model plus a map to look at." It is a
**refinement-ready package in a single coordinate frame** that must satisfy two
consumers at once:

| Consumer | Space | Needs |
|---|---|---|
| **Coot / Moorhen** | real (interactive) | full-res, *trimmed* (renderable) map, registered with the model |
| **servalcat** | real map in → SFs internally | the *same* frame's map(s) + model; **prefers two half-maps** (FSC cross-validation), does its own trim/mask/FFT |

Both consumers take a **real-space map** in the **model's frame** — not structure
factors. That is the key simplification the port makes over the old draft: **no
map-coefficient conversion, no spoofed `F`/`σ`.** servalcat's `refine_spa` takes
`--map`/`--halfmaps` directly and FFTs internally with its own weighting (§0's
help dump confirms it), and Moorhen contours the map on the GPU. So the emitted
map must be, simultaneously:

- **full resolution** (down-sampling is placement-only, never the deliverable);
- **trimmed** to the molecule + a **solvent border** — tight enough that Moorhen
  can contour it (the failure mode is *extent*: a 400³–512³ box is millions of
  voxels and an isosurface to match, and WebGL falls over), loose enough that a
  clean FFT is still possible and the mask has room to fall off. (servalcat will
  *also* trim internally, so this trim is driven by Moorhen's render budget; a box
  that satisfies Moorhen satisfies servalcat.);
- a **valid P1 cell**: `CELLA = grid × voxel`, angles 90°, `MX/MY/MZ` = the new
  counts, **voxel unchanged**, on **composite grid dimensions** (products of small
  primes) so any FFT is efficient and exact;
- **axis-order normalised** (`gemmi.Ccp4.setup`) to a standard axis order;
- **grid-origin-consistent with the model — via `nxstart`, never the Å `ORIGIN`**
  (§4, confirmed from Moorhen/coot source).

### 2a. Half-maps: optional in, for cross-validation out

servalcat *prefers two half-maps* and only grudgingly accepts a single `--map`
("Use this only if you really do not have half maps") — because the two
independent half-maps give it **FSC-based cross-validation**: the resolution-
dependent weighting and the guard against over-fitting both come from the
half-map FSC. A single map has no such error model.

So the task **optionally accepts two half-maps** (`HALFMAP1`, `HALFMAP2`) beside
the required full `MAPIN`:

- **Placement** always runs on the **full map** (best signal for MR).
- The chosen hand's **flip + trim + frame transform is applied *identically* to
  the full map and to both half-maps** — one transform, three maps, all landing
  in the same box/frame. (If no half-maps were given, `MAPIN` can be taken *as*
  the full map and no half-maps are emitted.)
- **Emit**: the trimmed full map (`MAPOUT`, for molrep-scoring provenance + Moorhen
  display), the two trimmed half-maps (`HALFMAPOUT1/2`, for servalcat's preferred
  `--halfmaps` cross-validated path) when supplied, the mask, and the model.

When half-maps are absent the package still refines — servalcat just falls back to
single-`--map` mode. Emitting them when present is what makes the downstream
refinement *properly* cross-validated rather than merely functional.

Because trimming for rendering and trimming for FFT are the *same* origin-aware
crop, getting it wrong mis-registers the *displayed* map against the *placed*
model — the most visible way this can break.

**A mask is part of the package, not an afterthought.** servalcat's single-
particle likelihood is masked (it refines against the molecule, not the empty
box). The same bounding-box work that trims the map defines that mask, so
`preprocess_map()` emits it alongside.

## 3. Pipeline: place → choose hand → prepare → emit

```
inputs:  XYZIN (model)   MAPIN (real map)   [HALFMAP1 HALFMAP2 (optional, §2a)]

place :  for hand in {original, inverted}:            # inversion = §4 flip
             m = downsample(trim_loose(MAPIN), factor)  # disposable, fast
             placed[hand], score[hand] = ENGINE.place(m, XYZIN, phased=True)  # §5,§6
choose:  hand* = argmax(score)                        # better-fitting hand wins
prepare: box   = molecule_bbox(placed[hand*]) + solvent_border    # §2
         T     = flip[hand*] ∘ trim(box) ∘ frame       # ONE transform (§2a, §4)
         MAPOUT       = apply(T, MAPIN)  → valid P1 cell, composite grid, axis-norm, nxstart
         HALFMAPOUT1/2 = apply(T, HALFMAP1/2)  if supplied   # identical transform
         MASKOUT      = mask(box, placed[hand*])
         XYZOUT       = placed[hand*]         # already in the box frame
emit  :  (MAPOUT, HALFMAPOUT1/2?, MASKOUT, XYZOUT, Δ, hand*, axis_perm)
```

`Δ` is the single real-space offset (Å) that carries the whole package back to the
**original deposited map's frame** (§5), so a later step can relate the placement
to the untrimmed map or to other maps of the same particle.

`molrep_map` **does not refine** — it prepares. The refinement is the *next* task
the user runs (`servalcat_pipe` / a Coot session) on the emitted package. This
keeps the task self-contained and the responsibilities clean.

## 4. Map surgery, in gemmi

All of this is execution-time gemmi/numpy — the slim server never runs it, so it
stays off the CCP4-free request/validation/digest path (see §8).

- **Hand flip = inversion through the origin.** `chapi.flip_hand` becomes a gemmi
  point inversion: `map(x,y,z) → map(-x,-y,-z)`, i.e. reverse the grid on all
  three axes with correct origin handling. **Inverting the hand also negates the
  offset**, so the flip lives inside the same Δ bookkeeping — doing it separately
  is a classic sign error. Pin it with a parity test against
  `coot_headless_api.flip_hand` on a CCP4-aware box (§9).
- **Dual-origin bookkeeping (the trap) — and which origin our consumers read.**
  An MRC/CCP4 map pins its density in space two ways: `NXSTART/NYSTART/NZSTART`
  (grid units, CCP4-native) **and** `ORIGIN` words 50–52 (Å, MRC-2000, what
  RELION/cryoSPARC write). Crystallographic tools read `nxstart` and are **blind to
  the Å `ORIGIN`** — the classic cryo-EM registration bug. This is not a guess for
  our stack: **Moorhen/coot confirm it in source.** Moorhen's MRC parser
  (`baby-gru/src/utils/mapHeaders.ts`) reads `nxstart/nystart/nzstart` and its
  header-word loop stops at word 49 — it *never reads* the Å `ORIGIN` (words
  50–52); and the density itself is placed by coot's `shim_read_ccp4_map` →
  **clipper**, which lives in the grid/cell frame and uses `nxstart`. (Moorhen's
  `originShift`/`drawOrigin` for EM maps are `−centre-of-mass` *camera centring*,
  not registration.) **Therefore the output must express its position through the
  grid (`nxstart`) and keep map + model in one shared frame — writing only the Å
  `ORIGIN` would render the map off its model in Moorhen.** Concretely: fold *both*
  input origin conventions into one Å offset via gemmi's fractional↔orthogonal
  maths, apply the trim as an explicit offset, and emit the box so that `nxstart`
  (or a zero-origin box with the model shifted to match) places the density on the
  model, recording `Δ` back to the deposited frame.
- **Valid-cell reconstruction, composite grid, solvent border, axis-norm** — as
  §2. gemmi's grid/cell setup and `Ccp4.setup(...)` do the cell and axis work; the
  border and composite-dimension snapping are a few lines of numpy on the grid
  extents.
- **Mask** — box/soft mask around the placed molecule for servalcat.

## 5. P1 has no origin → an unphased translation function is meaningless → the target must be *phased*

This is the crystallographic crux, and it constrains **both** engines.

A cryo-EM box is **P1**, which has **continuous origin freedom** — every origin is
equivalent. A conventional (unphased) **translation function searches for the
origin-relative position of the model, and in P1 that is undefined**: there is no
FTF to run. Porting the task with a standard TF target would silently produce
nonsense.

The resolution is that **the map is phased**: a real-space map FFTs to amplitudes
*and phases*, and those phases **fix the origin**. So a **phased** search — phased
rotation function + phased translation against the map's own density — *is*
well-defined in P1. The map is the fixed, phased reference; the model is placed to
match it.

**Therefore the MR engine must be configured with a phased, map-referenced
target, never a bare P1 TF.** Concretely:

- **molrep**: `-f <map>` with the phased rotation function (`PRF`) — the source
  `molrep_map.py` already does this; carry the keyword across and confirm the
  target is the phased/map-fit one, not a default TF.
- **phaser**: see §6 — must likewise use a phased placement against the map, not
  an unphased TF.

Say this loudly in the plugin: **choosing the wrong molrep target here is the
single most likely way to ship a task that "runs" and returns garbage.**

## 6. MR engine is pluggable: molrep (default) or phaser

Abstract the placement behind a small interface —
`ENGINE.place(map_or_sfs, model, phased=True, hand, time_limit) -> (placed_model, score)` —
with two implementations:

- **`molrep`** — the proven engine, ported from the `chapi` branch. Runs `molrep`
  as a sub-step (subprocess as in the source, or a sub-plugin), phased map-fit
  target (§5), time-limited, on the down-sampled map. **This is the default and
  the reference implementation; everything else is validated against it.**
- **`phaser`** — an alternative, motivated by the Phaser stack now running
  **direct from PHIL** (the `phaser_*` PHIL tasks in `core/tasks.py`). Phaser is
  reciprocal-space MR, so the shared preprocessing would hand it the map as
  **phased structure factors in the P1 box** (FFT of the prepared box, amplitudes
  + phases) and drive it via its PHIL tree. **Open design question to settle
  before wiring:** the exact Phaser mode for placing a model against a *phased P1
  map* (its phased TF / EM path) — Phaser must, like molrep, use the phases, not
  an unphased P1 TF. Until that mode is confirmed with Phaser's own docs/PHIL,
  **implement `molrep` fully and leave `phaser` as a declared engine option with a
  clear "not yet wired — pending phased-map mode confirmation" guard**, so the
  abstraction is in place without shipping an unproven path. The PHIL-driven
  Phaser stack makes the *wiring* cheap once the mode is known; the crystallography
  is the gate.

The engine choice is a control parameter (`ENGINE = molrep | phaser`, default
`molrep`).

## 7. Registry impact (small, mostly local)

The work is almost entirely inside the new plugin directory, with three registry
touch-points:

1. **The plugin dir** — `server/ccp4i2/pipelines/molrep_map/` (pipeline, since it
   orchestrates two placement runs + preprocessing; a wrapper is fine if it stays
   self-contained): `molrep_map.py`, `molrep_map.def.xml`, `molrep_map_report.py`,
   `__init__.py`s, plus a `preprocess_map.py` (the gemmi surgery) and an
   `engines.py` (molrep/phaser interface).
2. **`server/ccp4i2/core/tasks.py`** — one `TASKS` entry (`pluginPath`,
   `defXmlPath`, `reportPath`).
3. **Frontend** — a category entry in `task-chooser.tsx`; `GenericInterface`
   (from the def.xml) is an acceptable first cut, so no bespoke `task-container`
   interface is required to ship.

`def.xml` sketch — inputs: `XYZIN` (CPdbDataFile), `MAPIN` (CMapDataFile),
`HALFMAP1`/`HALFMAP2` (CMapDataFile, *optional* — §2a); outputs: `XYZOUT`
(CPdbDataFile), `MAPOUT` (CMapDataFile), `HALFMAPOUT1`/`HALFMAPOUT2` (CMapDataFile,
emitted only when half-maps were supplied), `MASKOUT` (CMapDataFile);
controlParameters: `ENGINE`, `DOWNSAMPLE` (or `SEARCH_RESOLUTION`), `TIME_LIMIT`,
`SOLVENT_BORDER`, molrep knobs (`NMON`, `BADD`, …), and the recorded `Δ`/`hand` for
the report.

## 8. Slim-server safety

- The **request/validation/digest path stays CCP4-free.** All map surgery and the
  engine calls happen in `process()` (execution, under ccp4-python), never in
  `validity()` or the digest. `chapi`/`coot_headless_api` is **not** imported
  anywhere; the flip is gemmi.
- `validity()` may cheaply check that `MAPIN` is set and (optionally, via the
  cached diagnosis pattern) that it is a real map — but must not read the map on
  every poll or import anything CCP4.
- Execution-time binaries (`molrep`, optionally `phaser`, and the downstream
  `servalcat`) are the full-CCP4 concern, exactly like every other MR/refinement
  task — not a slim-boundary violation.

## 9. Tests — the round-trip is the gate

- **Two-sided round-trip (the acceptance test).** Place a known model into a map,
  then check both consumers: (a) real space — `XYZOUT` overlays the *original*
  (untrimmed) map to sub-voxel accuracy after applying `Δ`, **and loads registered
  in Moorhen via `nxstart`** (not the Å `ORIGIN`); (b) servalcat — feeding
  `MAPOUT` (or `HALFMAPOUT1/2` when present) + `MASKOUT` + `XYZOUT` to `refine_spa`
  moves R/FSC the right way. servalcat FFTs internally, so this needs no SF
  synthesis. One test, both failure surfaces.
- **`flip_hand` parity** — gemmi origin-inversion vs `coot_headless_api.flip_hand`
  on a sample map, on a CCP4-aware box (skip on slim), in `tests/parity/`.
- **i2run E2E** — a small demo cryo-EM `(map, model)` through the task, asserting
  the emitted package is a valid P1 cell, renderable extent, and correct frame.
- **Hand discrimination** — feed a map of known hand and assert the task picks it.

## 10. Sequence to implement

1. `preprocess_map.py` (gemmi): dual-origin → Δ, flip, trim+border, valid-cell,
   composite grid, axis-norm, mask. Unit-test each with synthetic gemmi maps
   (slim-safe) + the parity test for the flip.
2. `engines.py`: the `place()` interface + the `molrep` engine (from source).
   `phaser` declared, guarded (§6).
3. `molrep_map.py` pipeline: place-both-hands → choose → prepare → emit; report.
4. `def.xml` + `tasks.py` + task-chooser entry.
5. i2run E2E + the two-sided round-trip test.
6. Open the PR against `django` at an appropriate time.

---

*Design note, grounded in the Qt-era `molrep_map` (`martinemnoble1/ccp4i2@chapi`),
the slim-server boundary, gemmi 0.7.5's map/grid/cell API, and the PHIL-driven
Phaser stack now in `core/tasks.py`. The two decisions to confirm before or
during implementation: the exact molrep phased target (§5) and, if pursued, the
Phaser phased-map mode (§6).*
