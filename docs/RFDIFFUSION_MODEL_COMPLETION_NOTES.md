# RFDiffusion for Crystallographic Model Completion — Design Notes

**Status:** exploratory design notes / discussion record. Nothing built.
**Provenance:** distilled from a conversation (2026-05-21) plus Randy Read's response
on the CCP4 forum, the reasoned replies to it, and a close read of the ROCKET paper
(2026-06-26).
**Owner:** Martin Noble (martin.noble@ncl.ac.uk)

**Key paper read:** Fadini, Li, McCoy, … Jovine, Terwilliger, Read, Hekstra &
AlQuraishi, *"AlphaFold as a prior: experimental structure determination conditioned
on a pretrained neural network,"* Nat. Methods 2026
(https://doi.org/10.1038/s41592-026-03047-4; OA preprint bioRxiv 2025.02.18.638828).
Quotations below are from this paper.

---

## 1. The idea in one paragraph

Use a denoising-diffusion structural prior (RFDiffusion) as a **density-guided
main-chain builder / loop completer**. The experimental map enters sampling not as
a hard constraint but through a **potential** — a custom auxiliary term whose
gradient is added to the diffusion score at each reverse step — that scores
agreement between a synthetic density rendered from the current backbone and the
observed map (real-space CC or a likelihood-style target). The already-built part
of the model is fixed as a motif (partial diffusion); the missing/incorrect
stretches are diffused into the density.

This is **guidance of a generative prior**, not reinforcement learning. RFDiffusion
is a diffusion model; the "policy you can shape" analogue is its `potentials/`
framework, which is genuinely designed for users to drop in new terms.

> ### Core thesis — the focusing rationale: *excisability*
>
> The single reason to reach for RFDiffusion **specifically** (rather than diffusion
> in general, or ROCKET/OpenFold) is that the loop-completion problem is naturally
> **local**, and RFDiffusion is the one tool whose architecture lets you *act*
> locally:
>
> - **Conditioning is by coordinates, not coevolution.** Motif scaffolding / partial
>   diffusion *is* "fix this local context, build into the gap." You hand it the
>   neighbourhood you already trust as fixed coordinates — exactly what a partial
>   crystallographic model gives you.
> - **The network is already partly spatial-local.** Its SE(3) structure track
>   typically runs on a **kNN neighbour graph**; only the 2D pair track is global, and
>   even that shrinks when you excise (smaller N). So a local input is not a hack —
>   it's close to how the model already attends.
> - **ROCKET structurally cannot do this.** Its lever is the **MSA cluster profile, a
>   whole-chain coevolutionary object**. There is no meaningful MSA of a spatial
>   sub-domain sphere; cut below domain level and the signal it optimises evaporates.
>   ROCKET is therefore pinned at whole-chain/whole-domain granularity (and its
>   ~500-residue/40 GB memory wall, §5).
>
> So *excisability is an architectural property of RFDiffusion, not a generic property
> of "using a density potential."* It is simultaneously the **scientific** argument
> (you can work at the natural scale of the problem) and the **cost** argument (small
> on-manifold inputs, forward-only, no MSA axis — §5).
>
> **Honesty guard, so it survives scrutiny:** the advantage is *relative*, not
> absolute. RFDiffusion was still trained on whole chains, so an excision wants a
> buffer (~loop + the secondary structure it packs against + a shell), should respect
> chain-continuity / secondary-structure boundaries, and should let the *local density*
> — not a larger prior input — carry long-range and crystal-contact information (§5).
> The claim is "**far more excisable than ROCKET**," not "freely excisable."

---

## 2. Architecture, and why the differentiability worry dissolves

- **Per residue the diffused quantity is an SE(3) frame** — a Cα position `t` plus
  a rotation `R` for the local N–Cα–C frame. N, C, O and Cβ are then fixed rigid
  offsets in that frame: `pos = t + R · ideal_local_xyz`. Ideal bond
  lengths/angles, peptide planarity at ω ≈ 180°.
- That reconstruction is **fully continuous and differentiable** in both `t` and
  `R`. There is *no* peptide-flip / Cβ-basin discontinuity, because the orientation
  is part of what is diffused — it is never *inferred* from Cα–Cα–Cα–Cα pseudo-
  dihedrals (which is where a pure-Cα reconstructor, PULCHRA-style, would have
  bimodal basins and gradient jumps). The diffusion already pays the cost of
  carrying orientation; we just read backbone+Cβ out for free.
- **Cost argument (the reason to prefer Cα-frame RFDiffusion over all-atom):**
  Cα → mainchain+Cβ is essentially deterministic rigid geometry, so generate
  main chain only with the cheap Cα-frame network; render synthetic density from
  backbone+Cβ; score; backprop to the frames. First-order sidechain placement is a
  fast **post-hoc scoring** step, kept *out* of the gradient (rotamers are discrete).
- **Discreteness that remains, and where to park it:** sidechain χ-basins (post-hoc),
  cis/trans peptides (assume trans), and sequence→trace register assignment (a
  separate register-finding step, not part of the sampling gradient).

### Competitive framing (important — these are *different* problems)

| Approach | What it does | Competes with |
|----------|--------------|---------------|
| RFDiffusion as **density-guided tracer** | builds backbone into density, little/no starting model | Buccaneer, ARP/wARP autobuild |
| Rocket / OpenFold-refinement | bends an existing *prediction* toward the data | refinement / remodelling |

Conflating these is the main source of confusion in the "do we need a new tool"
debate (§4). They overlap in the loop-completion case but start from different
objects: a generative prior sampling conformations vs. a prediction prior being
optimised toward data.

---

## 3. Prior art

- **ROCKET** = *Refining OpenFold with Crystallographic/cryo-EM likelihood targets*
  (Fadini et al., Nat. Methods 2026) — **the most relevant prior art, and the
  incumbent.** Precisely what it does, from the paper:
  - **Optimises MSA cluster-profile embeddings, not coordinates.** At inference time
    (no retraining) it learns *"multiplicative and additive adjustments to MSA
    cluster profiles that maximize their agreement with experimental data,"* by
    gradient descent — *"structure optimization in the space of coevolutionary
    embeddings rather than Cartesian coordinates."* This is the lever because MSA
    depth/statistics determine AF2's geometry.
  - **Differentiable reciprocal-space likelihoods.** X-ray target `Lxtal` depends on
    observed intensities and their errors `(Io, σI)` and model structure-factor
    amplitudes `Fc`; cryo-EM target `Lcryo` uses complex Fourier terms from the two
    half-maps `(F1, F2)` and model `Fc`. `σA` refined in resolution bins each
    iteration. Initial MR/cryo-EM docking to place the model; a final `phenix.refine`
    polish for local geometry.
  - **Resolutions demonstrated:** 27 high-res (<3 Å) X-ray structures (matches
    phenix.refine / PredictAndBuild / ModelCraft); reduced-resolution cryo-EM series
    at 6/8/10 Å (made by adding noise to half-maps until FSC matches target); frontier
    cases including a 9.6 Å subtomogram average, a 3.82 Å X-ray multidomain case, a
    3.4 Å preferred-orientation cryo-EM complex, and **Jovine's 8.6 Å ZPD egg-coat
    filament, later validated against an independent 4.6 Å map**.
  - **It is a point estimate, not an ensemble** — gradient descent to one answer
    (multiple MSA-subsampled starts / independent 100-iteration traces help escape
    bad starts, but the output is a single model).
  - **Documented failure modes we should care about:** gradient descent gets
    *"trapped in a local minimum"* (GroEL at 6.8 Å — the correct conformation scores
    a *higher* likelihood, LLG 1583 vs 1448, but is never reached); it *"can fail to
    flip small loops (3–4 residues in length) that contain bulky side chains"*; one
    chain at a time; not aware of crystal contacts.
  - Available without licensing friction.
- **OpenFold3 + Alisia Fadini (AlQuraishi lab)** — Fadini is **ROCKET's first author**
  and has moved to AlQuraishi; the ROCKET paper itself names the next step:
  *"integrating generative models to account for conformational ensembles"* and
  learning an amortised *"mapping from experimental observables to a profile bias
  matrix."* OpenFold3 is diffusion-based. This is the single most important fact for
  the "why RFDiffusion?" question (§4): it neutralises "because diffusion" as a
  differentiator **and** means the diffusion-for-ensembles direction is the ROCKET
  team's own roadmap, not a rival programme.
- **DiffModeler** (Kihara lab) — diffusion-based model building for cryo-EM maps;
  closest existing template for "diffusion + density."
- **Chroma** (Generate Biomedicines) — generative protein model designed around
  classifier-style conditioners; arguably an easier substrate for density guidance
  than RFDiffusion proper.
- The mainstream X-ray trajectory remains "predict (AlphaFold) → dock → refine"
  (predicted-model-aware autobuild in Phenix; Buccaneer + AF templates). A
  diffusion-guided completer competes with that, not with autobuild-from-scratch.

---

## 4. Randy Read's questions, and reasoned responses

Randy framed three loop-modelling scenarios (all assume the experimental
conformation is one a predictor *won't* readily produce — otherwise you'd just
refine from the prediction) and pressed two hard questions: **do we need a new tool
given Rocket exists**, and **why expect RFDiffusion to beat OpenFold/AF2**. These
are fair, and the honest answers are mostly conceding his framing while isolating a
genuinely differentiated niche.

### Scenario 1 — reasonable resolution, loop well-ordered
Rocket and many building tools already handle this. **No new tool is justified
here.** A density-guided diffusion tracer would, at best, match existing methods.

### Scenario 2 — loop has no density (poorly ordered)
There is no single answer; the right object is an **ensemble of conformations
consistent with the data and with protein-structure priors**. Randy correctly says
this is "not really loop completion." **This is the one place a generative diffusion
prior is the architecturally natural tool:** guided reverse diffusion samples from
(roughly) prior × data-likelihood, so multiple trajectories give a *posterior
ensemble* directly. ROCKET is by construction a **point-estimate gradient-descent
optimiser** in MSA-embedding space, and its authors name *"integrating generative
models to account for conformational ensembles"* as the explicit next step. So the
differentiated value of a diffusion approach is ensembles — exactly the active area
Randy, Fadini and the upcoming CASP are all moving toward, and which ROCKET itself
does not yet address.

### Scenario 3 — low resolution (≲3.5 Å X-ray, ≲4 Å cryo-EM), loop not disordered
This was a **main emphasis of the Rocket paper**, and Rocket demonstrably works here.
The honest position: **for scenario 3, Rocket is the incumbent with published
results, and the burden of proof is on any diffusion-prior approach to show a case
Rocket cannot solve.** Two concrete openings the paper itself exposes: (a) ROCKET's
gradient descent demonstrably **gets trapped in local minima** even when the correct
conformation scores a higher likelihood (the GroEL 6.8 Å case, LLG 1583 vs 1448) — a
stochastic sampler is, in principle, less prone to this; and (b) it *"can fail to flip
small loops (3–4 residues) that contain bulky side chains"* — a specific loop-building
gap. Both are hypotheses for where a diffusion sampler might win, to be demonstrated,
not assumed.

### "Do we need a new tool, given Rocket?"
For scenarios 1 and 3 as stated: **probably not** — Rocket (and the OpenFold3
successor) targets them directly. The defensible reasons to still pursue a
diffusion-prior approach are narrow and specific:

1. **Ensembles (scenario 2)** — generative sampling is the more principled object
   for data-consistent conformational diversity.
2. **De novo designed proteins** (the actual day-job of Martin's group). ROCKET's
   entire mechanism is *reshaping coevolutionary signal* — its c-Abl analysis shows
   it works by altering MSA mutual-information patterns, and it explicitly leverages
   the fact that *"MSAs may implicitly encode multiple conformations that a protein
   can adopt."* For a designed protein, the **only** alignment available is a
   synthetic, single-conformation profile: the design pipeline (RFDiffusion backbone →
   **ProteinMPNN** → AF2 filter) does yield many scaffold-compatible sequences, but
   they are all sampled `p(sequence | one design backbone)` — effectively a fancy PSSM
   with some structure-induced coupling, **not** a deep natural MSA spanning the fold's
   conformational landscape. So an MSA-lever method gets almost none of the
   alternative-conformation signal it relies on, and what coupling it has points *at*
   the designed conformation — i.e. toward the state the crystal may be contradicting.
   *(Note: do not say "there is no MSA" — that overstates it; the precise claim is that
   the available MSA is synthetic and single-state, and a secondary point is that MPNN
   sequences are themselves off-distribution for AF2's natural-MSA channel.)* RFDiffusion
   sidesteps the question entirely: no MSA needed, and its prior is *exactly the one the
   protein was designed under*. This is the strongest, most concrete "why RFDiffusion
   not OpenFold" answer tied to a real, recurring use case.

The honest flip side, worth stating to Randy: an unconditioned generative prior that
"bends more freely" also carries *less* protein-specific knowledge to stay physical
when the density is weak — which is precisely scenario 3's risk. Freedom to move and
reliability under weak data are in tension.

### "Why expect RFDiffusion to beat OpenFold (AF2 reimplementation)?"
- **Not "because diffusion."** OpenFold3 already brings a diffusion framework, and
  Fadini is steering it at Rocket's goals. "Diffusion" is therefore not a
  differentiator.
- What *could* still differentiate RFDiffusion is that its prior is built for
  **generation under constraints** — partial diffusion, motif scaffolding, symmetry,
  and a mature `potentials/` interface for arbitrary guidance — rather than for
  prediction from sequence+MSA. For MSA-poor / designed cases this matched prior is
  a real advantage. But this is partly a *transient engineering* advantage (ease of
  dropping in a density potential today), which OpenFold3 may erase.
- Net: **be honest that the case for RFDiffusion specifically is narrow** — designed
  proteins, MSA-poor targets, and ensemble generation — and that for mainstream
  low-res remodelling Rocket/OpenFold3 are the sensible default.

### The benchmark / test-case problem
Randy's point — that you need both (a) a case where prediction fails without
experimental data and (b) genuine low-res data with a higher-res ground truth — is
**real, tool-agnostic, and the actual gating resource for this whole research line.**
It penalises Rocket, OpenFold3 and any RFDiffusion approach equally. Constructive
contributions ccp4i2 could make:

- **Reuse Randy's cryo-EM signal-degradation recipe for X-ray:** truncate a trusted
  high-res structure's data to low resolution and add realistic noise to the
  structure factors, giving controlled (if simulated) test cases with known ground
  truth. Not "genuine" low-res data, but reproducible and immediately available.
- **Mine the PDB for genuine pairs** instead of trawling Diamond/ESRF raw archives:
  look for same-construct entries where an early low-resolution deposition was later
  superseded by a higher-resolution one — cross-reference by UniProt/sequence +
  deposition date + "supersedes"/related-entry records. This is more tractable than
  asking beamlines to mine raw-data troves, and ccp4i2 already has PDB/PDB-REDO
  fetch infrastructure in its test fixtures to assemble a small curated set.
- The de novo niche partly sidesteps the problem: for designed proteins we often
  *have* the intended design model as a reference and know predictors struggle, so
  constructing honest test cases is easier.

### "Just run Rocket on your candidate cases"
This is the **correct first experiment**, and the collegial answer is yes. Before
building anything: take the candidate structures where RFDiffusion is expected to
help, run Rocket (and, when available, the Fadini/OpenFold3 tool) on them. If Rocket
solves them, the RFDiffusion tool is unjustified for that regime. Pursue the
diffusion-prior path only where a case is demonstrated in which the
prediction-refinement framing fails but generative-prior + density guidance
plausibly would not — most likely the designed-protein and ensemble regimes above.

---

## 5. Computational cost — does diffusion buy you anything?

A recurring intuition: a diffusion loop-builder feels "local to the loop and crystal"
and therefore cheaper / lower hardware demand. The intuition lands on a real saving,
but **not for the reason "locality" suggests** — so it's worth being precise, because
Randy will probe it.

**Where "local → cheaper" does *not* hold.** "Local to the loop" does not mean the
network only sees the loop. Both RFDiffusion and OpenFold are global-attention
transformers; a loop's conformation depends on global context and the prior is a
whole-structure network, so the 2D/pair track scales ~O(N²) in *total* residues, not
loop length. You cannot amputate the protein to the loop and run cheaply. The density
term *is* local and cheap (rendering synthetic density for a few residues is trivial),
but that's true of every method's data term — it isn't where the cost lives. The
hardware *floor* (one decent GPU, single domain) is the same for both.

**Can you excise a crystal volume to make the prior input small?** Yes — and this is
the right instinct — but as a *spectrum*, governed by "how far off the prior's
manifold does the cut take me?" Two senses, opposite answers:
- *The data term:* excise freely. You only ever render/score density in a box around
  the loop; already local and cheap, and not where cost lives.
- *The prior input:* damage scales with how violent the cut is. **Excise a whole
  domain** → fine; domains are on-manifold and ROCKET *already* does this (it works one
  chain/domain at a time because of the 500-residue wall). **Excise loop + the
  secondary structure it packs against + a ~12–15 Å buffer, fixed as a motif** → the
  workable sweet spot: 50–150 residues, well under any memory wall, fast per sample.
  **Excise a tight sphere that severs secondary structure / drops long-range
  contacts** → too far off-manifold (fake severed termini, exposed core, and you may
  have deleted the residues that *determine* the loop).
Two reasons this favours the diffusion route specifically: (1) excision is
**architecturally native to RFDiffusion** (motif scaffolding *is* "fix local context,
build the gap"; its SE(3) structure track is typically a kNN neighbour graph, so it is
already partly spatial-local — only the pair track is global), whereas ROCKET *cannot*
excise below domain level because its lever, the MSA cluster profile, is a whole-chain
coevolutionary object. (2) For crystals, the loop is often fixed by **crystal contacts**
and the paper notes OpenFold *"is not explicitly aware of crystal contacts."* Clean
resolution: excise loop + buffer for the *prior* (small, on-manifold), but score
against the *full local experimental density*, which already encodes the symmetry-mate
contacts — the data term supplies what the excision removed. So the accessibility win is
real, provided "crystal volume" means *a coherent local neighbourhood with a buffer,
fixed as a motif*, not a raw sphere.

**Where the saving is real — and ROCKET's own paper quantifies it.** ROCKET's cost is
dominated by **backpropagation through OpenFold every iteration**, and the authors
state the consequence plainly: *"iterative backpropagation through OpenFold is memory
intensive and limits the maximum size of the protein or domain that can be refined at
once — about 500 residues on a 40-GB A100 GPU."* That ceiling comes from holding
activations for the backward pass. A diffusion sampler runs the network **forward
only** — the density-potential gradient is autodiffed through the cheap
reconstruction+render step, *not* through the network — so it has no backward-pass
activation store and a much higher residue ceiling for the same VRAM. This is the
honest, evidenced version of your intuition: the win is **forward-only vs
backprop-through-the-predictor**, not loop-locality. Three components:

1. **Forward-only.** No backward pass through the prior → lower peak memory (ROCKET's
   actual wall) and ~half the per-step compute.
2. **Partial diffusion starts near the answer.** Fix the built part, noise only the
   loop to a low level, run few reverse steps — fewer network passes per result.
3. **No MSA axis.** RFDiffusion conditions on structure, not an MSA cluster profile,
   so it deletes the dimension ROCKET both depends on and pays for.

**Where it compounds: ensembles.** Forward-only sampling makes N data-consistent
conformations N cheap, batchable runs; getting N diverse answers from ROCKET means N
backprop-heavy optimisation runs. So the per-sample saving multiplies exactly in the
scenario-2 regime that is the differentiated use case.

**The deeper accessibility point.** If "local and cheap" is the real driver, the
honest competitor isn't ROCKET at all — it's **classical density-guided loop building**
(CCD / kinematic closure / Rosetta loop modelling, Coot's loop fitting): genuinely
local, CPU-only, near-zero hardware demand. A learned global prior — diffusion or
OpenFold — only earns its GPU cost when the loop is long/ambiguous, you need a
data-consistent *ensemble*, or there's no MSA (designed proteins). The ladder:

| Rung | Tool | Hardware | Earns its place when |
|------|------|----------|----------------------|
| Cheapest, most local | Coot / Rosetta density-guided loop closure | CPU | short/medium loops, decent density |
| Middle | RFDiffusion + density potential | 1 GPU, **forward-only** | long/ambiguous loops, **ensembles**, **no MSA** |
| Heaviest | ROCKET / OpenFold(3) refinement | 1 GPU, **backprop + MSA**, ≤~500 res / 40 GB | low-res remodelling where coevolutionary signal helps |

---

## 6. Proposed protocol — density-conditioned inpainting of a loop in a fixed motif

The concrete proposal: treat loop completion as **density-guided partial diffusion of
a loop within an excised local neighbourhood that is fixed as a motif.** One sentence:
*excise loop + the structure it packs against + a ~12–15 Å buffer, fix it as a motif,
and diffuse the loop into the local density.* Small enough to be cheap and fast, large
enough to stay on-manifold (§1 thesis, §5 cost).

Step by step:

1. **Target.** Identify the loop residues to rebuild (wrong/poor-fitting stretch).
   Sequence is known (crystallographic completion).
2. **Excise the neighbourhood.** Take every residue with an atom within ~12–15 Å of
   any loop atom — this captures the secondary structure the loop packs against. Then
   **snap to structural boundaries**: never cut mid-helix/strand; include the whole SSE
   if it is partially in the shell. Result is typically ~50–150 residues.
3. **Fix as motif, mark loop as diffused.** Partial diffusion: the neighbourhood is
   held as fixed coordinates (the motif); only the loop (plus ~1 stem residue of
   overlap each side for closure) is noised. The gap is built to connect the two fixed
   anchor frames — this is RFDiffusion's native scaffolding/inpainting mode, so
   loop-closure to fixed anchors is handled by the diffusion process, not a separate
   closure step.
4. **Noise schedule.** Start from a low/moderate noise level (you begin near the
   answer) → few reverse steps. Tune.
5. **Density potential, added to the score each reverse step.** Reconstruct
   backbone+Cβ from the SE(3) frames (deterministic, differentiable); render synthetic
   density at the map resolution (gemmi); score against the **full local observed
   density** in a box around the loop — real-space CC or a likelihood target. Use the
   full local density even though the prior *input* was excised: it implicitly carries
   crystal-contact / symmetry-mate information the excision dropped. Backprop the
   scalar to the frames through the reconstruction+render only — **not** through the
   network (forward-only network passes; this is the §5 cost win).
6. **Sample an ensemble.** Run N independent trajectories → a posterior ensemble of
   loop conformations consistent with the local density. This is the scenario-2
   value-add ROCKET does not provide.
7. **Post-hoc per sample.** First-order rotamer placement (sequence known) kept *out*
   of the gradient; score each conformer against density; rank. Optional final
   geometry polish in the full model (cf. ROCKET's `phenix.refine` step).
8. **Stitch back.** Graft the chosen conformer(s) onto the fixed stems in the full
   structure; refine in the complete crystal context.

**Design decisions still open (flag these to Randy / decide empirically):**
- **Buffer radius** (12–15 Å is a starting guess) — trades on-manifold fidelity
  against cost; tune on the validation set.
- **Crystal contacts:** rely on the local density to encode them (cleaner, preferred)
  vs. explicitly including contacting symmetry-mate residues in the motif (more
  faithful prior input, but those are separate chains and themselves somewhat
  off-distribution). Default to density-only first.
- **Cα-frame vs all-atom** for the motif context — Cα-frame for the diffused loop is
  settled (cheap, backbone+Cβ is enough density signal); whether the *fixed* motif
  benefits from all-atom context for a richer density comparison is worth a test.
- **Stem closure quality** — verify the built loop closes cleanly to both anchors
  without strained geometry; this is the classic loop-closure failure point and the
  one thing to watch in early runs.

**Validation (see §7 next steps):** take a trusted high-res structure,
blank a loop, excise the neighbourhood, condition on (optionally degraded) density,
and check the rebuild against ground truth; and run the de novo cases where the prior
is matched.

---

## 7. Where this leaves us — recommended next steps

1. **Assemble candidate cases** (ideally de novo designs solved crystallographically,
   and any low-res-then-high-res pairs) and **run Rocket on them first.** Let the
   incumbent define the bar.
2. **Build the simulated low-res benchmark** via resolution truncation + structure-
   factor noise from trusted high-res structures (Randy's degradation idea, applied
   to X-ray), using ccp4i2's existing fetch infra.
3. **Only if (1) leaves a gap:** prototype the density potential against RFDiffusion
   in the designed-protein / ensemble regime, where the prior is matched and the
   value (a posterior ensemble of data-consistent conformations) is something Rocket
   does not naturally produce.
4. Keep talking to Randy/Fadini — the OpenFold3 diffusion + ensemble work may simply
   be the better vehicle, in which case the contribution is the **benchmark and the
   crystallographic-integration / Moorhen-facing UX**, not a competing core method.

---

## 8. Environment to prototype (when it comes to that)

Ignore the Container Apps environment for the experimentation phase — wrong tool for
the "edit a potential, re-run, look at output" loop. Fastest path to results:

- **Needs:** one NVIDIA GPU ≥16 GB VRAM (A100 40/80 GB comfortable; T4 works for
  small jobs), ~5 GB disk for code+weights, inference only (no multi-GPU). Use the
  **Baker lab's official RFDiffusion Docker image** — the SE(3)-Transformer custom
  kernels are notoriously fragile to install by hand.
- **Ranked options:**
  1. RunPod / Lambda / Vast.ai notebook with an A100 (~$1–2/hr) — running within an
     hour, cheapest in *your* time for a 50–200 GPU-hr prototype. Do this first.
  2. Azure ML compute instance, `Standard_NC24ads_A100_v4` (single A100 80 GB) in
     UK South — stays in-tenancy; ~$3–4/hr, shut down when idle.
  3. Plain Azure NC/ND VM + the Docker image — same, more DIY, good for a long-lived
     dev box.
  4. Azure Container Apps serverless GPU (T4/A100 consumption) — the eventual
     *serving* destination if a working potential becomes a Moorhen-callable
     endpoint, but a poor fit for interactive prototyping (cold-start image pulls,
     stateless abstraction). Verify UK South availability and NC/ND_A100 quota
     before committing budget.
- **Effort gradient:** Day 0–1, run unmodified RFDiffusion inference to confirm the
  stack. Day 1–N, the 90% of real work: write the density-agreement potential in
  `rfdiffusion/potentials/` — render synthetic density from current backbone+Cβ
  (gemmi), score against the observed map, return a scalar torch can differentiate.
  Validation: take a trusted refined structure, blank a loop, fix the rest as a
  motif scaffold, condition on the experimental map, and see if it rebuilds the loop.
