---
marp: true
theme: default
paginate: true
size: 16:9
header: 'The new CCP4i2 · CCP4 Developers Meeting, Cosener\'s House'
footer: 'Martin Noble · Newcastle University'
---

<!--
SPEAKER NOTES — how to use this deck
- This is a Marp deck. Preview/export with the "Marp for VS Code" extension,
  or:  npx @marp-team/marp-cli docs/presentations/ccp4i2-coseners-2026.md --pdf
- Notes for each slide live in HTML comments like this one — they don't render.
- Target: 20 min. ~18 min talk + 2 min Q&A. Slide timings are in the notes.
- Audience: CCP4 developers. They wrote/maintain the wrappers, viewers and
  libraries you're building on. Answer their two instinctive questions early:
  "did you break my task?" (no) and "can I still contribute my way?" (yes, easier).
-->

# The new CCP4i2
## A Qt-free, web-native rebuild — that keeps the crystallography

**Martin Noble**, Newcastle University
CCP4 Developers Meeting · Cosener's House

<!--
30 sec. One-liner while the title is up:
"I want to show you that we rebuilt the *plumbing* of CCP4i2 — not the science —
and that in doing so we've opened a path to retire not one but two legacy layers."
-->

---

# The thesis

> We replaced the **plumbing** — Qt, the dependency thicket, the plugin
> registry — **without touching the crystallography**.
>
> And in doing so, we can also retire the *other* legacy: **classic CCP4i**.

**Two maintenance burdens → one.**

<!--
1 min. This is the whole talk in one slide. Say the number out loud:
the win for CCP4 as a project is FEWER things to maintain, not more.
Everything that follows is evidence for this claim.
Don't frame it as "a shiny new rewrite" — this room is wary of rewrites.
-->

---

# Why we had to move

- CCP4i2 was welded to **Qt / PySide2**
  - packaging & deployment fragility across platforms
  - desktop-only; no browser, no shared access
- A heavy stack of **native dependencies** (clipper, mmdb, cctbx, CCP4 binaries)
  entangled with the *interactive* layer
- Adding a task meant a **generated registry** you had to keep in sync

<!--
2 min mark (so ~1.5 min here). This room knows the pain — don't wallow.
The point of this slide is to earn the right to the rest. Land the three
categories (UI toolkit / dependencies / registry) because slides 4, 5, 6
each kill one of them.
-->

---

# The architecture, in one picture

```
   React / Electron client            Django REST server              Worker
   (thin — no business logic)   →   (validation · execution · DB)  →  (CCP4 jobs)
        browser OR desktop              single source of truth        out-of-process
```

- **Server is the sole authority**: validation, job execution, persistence
- Client just **renders what the server reports** — same code, browser *or* desktop
- Crystallographic jobs run **out-of-process** in a worker

<!--
2 min. The developer-facing point: business logic lives in ONE place.
The client is deliberately dumb. This is what makes "browser or desktop from
one codebase" possible, and it's what makes the slim server (slide 4) possible.
-->

---

# Your investment is safe

Preserved and running **unchanged**:

- `wrappers/` · `pipelines/` · **`def.xml`** · `pimple` · `smartie` · reports

**The kicker:**

> A wrapped task with **zero frontend code** still gets a working UI —
> `GenericInterface` auto-renders it from the `def.xml`.

<!--
2 min. THE most important slide for this audience. Make eye contact here.
The fear is "my task got rewritten / abandoned." The answer is no — the
crystallographic logic is untouched, and it's now *cheaper* to surface because
the def.xml drives the UI for free. If you demo nothing else, demo a task that
has no bespoke React and still renders.
-->

---

# We swept away the dependency thicket

Rebased crystallographic utilities on **gemmi** — and where gemmi lacked it,
**ported the algorithm and deleted the dependency** (with parity tests):

| Operation | Was | Now |
|---|---|---|
| MTZ split / merge / columns | `sftools`, `splitHklout` | **in-process gemmi** |
| f′ / f″ anomalous scattering | `crossec` binary | **`gemmi.cromer_liberman()`** (~1e-2) |
| HL ↔ PHI/FOM phase conv. | clipper C++ | **gemmi + numpy** (parity-tested) |
| Unit-cell compatibility | clipper `Cell::equals` | **pure numpy** (bit-for-bit) |

**Payoff:** the interactive API runs on **stock CPython + pip `gemmi`** — *no CCP4
install*. A CI test *bans* clipper/cctbx/mmdb imports from the server surface.

<!--
3 min. This is the maintainer's-love slide — expect quiet nods.
Emphasise the discipline: we didn't just swap a library, we ported algorithms
AND wrote parity tests proving numerical agreement with clipper/crossec/chltofom.
The banned-imports CI test is the punchline — it's an enforced architectural
invariant, not an aspiration.
Files if asked: server/ccp4i2/core/conversions/{phase_data_converter,model_converter}.py,
lib/utils/formats/gemmi_split_mtz.py, tests/unit/slim/test_no_ccp4_native_imports.py
-->

---

# Registry: generated JSON → one dict entry

**Before:** scan + generate `plugin_lookup.json`, keep it in sync.

**Now:** add one entry to `core/tasks.py`.

```python
"acorn": Task(
    title="ACORN density refinement",
    pluginPath="ccp4i2.wrappers.acorn.script.acorn:acorn",
    defXmlPath="wrappers/acorn/script/acorn.def.xml",
    reportPath="...:acorn_report",   # optional
),
```

No regeneration step. Fewer moving parts.

<!--
1 min. Quick win. "Adding your program is one def.xml + one dict entry."
This directly answers "can I still contribute my way?" — yes, and it's simpler.
-->

---

# Two ways to declare a task — including the tool's own PHIL

A task's parameters can be declared by **`def.xml`** — *or* by ingesting a
tool's native **PHIL**. Both are first-class.

- `PhilPluginScript` reads a tool's `master_phil` and builds the **same**
  parameter model at runtime — rendered by the **same** UI, expert-level
  filtering for free
- **~400 PhaserTNG parameters:** ~1000 lines of hand-written def.xml
  → **~200 lines of Python** (def.xml carries only input/output files)
- Upstream adds a parameter or changes a default? The wrapper **auto-syncs** —
  no rebuild, no drift

> Wrapping **Phenix / Phaser / DIALS** stops being a maintenance tax.

<!--
1 min. The maintainer's-eyes-light-up slide for anyone who's hand-maintained
Phenix parameters. Key points to land:
- It's NOT a hack/one-off: same CData model, same React widgets, real test
  suite (tests/unit/phil/), two production wrappers (phasertng_picard/riker).
- Runtime conversion => auto-sync with upstream. That's the killer property.
- Rich CCP4i2 file types survive via "shims" (MtzFileShim, PdbFileListShim,
  AsuContentShim...) so you keep file browsers/column pickers AND get PHIL.
Files if asked: core/PhilPluginScript.py, utils/phil_to_cdata.py, utils/phil_shims.py,
wrappers/PHIL_TASK_GUIDE.md.
CUT/MERGE NOTE: if running long, fold the headline bullet into the dependency
slide (6) and drop this one.
-->

---

# The strategic move: absorb classic CCP4i

The classic **CCP4i** (Tcl/Tk/X11) is the last home of the *files-in / files-out*
idiom. We can reproduce that idiom inside CCP4i2 — and retire it.

**Why it's tractable, not hand-waving:**

- ~**181** legacy `.def` + `.tcl` task pairs
- The program contract lives in **declarative `.com` templates**, *not* in Tcl
- So: a `.com` interpreter + `.def`→`def.xml` generator **auto-ports the bulk**

> *"It can't be retired by being different from it, only by reproducing the
> idiom well enough that there is no reason to stay."*

<!--
2.5 min. Frame as a PROPOSAL to the community, not a fait accompli — this is a
dev meeting, invite discussion. The credibility comes from the mechanism: the
contract is externalised as data, so translation is mechanical for most tasks.
Source: CCP4I_CLASSIC_MODE_THINKING.md (phased plan, Phase 0 = one LABIN widget
+ 3-4 tasks to prove it end to end).
This is the slide that separates you from "just a UI refresh."
-->

---

# And it's actually richer to use

All in-browser, no Qt:

- **MTZ reciprocal-space previewer** — reflections on hk/hl/kl planes, coloured
  by intensity, resolution rings, Miller-index tooltips (Canvas)
- **Moorhen** 3D — models + electron density (WebGL/WASM)
- **MolRep self-rotation** — native **PostScript → SVG**, live in the browser
- Report graphs (Chart.js), Clustal alignment viewer, pipeline **DAG**, validation viewer

<!--
2 min (mostly the demo on the next slide). The MolRep point lands with this crowd:
we parse MolRep's PostScript self-rotation plot and re-emit it as SVG — NCS peaks
viewable in the browser with no PostScript viewer. Custom regex parser, contour
levels coloured blue→red. (server/.../wrappers/molrep_selfrot/script/ps2svg.py)
The MTZ previewer is the "ooh" — it's not a table of numbers, it's reciprocal space.
-->

---

# Demo

1. Run a task
2. Preview its MTZ output — **reciprocal space in the browser**
3. Open the result in **Moorhen** — model + map
4. Same thing, **browser tab and desktop** — one codebase

<!--
4 min. This is worth more than any slide to this audience.
HAVE A RECORDED FALLBACK ready in case of network/venue issues.
The money shot is step 4: identical result desktop and browser, same server.
Keep narration minimal — let it run.
-->

---

# Deployment — deliberately modest

> Not a cloud service. A **containerised, self-hosted option** for a
> medium-sized lab that already prefers CCP4i2 and wants shared, browser-based
> access on its own hardware.

**Complementary to CCP4Cloud — not competing.**

<!--
1 min. One slide, one sentence, said plainly. Do NOT stray into Azure/Bicep
internals — reads as someone else's problem and risks a turf reaction re CCP4Cloud.
The niche is: self-hosted, lab-scale, CCP4i2-native. Humble and useful.
-->

---

# Where this leaves CCP4

- One engine, **two front doors**: modern typed + file-oriented
- Retire **Qt**, retire (eventually) **classic CCP4i / X11**
- Interactive layer with **no native CCP4 dependency**
- Your wrappers & def.xml: **unchanged, now web-native and scriptable**

## Get involved
- Add a program = **one `def.xml` + one dict entry**
- The **CCP4i-retirement path** is open for discussion — talk to me

<!--
1 min. Close on the project-level win (fewer things to maintain) and a concrete,
low-friction invitation. Then hand to Q&A.
Repo / doc pointers to have ready: CCP4I_CLASSIC_MODE_THINKING.md, core/tasks.py.
-->

---

# Thank you

**Questions?**

Martin Noble · martin.noble@ncl.ac.uk

<!--
Likely questions to pre-load:
- "How does this relate to CCP4Cloud?" → complementary; self-hosted lab niche.
- "What about Windows?" → cross-platform; no-emoji-in-print discipline, tested.
- "Do old i2 projects load?" → (answer honestly re: DB/project import status).
- "Who maintains the gemmi ports?" → parity tests pin them to clipper/crossec behaviour.
- "Timeline for CCP4i retirement?" → it's a proposal; Phase 0 is a small proof.
-->
