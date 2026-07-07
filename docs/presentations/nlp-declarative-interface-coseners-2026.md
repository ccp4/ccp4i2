---
marp: true
theme: default
paginate: true
size: 16:9
header: 'Natural language as an interface to complex capability · Cosener''s House break-out'
footer: 'Martin Noble · Newcastle University'
style: |
  section { font-size: 25px; line-height: 1.3; }
  section h1 { font-size: 1.5em; line-height: 1.15; }
  section li { margin: 0.1em 0; }
  section blockquote { font-size: 0.95em; }
  section table { font-size: 0.82em; }
---

<!--
SPEAKER NOTES — break-out contribution
- Three-slide worked example for the "NL as an interface to complex capability"
  break-out. Marp deck; preview/export with "Marp for VS Code" or:
    npx @marp-team/marp-cli docs/presentations/nlp-declarative-interface-coseners-2026.md --pdf
- The claim isn't "we bolted an LLM on". It's a REUSABLE RECIPE: put a narrow,
  declarative, versioned model at the waist; let the AI emit DATA against it;
  keep all semantics deterministic on the backend. Molecular graphics is the
  proof-of-concept; crystallographic job construction is the same recipe at a
  harder point on the risk/state axis (slide 3).
- Source docs if pressed: client/renderer/MOORHEN_SCENES_SCHEMA_V1_DESIGN.md
  (esp. §3 stratification, §4 conformance rule, §7a/§12 AI authoring),
  docs/NLP_JOB_CONSTRUCTION_DISCUSSION.md.
-->

# Natural language over a declarative waist
## A worked example: teaching AIs to compose molecular graphics — and where it goes next

**Martin Noble**, Newcastle University
Break-out: *Natural language as an interface to complex capability*

<!--
30 sec. Framing line while the title is up: "The interesting bit isn't the model.
It's the SHAPE we gave the problem so a model could touch it safely — a narrow
declarative interface with all the meaning kept deterministic behind it."
-->

---

# (a) The approach: a declarative model, then two ways to teach an AI to write it

We did **not** teach an AI to *drive* Moorhen (imperative API, huge fragile surface).
We built a **tight declarative model of a scene** — *what is shown*, not *how to draw it* — and taught AIs to **emit that data**.

- **One source of truth** (a Zod schema) → *generates* every artifact, drift-tested to byte-equality:
  machine **JSON Schema** · human **grammar** · compact **LLM brief** · strict **structured-output schema**
- **Two delivery paths, one source:**
  - **Capable model** → feed the schema to the **decoder** (constrained decoding): structure *guaranteed*, ~zero prompt tokens on shape, tokens spent on *intent*
  - **Small / on-device model** → a compact in-prompt **mini-grammar** (~1k tokens)
- A **deterministic resolver** applies the scene — validate / repair / render. **The AI emits data describing a state; it never issues render commands.**

> The trust boundary is the schema + resolver, not the model.

<!--
2 min. This is the corrected version of the one-liner ("declarative model + tight
prompt"). Two refinements to land out loud:
(1) the prompt is ONE generated artifact of several — the same Zod source also
    emits the JSON Schema, the human docs and the strict structured-output schema,
    all guarded by a byte-equality drift test. Single source => they can't diverge.
(2) for capable models the schema goes to the DECODER not the prompt (constrained
    decoding) — shape is guaranteed, so the context budget buys intent, not syntax.
Declarative-not-imperative is the load-bearing choice: a scene is idempotent state
with no side effects, so a wrong scene is cheap and instantly visible. Contrast
with "AI drives the viewer", which has no ground truth and no safe undo.
-->

---

# (b) Why this shape — its strengths, and its edges

**Why we arrived here.** An imperative "AI drives the viewer" interface has no ground truth, no undo, and a boundless surface. A **declarative waist** inverts that: the model emits *verifiable data*, and every hard guarantee lives in code we own.

| Strengths | Edges / weaknesses |
|---|---|
| **Verifiable** — schema validates; bad output fails loudly, not silently | **Expressiveness ceiling** — only what the schema models; novel viz needs a schema change |
| **Model-agnostic & drift-free** — one source feeds every model | **Schema must evolve carefully** once published/upstreamed |
| **Token-economical** — shape via decoder, context spent on intent | **Constrained-decoding caps are real** — Azure's 100-property strict limit forced a pruned 83-property authoring subset |
| **Graceful degradation** — *core* must-honour vs *hints* advisory (the "wrong vs merely plainer" test) | **"Fiction" risk** — a hint that maps to nothing; we defer such fields rather than fake them |

**Desktop / laptop fit.** Compactness is a *design constraint*, not an afterthought: the core grammar is small enough to *aspire* to an **on-device model** (Apple Intelligence, Windows Copilot) — **tiered prompting** escalates only complex scenes to the cloud. Honest status: **capable-model-via-constrained-decoding is the proven path today; on-device is measured-and-plausible, not yet shipped.**

<!--
2.5 min. The conformance split (slide's "core vs hints") is the elegant bit worth
dwelling on: the dividing test is "would a renderer that ignored this produce a
WRONG image or merely a PLAINER one?" Wrong => core, must-honour, carries a unit
(Å). Plainer => hint, advisory, safely ignorable. That single rule is what lets
the same scene target a full renderer and a weak one and degrade sensibly.
On desktop suitability — be candid: the compactness work (terse keys, strong
defaults, an explicit small authoring subset) is REAL and the brief's token count
is MEASURED, so "fits on-device" is a number not a hope. But we prove it on capable
cloud models via constrained decoding; the on-device path is designed-for, not yet
demonstrated end to end. Don't oversell.
-->

---

# (c) Mapping it onto CCP4i2's crystallographic capability

The same recipe, moved to a harder point on the risk axis. A **scene is declarative state**; a **job is a plan with side effects** — so the pattern carries, but the guarantees must be stronger.

**What transfers unchanged** — LLM-as-parser (emits a structured `JobPlan`, never code, never a task name it invents) · backend owns the canonical catalogue · single-source strict schema · **confirm-before-act**.

| | Molecular graphics scene | Crystallographic job |
|---|---|---|
| AI emits | declarative **state** | a **plan / DAG** of jobs + bindings |
| Cost of being wrong | cheap, instantly visible | wasted compute, *quietly* wrong science |
| Side effects | none (idempotent) | queue, filesystem, external fetches |
| Confirmation | optional | **non-negotiable** — preview every binding, then Submit |

- **The bridge:** the LLM types the user's words (*"the last MR job's coordinates"*, *"a dictionary from SMILES c1ccccc1"*); a **deterministic resolver** maps intent → real task names, job outputs, sub-jobs — exactly as the scene resolver maps a scene to renderer calls.
- **v1 is deliberately safe:** resolver-first (slices 1–6 have *no* LLM), strict clarify over silent inference, an audit row per submission. *(docs/NLP_JOB_CONSTRUCTION_DISCUSSION.md)*

> One recipe: narrow declarative interface · backend owns the semantics · AI only parses intent · human confirms. Scenes proved it cheap; jobs apply it where it pays.

<!--
2.5 min. This is the "so what for CCP4" payoff. The key insight to voice: the scene
work isn't a one-off toy — it's the low-risk end of a SPECTRUM. Same architecture
(declarative waist, backend owns semantics, LLM parses intent, deterministic
resolver, confirm loop) reappears for job construction, where the stakes are higher
so the confirm step becomes mandatory and the resolver does more.
The worked crystallographer prompt to say aloud if there's time: "Refine the output
coordinates of the last MR job using servalcat and a dictionary made from SMILES
c1ccccc1" — that one sentence names a task, an output reference, a synthesised
sub-job (acedrg on the SMILES) and a literal; the resolver assembles the DAG and
shows it for confirmation. ~15 clicks -> one sentence + one Submit.
Close on the reusable-recipe line; that's the transferable idea the break-out is
actually about.
-->
