---
marp: true
theme: default
paginate: true
size: 16:9
header: 'Natural language as an interface to complex capability · Cosener''s House break-out'
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

# Natural language over a declarative interface

## Teaching a model to compose molecular graphics, and where the idea goes next

Martin Noble, Newcastle University
Break-out: *Natural language as an interface to complex capability*

<!--
30 sec. Framing line while the title is up: "The interesting bit isn't the model.
It's the SHAPE we gave the problem so a model could touch it safely — a narrow
declarative interface with all the meaning kept deterministic behind it."
-->

---

<style scoped>
section { padding-top: 40px; }
.worked { display: flex; gap: 0.6em; align-items: stretch; margin-top: 0.1em; }
.worked .col { flex: 1 1 0; display: flex; flex-direction: column; min-width: 0; }
.worked .yamlcol { flex: 1.15 1 0; }
.worked .col h3 { font-size: 0.58em; margin: 0 0 0.3em 0; letter-spacing: 0.03em;
  text-transform: uppercase; color: #666; font-weight: 700; }
.worked .prompt { font-size: 0.74em; line-height: 1.4; background: #f4f6f8;
  border-left: 4px solid #4a7cbf; border-radius: 4px; padding: 0.75em 0.85em; }
.worked .prompt b { color: #1f4e79; }
.worked pre { font-size: 8.5px; line-height: 1.32; background: #1e1e28; color: #cfcfe0;
  border-radius: 4px; padding: 0.6em 0.8em; margin: 0; flex: 1; white-space: pre;
  overflow: hidden; font-family: "SFMono-Regular", Menlo, Consolas, monospace; }
.worked pre b { color: #7fb0e6; font-weight: 600; }   /* keys */
.worked pre i { color: #c8e6a0; font-style: normal; } /* values */
.worked pre u { color: #ff8f6b; text-decoration: none; } /* colour hexes */
.worked .shot img { width: 100%; border-radius: 4px; box-shadow: 0 1px 6px rgba(0,0,0,0.25); }
.worked .arrow { align-self: center; font-size: 1.3em; color: #4a7cbf; }
.worked-cap { text-align: center; font-size: 0.7em; color: #555; margin-top: 0.55em; }
</style>

# One sentence in, a rendered scene out

<div class="worked">
<div class="col">
<h3>1 · Plain-English prompt</h3>
<div class="prompt">
"Using PDB file <b>1HCK</b>, draw <b>chain A as ribbons</b> and the <b>ATP as metaballs</b>. Residues <b>1–84</b> should be white except for the <b>C-helix (45–58)</b> which should be <b>red</b>. Residues <b>85–300</b> should be <b>golden</b>."
</div>
</div>
<div class="arrow">→</div>
<div class="col yamlcol">
<h3>2 · Generated scene (validated YAML)</h3>
<pre><b>scene:</b> <i>1HCK coloured kinase domains</i>
<b>version:</b> 1
<b>files:</b>
  - {<b>name:</b> <i>1hck</i>, <b>pdb:</b> <i>1HCK</i>}
<b>elements:</b>
  - <b>file:</b> <i>1hck</i>
    <b>representations:</b>
      - <b>style:</b> <i>CRs</i>
        <b>selection:</b> <i>//A</i>
        <b>colour:</b>
          - {<b>selection:</b> <i>//A/1-44</i>,   <b>colour:</b> <u>"#ffffff"</u>}
          - {<b>selection:</b> <i>//A/45-58</i>,  <b>colour:</b> <u>"#ff0000"</u>}
          - {<b>selection:</b> <i>//A/59-84</i>,  <b>colour:</b> <u>"#ffffff"</u>}
          - {<b>selection:</b> <i>//A/85-300</i>, <b>colour:</b> <u>"#d4af37"</u>}
      - <b>style:</b> <i>MetaBalls</i>
        <b>selection:</b> <i>//A/400</i>
<b>view:</b>
  <b>centre:</b> {<b>file:</b> <i>1hck</i>, <b>selection:</b> <i>//A/400</i>}
  <b>slab:</b>   {<b>file:</b> <i>1hck</i>, <b>selection:</b> <i>//A</i>}
  <b>background:</b> <u>"#ffffff"</u></pre>
</div>
<div class="arrow">→</div>
<div class="col shot">
<h3>3 · Rendered in Moorhen</h3>
<img src="images/nlp-declarative-moorhen-scenes.png" alt="1HCK rendered from the generated scene: chain A as ribbons with white/red/golden domain colouring and ATP as metaballs" />
</div>
</div>

<div class="worked-cap">The model writes data that fits the schema. A resolver we control turns that data into the picture. The model never issues a drawing command.</div>

<!--
1 min. Open the deck on the CONCRETE thing before any theory. Read the prompt aloud,
then: "That English becomes this YAML — and note what the YAML is: it names a PDB,
some selections, some colours. It is DATA, not instructions. The AI never told Moorhen
to draw anything; it described a state, and our resolver rendered it." Point at the
C-helix red stripe in the image as the 'it understood 45-58 is a sub-range of 1-84'
tell. This slide earns the right to the abstract 'declarative waist' framing on the
next slide — show it works, then explain why the shape matters.
-->

---

# (a) The approach: a declarative model, and two ways to write it

Rather than drive Moorhen through its imperative API (a large, fragile surface), we defined a small declarative model of a scene — a description of what is shown — and taught the model to produce that description.

- The Zod schema is the single source. It generates every other artifact, and a drift test holds them byte-identical: the JSON Schema, the human grammar, the LLM brief, and the strict structured-output schema.
- Two ways to deliver it:
  - A capable model gets the schema at the decoder (constrained decoding). The structure is then guaranteed, so almost no prompt tokens go on shape and the budget goes on intent.
  - A small or on-device model gets a compact grammar in the prompt (~1k tokens).
- A resolver we control applies the scene: validate, repair, render. The model produces a description of a state; it never issues a render command.

> What we trust is the schema and the resolver. We don't have to trust the model.

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

# (b) Why this shape — what it buys, and where it strains

An "AI drives the viewer" interface has no ground truth, no undo, and an open-ended surface. A declarative interface avoids all three: the model produces data we can check, and the guarantees live in code we own.

| Strengths | Costs |
|---|---|
| Checkable — the schema validates; bad output fails loudly | Only expresses what the schema models; a new kind of view needs a schema change |
| One source feeds every model, and can't drift | The schema has to change carefully once it is published |
| Cheap on tokens — shape at the decoder, context spent on intent | Decoder limits are real: Azure caps a strict schema at 100 properties, so we author against a pruned 83-property subset |
| Degrades gracefully — *core* fields must be honoured, *hints* are advisory | A hint can map to nothing; where it would, we leave the field out rather than fake it |

The core grammar is small on purpose — small enough that an on-device model (Apple Intelligence, Windows Copilot) is plausible, with only the harder scenes escalated to the cloud. To be clear about status: the constrained-decoding path works today on cloud models; the on-device path is designed for but not yet demonstrated.

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

<style scoped>
section { font-size: 22px; }
h1 { margin-bottom: 0.25em; }
table { font-size: 0.8em; margin: 0.35em 0; }
table th, table td { padding: 0.25em 0.55em; line-height: 1.2; }
ul { margin-top: 0.3em; }
blockquote { margin-top: 0.4em; }
</style>

# (c) Applying the same idea to CCP4i2 jobs

The same approach, at higher stakes. A scene is inert state; a job runs and has side effects. The pattern carries over, but the guarantees have to be stronger.

What stays the same: the model parses the request into a structured `JobPlan` (never code, never a task name it invents), the backend owns the catalogue of real tasks, one strict schema, and a confirmation step before anything runs.

| | Molecular graphics scene | Crystallographic job |
|---|---|---|
| Model produces | a description of state | a plan (DAG) of jobs and their inputs |
| Cost of being wrong | small, and visible at once | wasted compute, or quietly wrong science |
| Side effects | none | queue, filesystem, external fetches |
| Confirmation | optional | required — every input shown, then Submit |

- The bridge is the resolver. The model turns the user's words (*"the last MR job's coordinates"*, *"a dictionary from SMILES c1ccccc1"*) into real task names, job outputs and sub-jobs — the same job the scene resolver does for rendering.
- v1 stays cautious: the first six slices use no model at all, it asks rather than guesses, and every submission leaves an audit record. *(docs/NLP_JOB_CONSTRUCTION_DISCUSSION.md)*

> The same idea throughout: a narrow declarative interface, the backend owning the meaning, the model only reading intent, a person confirming. Scenes showed it was cheap; jobs are where it pays off.

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
