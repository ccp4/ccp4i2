# NLP-driven CCP4i2 job construction — discussion document

> **Status.** Discussion-stage. This document is the seed for a dedicated
> implementation conversation; it is **not** a locked proposal yet. Builds
> on the architectural principles of `apps/compounds/docs/NLP_QUERY_PROPOSAL.md`
> ("Compounds NLP") and reapplies them to a much richer problem space:
> describing a CCP4i2 *job* in English and having the system construct,
> validate, and submit it.
>
> Locked decisions: 0. Open questions: many. The point of the document
> is to surface the design space cleanly, name the trade-offs, and flag
> what reuses from Compounds NLP and what is genuinely new.

---

## 1. Motivation

The Compounds NLP work proved out an architecture chemists like:

> *"CDK4 pyrimidines registered by Alice in 2025 with HTRF IC50 < 100 nM"*

→ deterministic resolution, click-to-confirm, redirect to the rich display
surface that already exists. The chemist saves a few minutes of clicking
and gets an audit-friendly artefact (the redirect URL).

CCP4i2 has a structurally similar pain point at a **different layer**:
crystallographers spend a lot of time clicking down task trees and
filling parameter forms to launch jobs whose intent is short and
explicit. A worked example from a real crystallographer prompt:

> *"Refine the output coordinates of the last molecular replacement job
> using servalcat and a dictionary made from the SMILES `c1ccccc1`."*

That sentence carries enough information to:

- pick the right task (`servalcat` family),
- locate the right input coordinates (the most-recent molecular-replacement
  job's output XYZ),
- synthesise a missing input (run `LidiaAcedrgNew` on the SMILES first to
  produce a dictionary),
- assemble it all into a runnable job,
- present it for confirm-and-submit.

Today the same intent takes ~15 clicks across two task creations. If we
can do for jobs what the Compounds NLP work did for queries, the
productivity win is real.

The hard part is that **a job is much richer than a query**. A query is
a structured *filter*; a job is a structured *workflow with side
effects*. The architectural principles carry over, but the surface area
is far larger and the failure modes are more expensive.

---

## 2. Sample prompts

A non-exhaustive sketch of what we want to handle. Each prompt's *target
shape* — the structure the LLM would need to emit — is the heart of the
design problem.

| Prompt | Notes |
|---|---|
| *"Refine the output coordinates of the last molecular replacement job using servalcat and a dictionary made from the SMILES `c1ccccc1`"* | Worked example. Task class + output reference + sub-job + literal value. |
| *"Run phaser MR on `<the imported MTZ>` with `<the deposited PDB>` as the search model"* | Two output references (or freshly-imported file references). |
| *"Solve this with phasertng using ensembles built from PDB IDs 1ABC, 2DEF, 3GHI"* | Sub-job synthesis (ensembler) with literal PDB-ID list. |
| *"Re-refine the last refmac job but with TLS turned on"* | Job-as-template — clone parameters of an existing job, override one. |
| *"Build the model from my last servalcat output"* | Output reference + downstream pipeline. |
| *"Extract a chain A-only PDB from the latest deposit and use it to MR my new dataset"* | Sub-job (PDB-tools chain extraction) + main task. |
| *"Run aimless on imported run 3"* | Recency-by-position reference into a specific job class. |

These hint at a recurring vocabulary:

- **Task class** ("molecular replacement", "servalcat", "refinement")
- **Output reference** ("the output coordinates", "the last MR job's MTZ")
- **Sub-job** ("a dictionary made from SMILES", "an ensemble from these PDB IDs")
- **Literal** ("`c1ccccc1`", "PDB IDs 1ABC, 2DEF")
- **Template-modification** ("the last refmac job but with TLS on")

That vocabulary plus a confirmation step is roughly the v1 surface.

---

## 3. What carries over from Compounds NLP

The architectural commitments of the compounds work — distilled in
`NLP_QUERY_PROPOSAL.md` §3 and §13 — apply unchanged:

| Principle | Same here? | Notes |
|---|---|---|
| LLM-as-parser, never query engine | **Same** | LLM emits a structured `JobPlan`; never invokes a job, never sees task internals beyond a curated catalogue. |
| Backend owns canonical names | **Same** | LLM types verbatim ("MR" / "refmac" / "the last MR"); resolver normalises to plugin-registry task names. |
| Discriminated-union response: Resolved / Clarify / Miss | **Same** | New entity types (task class, output ref, sub-job recipe) but the same triad. |
| Strict-mode JSON schema, every property required | **Same** | Azure OpenAI structured output, deferred imports, daily-cap, etc. all transfer wholesale. |
| Today injection on the user message | **Same** | "the last X" / "yesterday's job" still needs a Today anchor. |
| Click-to-confirm, no auto-submit in v1 | **Same, sharpened** | Cost of a wrong job is much higher than a wrong query. **Confirmation is non-negotiable in v1.** |
| Audit trail (request / resolved spec / outcome) | **Same** | Every NLP-driven job records the prompt + the resolved plan so misfires are diagnosable post-hoc. |
| Daily cap as anti-runaway guard | **Same** | Per-user, per-day; 500 default. |
| Golden-set evaluation (offline + optional online) | **Same** | Format-identical YAML harness; entries describe (prompt, expected JobPlan). |

**Modules that transfer near-wholesale**: `azure_client.py`, the
`SYSTEM_PROMPT` discipline (different content), `parse_prompt`, the
strict-mode schema-drift guards, the `test_golden.py` round-trip
pattern, the `_to_parse_result` translator pattern.

---

## 4. What's genuinely new

| Concern | Compounds NLP | CCP4i2 jobs NLP |
|---|---|---|
| **Vocabulary size** | ~5 entity types (Target, Protocol, User, Scaffold, CompoundRef) | Hundreds of task wrappers + dozens of pipelines, each with its own def.xml schema. |
| **Output shape** | A *filter* (set predicate) | A *plan* (DAG of jobs with parameters and input bindings). |
| **Reference depth** | Flat (compound IDs are just lookups) | Hierarchical (output of job X feeds input of job Y; sub-jobs synthesise inputs). |
| **State** | Stateless (each query is independent) | Project-scoped (jobs live in projects; "the last X" depends on which project). |
| **Side effects** | None (selection → redirect URL) | Compute, queue slots, file system, possibly external services (e.g. PDBe fetches). |
| **Risk** | Wrong query → wrong page contents (cheap) | Wrong job → wasted compute / wrong scientific output (expensive). |
| **Failure modes** | Mostly Miss / Clarify | Plus: invalid def.xml params, missing prerequisites, version mismatches, server-side validity rejections. |

The big shift is **plan vs. filter**. The Compounds executor intersects
sets and emits IDs; a job executor must build, validate, and submit a
multi-step plan whose nodes reference each other.

---

## 5. Architectural principle (proposed)

The Compounds principle was *"LLM emits a structured filter; backend
applies it"*. The proposed analogue:

> **The LLM emits a structured `JobPlan`; the backend deterministically
> resolves references, validates against def.xml, and submits.**

Concretely:

- The LLM **never** writes Python, never picks file paths, never invents
  CCP4i2 task names. It emits a JSON object whose enums are validated
  against a registry the backend owns.
- The backend **never** trusts a parameter blindly. Every job spec runs
  through the existing `validity()` and `runTimeValidity()` (per
  `CLAUDE.md`'s task-validation architecture) before submission.
- The user **always** sees the resolved plan before it submits — task
  class, every input binding, every prereq sub-job — and can edit or
  cancel.

The LLM-as-parser line is unchanged; what's new is that "parse intent"
now means "construct a small DAG", not "describe a filter".

---

## 6. Sketched JobPlan schema

For discussion only — the eventual schema will be hammered out in the
dedicated conversation. This is the simplest shape that covers the
worked example.

```python
@dataclass
class JobPlan:
    root: JobNode

@dataclass
class JobNode:
    task_class_as_typed: str            # verbatim user phrasing
    parameters: List[ParameterBinding]  # explicit param overrides
    input_bindings: List[InputBinding]  # how to fill each required input

@dataclass
class ParameterBinding:
    path: str                           # def.xml dot-path, e.g. "controlParameters.TLS"
    value_as_typed: str                 # LLM-typed value, backend coerces

@dataclass
class InputBinding:
    slot_as_typed: str                  # user's words, e.g. "coordinates"
    binding: BindingSpec

# Discriminated union — exactly one of these is set per binding.
BindingSpec = Union[
    OutputReference,    # "the last MR job's coordinates"
    Subjob,             # "a dictionary from SMILES c1ccccc1"
    Literal,            # explicit file path / SMILES / number
    Implicit,           # backend infers or clarifies
]

@dataclass
class OutputReference:
    source_task_class_as_typed: str       # "molecular replacement"
    position: str                         # "last" | "specific:<job_uuid>"
    output_slot_as_typed: str             # "coordinates" | "MTZ"

@dataclass
class Subjob:
    task_class_as_typed: str              # "dictionary from SMILES"
    inputs: List[InputBinding]            # nested
    parameters: List[ParameterBinding]

@dataclass
class Literal:
    value: str
    kind: str                             # "smiles" | "pdb_id" | "file_path" | "number"
```

Open shape questions are listed in §13.

---

## 7. Resolution rules (sketch)

These mirror the §6 rules of `NLP_QUERY_PROPOSAL.md` but adapted.

### 7.1 Task class
- LLM emits the user's phrase; resolver matches against the plugin
  registry (`core/task_manager/plugin_registry.py` already exists and
  enumerates ~176 plugins per `CLAUDE.md`).
- Aliases come from a curated table — *"MR"* / *"molecular replacement"*
  → `{molrep_mr, phaser_simple, phaser_pipeline, phasertng_picard, …}`.
  Single match → Resolved; multiple → Clarify chip picker; none → Miss
  with closest-name suggestions.
- The aliases table is the new analogue of the scaffold catalog: a
  small curated Python module the maintainer edits and commits.

### 7.2 Output reference
*"the last MR job"*:
1. Identify the user's project context (request-scoped — see §10).
2. Filter `Job` rows by task-class set, status=`finished`, ordered by
   `created_at` desc.
3. Recency: position=`"last"` → first row.
4. Map `output_slot_as_typed` (*"coordinates"*) onto the resolved job's
   output container — backend reads the def.xml to find candidate slots
   and matches by name / role. Multiple plausible slots → Clarify.

### 7.3 Sub-job synthesis
A small curated **recipe table** maps user-typed sub-job phrasings to
known prerequisite tasks:

| User phrasing | Recipe |
|---|---|
| *"dictionary from SMILES X"* | `LidiaAcedrgNew` with `inputData.SMILES = X` |
| *"ensemble from PDB IDs A,B,C"* | `phaser_ensembler` with PDB-ID list |
| *"chain A-only PDB from job J"* | `pdb_tools` chain extraction |

The LLM emits *the user's phrasing*, not the recipe name. The resolver
maps phrasing → recipe → concrete sub-job. Recipes that match nothing
become a Miss with suggestions of available recipes; ambiguous ones
clarify.

### 7.4 Literal
SMILES, PDB IDs, explicit file paths, numbers. Validated server-side
(SMILES via gemmi/RDKit, PDB IDs via regex + optional pdbe-fetch ping,
numbers per def.xml type). The LLM emits the typed value verbatim.

### 7.5 Implicit / defaulting
When a required input is not bound by the user's prompt:
- **Strict mode (proposed v1):** clarify always — never silently fill.
  Avoids the canonical failure mode of "I ran the refinement against
  the wrong MTZ because the system inferred."
- **Smart mode (future):** infer from project-context conventions
  (e.g. *the* MTZ if there's only one), still surface in the
  preview-and-confirm so the user can see the inference.

---

## 8. LLM contract

What the LLM sees at parse time:

1. A static system prompt (Compounds NLP discipline, restated for jobs).
2. The current project context — project name, recent jobs (e.g. last
   10 finished jobs with task class + datetime + truncated title).
   This is **non-static**, so it lives in the user message, not the
   system prompt — same trick as Today injection; system prompt stays
   cache-friendly.
3. `[Today: YYYY-MM-DD]` line.
4. The user's prompt.

What it emits:

- A `JobPlan` (with `compound_refs_as_typed`-style verbatim strings
  throughout) **or** `not_a_job` with a reason.

Discipline (carried over verbatim):

- Never canonicalise — emit user's typed names for tasks / slots / files.
- Never invent a task or slot name.
- Never emit code, parameters the user didn't mention, or guesses.
- Filler nouns ("rerun", "redo", "kick off") carry no parameter meaning.
- Display preferences ("show me the log", "open in viewer") are
  ignored — the LLM's job ends at the plan.

### 8.1 The vocabulary problem

The hardest LLM-contract question: *what is the LLM's view of the task
catalogue?*

Two extremes:

- **Full enum.** The schema includes every plugin name as an enum value.
  Hundreds of entries; bloats the schema; every catalogue change
  requires a redeploy. But the LLM can never emit an unknown task name.
- **Free-form.** The LLM emits any string; the resolver fuzzily matches.
  Schema stays small; LLM can hallucinate task names that don't exist.

**Probable v1 shape (to be debated):** *neither extreme*. Use an
**aliased class catalogue** — the LLM sees a curated list of ~30 task
**classes** (MR, refinement, model-building, density-modification, etc.)
plus their canonical names; emits a class-or-name; resolver does the
fan-out. Same trick as Compounds-NLP's Target / Scaffold resolvers —
LLM emits intent, backend owns the catalogue.

This is the open question with the largest design impact. Flagged for
the dedicated conversation.

---

## 9. Clarification loop

Same shape as Compounds NLP's Clarify continuation. New axes:

- **Task-class clarify** — *"MR"* matches multiple plugins → chip picker.
- **Output-slot clarify** — *"coordinates"* matches multiple slots in
  the resolved job's output container → chip picker.
- **Recency clarify** — *"the last MR"* with two MR jobs finished
  within minutes of each other → chip picker showing job titles + times.
- **Sub-job recipe clarify** — *"dictionary"* matches multiple recipes
  → chip picker.
- **Default-input clarify** — required input not specified, no obvious
  default → text/file picker.

The continuation pattern (`partial_plan` round-trips back, picker
choice pins via an `_id`) is identical to Compounds NLP. The pinning
fields are populated by the view, never by the LLM.

---

## 10. UI surface

Two anchors, both naturally project-scoped:

1. **In-project command bar.** Lives on the project page; submissions
   land in *that* project. This is the primary surface — it makes
   "the last X" trivially scoped.
2. **Global command bar (later).** Same UI but requires the prompt to
   name a project. Defer to v2 unless eval data justifies it earlier.

Result rendering:

- **Plan preview.** Editable form showing the resolved task class,
  every input binding (resolved file or pinned source), and every
  prereq sub-job. Submit / Cancel / Clarify buttons.
- **Submission.** "Submit" sends the plan through the existing
  job-creation API endpoint. The originating prompt is stored on the
  Job row for audit.

The aggregation-page-handoff trick from Compounds NLP doesn't apply
here — there isn't a richer downstream display surface to redirect to.
The plan preview *is* the surface.

---

## 11. Confirmation discipline

Mandatory in v1, no opt-out. Reasons:

- Cost is real (CPU minutes, occasionally hours).
- Wrong-spec failure modes can take a long time to surface (a refinement
  pointed at the wrong MTZ may run to completion and quietly produce
  garbage statistics).
- The LLM's reliability against real crystallographer prompts is unproven.

Mechanism:

- The view returns `status: "plan"` with the resolved plan.
- The frontend renders an editable preview.
- The user clicks "Submit" — only then does the job-creation endpoint
  fire. Until then, no DB write, no queue activity.

Auto-submit may become a setting once eval data shows reliability is
adequate for low-risk task classes. Locked-out for v1.

---

## 12. Audit and observability

Every NLP-driven submission writes a small `NLPJobAudit` row alongside
the Job:

- `prompt` (verbatim user input)
- `resolved_plan_json` (the canonicalised plan after resolution)
- `model_name` (e.g. `"gpt-4o"`)
- `prompt_template_version` (a hash of the system prompt at parse time)
- `n_clarify_rounds` (how many round-trips to convergence)
- `submitted_at`, `user_id`, `project_id`

This is what makes the eventual *"the LLM got this wrong"* failures
diagnosable — without it, we are debugging from prompts alone, which
is brittle.

---

## 13. Open questions

The design space is wide. The dedicated conversation should pick a
position on each before any code lands.

| # | Question | Why it matters |
|---|---|---|
| Q1 | Which task classes are in v1? Curated short list or full plugin registry? | Drives prompt-engineering effort, eval scope, and aliasing-table size. |
| Q2 | What does the LLM see of the task catalogue (§8.1)? | Largest schema-design lever. |
| Q3 | Sub-job depth — limit to one? Allow nesting? | Two-level plans cover most prompts; nested plans risk LLM rambling. |
| Q4 | Strict-clarify-everything vs. infer-and-show-in-preview? | Strict is safer; smart-defaults are ergonomically nicer. |
| Q5 | "Re-run last refmac with TLS on" — template-modification a v1 feature? | Adds significant resolver complexity (clone parameters, override one). Defer? |
| Q6 | Cross-project references — *"refine MR from project X using data from Y"* — defer? | Probably defer; project-scoped command bar makes most cases unnecessary. |
| Q7 | When the resolved "last MR job" failed — use its outputs anyway, or refuse? | Failed jobs may have partial outputs; needs explicit policy. |
| Q8 | Server-side validity — surface the existing `validity()` errors as Clarify entries? | Likely yes — same UX as resolver clarifies. |
| Q9 | Cost gating beyond the daily cap? | Some tasks (model-building pipelines) are far more expensive than others. Per-task-class budget? |
| Q10 | Where does this code live (apps/ vs. server/) given the proposed Compounds/CCP4i2 repo split? | See §14. |
| Q11 | Schema versioning — when the JobPlan grows a field, how do continuations from older clients behave? | Same problem as Compounds-NLP's `partial_selector` shape change. Solved there; same solution here. |
| Q12 | Failure-mode surfacing — when a sub-job recipe doesn't exist for the user's request, the Miss should suggest neighbouring recipes. Curate that table or compute? | Curated catalogue mirroring `substructures.py` is the obvious shape. |

---

## 14. Reuse from Compounds NLP

Given the proposed Compounds/CCP4i2 repo split (per
`apps/compounds/docs/CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md`), this
matters operationally — what *physical code* moves across the boundary?

### 14.1 Wholesale-reusable
- `azure_client.py` — Azure OpenAI factory + managed-identity wiring.
- The structured-output / strict-mode schema discipline.
- The Today-injection helper.
- The daily-cap (`LocMemCache`-backed) middleware.
- The eval-harness shape (offline round-trip + optional online).
- The discriminated-union response convention (Selection / Clarify /
  Miss / NotA{Query,Job}).

### 14.2 Same-shape, different content
- The system prompt — same discipline, different vocabulary.
- The frontend command-bar component — same shape, different result
  rendering (`PlanPreview` vs. `NLPResults`).
- The Resolved/Clarify/Miss triad on each entity — same dataclass shape,
  different entities (TaskClass, OutputRef, Recipe vs. Target,
  Protocol, User, Scaffold, CompoundRef).

### 14.3 Genuinely new
- The plugin-registry-backed task resolver.
- The output-reference resolver (job + slot together).
- The sub-job recipe catalogue.
- The plan validator (def.xml-aware).
- The plan-preview UI component.
- The submission audit row.

### 14.4 Where shared code should live

In light of the repo split:

- **Pure-infrastructure modules** (Azure client, daily-cap, eval shape,
  Today helper, structured-output guards) → a small **shared library**
  consumed by both repos. Sketch name: `nlp-core` or
  `llm_parser_kit`. Tiny surface, no domain content.
- **Domain modules** (task resolver, recipe catalogue, plan executor)
  → entirely on the CCP4i2 side post-split.
- **Compounds NLP** → entirely on the Compounds side post-split.

The two domains share *patterns*, not *vocabulary*. Trying to abstract
a unified "NLP framework" would couple the repos in a way the split
explicitly aims to undo. Better to copy the small shared bits and let
the domain code diverge.

---

## 15. Sketched slice plan

Same staged-vertical-slice approach as `NLP_QUERY_PROPOSAL.md` §18.1.
Slices land deterministic-resolver work first, defer the LLM until last,
keep eval honest from the start.

| Slice | Scope |
|---|---|
| 1 | Task-class resolver (alias catalogue + plugin-registry lookup), unit-tested against fixtures. |
| 2 | Output-reference resolver (project-scoped job lookup + def.xml output-slot mapping). |
| 3 | Sub-job recipe catalogue + recipe resolver. |
| 4 | `JobPlan` dataclass + plan executor (no sub-jobs yet, no LLM yet — pure deterministic plan-to-submission path). |
| 5 | Plan validator: surface `validity()` / `runTimeValidity()` errors as plan-level clarifies. |
| 6 | Sub-job synthesis end-to-end (depth-1). |
| 7 | Azure OpenAI client + system prompt + structured output. |
| 8 | DRF view, plan preview round-trip, audit row. |
| 9 | Frontend command bar + plan preview component. |
| 10 | Eval harness + golden set. |

Each slice is independently shippable behind a feature flag. The first
five have no LLM dependency at all — they harden the deterministic
substrate before the model gets near it.

---

## 16. Deliberately out of scope (v1)

To protect against scope creep:

- **Generative analysis.** *"Did this refinement converge?"* / *"Why is
  R-free not dropping?"* — different architecture, different data-
  visibility model for the LLM, separate proposal.
- **Workflow modification beyond plan-and-submit.** *"Re-queue all my
  failed MRs with looser tolerances"* is a multi-step workflow, not a
  single job spec.
- **Cross-project references.** Anchor everything to the current
  project until eval data justifies otherwise.
- **Fine-grained parameter tweaking via NLP** — *"Set `refmac.weight =
  0.3`"* — better solved by the existing parameter forms; NLP for that
  is uncanny-valley.
- **Auto-submit** — locked off (§11).

---

## 17. Risks and failure modes

| Risk | Mitigation |
|---|---|
| LLM hallucinates a task or slot name | Strict-mode enum on resolved task-class IDs; resolver Misses cleanly; plan preview shows real names. |
| Wrong file going into wrong slot | The plan preview surfaces every binding for confirmation; the user must click Submit. |
| Project-context confusion | Command bar is project-scoped; the prompt's *"the last X"* always resolves against *that* project's history. |
| Server-side validity rejects a synthesised plan | `runTimeValidity()` surfaces errors as plan-level Clarify entries; user fixes in the preview. |
| LLM-driven cost runaway | Daily cap (anti-runaway) + per-task-class budget under discussion (Q9). |
| Sub-job recipe mismatch | Curated recipe table, Miss surfaces near-neighbours, no synthesis without an explicit recipe. |
| LLM versions drift behind eval | Audit row records `model_name` + `prompt_template_version`; eval re-runs detect regressions. |
| Compute spent on a misconstrued prompt | Confirmation discipline (§11) — no submission without explicit Submit click. |

---

## 18. What "done for v1" might look like

A crystallographer types into a project's command bar:

> *"Refine the output coordinates of the last MR job using servalcat
> and a dictionary made from SMILES `c1ccccc1`."*

The system responds in &lt;3 s with a plan preview:

```
Plan: servalcat refinement
─────────────────────────
Inputs:
  XYZIN  ← Job 47 (phaser_simple, finished 2026-04-26 09:14)
            output: XYZOUT
  HKLIN  ← <not specified — pick one>
  DICT   ← Sub-job: LidiaAcedrgNew (acedrg)
            input: SMILES = "c1ccccc1"
Parameters:
  (defaults)

[Edit]   [Submit]   [Cancel]
```

The user picks `HKLIN` from a chip list of recent imported MTZs, clicks
Submit, and gets dropped on the running-job page. Audit row written.
Done.

That is the bar. Hitting it requires the substrate of slices 1–6 above;
the LLM is the last step, not the first.

---

## 19. Appendix — links to related docs

- Compounds NLP (canonical reference for the architecture): `apps/compounds/docs/NLP_QUERY_PROPOSAL.md`
- Task validation discipline: `CLAUDE.md` §"Task Validation: `validity()` and `runTimeValidity()`"
- Plugin registry: `core/task_manager/plugin_registry.py` (per CLAUDE.md)
- Repo split context: `apps/compounds/docs/CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md`
- Azure OpenAI infrastructure (already deployed for Compounds NLP): `NLP_QUERY_PROPOSAL.md` §11 + §16.2
