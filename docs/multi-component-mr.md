# Molecular replacement of a complex

> The one statement of how CCP4i2 places an asymmetric unit that holds
> several kinds of chain. The Phaser tasks' validity messages, the
> judgement files (`*.agent.yaml`) for mrparse, MrBUMP and the Phaser
> tasks, and the help pages all refer to this; none restates it.
> Settled 2026-10-09 after the CDK4/cyclin D1 trials (agent-evaluation.md
> §1).

## The question

The AU contents (`ProvideAsuContents`) name K kinds of chain, each with a
number of copies. Molecular replacement must account for every kind, in
enough copies. There are three ways, and a job is judged by what its
models **cover**, not by how many models it has.

## The three routes

**A. A complex template: one model holding every component, searched as
one rigid body.** A PDB entry of a homologous complex (6p8e for CDK4/cyclin
D1), or of the same complex in another crystal form, cut to the matching
chains and kept in the entry's own frame. The number of copies is the
number of copies of the *complex*. In CCP4i2: `phaser_simple_phil` with
the template as the search model and `NCOPIES`, or one `ENSEMBLES` entry
in `phaser_pipeline_phil`. First choice when such an entry exists and the
subunits are expected to sit as they do there: the search has the whole
complex's scattering behind it, and placing the components separately can
only lose signal. CDK4/cyclin D1 was solved so (`CDK4_CyclinD1` jobs 3 and
9: TFZ 27.5; two copies, TFZ 22.1 and 32.0; R-free 0.404 after ten
jelly-body cycles). mrparse writes such a template whenever hits from one
entry match more than one of the sequences it was given ("PDB complex
template: 6p8e chains B (CDK4 92%), A (CyclinD1 100%)").

**B. Separate component models, placed together in one run.** One search
model (or ensemble) per kind of chain, each with its own number of copies;
Phaser places them in the order of their expected signal, each in the
context of those already placed. In CCP4i2: `phaser_pipeline_phil`, one
`ENSEMBLES` entry per kind; or `phasertng_picard`, which takes the models
and the AU contents and chooses its own strategy (registered, runnable,
not yet judged: no judgement file, not in the chooser). The route when no
entry holds the complex, or the subunits may be arranged differently in
this crystal. BetaBlip job 4 (beta-lactamase and BLIP) is the worked
example: search TFZs 10.4 then 18.6, R-free 0.343.

**C. The components in turn.** Place one component; then search for the
next with the placed one fixed; and so on. In CCP4i2: `phaser_simple_phil`
with `INPUT_FIXED` and `XYZIN_FIXED` = the *earlier job's* `XYZOUT` (the
structure as placed in this crystal), or `phaser_pipeline_phil` with
`SOLIN` = the earlier job's `SOLOUT` and copies 0 on what is placed. The
route when B fails for one component, or when a person wants to judge each
placement before the next. It is what an unguided search of the largest
component first amounts to, and the least likely to work for a component
whose signal alone is weak.

## The rules

1. **Only a structure placed in THIS crystal is ever "already placed".** A
   model from mrparse, a database or a prediction is in its own crystal's
   frame or none, and is *searched for*. `phaser_simple_phil` refuses a
   fixed structure whose cell or point group differs from the data's, or
   that has no cell (codes 210, 221); the fixed slot is never filled from
   a previous job's files (`fromPreviousJob` False): that something is
   placed is a decision, by a person or an agent.
2. **A kind of chain is covered when some search model, fixed structure or
   earlier solution holds a chain that aligns to its sequence**
   (`lib/utils/formats/model_coverage.py`: the best gemmi alignment, with
   a positive score and at least 30 matched residues; a wrong pairing
   scores negative). The Phaser tasks warn (code 222) only for a kind
   nothing covers, say which model covers what, and name the three routes.
   They never block: placing one component now and the next later is
   route C.
3. **Items of one ensemble are alternatives for ONE component**, superposed
   beforehand; different components go in different ensembles (code 115).
   A complex template is one item, not two.
4. **Give the whole AU contents** whatever the route (`COMP_BY` ASU): the
   components not searched for still scatter.
5. **Copies are copies of what the model holds.** A complex template
   counted twice is two complexes; a single-chain model's copies are that
   chain's.

## How the tasks and their judgements express this

| Where | What it says |
|---|---|
| `phaser_pipeline_phil._check_components_searched` | rule 2's warning, naming the routes |
| `phaser_simple_phil.def.xml` `XYZIN_FIXED` | rule 1: `sameCrystalAs` and `fromPreviousJob` False |
| `phaser_simple_phil.agent.yaml` | route A with a template, route C's fixed slot; one kind of molecule otherwise |
| `phaser_pipeline_phil.agent.yaml` | route B, and A as one ensemble |
| `mrparse.agent.yaml` | templates first (A), else one ensemble per component (B); C last |
| `mrbump_basic.agent.yaml` | not for a complex to be searched together |
| help: `phaser_simple_phil`, `phaser_pipeline`, `mrparse` | the same in prose, for people |

## Known misuses, and what catches them

| Misuse | Caught by |
|---|---|
| A homologue, cut from another entry, given as the structure already placed (Haiku, 6p8e_B) | 210/221 at run time; the slot no longer fills itself |
| Two components as items of one ensemble | 115 |
| A complex template told it searches for one kind of two (the count, before 2026-10-09) | rule 2: coverage |
| One component searched, the other forgotten | 222 names it |

## Open

- A judgement for `phasertng_picard`, from evidence runs on BetaBlip and
  CDK4/cyclin D1, and its place in the chooser: the route-B variant that
  needs no decision about the search order.
- Route A when the AU holds *more* copies of one component than of
  another (CDK4 twice, cyclin once): the template covers one copy of each,
  and the extra copy is a route-B or route-C search afterwards.
  `_check_components_searched` judges kinds, not copies.
