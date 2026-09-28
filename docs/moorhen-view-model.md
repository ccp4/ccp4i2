# What a Moorhen session is about, and what saving it means

**Status:** design note, for discussion. Describes existing behaviour and
proposes a way to stop it multiplying.

## The observation

There are now six ways into Moorhen:

1. campaign summary — every hit in a campaign, unsuperposed
2. campaign site — one curated site, every hit there, superposed
3. campaign run-site — one PanDDA site, every pose there, superposed
4. campaign dataset — one member selected, with its maps and its ligand
5. generic job — whatever a job produced, honouring the job's own scene if it wrote one
6. Moorhen as a task — a job whose *purpose* is the model you build in it

Six is not the problem. The problem is that "which of these am I?" is currently
spread across independent booleans — `summaryMode`, `siteParam`, `jobParam`,
`viewParam`, `selectedMemberProjectId`, `siteSummary` — so the affordances are
gated ad hoc:

```tsx
{selectedMemberProjectId && onSetVerdict && (<SiteVerdictControl … />)}   // campaign-control-panel.tsx:780
{siteSummary ? "Loaded datasets" : "Push to CCP4i2"}                       // :843
```

Six booleans describe sixty-four states, of which about seven are real. Every new
window doubles the space again, and the checks that keep the meaningless states
off the screen are scattered across the panel.

## Two axes, not one flag

A session has a **subject** — what it is about — and a **commit target** — where a
save goes and what it means. Affordances follow from those two, and every
existing window is a combination of them:

| window | subject | commit target | what a commit means |
|---|---|---|---|
| campaign summary | the campaign | none | — |
| curated site | a site × its hits | verdicts only | a judgement row |
| run site | a run's site × its poses | origin only | which pose was taken, and where it goes next |
| campaign dataset | one dataset | the campaign head for that dataset | new head, and the per-site record reconciles |
| generic job | one job's outputs | that project | an ordinary new job |
| Moorhen task | the task's inputs | **this job's outputs** | the job's result |

Two things fall straight out.

**"Add ligand here" is absent from the site views because the subject cannot carry
it, not because we withheld it.** *Which* model, and *here* relative to which
frame, are both unanswerable in a 22-pose overlay.

**The campaign dataset view and the generic job view differ only in commit
target.** Same subject, same editing affordances, different meaning of save. They
should be one viewer parameterised by where the save lands, not two pages.

## Editing is local; only commits differ

Every session that holds a molecule can edit it — that is Coot, and it is
instant. So "read-only vs editable" is the wrong cut: it is not editing that
differs between windows, it is **what a commit does**. A session with no commit
target still lets you pull a side chain around; it simply has nowhere to put it,
and should say so rather than pretend to be read-only.

## Site views judge and dispatch

The site views' primary action is a **subject switch**: "open this member in its
dataset view, carrying this event". They triage; the dataset view models.

This is worth defending on more than tidiness. In a site view every pose is drawn
in the exemplar's frame, so a modelling commit there would have to invert the fit
before storing, and getting it wrong would look perfectly correct in the view that
produced it. Keeping the modelling commit in the dataset view means **no
coordinate crosses a frame boundary on the write path at all**. See
[`pandda-pose-acceptance-design.md`](pandda-pose-acceptance-design.md) §8.

## Representation belongs to the scene

A scene authors its representations: the site scene deliberately makes the
reference a ribbon and the poses sticks. A global BONDS / RIBBONS / SURFACE
control silently destroys that authorship, and per-molecule × three buttons × 22
molecules is worse.

The list of loaded molecules already exists in the panel. Put one compact
representation control on each row there, and change the global control from
"make everything ribbons" to **"reset to the scene's representation"**.

## The panel is the subject switcher

A Moorhen page is full-window, so a panel that shows only its own subject strands
the person. Every session's panel must either embed the project browser or carry
an affordance that opens it as a modal. With subject as an explicit value, that
browser *is* the subject switcher rather than a separate navigation concept.

## What this asks for

One discriminated union for the subject, one value for the commit target, and
panel sections derived from them instead of from booleans. Then a new window —
the run-site panel is the next one — is a **variant**, not a seventh page and an
eighth boolean.

It is a bigger diff than adding the run-site panel alongside the others. It is
proposed anyway, because the alternative makes the following view worse again, and
because the frame argument above means the structure is load-bearing rather than
cosmetic.
