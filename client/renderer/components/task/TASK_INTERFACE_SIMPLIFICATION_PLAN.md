# Task Interface Simplification — Plan

Working doc for the refactor that reduces boilerplate in task-interface files
so human authors (not just Claude) can curate and create them comfortably.

## Why

~100 task-interface files repeat the same plumbing: destructure `{ job }` from
props, call `useJob(job.id)`, then spread `{...props}` (or `job={job}`) on every
`CCP4i2TaskElement` / `CCP4i2ContainerElement`. Interfaces with several
`CBoolean` toggles also boilerplate ~6 lines of `useBoolToggle(useTaskItem, …)`
per file. None of this is complex — but it's a tax, and a confusing one for
newcomers.

A survey of ~10 interfaces identified several candidate abstractions. After
discussion we rejected "standard section components" (would impose a naming
convention on def.xml authors) and settled on **structural, name-agnostic**
abstractions only.

## Already Landed (commit 22d13c959)

| Piece | Path |
|-------|------|
| `TaskInterfaceProvider` + `useTaskInterface` | [task-interface-context.tsx](task-elements/task-interface-context.tsx) |
| `useTaskToggles([...] as const)` | [shared-hooks.ts](task-elements/shared-hooks.ts) |
| vitest + jsdom + Testing Library | [client/vitest.config.ts](../../../vitest.config.ts), [setup.ts](../../__tests__/setup.ts) |
| Narrow context-wiring test | [task-interface-context.test.tsx](../../__tests__/task-interface-context.test.tsx) |
| `test` / `test:watch` scripts | [client/package.json](../../../package.json) |

The provider + hook and `useTaskToggles` are unused in production code
at this point — no interface has been migrated yet.

## Next: The "Rip" (chunk 2)

Make the existing element primitives (`CCP4i2TaskElement`,
`CCP4i2ContainerElement`, `CListElement`, …) read `job` from context instead
of taking it as a prop, wrap every interface at the `TaskContainer` level once,
and codemod interface files to drop the now-redundant `{...props}` and
`job={job}`.

Shape:

1. **Element primitives: context-only.** In each file under [task-elements/](task-elements/)
   that declares or consumes `job: Job` (27 files identified by
   `grep "job: Job|props\.job|\{\s*job\s*[,}]"`), replace the `job` prop with
   `const { job } = useTaskInterface();`. Remove `job` from the props type
   (`CCP4i2TaskElementProps` etc.).

2. **Wrap once in `TaskContainer`.** At
   [task-interfaces/task-container.tsx:273](task-interfaces/task-container.tsx#L273),
   change `<Component job={job} />` to
   `<TaskInterfaceProvider job={job}><Component job={job} /></TaskInterfaceProvider>`.
   (Leave `job` passed to `Component` for now — interfaces still accept it in
   their props type. Can remove later when all interfaces ignore it.)

3. **Codemod interface files.** Across
   [task-interfaces/*.tsx](task-interfaces/) (~119 files):
   - Remove `job={job}` from element JSX
   - Remove `{...props}` from element JSX
   - Leave the interface component's outer signature alone — still takes
     `CCP4i2TaskInterfaceProps` (i.e. `{ job: Job }`)

   Two mechanical rules; a sed / jscodeshift pass should handle it. Spot-check
   for false positives (e.g. `{...props}` on non-element MUI components).

4. **Verify.**
   - `cd client && npx tsc --noEmit -p renderer/tsconfig.json` — TS catches any
     stale `job={job}` / props spread on an element that no longer accepts
     `job`. This is the primary safety net.
   - `cd client && npm test` — the narrow provider-wiring test must still pass.
     If the provider is mis-wired in step 2, the test fails.

### Why this ordering

Steps 1 and 2 together are "turn on context". Step 3 is cleanup. If step 3
landed first without step 1, nothing would break, but no savings materialise
either. If step 1 landed first without step 3, `{...props}` spreading an
object that includes `job` onto an element that no longer accepts `job` would
TS-error everywhere — which is actually a fine forcing function. Pick either
order.

## Decisions Locked In (don't re-litigate)

- **No "standard section" components** (e.g. `<ReflectionDataSection>`). Those
  would require plugin authors to use conventional item names in def.xml,
  which is not a guarantee the backend provides. Structural abstractions only.
- **Full rip over opt-in wrappers.** The opt-in approach
  (parallel `TaskElement` / `TaskContainer` wrappers that read from context)
  was prototyped and rejected — it doesn't realise savings until every
  interface is migrated, which defeats the point. The direct approach (change
  the primitives, codemod callers) is mechanically bigger but lands the
  benefit immediately and is fully type-checked.
- **Opt to break rather than deprecate `job?: Job`.** We considered keeping
  `job` as an optional, ignored prop during migration. Rejected: leaves a
  confusing "this prop does nothing" API around for an indefinite period.

## Known Testing Limitation

The renderer's import graph has deep fan-out and module-load side effects
(`@monaco-editor/react` calls `loader.init()` at module top;
`moorhen` / `rdkit` / `seqviz` / vis-network are heavy). Mounting a real task
interface in vitest hangs or is extremely slow because the provider/viewer
chain loads transitively through `cdatafile` → `providers/file-preview-context`.

We made the pragmatic call to **not** invest in a full "mount-every-interface"
smoke test for this refactor. The narrow context test + `tsc --noEmit` cover
the refactor's specific risks:

- TS: every stale `{...props}` / `job={job}` compile-errors
- Narrow test: catches provider wiring regressions at `TaskContainer`

Broader smoke coverage is a worthwhile follow-up, separately, and will
probably require a `__mocks__` folder convention for `providers/` and for
the heavy external deps. Not blocking this refactor.

## Candidate Follow-ups (not in scope here)

- `<ConditionalSection>` helper for the `{flag && <Box sx={{pl:3}}>…</Box>}`
  pattern that appears throughout
- Scaffold templates (a simple and a tabbed starting file) so a new interface
  doesn't start from blank
- A file-metadata extractor hook generalising the
  `fetchDigest → forceUpdate(WAVELENGTH/CELL/SG)` pattern used by ~6
  complex pipelines
- Split the two oversized interfaces
  ([crank2.tsx](task-interfaces/crank2.tsx) ~1350 LOC,
  [shelx.tsx](task-interfaces/shelx.tsx) ~1000 LOC) into folder-modules
- Properly plumbed `mount-every-interface` smoke test (see above)

## Context to Re-hydrate a Fresh Session

- Start by reading [TASK_INTERFACE_IMPLEMENTATION_GUIDE.md](task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md)
  for the current author-facing model.
- Then read the 3 files landed in `22d13c959` listed above.
- Look at any medium-complexity interface (e.g.
  [acorn.tsx](task-interfaces/acorn.tsx),
  [modelcraft.tsx](task-interfaces/modelcraft.tsx)) to see the boilerplate
  being targeted.
- The `TaskInterfaceProvider` + `useTaskInterface` + `useTaskToggles` live
  at [task-interface-context.tsx](task-elements/task-interface-context.tsx) and
  [shared-hooks.ts](task-elements/shared-hooks.ts).
