# Task Interface Refactoring Plan

This document summarises the approach for simplifying and compacting each task interface, using the patterns established in the `aimless-pipe/` reference implementation.

## Toolkit

| Tool | Purpose |
|------|---------|
| `InlineField` | Replace verbose `Box`+`Typography` flex rows with `[label] [widget] [hint] [after]` |
| `useBoolToggle` | Replace `useState`/`useEffect`/`useCallback` CBoolean triads with a single hook call |
| `isTruthy` (shared) | Replace per-file `isTruthy` definitions with the shared import |
| Folder modularisation | Split 500+ line multi-tab files into `folder/index.tsx` + per-tab components |

## Priority Tiers

### Tier 1 — Large, high-impact (folder modularisation + InlineField + useBoolToggle)

#### 1. `servalcat_pipe.tsx` (1765 lines, 4 tabs)
- **Tabs**: Input Data · Parameterisation · Restraints · Advanced
- **InlineField**: ~17 flex rows with labels (weight, cycles, resolution, B-factor, etc.)
- **useBoolToggle**: Not applicable — uses `isTruthy()` inline on `useTaskItem` values (no useState/useEffect pattern). However, the local `isTruthy` definition should be replaced with the shared import.
- **Modularise**: Yes — split into 4 tab files. The orchestrator will hold ~80 `useTaskItem` calls and a large `useMemo` visibility block.
- **Special**: Uses `useFreeRWarning`, `useFileDigest`, `CChainSelectElement`. Has a complex `onFSIGFChange` callback. The Restraints tab alone is ~850 lines.
- **Estimated saving**: ~80 lines from InlineField, ~3 lines from shared isTruthy. Main benefit is maintainability from modularisation.

#### 2. `prosmart_refmac.tsx` (1413 lines, 5 tabs)
- **Tabs**: Input data · Parameterisation · Restraints · Output · Advanced
- **InlineField**: ~15 flex rows
- **useBoolToggle**: Not applicable — same pattern as servalcat (inline isTruthy). Replace local `isTruthy` with shared import.
- **Modularise**: Yes — split into 5 tab files. Very similar structure to servalcat_pipe (they share Restraints and Parameterisation sections).
- **Special**: Uses `useFreeRWarning`, `useFileDigest`, `CChainSelectElement`. Has `onFSIGFChange` callback with digest fetching.
- **Note**: servalcat_pipe and prosmart_refmac share a great deal of structure (Restraints, Parameterisation). Consider whether shared sub-components could be extracted, or at minimum ensure consistent naming.

#### 3. `crank2.tsx` (1377 lines, 3 tabs)
- **Tabs**: Input Data · Important Options · Advanced Options
- **InlineField**: ~3 flex rows (scattering coefficient rows via `ScatteringRow` component)
- **useBoolToggle**: **8 useState + 11 useEffect + 8 useCallback** — prime candidate. Items: `inputPartial`, `inputSequence`, `inputPhases`, `nonMtz`, `mad2`, `mad3`, `mad4`, `native`, `shelxcde`, `doHanddet`, `useComb`.
- **Modularise**: Yes — 3 tabs. The Input Data tab is complex with multi-dataset anomalous scattering rows.
- **Special**: Has pipeline step logic (`getBaseSteps`, `checkStartEnd`) and a custom `ScatteringRow` component. These helper functions/components should stay in the folder as shared utilities.
- **Estimated saving**: ~100+ lines from useBoolToggle alone.

### Tier 2 — Medium, good candidates (InlineField + useBoolToggle, optional modularisation)

#### 4. `acorn.tsx` (586 lines, 2 tabs)
- **Tabs**: Input Data · Advanced Acorn Parameters
- **InlineField**: ~14 flex rows — highest density of InlineField opportunities
- **useBoolToggle**: **9 useState + 9 useEffect + 9 useCallback** — major candidate. Items: `extend`, `bgrid`, `bseed`, `bresol`, `bexclude`, `becut`, `customPhase`, `custddm`, `peaksearch`.
- **Modularise**: Optional — 2 tabs, 586 lines. Could stay as single file after simplification (would shrink to ~350 lines).
- **Estimated saving**: ~120 lines from useBoolToggle, ~50 lines from InlineField.

#### 5. `arcimboldo.tsx` (498 lines, 3 tabs)
- **InlineField**: ~8 flex rows
- **useBoolToggle**: 1 pair (`litePartial`)
- **Modularise**: Optional — 3 tabs but moderate size. After InlineField it would shrink enough to stay as one file.
- **Estimated saving**: ~30 lines from InlineField, ~8 lines from useBoolToggle.

#### 6. `molrep_pipe.tsx` (478 lines, 3 tabs)
- **InlineField**: ~4 flex rows
- **useBoolToggle**: None needed
- **Modularise**: No — will shrink to ~460 lines after InlineField.
- **Estimated saving**: ~15 lines from InlineField.

#### 7. `mrbump_basic.tsx` (320 lines, 2 tabs)
- **InlineField**: ~6 flex rows
- **useBoolToggle**: **2 useState + 3 useEffect + useCallback** — items: `searchPdb`, `searchAfdb`, `includeLocal`.
- **Modularise**: No — too small.
- **Estimated saving**: ~25 lines from useBoolToggle, ~20 lines from InlineField.

### Tier 3 — Special cases (limited applicability)

#### 8. `import_merged.tsx` (1168 lines, 1 tab)
- **Single tab** ("Main inputs") — not a modularisation candidate by tabs.
- **InlineField**: ~2 flex rows — minimal benefit.
- **useBoolToggle**: Not applicable — its `useEffect` calls handle data-loading logic, not CBoolean visibility.
- **Nature**: Mostly custom UI (column selection dialogs, digest parsing, file format detection). Very domain-specific; the bulk of code is bespoke logic not reducible by our toolkit.
- **Action**: Replace local `isTruthy` with shared import. Otherwise leave as-is.

#### 9. `splitMtz.tsx` (1068 lines, 1 tab)
- **Single tab** ("Input") — not a modularisation candidate by tabs.
- **InlineField**: ~5 flex rows.
- **useBoolToggle**: 1 useState (UI toggle, not CBoolean pattern).
- **Nature**: Mostly custom UI (column transfer list, split configuration dialogs). Similar to import_merged — bespoke domain logic.
- **Action**: Apply InlineField where applicable (~20 lines saved). Otherwise leave as-is.

### Tier 4 — Already clean (no action needed)

These files are small enough or clean enough that refactoring would add complexity without meaningful benefit:

| File | Lines | Notes |
|------|-------|-------|
| `phaser_simple.tsx` | 431 | No InlineField patterns, no CBoolean pairs |
| `shelx.tsx` | 379 | Complex param sync, not CBoolean pattern |
| `phaser_pipeline.tsx` | 322 | Clean, job-change reset logic only |
| `phaser_rnp_pipeline.tsx` | 320 | Same as phaser_pipeline |
| `ProvideAsuContents.tsx` | 291 | Clean |
| `phaser_EP_LLG.tsx` | 231 | Clean |
| `modelcraft.tsx` | 231 | 1 useBoolToggle candidate but marginal |
| `phaser_simple_mine.tsx` | 228 | Clean |
| `phaser_EP.tsx` | 226 | Clean |
| `ample.tsx` | 208 | Clean |
| `molrep_selfrot.tsx` | 180 | Clean |
| `parrot.tsx` | 151 | Clean |
| `freerflag.tsx` | 124 | 1 useBoolToggle candidate but marginal |
| `LidiaAcedrgNew.tsx` | 123 | Clean |
| `xia2_multiplex.tsx` | 114 | Clean |
| `gesamt.tsx` | 94 | No tabs |
| `SIMBAD.tsx` | 93 | Clean |
| `SubstituteLigand.tsx` | 92 | No tabs |
| `auspex.tsx` | 77 | Clean |
| `ProvideSequence.tsx` | 62 | Clean |
| `csymmatch.tsx` | 61 | No tabs |
| `i2Dimple.tsx` | 59 | No tabs |
| `clustalw.tsx` | 56 | Clean |
| `generic.tsx` | 46 | Fallback renderer |

## Execution Order

1. **acorn.tsx** — Highest density of both useBoolToggle and InlineField patterns. Good warm-up before the large files. Demonstrates value clearly.
2. **crank2.tsx** — 8 useBoolToggle conversions + folder modularisation (3 tabs). Complex pipeline logic makes modularisation especially valuable.
3. **mrbump_basic.tsx** — Quick win: 3 useBoolToggle + 6 InlineField replacements.
4. **servalcat_pipe.tsx** — Largest file. Folder modularisation into 4 tabs. InlineField throughout. Shared isTruthy import.
5. **prosmart_refmac.tsx** — Similar to servalcat. 5 tabs. Consider shared sub-components with servalcat after both are modularised.
6. **arcimboldo.tsx** — Light touch: InlineField + 1 useBoolToggle.
7. **molrep_pipe.tsx** — Light touch: InlineField only.
8. **splitMtz.tsx** — Light touch: InlineField where applicable.
9. **import_merged.tsx** — Minimal: shared isTruthy import only.

## Shared isTruthy Cleanup

The following files define their own local `isTruthy`:
- `servalcat_pipe.tsx`
- `prosmart_refmac.tsx`
- `crank2.tsx`
- `acorn.tsx`

All should be updated to `import { isTruthy } from "../../task-elements/shared-hooks"`.

## Rules of Engagement

1. **Build after each file** — Run `npm run build` to verify no regressions.
2. **Preserve behaviour** — No functional changes. The rendered UI must be identical.
3. **Tab components return fragments** — `<>...</>`, never `<CCP4i2Tab>`. The orchestrator owns tab wrappers.
4. **Visibility functions in orchestrator** — All `useMemo` visibility helpers stay in `index.tsx`, passed as typed props.
5. **Don't over-modularise** — Files under ~400 lines post-simplification can stay as single files.
