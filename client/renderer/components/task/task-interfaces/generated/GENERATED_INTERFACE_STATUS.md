# Auto-Generated Task Interface Status

Generated 2026-03-08. Updated 2026-03-08 after manual fixes.

Covers the 33 interfaces originally produced by `scripts/generate_task_interfaces.py` from legacy Qt `CTaskWidget` GUI files on the `main` branch.

The generator parses `drawContents()` methods and maps `createLine(['widget', ...])` calls to `<CCP4i2TaskElement>` components. It handles `openFolder` → `CCP4i2ContainerElement` and basic `toggle` conditions. It cannot parse raw Qt widget construction (`QtWidgets.QRadioButton`, `QButtonGroup`, etc.), dynamic slot logic, or complex conditional visibility chains.

## Fixes Applied (2026-03-08)

### Enumerator toggle fixes
Replaced incorrect `useBoolToggle` with `useTaskItem` + string comparison for:
- **sculptor** — `ALIGNMENTORSEQUENCEIN` enum (`ALIGNMENT`/`SEQUENCE`): shows ALIGNIN+TARGETINDEX vs SEQUENCEIN+CHAINIDS correctly
- **buster** — `WAT` enum (`OFF`/`ON`/`MAN`): WATCYC only shown when `MAN`
- **pairef** — `SH_TYPE` enum (`semi`/`manual`/`auto`): NSHELL+WSHELL for `semi`, MANSHELL for `manual`. Also added proper boolean toggles for USE_SHAKE, AUTO_WGT, FIXED_TLS
- **mrparse** — `DATABASE` enum (`PDB`/`AFDB`/`All`): PDB settings hidden when `AFDB`, AFDB settings hidden when `PDB`, both shown for `All`. Added USEAPI boolean toggle for AFDBSEQDB
- **phaser_singleMR** — `COMP_BY` enum (`DEFAULT`/`MW`/`ASU`): ASUFILE for `ASU`, protein+nucleicacid MW for `MW`
- **phaser_mr** — `COMP_BY` + `SGALT_SELECT` enums: composition fields conditional on COMP_BY, SGALT_TEST shown only for `LIST`
- **pyphaser_mr** — Same COMP_BY + SGALT_SELECT fixes as phaser_mr

### Missing widgets and folder structure added
- **chainsaw** — Added XYZIN widget, MODE widget in new "Simple Options" folder
- **scaleit** — Added RESOLUTION_MAX widget with label, descriptive text
- **matthews** — Added labels for NRES ("Number of residues") and MOLWT ("Molecular weight"), section titles
- **edstats** — Added section titles ("Model to analyse", "Map coefficients"), descriptive labels for map types, OUTPUT_PDB_FILE label
- **pdbset_ui** — Added section subtitles ("Input Structure", "Keywords")
- **ProvideTLS** — Added descriptive text about inferring/copying TLS sets, TLSIN label
- **nautilus_build_refine** — Added XYZIN_MODE boolean toggle (XYZIN conditional), section subtitles, pipeline control warning, resolution limit label
- **MakeMonster** — Replaced broken `label` itemName with actual data type widgets (NObsObjects, Obs_1, NPhiObjects, etc.)

### Empty interfaces written from Qt source
- **import_xia2** — Added XIA2_DIRECTORY widget with advisory text
- **mtzutils** — Added Protocol (MODE), Input Data (HKLIN1, FSIGF, HKLIN2), Output Data (HKLOUT) in 3 folders

### Conditional visibility added
- **newProject_fromMerged** — MODE enum toggle: "Composition for MR" folder only shown when mode is `molrep`. Added HKLIN, NEED_USER_CELLSYMM conditional, FREER_IN_HKLIN warning, CONTAINS_SEMET

## Quality Tiers (Updated)

### Tier 1: Good — likely usable as-is (21 tasks)

These have correct toggle logic, all key widgets present, and proper folder structure.

| Task | Notes |
|------|-------|
| `buster` | Correct WAT enum toggle. All input widgets present. |
| `chainsaw` | XYZIN + ALIGNIN + MODE in proper folders. |
| `cif2mtz` | Complete. |
| `dials_image` | Simple image viewer launcher. |
| `dials_rlattice` | Reciprocal lattice viewer launcher. |
| `dui` | Minimal GUI — just launches DUI. |
| `edstats` | Full Input Data + Options folders with descriptive text. |
| `findmyseq` | Complete. |
| `imosflm` | Launches iMosflm — no input params. |
| `import_xia2` | XIA2_DIRECTORY with advisory text. |
| `matthews` | HKLIN + 3 composition methods with labels. |
| `mrparse` | Correct DATABASE enum toggle with nested USEAPI boolean. |
| `mtzutils` | Protocol + Input Data + Output Data folders. |
| `pdbset_ui` | XYZIN + multiline keywords with section titles. |
| `ProvideTLS` | XYZIN + TLSIN + TLSTEXT with descriptive text. |
| `prosmart` | Complete. |
| `scaleit` | MERGEDFILES + RESOLUTION_MAX with descriptive text. |
| `sculptor` | Correct ALIGNMENTORSEQUENCEIN enum toggle. |
| `tableone` | Complete. |
| `unique` | Complete. |
| `nautilus_build_refine` | XYZIN_MODE toggle, 3 folders including Advanced options. |

### Tier 2: Partial — functional but missing some refinements (6 tasks)

Core functionality works but some conditional visibility or advanced features not yet implemented.

| Task | What's missing |
|------|----------------|
| `MakeMonster` | Shows first data object per type only; Qt GUI dynamically shows/hides 1-6 objects per type based on count. |
| `newProject_fromMerged` | MR_MODE radio-stack (no homologs / directory / file list) not fully wired; `phasing` mode folder empty. |
| `pairef` | All 4 toggle types working (enum SH_TYPE + booleans USE_SHAKE, AUTO_WGT, FIXED_TLS). Missing USE_PREREF conditional subframe. |
| `phaser_singleMR` | Correct COMP_BY toggle. Missing LLG_COMPL toggle and some advanced options. |
| `phaser_mr` | Correct COMP_BY + SGALT_SELECT toggles. Missing MODE_TY toggle for RFILEIN/SOLIN, TARG_TRAN phases subframe. |
| `pyphaser_mr` | Correct COMP_BY + SGALT_SELECT toggles. Missing some ensemble configuration options. |

### Tier 3: Significantly incomplete — needs manual rewrite (6 tasks)

These tasks have complex Qt GUIs that the generator couldn't handle.

| Task | Why incomplete |
|------|----------------|
| `arp_warp_classic` | 48 `createLine` calls with raw Qt widgets and complex toggle chains. 0 elements extracted. |
| `buccaneer_build_refine_mr` | 64 `createLine` calls across 5 folders with 19 toggle conditions. ~15% coverage. |
| `import_serial` | Good element extraction but toggle conditions not all correctly mapped. |
| `import_serial_pipe` | Same as import_serial. |
| `molrep_den` | 40 `createLine` calls with 10 toggles across 3 folders. Multiline text missed. |
| `slicendice` | 27 `createLine` calls with 4 toggles. Multiline text missed. |

## Remaining Common Issues

### 1. Inline label+widget patterns
The generator partially handles `['label', 'text', 'widget', 'PARAM']` but misses multi-widget lines. These appear as separate elements rather than inline fields. Affected: **buster**, **pairef**, **buccaneer_build_refine_mr**.

### 2. Multiline text widgets
`createLine(['widget', '-guiMode', 'multiLine', 'PARAM'])` is parsed correctly as a widget but no hint is passed to the React component about multiline rendering. Affected: **pdbset_ui**, **ProvideTLS**, **slicendice**, **molrep_den**, **buccaneer_build_refine_mr**.

### 3. Dynamic widget counts
Some tasks dynamically create/hide widgets based on a count parameter (e.g., MakeMonster shows 1-6 data objects per type). The static React interface shows only the first. Affected: **MakeMonster**.

## Priority Recommendations for Alpha Testing

1. **Immediate use**: Tier 1 tasks (21) — deploy as-is for testing
2. **Minor polish** (< 30 min each): Tier 2 tasks (6) — mostly adding missing advanced options
3. **Schedule for manual rewrite**: Tier 3 tasks (6) — especially arp_warp_classic, buccaneer_build_refine_mr
4. **Consider removing**: Tasks that may be obsolete (dui, imosflm, dials_image, dials_rlattice are GUI-launcher wrappers that may not apply in a web context)
