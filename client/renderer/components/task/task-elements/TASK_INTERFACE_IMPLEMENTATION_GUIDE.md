# CCP4i2 Task Interface — Definitive Guide

This is the single reference for building task interfaces in CCP4i2. It covers everything from understanding a task's narrative description to producing a fully working, reactive interface.

---

## Table of Contents

1. [Workflow Overview](#workflow-overview)
2. [Anatomy of a Task Interface](#anatomy-of-a-task-interface)
3. [Step 1 — Understand the Task](#step-1--understand-the-task)
4. [Step 2 — Scaffold the Interface](#step-2--scaffold-the-interface)
5. [Step 3 — Layout and Grouping](#step-3--layout-and-grouping)
6. [Step 4 — Conditional Visibility](#step-4--conditional-visibility)
7. [Step 5 — Reactive Behaviour (onChange, auto-sync)](#step-5--reactive-behaviour-onchange-auto-sync)
8. [Step 6 — Register and Test](#step-6--register-and-test)
9. [Props Reference](#props-reference)
10. [Layout Components](#layout-components)
11. [Common Patterns Cookbook](#common-patterns-cookbook)
12. [Pitfalls and Hard-Won Lessons](#pitfalls-and-hard-won-lessons)
13. [Complete Worked Example — ModelCraft](#complete-worked-example--modelcraft)
14. [Complete Worked Example — ProSMART-Refmac (Multi-Tab, Digest-Driven)](#complete-worked-example--prosmart-refmac-multi-tab-digest-driven)

---

## Workflow Overview

```
  Input source            (narrative description, screenshots, or both)
         │
         ▼
  Read the .def.xml      (parameter names, types, defaults, qualifiers)
         │
         ▼
  Look at generated/     (auto-generated baseline — useful for parameter names)
         │
         ▼
  Write the interface    (group, label, show/hide, react)
         │
         ▼
  Register in task-container.tsx
         │
         ▼
  Test: toggle booleans, load files, check validation borders
```

---

## Anatomy of a Task Interface

Every custom interface is a React functional component that receives a single prop:

```tsx
import { CCP4i2TaskInterfaceProps } from "./task-container";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  // ... build the UI
  return <Paper>...</Paper>;
};

export default TaskInterface;
```

`CCP4i2TaskInterfaceProps` contains one field: `job: Job` (with `job.id`, `job.task_name`, etc.).

---

## Step 1 — Understand the Task

Before writing any code:

1. **Read the `.def.xml`** file for the task (in `wrappers/`, `wrappers2/`, or `pipelines/`). This defines every parameter: its name, type (`CBoolean`, `CObsDataFile`, `CInt`, etc.), default value, and qualifiers (`guiLabel`, `toolTip`, `min`, `max`, enumerators).

2. **Open the generated interface** (if one exists) in `task-interfaces/generated/<taskname>.tsx`. This shows which parameters the legacy GUI exposed and in what order. It's a useful starting point but will be flat and unstyled.

3. **Run the task in the generic interface** to see all the parameters rendered automatically. Note which groups of parameters belong together logically.

### Working from Screenshots

Claude can work from **screenshots of the legacy Qt interface** to phenocopy the layout, groupings, labels, and conditional visibility logic. This is often the most efficient approach for complex interfaces:

1. **Provide one screenshot per tab** — Claude will identify all visible widgets, labels, checkboxes, dropdowns, and inline text patterns.
2. **Provide screenshots with dependent fields revealed** — Toggle checkboxes and expand sections in the old interface to show conditional elements, then screenshot those states. Claude needs to see both the collapsed and expanded states to implement visibility logic correctly.
3. **Note any dynamic content** — If a section's visibility depends on the input data (e.g. "no nucleotide chains" messages that depend on the coordinate file composition), describe this relationship explicitly.

Claude cross-references the screenshots against the `.def.xml` to map visual elements to parameter names, types, and containers. This is how `prosmart_refmac.tsx` was built — from 5 tab screenshots plus an additional screenshot showing the Restraints tab with dependent elements revealed.

---

## Step 2 — Scaffold the Interface

Create a new file: `task-interfaces/<taskname>.tsx`

Minimal scaffold:

```tsx
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2TaskElement itemName="F_SIGF" {...props} />
      <CCP4i2TaskElement itemName="XYZIN" {...props} />
    </Paper>
  );
};

export default TaskInterface;
```

**Key rule:** always spread `{...props}` onto every `CCP4i2TaskElement` and `CCP4i2ContainerElement`. This passes the `job` (and any future props) down.

---

## Step 3 — Layout and Grouping

### Containers

Use `CCP4i2ContainerElement` with an empty `itemName=""` to create visual grouping without tying to a real container in the data model:

```tsx
<CCP4i2ContainerElement
  {...props}
  itemName=""
  qualifiers={{ guiLabel: "Reflection data" }}
  containerHint="FolderLevel"    // Collapsible card with header
>
  <CCP4i2TaskElement itemName="F_SIGF" {...props} />
  <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
</CCP4i2ContainerElement>
```

| `containerHint` | Appearance |
|-----------------|------------|
| `FolderLevel`   | Card with collapsible header (default for top-level sections) |
| `BlockLevel`    | Subtle bordered box with a label (for sub-sections) |
| `RowLevel`      | Horizontal flow layout |

### Tabs

For interfaces with many parameters, use tabs:

```tsx
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";

<CCP4i2Tabs {...props}>
  <CCP4i2Tab label="Input data">
    {/* Input fields */}
  </CCP4i2Tab>
  <CCP4i2Tab label="Options">
    {/* Options fields */}
  </CCP4i2Tab>
</CCP4i2Tabs>
```

### Size Constants

Defined in `field-sizes.ts`:

| Size | Value | Use Case |
|------|-------|----------|
| `xs` | 8rem  | Single digits, cell parameters |
| `sm` | 12rem | Short numbers, small enums |
| `md` | 20rem | Default for standalone fields |
| `lg` | 32rem | File paths, longer text |
| `full` | 100% | Full container width |

### FieldRow — Horizontal Layout for Fields

```tsx
import { FieldRow } from "../task-elements/field-row";

// Equal-width fields side-by-side
<FieldRow>
  <CCP4i2TaskElement itemName="PARAM_A" {...props} />
  <CCP4i2TaskElement itemName="PARAM_B" {...props} />
</FieldRow>

// Constrained-width fields
<FieldRow equalWidth={false} size="xs">
  <CCP4i2TaskElement itemName="CELL_A" {...props} qualifiers={{ guiLabel: "a" }} />
  <CCP4i2TaskElement itemName="CELL_B" {...props} qualifiers={{ guiLabel: "b" }} />
  <CCP4i2TaskElement itemName="CELL_C" {...props} qualifiers={{ guiLabel: "c" }} />
</FieldRow>
```

### Grid2 — Precise Column Control

```tsx
import { Grid2 } from "@mui/material";

<Grid2 container spacing={2}>
  <Grid2 size={{ xs: 12, md: 4 }}>
    <CCP4i2TaskElement itemName="SPACEGROUP" {...props} qualifiers={{ guiLabel: "Space group" }} />
  </Grid2>
  <Grid2 size={{ xs: 12, md: 8 }}>
    <CCP4i2TaskElement itemName="UNITCELL" {...props} />
  </Grid2>
</Grid2>
```

### Inline Text + Field Pattern

For "Run for ___ cycles" style layouts:

```tsx
<Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
  <Typography variant="body1">Run for</Typography>
  <Box sx={{ width: "8rem" }}>
    <CCP4i2TaskElement itemName="CYCLES" {...props} qualifiers={{ guiLabel: " " }} />
  </Box>
  <Typography variant="body1">cycles</Typography>
</Box>
```

Use `guiLabel: " "` (a space) to suppress the label above the field.

### What NOT to Use for Field Layout

```tsx
// DON'T: Stack doesn't size fields properly
<Stack direction="row">
  <CCP4i2TaskElement itemName="FIELD1" {...props} />
</Stack>

// DON'T: Add minWidth directly to fields
<CCP4i2TaskElement sx={{ minWidth: "10rem" }} itemName="FIELD1" {...props} />

// DON'T: Use the deprecated elementSx on containers
<CCP4i2ContainerElement elementSx={{ width: "8rem" }} />
```

Use `FieldRow` or `Grid2` instead.

---

## Step 4 — Conditional Visibility

### Simple Visibility (via `visibility` prop)

For showing/hiding fields based on enum or boolean values. The `visibility` prop takes a function returning `true` (show) or `false` (hide):

```tsx
const { value: MODE } = useTaskItem("MODE");

<CCP4i2TaskElement
  itemName="REFERENCE_FILE"
  {...props}
  visibility={() => MODE === "MATCH"}
/>
```

### Boolean Conditional Rendering (via onChange + local state)

**This is the recommended pattern for CBoolean-driven visibility.** It provides immediate UI response:

```tsx
import { useCallback, useEffect, useState } from "react";

const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { value: MY_TOGGLE_RAW } = useTaskItem("MY_TOGGLE");

  // Local state mirrors the server value
  const [myToggle, setMyToggle] = useState(() => isTruthy(MY_TOGGLE_RAW));

  // Sync from server for programmatic changes (initial load, parameter file import)
  useEffect(() => {
    setMyToggle(isTruthy(MY_TOGGLE_RAW));
  }, [MY_TOGGLE_RAW]);

  // onChange fires when the user clicks the checkbox
  const handleToggle = useCallback(async (updatedItem: any) => {
    setMyToggle(isTruthy(updatedItem._value));
  }, []);

  return (
    <>
      <CCP4i2TaskElement
        itemName="MY_TOGGLE"
        {...props}
        onChange={handleToggle}
      />
      {myToggle && (
        <CCP4i2TaskElement itemName="DEPENDENT_FIELD" {...props} />
      )}
    </>
  );
};
```

**Why local state?** The SWR cache patch path (`mutateContainer` -> `patchContainer`) is reliable for data but React may not always re-render every subscriber. Local state gives deterministic, immediate visibility toggling.

**Why `isTruthy()`?** CBoolean values arrive as JS `true`/`false` from the server, but some code paths may return the strings `"True"`/`"False"`. Always normalize.

### Organizing Visibility for Complex Interfaces

```tsx
const { value: mode } = useTaskItem("MODE");
const { value: mapSharp } = useTaskItem("MAP_SHARP");
const { value: mapSharpCustom } = useTaskItem("MAP_SHARP_CUSTOM");

const visibility = useMemo(() => ({
  showMapSharpening: () => mapSharp === true,
  showCustomBFactor: () => mapSharp === true && mapSharpCustom === true,
  isMatchMode: () => mode === "MATCH",
}), [mode, mapSharp, mapSharpCustom]);

// Usage
<CCP4i2TaskElement
  itemName="BSHARP"
  {...props}
  visibility={visibility.showCustomBFactor}
/>
```

### Visibility on Containers

Containers also accept `visibility`:

```tsx
<CCP4i2ContainerElement
  {...props}
  itemName=""
  qualifiers={{ guiLabel: "NCS model" }}
  containerHint="BlockLevel"
  visibility={() => xyzinMode !== "no"}
>
  {/* children */}
</CCP4i2ContainerElement>
```

### Disabled vs Visibility

| Use Case | Pattern |
|----------|---------|
| Closely related option (user should see it exists) | `disabled` |
| Completely different mode/section | `visibility` |
| Cascading: immediate children | `disabled`; deeper descendants -> `visibility` |

```tsx
<CCP4i2TaskElement
  itemName="REFERENCE_TYPE"
  {...props}
  disabled={() => mode !== "MATCH"}
/>
```

---

## Step 5 — Reactive Behaviour (onChange, auto-sync)

### onChange — Responding to User Actions

The `onChange` callback fires when a parameter is successfully set on the server. It receives the full serialized item:

```tsx
<CCP4i2TaskElement
  itemName="F_SIGF"
  {...props}
  onChange={handleFileChange}
/>
```

Common uses:
- **Toggle visibility** (see Step 4 above)
- **Extract metadata from files** (see below)
- **Clear dependent fields** when a parent changes

### Extracting File Metadata with fetchDigest

When the user selects a file, extract metadata (wavelength, cell, spacegroup):

```tsx
const { useTaskItem, fetchDigest } = useJob(job.id);
const { item: F_SIGFItem } = useTaskItem("F_SIGF");
const { forceUpdate: setWavelength } = useTaskItem("WAVELENGTH");
const { forceUpdate: setSpaceGroup } = useTaskItem("SPACEGROUP");

const handleFileChange = useCallback(async () => {
  if (!F_SIGFItem?._objectPath) return;
  const digest = await fetchDigest(F_SIGFItem._objectPath);
  if (!digest) return;

  if (digest.spaceGroup) await setSpaceGroup(digest.spaceGroup.replace(/\s+/g, ""));
  if (digest.wavelength > 0) await setWavelength(digest.wavelength);
}, [F_SIGFItem?._objectPath, fetchDigest, setSpaceGroup, setWavelength]);
```

Available digest fields (reflection data):

| Field | Type | Description |
|-------|------|-------------|
| `spaceGroup` | `string` | e.g. "P212121" |
| `cell` | `{ a, b, c, alpha, beta, gamma }` | Unit cell |
| `wavelength` | `number` | Single wavelength |
| `wavelengths` | `number[]` | Multi-wavelength |
| `resolution` | `{ low, high }` | Resolution range |
| `contentFlag` | `number` | 1=anom I, 2=anom F, 3=mean I, 4=mean F |
| `hasFreeR` | `boolean` | FreeR flags present |

**Two ways to access digests:**

| Method | When to use |
|--------|------------|
| `fetchDigest(objectPath)` | Inside `onChange` handlers — imperative, one-shot (e.g. extract wavelength when file is selected) |
| `useFileDigest(objectPath)` | For render-time data — reactive SWR hook, re-renders when digest changes (e.g. chain composition for visibility) |

Both are returned from `useJob(job.id)`. See the [Digest-Driven Composition Visibility](#digest-driven-composition-visibility) cookbook pattern for `useFileDigest` usage with coordinate files.

### Auto-Syncing Parameters with `syncTo`

When one parameter's value should be automatically derived from another:

```tsx
const { syncTo: syncToAimlessRef } = useTaskItem("REFERENCE_FOR_AIMLESS");
const { value: modeValue } = useTaskItem("MODE");

// REFERENCE_FOR_AIMLESS = true when MODE is "MATCH", false otherwise
useEffect(() => {
  syncToAimlessRef(modeValue === "MATCH");
}, [modeValue, syncToAimlessRef]);

// Don't render REFERENCE_FOR_AIMLESS in the UI — it's derived
```

`syncTo` handles bounce prevention internally (pending-update tracking). Always prefer it over raw `update` for auto-sync effects.

### forceUpdate vs update vs syncTo

| Function | Use Case |
|----------|----------|
| `forceUpdate` | Programmatic updates from code (file change handlers, digest processing). Updates server + patches SWR cache + triggers validation. |
| `update` | Raw update function. Rarely needed directly. |
| `syncTo` | Auto-sync one parameter from another. Handles bounce prevention. |

---

## Step 6 — Register and Test

### 1. Add to task-container.tsx

```tsx
// At the top — import
import MyTaskInterface from "./mytask";

// In the switch statement
case "mytask":
  return <MyTaskInterface job={job} />;
```

### 2. Test Checklist

- [ ] All fields render with correct labels
- [ ] File selectors accept the right file types
- [ ] Boolean toggles show/hide dependent fields **in both directions**
- [ ] Validation borders appear (red = error, green = valid)
- [ ] Enum/autocomplete fields show correct options
- [ ] `forceUpdate` calls propagate (e.g. after file metadata extraction)
- [ ] TypeScript compiles: `cd client && npx tsc --noEmit -p renderer/tsconfig.json`

---

## Props Reference

### CCP4i2TaskElementProps

```tsx
interface CCP4i2TaskElementProps {
  job: Job;                              // Required — the job object
  itemName: string;                      // Parameter name (e.g. "F_SIGF", "CYCLES")
  sx?: SxProps<Theme>;                   // MUI sx styling
  visibility?: boolean | (() => boolean); // Show/hide control
  disabled?: boolean | (() => boolean);   // Enable/disable control
  qualifiers?: any;                       // Override guiLabel, toolTip, etc.
  onChange?: (updatedItem: any) => void;  // Fires after successful server update
  suppressMutations?: boolean;           // Skip SWR cache mutations (for bulk operations)
}
```

### CCP4i2ContainerElementProps

```tsx
interface CCP4i2ContainerElementProps extends CCP4i2TaskElementProps {
  initiallyOpen?: boolean;                          // Start collapsed/expanded (default: true)
  containerHint?: "FolderLevel" | "BlockLevel" | "RowLevel";
  excludeItems?: string[];                          // Hide specific child items
}
```

### Key qualifiers

| Qualifier | Effect |
|-----------|--------|
| `guiLabel` | Label shown above the field |
| `toolTip` | Tooltip text on hover |
| `guiMode: "radio"` | Render enum as radio buttons instead of dropdown |
| `min`, `max` | Numeric bounds |
| `enumerators` | Override available choices |

---

## Layout Components

### Quick Reference

| Want to... | Use |
|------------|-----|
| Collapsible section | `<CCP4i2ContainerElement containerHint="FolderLevel">` |
| Sub-section with border | `<CCP4i2ContainerElement containerHint="BlockLevel">` |
| Two fields side-by-side, equal width | `<FieldRow>` |
| Multiple small fields in a row | `<FieldRow equalWidth={false} size="xs">` |
| Specific column proportions | `<Grid2 size={{ xs: 4 }}>` |
| Responsive breakpoints | `<Grid2 size={{ xs: 12, md: 6 }}>` |
| Display text + icons in a row | `<Stack direction="row">` |
| Full-width field | Just use `<CCP4i2TaskElement>` (default is full-width) |
| Tabs | `<CCP4i2Tabs>` + `<CCP4i2Tab>` |

### Composite Widgets (Pre-Built)

These render automatically based on the parameter's CData type:

| Data Type | Widget | Layout |
|-----------|--------|--------|
| `CCell` | `CCellElement` | 6 fields (a,b,c,alpha,beta,gamma) in a row |
| `CReindexOperator` | `CReindexOperatorElement` | Matrix values in a row |
| `CEnsemble` | `CEnsembleElement` | Grid with copies, label, use fields |
| `CImportUnmerged` | `CImportUnmergedElement` | File + metadata grid |

### Standalone Widgets (Explicitly Used)

These are used directly in task interfaces, not auto-dispatched:

| Widget | Import | Use Case |
|--------|--------|----------|
| `CChainSelectElement` | `../task-elements/cchainselect` | Multi-select chain IDs as chips, backed by a comma-separated CString |

---

## Common Patterns Cookbook

### Checkbox with Inline Text

```tsx
<Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
  <CCP4i2TaskElement
    itemName="AUTO_STOP"
    {...props}
    qualifiers={{ guiLabel: " " }}
    sx={{ width: "auto" }}
  />
  <Typography variant="body1">
    Stop automatically if R-free does not improve in
  </Typography>
  <Box sx={{ width: "8rem" }}>
    <CCP4i2TaskElement itemName="STOP_CYCLES" {...props} qualifiers={{ guiLabel: " " }} />
  </Box>
  <Typography variant="body1">cycles</Typography>
</Box>
```

### Conditional Section Based on Enum

```tsx
const { value: XYZIN_MODE } = useTaskItem("XYZIN_MODE");

<CCP4i2TaskElement itemName="XYZIN_MODE" {...props} qualifiers={{ guiLabel: "Use of NCS" }} />

<CCP4i2ContainerElement
  {...props}
  itemName=""
  qualifiers={{ guiLabel: "Model from which to infer NCS" }}
  containerHint="BlockLevel"
  visibility={() => XYZIN_MODE !== "no"}
>
  <CCP4i2TaskElement
    itemName="XYZIN_HA"
    {...props}
    visibility={() => XYZIN_MODE === "ha"}
  />
  <CCP4i2TaskElement
    itemName="XYZIN_MR"
    {...props}
    visibility={() => XYZIN_MODE === "mr"}
  />
</CCP4i2ContainerElement>
```

### File Selection with Metadata Extraction

```tsx
const { useTaskItem, fetchDigest, mutateValidation } = useJob(job.id);
const { item: HKLINItem } = useTaskItem("HKLIN");
const { forceUpdate: setSpaceGroup } = useTaskItem("SPACEGROUP");
const { forceUpdate: setCell } = useTaskItem("UNITCELL");

const handleHKLIN = useCallback(async () => {
  if (!HKLINItem?._objectPath) return;
  const digest = await fetchDigest(HKLINItem._objectPath);
  if (!digest) return;
  if (digest.spaceGroup) await setSpaceGroup(digest.spaceGroup.replace(/\s+/g, ""));
  if (digest.cell) await setCell(digest.cell);
  await mutateValidation();
}, [HKLINItem?._objectPath, fetchDigest, setSpaceGroup, setCell, mutateValidation]);

<CCP4i2TaskElement itemName="HKLIN" {...props} onChange={handleHKLIN} />
```

### Side-by-Side with Disabled State

```tsx
<FieldRow>
  <CCP4i2TaskElement
    {...props}
    itemName="MODE"
    qualifiers={{ guiLabel: "Pipeline mode" }}
  />
  <CCP4i2TaskElement
    {...props}
    itemName="REFERENCE_DATASET"
    qualifiers={{ guiLabel: "Reference type" }}
    disabled={() => mode !== "MATCH"}
  />
</FieldRow>
```

### Push Default Values on Boolean Toggle

When checking/unchecking a boolean should set a sensible default on a sibling field (but still allow the user to override afterward):

```tsx
const { forceUpdate: forceSetCYCLES } = useTaskItem("CYCLES");

const handleBASIC = useCallback(
  async (updatedItem: any) => {
    await forceSetCYCLES(isTruthy(updatedItem._value) ? 5 : 25);
  },
  [forceSetCYCLES]
);

<CCP4i2TaskElement
  itemName="BASIC"
  {...props}
  qualifiers={{ guiLabel: "Run a quicker basic pipeline" }}
  onChange={handleBASIC}
/>
```

This is a one-shot default — the user can still manually change CYCLES afterward. Use this instead of `syncTo` when the relationship is "suggest a default" rather than "always derive."

### Digest-Driven Composition Visibility

When a section's visibility depends on the composition of an input coordinate file (e.g. "show ProSMART protein options only if the model has protein chains"), use `useFileDigest` to reactively read the XYZIN digest:

```tsx
const { useTaskItem, useFileDigest } = useJob(job.id);
const { item: XYZINItem, value: XYZINValue } = useTaskItem("XYZIN");

// Only fetch digest when a file is loaded
const xyzinDigestPath =
  XYZINValue?.dbFileId && XYZINItem?._objectPath ? XYZINItem._objectPath : "";
const { data: xyzinDigest } = useFileDigest(xyzinDigestPath);
const composition = xyzinDigest?.composition;

const hasProteinChains =
  Array.isArray(composition?.peptides) && composition.peptides.length > 0;
const hasNucleotideChains =
  Array.isArray(composition?.nucleics) && composition.nucleics.length > 0;

// In JSX:
{!hasProteinChains && (
  <Typography variant="body2" sx={{ fontStyle: "italic" }}>
    Input atomic model contains no protein chains
  </Typography>
)}
{hasProteinChains && (
  <>{/* Show protein-specific options */}</>
)}
```

Available composition fields from a CPdbDataFile digest:

| Field | Type | Description |
|-------|------|-------------|
| `composition.chains` | `string[]` | All chain IDs |
| `composition.peptides` | `string[]` | Protein chain IDs |
| `composition.nucleics` | `string[]` | Nucleic acid chain IDs |
| `composition.solventChains` | `string[]` | Solvent chain IDs |
| `composition.saccharides` | `string[]` | Sugar chain IDs |
| `composition.nChains` | `number` | Total chain count |
| `composition.nResidues` | `number` | Total residue count |
| `composition.nAtoms` | `number` | Total atom count |
| `composition.ligands` | `LigandInfo[]` | Ligand details |
| `composition.chainDetails` | `ChainDetail[]` | Per-chain composition |

**Key distinction:** `useFileDigest` is reactive (SWR-based, updates when the file changes). Use this for render-time data. `fetchDigest` is imperative (async call) — use this inside `onChange` handlers for one-time metadata extraction like wavelength.

### Multi-Select Chain Selector (`CChainSelectElement`)

For parameters that store a comma-separated list of chain IDs (e.g. `prosmartProtein.CHAINLIST_1`), use the dedicated multi-select widget:

```tsx
import { CChainSelectElement } from "../task-elements/cchainselect";

<CChainSelectElement
  job={job}
  itemName="prosmartProtein.CHAINLIST_1"
  options={composition?.peptides || []}    // Available chains from digest
  label=" "
  visibility={() => isTruthy(prosmartProteinToggle)}
/>
```

This renders an MUI Autocomplete with `multiple` and chip tags. It reads a comma-separated CString (e.g. `"A,B,C"`) and writes back as a comma-joined string. The `options` prop should be populated from the XYZIN digest composition.

### Info/Warning Text

```tsx
<Typography variant="body2" color="warning.main" sx={{ fontStyle: "italic", fontWeight: "bold" }}>
  You should normally let Parrot choose reference structures
</Typography>
```

---

## Pitfalls and Hard-Won Lessons

### 1. CBoolean `__bool__` Trap (Server-Side)

Python's CBoolean class has a `__bool__` method that returns `False` when the boolean value is `False`. This means any `if obj:` check on a CBoolean set to False will skip that code path. **Always use `if obj is not None:` in Python code that handles CData objects.**

This caused a real bug where unchecking a checkbox worked on the server, but the API response omitted the `updated_item` field — so the frontend's `onChange` callback never fired and the UI never toggled.

### 2. Always Normalize CBoolean Values

The server returns `_value: true/false` (JS boolean) for CBoolean items, but some code paths may return the strings `"True"`/`"False"`. Always use:

```tsx
const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";
```

### 3. Use Bound `mutateContainer`, Not Global `mutate`

When patching the SWR cache, use the bound `mutateContainer` from `useSWR` rather than the global `mutate(key, fn, opts)`. Global mutate can miss subscriber notifications in some SWR edge cases.

### 4. onChange Only Fires When `updated_item` Exists

The `onChange` callback in `CSimpleTextFieldElement` only fires when the server response includes `result.data.updated_item`. If the server doesn't return it (e.g. due to the CBoolean truthiness bug above), onChange silently doesn't fire. There's no error — it just doesn't call your handler.

### 5. Spread `{...props}` on Every Element

Every `CCP4i2TaskElement` and `CCP4i2ContainerElement` needs `{...props}` to receive the `job` object. Forgetting this causes "job.id is undefined" errors.

### 6. `qualifiers` Override, They Don't Replace

When you pass `qualifiers={{ guiLabel: "My label" }}`, this is **merged** with the item's own qualifiers from the `.def.xml`. Your values take precedence for any keys you specify, but other qualifiers (like `min`, `max`, `toolTip`) are preserved.

### 7. useTaskItem Returns Undefined Initially

The first render will have `item: undefined`, `value: undefined`. Always guard against this:

```tsx
const { value: MODE } = useTaskItem("MODE");
// MODE may be undefined on first render — visibility functions handle this gracefully
```

### 8. Don't Show Derived Parameters

If a parameter is always computed from another (via `syncTo`), don't render it in the UI. The user can't meaningfully edit it, and showing it creates confusion.

### 9. Include All Dependencies in useCallback/useEffect/useMemo

React hooks rules apply. Missing dependencies cause stale closures and subtle bugs:

```tsx
// BAD: forceSetXYZIN missing from deps
const handleToggle = useCallback(async (item: any) => {
  if (!isTruthy(item._value) && XYZIN?.dbFileId) {
    forceSetXYZIN({});  // stale reference!
  }
}, [XYZIN]);

// GOOD
const handleToggle = useCallback(async (item: any) => {
  if (!isTruthy(item._value) && XYZIN?.dbFileId) {
    forceSetXYZIN({});
  }
}, [XYZIN, forceSetXYZIN]);
```

### 10. Test Boolean Toggles in Both Directions

After implementing CBoolean-driven visibility, always test:
1. Check -> Uncheck (fields should appear/disappear depending on your logic)
2. Uncheck -> Check (must work symmetrically)
3. Page reload while unchecked (state must restore from server)

The most common failure mode is one-way toggling, where only one direction works.

---

## Complete Worked Example — ModelCraft

### The Narrative (What the Scientist Wrote)

> The modelcraft interface has at the top a container with label "Reflection data". That container contains the reflection data widget, the Free R set widget, and a tick box labelled "Get initial phases from refining the starting model (uncheck to specify starting phases, e.g. from experimental phasing)". If the element is unticked, an additional widget is revealed for input phases, together with an additional check box labelled "Phases are unbiased and should be used as refinement restraints when the model is poor".
>
> Below that is a container labelled "Asymmetric unit contents". The only widget in this container is the widget to specify an asymmetric unit content file.
>
> Below that is a container labelled "Starting model". That container contains the widget for providing a starting model from which to build.
>
> Below that is a container labelled "Options". This container contains rows with:
> 1. A tick box with label "Run a quicker basic pipeline"
> 2. Elements which place the text "Run for" in front of a widget to specify the number of cycles, and then the word "cycles"
> 3. Elements that provide a tickbox, followed by the text "Stop automatically if R-free does not improve in", then a number field to specify some sort of convergence criteria, followed by the word "cycles"
> 4. A tick box with label "Build selenomethionine (MSE) instead of methionine (MET)"
> 5. A tick box with label "Use twinned refinement"
>
> Below that is a container labelled "Optional pipeline steps" which contains rows that are all tickboxes with labels:
> 1. "Preliminary low-resolution refinement with Sheetbend"
> 2. "Residue and chain pruning"
> 3. "Classical density modification with Parrot"
> 4. "Phase improvement through addition and refinement of dummy atoms"
> 5. "Addition of waters"
> 6. "Final side-chain fixing"

### How the Narrative Maps to Code

| Narrative phrase | Implementation decision |
|-----------------|----------------------|
| "container with label" | `CCP4i2ContainerElement` with `containerHint="FolderLevel"` and `qualifiers={{ guiLabel: "..." }}` |
| "reflection data widget" | `<CCP4i2TaskElement itemName="F_SIGF" />` — inferred from the task's `.def.xml` |
| "tick box labelled ..." | `<CCP4i2TaskElement itemName="USE_MODEL_PHASES" qualifiers={{ guiLabel: "..." }} />` |
| "If unticked, an additional widget is revealed" | `onChange` + local state + `{!useModelPhases && (...)}` conditional rendering |
| "Run a quicker basic pipeline" (checkbox) | `onChange={handleBASIC}` pushes default CYCLES value (5 or 25) |
| "Run for ___ cycles" | `Box` with flex layout, `Typography` text, constrained-width field with `guiLabel: " "` |
| "tickbox, followed by text, then number field" | Same inline flex pattern with checkbox `sx={{ width: "auto" }}` |
| "rows that are all tickboxes" | Simple `CCP4i2TaskElement` items with `qualifiers={{ guiLabel: "..." }}` — one per line |

### The Resulting Code

This shows the real interface with grouped sections, conditional visibility driven by a CBoolean, inline field-with-text layouts, and reactive default-value pushing.

```tsx
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useCallback, useEffect, useState } from "react";

/** Normalize CBoolean values - server may return boolean or string */
const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { value: USE_MODEL_PHASES_RAW } = useTaskItem("USE_MODEL_PHASES");
  const { value: XYZIN, forceUpdate: forceSetXYZIN } = useTaskItem("XYZIN");
  const { forceUpdate: forceSetCYCLES } = useTaskItem("CYCLES");

  // Local state for immediate UI toggle
  const [useModelPhases, setUseModelPhases] = useState(() =>
    isTruthy(USE_MODEL_PHASES_RAW)
  );

  // Sync from container for programmatic changes (initial load, parameter file import)
  useEffect(() => {
    setUseModelPhases(isTruthy(USE_MODEL_PHASES_RAW));
  }, [USE_MODEL_PHASES_RAW]);

  const handleUSE_MODEL_PHASES = useCallback(
    async (new_USE_MODEL_PHASES: any) => {
      const newValue = isTruthy(new_USE_MODEL_PHASES._value);
      setUseModelPhases(newValue);
      // Clear XYZIN if unchecking and a model file is loaded
      if (!newValue && XYZIN?.dbFileId) {
        forceSetXYZIN({});
      }
    },
    [XYZIN, forceSetXYZIN]
  );

  // Push sensible default for CYCLES when BASIC toggle changes
  const handleBASIC = useCallback(
    async (updatedItem: any) => {
      await forceSetCYCLES(isTruthy(updatedItem._value) ? 5 : 25);
    },
    [forceSetCYCLES]
  );

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Reflection data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Reflection data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
        <CCP4i2TaskElement
          itemName="USE_MODEL_PHASES"
          {...props}
          qualifiers={{
            guiLabel:
              "Get initial phases from refining the starting model (uncheck to specify starting phases, e.g. from experimental phasing)",
          }}
          onChange={handleUSE_MODEL_PHASES}
        />
        {!useModelPhases && (
          <>
            <CCP4i2TaskElement itemName="PHASES" {...props} />
            <CCP4i2TaskElement
              itemName="UNBIASED"
              {...props}
              qualifiers={{
                guiLabel:
                  "Phases are unbiased and should be used as refinement restraints when the model is poor",
              }}
            />
          </>
        )}
      </CCP4i2ContainerElement>

      {/* Asymmetric unit contents */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Asymmetric unit contents" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="ASUIN" {...props} />
      </CCP4i2ContainerElement>

      {/* Starting model */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Starting model" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="BASIC"
          {...props}
          qualifiers={{ guiLabel: "Run a quicker basic pipeline" }}
          onChange={handleBASIC}
        />

        <Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
          <Typography variant="body1">Run for</Typography>
          <Box sx={{ width: "8rem" }}>
            <CCP4i2TaskElement itemName="CYCLES" {...props} qualifiers={{ guiLabel: " " }} />
          </Box>
          <Typography variant="body1">cycles</Typography>
        </Box>

        <CCP4i2TaskElement
          itemName="SELENOMET"
          {...props}
          qualifiers={{
            guiLabel: "Build selenomethionine (MSE) instead of methionine (MET)",
          }}
        />
        <CCP4i2TaskElement
          itemName="TWINNED"
          {...props}
          qualifiers={{ guiLabel: "Use twinned refinement" }}
        />
      </CCP4i2ContainerElement>

      {/* Optional pipeline steps */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Optional pipeline steps" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="SHEETBEND"
          {...props}
          qualifiers={{ guiLabel: "Preliminary low-resolution refinement with Sheetbend" }}
        />
        <CCP4i2TaskElement itemName="PRUNING" {...props} qualifiers={{ guiLabel: "Residue and chain pruning" }} />
        <CCP4i2TaskElement itemName="PARROT" {...props} qualifiers={{ guiLabel: "Classical density modification with Parrot" }} />
        <CCP4i2TaskElement itemName="WATERS" {...props} qualifiers={{ guiLabel: "Addition of waters" }} />
        <CCP4i2TaskElement itemName="SIDE_CHAIN_FIXING" {...props} qualifiers={{ guiLabel: "Final side-chain fixing" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
```

### What this example demonstrates:

| Technique | Where |
|-----------|-------|
| Grouped sections with `FolderLevel` containers | All five `CCP4i2ContainerElement` blocks |
| CBoolean conditional visibility via onChange + local state | `USE_MODEL_PHASES` -> PHASES/UNBIASED |
| `isTruthy()` normalization | `useState` init, `useEffect` sync, `handleUSE_MODEL_PHASES` |
| Clear dependent field on toggle | `forceSetXYZIN({})` when unchecking |
| Push default value on toggle | `BASIC` onChange → `forceSetCYCLES(5 or 25)` |
| Inline text + field layout | "Run for ___ cycles" |
| Qualifier override | Custom `guiLabel` on most fields |
| `useEffect` sync from server | Handles initial load + parameter file import |

---

## Complete Worked Example — ProSMART-Refmac (Multi-Tab, Digest-Driven)

ProSMART-Refmac is the most complex task interface in the codebase. It demonstrates patterns beyond what ModelCraft covers: tabbed layout, file digest-driven composition visibility, sub-container dotted paths, and custom widgets.

### How It Was Built

The interface was built from **screenshots of the legacy Qt GUI** — one per tab (Input data, Parameterisation, Restraints, Output, Advanced) plus a screenshot showing the Restraints tab with dependent elements revealed. Claude cross-referenced the screenshots against `prosmart_refmac.def.xml` and `refmac.def.xml` to map every visible widget to its parameter name and type.

### Key Techniques Demonstrated

| Technique | Where in prosmart_refmac.tsx |
|-----------|----------------------------|
| Tabbed layout | `<CCP4i2Tabs>` with 5 `<CCP4i2Tab>` children |
| Digest-driven visibility | XYZIN digest `composition.peptides`/`nucleics` for ProSMART sections |
| Sub-container dotted paths | `prosmartProtein.TOGGLE`, `prosmartProtein.CHAINLIST_1`, `platonyzer.MODE`, `controlParameters.WEIGHT` |
| Custom multi-select widget | `<CChainSelectElement>` for chain selection |
| File onChange with side effects | F_SIGF onChange → wavelength, anomalous, twinning |
| FreeR warning hook | `useFreeRWarning()` for cross-validation checks |
| Conditional "no data" messages | Italic text when model has no protein/nucleotide chains |
| Nested visibility conditions | Occupancy groups visible only when occupancy refinement is enabled |
| Inline natural-language layouts | "Refine [isotropic] B-factors", "with sigma: [0.01] and max distance: [4.2]" |

### Sub-Container Dotted Paths

When the `.def.xml` defines nested containers (e.g. `prosmartProtein` containing `TOGGLE`, `CHAINLIST_1`, `WEIGHT`), access them with dotted notation:

```tsx
const { value: prosmartProteinToggle } = useTaskItem("prosmartProtein.TOGGLE");

<CCP4i2TaskElement
  {...props}
  itemName="prosmartProtein.WEIGHT"
  qualifiers={{ guiLabel: " " }}
/>
```

### Digest-Driven Section (ProSMART Pattern)

```tsx
// Get XYZIN digest for chain composition
const { item: XYZINItem, value: XYZINValue } = useTaskItem("XYZIN");
const xyzinDigestPath =
  XYZINValue?.dbFileId && XYZINItem?._objectPath ? XYZINItem._objectPath : "";
const { data: xyzinDigest } = useFileDigest(xyzinDigestPath);
const xyzinComposition = xyzinDigest?.composition;
const hasProteinChains =
  Array.isArray(xyzinComposition?.peptides) && xyzinComposition.peptides.length > 0;

// In JSX:
<CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "ProSMART - protein" }} containerHint="FolderLevel">
  {!hasProteinChains && (
    <Typography variant="body2" sx={{ fontStyle: "italic" }}>
      Input atomic model contains no protein chains
    </Typography>
  )}
  {hasProteinChains && (
    <>
      <CChainSelectElement
        job={job}
        itemName="prosmartProtein.CHAINLIST_1"
        options={xyzinComposition?.peptides || []}
        label=" "
      />
      {/* ... more ProSMART protein options */}
    </>
  )}
</CCP4i2ContainerElement>
```

### Multi-Level Inline Layout (Restraints Tab)

The Restraints tab demonstrates complex inline patterns with multiple interleaved text and field elements:

```tsx
{/* Jelly-body: checkbox + conditional inline parameters */}
<CCP4i2TaskElement {...props} itemName="USE_JELLY" qualifiers={{ guiLabel: "Use jelly-body restraints" }} />
{isTruthy(useJelly) && (
  <Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
    <Typography variant="body1">with sigma:</Typography>
    <Box sx={{ width: "8rem" }}>
      <CCP4i2TaskElement {...props} itemName="JELLY_SIGMA" qualifiers={{ guiLabel: " " }} />
    </Box>
    <Typography variant="body1">and max distance:</Typography>
    <Box sx={{ width: "8rem" }}>
      <CCP4i2TaskElement {...props} itemName="JELLY_DIST" qualifiers={{ guiLabel: " " }} />
    </Box>
  </Box>
)}

{/* ProSMART: interleaved text + dropdown + text + field + text */}
<Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
  <Typography variant="body1">Use</Typography>
  <Box sx={{ width: "14rem" }}>
    <CCP4i2TaskElement {...props} itemName="prosmartProtein.ALL_BEST" qualifiers={{ guiLabel: " " }} />
  </Box>
  <Typography variant="body1">chain(s) from reference model(s). Minimum sequence identity:</Typography>
  <Box sx={{ width: "6rem" }}>
    <CCP4i2TaskElement {...props} itemName="prosmartProtein.SEQID" qualifiers={{ guiLabel: " " }} />
  </Box>
  <Typography variant="body1">%</Typography>
</Box>
```

---

## Imports Cheatsheet

```tsx
// Core interface
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

// Tabs
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";

// Layout
import { FieldRow } from "../task-elements/field-row";
import { Box, Grid2, Stack, Paper, Typography, Card, CardHeader, CardContent } from "@mui/material";

// Standalone widgets
import { CChainSelectElement } from "../task-elements/cchainselect";

// Validation hooks
import { useFreeRWarning } from "../../../providers/run-check-provider";

// React hooks
import { useCallback, useEffect, useMemo, useState } from "react";
```
