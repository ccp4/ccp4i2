# Task Interface Layout Guide

This guide explains how to control the size and layout of form fields in CCP4i2 task interfaces.

## Core Design Principle

**Base widgets are full-width by default. Containers control sizing.**

This separation of concerns means:
- Base widgets (TextField, Autocomplete) just render content — no width logic
- Container components (FieldRow, Grid2) are the single source of truth for sizing
- No duplication of width constraints across widgets

## Available Size Constants

Defined in `field-sizes.ts`:

| Size | Value | Use Case |
|------|-------|----------|
| `xs` | 8rem (128px) | Single digits, cell parameters |
| `sm` | 12rem (192px) | Short numbers, small enums |
| `md` | 20rem (320px) | Default for standalone fields |
| `lg` | 32rem (512px) | File paths, longer text |
| `full` | 100% | Full container width |

---

## Layout Components

### 1. FieldRow — Equal-Width Horizontal Layout

Use `FieldRow` when you want multiple fields side-by-side with equal widths.

```tsx
import { FieldRow } from "../task-elements/field-row";

// Two fields sharing space equally (default behavior)
<FieldRow>
  <CCP4i2TaskElement itemName="CRYSTALNAME" qualifiers={{ guiLabel: "Crystal name" }} />
  <CCP4i2TaskElement itemName="DATASETNAME" qualifiers={{ guiLabel: "Dataset name" }} />
</FieldRow>
```

**Result:** Each field gets 50% width, wraps on narrow screens.

### 2. FieldRow with Size — Constrained Widths

Use `FieldRow` with `equalWidth={false}` and `size` to constrain each field to a specific width.

```tsx
// Cell parameters - each field constrained to 8rem
<FieldRow equalWidth={false} size="xs">
  <CCP4i2TaskElement itemName="CELL_A" qualifiers={{ guiLabel: "a" }} />
  <CCP4i2TaskElement itemName="CELL_B" qualifiers={{ guiLabel: "b" }} />
  <CCP4i2TaskElement itemName="CELL_C" qualifiers={{ guiLabel: "c" }} />
</FieldRow>
```

**Result:** Each field is max 8rem wide, fields wrap when container is narrow.

### 3. Grid2 — Precise Column Control

Use MUI's `Grid2` when you need specific column proportions or responsive breakpoints.

```tsx
import { Grid2 } from "@mui/material";

<Grid2 container spacing={2}>
  <Grid2 size={{ xs: 12, md: 4 }}>
    <CCP4i2TaskElement itemName="SPACEGROUP" qualifiers={{ guiLabel: "Space group" }} />
  </Grid2>
  <Grid2 size={{ xs: 12, md: 8 }}>
    <CCP4i2TaskElement itemName="UNITCELL" />
  </Grid2>
</Grid2>
```

**Result:** On mobile (xs), each field is full-width. On desktop (md+), space group is 1/3, unit cell is 2/3.

### 4. Stack — For Non-Field Content Only

Use `Stack direction="row"` **only** for display elements (Typography, Chip, Icon, Button).

```tsx
// CORRECT: Stack for display elements
<Stack direction="row" spacing={1} alignItems="center">
  <Typography variant="body2">Resolution:</Typography>
  <Chip label="2.5 Å" size="small" />
</Stack>

// WRONG: Don't use Stack for form fields
<Stack direction="row">
  <CCP4i2TaskElement itemName="FIELD1" />  {/* Fields won't size properly! */}
  <CCP4i2TaskElement itemName="FIELD2" />
</Stack>
```

---

## Common Patterns

### Pattern 1: Label + Field + Help Text

```tsx
<Box>
  <Typography variant="subtitle2" gutterBottom>
    Resolution Range
  </Typography>
  <FieldRow equalWidth={false} size="sm">
    <CCP4i2TaskElement itemName="RES_LOW" qualifiers={{ guiLabel: "Low" }} />
    <CCP4i2TaskElement itemName="RES_HIGH" qualifiers={{ guiLabel: "High" }} />
  </FieldRow>
  <Typography variant="caption" color="text.secondary">
    Enter resolution limits in Ångstroms
  </Typography>
</Box>
```

### Pattern 2: Mixed Layout (Some Fields Wide, Some Narrow)

```tsx
<Grid2 container spacing={2}>
  {/* Full-width file selector */}
  <Grid2 size={{ xs: 12 }}>
    <CCP4i2TaskElement itemName="INPUT_FILE" />
  </Grid2>

  {/* Three narrow fields in a row */}
  <Grid2 size={{ xs: 12 }}>
    <FieldRow equalWidth={false} size="sm">
      <CCP4i2TaskElement itemName="PARAM1" qualifiers={{ guiLabel: "Param 1" }} />
      <CCP4i2TaskElement itemName="PARAM2" qualifiers={{ guiLabel: "Param 2" }} />
      <CCP4i2TaskElement itemName="PARAM3" qualifiers={{ guiLabel: "Param 3" }} />
    </FieldRow>
  </Grid2>

  {/* Two equal-width fields */}
  <Grid2 size={{ xs: 12 }}>
    <FieldRow>
      <CCP4i2TaskElement itemName="OPTION_A" qualifiers={{ guiLabel: "Option A" }} />
      <CCP4i2TaskElement itemName="OPTION_B" qualifiers={{ guiLabel: "Option B" }} />
    </FieldRow>
  </Grid2>
</Grid2>
```

### Pattern 3: Responsive 2-Column Form

```tsx
<Grid2 container spacing={2}>
  <Grid2 size={{ xs: 12, sm: 6 }}>
    <CCP4i2TaskElement itemName="CRYSTAL_NAME" qualifiers={{ guiLabel: "Crystal" }} />
  </Grid2>
  <Grid2 size={{ xs: 12, sm: 6 }}>
    <CCP4i2TaskElement itemName="DATASET_NAME" qualifiers={{ guiLabel: "Dataset" }} />
  </Grid2>
  <Grid2 size={{ xs: 12, sm: 6 }}>
    <CCP4i2TaskElement itemName="WAVELENGTH" qualifiers={{ guiLabel: "Wavelength" }} />
  </Grid2>
  <Grid2 size={{ xs: 12, sm: 6 }}>
    <CCP4i2TaskElement itemName="RESOLUTION" qualifiers={{ guiLabel: "Resolution" }} />
  </Grid2>
</Grid2>
```

### Pattern 4: Field with Inline Suffix

```tsx
<Stack direction="row" spacing={1} alignItems="center">
  <Box sx={{ width: "12rem" }}>
    <CCP4i2TaskElement itemName="WAVELENGTH" qualifiers={{ guiLabel: "Wavelength" }} />
  </Box>
  <Typography variant="body2" color="text.secondary">Å</Typography>
</Stack>
```

### Pattern 5: Card with Grouped Fields

```tsx
import { Card, CardHeader, CardContent } from "@mui/material";

<Card>
  <CardHeader title="Crystal Information" />
  <CardContent>
    <Grid2 container spacing={2}>
      <Grid2 size={{ xs: 6, md: 4 }}>
        <CCP4i2TaskElement itemName="SPACEGROUP" qualifiers={{ guiLabel: "Space group" }} />
      </Grid2>
      <Grid2 size={{ xs: 6, md: 8 }}>
        <CCP4i2TaskElement itemName="UNITCELL" />
      </Grid2>
      <Grid2 size={{ xs: 12 }}>
        <FieldRow>
          <CCP4i2TaskElement itemName="CRYSTALNAME" qualifiers={{ guiLabel: "Crystal" }} />
          <CCP4i2TaskElement itemName="DATASETNAME" qualifiers={{ guiLabel: "Dataset" }} />
        </FieldRow>
      </Grid2>
    </Grid2>
  </CardContent>
</Card>
```

---

## Composite Widgets (Pre-Built Layouts)

Some data types have dedicated composite widgets with built-in layout:

| Data Type | Widget | Layout |
|-----------|--------|--------|
| `CCell` | `CCellElement` | 6 fields (a,b,c,α,β,γ) in a row, `xs` size each |
| `CReindexOperator` | `CReindexOperatorElement` | Matrix values in a row, `xs` size each |
| `CEnsemble` | `CEnsembleElement` | Grid layout with copies, label, use fields |
| `CImportUnmerged` | `CImportUnmergedElement` | File + metadata grid |

These are automatically used when rendering the corresponding data types.

---

## Quick Reference

| Want to... | Use |
|------------|-----|
| Two fields side-by-side, equal width | `<FieldRow>` |
| Multiple small fields in a row | `<FieldRow equalWidth={false} size="xs">` |
| Specific column proportions | `<Grid2 size={{ xs: 4 }}>` |
| Responsive breakpoints | `<Grid2 size={{ xs: 12, md: 6 }}>` |
| Display text + icons in a row | `<Stack direction="row">` |
| Full-width field | Just use `<CCP4i2TaskElement>` (default is full-width) |

---

## What NOT to Do

```tsx
// DON'T: Use Stack for form fields
<Stack direction="row">
  <CCP4i2TaskElement itemName="FIELD1" />
</Stack>

// DON'T: Add minWidth directly to fields
<CCP4i2TaskElement sx={{ minWidth: "10rem" }} itemName="FIELD1" />

// DON'T: Use elementSx on containers (deprecated)
<CCP4i2ContainerElement elementSx={{ width: "8rem" }} />
```

Instead:
```tsx
// DO: Use FieldRow to control width
<FieldRow equalWidth={false} size="sm">
  <CCP4i2TaskElement itemName="FIELD1" />
</FieldRow>

// DO: Use Grid2 for precise control
<Grid2 size={{ xs: 4 }}>
  <CCP4i2TaskElement itemName="FIELD1" />
</Grid2>
```

---

## Imports

```tsx
// Layout components
import { FieldRow } from "../task-elements/field-row";
import { Grid2, Stack, Box, Card, CardHeader, CardContent, Typography } from "@mui/material";

// Task elements
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
```

---

## Visibility Control

Fields can be shown or hidden based on the values of other parameters using the `visibility` prop.

### Basic Visibility

Pass a function that returns `true` (show) or `false` (hide):

```tsx
const { value: INPUT_FIXEDValue } = useTaskItem("INPUT_FIXED");

<CCP4i2TaskElement
  {...props}
  itemName="XYZIN_FIXED"
  qualifiers={{ guiLabel: "Known partial model" }}
  visibility={() => INPUT_FIXEDValue === true}
/>
```

### Visibility Based on Enum Values

```tsx
const { value: COMP_BYValue } = useTaskItem("COMP_BY");

// Show ASU file selector only when "ASU" is selected
<CCP4i2TaskElement
  {...props}
  itemName="ASUFILE"
  qualifiers={{ guiLabel: "CCP4i2 ASU file" }}
  visibility={() => COMP_BYValue === "ASU"}
/>

// Show molecular weight fields only when "MW" is selected
<FieldRow>
  <CCP4i2TaskElement
    {...props}
    itemName="ASU_NUCLEICACID_MW"
    qualifiers={{ guiLabel: "nucleic acid (Da)" }}
    visibility={() => COMP_BYValue === "MW"}
  />
  <CCP4i2TaskElement
    {...props}
    itemName="ASU_PROTEIN_MW"
    qualifiers={{ guiLabel: "protein (Da)" }}
    visibility={() => COMP_BYValue === "MW"}
  />
</FieldRow>
```

### Visibility Based on File Content

```tsx
const { value: F_SIGFValue } = useTaskItem("F_SIGF");

// Show F/I selector only when file contains both F and I data
<CCP4i2TaskElement
  {...props}
  itemName="F_OR_I"
  qualifiers={{ guiLabel: "Use Fs or Is" }}
  visibility={() => [1, 3].includes(F_SIGFValue?.contentFlag)}
/>
```

### Organizing Visibility Functions

For complex interfaces, group visibility conditions in a `useMemo`:

```tsx
const { value: partialModeOrMap } = useTaskItem("PARTIALMODELORMAP");
const { value: compBy } = useTaskItem("COMP_BY");
const { value: runModelCraft } = useTaskItem("RUNMODELCRAFT");

const visibility = useMemo(() => ({
  isPartialModel: () => partialModeOrMap === "MODEL",
  isNoPartial: () => partialModeOrMap === "NONE",
  isSearch: () => partialModeOrMap === "SEARCH",
  isAsuFile: () => compBy === "ASU",
  isMolecularWeight: () => compBy === "MW",
  hasModelCraft: () => runModelCraft === true,
}), [partialModeOrMap, compBy, runModelCraft]);

// Usage
<CCP4i2TaskElement
  {...props}
  itemName="XYZIN_PARTIAL"
  qualifiers={{ guiLabel: "Partial model coordinates" }}
  visibility={visibility.isPartialModel}
/>
```

---

## Reactive Interfaces — Responding to Value Changes

Task interfaces can respond to parameter changes to auto-populate related fields.

### Reading Task Values

Use `useTaskItem` to access current values:

```tsx
const { useTaskItem } = useJob(job.id);

// Read-only access to value
const { value: wavelengthValue } = useTaskItem("WAVELENGTH");

// Access to both value and item metadata
const { item: F_SIGFItem, value: F_SIGFValue } = useTaskItem("F_SIGF");
```

### Updating Task Values

Use `forceUpdate` to programmatically update parameters:

```tsx
const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");
const { forceUpdate: forceUpdateSPACEGROUP } = useTaskItem("SPACEGROUP");

// Update wavelength
await forceUpdateWAVELENGTH(1.54);

// Update space group
await forceUpdateSPACEGROUP("P212121");
```

### Auto-Syncing Parameter Values

Sometimes a parameter should be automatically derived from another parameter's value, removing the need for the user to set it manually. Use `useEffect` to sync values and hide the derived control from the UI.

**Example: Auto-set a toggle based on mode selection**

In the aimless_pipe interface, when the user selects "MATCH" mode, a reference is always required. Rather than showing a toggle asking "Provide reference?", we auto-set it based on the mode:

```tsx
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Get update function for the derived parameter
  const { update: updateAimlessRef } = useTaskItem("REFERENCE_FOR_AIMLESS");

  // Get values needed to compute the derived value
  const { value: modeValue } = useTaskItem("MODE");
  const { value: aimlessRefValue } = useTaskItem("REFERENCE_FOR_AIMLESS");
  const { value: referenceDataset } = useTaskItem("REFERENCE_DATASET");

  // Track pending updates to prevent bouncing loops
  const pendingAimlessRef = useRef<boolean | null>(null);

  // Auto-sync REFERENCE_FOR_AIMLESS based on MODE:
  // - true when MODE is "MATCH" (always provide reference when matching)
  // - false otherwise
  useEffect(() => {
    const isMatchMode = modeValue === "MATCH";
    const currentRef = aimlessRefValue;
    const targetRef = isMatchMode;

    // Check if our pending update has been applied
    if (pendingAimlessRef.current !== null) {
      if (pendingAimlessRef.current === currentRef) {
        pendingAimlessRef.current = null; // Update applied
      } else {
        return; // Still waiting, don't trigger another update
      }
    }

    // Only update if value needs to change
    if (targetRef !== currentRef) {
      pendingAimlessRef.current = targetRef;
      updateAimlessRef(targetRef);
    }
  }, [modeValue, aimlessRefValue, updateAimlessRef]);

  // Visibility - since REFERENCE_FOR_AIMLESS is auto-synced, we don't show it
  // but downstream fields that depend on it can assume it's true in MATCH mode
  const visibility = useMemo(() => ({
    isMatchMode: () => modeValue === "MATCH",
    // No need to check aimlessRefValue - it's always true when mode is "MATCH"
    isHklReference: () => modeValue === "MATCH" && referenceDataset === "HKL",
  }), [modeValue, referenceDataset]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Options">
        {/* User selects mode */}
        <CCP4i2TaskElement
          {...props}
          itemName="MODE"
          qualifiers={{ guiLabel: "Pipeline mode" }}
        />

        {/* REFERENCE_FOR_AIMLESS is NOT shown - it's auto-synced with MODE */}

        {/* Reference type selector - shown when mode is MATCH */}
        <CCP4i2TaskElement
          {...props}
          itemName="REFERENCE_DATASET"
          qualifiers={{ guiLabel: "Reference type" }}
          visibility={visibility.isMatchMode}
        />

        {/* Reference file - shown based on reference type */}
        <CCP4i2TaskElement
          {...props}
          itemName="HKLIN_REF"
          qualifiers={{ guiLabel: "Reference reflections" }}
          visibility={visibility.isHklReference}
        />
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
```

**Key principles for auto-syncing:**

1. **Check before updating** - Always compare current value before calling update to avoid infinite loops
2. **Remove from UI** - Don't show the derived parameter to the user; it's computed, not chosen
3. **Simplify downstream visibility** - Once a parameter is auto-synced, downstream visibility checks can be simplified (no need to check the derived value)
4. **Use `syncTo` for auto-sync** - The `syncTo` method handles bounce prevention internally

### Using `syncTo` for Auto-Sync (Recommended)

The `syncTo` method is specifically designed for auto-syncing one parameter based on another. It handles bounce prevention internally, making your code much simpler:

```tsx
// Get syncTo for the parameter you want to auto-set
const { syncTo: syncToAimlessRef } = useTaskItem("REFERENCE_FOR_AIMLESS");

// Get the driving value
const { value: modeValue } = useTaskItem("MODE");

// Auto-sync: REFERENCE_FOR_AIMLESS = true when MODE is "MATCH", false otherwise
useEffect(() => {
  syncToAimlessRef(modeValue === "MATCH");
}, [modeValue, syncToAimlessRef]);
```

That's it! The `syncTo` method:
- Tracks pending updates internally to prevent bouncing
- Only triggers an update if the value actually needs to change
- Handles all the edge cases around re-renders and stale values

### How `syncTo` Prevents Bouncing (Under the Hood)

When `update` is called, the container re-renders before the new value has fully propagated. This can cause the useEffect to fire again with stale values, creating an infinite loop. The `syncTo` method prevents this by:

1. Storing the target value in a module-level Map before triggering the update
2. On subsequent calls (while waiting for propagation):
   - If the stored value matches current value → update complete, clear pending state
   - If stored value differs from current → still propagating, skip this call
3. Once pending state is cleared, normal comparison determines if a new update is needed

**Manual implementation (for reference only):**
```tsx
// You don't need to write this - use syncTo instead!
const pendingRef = useRef<boolean | null>(null);

useEffect(() => {
  const targetRef = modeValue === "MATCH";

  // Check if pending update has been applied
  if (pendingRef.current !== null) {
    if (pendingRef.current === currentRef) {
      pendingRef.current = null;  // Update applied
    } else {
      return;  // Still waiting
    }
  }

  if (targetRef !== currentRef) {
    pendingRef.current = targetRef;
    updateAimlessRef(targetRef);
  }
}, [modeValue, currentRef, updateAimlessRef]);
```

### Responding to File Selection — Extract Metadata

When a user selects a file, extract metadata (wavelength, cell, spacegroup) from the file digest:

```tsx
const { useTaskItem, fetchDigest } = useJob(job.id);

const { item: F_SIGFItem } = useTaskItem("F_SIGF");
const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");

// Handler called when file changes
const handleFileChange = useCallback(async () => {
  if (!F_SIGFItem?._objectPath) return;

  // Fetch file digest (contains metadata extracted from file)
  const digestData = await fetchDigest(F_SIGFItem._objectPath);

  // Extract wavelength from digest
  const wavelength = digestData?.wavelengths?.at(-1);
  if (wavelength && wavelength > 0 && wavelength < 9) {
    await forceUpdateWAVELENGTH(wavelength);
  }
}, [F_SIGFItem?._objectPath, fetchDigest, forceUpdateWAVELENGTH]);

// Pass handler to the file element
<CCP4i2TaskElement
  {...props}
  itemName="F_SIGF"
  qualifiers={{ guiLabel: "Reflections" }}
  onChange={handleFileChange}
/>
```

### Full Example: Auto-populate from Reflection File

```tsx
import { useCallback } from "react";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest, mutateValidation } = useJob(job.id);

  // Items for reading/writing
  const { item: HKLINItem } = useTaskItem("HKLIN");
  const { forceUpdate: forceUpdateSPACEGROUP } = useTaskItem("SPACEGROUP");
  const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");
  const { forceUpdate: forceUpdateUNITCELL } = useTaskItem("UNITCELL");

  // Handle reflection file change
  const handleHKLINChange = useCallback(async () => {
    if (!HKLINItem?._objectPath) return;

    const digest = await fetchDigest(HKLINItem._objectPath);
    if (!digest) return;

    // Update space group
    if (digest.spaceGroup) {
      await forceUpdateSPACEGROUP(digest.spaceGroup.replace(/\s+/g, ""));
    }

    // Update wavelength
    if (digest.wavelength && digest.wavelength > 0) {
      await forceUpdateWAVELENGTH(digest.wavelength);
    }

    // Update unit cell (object with a, b, c, alpha, beta, gamma)
    if (digest.cell) {
      await forceUpdateUNITCELL(digest.cell);
    }

    // Refresh validation after updates
    await mutateValidation();
  }, [
    HKLINItem?._objectPath,
    fetchDigest,
    forceUpdateSPACEGROUP,
    forceUpdateWAVELENGTH,
    forceUpdateUNITCELL,
    mutateValidation,
  ]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Input">
        <CCP4i2TaskElement
          {...props}
          itemName="HKLIN"
          qualifiers={{ guiLabel: "Reflections" }}
          onChange={handleHKLINChange}
        />
        <CCP4i2TaskElement
          {...props}
          itemName="SPACEGROUP"
          qualifiers={{ guiLabel: "Space group" }}
        />
        <CCP4i2TaskElement {...props} itemName="UNITCELL" />
        <CCP4i2TaskElement
          {...props}
          itemName="WAVELENGTH"
          qualifiers={{ guiLabel: "Wavelength" }}
        />
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
```

### Available Digest Fields

When you call `fetchDigest()`, the returned object may contain:

| Field | Type | Description |
|-------|------|-------------|
| `spaceGroup` | `string` | Space group name (e.g., "P212121") |
| `cell` | `{ a, b, c, alpha, beta, gamma }` | Unit cell parameters |
| `wavelength` | `number` | Single wavelength value |
| `wavelengths` | `number[]` | Array of wavelengths (multi-wavelength data) |
| `resolution` | `{ low, high }` | Resolution range |
| `contentFlag` | `number` | Data type: 1=anom I, 2=anom F, 3=mean I, 4=mean F |
| `hasFreeR` | `boolean` | Whether file contains FreeR flags |
| `freerValid` | `boolean` | Whether FreeR flags are valid |
| `crystalNames` | `string[]` | Crystal names from file |
| `datasets` | `string[]` | Dataset names from file |

### When to Use forceUpdate vs update

| Function | Use Case |
|----------|----------|
| `forceUpdate` | Programmatic updates from code (file change handlers, digest processing) |
| `update` | Less common; used when you need the raw update function |

Both functions:
- Update the server
- Patch the local SWR cache
- Trigger validation refresh

---

## Combining Visibility and Reactivity

A complete example showing both patterns together:

```tsx
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest } = useJob(job.id);

  // Values for visibility
  const { value: obsTypeValue } = useTaskItem("OBS_TYPE");
  const { value: F_SIGFValue } = useTaskItem("F_SIGF");

  // Update functions for reactivity
  const { item: F_SIGFItem } = useTaskItem("F_SIGF");
  const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");

  // Reactive handler
  const handleFileChange = useCallback(async () => {
    if (!F_SIGFItem?._objectPath) return;
    const digest = await fetchDigest(F_SIGFItem._objectPath);
    if (digest?.wavelength) {
      await forceUpdateWAVELENGTH(digest.wavelength);
    }
  }, [F_SIGFItem?._objectPath, fetchDigest, forceUpdateWAVELENGTH]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Data">
        {/* File selector with onChange handler */}
        <CCP4i2TaskElement
          {...props}
          itemName="F_SIGF"
          qualifiers={{ guiLabel: "Reflections" }}
          onChange={handleFileChange}
        />

        {/* Conditionally visible based on file content */}
        <CCP4i2TaskElement
          {...props}
          itemName="F_OR_I"
          qualifiers={{ guiLabel: "Use Fs or Is" }}
          visibility={() => [1, 3].includes(F_SIGFValue?.contentFlag)}
        />

        {/* Conditionally visible based on enum selection */}
        <CCP4i2TaskElement
          {...props}
          itemName="ANOM_WAVELENGTH"
          qualifiers={{ guiLabel: "Anomalous wavelength" }}
          visibility={() => obsTypeValue === "ANOMALOUS"}
        />

        {/* Always visible, auto-populated from file */}
        <CCP4i2TaskElement
          {...props}
          itemName="WAVELENGTH"
          qualifiers={{ guiLabel: "Wavelength" }}
        />
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
```

---

## Decision Guide: Which Pattern to Use?

| Scenario | Pattern | Example |
|----------|---------|---------|
| Show/hide based on another parameter | **Visibility** | Show "anomalous wavelength" only when obs type is "ANOMALOUS" |
| Show/hide based on file content | **Visibility** | Show "Use Fs or Is" only when file contains both |
| Auto-populate from file metadata | **onChange + fetchDigest** | Extract wavelength when reflection file is selected |
| Parameter value determined by another | **useEffect auto-sync** | Toggle always true when mode is "MATCH" |
| Derived parameter user shouldn't control | **Auto-sync + hide from UI** | Don't show toggle, just sync it with mode |
| Cascade of dependent fields | **Combine visibility + auto-sync** | Mode → auto-sync toggle → visibility for downstream fields |
| Preview dependent option before enabling | **Disabled** | Show "Reference type" greyed out until "MATCH" mode selected |
| Related options side-by-side | **FieldRow + disabled** | MODE and REFERENCE_TYPE on same row, type disabled until MATCH |

### Pattern Comparison

**Visibility (show/hide)**
```tsx
// User sees field only when condition is met
<CCP4i2TaskElement
  itemName="SPECIAL_OPTION"
  visibility={() => modeValue === "ADVANCED"}
/>
```
- User still controls the value when visible
- Field state persists when hidden

**Auto-sync (computed value)**
```tsx
// Track pending to prevent bouncing
const pendingRef = useRef<boolean | null>(null);

useEffect(() => {
  const targetRef = modeValue === "MATCH";

  // Wait for pending update to propagate
  if (pendingRef.current !== null) {
    if (pendingRef.current === refValue) pendingRef.current = null;
    else return;
  }

  if (targetRef !== refValue) {
    pendingRef.current = targetRef;
    updateRef(targetRef);
  }
}, [modeValue, refValue, updateRef]);

// Don't render the field at all - it's derived
```
- Value is always derived from other state
- Remove from UI entirely
- **Must use ref to prevent bouncing loops**

**onChange handler (reactive to user action)**
```tsx
// Respond to explicit user action (file selection)
<CCP4i2TaskElement
  itemName="HKLIN"
  onChange={async () => {
    const digest = await fetchDigest(item._objectPath);
    if (digest?.wavelength) {
      await forceUpdateWAVELENGTH(digest.wavelength);
    }
  }}
/>
```
- Triggered by user interaction
- Use for extracting metadata from files

**Disabled state (visible but inactive)**
```tsx
// Field is always visible but disabled based on condition
const isReferenceTypeDisabled = useMemo(
  () => () => modeValue !== "MATCH",
  [modeValue]
);

<CCP4i2TaskElement
  itemName="REFERENCE_TYPE"
  qualifiers={{ guiLabel: "Reference type" }}
  disabled={isReferenceTypeDisabled}
/>
```
- User can see what options exist before enabling them
- Use when you want to preview dependent options
- Better UX than hiding when the relationship should be obvious

### When to Use Disabled vs Visibility

| Scenario | Use | Reason |
|----------|-----|--------|
| Related options on same row | **disabled** | Shows relationship, user can preview |
| Completely different mode/section | **visibility** | Reduces clutter, unrelated to current mode |
| Cascading selections | **disabled** for immediate children, **visibility** for deeper | Balance between preview and clutter |

**Example: Side-by-side with disabled state**

Place a mode selector next to its dependent field, with the dependent field disabled until the mode enables it:

```tsx
// Disabled helper - reference type is visible but disabled unless mode is MATCH
const isReferenceTypeDisabled = useMemo(
  () => () => taskValues.mode !== "MATCH",
  [taskValues.mode]
);

// Side-by-side layout: MODE always editable, REFERENCE_DATASET disabled unless MATCH
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
    disabled={isReferenceTypeDisabled}
  />
</FieldRow>

{/* File selector only shown when mode is MATCH and type is selected */}
<CCP4i2TaskElement
  {...props}
  itemName="HKLIN_REF"
  qualifiers={{ guiLabel: "Reference reflections" }}
  visibility={() => taskValues.mode === "MATCH" && taskValues.referenceDataset === "HKL"}
/>
```

This pattern:
- Shows users what options exist when they select "MATCH" mode
- Keeps the UI compact with side-by-side layout
- Uses visibility for the file selector since it depends on two conditions

### Common Mistakes to Avoid

1. **Showing derived parameters** - If a value is always computed from another, hide it from the UI
2. **Missing change guards in useEffect** - Always check if value needs to change before calling update
3. **Using visibility for derived values** - If a parameter should always have a specific value based on state, use auto-sync, not visibility
4. **Forgetting dependencies** - Include all referenced values in useEffect/useMemo dependency arrays
5. **Using visibility when disabled is better** - If the dependent field is closely related and should preview, use disabled instead
6. **Not using `syncTo` for auto-sync effects** - Always use `syncTo` instead of `update` when auto-syncing parameters to prevent bouncing loops (see "Using `syncTo` for Auto-Sync" section)
