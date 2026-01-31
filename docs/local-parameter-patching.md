# Local Parameter Patching Design

## Overview

This document describes a refactoring of the parameter update architecture in CCP4i2's Django/React frontend. The goal is to replace the current "full container invalidation" approach with local patching, eliminating timing-based heuristics and simplifying the developer experience for task interface authors.

## Problem Statement

### Current Architecture

When a parameter is updated:

1. Frontend sends parameter change to server via `POST /api/jobs/{id}/set_parameter/`
2. Server loads plugin from `.def.xml`, overlays `input_params.xml`, applies new value, saves
3. Frontend invalidates entire container cache via `mutateContainer()`
4. Full container is refetched from server
5. Each field checks `wasRecentlyChanged(path, 10000)` to decide whether to accept server value

### Issues

| Problem | Manifestation |
|---------|---------------|
| **Race conditions** | Server response may not reflect other pending changes |
| **Arbitrary timing** | Intent system uses magic numbers (10s, 30s) with no principled basis |
| **Infinite loops** | Derived updates trigger container refetch, which triggers effects, which trigger updates |
| **Complex guards** | Each task interface invents its own de-duplication logic |
| **Poor developer experience** | Authors must understand intent system, effect dependencies, `forceUpdate` vs `update` vs `updateNoMutate` |

### Example: Current Complexity in import_merged.tsx

```typescript
// Manual de-duplication tracking
const [processedDigestKey, setProcessedDigestKey] = useState<string | null>(null);

useEffect(() => {
  const digestKey = HKLINDigest ? JSON.stringify({
    spaceGroup: HKLINDigest.spaceGroup,
    wavelength: HKLINDigest.wavelength,
    // ...
  }) : null;

  if (!digestKey || digestKey === processedDigestKey || job?.status !== 1) {
    return;
  }

  // Multiple forceUpdate calls, each potentially triggering refetch
  await forceUpdateSPACEGROUP(cleanedSG);
  await forceUpdateWAVELENGTH(HKLINDigest.wavelength);
  // ...

  setProcessedDigestKey(digestKey);
  // eslint-disable-next-line react-hooks/exhaustive-deps  <-- suppressing lint!
}, [...]);
```

## Solution: Phased Refactoring

### Phase 1: Local Patching with SWR Cache Mutation

**Goal**: Eliminate full container invalidation and the intent system.

#### API Changes

Modify `set_parameter` response to return the updated item with its full path:

```typescript
// Current response
{
  success: true,
  data: {
    path: "import_merged.container.inputData.HKLIN",
    value: "/path/to/file.mtz",
    object_type: "CDataFile",
    file_path: "/path/to/file.mtz"
  }
}

// New response
{
  success: true,
  data: {
    path: "import_merged.container.inputData.HKLIN",
    updated_item: {
      _objectPath: "import_merged.container.inputData.HKLIN",
      _baseClass: "CDataFile",
      _value: { ... },  // Full item structure as returned in container
      // ... all item properties
    },
    version: 43  // Optional: for future sync verification
  }
}
```

#### Frontend Changes

**1. Patch SWR cache instead of invalidating**

```typescript
// In useJob.setParameter
const setParameter = useCallback(async (arg: SetParameterArg) => {
  const result = await api.post<SetParameterResponse>(
    `jobs/${job.id}/set_parameter`,
    arg
  );

  if (result.success) {
    // Patch the container cache directly
    mutateContainer(
      (currentContainer) => patchContainer(currentContainer, result.data.updated_item),
      { revalidate: false }  // Don't refetch
    );
  }

  return result;
}, [...]);
```

**2. Remove intent system entirely**

Delete:
- `parameter-change-intent-provider.tsx`
- All `wasRecentlyChanged` checks in field components
- `setIntent`, `clearIntent`, `setIntentForPath`, `clearIntentForPath` calls

**3. Add event-driven sync**

Refetch full container at meaningful decision points, not arbitrary intervals:

```typescript
// Sync on tab/window focus
useEffect(() => {
  const handleFocus = () => {
    mutateContainer();  // Full refetch
  };
  window.addEventListener('focus', handleFocus);
  return () => window.removeEventListener('focus', handleFocus);
}, [mutateContainer]);

// Sync before job submission
const handleRun = async () => {
  await mutateContainer();  // Ensure we have latest state
  // Proceed with run
};
```

**4. Simplify field components**

```typescript
// csimple-textfield.tsx - BEFORE
useEffect(() => {
  if (objectPath && wasRecentlyChanged(objectPath, 10000)) {
    return;  // Skip sync
  }
  setValue(initialValue);
}, [initialValue, objectPath, wasRecentlyChanged]);

// csimple-textfield.tsx - AFTER
useEffect(() => {
  setValue(initialValue);
}, [initialValue]);
```

#### Container Patching Implementation

```typescript
// utils/container-patch.ts

/**
 * Patches a container with an updated item, using the lookup table for O(1) access.
 */
export function patchContainer(
  container: ContainerWithLookup,
  updatedItem: ContainerItem
): ContainerWithLookup {
  const path = updatedItem._objectPath;

  // Deep clone to avoid mutating cached data
  const newContainer = structuredClone(container);

  // Update the lookup table
  newContainer.lookup[path] = updatedItem;

  // Update the nested structure (walk the path)
  setNestedItem(newContainer, path, updatedItem);

  return newContainer;
}

/**
 * Sets an item at a nested path within the container structure.
 */
function setNestedItem(
  container: any,
  path: string,
  item: ContainerItem
): void {
  const segments = path.split('.');
  let current = container;

  // Navigate to parent, skipping first segment (task name) and last (item name)
  for (let i = 1; i < segments.length - 1; i++) {
    const segment = segments[i];
    if (segment === 'container') continue;  // Skip .container. segment

    if (current._value && typeof current._value === 'object') {
      current = current._value[segment];
    } else if (current[segment]) {
      current = current[segment];
    }
  }

  // Set the item at the final segment
  const finalSegment = segments[segments.length - 1];
  if (current._value && typeof current._value === 'object') {
    current._value[finalSegment] = item;
  }
}
```

### Phase 2: Server-Side Validity

**Goal**: Move validation logic from frontend effects to backend plugin methods.

#### Current State

Validity logic is scattered across task interfaces:

```typescript
// aimless_pipe.tsx
// TODO: This filtering should really be done in Python's validity() method
// to avoid race conditions between frontend filtering and server validation
const filteredErrors = useMemo(() => {
  return errors.filter(e => !shouldSuppressError(e));
}, [errors]);
```

#### Target State

Plugins define validity logic in Python:

```python
# wrappers/aimless_pipe/script/aimless_pipe.py
class aimless_pipe(CPluginScript):

    def validity(self, container):
        errors = []

        # Example: HKLIN required
        if not container.inputData.HKLIN.isSet():
            errors.append(CValidityError(
                path='inputData.HKLIN',
                severity=2,  # Error
                message='Reflection data file is required'
            ))

        # Example: Conditional requirement
        if container.controlParameters.SCALING_MODE.value == 'ANOMALOUS':
            if not container.inputData.WAVELENGTH.isSet():
                errors.append(CValidityError(
                    path='inputData.WAVELENGTH',
                    severity=1,  # Warning
                    message='Wavelength recommended for anomalous scaling'
                ))

        return errors
```

#### API Changes

Include validation in `set_parameter` response:

```typescript
{
  success: true,
  data: {
    path: "...",
    updated_item: { ... },
    validation: {
      "inputData.HKLIN": { maxSeverity: 0 },
      "inputData.XYZIN": { maxSeverity: 2, message: "Required" },
      // ... validation state for all parameters
    }
  }
}
```

#### Migration Strategy

1. For each task interface with frontend validity logic:
   - Identify the conditions being checked
   - Translate to equivalent Python in the plugin's `validity()` method
   - Remove frontend filtering/suppression logic

2. Frontend validation display becomes purely presentational:
   - Receives validation state from server
   - Renders appropriate indicators
   - No conditional logic

### Phase 3: Backend-Driven Extraction

**Goal**: Eliminate complex frontend logic for derived parameter updates (e.g., extracting metadata from uploaded files).

#### Current State

`CImportUnmerged` and similar widgets contain extraction logic:

```typescript
// When file uploaded, extract metadata and update multiple parameters
useEffect(() => {
  if (!digest) return;

  await forceUpdateSPACEGROUP(digest.spaceGroup);
  await forceUpdateWAVELENGTH(digest.wavelength);
  await forceUpdateUNITCELL(digest.cell);
  // ... more updates, each potentially causing issues
}, [digest]);
```

#### Target State

Declare extraction relationships in `.def.xml`:

```xml
<container name="import_merged">
  <inputData>
    <CDataFile name="HKLIN" mimeType="application/mtz">
      <extractsTo target="SPACEGROUP" source="spaceGroup" />
      <extractsTo target="WAVELENGTH" source="wavelength" />
      <extractsTo target="UNITCELL" source="cell" />
      <extractsTo target="CRYSTALNAME" source="crystalNames[0]" />
      <extractsTo target="DATASETNAME" source="datasets[0]" />
    </CDataFile>

    <CString name="SPACEGROUP" />
    <CFloat name="WAVELENGTH" />
    <CCell name="UNITCELL" />
    <!-- ... -->
  </inputData>
</container>
```

#### Backend Implementation

When `set_parameter` or `upload_file_param` is called on a `CDataFile` with `extractsTo` declarations:

1. Set the file parameter
2. Extract digest/metadata from the file
3. Apply `extractsTo` mappings to set derived parameters
4. Return all affected parameters in response:

```typescript
{
  success: true,
  data: {
    affected_items: [
      { path: "inputData.HKLIN", item: { ... } },
      { path: "inputData.SPACEGROUP", item: { ... } },
      { path: "inputData.WAVELENGTH", item: { ... } },
      { path: "inputData.UNITCELL", item: { ... } },
    ],
    validation: { ... }
  }
}
```

#### Frontend Changes

Patching function handles multiple items:

```typescript
const uploadFileParam = useCallback(async (arg) => {
  const result = await api.post(`jobs/${job.id}/upload_file_param`, formData);

  if (result.success) {
    // Patch all affected items
    mutateContainer(
      (current) => patchContainerMultiple(current, result.data.affected_items),
      { revalidate: false }
    );
  }

  return result;
}, [...]);
```

Complex widgets become simple:

```typescript
// CImportUnmerged - AFTER
const CImportUnmerged = ({ itemName }) => {
  const { uploadFileParam } = useJob(jobId);

  const handleFileSelect = async (file) => {
    await uploadFileParam({ objectPath: itemName, file });
    // That's it. Server handles extraction, response patches all affected fields.
  };

  return <FileDropZone onFileSelect={handleFileSelect} />;
};
```

## Event-Driven Sync Points

Rather than arbitrary polling intervals, sync (full container refetch) at meaningful events:

| Event | Rationale |
|-------|-----------|
| Window/tab gains focus | User returning from elsewhere; state may have changed |
| Before job submission | Ensure server has consistent view before run |
| Job status changes | Server-side processing may have modified parameters |
| Explicit user action | "Refresh" button if provided |

## Versioning (Optional Enhancement)

For robust sync verification without timing heuristics:

```typescript
// Server maintains version per job's parameter state
interface ParameterResponse {
  // ...
  version: number;  // Incremented on each parameter change
}

// Client tracks version
const [localVersion, setLocalVersion] = useState<number | null>(null);

// On sync, compare versions
const handleSync = async () => {
  const serverContainer = await fetchContainer();
  if (serverContainer.version !== localVersion) {
    // Versions differ - full reconciliation needed
    // Could prompt user if there are local changes
  }
};
```

This is non-arbitrary: version mismatch is a definitive signal, not a timing guess.

## Migration Path

### Phase 1 Deliverables
- [ ] Modify `set_parameter` API to return `updated_item`
- [ ] Implement `patchContainer` utility
- [ ] Update `useJob.setParameter` to patch cache
- [ ] Update `useJob.uploadFileParam` to patch cache
- [ ] Remove `ParameterChangeIntentProvider`
- [ ] Remove `wasRecentlyChanged` checks from field components
- [ ] Add event-driven sync (focus, pre-run)
- [ ] Update tests

### Phase 2 Deliverables
- [ ] Define `validity()` method signature in `CPluginScript`
- [ ] Migrate validity logic for key task interfaces
- [ ] Include validation in `set_parameter` response
- [ ] Remove frontend validity filtering
- [ ] Update tests

### Phase 3 Deliverables
- [ ] Define `extractsTo` schema in `.def.xml` format
- [ ] Implement extraction logic in parameter-setting code
- [ ] Return `affected_items` array in responses
- [ ] Update `patchContainer` to handle multiple items
- [ ] Simplify complex widgets
- [ ] Update tests

## Files Affected

### Phase 1
- `server/ccp4i2/lib/utils/parameters/set_param.py` - Return updated item
- `server/ccp4i2/api/JobViewSet.py` - Modify response format
- `client/renderer/utils.ts` - Patch cache instead of invalidate
- `client/renderer/providers/parameter-change-intent-provider.tsx` - Delete
- `client/renderer/components/task/task-elements/csimple-textfield.tsx` - Remove intent checks
- `client/renderer/components/task/task-elements/*.tsx` - Remove intent checks from other field types

### Phase 2
- `core/CCP4PluginScript.py` - Add `validity()` method
- `wrappers/*/script/*.py` - Implement validity per task
- `server/ccp4i2/api/JobViewSet.py` - Include validation in response
- `client/renderer/components/task/task-interfaces/*.tsx` - Remove validity filtering

### Phase 3
- `core/CCP4Container.py` - Parse `extractsTo` declarations
- `server/ccp4i2/lib/utils/parameters/set_param.py` - Apply extractions
- `server/ccp4i2/api/JobViewSet.py` - Return affected_items
- `client/renderer/utils.ts` - Handle multiple patches
- `client/renderer/components/task/task-elements/cimportunmerged.tsx` - Simplify

## Success Criteria

### Phase 1 Complete When
- Intent system fully removed
- No timing-based heuristics in parameter handling
- Field components sync directly from props
- No race conditions on rapid parameter changes
- Task interface code significantly simpler

### Phase 2 Complete When
- All validity logic in Python plugins
- Frontend displays validation, doesn't compute it
- No `eslint-disable` comments for effect dependencies in validity code

### Phase 3 Complete When
- File upload + metadata extraction is single round-trip
- No frontend effect chains for derived parameter updates
- Complex widgets (CImportUnmerged, etc.) reduced to simple file handling

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Breaking existing task interfaces | Phased rollout; thorough testing per phase |
| Performance of deep cloning containers | Profile; consider immutable data structures if needed |
| Edge cases in container patching | Comprehensive test coverage for nested structures, lists |
| Validity migration effort | Start with most problematic interfaces; others can migrate incrementally |

## Resolved Questions

### List Items (CList Patching)

**Decision**: Use index-based patching for `CList` items.

Both frontend and backend data structures preserve item order, so index is a reliable identifier. The `_objectPath` for list items already includes the index (e.g., `inputData.SEQUENCES[0]`, `inputData.SEQUENCES[1]`).

Patching implementation:
```typescript
// Path like "task.container.inputData.SEQUENCES[2].name"
// Extract index from path segment: SEQUENCES[2] -> index 2
const listIndexRegex = /^(\w+)\[(\d+)\]$/;
```

### Computed/Derived Values

**Decision**: Not a separate concern for this refactoring.

Analysis of computed values in the codebase shows they fall into existing categories:

| Type | Examples | Handling |
|------|----------|----------|
| **Heavy calculations** | Molecular weight, Matthews coefficient | Already use server API calls with SWR caching - no change needed |
| **Pure UI state** | Visibility functions, text formatting, element dispatch | Conventional React patterns - no change needed |
| **Validation workarounds** | Synthetic values to suppress warnings (e.g., `asuFileWarningValue`) | Eliminated by Phase 2 (server-side validity) |
| **Digest-derived values** | Column groups from file metadata | Handled by Phase 3 (backend-driven extraction) |

No additional infrastructure required for computed values.

## Open Questions

1. **Undo/redo**: Local patching makes undo more feasible. Worth considering as future enhancement?

2. **Concurrent editing**: If same job opened in multiple tabs, how to handle? Version-based conflict detection?
