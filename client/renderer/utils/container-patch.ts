/**
 * Container patching utilities for local SWR cache updates.
 *
 * These utilities enable patching the container cache directly with updated items
 * from API responses, eliminating the need for full container refetch and the
 * timing-based intent system.
 *
 * @see docs/local-parameter-patching.md for design documentation
 */

/**
 * Represents an item in the container structure as returned by CCP4i2JsonEncoder.
 */
export interface ContainerItem {
  _class: string;
  _value: any;
  _qualifiers: Record<string, any>;
  _CONTENTS_ORDER: string[];
  _objectPath: string;
  _baseClass: string;
  _subItem?: any;
}

/**
 * Container with lookup table for O(1) path-based access.
 * Note: The API returns { container, lookup } not { result, lookup }
 */
export interface ContainerWithLookup {
  container: ContainerItem;
  lookup: Record<string, ContainerItem>;
}

/**
 * Regex to match list index segments in paths.
 * Matches patterns like: SEQUENCES[0], REFERENCE_MODELS[5]
 */
const LIST_INDEX_REGEX = /^(\w+)\[(\d+)\]$/;

/**
 * Parses a path segment that may contain a list index.
 *
 * @param segment - Path segment like "SEQUENCES" or "SEQUENCES[0]"
 * @returns Object with name and optional index
 */
function parsePathSegment(segment: string): { name: string; index?: number } {
  const match = segment.match(LIST_INDEX_REGEX);
  if (match) {
    return { name: match[1], index: parseInt(match[2], 10) };
  }
  return { name: segment };
}

/**
 * Deep clones an object using structuredClone.
 * Falls back to JSON parse/stringify for environments without structuredClone.
 */
function deepClone<T>(obj: T): T {
  if (typeof structuredClone === "function") {
    return structuredClone(obj);
  }
  return JSON.parse(JSON.stringify(obj));
}

/**
 * Removes all lookup entries that are children of a given path prefix.
 * This is necessary when patching CLists because:
 * 1. Items may be deleted from start/middle, causing reindexing
 * 2. After deepClone, object identity checks won't work
 * 3. We need to remove ALL entries under the CList path, not just exact matches
 *
 * For example, when patching UNMERGEDFILES (a CList), this removes:
 * - UNMERGEDFILES[0], UNMERGEDFILES[1], etc.
 * - UNMERGEDFILES[0].dbFileId, UNMERGEDFILES[1].annotation, etc.
 * - Any nested children of list items
 *
 * @param pathPrefix - The path prefix to match (e.g., "task.container.inputData.UNMERGEDFILES")
 * @param lookup - The lookup table to clean (mutates in place)
 */
function removeChildLookupEntries(
  pathPrefix: string,
  lookup: Record<string, ContainerItem>
): void {
  // Build regex to match the path prefix followed by [ (for list items)
  // or nothing (for the list itself)
  // This catches: UNMERGEDFILES, UNMERGEDFILES[0], UNMERGEDFILES[0].dbFileId, etc.
  const prefixWithBracket = `${pathPrefix}[`;

  const keysToRemove: string[] = [];

  for (const key of Object.keys(lookup)) {
    // Remove if key equals the prefix exactly, or starts with prefix[
    // This handles both the list itself and all indexed children
    if (key === pathPrefix || key.startsWith(prefixWithBracket)) {
      keysToRemove.push(key);
    }
  }

  for (const key of keysToRemove) {
    delete lookup[key];
  }
}

/**
 * Recursively builds lookup table entries for an item and all its children.
 * This mirrors the buildLookup function in api.ts, adding suffix paths
 * for O(1) lookup by short name (e.g., 'UNMERGEDFILES') or full path.
 *
 * @param item - The item to add to the lookup (including all children)
 * @param lookup - The lookup table to populate (mutates in place)
 */
function addToLookup(
  item: ContainerItem,
  lookup: Record<string, ContainerItem>
): void {
  // Add this item if it has an _objectPath
  const objectPath = item._objectPath;
  if (objectPath) {
    const pathElements = objectPath.split(".");

    // Add all suffix paths to lookup (mirrors buildLookup in api.ts)
    // This enables lookup by short name like 'UNMERGEDFILES' as well as full path
    for (let i = 0; i < pathElements.length; i++) {
      const subPath = pathElements.slice(-i).join(".");
      if (subPath) {
        lookup[subPath] = item;
      }
    }
    // Also add the full path explicitly (slice(-0) returns full array, not empty)
    lookup[objectPath] = item;
  }

  // Recurse into children based on type
  if (item._baseClass === "CList" && Array.isArray(item._value)) {
    // CList: _value is an array of child items
    item._value.forEach((child: any) => {
      if (child && typeof child === "object" && child._objectPath) {
        addToLookup(child, lookup);
      }
    });
  } else if (item._value && item._value.constructor === Object) {
    // CContainer or CDataFile: _value is an object with named children
    Object.values(item._value).forEach((child: any) => {
      if (child && typeof child === "object" && child._objectPath) {
        addToLookup(child, lookup);
      }
    });
  }
}

/**
 * Sets an item at a nested path within the container structure.
 *
 * Handles both regular object properties and list indices.
 * Path format: "taskname.container.inputData.SEQUENCES[0].name"
 *
 * @param container - The container object to modify (mutates in place)
 * @param path - Full object path to the item
 * @param item - The item to set at the path
 */
function setNestedItem(
  container: ContainerItem,
  path: string,
  item: ContainerItem
): void {
  const segments = path.split(".");

  // Navigate to parent, handling special segments
  let current: any = container;

  // Skip first segment (task name) and last segment (item name)
  // Process segments 1 through n-2
  for (let i = 1; i < segments.length - 1; i++) {
    const segment = segments[i];

    // Skip .container. segment
    if (segment === "container") continue;

    const { name, index } = parsePathSegment(segment);

    // Navigate into _value for container-like objects
    if (current._value && typeof current._value === "object") {
      if (index !== undefined) {
        // List access: _value is an array
        if (Array.isArray(current._value[name]?._value)) {
          current = current._value[name]._value[index];
        } else if (Array.isArray(current._value)) {
          // Direct array
          current = current._value[index];
        } else {
          console.warn(
            `[patchContainer] Expected array at ${name}, got:`,
            current._value[name]
          );
          return;
        }
      } else {
        // Regular property access
        current = current._value[name];
      }
    } else {
      console.warn(
        `[patchContainer] Cannot navigate into non-object at segment ${segment}:`,
        current
      );
      return;
    }

    if (!current) {
      console.warn(
        `[patchContainer] Path segment not found: ${segment} in ${path}`
      );
      return;
    }
  }

  // Set the item at the final segment
  const finalSegment = segments[segments.length - 1];
  const { name, index } = parsePathSegment(finalSegment);

  if (current._value && typeof current._value === "object") {
    if (index !== undefined) {
      // Setting item in a list
      if (Array.isArray(current._value[name]?._value)) {
        current._value[name]._value[index] = item;
      } else if (Array.isArray(current._value)) {
        current._value[index] = item;
      } else {
        console.warn(
          `[patchContainer] Expected array for final segment ${finalSegment}`
        );
      }
    } else {
      // Setting regular property
      current._value[name] = item;
    }
  }
}

/**
 * Patches a container with an updated item.
 *
 * Creates a deep clone of the container and updates both:
 * 1. The lookup table (O(1) access by path)
 * 2. The nested structure (for traversal)
 *
 * @param container - The current container state
 * @param updatedItem - The updated item from the API response
 * @returns New container with the item patched in
 */
export function patchContainer(
  container: ContainerWithLookup,
  updatedItem: ContainerItem
): ContainerWithLookup {
  if (!container || !updatedItem) {
    console.warn("[patchContainer] Invalid arguments:", { container, updatedItem });
    return container;
  }

  const path = updatedItem._objectPath;
  if (!path) {
    console.warn("[patchContainer] Updated item has no _objectPath:", updatedItem);
    return container;
  }

  // Deep clone to avoid mutating cached data
  const newContainer = deepClone(container);

  // For CLists (and items with children), remove all existing lookup entries
  // under this path before adding the new ones. This is critical because:
  // 1. Delete from middle causes reindexing (ITEM[2] becomes ITEM[1])
  // 2. Old entries like ITEM[2].dbFileId would remain stale
  // 3. We can't use object identity after deepClone
  //
  // For simple scalar items, this is a no-op (removes just the item itself)
  removeChildLookupEntries(path, newContainer.lookup);

  // Update the lookup table - recursively add all children too
  // This is crucial for CList items where new child items need lookup entries
  addToLookup(updatedItem, newContainer.lookup);

  // Update the nested structure
  setNestedItem(newContainer.container, path, updatedItem);

  return newContainer;
}

/**
 * Patches a container with multiple updated items.
 *
 * Used when a single API call affects multiple parameters (e.g., file upload
 * with extraction that populates derived fields).
 *
 * @param container - The current container state
 * @param updatedItems - Array of updated items from the API response
 * @returns New container with all items patched in
 */
export function patchContainerMultiple(
  container: ContainerWithLookup,
  updatedItems: ContainerItem[]
): ContainerWithLookup {
  if (!container || !updatedItems?.length) {
    return container;
  }

  // Deep clone once, then mutate
  let newContainer = deepClone(container);

  for (const item of updatedItems) {
    if (!item._objectPath) {
      console.warn("[patchContainerMultiple] Item has no _objectPath:", item);
      continue;
    }

    // Remove old lookup entries before adding new ones (handles CList reindexing)
    removeChildLookupEntries(item._objectPath, newContainer.lookup);

    // Update lookup table - recursively add all children too
    addToLookup(item, newContainer.lookup);

    // Update nested structure
    setNestedItem(newContainer.container, item._objectPath, item);
  }

  return newContainer;
}
