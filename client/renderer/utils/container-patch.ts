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
 */
export interface ContainerWithLookup {
  result: ContainerItem;
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

  // Update the lookup table
  newContainer.lookup[path] = updatedItem;

  // Update the nested structure
  setNestedItem(newContainer.result, path, updatedItem);

  console.log(`[patchContainer] Patched item at path: ${path}`);
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

    // Update lookup table
    newContainer.lookup[item._objectPath] = item;

    // Update nested structure
    setNestedItem(newContainer.result, item._objectPath, item);

    console.log(`[patchContainerMultiple] Patched item at path: ${item._objectPath}`);
  }

  return newContainer;
}
