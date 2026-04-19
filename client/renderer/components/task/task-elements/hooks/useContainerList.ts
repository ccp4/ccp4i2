import { useCallback, useMemo } from "react";

import { Job, Project } from "../../../../types/models";
import { SetParameterResponse, valueOfItem } from "../../../../utils";
import { useContainerField } from "./useContainerField";

export interface UseContainerListOptions {
  job: Job;
  itemName: string;
  /** Current project — used to stamp CDataFile defaults when adding rows. */
  project?: Project;
  visibility?: boolean | (() => boolean);
  disabled?: boolean | (() => boolean);
  suppressMutations?: boolean;
  onChange?: (updatedItem: any) => void | Promise<void>;
}

export interface UseContainerListResult {
  item: any;
  items: any[];
  objectPath: string | null;
  isVisible: boolean;
  isEditable: boolean;
  isSubmitting: boolean;
  validationColor: string | undefined;
  hasValidationError: boolean;

  /**
   * Append a row using the item's `_subItem` template.  Defaults are
   * derived by class (CInt→0, CFloat→0.0, CString→"", CAltSpaceGroup→"P1"),
   * with child min-qualifiers respected.  CDataFile rows are stamped with
   * the project UUID and a placeholder baseName.
   */
  addItem: () => Promise<SetParameterResponse | undefined>;

  /** Delete the row at `index`. */
  deleteAt: (index: number) => Promise<SetParameterResponse | undefined>;

  /**
   * Delete the first row matching `predicate`.  Useful when the caller
   * has a reference to the row object but not its index.
   */
  deleteMatching: (
    predicate: (row: any, index: number) => boolean
  ) => Promise<SetParameterResponse | undefined>;

  /** Replace the entire array — escape hatch for bespoke mutations. */
  replaceArray: (
    newArray: any[]
  ) => Promise<SetParameterResponse | undefined>;
}

const DEFAULT_VALUES: Record<string, any> = {
  CAltSpaceGroup: "P1",
  CInt: 0,
  CFloat: 0.0,
  CString: "",
};

/** Floor a numeric value to the item's min qualifier if defined. */
const respectMin = (val: any, item: any): any => {
  const min = item?._qualifiers?.min;
  if (min !== undefined && typeof val === "number" && val < min) return min;
  return val;
};

/**
 * Build the default value for a new list row from the list's `_subItem`
 * template, mirroring server-side default construction.  Compound types
 * build child-by-child respecting min qualifiers; CDataFile rows get
 * stamped with the current project UUID.
 */
const createNewItemValue = (
  subItemTemplate: any,
  project: Project | undefined
): any => {
  let newItemValue = valueOfItem(subItemTemplate);

  if (subItemTemplate._baseClass === "CDataFile" && newItemValue && project) {
    return {
      ...newItemValue,
      project: project.uuid.replace(/-/g, ""),
      baseName: "UNDEFINED",
    };
  }

  if (
    subItemTemplate._value &&
    typeof subItemTemplate._value === "object" &&
    !Array.isArray(subItemTemplate._value)
  ) {
    const result: Record<string, any> = {};
    for (const [key, child] of Object.entries(
      subItemTemplate._value as Record<string, any>
    )) {
      if (!child) {
        result[key] = null;
        continue;
      }
      const childDefault =
        DEFAULT_VALUES[(child as any)._class] ??
        DEFAULT_VALUES[(child as any)._baseClass];
      const val =
        childDefault !== undefined ? childDefault : valueOfItem(child);
      result[key] = respectMin(val, child);
    }
    return result;
  }

  const defaultValue =
    DEFAULT_VALUES[subItemTemplate._class] ??
    DEFAULT_VALUES[subItemTemplate._baseClass];
  const result = defaultValue !== undefined ? defaultValue : newItemValue;
  return respectMin(result, subItemTemplate);
};

/** Rewrite the `[?]` placeholder in a _subItem's paths to a concrete index. */
const stampIndex = (element: any, index: number): any => {
  if (!element) return element;
  const stamped = { ...element };
  stamped._objectPath = element._objectPath?.replace("[?]", `[${index}]`);
  if (typeof stamped._value === "object" && stamped._value) {
    stamped._value = Object.keys(stamped._value).reduce((acc: any, key) => {
      const child = stamped._value[key];
      acc[key] = child?._objectPath
        ? {
            ...child,
            _objectPath: child._objectPath.replace("[?]", `[${index}]`),
          }
        : child;
      return acc;
    }, {});
  }
  return stamped;
};

/**
 * Shared list-management plumbing for task elements that render a CList.
 * Owns the array-rewrite pattern (`setParameter({object_path: list,
 * value: newArray})`) and the gnarly default-value construction for new
 * rows, so individual list widgets only need to render their UI and call
 * `addItem` / `deleteAt`.
 *
 * Cell-level edits inside a row should use `useContainerField` with the
 * row's nested path (e.g. `parent[i].field`) — the list hook does not
 * own those writes.
 */
export function useContainerList(
  options: UseContainerListOptions
): UseContainerListResult {
  const { job, itemName, project, visibility, disabled, suppressMutations, onChange } =
    options;

  const {
    item,
    objectPath,
    isVisible,
    isSubmitting,
    validationColor,
    hasValidationError,
    commit,
  } = useContainerField<any[]>({
    job,
    itemName,
    visibility,
    disabled,
    suppressMutations,
    onChange,
  });

  const isEditable = job.status === 1;

  const items: any[] = useMemo(
    () => (Array.isArray(item?._value) ? item._value : []),
    [item]
  );

  const currentArray = useCallback((): any[] => {
    const v = valueOfItem(item);
    return Array.isArray(v) ? [...v] : [];
  }, [item]);

  const replaceArray = useCallback(
    (newArray: any[]) => commit(newArray),
    [commit]
  );

  const addItem = useCallback(async () => {
    if (!item?._subItem) {
      console.error(
        `useContainerList: item at ${itemName} has no _subItem template`
      );
      return undefined;
    }
    const template = JSON.parse(JSON.stringify(item._subItem));
    const nextIndex = items.length;
    const stamped = stampIndex(template, nextIndex);
    const newRow = createNewItemValue(stamped, project);

    const next = currentArray();
    next.push(newRow);
    return commit(next);
  }, [item, itemName, items.length, project, currentArray, commit]);

  const deleteAt = useCallback(
    async (index: number) => {
      const next = currentArray();
      if (index < 0 || index >= next.length) return undefined;
      next.splice(index, 1);
      return commit(next);
    },
    [currentArray, commit]
  );

  const deleteMatching = useCallback(
    async (predicate: (row: any, index: number) => boolean) => {
      const next = currentArray();
      const idx = next.findIndex((row, i) => predicate(row, i));
      if (idx === -1) return undefined;
      next.splice(idx, 1);
      return commit(next);
    },
    [currentArray, commit]
  );

  return {
    item,
    items,
    objectPath,
    isVisible,
    isEditable,
    isSubmitting,
    validationColor,
    hasValidationError,
    addItem,
    deleteAt,
    deleteMatching,
    replaceArray,
  };
}
