import { useCallback, useEffect, useMemo, useState } from "react";

import { Job } from "../../../../types/models";
import {
  useJob,
  SetParameterArg,
  SetParameterResponse,
  valueOfItem,
} from "../../../../utils";
import { useTaskInterface } from "../../../../providers/task-provider";
import { usePopcorn } from "../../../../providers/popcorn-provider";

export interface UseContainerFieldOptions {
  job: Job;
  itemName: string;
  visibility?: boolean | (() => boolean);
  disabled?: boolean | (() => boolean);
  suppressMutations?: boolean;
  onChange?: (updatedItem: any) => void | Promise<void>;
}

export interface CommitOptions {
  /**
   * Append to the item's objectPath before writing (e.g. ".selection" or
   * ".selected"). Use when the write target is a child of the hook's item
   * rather than the item itself.
   */
  subPath?: string;
}

export interface UseContainerFieldResult<V = any> {
  item: any;
  objectPath: string | null;
  /**
   * Raw server-side `_value` as it sits on the item. For primitives this is
   * the actual value; for complex items (e.g. CDataFile) this is a nested
   * object of `{ _value, _class, ... }` children. Prefer `unwrappedValue`
   * when you want plain JS values.
   */
  serverValue: V;
  /**
   * `_value` recursively unwrapped — nested `_value` children are flattened
   * so a CDataFile becomes `{ dbFileId, baseName, ... }` rather than
   * `{ dbFileId: {_value, _class}, baseName: {_value, _class}, ... }`.
   */
  unwrappedValue: any;
  isVisible: boolean;
  isDisabled: boolean;
  isSubmitting: boolean;
  validationColor: string | undefined;
  hasValidationError: boolean;
  commit: (
    newValue: any,
    options?: CommitOptions
  ) => Promise<SetParameterResponse | undefined>;
}

/**
 * Shared field plumbing for task elements: item lookup, visibility,
 * disabled state, validation color, and the parameter-write flow
 * (submission state, error toast, revert on failure, optional onChange).
 *
 * Element components own their own local editing state (debounced text,
 * autocomplete input, etc.) and call `commit` when a value should be sent
 * to the server.
 */
export function useContainerField<V = any>(
  options: UseContainerFieldOptions
): UseContainerFieldResult<V> {
  const { job, itemName, visibility, disabled, suppressMutations, onChange } =
    options;

  const {
    useTaskItem,
    getValidationColor,
    setParameter,
    setParameterNoMutate,
  } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const { inFlight, setInFlight } = useTaskInterface();
  const { setMessage } = usePopcorn();

  const [isSubmitting, setIsSubmitting] = useState(false);

  const objectPath: string | null = useMemo(
    () => item?._objectPath || null,
    [item]
  );

  const serverValue = item?._value as V;
  const unwrappedValue = useMemo(() => valueOfItem(item), [item]);

  const isVisible = useMemo(() => {
    if (typeof visibility === "function") {
      try {
        return visibility();
      } catch {
        return true;
      }
    }
    return visibility !== false;
  }, [visibility]);

  const isDisabled = useMemo(() => {
    const base = inFlight || isSubmitting || job.status !== 1;
    if (typeof disabled === "function") return disabled() || base;
    return Boolean(disabled) || base;
  }, [disabled, inFlight, isSubmitting, job.status]);

  const validationColor = useMemo(
    () => getValidationColor(item),
    [getValidationColor, item]
  );

  const hasValidationError = validationColor === "error.light";

  const commit = useCallback(
    async (
      newValue: any,
      options?: CommitOptions
    ): Promise<SetParameterResponse | undefined> => {
      if (!objectPath) {
        console.error(`useContainerField: no objectPath for ${itemName}`);
        return undefined;
      }

      const writePath = options?.subPath
        ? `${objectPath}${options.subPath}`
        : objectPath;
      const arg: SetParameterArg = { object_path: writePath, value: newValue };
      const updateFn = suppressMutations ? setParameterNoMutate : setParameter;

      setInFlight(true);
      setIsSubmitting(true);
      try {
        const result = await updateFn(arg);

        if (result && !result.success) {
          setMessage(`Unacceptable value provided: "${newValue}"`);
          return result;
        }
        if (result?.success && result.data?.updated_item && onChange) {
          await onChange(result.data.updated_item);
        }
        return result;
      } catch (error) {
        const errorMessage =
          error instanceof Error ? error.message : String(error);
        setMessage(`Error updating parameter: ${errorMessage}`);
        console.error(`Parameter update failed for ${itemName}:`, error);
        return undefined;
      } finally {
        setInFlight(false);
        setIsSubmitting(false);
      }
    },
    [
      objectPath,
      itemName,
      suppressMutations,
      setParameter,
      setParameterNoMutate,
      setInFlight,
      setMessage,
      onChange,
    ]
  );

  return {
    item,
    objectPath,
    serverValue,
    unwrappedValue,
    isVisible,
    isDisabled,
    isSubmitting,
    validationColor,
    hasValidationError,
    commit,
  };
}

/**
 * Small helper for elements that keep local editing state: whenever the
 * server value changes, adopt it.  Returns [local, setLocal].
 */
export function useSyncedLocalValue<V>(serverValue: V): [V, (v: V) => void] {
  const [local, setLocal] = useState<V>(serverValue);
  useEffect(() => {
    setLocal(serverValue);
  }, [serverValue]);
  return [local, setLocal];
}
