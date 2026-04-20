import { useCallback, useEffect, useState } from "react";
import { useTaskInterface } from "./task-interface-context";

/** Normalize CBoolean values - server may return boolean or string */
export const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

/** Return type from useBoolToggle */
export interface BoolToggle {
  value: boolean;
  onChange: (item: any) => Promise<void>;
}

/**
 * Hook to manage a CBoolean toggle with local state for immediate UI response.
 * Standard pattern: local state + useEffect sync from server + onChange callback.
 */
export function useBoolToggle(
  useTaskItemFn: (name: string) => { value: any },
  itemName: string
): BoolToggle {
  const { value: raw } = useTaskItemFn(itemName);
  const [value, setValue] = useState(() => isTruthy(raw));
  useEffect(() => setValue(isTruthy(raw)), [raw]);
  const onChange = useCallback(async (item: any) => {
    setValue(isTruthy(item._value));
  }, []);
  return { value, onChange };
}

/**
 * Batch version of useBoolToggle. Must be called inside <TaskInterfaceProvider>.
 *
 * Pass a readonly tuple of item names — the length must be stable across renders
 * (normal for a task interface, since names are literals). Returns a keyed record
 * so callers can write `toggles.ACORN_BECUT.value` / `toggles.ACORN_BECUT.onChange`.
 *
 * Example:
 *   const toggles = useTaskToggles([
 *     "ACOPH_CUSTOM", "ACORN_BECUT", "ACORN_BRESOL",
 *   ] as const);
 */
export function useTaskToggles<const Names extends readonly string[]>(
  names: Names
): { [K in Names[number]]: BoolToggle } {
  const { useTaskItem } = useTaskInterface();
  const result = {} as { [K in Names[number]]: BoolToggle };
  for (const name of names) {
    // Hook order is stable: `names` is a const tuple whose length must not change.
    // eslint-disable-next-line react-hooks/rules-of-hooks
    (result as any)[name] = useBoolToggle(useTaskItem, name);
  }
  return result;
}
