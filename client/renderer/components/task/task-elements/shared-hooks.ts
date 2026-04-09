import { useCallback, useEffect, useState } from "react";

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
