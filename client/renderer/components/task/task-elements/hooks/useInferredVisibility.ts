import { useMemo } from "react";

/**
 * Infer visibility from a prop that can be a boolean, function, or undefined.
 *
 * This hook consolidates the common visibility inference pattern used across
 * task element components. It handles three cases:
 * - undefined/null: defaults to visible (true)
 * - boolean: uses the boolean value directly
 * - function: calls the function to get the boolean value
 *
 * @param visibility - The visibility prop (boolean, function, or undefined)
 * @returns boolean - Whether the component should be visible
 */
export const useInferredVisibility = (
  visibility: boolean | (() => boolean) | undefined
): boolean => {
  return useMemo(() => {
    if (visibility === undefined || visibility === null) {
      return true;
    }
    if (typeof visibility === "function") {
      try {
        return visibility();
      } catch (error) {
        console.error("Error evaluating visibility function:", error);
        return true; // Default to visible on error
      }
    }
    return Boolean(visibility);
  }, [visibility]);
};

export default useInferredVisibility;
