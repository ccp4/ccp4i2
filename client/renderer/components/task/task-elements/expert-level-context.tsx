import { createContext, useContext } from "react";

/**
 * Context for PHIL expert level filtering.
 *
 * When set, CContainer auto-rendering will hide children whose
 * expertLevel qualifier exceeds this value.  undefined means no
 * filtering (all children visible).
 */
export const ExpertLevelContext = createContext<number | undefined>(undefined);

export const useExpertLevel = (): number | undefined => {
  return useContext(ExpertLevelContext);
};
