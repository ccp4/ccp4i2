/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
