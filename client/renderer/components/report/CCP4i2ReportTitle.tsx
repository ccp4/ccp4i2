/*
 * Copyright (C) 2025 Newcastle University
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
import React from "react";
import { CardHeader, Skeleton } from "@mui/material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";

export const CCP4i2ReportTitle: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  return false ? (
    <CardHeader
      title={$(props.item).attr("title1")}
      subtitle={$(props.item).attr("title2")}
    />
  ) : (
    <Skeleton />
  );
};
