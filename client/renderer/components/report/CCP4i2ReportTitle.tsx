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
