import { useApi } from "../../../api";
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Button, Card, CardContent, CardHeader, Grid2 } from "@mui/material";
import { ErrorInfo } from "./error-info";
import { useMemo, useState } from "react";
import { ExpandLess, ExpandMore } from "@mui/icons-material";
import { CCP4i2ContainerElement } from "./ccontainer";

export const CReindexOperatorElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => (
  <CCP4i2ContainerElement
    {...props}
    qualifiers={props.qualifiers}
    containerHint="RowLevel"
    elementSx={{ width: "5rem" }}
  />
);
