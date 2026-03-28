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
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { CCP4i2ContainerElement } from "./ccontainer";
import { Typography } from "@mui/material";

export const CRangeElement: React.FC<CCP4i2TaskElementProps> = (props) => (
  <CCP4i2ContainerElement
    {...props}
    itemName=""
    qualifiers={props.qualifiers}
    containerHint="RowLevel"
  >
    <Typography variant="body2">Resolution</Typography>
    <CCP4i2TaskElement job={props.job} itemName={`${props.itemName}.start`} />
    <CCP4i2TaskElement job={props.job} itemName={`${props.itemName}.end`} />
  </CCP4i2ContainerElement>
);
