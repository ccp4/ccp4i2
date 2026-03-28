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
import React from "react";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { isTruthy } from "../../task-elements/shared-hooks";

interface OutputTabProps extends CCP4i2TaskInterfaceProps {
  mapSharp: any;
  mapSharpCustom: any;
}

export const OutputTab: React.FC<OutputTabProps> = (props) => {
  const { mapSharp, mapSharpCustom, ...taskProps } = props;

  return (
    <>
      {/* Output options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Output options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="OUTPUT_HYDROGENS"
          qualifiers={{
            guiLabel: "Output calculated riding hydrogens to file",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Map calculation */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Map calculation" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="MAP_SHARP"
          qualifiers={{
            guiLabel: "Perform map sharpening when calculating maps",
          }}
        />
        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="MAP_SHARP_CUSTOM"
            qualifiers={{
              guiLabel: "Use custom sharpening parameter (B-factor)",
            }}
            visibility={() => isTruthy(mapSharp)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="BSHARP"
            qualifiers={{ guiLabel: "B-factor" }}
            visibility={() =>
              isTruthy(mapSharp) && isTruthy(mapSharpCustom)
            }
          />
        </FieldRow>
      </CCP4i2ContainerElement>

      {/* Validation and analysis */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Validation and analysis" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="VALIDATE_IRIS"
          qualifiers={{ guiLabel: "Run IRIS for validation" }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="VALIDATE_BAVERAGE"
          qualifiers={{ guiLabel: "Analyse B-factor distributions" }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="VALIDATE_RAMACHANDRAN"
          qualifiers={{ guiLabel: "Calculate Ramachandran plots" }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="VALIDATE_MOLPROBITY"
          qualifiers={{
            guiLabel: "Run MolProbity to analyse geometry",
          }}
        />
      </CCP4i2ContainerElement>
    </>
  );
};
