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
import { Box, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { isTruthy } from "../../task-elements/shared-hooks";

interface AdvancedTabProps extends CCP4i2TaskInterfaceProps {
  scatteringFactors: string;
  hydrUse: any;
  hdInitToggle: any;
  bfacSetUse: any;
  resCustom: any;
}

export const AdvancedTab: React.FC<AdvancedTabProps> = (props) => {
  const {
    scatteringFactors,
    hydrUse,
    hdInitToggle,
    bfacSetUse,
    resCustom,
    ...taskProps
  } = props;

  return (
    <>
      {/* Experiment */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Experiment" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="SCATTERING_FACTORS"
          qualifiers={{ guiLabel: "Diffraction experiment type:" }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="SCATTERING_ELECTRON"
          qualifiers={{ guiLabel: "Form factor calculation:" }}
          visibility={() => scatteringFactors === "ELECTRON"}
        />

        {/* Neutron refinement options */}
        <CCP4i2ContainerElement
          {...taskProps}
          itemName=""
          qualifiers={{ guiLabel: "Neutron refinement options" }}
          containerHint="BlockLevel"
          visibility={() => scatteringFactors === "NEUTRON"}
        >
          <FieldRow>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="HYDR_USE"
              qualifiers={{
                guiLabel: "Use hydrogens during refinement",
              }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="HYDR_ALL"
              visibility={() => isTruthy(hydrUse)}
            />
          </FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HD_INIT_TOGGLE"
            qualifiers={{
              guiLabel: "Initialise hydrogen/deuterium fractions",
            }}
            visibility={() => isTruthy(hydrUse)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HD_INIT"
            visibility={() =>
              isTruthy(hydrUse) && isTruthy(hdInitToggle)
            }
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HD_FRACTION"
            qualifiers={{
              guiLabel: "Refine hydrogen/deuterium fractions",
            }}
            visibility={() => isTruthy(hydrUse)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HD_FRACTION_TYPE"
            qualifiers={{ guiLabel: "for" }}
            visibility={() => isTruthy(hydrUse)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="H_REFINE"
            qualifiers={{ guiLabel: "Refine hydrogen positions" }}
            visibility={() => isTruthy(hydrUse)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="H_REFINE_SELECT"
            qualifiers={{ guiLabel: "for" }}
            visibility={() => isTruthy(hydrUse)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="H_TORSION"
            qualifiers={{
              guiLabel: "Use hydrogen torsion angle restraints",
            }}
            visibility={() => isTruthy(hydrUse)}
          />
        </CCP4i2ContainerElement>
      </CCP4i2ContainerElement>

      {/* Resolution */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Resolution" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="RES_CUSTOM"
          qualifiers={{ guiLabel: "Use custom resolution limits" }}
        />
        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="RES_MIN"
            qualifiers={{ guiLabel: "Low (dmax):" }}
            visibility={() => isTruthy(resCustom)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="RES_MAX"
            qualifiers={{ guiLabel: "High (dmin):" }}
            visibility={() => isTruthy(resCustom)}
          />
        </FieldRow>
      </CCP4i2ContainerElement>

      {/* Reset B-factors */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Reset B-factors" }}
        containerHint="FolderLevel"
      >
        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            gap: 1,
            flexWrap: "wrap",
          }}
        >
          <CCP4i2TaskElement
            {...taskProps}
            itemName="BFACSETUSE"
            qualifiers={{ guiLabel: "Reset all B-factors at start" }}
            sx={{ width: "auto" }}
          />
          {isTruthy(bfacSetUse) && (
            <>
              <Typography variant="body1">to fixed value:</Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="BFACSET"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </>
          )}
        </Box>
      </CCP4i2ContainerElement>

      {/* Sequence and structure */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Sequence and structure" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement {...taskProps} itemName="ASUIN" />
      </CCP4i2ContainerElement>

      {/* Extra keywords */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Extra keywords" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="EXTRAREFMACKEYWORDS"
          qualifiers={{ guiLabel: " ", guiMode: "multiLine" }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="REFMAC_KEYWORD_FILE"
        />
      </CCP4i2ContainerElement>

      {/* Cleanup */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Cleanup" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="REFMAC_CLEANUP"
          qualifiers={{
            guiLabel: "Clean up intermediate files at end of job",
          }}
        />
      </CCP4i2ContainerElement>
    </>
  );
};
