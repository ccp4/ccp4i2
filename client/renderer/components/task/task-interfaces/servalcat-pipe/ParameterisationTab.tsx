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
import { isTruthy } from "../../task-elements/shared-hooks";

interface ParameterisationTabProps extends CCP4i2TaskInterfaceProps {
  occupancyGroups: any;
  occupancyComplete: any;
  occupancyIncomplete: any;
  occupancyRefinement: any;
}

export const ParameterisationTab: React.FC<ParameterisationTabProps> = (
  props
) => {
  const {
    occupancyGroups,
    occupancyComplete,
    occupancyIncomplete,
    occupancyRefinement,
    ...taskProps
  } = props;

  return (
    <>
      {/* Atomic displacement parameters (ADPs) */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Atomic displacement parameters (ADPs)",
        }}
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
          <Box sx={{ width: "10rem" }}>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="B_REFINEMENT_MODE"
              qualifiers={{ guiLabel: " " }}
            />
          </Box>
          <Typography variant="body1">ADPs</Typography>
        </Box>
      </CCP4i2ContainerElement>

      {/* Scaling */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Scaling" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="NO_SOLVENT"
          qualifiers={{
            guiLabel: "Do not consider bulk solvent contribution",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Conformer groups and occupancy refinement */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Conformer groups and occupancy refinement",
        }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="OCCUPANCY_GROUPS"
          qualifiers={{
            guiLabel:
              "Specify partial occupancy groups (alternative conformers)",
          }}
        />
        {isTruthy(occupancyGroups) && (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_SELECTION"
            />

            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_COMPLETE"
              qualifiers={{
                guiLabel:
                  "Specify overlapping alternative conformer groups (constrain occupancies to sum to one)",
              }}
            />
            {isTruthy(occupancyComplete) && (
              <CCP4i2TaskElement
                {...taskProps}
                itemName="OCCUPANCY_COMPLETE_TABLE"
              />
            )}

            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_INCOMPLETE"
              qualifiers={{
                guiLabel:
                  "Specify overlapping alternative conformer groups (occupancies sum to less than one)",
              }}
            />
            {isTruthy(occupancyIncomplete) && (
              <CCP4i2TaskElement
                {...taskProps}
                itemName="OCCUPANCY_INCOMPLETE_TABLE"
              />
            )}
          </>
        )}
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
            itemName="OCCUPANCY_REFINEMENT"
            qualifiers={{
              guiLabel:
                "Perform refinement of atomic occupancies every",
            }}
            sx={{ width: "auto" }}
          />
          {isTruthy(occupancyRefinement) && (
            <>
              <Box sx={{ width: "6rem" }}>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="OCCUPANCY_NCYCLE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">cycle.</Typography>
            </>
          )}
        </Box>
      </CCP4i2ContainerElement>
    </>
  );
};
