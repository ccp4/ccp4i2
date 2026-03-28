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
import { InlineField } from "../../task-elements/inline-field";
import { isTruthy } from "../../task-elements/shared-hooks";

interface ParameterisationTabProps extends CCP4i2TaskInterfaceProps {
  refinementMode: string;
  solventMaskType: string;
  solventAdvanced: any;
  tlsMode: string;
  bfacSetUse: any;
  weightOpt: string;
  occupancyRefinement: any;
  occupancyGroups: any;
  occupancyComplete: any;
  occupancyIncomplete: any;
}

const NotAvailableMessage = () => (
  <Typography variant="body2" sx={{ fontStyle: "italic" }}>
    Not available in Rigid Body mode.
  </Typography>
);

export const ParameterisationTab: React.FC<ParameterisationTabProps> = (
  props
) => {
  const {
    refinementMode,
    solventMaskType,
    solventAdvanced,
    tlsMode,
    bfacSetUse,
    weightOpt,
    occupancyRefinement,
    occupancyGroups,
    occupancyComplete,
    occupancyIncomplete,
    ...taskProps
  } = props;

  const isRigidBody = refinementMode === "RIGID";

  return (
    <>
      {/* B-factors */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "B-factors" }}
        containerHint="FolderLevel"
      >
        {isRigidBody ? (
          <NotAvailableMessage />
        ) : (
          <InlineField
            label="Refine"
            width="16rem"
            after={<Typography variant="body1">B-factors</Typography>}
          >
            <CCP4i2TaskElement
              {...taskProps}
              itemName="B_REFINEMENT_MODE"
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>
        )}
      </CCP4i2ContainerElement>

      {/* Scaling */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Scaling" }}
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
          <Typography variant="body1">Use</Typography>
          <Box sx={{ width: "10rem" }}>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="SCALE_TYPE"
              qualifiers={{ guiLabel: " " }}
            />
          </Box>
          <Typography variant="body1">scaling, with</Typography>
          <Box sx={{ width: "10rem" }}>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="SOLVENT_MASK_TYPE"
              qualifiers={{ guiLabel: " " }}
            />
          </Box>
          <Typography variant="body1">solvent mask</Typography>
        </Box>

        <CCP4i2TaskElement
          {...taskProps}
          itemName="SOLVENT_ADVANCED"
          qualifiers={{ guiLabel: "Use custom solvent mask parameters" }}
          visibility={() => solventMaskType === "EXPLICIT"}
        />

        <CCP4i2ContainerElement
          {...taskProps}
          itemName=""
          qualifiers={{ guiLabel: "Custom solvent mask parameters" }}
          containerHint="BlockLevel"
          visibility={() =>
            solventMaskType === "EXPLICIT" && isTruthy(solventAdvanced)
          }
        >
          <CCP4i2TaskElement
            {...taskProps}
            itemName="SOLVENT_VDW_RADIUS"
            qualifiers={{
              guiLabel: "Increase VDW Radius of non-ion atoms by",
            }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="SOLVENT_IONIC_RADIUS"
            qualifiers={{
              guiLabel: "Increase VDW Radius of potential ion atoms by",
            }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="SOLVENT_SHRINK"
            qualifiers={{
              guiLabel: "Shrink the mask area by a factor of",
            }}
          />
        </CCP4i2ContainerElement>
      </CCP4i2ContainerElement>

      {/* Translation-Libration-Screw (TLS) */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Translation-Libration-Screw (TLS)" }}
        containerHint="FolderLevel"
      >
        {isRigidBody ? (
          <NotAvailableMessage />
        ) : (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="TLSMODE"
              qualifiers={{ guiLabel: "TLS parameters:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="NTLSCYCLES"
              qualifiers={{ guiLabel: "Number of TLS cycles:" }}
              visibility={() => tlsMode !== "NONE"}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="TLSIN"
              qualifiers={{ guiLabel: "TLS group definitions:" }}
              visibility={() => tlsMode === "FILE"}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="TLSTEXT"
              qualifiers={{ guiLabel: "TLS group definitions:" }}
              visibility={() => tlsMode !== "FILE" && tlsMode !== "NONE"}
            />

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
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="BFACSET"
                  qualifiers={{ guiLabel: "to fixed value:" }}
                  visibility={() => isTruthy(bfacSetUse)}
                />
              </Box>
            </Box>

            <CCP4i2TaskElement
              {...taskProps}
              itemName="TLSOUT_ADDU"
              qualifiers={{
                guiLabel:
                  "Add TLS contribution to output B-factors (only for analysis and deposition)",
              }}
              visibility={() => tlsMode !== "NONE"}
            />
          </>
        )}
      </CCP4i2ContainerElement>

      {/* Rigid body groups — only in Rigid Body mode */}
      {isRigidBody && (
        <CCP4i2ContainerElement
          {...taskProps}
          itemName=""
          qualifiers={{ guiLabel: "Rigid body groups" }}
          containerHint="FolderLevel"
        >
          <CCP4i2TaskElement
            {...taskProps}
            itemName="RIGID_BODY_SELECTION"
          />
        </CCP4i2ContainerElement>
      )}

      {/* Weights */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Weights" }}
        containerHint="FolderLevel"
      >
        <InlineField
          label="Weight restraints versus experimental data using"
          width="10rem"
          after={<Typography variant="body1">weighting</Typography>}
        >
          <CCP4i2TaskElement
            {...taskProps}
            itemName="WEIGHT_OPT"
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
        <CCP4i2TaskElement
          {...taskProps}
          itemName="controlParameters.WEIGHT"
          qualifiers={{ guiLabel: "Weight:" }}
          visibility={() => weightOpt === "MANUAL"}
        />
      </CCP4i2ContainerElement>

      {/* Conformer groups and occupancy refinement */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Conformer groups and occupancy refinement" }}
        containerHint="FolderLevel"
      >
        {isRigidBody ? (
          <NotAvailableMessage />
        ) : (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_REFINEMENT"
              qualifiers={{ guiLabel: "Refine conformer occupancies" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_GROUPS"
              qualifiers={{ guiLabel: "Specify partial occupancy groups" }}
              visibility={() => isTruthy(occupancyRefinement)}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_SELECTION"
              visibility={() =>
                isTruthy(occupancyRefinement) && isTruthy(occupancyGroups)
              }
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_COMPLETE"
              qualifiers={{ guiLabel: "Complete groups (sum to 1.0)" }}
              visibility={() =>
                isTruthy(occupancyRefinement) && isTruthy(occupancyGroups)
              }
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_COMPLETE_TABLE"
              visibility={() =>
                isTruthy(occupancyRefinement) &&
                isTruthy(occupancyGroups) &&
                isTruthy(occupancyComplete)
              }
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_INCOMPLETE"
              qualifiers={{ guiLabel: "Incomplete groups (sum to < 1.0)" }}
              visibility={() =>
                isTruthy(occupancyRefinement) && isTruthy(occupancyGroups)
              }
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="OCCUPANCY_INCOMPLETE_TABLE"
              visibility={() =>
                isTruthy(occupancyRefinement) &&
                isTruthy(occupancyGroups) &&
                isTruthy(occupancyIncomplete)
              }
            />
          </>
        )}
      </CCP4i2ContainerElement>
    </>
  );
};
