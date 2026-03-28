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

interface InputDataTabProps extends CCP4i2TaskInterfaceProps {
  isXtal: boolean;
  isSpa: boolean;
  isMerged: boolean;
  intensitiesAvailable: boolean;
  hydrUse: any;
  addWaters: any;
  useAnomalous: any;
  HKLINValue: any;
}

export const InputDataTab: React.FC<InputDataTabProps> = (props) => {
  const {
    isXtal,
    isSpa,
    isMerged,
    intensitiesAvailable,
    hydrUse,
    addWaters,
    useAnomalous,
    HKLINValue,
    ...taskProps
  } = props;

  return (
    <>
      {/* Data source */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Data source" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="DATA_METHOD"
          qualifiers={{ guiLabel: "Data type:" }}
        />
        {isXtal && (
          <CCP4i2TaskElement
            {...taskProps}
            itemName="MERGED_OR_UNMERGED"
            qualifiers={{ guiLabel: "Diffraction data are" }}
          />
        )}
      </CCP4i2ContainerElement>

      {/* Main inputs */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Main inputs" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement {...taskProps} itemName="container.inputData.XYZIN" />

        {/* X-ray merged inputs */}
        {isXtal && isMerged && (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="HKLIN"
              qualifiers={{ guiLabel: "Reflections" }}
            />
            {intensitiesAvailable ? (
              <CCP4i2TaskElement
                {...taskProps}
                itemName="F_SIGF_OR_I_SIGI"
                qualifiers={{ guiLabel: "Refinement against" }}
              />
            ) : (
              <Typography
                variant="body1"
                sx={{ fontWeight: "medium", mt: 1 }}
              >
                Using <strong>amplitudes</strong>
              </Typography>
            )}
            <CCP4i2TaskElement {...taskProps} itemName="FREERFLAG" />
          </>
        )}

        {/* X-ray unmerged inputs */}
        {isXtal && !isMerged && (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="HKLIN_UNMERGED"
              qualifiers={{ guiLabel: "Unmerged reflection data" }}
            />
            <CCP4i2TaskElement {...taskProps} itemName="FREERFLAG" />
          </>
        )}

        {/* SPA map inputs */}
        {isSpa && (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="MAPIN1"
              qualifiers={{ guiLabel: "Half map 1" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="MAPIN2"
              qualifiers={{ guiLabel: "Half map 2" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="MAPMASK"
              qualifiers={{ guiLabel: "Mask" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="RES_MIN"
              qualifiers={{ guiLabel: "Resolution:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="MASK_RADIUS"
              qualifiers={{ guiLabel: "Mask radius:" }}
            />
          </>
        )}

        {isXtal && (
          <CCP4i2TaskElement
            {...taskProps}
            itemName="USE_TWIN"
            qualifiers={{ guiLabel: "Twin refinement" }}
          />
        )}
      </CCP4i2ContainerElement>

      {/* Additional geometry dictionaries */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Additional geometry dictionaries" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement {...taskProps} itemName="DICT_LIST" />
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="NCYCLES"
          qualifiers={{ guiLabel: "Number of refinement cycles:" }}
        />

        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HYDR_USE"
            qualifiers={{
              guiLabel: "Use riding hydrogens during refinement",
            }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HYDR_ALL"
            visibility={() => isTruthy(hydrUse)}
          />
        </FieldRow>

        {isXtal && (
          <>
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
                itemName="ADD_WATERS"
                qualifiers={{
                  guiLabel: "Add waters and then perform",
                }}
                sx={{ width: "auto" }}
              />
              {isTruthy(addWaters) && (
                <>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="NCYCLES_AFTER_ADD_WATERS"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    further refinement cycles.
                  </Typography>
                </>
              )}
            </Box>
          </>
        )}
      </CCP4i2ContainerElement>

      {/* Anomalous signal (xtal only, when anomalous data available) */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Anomalous signal" }}
        containerHint="FolderLevel"
        visibility={() =>
          isXtal &&
          HKLINValue?.contentFlag !== undefined &&
          [1, 2].includes(HKLINValue.contentFlag)
        }
      >
        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="USEANOMALOUS"
            qualifiers={{ guiLabel: "Use anomalous" }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="USEANOMALOUSFOR"
            qualifiers={{ guiLabel: "Use for:" }}
            visibility={() => isTruthy(useAnomalous)}
          />
        </FieldRow>
      </CCP4i2ContainerElement>
    </>
  );
};
