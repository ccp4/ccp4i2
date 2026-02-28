import React from "react";
import { Box, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";

interface InputDataTabProps extends CCP4i2TaskInterfaceProps {
  mode: string;
  chooseMode: string;
  referenceDataset: string;
  vis: {
    isChooseMode: () => boolean;
    isChooseSolution: () => boolean;
    isChooseSpacegroup: () => boolean;
    isReindexSpace: () => boolean;
    isChooseLauegroup: () => boolean;
  };
}

export const InputDataTab: React.FC<InputDataTabProps> = (props) => {
  const { mode, referenceDataset, vis, ...taskProps } = props;

  return (
    <>
      {/* Select unmerged data files */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Select unmerged data files" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="UNMERGEDFILES" {...taskProps} />
      </CCP4i2ContainerElement>

      {/* Resolution */}
      <CCP4i2TaskElement
        itemName="RESOLUTION_RANGE"
        {...taskProps}
        qualifiers={{ guiLabel: "Resolution range (Å)" }}
      />
      <CCP4i2TaskElement
        itemName="POINTLESS_USE_RESOLUTION_RANGE"
        {...taskProps}
        qualifiers={{
          guiLabel:
            "use explicit resolution range in symmetry determination as well as in scaling",
        }}
      />
      <CCP4i2TaskElement
        itemName="AUTOCUTOFF"
        {...taskProps}
        qualifiers={{
          guiLabel:
            "automatically cut resolution range based on a first Aimless run",
        }}
      />

      {/* Symmetry determination */}
      <CCP4i2TaskElement
        itemName="MODE"
        {...taskProps}
        qualifiers={{ guiLabel: "Options for symmetry determination" }}
      />

      {/* Choose mode options (visible when MODE === "CHOOSE") */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Options for choice of space group or Laue group",
        }}
        containerHint="BlockLevel"
        visibility={vis.isChooseMode}
      >
        <FieldRow>
          <CCP4i2TaskElement
            itemName="CHOOSE_MODE"
            {...taskProps}
            qualifiers={{ guiLabel: " " }}
          />
          <CCP4i2TaskElement
            itemName="CHOOSE_SPACEGROUP"
            {...taskProps}
            visibility={vis.isChooseSpacegroup}
          />
          <CCP4i2TaskElement
            itemName="CHOOSE_SOLUTION_NO"
            {...taskProps}
            visibility={vis.isChooseSolution}
          />
          <CCP4i2TaskElement
            itemName="CHOOSE_LAUEGROUP"
            {...taskProps}
            visibility={vis.isChooseLauegroup}
          />
        </FieldRow>
        <CCP4i2TaskElement
          itemName="REINDEX_OPERATOR"
          {...taskProps}
          qualifiers={{ guiLabel: "Reindex operator" }}
          visibility={vis.isReindexSpace}
        />
      </CCP4i2ContainerElement>

      {/* Optional reference data */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel:
            "Optional reference data to resolve indexing ambiguity and space group (also later option to scale against this reference set)",
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
          <Typography variant="body1">Reference data are</Typography>
          <Box sx={{ width: "12rem" }}>
            <CCP4i2TaskElement
              itemName="REFERENCE_DATASET"
              {...taskProps}
              qualifiers={{ guiLabel: " " }}
            />
          </Box>
          <Typography variant="body2" sx={{ fontStyle: "italic" }}>
            {mode === "MATCH"
              ? "and MUST be defined in next line"
              : "and is optionally defined in next line"}
          </Typography>
        </Box>

        {referenceDataset === "XYZ" ? (
          <CCP4i2TaskElement
            itemName="XYZIN_REF"
            {...taskProps}
            qualifiers={{ guiLabel: "Atomic model" }}
          />
        ) : (
          <CCP4i2TaskElement
            itemName="HKLIN_REF"
            {...taskProps}
            qualifiers={{ guiLabel: "Reflections" }}
          />
        )}
      </CCP4i2ContainerElement>

      {/* Optional existing FreeR set */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel:
            "Optional existing FreeR set, define to copy or extend if necessary",
        }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="FREERFLAG"
          {...taskProps}
          qualifiers={{ guiLabel: "Free R set" }}
        />
        <CCP4i2TaskElement
          itemName="CUTRESOLUTION"
          {...taskProps}
          qualifiers={{
            guiLabel:
              "Cut resolution of FreeR set if necessary to match the data",
          }}
        />
      </CCP4i2ContainerElement>
    </>
  );
};
