import React, { useMemo } from "react";
import { Box, Chip, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { useJob } from "../../../../utils";

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

/**
 * Hook to compute the maximum resolution across all files in the UNMERGEDFILES list.
 * Reads highRes from each list item's file digest (CUnmergedDataContent.highRes).
 */
const useMaxResolution = (jobId: number) => {
  const { useTaskItem, useFileDigest } = useJob(jobId);
  const { item: listItem } = useTaskItem("UNMERGEDFILES");
  const listItems: any[] = listItem?._value || [];

  // Build file object paths for all list items that have a file uploaded
  const filePaths = useMemo(() => {
    return listItems
      .map((entry: any) => {
        const fileValue = entry?.file;
        if (!fileValue?.dbFileId) return "";
        // Construct the objectPath for the file within this list entry
        return entry?._objectPath ? `${entry._objectPath}.file` : "";
      })
      .filter(Boolean);
  }, [listItems]);

  // Fetch digest for first file (most common case is single file)
  const { data: digest0 } = useFileDigest(filePaths[0] || "");
  const { data: digest1 } = useFileDigest(filePaths[1] || "");
  const { data: digest2 } = useFileDigest(filePaths[2] || "");

  return useMemo(() => {
    const digests = [digest0, digest1, digest2].filter(Boolean);
    if (digests.length === 0) return null;

    let best: number | null = null;
    for (const d of digests) {
      const hr = d?.highRes ?? d?.resolutionRange?.high;
      if (typeof hr === "number" && hr > 0) {
        if (best === null || hr < best) best = hr;
      }
    }
    return best;
  }, [digest0, digest1, digest2]);
};

export const InputDataTab: React.FC<InputDataTabProps> = (props) => {
  const { mode, referenceDataset, vis, ...taskProps } = props;
  const maxResolution = useMaxResolution(props.job.id);

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
      <FieldRow>
        <CCP4i2TaskElement
          itemName="RESOLUTION_RANGE"
          {...taskProps}
          qualifiers={{ guiLabel: "Resolution range (Å)" }}
        />
        {maxResolution != null && (
          <Box sx={{ display: "flex", alignItems: "center", gap: 1, ml: 2 }}>
            <Typography variant="body2" color="text.secondary">
              Maximum resolution in files
            </Typography>
            <Chip
              label={`${maxResolution.toFixed(2)}Å`}
              size="small"
              variant="outlined"
              color="primary"
            />
          </Box>
        )}
      </FieldRow>
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
