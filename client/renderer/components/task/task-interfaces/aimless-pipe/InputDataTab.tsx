import React, { useCallback, useEffect, useMemo } from "react";
import {
  Alert,
  Box,
  Chip,
  List,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Radio,
  Typography,
} from "@mui/material";
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

interface RBlockInfo {
  bname: string;
  entry_id?: string;
  cell?: number[];
  spacegroup_name?: string;
  wavelength?: number;
  highres?: number;
  unmerged?: boolean;
  info?: string;
  hklcheckformat?: string;
  suitableForMerging?: boolean;
  columnSetsText?: string[];
  hasFreeR?: boolean;
  freerValid?: boolean;
  freerWarnings?: string[];
  details?: string;
}

/** Build a multi-line description of a reflection block, mirroring
 *  the legacy Qt interface display. */
function formatBlockDescription(b: RBlockInfo): string {
  const lines: string[] = [];
  // Column content (e.g. "unmerged I")
  if (b.info) lines.push(`Column content type: ${b.info}`);
  // hkl list assessment (e.g. "Unmerged data")
  if (b.hklcheckformat) lines.push(`hkl list type: ${b.hklcheckformat}`);
  // FreeR status
  if (b.hasFreeR) {
    const status = b.freerValid ? "present" : "present but invalid";
    lines.push(`FreeR status: ${status}`);
    if (b.freerWarnings?.length) {
      for (const w of b.freerWarnings) lines.push(w);
    }
  } else {
    lines.push("FreeR status: absent");
  }
  // Cell & spacegroup
  const meta: string[] = [];
  if (b.spacegroup_name) meta.push(b.spacegroup_name);
  if (b.highres !== undefined) meta.push(`${b.highres.toFixed(2)} Å`);
  if (meta.length) lines.push(meta.join(", "));
  // Diffrn details if present
  if (b.details) lines.push(b.details);
  return lines.join("\n");
}

/**
 * mmCIF block selector for unmerged files.
 *
 * Uses the container-level MMCIF_SELECTED_BLOCK parameter, so only a
 * single mmCIF input file is supported at present.
 *
 * TODO: extend or specialise CImportUnmerged to carry a per-entry
 * selectedBlock field so that multiple mmCIF inputs can each specify
 * their own block independently.
 */
const MmcifBlockSelector: React.FC<{ jobId: number }> = ({ jobId }) => {
  const { useTaskItem, useFileDigest, mutateValidation } = useJob(jobId);
  const { item: listItem } = useTaskItem("UNMERGEDFILES");
  const { value: selectedBlock, forceUpdate: setSelectedBlock } =
    useTaskItem("MMCIF_SELECTED_BLOCK");
  const { forceUpdate: setSelectedInfo } =
    useTaskItem("MMCIF_SELECTED_INFO");
  const { forceUpdate: setSelectedColumns } =
    useTaskItem("MMCIF_SELECTED_COLUMNS");
  const { forceUpdate: setSelectedDetails } =
    useTaskItem("MMCIF_SELECTED_DETAILS");

  const listItems: any[] = listItem?._value || [];

  // The _value list entries are shallow stubs — the file sub-object
  // is not hydrated there. Instead, use useTaskItem to read the file
  // value from the first list entry via its known path.
  const firstEntryPath = listItems[0]?._objectPath || "";
  const { value: firstFileValue } = useTaskItem(
    firstEntryPath ? `${firstEntryPath}.file` : ""
  );

  // Build the digest path when the file is set
  const firstFilePath =
    firstFileValue?.dbFileId && firstEntryPath
      ? `${firstEntryPath}.file`
      : "";

  const { data: digest } = useFileDigest(firstFilePath) as {
    data: { format?: string; rblock_infos?: RBlockInfo[] } | null;
  };

  // Only show the block selector when the digest confirms mmCIF format
  const isMMCIF = digest?.format?.toLowerCase() === "mmcif";
  const blocks = digest?.rblock_infos || [];
  const suitable = useMemo(
    () => blocks.filter((b) => b.suitableForMerging),
    [blocks]
  );
  const unsuitable = useMemo(
    () => blocks.filter((b) => !b.suitableForMerging),
    [blocks]
  );

  // Auto-select when there is exactly one suitable block
  useEffect(() => {
    if (suitable.length === 1 && !selectedBlock && setSelectedBlock) {
      const b = suitable[0];
      setSelectedBlock(b.bname);
      if (setSelectedInfo) setSelectedInfo(b.info || "");
      if (setSelectedColumns)
        setSelectedColumns(b.columnSetsText?.join(", ") || "");
      if (setSelectedDetails) setSelectedDetails(b.hklcheckformat || "");
      mutateValidation();
    }
  }, [
    suitable,
    selectedBlock,
    setSelectedBlock,
    setSelectedInfo,
    setSelectedColumns,
    setSelectedDetails,
    mutateValidation,
  ]);

  const handleSelect = useCallback(
    async (block: RBlockInfo) => {
      if (!setSelectedBlock) return;
      await setSelectedBlock(block.bname);
      if (setSelectedInfo) await setSelectedInfo(block.info || "");
      if (setSelectedColumns)
        await setSelectedColumns(block.columnSetsText?.join(", ") || "");
      if (setSelectedDetails)
        await setSelectedDetails(block.hklcheckformat || "");
      await mutateValidation();
    },
    [
      setSelectedBlock,
      setSelectedInfo,
      setSelectedColumns,
      setSelectedDetails,
      mutateValidation,
    ]
  );

  if (!isMMCIF || !blocks.length) return null;

  return (
    <Box sx={{ mt: 1 }}>
      <Typography variant="body2" sx={{ fontWeight: 600, mb: 0.5 }}>
        Select mmCIF unmerged reflection block to use
        {suitable.length > 0 && " (first is default)"}
      </Typography>

      {suitable.length === 0 && (
        <Alert severity="error" sx={{ mb: 1 }}>
          No unmerged reflection blocks suitable for merging were found in
          this file.
        </Alert>
      )}

      {suitable.length > 0 && (
        <List dense disablePadding>
          {suitable.map((b, i) => (
            <ListItemButton
              key={`s-${b.bname}-${i}`}
              selected={selectedBlock === b.bname}
              onClick={() => handleSelect(b)}
              sx={{ py: 0.5 }}
            >
              <ListItemIcon sx={{ minWidth: 36 }}>
                <Radio
                  checked={selectedBlock === b.bname}
                  size="small"
                  tabIndex={-1}
                />
              </ListItemIcon>
              <ListItemText
                primary={b.bname}
                secondary={formatBlockDescription(b)}
                secondaryTypographyProps={{
                  component: "pre",
                  sx: { fontFamily: "inherit", whiteSpace: "pre-wrap", m: 0 },
                }}
              />
            </ListItemButton>
          ))}
        </List>
      )}

      {unsuitable.length > 0 && (
        <>
          <Typography
            variant="body2"
            color="text.secondary"
            sx={{ mt: 1, fontWeight: 600 }}
          >
            Other blocks (not suitable for merging)
          </Typography>
          <List dense disablePadding>
            {unsuitable.map((b, i) => (
              <ListItemButton
                key={`u-${b.bname}-${i}`}
                disabled
                sx={{ py: 0.5 }}
              >
                <ListItemText
                  primary={b.bname}
                  secondary={formatBlockDescription(b)}
                  primaryTypographyProps={{ color: "text.secondary" }}
                  secondaryTypographyProps={{
                    component: "pre",
                    sx: { fontFamily: "inherit", whiteSpace: "pre-wrap", m: 0 },
                  }}
                />
              </ListItemButton>
            ))}
          </List>
        </>
      )}
    </Box>
  );
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
        <MmcifBlockSelector jobId={props.job.id} />
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

      {/* Reference data — only relevant when MODE === "MATCH" */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel:
            "Reference data to match indexing and space group against (also later option to scale against this reference set)",
        }}
        containerHint="FolderLevel"
        visibility={() => mode === "MATCH"}
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
