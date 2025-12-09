import { CCP4i2TaskInterfaceProps } from "./task-container";
import {
  CCP4i2TaskElement,
  CCP4i2TaskElementProps,
} from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { doRetrieve, useApi } from "../../../api";
import { useJob, usePrevious } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback, useEffect, useMemo, useState } from "react";
import { showMtzColumnDialog, parseMtzColumns } from "../task-elements/mtz-column-dialog";
import { Job } from "../../../types/models";
import { useCCP4i2Window } from "../../../app-context";
import {
  Alert,
  Box,
  Card,
  CardContent,
  CardHeader,
  Checkbox,
  Chip,
  FormControlLabel,
  Grid2,
  List,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Paper,
  Radio,
  RadioGroup,
  Stack,
  Typography,
} from "@mui/material";
import TableChartIcon from "@mui/icons-material/TableChart";
import WarningIcon from "@mui/icons-material/Warning";
import CheckCircleIcon from "@mui/icons-material/CheckCircle";

// Type definitions matching digest output
interface MtzColumn {
  columnLabel: string;
  columnType: string;
  dataset: string;
  groupIndex?: number;
}

interface ColumnGroup {
  columnGroupType: string;
  contentFlag: number;
  dataset: string;
  columnList: MtzColumn[];
}

interface RBlockInfo {
  bname: string;
  entry_id: string;
  cell: number[];
  spacegroup_name: string;
  wavelength: number;
  highres: number;
  info: string;
  hklcheckformat: string;
  hasFreeR: boolean;
  freerValid: boolean;
  freerWarnings: string[];
  typeCodes?: number[];
  columnSetsText?: string[];
  columnnames?: Record<string, string[]>;
}

interface GenericReflDigest {
  format: string;
  merged: boolean;
  spaceGroup?: string;
  cell?: { a: number; b: number; c: number; alpha: number; beta: number; gamma: number };
  wavelength?: number;
  wavelengths?: number[];
  datasets?: string[];
  crystalNames?: string[];
  listOfColumns?: MtzColumn[];
  columnGroups?: ColumnGroup[];
  rblock_infos?: RBlockInfo[];
  hasFreeR: boolean;
  freerValid: boolean;
  freerWarnings: string[];
  freerColumnLabel?: string;
}

/**
 * Pattern definitions for column type signatures (same as splitMtz.tsx).
 * contentFlag values match backend CONTENT_SIGNATURE_LIST order:
 * - Obs: 1=Anom I (KMKM), 2=Anom SF (GLGL), 3=Mean I (JQ), 4=Mean SF (FQ)
 */
const COLUMN_PATTERNS: [string, string, number][] = [
  ["KMKM", "Obs", 1], // I+, sigI+, I-, sigI- (Anomalous Intensities)
  ["GLGL", "Obs", 2], // F+, sigF+, F-, sigF- (Anomalous SFs)
  ["JQ", "Obs", 3], // I, sigI (Mean Intensities)
  ["FQ", "Obs", 4], // F, sigF (Mean SFs)
  ["I", "FreeR", 1], // Integer flag
];

const OBS_CONTENT_LABELS: Record<number, string> = {
  1: "Anomalous Intensities (I+, σI+, I-, σI-)",
  2: "Anomalous Structure Factors (F+, σF+, F-, σF-)",
  3: "Mean Intensities (I, σI)",
  4: "Mean Structure Factors (F, σF)",
};

/**
 * mmCIF type code labels (from mmcifutils.py CIFLabelSets.ACCEPTED_SETS)
 * typeCode > 0 = observation data, = 0 FreeR, < 0 other
 */
const MMCIF_TYPE_LABELS: Record<number, string> = {
  1: "Mean Intensities (Imean)",
  2: "Anomalous Intensities (I±)",
  3: "Mean Structure Factors (Fmean)",
  4: "Anomalous Structure Factors (F±)",
  0: "FreeR Flag",
  [-1]: "Anomalous Difference (ΔF)",
  [-2]: "Phases (HL)",
  [-3]: "Map Coefficients",
};

/**
 * Group columns using pattern matching (simplified from splitMtz.tsx)
 */
function groupColumnsByPattern(columns: MtzColumn[]): ColumnGroup[] {
  if (!columns || columns.length === 0) return [];

  const groups: ColumnGroup[] = [];
  const used = new Array(columns.length).fill(false);
  const typeString = columns.map((c) => c.columnType).join("");

  let i = 0;
  while (i < columns.length) {
    if (used[i]) {
      i++;
      continue;
    }

    let matched = false;

    for (const [pattern, groupType, contentFlag] of COLUMN_PATTERNS) {
      const patternLen = pattern.length;
      if (i + patternLen > columns.length) continue;

      const candidate = typeString.slice(i, i + patternLen);
      if (candidate !== pattern) continue;

      // Special check for FreeR: label must contain 'free'
      if (groupType === "FreeR") {
        const label = columns[i].columnLabel;
        if (!label.toLowerCase().includes("free")) continue;
      }

      // For Obs data, columns must come from the same dataset
      if (groupType === "Obs") {
        const firstDataset = columns[i].dataset;
        if (firstDataset) {
          const sameDataset = Array.from({ length: patternLen }).every(
            (_, j) => columns[i + j].dataset === firstDataset
          );
          if (!sameDataset) continue;
        }
      }

      groups.push({
        columnGroupType: groupType,
        contentFlag,
        dataset: columns[i].dataset || "",
        columnList: columns.slice(i, i + patternLen),
      });

      for (let j = 0; j < patternLen; j++) {
        used[i + j] = true;
      }

      i += patternLen;
      matched = true;
      break;
    }

    if (!matched) {
      i++;
    }
  }

  return groups;
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { job } = props;
  const { useFileDigest, useTaskItem, mutateContainer, mutateValidation, uploadFileParam } =
    useJob(job.id);
  const { cootModule } = useCCP4i2Window();

  const { item: HKLINItem, value: HKLINValue } = useTaskItem("HKLIN");
  const oldHKLINValue = usePrevious(HKLINValue);

  // Use task-qualified path for digest API - only fetch when file has been uploaded (has dbFileId)
  const hasUploadedFile = Boolean(HKLINValue?.dbFileId);
  const digestObjectPath = hasUploadedFile ? "import_merged.inputData.HKLIN" : "";
  const { data: HKLINDigest, isLoading: digestLoading, error: digestError } = useFileDigest(digestObjectPath) as {
    data: GenericReflDigest | null;
    isLoading: boolean;
    error: Error | null;
  };

  // Debug: log digest fetching status
  useEffect(() => {
    console.log("[import_merged] Digest status:", {
      hasUploadedFile,
      digestObjectPath,
      HKLINValue,
      HKLINDigest,
      digestLoading,
      digestError,
    });
  }, [hasUploadedFile, digestObjectPath, HKLINValue, HKLINDigest, digestLoading, digestError]);

  const { item: HKLIN_OBSItem } = useTaskItem("HKLIN_OBS");
  const { update: updateSPACEGROUP } = useTaskItem("SPACEGROUP");
  const { update: updateUNITCELL } = useTaskItem("UNITCELL");
  const { update: updateWAVELENGTH } = useTaskItem("WAVELENGTH");
  const { update: updateCRYSTALNAME } = useTaskItem("CRYSTALNAME");
  const { update: updateDATASETNAME } = useTaskItem("DATASETNAME");
  const { update: setHKLIN_OBS_COLUMNS } = useTaskItem("HKLIN_OBS_COLUMNS");
  const { update: setHKLIN_OBS_CONTENT_FLAG } = useTaskItem("HKLIN_OBS_CONTENT_FLAG");
  const { value: HKLIN_FORMATValue, update: updateHKLIN_FORMAT } = useTaskItem("HKLIN_FORMAT");

  // mmCIF specific
  const { update: setMMCIF_SELECTED_BLOCK } = useTaskItem("MMCIF_SELECTED_BLOCK");
  const { update: setMMCIF_SELECTED_COLUMNS } = useTaskItem("MMCIF_SELECTED_COLUMNS");
  const { update: setMMCIF_SELECTED_ISINTENSITY } = useTaskItem("MMCIF_SELECTED_ISINTENSITY");
  const { update: setMMCIF_SELECTED_INFO } = useTaskItem("MMCIF_SELECTED_INFO");
  const { update: setMMCIF_SELECTED_CONTENT } = useTaskItem("MMCIF_SELECTED_CONTENT");
  const { update: setHASFREER } = useTaskItem("HASFREER");

  // Local state for UI
  const [selectedObsGroup, setSelectedObsGroup] = useState<ColumnGroup | null>(null);
  const [selectedMmcifBlock, setSelectedMmcifBlock] = useState<string | null>(null);
  const [selectedMmcifContentType, setSelectedMmcifContentType] = useState<number | null>(null);

  // Compute column groups from digest
  const columnGroups = useMemo(() => {
    if (!HKLINDigest) return [];

    // Use pre-computed column groups if available
    if (HKLINDigest.columnGroups && HKLINDigest.columnGroups.length > 0) {
      return HKLINDigest.columnGroups;
    }

    // Otherwise compute from listOfColumns
    if (HKLINDigest.listOfColumns) {
      return groupColumnsByPattern(HKLINDigest.listOfColumns);
    }

    return [];
  }, [HKLINDigest]);

  // Filter observation groups (exclude FreeR, phases, etc.)
  const obsGroups = useMemo(() => {
    return columnGroups.filter((g) => g.columnGroupType === "Obs");
  }, [columnGroups]);

  // Find FreeR groups
  const freerGroups = useMemo(() => {
    return columnGroups.filter((g) => g.columnGroupType === "FreeR");
  }, [columnGroups]);

  // Track if we've already processed the current digest to prevent re-processing
  const [processedDigestKey, setProcessedDigestKey] = useState<string | null>(null);

  // Effect: Handle digest changes - auto-populate metadata
  // Note: We intentionally exclude update functions from dependencies to prevent infinite loops.
  // The update functions change on every render, but we only want to run this effect
  // when the digest data actually changes.
  useEffect(() => {
    const digestKey = HKLINDigest ? JSON.stringify({
      spaceGroup: HKLINDigest.spaceGroup,
      wavelength: HKLINDigest.wavelength,
      format: HKLINDigest.format,
      cell: HKLINDigest.cell,
      hasFreeR: HKLINDigest.hasFreeR,
    }) : null;

    // Skip if we've already processed this digest or if nothing to process
    if (!digestKey || digestKey === processedDigestKey || job?.status !== 1) {
      return;
    }

    console.log("[import_merged] Processing new digest:", HKLINDigest);
    console.log("[import_merged] Update functions available:", {
      updateSPACEGROUP: !!updateSPACEGROUP,
      updateUNITCELL: !!updateUNITCELL,
      updateWAVELENGTH: !!updateWAVELENGTH,
      updateHKLIN_FORMAT: !!updateHKLIN_FORMAT,
      updateCRYSTALNAME: !!updateCRYSTALNAME,
      updateDATASETNAME: !!updateDATASETNAME,
    });

    const processDigest = async () => {
      if (!HKLINDigest) return;

      let parametersChanged = false;

      // Update space group
      if (HKLINDigest.spaceGroup) {
        const cleanedSG = String(HKLINDigest.spaceGroup).replace(/\s+/g, "");
        console.log("[import_merged] Setting spaceGroup to:", cleanedSG);
        if (updateSPACEGROUP) {
          try {
            const result = await updateSPACEGROUP(cleanedSG);
            console.log("[import_merged] updateSPACEGROUP result:", result);
            parametersChanged = parametersChanged || Boolean(result);
          } catch (e) {
            console.error("[import_merged] updateSPACEGROUP error:", e);
          }
        } else {
          console.warn("[import_merged] updateSPACEGROUP is not available!");
        }
      }

      // Update wavelength
      if (HKLINDigest.wavelength) {
        console.log("[import_merged] Setting wavelength to:", HKLINDigest.wavelength);
        if (updateWAVELENGTH) {
          try {
            const result = await updateWAVELENGTH(HKLINDigest.wavelength);
            console.log("[import_merged] updateWAVELENGTH result:", result);
            parametersChanged = parametersChanged || Boolean(result);
          } catch (e) {
            console.error("[import_merged] updateWAVELENGTH error:", e);
          }
        }
      }

      // Update format
      if (HKLINDigest.format) {
        console.log("[import_merged] Setting format to:", HKLINDigest.format.toUpperCase());
        if (updateHKLIN_FORMAT) {
          try {
            const result = await updateHKLIN_FORMAT(HKLINDigest.format.toUpperCase());
            console.log("[import_merged] updateHKLIN_FORMAT result:", result);
            parametersChanged = parametersChanged || Boolean(result);
          } catch (e) {
            console.error("[import_merged] updateHKLIN_FORMAT error:", e);
          }
        }
      }

      // Update unit cell
      if (HKLINDigest.cell) {
        console.log("[import_merged] Setting cell to:", HKLINDigest.cell);
        if (updateUNITCELL) {
          try {
            const result = await updateUNITCELL(HKLINDigest.cell);
            console.log("[import_merged] updateUNITCELL result:", result);
            parametersChanged = parametersChanged || Boolean(result);
          } catch (e) {
            console.error("[import_merged] updateUNITCELL error:", e);
          }
        }
      }

      // Update crystal name and dataset name from digest
      if (HKLINDigest.crystalNames && HKLINDigest.crystalNames.length > 0 ) {
        console.log("[import_merged] Setting crystalName to:", HKLINDigest.crystalNames);
        if (updateCRYSTALNAME) {
          try {
            const result = await updateCRYSTALNAME(HKLINDigest.crystalNames[0]);
            console.log("[import_merged] updateCRYSTALNAME result:", result);
            parametersChanged = parametersChanged || Boolean(result);
          } catch (e) {
            console.error("[import_merged] updateCRYSTALNAME error:", e);
          }
        }
      }

      if (HKLINDigest.datasets && HKLINDigest.datasets.length > 0 ) {
        console.log("[import_merged] Setting datasetName to:", HKLINDigest.datasets[0]);
        if (updateDATASETNAME) {
          try {
            const result = await updateDATASETNAME(HKLINDigest.datasets[0]);
            console.log("[import_merged] updateDATASETNAME result:", result);
            parametersChanged = parametersChanged || Boolean(result);
          } catch (e) {
            console.error("[import_merged] updateDATASETNAME error:", e);
          }
        }
      }

      // Update HASFREER
      if (setHASFREER) {
        await setHASFREER(Boolean(HKLINDigest.hasFreeR && HKLINDigest.freerValid));
      }

      console.log("[import_merged] parametersChanged:", parametersChanged);
      if (parametersChanged) {
        await mutateContainer();
        await mutateValidation();
      }

      // Mark this digest as processed
      setProcessedDigestKey(digestKey);
    };

    processDigest();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [HKLINDigest, job?.status, processedDigestKey]);

  // Handle observation group selection (MTZ)
  const handleObsGroupSelect = useCallback(
    async (group: ColumnGroup) => {
      if (!setHKLIN_OBS_CONTENT_FLAG || !setHKLIN_OBS_COLUMNS || job?.status !== 1)
        return;

      setSelectedObsGroup(group);

      const columnLabels = group.columnList.map((c) => c.columnLabel).join(",");
      await setHKLIN_OBS_COLUMNS(columnLabels);
      await setHKLIN_OBS_CONTENT_FLAG(group.contentFlag);

      // Update dataset and crystal name
      if (group.dataset && updateDATASETNAME) {
        await updateDATASETNAME(group.dataset);
      }

      // Try to find crystal name from digest
      if (HKLINDigest?.crystalNames && HKLINDigest.datasets && updateCRYSTALNAME) {
        const datasetIndex = HKLINDigest.datasets.indexOf(group.dataset);
        if (datasetIndex >= 0 && datasetIndex < HKLINDigest.crystalNames.length) {
          await updateCRYSTALNAME(HKLINDigest.crystalNames[datasetIndex]);
        }
      }

      await mutateContainer();
      await mutateValidation();
    },
    [
      job?.status,
      setHKLIN_OBS_COLUMNS,
      setHKLIN_OBS_CONTENT_FLAG,
      updateDATASETNAME,
      updateCRYSTALNAME,
      HKLINDigest,
      mutateContainer,
      mutateValidation,
    ]
  );

  // Handle mmCIF block selection
  const handleMmcifBlockSelect = useCallback(
    async (blockName: string) => {
      if (!HKLINDigest?.rblock_infos || job?.status !== 1) return;

      const block = HKLINDigest.rblock_infos.find((b) => b.bname === blockName);
      if (!block) return;

      setSelectedMmcifBlock(blockName);

      // Update parameters from block
      if (setMMCIF_SELECTED_BLOCK) await setMMCIF_SELECTED_BLOCK(blockName);
      if (updateDATASETNAME) await updateDATASETNAME(blockName);
      if (updateCRYSTALNAME) await updateCRYSTALNAME(blockName);
      if (updateSPACEGROUP) await updateSPACEGROUP(block.spacegroup_name);
      if (updateWAVELENGTH && block.wavelength) await updateWAVELENGTH(block.wavelength);
      if (updateUNITCELL && block.cell) {
        await updateUNITCELL({
          a: block.cell[0],
          b: block.cell[1],
          c: block.cell[2],
          alpha: block.cell[3],
          beta: block.cell[4],
          gamma: block.cell[5],
        });
      }

      // Set info
      if (setMMCIF_SELECTED_INFO) {
        const info = `${block.info}\nhkl list type: ${block.hklcheckformat}\nHighest resolution: ${block.highres}`;
        await setMMCIF_SELECTED_INFO(info);
      }

      // Update HASFREER based on block
      if (setHASFREER) {
        await setHASFREER(block.hasFreeR && block.freerValid);
      }

      // Find the best observation data type
      if (block.typeCodes && block.typeCodes.length > 0) {
        // Find first positive type code (observation data)
        const obsTypeCode = block.typeCodes.find((tc) => tc > 0);
        if (obsTypeCode !== undefined) {
          setSelectedMmcifContentType(obsTypeCode);
          await handleMmcifContentSelect(block, obsTypeCode);
        }
      }

      await mutateContainer();
      await mutateValidation();
    },
    [
      HKLINDigest,
      job?.status,
      setMMCIF_SELECTED_BLOCK,
      setMMCIF_SELECTED_INFO,
      updateDATASETNAME,
      updateCRYSTALNAME,
      updateSPACEGROUP,
      updateWAVELENGTH,
      updateUNITCELL,
      setHASFREER,
      mutateContainer,
      mutateValidation,
    ]
  );

  // Handle mmCIF content type selection within a block
  const handleMmcifContentSelect = useCallback(
    async (block: RBlockInfo, typeCode: number) => {
      if (job?.status !== 1) return;

      setSelectedMmcifContentType(typeCode);

      // Set content type label
      const contentLabel = block.columnSetsText?.[block.typeCodes?.indexOf(typeCode) ?? 0] || "";
      if (setMMCIF_SELECTED_CONTENT) await setMMCIF_SELECTED_CONTENT(contentLabel);

      // Set columns based on content type
      if (block.columnnames && setMMCIF_SELECTED_COLUMNS) {
        const columns = block.columnnames[contentLabel];
        if (columns) {
          await setMMCIF_SELECTED_COLUMNS(columns.join(", "));
        }
      }

      // Set isIntensity: +1 for I types (1, 2), -1 for F types (3, 4)
      if (setMMCIF_SELECTED_ISINTENSITY) {
        const isIntensity = typeCode === 1 || typeCode === 2 ? 1 : -1;
        await setMMCIF_SELECTED_ISINTENSITY(isIntensity);
      }

      await mutateContainer();
      await mutateValidation();
    },
    [job?.status, setMMCIF_SELECTED_CONTENT, setMMCIF_SELECTED_COLUMNS, setMMCIF_SELECTED_ISINTENSITY, mutateContainer, mutateValidation]
  );

  // Process column selection from MTZ dialog (legacy path)
  const processColumnSelection = useCallback(
    async (columnPath: string, file: File) => {
      if (!setHKLIN_OBS_CONTENT_FLAG || !setHKLIN_OBS_COLUMNS) return;

      const match = columnPath.match(/\[([^\]]+)\]/);
      if (match) {
        await setHKLIN_OBS_COLUMNS(match[1]);
        const columnNames = match[1].split(",").map((name) => name.trim());

        // Determine content flag from column types
        if (HKLINDigest?.listOfColumns) {
          const columnTypes = columnNames.map(
            (name) =>
              HKLINDigest.listOfColumns?.find(
                (col) => col.columnLabel === name
              )?.columnType
          );
          const signature = columnTypes.join("");
          const contentFlag = ["KMKM", "GLGL", "JQ", "FQ"].indexOf(signature);
          if (contentFlag > -1) {
            await setHKLIN_OBS_CONTENT_FLAG(contentFlag + 1);
          }
        }
      }

      // Upload the file
      if (columnPath && columnPath.trim().length > 0 && HKLIN_OBSItem) {
        await uploadFileParam({
          objectPath: HKLIN_OBSItem._objectPath,
          file: file,
          fileName: file.name,
          columnSelector: columnPath,
        });
      }
    },
    [HKLINDigest, HKLIN_OBSItem, setHKLIN_OBS_COLUMNS, setHKLIN_OBS_CONTENT_FLAG, uploadFileParam]
  );

  // Handle HKLIN file change (trigger column dialog for MTZ)
  const handleHKLINFileChange = useCallback(
    async (hklinValue: any) => {
      if (
        !hklinValue?.dbFileId ||
        !hklinValue?.baseName ||
        !oldHKLINValue ||
        !cootModule ||
        job?.status !== 1
      )
        return;
      if (JSON.stringify(hklinValue) === JSON.stringify(oldHKLINValue)) return;

      const isMtzFile = hklinValue.baseName.toLowerCase().endsWith(".mtz");
      if (!isMtzFile) return;

      // Download and parse MTZ
      const downloadURL = `files/${hklinValue.dbFileId}/download_by_uuid`;
      const arrayBuffer = await doRetrieve(downloadURL, hklinValue.baseName);
      const blob = new Blob([arrayBuffer], { type: "application/CCP4-mtz-file" });
      const file = new File([blob], hklinValue.baseName, { type: "application/CCP4-mtz-file" });

      const columnNames = await parseMtzColumns(file, cootModule);
      if (!columnNames) return;

      const columnPath = await showMtzColumnDialog(columnNames, HKLIN_OBSItem);
      if (!columnPath) return;

      await processColumnSelection(columnPath, file);
    },
    [oldHKLINValue, cootModule, job?.status, HKLIN_OBSItem, processColumnSelection]
  );

  // Effect: Handle HKLIN value changes
  useEffect(() => {
    if (HKLINValue) {
      handleHKLINFileChange(HKLINValue);
    }
  }, [HKLINValue, handleHKLINFileChange]);

  return (
    <>
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Main inputs" key="1">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input data", initiallyOpen: true }}
            key="Input data"
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              key="HKLIN"
              itemName="HKLIN"
              qualifiers={{ guiLabel: "Reflections" }}
            />

            {/* Metadata Card */}
            <Card sx={{ mb: 2 }}>
              <CardHeader title="Crystal Information" />
              <CardContent>
                <Grid2 container direction="row" sx={{ mb: 2 }} spacing={2}>
                  <Grid2 size={{ xs: 6, md: 4 }}>
                    <CCP4i2TaskElement
                      {...props}
                      key="SPACEGROUP"
                      itemName="SPACEGROUP"
                      qualifiers={{ guiLabel: "Space group" }}
                    />
                  </Grid2>
                  <Grid2 size={{ xs: 6, md: 8 }}>
                    <CCP4i2TaskElement {...props} key="UNITCELL" itemName="UNITCELL" />
                  </Grid2>
                  <Stack direction="row" spacing={2}>
                    <CCP4i2TaskElement
                      {...props}
                      key="CRYSTALNAME"
                      itemName="CRYSTALNAME"
                      qualifiers={{ guiLabel: "Crystal name" }}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      key="DATASETNAME"
                      itemName="DATASETNAME"
                      qualifiers={{ guiLabel: "Dataset name" }}
                    />
                  </Stack>
                </Grid2>
                <CCP4i2TaskElement
                  {...props}
                  key="WAVELENGTH"
                  itemName="WAVELENGTH"
                  qualifiers={{ guiLabel: "Wavelength" }}
                />
              </CardContent>
            </Card>

            {/* Loading state while fetching digest */}
            {hasUploadedFile && digestLoading && (
              <Card sx={{ mb: 2 }}>
                <CardContent>
                  <Typography color="text.secondary">
                    Loading reflection file information...
                  </Typography>
                </CardContent>
              </Card>
            )}

            {/* Error state */}
            {hasUploadedFile && digestError && (
              <Alert severity="error" sx={{ mb: 2 }}>
                Error loading reflection file: {digestError.message}
              </Alert>
            )}

            {/* Format-specific panels - use digest format directly, not the task parameter */}
            {HKLINDigest?.format?.toUpperCase() === "MTZ" && (
              <MtzReflectionPanel
                {...props}
                digest={HKLINDigest}
                obsGroups={obsGroups}
                freerGroups={freerGroups}
                selectedObsGroup={selectedObsGroup}
                onObsGroupSelect={handleObsGroupSelect}
              />
            )}

            {HKLINDigest?.format?.toUpperCase() === "MMCIF" && (
              <MmcifReflectionPanel
                {...props}
                digest={HKLINDigest}
                selectedBlock={selectedMmcifBlock}
                selectedContentType={selectedMmcifContentType}
                onBlockSelect={handleMmcifBlockSelect}
                onContentSelect={handleMmcifContentSelect}
              />
            )}

            {/* For other formats (Scalepack, XDS) show basic info from digest */}
            {HKLINDigest && !["MTZ", "MMCIF"].includes(HKLINDigest.format?.toUpperCase() || "") && (
              <Card sx={{ mb: 2 }}>
                <CardHeader title={`${HKLINDigest.format} Reflection Data`} />
                <CardContent>
                  <Typography variant="body2" color="text.secondary">
                    Format: {HKLINDigest.format}
                    {HKLINDigest.merged !== undefined && (
                      <> | {HKLINDigest.merged ? "Merged" : "Unmerged"}</>
                    )}
                  </Typography>
                  <FreeRStatusDisplay
                    hasFreeR={HKLINDigest.hasFreeR}
                    freerValid={HKLINDigest.freerValid}
                    freerWarnings={HKLINDigest.freerWarnings || []}
                  />
                </CardContent>
              </Card>
            )}

            {/* FreeR input */}
            <CCP4i2TaskElement {...props} itemName="FREERFLAG" />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </>
  );
};

// ============================================================================
// MTZ Reflection Panel
// ============================================================================

interface MtzReflectionPanelProps {
  digest: GenericReflDigest;
  obsGroups: ColumnGroup[];
  freerGroups: ColumnGroup[];
  selectedObsGroup: ColumnGroup | null;
  onObsGroupSelect: (group: ColumnGroup) => void;
}

const MtzReflectionPanel: React.FC<MtzReflectionPanelProps> = ({
  digest,
  obsGroups,
  freerGroups,
  selectedObsGroup,
  onObsGroupSelect,
}) => {
  return (
    <Card sx={{ mb: 2 }}>
      <CardHeader title="MTZ Reflection Data" />
      <CardContent>
        {/* Observation Data Selection */}
        <Typography variant="subtitle2" gutterBottom>
          Select Observation Data
        </Typography>

        {obsGroups.length === 0 ? (
          <Alert severity="warning" sx={{ mb: 2 }}>
            No observation data groups found in this MTZ file
          </Alert>
        ) : (
          <List dense sx={{ mb: 2 }}>
            {obsGroups.map((group, idx) => {
              const isSelected = Boolean(
                selectedObsGroup &&
                group.columnList.map((c) => c.columnLabel).join(",") ===
                  selectedObsGroup.columnList.map((c) => c.columnLabel).join(",")
              );
              const labels = group.columnList.map((c) => c.columnLabel).join(", ");

              return (
                <ListItemButton
                  key={idx}
                  selected={isSelected}
                  onClick={() => onObsGroupSelect(group)}
                  sx={{
                    border: 1,
                    borderColor: isSelected ? "primary.main" : "divider",
                    borderRadius: 1,
                    mb: 0.5,
                  }}
                >
                  <ListItemIcon>
                    <TableChartIcon color={isSelected ? "primary" : "inherit"} />
                  </ListItemIcon>
                  <ListItemText
                    primary={
                      <Stack direction="row" spacing={1} alignItems="center">
                        <Typography variant="body2">{labels}</Typography>
                        <Chip
                          label={OBS_CONTENT_LABELS[group.contentFlag] || `Type ${group.contentFlag}`}
                          size="small"
                          color={isSelected ? "primary" : "default"}
                        />
                      </Stack>
                    }
                    secondary={`Dataset: ${group.dataset || "default"}`}
                  />
                </ListItemButton>
              );
            })}
          </List>
        )}

        {/* FreeR Status */}
        <FreeRStatusDisplay
          hasFreeR={digest.hasFreeR}
          freerValid={digest.freerValid}
          freerWarnings={digest.freerWarnings}
          freerColumnLabel={digest.freerColumnLabel}
        />

        {/* Show selected columns */}
        {selectedObsGroup && (
          <Box sx={{ mt: 2, p: 1, bgcolor: "action.hover", borderRadius: 1 }}>
            <Typography variant="caption" color="text.secondary">
              Selected: {selectedObsGroup.columnList.map((c) => c.columnLabel).join(", ")}
            </Typography>
          </Box>
        )}
      </CardContent>
    </Card>
  );
};

// ============================================================================
// mmCIF Reflection Panel
// ============================================================================

interface MmcifReflectionPanelProps {
  digest: GenericReflDigest;
  selectedBlock: string | null;
  selectedContentType: number | null;
  onBlockSelect: (blockName: string) => void;
  onContentSelect: (block: RBlockInfo, typeCode: number) => void;
}

const MmcifReflectionPanel: React.FC<MmcifReflectionPanelProps> = ({
  digest,
  selectedBlock,
  selectedContentType,
  onBlockSelect,
  onContentSelect,
}) => {
  const rblocks = digest.rblock_infos || [];

  // Filter to blocks with merged diffraction data (positive typeCodes)
  const validBlocks = useMemo(() => {
    return rblocks.filter((block) => {
      if (!block.typeCodes) return false;
      return block.typeCodes.some((tc) => tc > 0);
    });
  }, [rblocks]);

  const currentBlock = selectedBlock
    ? rblocks.find((b) => b.bname === selectedBlock)
    : null;

  return (
    <Card sx={{ mb: 2 }}>
      <CardHeader title="mmCIF Reflection Data" />
      <CardContent>
        {validBlocks.length === 0 ? (
          <Alert severity="warning">
            No suitable reflection blocks found in this mmCIF file
          </Alert>
        ) : (
          <>
            {/* Block Selection */}
            <Typography variant="subtitle2" gutterBottom>
              {validBlocks.length === 1
                ? "Reflection Block"
                : "Select Reflection Block"}
            </Typography>

            <List dense sx={{ mb: 2 }}>
              {validBlocks.map((block, idx) => {
                const isSelected = selectedBlock === block.bname;

                return (
                  <ListItemButton
                    key={idx}
                    selected={isSelected}
                    onClick={() => onBlockSelect(block.bname)}
                    sx={{
                      border: 1,
                      borderColor: isSelected ? "primary.main" : "divider",
                      borderRadius: 1,
                      mb: 0.5,
                    }}
                  >
                    <ListItemIcon>
                      <TableChartIcon color={isSelected ? "primary" : "inherit"} />
                    </ListItemIcon>
                    <ListItemText
                      primary={
                        <Stack direction="row" spacing={1} alignItems="center">
                          <Typography variant="body2" fontWeight="medium">
                            {block.bname}
                          </Typography>
                          <Chip
                            label={`${block.highres.toFixed(2)}Å`}
                            size="small"
                            variant="outlined"
                          />
                        </Stack>
                      }
                      secondary={
                        <Typography variant="caption" component="div">
                          {block.info}
                          <br />
                          {block.hklcheckformat}
                        </Typography>
                      }
                    />
                  </ListItemButton>
                );
              })}
            </List>

            {/* Content Type Selection within selected block */}
            {currentBlock && currentBlock.typeCodes && (
              <>
                <Typography variant="subtitle2" gutterBottom>
                  Data Type
                </Typography>
                <RadioGroup
                  value={selectedContentType ?? ""}
                  onChange={(e) => {
                    const typeCode = parseInt(e.target.value, 10);
                    onContentSelect(currentBlock, typeCode);
                  }}
                >
                  {currentBlock.typeCodes
                    .filter((tc) => tc > 0) // Only show observation data types
                    .map((typeCode, idx) => (
                      <FormControlLabel
                        key={idx}
                        value={typeCode}
                        control={<Radio size="small" />}
                        label={
                          <Typography variant="body2">
                            {currentBlock.columnSetsText?.[
                              currentBlock.typeCodes?.indexOf(typeCode) ?? idx
                            ] || MMCIF_TYPE_LABELS[typeCode] || `Type ${typeCode}`}
                          </Typography>
                        }
                      />
                    ))}
                </RadioGroup>

                {/* FreeR Status */}
                <Box sx={{ mt: 2 }}>
                  <FreeRStatusDisplay
                    hasFreeR={currentBlock.hasFreeR}
                    freerValid={currentBlock.freerValid}
                    freerWarnings={currentBlock.freerWarnings}
                  />
                </Box>
              </>
            )}
          </>
        )}
      </CardContent>
    </Card>
  );
};

// ============================================================================
// FreeR Status Display Component
// ============================================================================

interface FreeRStatusDisplayProps {
  hasFreeR: boolean;
  freerValid: boolean;
  freerWarnings: string[];
  freerColumnLabel?: string;
}

const FreeRStatusDisplay: React.FC<FreeRStatusDisplayProps> = ({
  hasFreeR,
  freerValid,
  freerWarnings,
  freerColumnLabel,
}) => {
  if (!hasFreeR) {
    return (
      <Alert severity="info" sx={{ mt: 1 }}>
        No FreeR flag found in input data. A new FreeR set will be generated.
      </Alert>
    );
  }

  if (!freerValid) {
    return (
      <Alert severity="warning" icon={<WarningIcon />} sx={{ mt: 1 }}>
        <Typography variant="body2" fontWeight="medium">
          FreeR flag found but may be invalid
        </Typography>
        {freerWarnings && freerWarnings.length > 0 && (
          <Box component="ul" sx={{ m: 0, pl: 2 }}>
            {freerWarnings.map((warning, idx) => (
              <li key={idx}>
                <Typography variant="caption">{warning}</Typography>
              </li>
            ))}
          </Box>
        )}
      </Alert>
    );
  }

  return (
    <Alert severity="success" icon={<CheckCircleIcon />} sx={{ mt: 1 }}>
      <Typography variant="body2">
        Valid FreeR flag found{freerColumnLabel ? `: ${freerColumnLabel}` : ""}
      </Typography>
    </Alert>
  );
};

export default TaskInterface;
