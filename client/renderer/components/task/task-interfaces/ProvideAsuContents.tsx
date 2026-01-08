import {
  Alert,
  Box,
  Grid2,
  LinearProgress,
  Stack,
  Typography,
} from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob, useProjectFiles } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback, useMemo } from "react";
import { useApi } from "../../../api";
import { BaseSpacegroupCellElement } from "../task-elements/base-spacegroup-cell-element";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const api = useApi();
  const { useTaskItem, useFileDigest, fetchDigest, getErrors, mutateValidation } = useJob(job.id);
  const { files: projectFiles } = useProjectFiles(job.project);
  const { forceUpdate: forceSetAsuContent } = useTaskItem("ASU_CONTENT");
  const { item: asuContentInItem } = useTaskItem("ASUCONTENTIN");
  const { item: asuContentItem } = useTaskItem("ASU_CONTENT");
  const { value: HKLINValue } = useTaskItem("HKLIN");

  // Check if there are any existing CAsuDataFiles in the project
  // Hide the "load existing" section if none exist (prevents recursive create button)
  const hasExistingAsuFiles = useMemo(() => {
    if (!projectFiles) return false;
    return projectFiles.some((file) => file.type === "application/CCP4-asu-content");
  }, [projectFiles]);

  // File digest for HKLIN (used for Matthews calculation)
  // Only fetch when a file has been uploaded (has dbFileId) - otherwise digest endpoint fails
  const hasHKLINFile = Boolean(HKLINValue?.dbFileId);
  const hklinDigestPath = hasHKLINFile ? "ProvideAsuContents.inputData.HKLIN" : "";
  const { data: HKLINDigest, mutate: mutateHKLINDigest } = useFileDigest(hklinDigestPath);

  // ASU content is valid when there are no validation errors for it
  // Uses the existing validation infrastructure from useJob().getErrors()
  // which fetches from /api/jobs/{id}/validation/ endpoint
  const asuContentErrors = getErrors(asuContentItem);
  const isAsuContentValid = asuContentErrors.length === 0 && (asuContentItem?._value?.length ?? 0) > 0;

  /**
   * Fetches the molecular weight for the current job's ASU content using SWR.
   * Only fetches when ASU_CONTENT is valid (passes all validation checks).
   */
  const { data: molWeight, mutate: mutateMolWeight } = api.objectMethod<any>(
    job.id,
    "ProvideAsuContents.inputData.ASU_CONTENT",
    "molecularWeight",
    {},
    [],
    isAsuContentValid
  );

  /**
   * Fetches and caches the Matthews coefficient analysis for the current job using SWR.
   * Only fetches when we have both a valid molecular weight result AND an HKLIN file.
   */
  const { data: matthewsAnalysis, mutate: mutateMatthews } = api.objectMethod<any>(
    job.id,
    "ProvideAsuContents.inputData.HKLIN.fileContent",
    "matthewsCoeff",
    { molWt: molWeight?.data?.result },
    [molWeight?.data?.result, HKLINDigest],
    !!(molWeight?.data?.result && HKLINDigest)
  );

  /**
   * Handle ASUCONTENTIN file change - explicitly fetch digest and populate ASU_CONTENT.
   * Uses imperative fetchDigest for deterministic, race-condition-free behavior.
   *
   * @param updatedItem - The updated file item passed from onChange (has fresh _objectPath)
   */
  const handleAsuContentInChange = useCallback(async (updatedItem?: any) => {
    // Use the updated item's path if provided (from onChange), otherwise fall back to current state
    // This is important because when selecting from pulldown, the container hasn't mutated yet
    const objectPath = updatedItem?._objectPath || asuContentInItem?._objectPath;
    if (!objectPath) return;

    // Fetch the digest for the newly uploaded/selected file
    const digestData = await fetchDigest(objectPath);

    // Extract seqList and populate ASU_CONTENT
    if (digestData?.seqList && Array.isArray(digestData.seqList)) {
      const seqList = digestData.seqList.map((seq: any) => ({
        name: seq.name,
        sequence: seq.sequence,
        polymerType: seq.polymerType,
        description: seq.description,
        nCopies: seq.nCopies ?? 1,  // Default to 1 if not provided
      }));
      await forceSetAsuContent(seqList);
      // Refresh validation, molecular weight, and Matthews after updating ASU content
      // Must await molWeight mutation so the new value is available for Matthews calculation
      mutateValidation();
      await mutateMolWeight();
      mutateMatthews();
    }
  }, [asuContentInItem?._objectPath, fetchDigest, forceSetAsuContent, mutateValidation, mutateMolWeight, mutateMatthews]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        {/* Only show "load existing" option if there are existing ASU files in project */}
        {hasExistingAsuFiles && (
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Optionally load existing ASU content file to edit",
            }}
            containerHint="BlockLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="ASUCONTENTIN"
              qualifiers={{ guiLabel: "ASU contents" }}
              onChange={handleAsuContentInChange}
            />
          </CCP4i2ContainerElement>
        )}
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{
            guiLabel:
              "Specify the protein/nucleic acid sequences in the crystal",
          }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="ASU_CONTENT"
            qualifiers={{ guiLabel: "ASU contents" }}
            onChange={async () => {
              // Re-check validity, molecular weight, and Matthews when ASU content changes
              mutateValidation();
              await mutateMolWeight();
              mutateMatthews();
            }}
          />
        </CCP4i2ContainerElement>
        <Typography>
          Molecular weight:{" "}
          {molWeight?.data?.result
            ? molWeight.data.result?.toFixed(2)
            : isAsuContentValid
              ? "(calculating...)"
              : asuContentErrors.length > 0
                ? `(${asuContentErrors.length} validation error(s))`
                : "(no valid sequences)"}
        </Typography>
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Solvent analysis" }}
          containerHint="BlockLevel"
          initiallyOpen={true}
        >
          <Grid2 container spacing={2}>
            <Grid2 size={{ xs: 12, sm: 8 }}>
              <CCP4i2TaskElement
                {...props}
                itemName="HKLIN"
                qualifiers={{ guiLabel: "MTZFile (for Matthews volume calc)" }}
                onChange={() => mutateHKLINDigest()}
              />
              {/* Show MTZ file info when digest is available */}
              {HKLINDigest && (
                <Stack spacing={1} sx={{ mt: 1 }}>
                  <BaseSpacegroupCellElement data={HKLINDigest} />
                  {/* Warning if no wavelengths (e.g., FreeR-only file) */}
                  {HKLINDigest.cell && (!HKLINDigest.wavelengths || HKLINDigest.wavelengths.length === 0) && (
                    <Alert severity="info" sx={{ py: 0 }}>
                      No wavelength information in this MTZ file
                    </Alert>
                  )}
                </Stack>
              )}
            </Grid2>
            <Grid2 size={{ xs: 12, sm: 4 }}>
              {matthewsAnalysis?.success && matthewsAnalysis?.data?.result ? (
                <Box
                  sx={{
                    p: 1.5,
                    borderRadius: 1,
                    bgcolor: "action.hover",
                    border: 1,
                    borderColor: "divider",
                  }}
                >
                  <Typography
                    variant="caption"
                    color="text.secondary"
                    sx={{ mb: 1, display: "block" }}
                  >
                    Matthews Analysis
                  </Typography>
                  <Stack spacing={1}>
                    {matthewsAnalysis?.data?.result?.results.map(
                      (result: {
                        nmol_in_asu: number;
                        percent_solvent: number;
                        prob_matth: number;
                      }) => {
                        const probability = result.prob_matth;
                        const isLikely = probability > 0.5;
                        return (
                          <Box
                            key={result.nmol_in_asu}
                            sx={{
                              p: 1,
                              borderRadius: 0.5,
                              bgcolor: isLikely
                                ? "success.main"
                                : "background.paper",
                              color: isLikely
                                ? "success.contrastText"
                                : "text.primary",
                              border: 1,
                              borderColor: isLikely
                                ? "success.main"
                                : "divider",
                            }}
                          >
                            <Stack
                              direction="row"
                              justifyContent="space-between"
                              alignItems="center"
                            >
                              <Typography variant="body2" fontWeight="medium">
                                {result.nmol_in_asu} mol/ASU
                              </Typography>
                              <Typography variant="body2" fontWeight="bold">
                                {(probability * 100).toFixed(0)}%
                              </Typography>
                            </Stack>
                            <Stack direction="row" spacing={2} sx={{ mt: 0.5 }}>
                              <Typography
                                variant="caption"
                                sx={{
                                  opacity: isLikely ? 0.9 : 0.7,
                                }}
                              >
                                {result.percent_solvent.toFixed(1)}% solvent
                              </Typography>
                              <LinearProgress
                                variant="determinate"
                                value={probability * 100}
                                sx={{
                                  flex: 1,
                                  alignSelf: "center",
                                  height: 4,
                                  borderRadius: 2,
                                  bgcolor: isLikely
                                    ? "success.light"
                                    : "action.disabledBackground",
                                  "& .MuiLinearProgress-bar": {
                                    bgcolor: isLikely
                                      ? "success.contrastText"
                                      : "primary.main",
                                  },
                                }}
                              />
                            </Stack>
                          </Box>
                        );
                      }
                    )}
                  </Stack>
                </Box>
              ) : (
                <Typography variant="body2" color="text.secondary">
                  Provide MTZ file to calculate Matthews coefficient
                </Typography>
              )}
            </Grid2>
          </Grid2>
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
export default TaskInterface;
