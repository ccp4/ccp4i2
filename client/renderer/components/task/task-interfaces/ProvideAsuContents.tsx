import {
  Alert,
  Grid2,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Typography,
} from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback } from "react";
import useSWR from "swr";
import { apiPost } from "../../../api-fetch";
import { BaseSpacegroupCellElement } from "../task-elements/base-spacegroup-cell-element";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, useFileDigest, fetchDigest, getErrors, mutateValidation } = useJob(job.id);
  const { update: setAsuContent } = useTaskItem("ASU_CONTENT");
  const { item: asuContentInItem } = useTaskItem("ASUCONTENTIN");
  const { item: asuContentItem } = useTaskItem("ASU_CONTENT");

  // File digest for HKLIN (used for Matthews calculation)
  // Uses SWR for automatic caching and revalidation - works on initial load
  // and when returning to a previously configured interface
  const { data: HKLINDigest, mutate: mutateHKLINDigest } = useFileDigest(
    `ProvideAsuContents.inputData.HKLIN`
  );

  // ASU content is valid when there are no validation errors for it
  // Uses the existing validation infrastructure from useJob().getErrors()
  // which fetches from /api/jobs/{id}/validation/ endpoint
  const asuContentErrors = getErrors(asuContentItem);
  const isAsuContentValid = asuContentErrors.length === 0 && (asuContentItem?._value?.length ?? 0) > 0;

  /**
   * Fetches the molecular weight for the current job's ASU content using SWR.
   * Only fetches when ASU_CONTENT is valid (passes all validation checks).
   */
  const { data: molWeight, mutate: mutateMolWeight } = useSWR(
    isAsuContentValid ? [`jobs/${job.id}/object_method`, "molecularWeight"] : null,
    ([url]) =>
      apiPost(url, {
        object_path: "ProvideAsuContents.inputData.ASU_CONTENT",
        method_name: "molecularWeight",
      })
  );

  /**
   * Fetches and caches the Matthews coefficient analysis for the current job using SWR.
   * Only fetches when we have both a valid molecular weight result AND an HKLIN file.
   */
  const { data: matthewsAnalysis, mutate: mutateMatthews } = useSWR(
    // Only fetch when we have molecular weight AND HKLIN digest
    // API response format: {success: true, data: {result: <value>}}
    molWeight?.data?.result && HKLINDigest
      ? [
          `jobs/${job.id}/object_method`,
          "matthewsCoeff",
          molWeight.data.result,
          HKLINDigest,
        ]
      : null,
    ([url, , molWeightResult]) =>
      apiPost(url, {
        object_path: "ProvideAsuContents.inputData.HKLIN.fileContent",
        method_name: "matthewsCoeff",
        kwargs: { molWt: molWeightResult },
      }),
    { keepPreviousData: true }
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
      await setAsuContent(seqList);
      // Refresh validation, molecular weight, and Matthews after updating ASU content
      // Must await molWeight mutation so the new value is available for Matthews calculation
      mutateValidation();
      await mutateMolWeight();
      mutateMatthews();
    }
  }, [asuContentInItem?._objectPath, fetchDigest, setAsuContent, mutateValidation, mutateMolWeight, mutateMatthews]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{
            guiLabel: "Optionally load existing AU content file to edit",
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
                <Table>
                  <TableHead>
                    <TableRow>
                      <TableCell>Multiplier</TableCell>
                      <TableCell>%Solvent</TableCell>
                      <TableCell>Probability</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {matthewsAnalysis?.data?.result?.results.map((result) => (
                      <TableRow key={result.nmol_in_asu}>
                        <TableCell>{result.nmol_in_asu}</TableCell>
                        <TableCell>
                          {result.percent_solvent.toFixed(2)}
                        </TableCell>
                        <TableCell>{result.prob_matth.toFixed(2)}</TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              ) : (
                "Provide MTZ file to calculate Matthews coefficient"
              )}
            </Grid2>
          </Grid2>
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
export default TaskInterface;
