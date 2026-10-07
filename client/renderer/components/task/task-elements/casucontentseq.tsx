import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import {
  Box,
  Chip,
  Divider,
  Grid2,
  Stack,
  Typography,
} from "@mui/material";
import { useJob, valueOfItem } from "../../../utils";
import { apiGet } from "../../../api-fetch";
import { useCallback } from "react";
import { Science } from "@mui/icons-material";
import { polymerTypeFromMolecule } from "./mmcif-sequence-parser";

/** Get color for polymer type */
const getPolymerTypeColor = (
  type: string
): "primary" | "secondary" | "success" | "warning" | "info" => {
  switch (type?.toUpperCase()) {
    case "PROTEIN":
      return "primary";
    case "DNA":
      return "info";
    case "RNA":
      return "success";
    default:
      return "warning";
  }
};

export const CAsuContentSeqElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { itemName, job } = props;
  const { useTaskItem, getValidationColor, mutateContainer } = useJob(job.id);

  const { item } = useTaskItem(itemName);
  const { update: setPolymerType } = useTaskItem(
    `${item._objectPath}.polymerType`
  );
  const { update: setName } = useTaskItem(`${item._objectPath}.name`);
  const { update: setSequence } = useTaskItem(`${item._objectPath}.sequence`);
  const { update: setDescription } = useTaskItem(
    `${item._objectPath}.description`
  );

  // Get current values for display
  const polymerType = item?._value?.polymerType?._value || "";
  const sequenceValue = item?._value?.sequence?._value || "";
  const seqLength = sequenceValue.replace(/\s/g, "").length;

  /**
   * A sequence file was chosen for THIS row: replace the row's sequence with
   * it. The slot is a sequence-file slot, and a PDB entry fetched into it
   * arrives already as a FASTA of one chain, so the digest is always
   * { name, moleculeType, sequence }. Whole models and entries belong to the
   * page's model slot, which adds every chain.
   */
  const setSEQUENCEFromSEQIN = useCallback(
    async (seqinDigestResponse: any, annotation: string) => {
      const seqinDigest = seqinDigestResponse?.data;
      if (
        !setSequence ||
        !setName ||
        !setPolymerType ||
        !setDescription ||
        !item ||
        !seqinDigest?.moleculeType ||
        job?.status != 1
      ) {
        return;
      }
      const { name, moleculeType, sequence } = seqinDigest;
      const sanitizedName = name.replace(/[^a-zA-Z0-9]/g, "_");
      // CSequence says PROTEIN/NUCLEIC; the row wants PROTEIN/DNA/RNA.
      await setPolymerType(polymerTypeFromMolecule(moleculeType, sequence));
      await setName(sanitizedName);
      await setSequence(sequence);
      await setDescription(annotation);
      // The number of copies is the user's to say, not the file's: keep it.
      props.onChange?.({ name, moleculeType, sequence });
    },
    [setSequence, setName, setPolymerType, setDescription, job, item, props]
  );

  const validationColor = getValidationColor(item);

  return (
    <>
      <Box
        sx={{
          borderRadius: 2,
          border: 1,
          borderColor: validationColor !== "inherit" ? validationColor : "divider",
          borderLeftWidth: validationColor !== "inherit" ? 4 : 1,
          bgcolor: "background.paper",
          overflow: "hidden",
        }}
      >
        {/* Header */}
        <Box
          sx={{
            px: 2,
            py: 1.5,
            bgcolor: "action.hover",
            borderBottom: 1,
            borderColor: "divider",
          }}
        >
          <Stack
            direction="row"
            justifyContent="space-between"
            alignItems="center"
          >
            <Stack direction="row" alignItems="center" spacing={1.5}>
              <Science color="primary" />
              <Typography variant="subtitle1" fontWeight="bold">
                {item._qualifiers.guiLabel || "Sequence Details"}
              </Typography>
              {polymerType && (
                <Chip
                  label={polymerType}
                  size="small"
                  color={getPolymerTypeColor(polymerType)}
                  variant="outlined"
                />
              )}
              {seqLength > 0 && (
                <Typography variant="body2" color="text.secondary">
                  {seqLength} residues
                </Typography>
              )}
            </Stack>
          </Stack>
        </Box>

        {/* Content */}
        <Box sx={{ p: 2 }}>
          <Stack spacing={2}>
            {/* Basic info row */}
            <Grid2 container spacing={2}>
              {item && (
                <Grid2 size={{ xs: 12, sm: 4 }}>
                  <CCP4i2TaskElement
                    {...props}
                    itemName={`${item._objectPath}.nCopies`}
                    qualifiers={{
                      guiLabel: "Copies in ASU",
                      enumerators: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                    }}
                  />
                </Grid2>
              )}
              {item && (
                <Grid2 size={{ xs: 12, sm: 4 }}>
                  <CCP4i2TaskElement
                    {...props}
                    itemName={`${item._objectPath}.polymerType`}
                    qualifiers={{
                      guiLabel: "Polymer Type",
                    }}
                  />
                </Grid2>
              )}
              {item && (
                <Grid2 size={{ xs: 12, sm: 4 }}>
                  <CCP4i2TaskElement
                    {...props}
                    itemName={`${item._objectPath}.name`}
                    qualifiers={{
                      guiLabel: "Name",
                    }}
                  />
                </Grid2>
              )}
            </Grid2>

            {/* Description */}
            <CCP4i2TaskElement
              {...props}
              itemName={`${item._objectPath}.description`}
              qualifiers={{
                guiLabel: "Description",
                guiMode: "multiLine",
              }}
            />

            {/* Sequence */}
            <Box
              sx={{
                "& textarea": {
                  fontFamily: "monospace !important",
                  fontSize: "0.85rem !important",
                  letterSpacing: "0.05em",
                  lineHeight: 1.6,
                },
              }}
            >
              <CCP4i2TaskElement
                {...props}
                itemName={`${item._objectPath}.sequence`}
                qualifiers={{
                  guiLabel: "Sequence",
                  guiMode: "multiLine",
                }}
              />
            </Box>

            <Divider />

            {/* Source file */}
            <Box>
              <Typography variant="caption" color="text.secondary" sx={{ mb: 1, display: "block" }}>
                Replace this sequence from a file, or from UniProt (a residue range trims it to the crystallised construct). Whole models and PDB entries are added from the page, not here.
              </Typography>
              <CCP4i2TaskElement
                {...props}
                itemName={`${item._objectPath}.source`}
                qualifiers={{
                  guiLabel: "Sequence file or UniProt",
                  guiMode: "multiLine",
                  mimeTypeName: "application/CCP4-seq",
                  downloadModes: ["uniprotFasta"],
                }}
                onChange={async (updatedItem: any) => {
                  const { dbFileId, annotation } = valueOfItem(updatedItem);
                  const digest = await apiGet(`files_by_uuid/${dbFileId}/digest/`);
                  setSEQUENCEFromSEQIN(digest, annotation);
                }}
                suppressMutations={true}
              />
            </Box>
          </Stack>
        </Box>
      </Box>

    </>
  );
};
