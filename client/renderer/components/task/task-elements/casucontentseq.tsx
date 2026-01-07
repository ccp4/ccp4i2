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
import { ErrorInfo } from "./error-info";
import { apiGet } from "../../../api-fetch";
import { useCallback } from "react";
import { useParameterChangeIntent } from "../../../providers/parameter-change-intent-provider";
import { Science } from "@mui/icons-material";

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
  const { clearIntentForPath } = useParameterChangeIntent();

  const { item } = useTaskItem(itemName);
  const { update: setPolymerType } = useTaskItem(
    `${item._objectPath}.polymerType`
  );
  const { update: setName } = useTaskItem(`${item._objectPath}.name`);
  const { update: setSequence } = useTaskItem(`${item._objectPath}.sequence`);
  const { update: setDescription } = useTaskItem(
    `${item._objectPath}.description`
  );
  const { update: setNCopies } = useTaskItem(`${item._objectPath}.nCopies`);

  // Get current values for display
  const polymerType = item?._value?.polymerType?._value || "";
  const sequenceValue = item?._value?.sequence?._value || "";
  const seqLength = sequenceValue.replace(/\s/g, "").length;

  const setSEQUENCEFromSEQIN = useCallback(
    async (seqinDigestResponse: any, annotation: string) => {
      const seqinDigest = seqinDigestResponse?.data;
      console.log("Setting SEQUENCE from SEQIN digest", seqinDigest);
      if (
        !setSequence ||
        !setName ||
        !setPolymerType ||
        !setDescription ||
        !setNCopies ||
        !item ||
        !seqinDigest ||
        job?.status != 1
      ) {
        console.log(
          "Cannot set SEQUENCE from SEQIN - missing data",
          seqinDigest
        );
        return;
      }
      if (seqinDigest?.moleculeType) {
        console.log("Seqin digest was a sequence file");
        const { name, moleculeType, sequence } = seqinDigest || {};
        const sanitizedName = name.replace(/[^a-zA-Z0-9]/g, "_");
        await setPolymerType(moleculeType);
        await setName(sanitizedName);
        await setSequence(sequence);
        await setDescription(annotation);
        await setNCopies(1);
        clearIntentForPath(`${item._objectPath}.polymerType`);
        clearIntentForPath(`${item._objectPath}.name`);
        clearIntentForPath(`${item._objectPath}.sequence`);
        clearIntentForPath(`${item._objectPath}.description`);
        clearIntentForPath(`${item._objectPath}.nCopies`);
        await mutateContainer();
        props.onChange?.({ name, moleculeType, sequence });
      } else if (seqinDigest?.composition) {
        console.log("Seqin digest was a coordinate file");
        const { name, moleculeType, sequence } = {
          name: `Chain_${seqinDigest.composition.peptides[0]}`,
          moleculeType: "PROTEIN",
          sequence:
            seqinDigest.sequences[seqinDigest.composition.peptides[0]] || "",
        };
        const sanitizedName = name.replace(/[^a-zA-Z0-9]/g, "_");
        await setPolymerType(moleculeType);
        await setName(sanitizedName);
        await setSequence(sequence);
        await setDescription(annotation);
        await setNCopies(1);
        clearIntentForPath(`${item._objectPath}.polymerType`);
        clearIntentForPath(`${item._objectPath}.name`);
        clearIntentForPath(`${item._objectPath}.sequence`);
        clearIntentForPath(`${item._objectPath}.description`);
        clearIntentForPath(`${item._objectPath}.nCopies`);
        await mutateContainer();
        props.onChange?.({ name, moleculeType, sequence });
      }
    },
    [
      setSequence,
      setName,
      setPolymerType,
      setDescription,
      setNCopies,
      job,
      item,
      mutateContainer,
      clearIntentForPath,
      props,
    ]
  );

  const validationColor = getValidationColor(item);

  return (
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
          <ErrorInfo {...props} />
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
              Or import from file / fetch from UniProt or PDB
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName={`${item._objectPath}.source`}
              qualifiers={{
                guiLabel: "Source File",
                guiMode: "multiLine",
                mimeTypeName: "application/CCP4-seq",
                downloadModes: ["uniprotFasta", "ebiPdb"],
              }}
              onChange={async (updatedItem: any) => {
                console.log("Fetch file for param", updatedItem);
                const { dbFileId, annotation } = valueOfItem(updatedItem);
                const digest = await apiGet(`files/${dbFileId}/digest_by_uuid`);
                console.log({ digest, annotation });
                setSEQUENCEFromSEQIN(digest, annotation);
              }}
              suppressMutations={true}
            />
          </Box>
        </Stack>
      </Box>
    </Box>
  );
};
