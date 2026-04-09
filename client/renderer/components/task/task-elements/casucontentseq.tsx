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
import { useCallback, useState } from "react";
import { Science } from "@mui/icons-material";
import { ChainPickerDialog } from "./chain-picker-dialog";
import type { ChainSequenceInfo } from "./mmcif-sequence-parser";

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

/** Classify chain type from composition lists */
function classifyChain(
  chainId: string,
  composition: any
): "PROTEIN" | "DNA" | "RNA" | "OTHER" {
  if (composition.peptides?.includes(chainId)) return "PROTEIN";
  if (composition.nucleics?.includes(chainId)) return "PROTEIN"; // nucleic → check chainDetails
  // Use chainDetails if available for finer classification
  const detail = composition.chainDetails?.find((d: any) => d.id === chainId);
  if (detail) {
    if (detail.type === "protein") return "PROTEIN";
    if (detail.type === "nucleic") return "DNA"; // Could be RNA, but chainDetails doesn't distinguish
  }
  return "OTHER";
}

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
  const { update: setNCopies } = useTaskItem(`${item._objectPath}.nCopies`);

  // Get current values for display
  const polymerType = item?._value?.polymerType?._value || "";
  const sequenceValue = item?._value?.sequence?._value || "";
  const seqLength = sequenceValue.replace(/\s/g, "").length;

  // Chain picker state for multi-chain coordinate file disambiguation
  const [chainPickerOpen, setChainPickerOpen] = useState(false);
  const [pendingChains, setPendingChains] = useState<ChainSequenceInfo[]>([]);
  const [pendingPdbId, setPendingPdbId] = useState("");
  const [pendingAnnotation, setPendingAnnotation] = useState("");

  /** Apply a selected chain's sequence to the ASU content fields */
  const applyChainSequence = useCallback(
    async (
      chainId: string,
      moleculeType: string,
      sequence: string,
      annotation: string
    ) => {
      if (
        !setSequence || !setName || !setPolymerType ||
        !setDescription || !setNCopies || !item || job?.status != 1
      ) return;

      const sanitizedName = `Chain_${chainId}`.replace(/[^a-zA-Z0-9_]/g, "_");
      await setPolymerType(moleculeType);
      await setName(sanitizedName);
      await setSequence(sequence);
      await setDescription(annotation);
      await setNCopies(1);
      props.onChange?.({ name: sanitizedName, moleculeType, sequence });
    },
    [setSequence, setName, setPolymerType, setDescription, setNCopies, job, item, props]
  );

  const setSEQUENCEFromSEQIN = useCallback(
    async (seqinDigestResponse: any, annotation: string) => {
      const seqinDigest = seqinDigestResponse?.data;
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
        return;
      }
      if (seqinDigest?.moleculeType) {
        const { name, moleculeType, sequence } = seqinDigest || {};
        const sanitizedName = name.replace(/[^a-zA-Z0-9]/g, "_");
        await setPolymerType(moleculeType);
        await setName(sanitizedName);
        await setSequence(sequence);
        await setDescription(annotation);
        await setNCopies(1);
        // With local patching, each update patches the cache immediately
        // No need to call mutateContainer or clear intents
        props.onChange?.({ name, moleculeType, sequence });
      } else if (seqinDigest?.composition && seqinDigest?.sequences) {
        const { composition, sequences } = seqinDigest;

        // Collect all polymer chains that have sequences
        const allPolymerChains = [
          ...(composition.peptides || []),
          ...(composition.nucleics || []),
        ];
        const chainsWithSeq = allPolymerChains.filter(
          (chainId: string) => sequences[chainId]
        );

        if (chainsWithSeq.length === 0) {
          return;
        }

        if (chainsWithSeq.length === 1) {
          // Single chain — apply directly
          const chainId = chainsWithSeq[0];
          const polyType = classifyChain(chainId, composition);
          await applyChainSequence(chainId, polyType, sequences[chainId], annotation);
        } else {
          // Multiple chains — show picker dialog
          const chainInfos: ChainSequenceInfo[] = chainsWithSeq.map(
            (chainId: string) => {
              const seq = sequences[chainId] || "";
              const polyType = classifyChain(chainId, composition);
              return {
                chainId,
                entityId: "",
                sequence: seq,
                polymerType: polyType,
                length: seq.length,
                description: `Chain ${chainId}`,
              };
            }
          );
          setPendingChains(chainInfos);
          setPendingPdbId(annotation || "structure");
          setPendingAnnotation(annotation);
          setChainPickerOpen(true);
        }
      } else if (seqinDigest?.composition) {
        // Legacy fallback: composition but no sequences dict (shouldn't happen with updated server)
        const chainId = seqinDigest.composition.peptides?.[0];
        if (chainId) {
          await applyChainSequence(chainId, "PROTEIN", "", annotation);
        }
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
      props,
      applyChainSequence,
    ]
  );

  /** Handle chain picker selection */
  const handleChainSelected = useCallback(
    (chain: ChainSequenceInfo) => {
      applyChainSequence(
        chain.chainId,
        chain.polymerType,
        chain.sequence,
        pendingAnnotation
      );
    },
    [applyChainSequence, pendingAnnotation]
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
                  downloadModes: ["uniprotFasta", "ebiPdb", "rcsbPdb"],
                }}
                onChange={async (updatedItem: any) => {
                  const { dbFileId, annotation } = valueOfItem(updatedItem);
                  const digest = await apiGet(`files/${dbFileId}/digest_by_uuid`);
                  setSEQUENCEFromSEQIN(digest, annotation);
                }}
                suppressMutations={true}
              />
            </Box>
          </Stack>
        </Box>
      </Box>

      <ChainPickerDialog
        open={chainPickerOpen}
        onClose={() => setChainPickerOpen(false)}
        onSelect={handleChainSelected}
        chains={pendingChains}
        pdbId={pendingPdbId}
      />
    </>
  );
};
