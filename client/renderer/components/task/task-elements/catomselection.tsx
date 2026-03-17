import React, { useMemo } from "react";
import { Chip, Stack, Typography } from "@mui/material";

import { CCP4i2TaskElementProps } from "./task-element";
import { CSimpleElement } from "./csimple";
import { useJob } from "../../../utils";
import { useApi } from "../../../api";
import { useInferredVisibility } from "./hooks/useInferredVisibility";

/**
 * CAtomSelectionElement renders the `text` sub-item of a CAtomSelection
 * as an editable text field, with a helper label showing the parent PDB
 * file and a live atom count resolved via the object_method API.
 */
export const CAtomSelectionElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const { job, itemName, visibility, qualifiers, ...restProps } = props;
  const { useTaskItem } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const api = useApi();

  // Look up the parent PDB file referenced by the pdbFileKey qualifier
  const pdbFileKey = item?._qualifiers?.pdbFileKey || qualifiers?.pdbFileKey;
  const { item: pdbItem } = useTaskItem(pdbFileKey || "");

  const isVisible = useInferredVisibility(visibility);

  // Extract values needed for display and the atom count query
  const selectionText = item?._value?.text?._value || "";
  const dbFileId = pdbItem?._value?.dbFileId?._value?.trim() || null;
  // object_method expects "plugin.inputData.ITEM" but _objectPath includes
  // "container." (e.g. "phaser_rnp_pipeline.container.inputData.XYZIN_PARENT").
  // Strip the "container." segment to match the expected format.
  const pdbObjectPath = (pdbItem?._objectPath || "").replace(".container.", ".");

  // Build a short description of the parent PDB
  const parentLabel = useMemo(() => {
    if (!pdbFileKey) return null;
    const annotation = pdbItem?._value?.annotation?._value;
    const baseName = pdbItem?._value?.baseName?._value;
    if (dbFileId && (annotation || baseName)) {
      return annotation || baseName;
    }
    return null;
  }, [pdbFileKey, pdbItem, dbFileId]);

  // Call atomCount on the parent CPdbDataFile via object_method endpoint.
  // Re-fetches when the selection string or the parent PDB file changes.
  const hasValidQuery = Boolean(
    pdbObjectPath && dbFileId && selectionText.trim()
  );
  const { data: atomCountResult } = api.objectMethod<any>(
    job.id,
    pdbObjectPath,
    "atomCount",
    {},
    [selectionText, dbFileId],
    hasValidQuery,
    [selectionText],
  );

  const atomCount: number | null = useMemo(() => {
    const n = atomCountResult?.data?.result;
    if (typeof n === "number" && n >= 0) return n;
    return null;
  }, [atomCountResult]);

  if (!isVisible || !item) return null;

  // The text child lives at item._value.text
  const textChild = item._value?.text;
  if (!textChild?._objectPath) return null;

  return (
    <Stack direction="row" spacing={1} alignItems="center" sx={{ width: "100%" }}>
      <CSimpleElement
        {...restProps}
        job={job}
        itemName={textChild._objectPath}
        qualifiers={qualifiers}
        visibility={visibility}
        type="text"
      />
      {atomCount !== null && (
        <Chip
          label={`${atomCount} atoms`}
          size="small"
          variant="outlined"
          color="default"
        />
      )}
      {pdbFileKey && (
        <Typography
          variant="caption"
          color="text.secondary"
          sx={{ whiteSpace: "nowrap", flexShrink: 0 }}
        >
          {parentLabel
            ? `on ${parentLabel}`
            : "Select parent coordinates above"}
        </Typography>
      )}
    </Stack>
  );
};
