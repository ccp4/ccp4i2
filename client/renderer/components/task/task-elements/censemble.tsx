/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Box, Stack } from "@mui/material";
import { useMemo } from "react";
import { useInferredVisibility } from "./hooks/useInferredVisibility";
import { FieldRow } from "./field-row";

export const CEnsembleElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const { job, itemName } = props;
  const { useTaskItem, getValidationColor } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  // Fetch all child items at top level to comply with React hooks rules
  const numberPath = item ? `${item._objectPath}.number` : "";
  const labelPath = item ? `${item._objectPath}.label` : "";
  const usePath = item ? `${item._objectPath}.use` : "";
  const pdbItemListPath = item ? `${item._objectPath}.pdbItemList` : "";

  const { item: numberItem } = useTaskItem(numberPath);
  const { item: useItem } = useTaskItem(usePath);
  const { item: pdbItemListItem } = useTaskItem(pdbItemListPath);

  const isVisible = useInferredVisibility(props.visibility);

  // Memoize qualifiers to avoid creating new objects on each render
  const numberQualifiers = useMemo(
    () => ({
      ...numberItem?._qualifiers,
      guiLabel: "copies",
    }),
    [numberItem?._qualifiers]
  );

  const labelQualifiers = useMemo(
    () => ({
      ...props.qualifiers,
      guiLabel: "label",
    }),
    [props.qualifiers]
  );

  const useQualifiers = useMemo(
    () => ({
      ...useItem?._qualifiers,
      guiLabel: "use",
    }),
    [useItem?._qualifiers]
  );

  const pdbItemListQualifiers = useMemo(
    () => ({
      ...pdbItemListItem?._qualifiers,
      guiLabel: "PDBs",
    }),
    [pdbItemListItem?._qualifiers]
  );

  const validationColor = getValidationColor(item);

  if (!isVisible) return null;

  return (
    <Box
      sx={{
        borderLeft: "3px solid",
        borderColor: validationColor,
        borderRadius: 1,
        bgcolor: "action.hover",
        py: 1,
        pl: 1,
      }}
    >
      {item && (
        <Stack spacing={0.5}>
          {/* Ensemble metadata: copies, label, use — in a compact row */}
          <FieldRow equalWidth={false} size="sm" sx={{ alignItems: "center" }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0 }}
              itemName={numberPath}
              qualifiers={numberQualifiers}
            />
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0 }}
              itemName={labelPath}
              qualifiers={labelQualifiers}
            />
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0 }}
              itemName={usePath}
              qualifiers={useQualifiers}
            />
          </FieldRow>

          {/* PDB list */}
          <CCP4i2TaskElement
            {...props}
            sx={{ my: 0, py: 0 }}
            itemName={pdbItemListPath}
            qualifiers={pdbItemListQualifiers}
          />
        </Stack>
      )}
    </Box>
  );
};
