import React, { useMemo } from "react";
import { Box, Typography } from "@mui/material";
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";

/**
 * CRunBatchRangeElement — renders a single batch range item inline:
 *   Run number [___]  for batch range [___]  to [___]  resolution limit [___]
 */
export const CRunBatchRangeElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName } = props;
  const { useTaskItem } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  const basePath = useMemo(() => item?._objectPath ?? "", [item]);

  if (!item) return null;

  return (
    <Box
      sx={{
        display: "flex",
        alignItems: "center",
        gap: 1,
        flexWrap: "wrap",
      }}
    >
      <Typography variant="body2">Run number</Typography>
      <Box sx={{ width: "6rem" }}>
        <CCP4i2TaskElement
          {...props}
          itemName={`${basePath}.runNumber`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>

      <Typography variant="body2">for batch range</Typography>
      <Box sx={{ width: "6rem" }}>
        <CCP4i2TaskElement
          {...props}
          itemName={`${basePath}.batchRange0`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>

      <Typography variant="body2">to</Typography>
      <Box sx={{ width: "6rem" }}>
        <CCP4i2TaskElement
          {...props}
          itemName={`${basePath}.batchRange1`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>

      <Typography variant="body2">resolution limit</Typography>
      <Box sx={{ width: "6rem" }}>
        <CCP4i2TaskElement
          {...props}
          itemName={`${basePath}.resolution`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>
    </Box>
  );
};
