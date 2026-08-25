import React, { useCallback, useMemo, useState } from "react";
import { Box, Stack, Typography } from "@mui/material";

import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { CPathElement } from "./cpath";
import { useJob } from "../../../utils";
import { apiJson } from "../../../api-fetch";

interface SweepDigest {
  template: string;
  directory: string;
  start: number | null;
  end: number | null;
  count: number;
}

/**
 * One entry of a CXia2ImageSelectionList: an image, and optionally the range of
 * the sweep to process.
 *
 * xia2 is given a single image and discovers the rest of the sweep itself, so
 * the interface asks for one image "from each dataset" rather than for a list
 * of frames. Having picked one, we ask the server which sweep it belongs to and
 * fill in the range, both to confirm to the user that the sweep was found and
 * to give them concrete numbers to narrow.
 */
export const CXia2ImageSelectionElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName } = props;
  const { useTaskItem, setParameter } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const [sweep, setSweep] = useState<SweepDigest | null>(null);

  const basePath = useMemo(() => item?._objectPath ?? "", [item]);

  /**
   * A newly picked image supersedes any range left over from the previous one,
   * so adopt the whole sweep; the user narrows it from there.
   */
  const handleImagePicked = useCallback(
    async (updatedItem: any) => {
      const baseName = updatedItem?._value?.baseName?._value ?? "";
      const relPath = updatedItem?._value?.relPath?._value ?? "";
      const path = relPath ? `${relPath}/${baseName}` : baseName;
      if (!path) {
        setSweep(null);
        return;
      }
      try {
        const response = await apiJson<{ success: boolean; data: SweepDigest }>(
          `image_sweep/?path=${encodeURIComponent(path)}`
        );
        if (!response?.success) {
          setSweep(null);
          return;
        }
        setSweep(response.data);
        if (response.data.start !== null && response.data.end !== null) {
          await setParameter({
            object_path: `${basePath}.imageStart`,
            value: response.data.start,
          });
          await setParameter({
            object_path: `${basePath}.imageEnd`,
            value: response.data.end,
          });
        }
      } catch {
        // The sweep summary is a convenience; a failed lookup must not stop
        // the user from setting the image and typing a range by hand.
        setSweep(null);
      }
    },
    [basePath, setParameter]
  );

  if (!item) return null;

  return (
    <Box sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}>
      <CPathElement
        {...props}
        itemName={`${basePath}.imageFile`}
        mode="file"
        qualifiers={{
          guiLabel: "Image file",
          toolTip:
            "xia2 will automatically find the other images in this dataset",
        }}
        onChange={handleImagePicked}
      />

      {sweep && (
        <Typography variant="caption" color="text.secondary" sx={{ pl: 6 }}>
          {sweep.count === 1
            ? `${sweep.template} — single image`
            : `${sweep.template} — ${sweep.count} images, ${sweep.start} to ${sweep.end}`}
        </Typography>
      )}

      <Stack direction="row" alignItems="center" spacing={1} sx={{ pl: 6 }}>
        <Typography variant="body2">Image range</Typography>
        <Box sx={{ width: "8rem" }}>
          <CCP4i2TaskElement
            {...props}
            itemName={`${basePath}.imageStart`}
            qualifiers={{ guiLabel: " " }}
          />
        </Box>
        <Typography variant="body2">to</Typography>
        <Box sx={{ width: "8rem" }}>
          <CCP4i2TaskElement
            {...props}
            itemName={`${basePath}.imageEnd`}
            qualifiers={{ guiLabel: " " }}
          />
        </Box>
      </Stack>
    </Box>
  );
};

CXia2ImageSelectionElement.displayName = "CXia2ImageSelectionElement";
