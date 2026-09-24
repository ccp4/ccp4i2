import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Box, Typography } from "@mui/material";
import { CSimpleDataFileElement } from "./csimpledatafile";
import { useInferredVisibility } from "./hooks/useInferredVisibility";
import { FIELD_SIZES } from "./field-sizes";

/**
 * CPanddaEventElement — one PanDDA event in a receipt's EVENTS list: its
 * scores on one row, then its event map and candidate pose as files.
 *
 * Every pose is a candidate merged onto the apo model, never the model of
 * record, and the optimal contour is in absolute map units, not sigma
 * (design note §7.3). Event numbers are ordinals within one run only (§7.4).
 */
export const CPanddaEventElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName } = props;
  const { useTaskItem } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const isVisible = useInferredVisibility(props.visibility);
  if (!isVisible || !item) return null;

  const path = item._objectPath;
  const scalar = (name: string, label: string, width: string) => (
    <>
      <Typography variant="body2" sx={{ flexShrink: 0 }}>
        {label}
      </Typography>
      <Box sx={{ width, flexShrink: 0 }}>
        <CCP4i2TaskElement
          job={job}
          itemName={`${path}.${name}`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>
    </>
  );

  return (
    <Box
      sx={{
        borderLeft: "2px solid",
        borderColor: "divider",
        borderRadius: 0.5,
        pl: 1,
        py: 0.5,
        display: "flex",
        flexDirection: "column",
        gap: 0.5,
      }}
    >
      <Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
        {scalar("EVENT_IDX", "event", FIELD_SIZES.xs)}
        {scalar("LIGAND_ID", "ligand", FIELD_SIZES.xs)}
        {scalar("SITE_IDX", "site", FIELD_SIZES.xs)}
        {scalar("BDC", "BDC", FIELD_SIZES.sm)}
        {scalar("SCORE", "score", FIELD_SIZES.sm)}
        {scalar("BUILD_SCORE", "build", FIELD_SIZES.sm)}
        {scalar("RSCC", "RSCC", FIELD_SIZES.sm)}
        {scalar("HIT_PROBABILITY", "hit prob.", FIELD_SIZES.sm)}
        {scalar("OPTIMAL_CONTOUR", "contour (abs.)", FIELD_SIZES.sm)}
      </Box>
      <CSimpleDataFileElement
        {...props}
        itemName={`${path}.EVENT_MAP`}
        qualifiers={{ guiLabel: "Event map" }}
      />
      <CSimpleDataFileElement
        {...props}
        itemName={`${path}.POSE`}
        qualifiers={{ guiLabel: "Candidate pose" }}
      />
    </Box>
  );
};
