import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Box, Typography } from "@mui/material";
import { CSimpleDataFileElement } from "./csimpledatafile";
import { useInferredVisibility } from "./hooks/useInferredVisibility";
import { FIELD_SIZES } from "./field-sizes";

/**
 * CPanddaDatasetElement — one crystal in a PanDDA run's DATASETS list:
 *
 *   label [ BAZ2BA-x425 ]
 *   Refined model     [ file ]
 *   Reflections       [ file ]
 *   Ligand dictionary [ file ]   (optional)
 *
 * The list is the whole input of the orchestrator (design note §3.1): it
 * never asks which datasets a campaign has, so what is shown here is exactly
 * what will be staged, and a dataset the campaign rule missed is added by
 * hand. The label is what Projects.csv records against the clean xtal-NNNN
 * name PanDDA itself sees.
 */
export const CPanddaDatasetElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName } = props;
  const { useTaskItem, getValidationColor } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const isVisible = useInferredVisibility(props.visibility);
  if (!isVisible || !item) return null;

  const colour = getValidationColor(item);
  const path = item._objectPath;

  return (
    <Box
      sx={{
        borderLeft: "2px solid",
        borderColor: colour,
        borderRadius: 0.5,
        pl: 1,
        py: 0.5,
        display: "flex",
        flexDirection: "column",
        gap: 0.5,
      }}
    >
      <Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
        <Typography variant="body2" sx={{ flexShrink: 0 }}>
          label
        </Typography>
        <Box sx={{ width: FIELD_SIZES.lg, flexShrink: 0 }}>
          <CCP4i2TaskElement
            job={job}
            itemName={`${path}.DTAG`}
            qualifiers={{ guiLabel: " " }}
          />
        </Box>
      </Box>
      <CSimpleDataFileElement
        {...props}
        itemName={`${path}.XYZIN`}
        qualifiers={{ guiLabel: "Refined model" }}
      />
      <CSimpleDataFileElement
        {...props}
        itemName={`${path}.HKLIN`}
        qualifiers={{ guiLabel: "Reflections" }}
      />
      <CSimpleDataFileElement
        {...props}
        itemName={`${path}.DICT`}
        qualifiers={{ guiLabel: "Ligand dictionary" }}
      />
    </Box>
  );
};
