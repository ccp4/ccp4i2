import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Grid2, Paper } from "@mui/material";
import { useMemo } from "react";
import { useInferredVisibility } from "./hooks/useInferredVisibility";

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

  if (!isVisible) return null;

  return (
    <Paper
      sx={{
        border: "3px solid",
        borderColor: getValidationColor(item),
        pt: 2,
      }}
    >
      {item && (
        <Grid2 container rowSpacing={0}>
          <Grid2 key="number" size={{ xs: 4 }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
              itemName={numberPath}
              qualifiers={numberQualifiers}
            />
          </Grid2>
          <Grid2 key="label" size={{ xs: 4 }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
              itemName={labelPath}
              qualifiers={labelQualifiers}
            />
          </Grid2>
          <Grid2 key="use" size={{ xs: 4 }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
              itemName={usePath}
              qualifiers={useQualifiers}
            />
          </Grid2>
        </Grid2>
      )}
      {item && (
        <CCP4i2TaskElement
          {...props}
          sx={{ my: 0, py: 0, minWidth: "10rem" }}
          itemName={pdbItemListPath}
          qualifiers={pdbItemListQualifiers}
        />
      )}
    </Paper>
  );
};
