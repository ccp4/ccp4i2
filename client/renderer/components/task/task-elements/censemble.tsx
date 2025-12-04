import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Card, CardContent, CardHeader, Grid2, Paper } from "@mui/material";
import { ErrorInfo } from "./error-info";
import { useMemo } from "react";

export const CEnsembleElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const { job, itemName } = props;
  const { useTaskItem, getValidationColor } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  const inferredVisibility = useMemo(() => {
    if (!props.visibility) return true;
    if (typeof props.visibility === "function") {
      return props.visibility();
    }
    return props.visibility;
  }, [props.visibility]);

  return inferredVisibility ? (
    <Paper
      sx={{
        border: "3px solid",
        borderColor: getValidationColor(item),
        pt: 2,
      }}
    >
      {item && (
        <Grid2 container rowSpacing={0}>
          <Grid2 key={"number"} size={{ xs: 4 }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
              itemName={`${item._objectPath}.number`}
              qualifiers={{
                ...useTaskItem(`${item._objectPath}.number`).item._qualifiers,
                guiLabel: "copies",
              }}
            />
          </Grid2>
          <Grid2 key={"label"} size={{ xs: 4 }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
              itemName={`${item._objectPath}.label`}
              qualifiers={{ ...props.qualifiers, guiLabel: "label" }}
            />
          </Grid2>
          <Grid2 key={"use"} size={{ xs: 4 }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
              itemName={`${item._objectPath}.use`}
              qualifiers={{
                ...useTaskItem(`${item._objectPath}.use`).item._qualifiers,
                guiLabel: "use",
              }}
            />
          </Grid2>
        </Grid2>
      )}
      <CCP4i2TaskElement
        {...props}
        sx={{ my: 0, py: 0, minWidth: "10rem" }}
        itemName={`${item._objectPath}.pdbItemList`}
        qualifiers={{
          ...useTaskItem(`${item._objectPath}.pdbItemList`).item._qualifiers,
          guiLabel: "PDBs",
        }}
      />
    </Paper>
  ) : null;
};
