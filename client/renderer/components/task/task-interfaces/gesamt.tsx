import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const GesamtInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { value: PAIRMULTI } = useTaskItem("PAIRMULTI");

  const isPairwise = PAIRMULTI === "PAIR" || PAIRMULTI === undefined;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Mode selection */}
      <Typography variant="subtitle1" color="primary" fontWeight="bold">
        Pairwise (2 structures) or multiple alignment
      </Typography>
      <CCP4i2TaskElement
        itemName="PAIRMULTI"
        {...props}
        qualifiers={{ guiLabel: " " }}
      />

      {/* Pairwise mode */}
      {isPairwise && (
        <>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Structure to move" }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              itemName="XYZIN_QUERY"
              {...props}
              qualifiers={{ guiLabel: "Atomic model" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Fixed structure" }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              itemName="XYZIN_TARGET"
              {...props}
              qualifiers={{ guiLabel: "Atomic model" }}
            />
          </CCP4i2ContainerElement>
        </>
      )}

      {/* Multiple mode */}
      {!isPairwise && (
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Structures to superpose" }}
          containerHint="BlockLevel"
        >
          <CCP4i2TaskElement
            itemName="XYZIN_LIST"
            {...props}
            qualifiers={{ guiLabel: "Atomic model" }}
          />
        </CCP4i2ContainerElement>
      )}

      {/* Quality mode */}
      <CCP4i2TaskElement
        itemName="MODE"
        {...props}
        qualifiers={{ guiLabel: " " }}
      />
    </Paper>
  );
};

export default GesamtInterface;
