import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: AMPLE_EXISTING_MODELS } = useTaskItem("AMPLE_EXISTING_MODELS");
  const { value: AMPLE_MODELS_SOURCE } = useTaskItem("AMPLE_MODELS_SOURCE");
  const { value: AMPLE_MODEL_TYPE } = useTaskItem("AMPLE_MODEL_TYPE");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Input sequence
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="AMPLE_SEQIN" {...props} qualifiers={{ toolTip: "Input sequence" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Input reflections
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="AMPLE_F_SIGF" {...props} qualifiers={{ toolTip: "Input reflections" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Class of protein:
          </Typography>
          <CCP4i2TaskElement itemName="AMPLE_PROTEIN_CLASS" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Do you have existing models?
          </Typography>
          <CCP4i2TaskElement itemName="AMPLE_EXISTING_MODELS" {...props} />
          {(AMPLE_EXISTING_MODELS === "True") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
                Models are from:
              </Typography>
              <CCP4i2TaskElement itemName="AMPLE_MODELS_SOURCE" {...props} />
              {(AMPLE_MODELS_SOURCE === "directory") && (
                <CCP4i2TaskElement itemName="AMPLE_MODELS_DIR" {...props} />
              )}
              {(AMPLE_MODELS_SOURCE === "file") && (
                <CCP4i2TaskElement itemName="AMPLE_MODELS_FILE" {...props} />
              )}
              <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
                What sort of models are these?
              </Typography>
              <CCP4i2TaskElement itemName="AMPLE_MODEL_TYPE" {...props} />
              <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
                NMR remodelling:
              </Typography>
              <CCP4i2TaskElement itemName="AMPLE_NMR_REMODEL" {...props} />
            </CCP4i2ContainerElement>
          )}
          {(AMPLE_EXISTING_MODELS === "True") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
                Model source:
              </Typography>
              <CCP4i2TaskElement itemName="AMPLE_MODEL_GENERATION" {...props} />
            </CCP4i2ContainerElement>
          )}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              ROSETTA paths
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              ROSETTA can be downloaded from: https://www.rosettacommons.org
            </Typography>
            <CCP4i2TaskElement itemName="AMPLE_ROSETTA_DIR" {...props} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Fragments should be created on the ROBETTA SERVER: http://robetta.bakerlab.org
            </Typography>
            <CCP4i2TaskElement itemName="AMPLE_ROSETTA_FRAGS3" {...props} />
            <CCP4i2TaskElement itemName="AMPLE_ROSETTA_FRAGS9" {...props} />
            <CCP4i2TaskElement itemName="AMPLE_CONTACT_FILE" {...props} />
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              What is the format of the contact file?
            </Typography>
            <CCP4i2TaskElement itemName="AMPLE_CONTACT_FORMAT" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              x
            </Typography>
            <CCP4i2TaskElement itemName="AMPLE_NPROC" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="inputData" label="Advanced Options">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Rebuilding options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="AMPLE_REFINE_REBUILD" {...props} qualifiers={{ guiLabel: "x" }} />
            <CCP4i2TaskElement itemName="AMPLE_USE_SHELXE" {...props} qualifiers={{ guiLabel: "x" }} />
            <CCP4i2TaskElement itemName="AMPLE_SHELXE_REBUILD" {...props} qualifiers={{ guiLabel: "x" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Ensembling options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="AMPLE_ENSEMBLING_TM" {...props} qualifiers={{ guiLabel: "x" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            x
          </Typography>
          <CCP4i2TaskElement itemName="AMPLE_EXTRA_FLAGS" {...props} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;