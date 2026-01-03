import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: workflow__steps } = useTaskItem("workflow__steps");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Search for xia2.ssx integrated files
          </Typography>
          <CCP4i2TaskElement itemName="SEARCH_PREFERENCE" {...props} qualifiers={{ guiLabel: "files" }} />
          <CCP4i2TaskElement itemName="SEARCH_ROOT_DIR" {...props} />
          <CCP4i2TaskElement itemName="DIALS_INTEGRATED" {...props} />
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Basic parameters
            </Typography>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced parameters">
          {(workflow__steps === "scale+merge") && (
            <CCP4i2TaskElement itemName="scaling__anomalous" {...props} qualifiers={{ guiLabel: "Keep anomalous pairs separate during scaling" }} />
          )}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Workflow options
            </Typography>
            <CCP4i2TaskElement itemName="symmetry__lattice_symmetry_max_delta" {...props} qualifiers={{ guiLabel: "Lattice tolerance for triggering mis-indexing assessment" }} />
            <CCP4i2TaskElement itemName="-guiMode" {...props} qualifiers={{ guiLabel: "Workflow:" }} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Reduction of batch size and number of processors prevents system memory issues.
            </Typography>
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Unit cell filtering
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              A filter is applied based on the median cell and the absolute tolerances.
            </Typography>
            <CCP4i2TaskElement itemName="MEDIAN_CELL" {...props} qualifiers={{ guiLabel: "Calculated median cell" }} />
            <CCP4i2TaskElement itemName="clustering__absolute_angle_tolerance" {...props} qualifiers={{ guiLabel: "Absolute angle tolerance (degrees)" }} />
            <CCP4i2TaskElement itemName="clustering__absolute_length_tolerance" {...props} qualifiers={{ guiLabel: "Absolute length tolerance (A)" }} />
            <CCP4i2TaskElement itemName="clustering__central_unit_cell" {...props} qualifiers={{ guiLabel: "Instead of using the median cell, use these central cell values for the cell filtering (Separate unit cell parameters using comma)" }} />
            <CCP4i2TaskElement itemName="clustering__threshold" {...props} qualifiers={{ guiLabel: "Instead of filtering based on cell values and tolerances, use a clustering approach to select a cell cluster. Clustering threshold (Andrews–Bernstein distance)" }} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="symmetry__space_group" {...props} qualifiers={{ guiLabel: "Space group for scaling and merging" }} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Reference model options (if generating reference from PDB model)
            </Typography>
            <CCP4i2TaskElement itemName="reference_model__k_sol" {...props} qualifiers={{ guiLabel: "Average solvent density (k-sol)" }} />
            <CCP4i2TaskElement itemName="reference_model__b_sol" {...props} qualifiers={{ guiLabel: "Average solvent B-factor (B-sol)" }} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Splitting mixed-condition data
            </Typography>
            <CCP4i2TaskElement itemName="dose_series_repeat" {...props} qualifiers={{ guiLabel: "Dose series - number of repeated measurements at each point" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;