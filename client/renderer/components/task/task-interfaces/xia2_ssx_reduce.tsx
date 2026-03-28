/*
 * Copyright (C) 2026 Newcastle University
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
import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: workflow__steps } = useTaskItem("workflow__steps");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Search for xia2.ssx integrated files" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Search for" hint="files">
              <CCP4i2TaskElement
                itemName="SEARCH_PREFERENCE"
                {...props}
                qualifiers={{ guiLabel: " ", guiMode: "radio" }}
              />
            </InlineField>
            <CCP4i2TaskElement itemName="SEARCH_ROOT_DIR" {...props} />
            <CCP4i2TaskElement itemName="inputData.DIALS_INTEGRATED" {...props} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Basic parameters" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="d_min" {...props} />
            <CCP4i2TaskElement itemName="reference" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Advanced parameters">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Scaling options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="partiality_threshold"
              {...props}
            />
            {workflow__steps === "scale+merge" && (
              <CCP4i2TaskElement
                itemName="scaling__anomalous"
                {...props}
                qualifiers={{
                  guiLabel:
                    "Keep anomalous pairs separate during scaling",
                }}
              />
            )}
            <CCP4i2TaskElement
              itemName="dials_cosym_phil_d_min"
              {...props}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Workflow options" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Lattice tolerance for triggering mis-indexing assessment">
              <CCP4i2TaskElement
                itemName="symmetry__lattice_symmetry_max_delta"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Workflow:">
              <CCP4i2TaskElement
                itemName="workflow__steps"
                {...props}
                qualifiers={{ guiLabel: " ", guiMode: "radio" }}
              />
            </InlineField>
            <Typography variant="body2" color="text.secondary">
              Reduction of batch size and number of processors prevents system
              memory issues.
            </Typography>
            <CCP4i2TaskElement
              itemName="reduction_batch_size"
              {...props}
            />
            <CCP4i2TaskElement
              itemName="multiprocessing__nproc"
              {...props}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Unit cell filtering" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" color="text.secondary">
              A filter is applied based on the median cell and the absolute
              tolerances.
            </Typography>
            <InlineField label="Calculated median cell">
              <CCP4i2TaskElement
                itemName="MEDIAN_CELL"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Absolute angle tolerance (degrees)">
              <CCP4i2TaskElement
                itemName="clustering__absolute_angle_tolerance"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Absolute length tolerance (A)">
              <CCP4i2TaskElement
                itemName="clustering__absolute_length_tolerance"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Instead of using the median cell, use these central cell values for the cell filtering (Separate unit cell parameters using comma)">
              <CCP4i2TaskElement
                itemName="clustering__central_unit_cell"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Instead of filtering based on cell values and tolerances, use a clustering approach to select a cell cluster. Clustering threshold (Andrews-Bernstein distance)">
              <CCP4i2TaskElement
                itemName="clustering__threshold"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Space group" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Space group for scaling and merging">
              <CCP4i2TaskElement
                itemName="symmetry__space_group"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Reference model options (if generating reference from PDB model)",
            }}
            containerHint="FolderLevel"
          >
            <InlineField label="Average solvent density (k-sol)">
              <CCP4i2TaskElement
                itemName="reference_model__k_sol"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Average solvent B-factor (B-sol)">
              <CCP4i2TaskElement
                itemName="reference_model__b_sol"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Splitting mixed-condition data" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Dose series - number of repeated measurements at each point">
              <CCP4i2TaskElement
                itemName="dose_series_repeat"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <CCP4i2TaskElement itemName="grouping" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
