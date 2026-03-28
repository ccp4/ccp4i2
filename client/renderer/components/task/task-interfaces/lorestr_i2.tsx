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
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Main inputs
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Optional additional inputs
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="TLSIN" {...props} />
            <CCP4i2TaskElement itemName="DICT" {...props} />
            <CCP4i2TaskElement itemName="REFERENCE_LIST" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Options">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Main options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <InlineField label="Automatically fetch homologues from following databases:">
              <CCP4i2TaskElement itemName="AUTO" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <CCP4i2TaskElement itemName="OVB" {...props} qualifiers={{ guiLabel: "Refine overall B-factors (for very low resolution)" }} />
            <CCP4i2TaskElement itemName="DNA" {...props} qualifiers={{ guiLabel: "Generate external restraints for DNA/RNA chains" }} />
            <InlineField label="Number of CPUs to use:">
              <CCP4i2TaskElement itemName="CPU" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            After Molecular Replacement
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MR" {...props} qualifiers={{ guiLabel: "Run 100-200 cycles of jelly body refinement first (if model is straight after MR)" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="advancedOptions" label="Advanced Options">
          <CCP4i2TaskElement itemName="SS" {...props} qualifiers={{ guiLabel: "Save disk space by removing excessive ProSMART output" }} />
          <InlineField label="Download and use homologues with resolution better than" hint="Angstrom">
            <CCP4i2TaskElement itemName="MINRES" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          <InlineField label="Use up to" hint="homologues for restraints">
            <CCP4i2TaskElement itemName="NH" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          <InlineField label="Use up to" hint="chains to generate restraints">
            <CCP4i2TaskElement itemName="NC" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
