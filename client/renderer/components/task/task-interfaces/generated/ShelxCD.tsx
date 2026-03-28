/*
 * Copyright (C) 2025-2026 Newcastle University
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
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MODE } = useTaskItem("MODE");

  if (!container) return <LinearProgress />;

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <InlineField label="Experiment type">
            <CCP4i2TaskElement itemName="MODE" {...props} qualifiers={{ guiLabel: " ", toolTip: "What sort of data are available for finding heavy atoms" }} />
          </InlineField>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <InlineField label="Atom type to seek">
              <CCP4i2TaskElement itemName="SFAC" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Number to find">
              <CCP4i2TaskElement itemName="FIND" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Number of trys">
              <CCP4i2TaskElement itemName="NTRY" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            {MODE === "SAD" && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Anomalous (SAD) dataset
                </Typography>
                <CCP4i2TaskElement itemName="SAD" {...props} />
              </>
            )}

            {MODE === "MAD" && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  High energy remote dataset
                </Typography>
                <CCP4i2TaskElement itemName="HREM" {...props} />
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Low energy remote dataset
                </Typography>
                <CCP4i2TaskElement itemName="LREM" {...props} />
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Inflection point dataset
                </Typography>
                <CCP4i2TaskElement itemName="INFL" {...props} />
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Peak dataset
                </Typography>
                <CCP4i2TaskElement itemName="PEAK" {...props} />
              </>
            )}

            {MODE === "SIR" && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Derivative dataset
                </Typography>
                <CCP4i2TaskElement itemName="SIR" {...props} />
              </>
            )}

            {MODE === "SIRAS" && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Anomalous derivative dataset
                </Typography>
                <CCP4i2TaskElement itemName="SIRA" {...props} />
              </>
            )}

            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Native dataset
            </Typography>
            <CCP4i2TaskElement itemName="NAT" {...props} />

            {MODE === "RIP" && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Radiation damaged dataset
                </Typography>
                <CCP4i2TaskElement itemName="RIP" {...props} />
              </>
            )}

            {MODE === "RIPAS" && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Radiation damaged anomalous dataset
                </Typography>
                <CCP4i2TaskElement itemName="RIPA" {...props} />
              </>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Keywords">
          {/* Auto-generated keywords excluding MODE, FIND, NTRY, SFAC */}
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
