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
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: BLUR_MODE } = useTaskItem("BLUR_MODE");
  const { value: FORM_FACTOR } = useTaskItem("FORM_FACTOR");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Single folder: inputData */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
          Structure
        </Typography>
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ initiallyOpen: true }}
          containerHint="BlockLevel"
        >
          <CCP4i2TaskElement itemName="XYZIN" {...props} />
        </CCP4i2ContainerElement>

        <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
          Options
        </Typography>
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ initiallyOpen: true }}
          containerHint="BlockLevel"
        >
          <InlineField label="Scattering form factor">
            <CCP4i2TaskElement itemName="FORM_FACTOR" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          <InlineField label="Resolution limit" hint={"\u00C5"}>
            <CCP4i2TaskElement itemName="D_MIN" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
        </CCP4i2ContainerElement>

        <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
          Advanced Options
        </Typography>
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ initiallyOpen: true }}
          containerHint="BlockLevel"
        >
          <InlineField label="Oversampling rate">
            <CCP4i2TaskElement itemName="RATE" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          <InlineField label="Blurring mode">
            <CCP4i2TaskElement itemName="BLUR_MODE" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          {BLUR_MODE === "custom" && (
            <InlineField label="Blurring B-factor" hint={"\u00C5\u00B2"}>
              <CCP4i2TaskElement itemName="BLUR" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          )}
          {BLUR_MODE !== "none" && (
            <CCP4i2TaskElement itemName="UNBLUR" {...props} qualifiers={{ guiLabel: "Unblur when calculating the reciprocal-space map coefficients" }} />
          )}
          <InlineField label="Density cutoff">
            <CCP4i2TaskElement itemName="CUTOFF" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          {FORM_FACTOR === "xray" && (
            <CCP4i2TaskElement itemName="MOTT_BETHE" {...props} qualifiers={{ guiLabel: "Approximate electron scattering factors using the Mott-Bethe formula" }} />
          )}
        </CCP4i2ContainerElement>
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
