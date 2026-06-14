import React from "react";
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { InlineField } from "../task-elements/inline-field";

/**
 * dm_multidomain - multi-domain NCS averaging via dm.
 *
 * The user supplies a model and per-domain residue ranges (the DOMAINS list of
 * CDmDomain items); the server derives the per-domain NCS operators and masks.
 * Different domains can follow different NCS transformations and be averaged,
 * refined, or excluded independently.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Input Data ===== */}
        <CCP4i2Tab label="Input Data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflection Data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props}
              qualifiers={{ guiLabel: "Reflections (F/SIGF or intensities)" }} />
            <CCP4i2TaskElement itemName="ABCD" {...props}
              qualifiers={{ guiLabel: "Starting phases" }} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props}
              qualifiers={{ guiLabel: "Free R set (optional)" }} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "NCS model" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props}
              qualifiers={{ guiLabel: "Model (one chain per NCS copy)" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Domains ===== */}
        <CCP4i2Tab label="Domains">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "NCS copies" }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" sx={{ mb: 1 }}>
              Leave blank to auto-detect the reference copy (first protein chain)
              and its NCS copies (the rest).
            </Typography>
            <InlineField label="Reference chain">
              <CCP4i2TaskElement itemName="REFERENCE_CHAIN" {...props}
                qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Copy chains (comma-separated)">
              <CCP4i2TaskElement itemName="COPY_CHAINS" {...props}
                qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Domains (residue range + averaging mode)" }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" sx={{ mb: 1 }}>
              One row per domain. Each domain is superposed independently across
              the NCS copies, so different domains may follow different NCS
              transformations. Mode: <b>average</b> (rigid NCS), <b>refine</b>
              (let dm refine the operators), or <b>exclude</b>.
            </Typography>
            <CCP4i2TaskElement itemName="DOMAINS" {...props}
              qualifiers={{ guiLabel: " " }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Parameters ===== */}
        <CCP4i2Tab label="Parameters">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Density modification" }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MODE_SOLVENT" {...props}
              qualifiers={{ guiLabel: "Solvent flattening" }} />
            <CCP4i2TaskElement itemName="MODE_HISTOGRAM" {...props}
              qualifiers={{ guiLabel: "Histogram matching" }} />
            <InlineField label="Number of cycles">
              <CCP4i2TaskElement itemName="NCYCLES" {...props}
                qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Solvent content"
              hint="blank = estimate from model and cell">
              <CCP4i2TaskElement itemName="SOLVENT_CONTENT" {...props}
                qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Mask radius" hint="Angstrom">
              <CCP4i2TaskElement itemName="MASK_RADIUS" {...props}
                qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
