import React from "react";
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";

/**
 * dm_multidomain - multi-domain NCS averaging via dm.
 *
 * The user supplies a model and per-domain residue ranges (the DOMAINS list of
 * CDmDomain items); the server derives the per-domain NCS operators and masks.
 * Different domains can follow different NCS transformations and be averaged,
 * refined, or excluded independently.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { value: PHASE_SOURCE } = useTaskItem("PHASE_SOURCE");
  const fromModel = PHASE_SOURCE === "model";

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
            <InlineField label="Starting phases from">
              <CCP4i2TaskElement itemName="PHASE_SOURCE" {...props}
                qualifiers={{ guiLabel: " " }} />
            </InlineField>
            {/* input phases shown only when not calculating from the model */}
            <CCP4i2TaskElement itemName="ABCD" {...props}
              qualifiers={{ guiLabel: "Starting phases" }}
              visibility={() => !fromModel} />
            {fromModel && (
              <Typography variant="body2" sx={{ pl: 1, color: "text.secondary" }}>
                Phases will be calculated from the model with servalcat sigmaa
                (bulk solvent + sigmaA weighting).
              </Typography>
            )}
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
            qualifiers={{ guiLabel: "NCS assembly" }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" sx={{ mb: 1 }}>
              One row per NCS copy; the first row is the reference. Each row maps
              entities to chains as space-separated <code>role=chain</code>{" "}
              tokens (e.g. <code>CDK=A cyclin=B</code>). For a homomer just give a
              bare chain id per row (<code>A</code>, <code>B</code>, …) — they
              share a single implicit role. Omit a role from a row for a partial
              copy (A<sub>m</sub>B<sub>n</sub>). Leave the list empty to
              auto-detect a homomer (first protein chain = reference, the rest
              its copies).
            </Typography>
            <CCP4i2TaskElement itemName="ASSEMBLY" {...props}
              qualifiers={{ guiLabel: " " }} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Rigid bodies (segments + averaging mode)" }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" sx={{ mb: 1 }}>
              One row per rigid body. <b>segments</b> is a comma-separated list of{" "}
              <code>role:first-last</code> residue ranges that move together —{" "}
              <code>340-485</code> for a homomer, or a cross-chain body such as{" "}
              <code>cyclin:10-95,CDK:45-60</code> (the CDK C-helix travelling with
              the cyclin N-lobe). Each body is superposed independently across the
              copies, so different bodies may follow different NCS transformations.
              Mode: <b>average</b> (rigid NCS), <b>refine</b> (let dm refine the
              operators), or <b>exclude</b>.
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
