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
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          {/* Atomic model */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Atomic model" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          </CCP4i2ContainerElement>

          {/* Options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Options" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Select B-factor treatment option - it is important this is set
              correctly
            </Typography>
            <CCP4i2TaskElement
              itemName="BTREATMENT"
              {...props}
              qualifiers={{ guiLabel: "" }}
            />
            <Typography variant="body2" color="text.secondary" sx={{ mt: 1, mb: 1 }}>
              Options for residues and regions
            </Typography>
            <CCP4i2TaskElement
              itemName="CONFCUT"
              {...props}
              qualifiers={{
                guiLabel:
                  "Remove low confidence residues (lddt or rmsd) (recommended)",
              }}
            />
            <CCP4i2TaskElement
              itemName="COMPACTREG"
              {...props}
              qualifiers={{
                guiLabel: "Split model into compact regions (recommended)",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Optional Files */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Optional Files" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              PAE file
            </Typography>
            <CCP4i2TaskElement itemName="PAEIN" {...props} />
            <Typography variant="body2" color="text.secondary" sx={{ mt: 1, mb: 1 }}>
              Distance model file
            </Typography>
            <CCP4i2TaskElement itemName="XYZDISTMOD" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="addSettings" label="Additional Settings">
          {/* Options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Options" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Maximum domains to obtain">
              <CCP4i2TaskElement
                itemName="MAXDOM"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Approximate size of domains to be found (Angstroms)">
              <CCP4i2TaskElement
                itemName="DOMAINSIZE"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Minimum domain length (residues)">
              <CCP4i2TaskElement
                itemName="MINDOML"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Maximum fraction close">
              <CCP4i2TaskElement
                itemName="MAXFRACCL"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Minimum sequential residues">
              <CCP4i2TaskElement
                itemName="MINSEQRESI"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Minimum remainder sequence length">
              <CCP4i2TaskElement
                itemName="MINREMSEQL"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Minimum LDDT for removing residues">
              <CCP4i2TaskElement
                itemName="MINLDDT"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Maximum RMSD for removing residues">
              <CCP4i2TaskElement
                itemName="MAXRMSD"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <CCP4i2TaskElement
              itemName="SUBMINB"
              {...props}
              qualifiers={{
                guiLabel: "Subtract the lowest B-value from all B-values",
              }}
            />
          </CCP4i2ContainerElement>

          {/* PAE File Options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "PAE File Options (if PAE matrix supplied)" }}
            containerHint="FolderLevel"
          >
            <InlineField label="PAE power">
              <CCP4i2TaskElement
                itemName="PAEPOWER"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="PAE cutoff">
              <CCP4i2TaskElement
                itemName="PAECUTOFF"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="PAE graph resolution">
              <CCP4i2TaskElement
                itemName="PAEGRAPHRES"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Distance Model Options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Distance Model Options (if suitable model provided)",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="WEIGHTCA"
              {...props}
              qualifiers={{
                guiLabel:
                  "Weight by CA-CA distance (if distance_model supplied)",
              }}
            />
            <InlineField label="Distance power (for weighting by CA-CA distance)">
              <CCP4i2TaskElement
                itemName="DISTPOW"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
