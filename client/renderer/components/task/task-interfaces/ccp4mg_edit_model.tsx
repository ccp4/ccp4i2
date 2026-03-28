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
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const searchPdb = useBoolToggle(useTaskItem, "SEARCH_PDB");
  const searchAfdb = useBoolToggle(useTaskItem, "SEARCH_AFDB");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Single folder — no tabs needed */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary">
          Sequences from AU content:
        </Typography>
        <CCP4i2TaskElement itemName="ASUIN" {...props} />
      </CCP4i2ContainerElement>

      {/* --- Model databases --- */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Model databases" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary">
          Databases to search for possible search models
        </Typography>

        <CCP4i2TaskElement itemName="SEARCH_PDB" {...props} qualifiers={{ guiLabel: "Search PDB for possible MR search models" }} />
        {searchPdb.value && (
          <>
            <Typography variant="body2" color="text.secondary" sx={{ pl: 3 }}>
              Non-redundancy level for homologue search:
            </Typography>
            <InlineField label="" sx={{ pl: 3 }}>
              <CCP4i2TaskElement itemName="REDUNDANCYLEVEL" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </>
        )}

        <CCP4i2TaskElement itemName="SEARCH_AFDB" {...props} qualifiers={{ guiLabel: "Search EBI-AFDB for possible MR search models" }} />
        {searchAfdb.value && (
          <>
            <Typography variant="body2" color="text.secondary" sx={{ pl: 3 }}>
              EBI-AFDB pLDDT residue score cut-off:
            </Typography>
            <InlineField label="" sx={{ pl: 3 }}>
              <CCP4i2TaskElement itemName="AFDBLEVEL" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </>
        )}

        <Typography variant="body2" color="text.secondary">
          Maximum no. of search models to create:
        </Typography>
        <CCP4i2TaskElement itemName="MRMAX" {...props} />
      </CCP4i2ContainerElement>

      {/* --- Optional Settings --- */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Optional Settings" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary">
          HHpred results and path to local PDB mirror
        </Typography>
        <CCP4i2TaskElement itemName="HHPREDIN" {...props} />
        <CCP4i2TaskElement itemName="PDBLOCAL" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
