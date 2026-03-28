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

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: REFERENCE } = useTaskItem("REFERENCE");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Reflection objects to manipulate" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Spacegroup and indexing" }}
        containerHint="FolderLevel"
      >
        {REFERENCE !== "ANALYSE" && REFERENCE !== "EXPAND" && (
          <InlineField label="Define new indexing and spacegroup using">
            <CCP4i2TaskElement
              itemName="REFERENCE"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>
        )}
        {REFERENCE === "ANALYSE" && (
          <InlineField label="Analyse data symmetry">
            <CCP4i2TaskElement
              itemName="REFERENCE"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>
        )}
        {REFERENCE === "EXPAND" && (
          <InlineField label="Expand to space group P1">
            <CCP4i2TaskElement
              itemName="REFERENCE"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>
        )}

        {REFERENCE === "HKLIN_FOBS_REF" && (
          <CCP4i2TaskElement itemName="HKLIN_FOBS_REF" {...props} />
        )}
        {REFERENCE === "HKLIN_FC_REF" && (
          <CCP4i2TaskElement itemName="HKLIN_FC_REF" {...props} />
        )}
        {REFERENCE === "HKLIN_FMAP_REF" && (
          <CCP4i2TaskElement itemName="HKLIN_FMAP_REF" {...props} />
        )}
        {REFERENCE === "XYZIN_REF" && (
          <CCP4i2TaskElement itemName="XYZIN_REF" {...props} />
        )}

        {REFERENCE === "SPECIFY" && (
          <>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Note that reindexing a merged file is only valid within the same
              point group
            </Typography>
            <CCP4i2TaskElement itemName="CHOOSE_SPACEGROUP" {...props} />
            <CCP4i2TaskElement itemName="REINDEX_OPERATOR" {...props} />
            <CCP4i2TaskElement
              itemName="USE_REINDEX"
              {...props}
              qualifiers={{ guiLabel: "use reindex operator" }}
            />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              You can give either a spacegroup, or a reindex operator (eg
              k,l,h), or both
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              If only one of these is given, Pointless is usually able to
              generate the other automatically, but you should check
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              The validity and consistency will be checked by Pointless: note
              warnings
            </Typography>
          </>
        )}

        {REFERENCE === "LATTICE" && (
          <>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Remove centred lattice absences: BEWARE dangerous
            </Typography>
            <InlineField label="Desired lattice centering type">
              <CCP4i2TaskElement
                itemName="LATTICE_CENTERING"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              This option should ONLY be used if you are sure that the wrong
              cell was used in integration
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Note that not all centred lattices are consistent with all Bravais
              lattices, check the result carefully
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Lattice type &quot;P&quot; is ignored
            </Typography>
          </>
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
