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
import React from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Enum values for conditional visibility
  const { value: EXISTING_MODELS } = useTaskItem("AMPLE_EXISTING_MODELS");
  const { value: MODELS_SOURCE } = useTaskItem("AMPLE_MODELS_SOURCE");
  const { value: MODEL_TYPE } = useTaskItem("AMPLE_MODEL_TYPE");
  const { value: MODEL_GENERATION } = useTaskItem("AMPLE_MODEL_GENERATION");

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input data ===== */}
        <CCP4i2Tab label="Input data">
          {/* Input sequence */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input sequence" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="AMPLE_SEQIN"
              {...props}
              qualifiers={{ guiLabel: "Sequence" }}
            />
          </CCP4i2ContainerElement>

          {/* Input reflections */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input reflections" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="AMPLE_F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections" }}
            />
          </CCP4i2ContainerElement>

          {/* Class of protein */}
          <CCP4i2TaskElement
            itemName="AMPLE_PROTEIN_CLASS"
            {...props}
            qualifiers={{ guiLabel: "Class of protein", guiMode: "radio" }}
          />

          {/* Do you have existing models? */}
          <CCP4i2TaskElement
            itemName="AMPLE_EXISTING_MODELS"
            {...props}
            qualifiers={{ guiLabel: "Do you have existing models?" }}
          />

          {/* When existing models = True: model source, dir/file, type */}
          {EXISTING_MODELS === "True" && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                itemName="AMPLE_MODELS_SOURCE"
                {...props}
                qualifiers={{ guiLabel: "Models are from", guiMode: "radio" }}
              />
              <CCP4i2TaskElement
                itemName="AMPLE_MODELS_DIR"
                {...props}
                qualifiers={{ guiLabel: "Models directory" }}
                visibility={() => MODELS_SOURCE === "directory"}
              />
              <CCP4i2TaskElement
                itemName="AMPLE_MODELS_FILE"
                {...props}
                qualifiers={{ guiLabel: "Models file" }}
                visibility={() => MODELS_SOURCE === "file"}
              />
              <CCP4i2TaskElement
                itemName="AMPLE_MODEL_TYPE"
                {...props}
                qualifiers={{ guiLabel: "What sort of models are these?" }}
              />
              <CCP4i2TaskElement
                itemName="AMPLE_NMR_REMODEL"
                {...props}
                qualifiers={{ guiLabel: "NMR remodelling" }}
                visibility={() => MODEL_TYPE === "nmr_ensemble"}
              />
            </CCP4i2ContainerElement>
          )}

          {/* When existing models = False: model generation */}
          {EXISTING_MODELS === "False" && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                itemName="AMPLE_MODEL_GENERATION"
                {...props}
                qualifiers={{ guiLabel: "Model generation" }}
              />
              {MODEL_GENERATION === "rosetta" && (
                <>
                  <CCP4i2TaskElement
                    itemName="AMPLE_ROSETTA_DIR"
                    {...props}
                    qualifiers={{ guiLabel: "ROSETTA directory" }}
                  />
                  <CCP4i2TaskElement
                    itemName="AMPLE_ROSETTA_FRAGS3"
                    {...props}
                    qualifiers={{ guiLabel: "ROSETTA 3-fragment database" }}
                  />
                  <CCP4i2TaskElement
                    itemName="AMPLE_ROSETTA_FRAGS9"
                    {...props}
                    qualifiers={{ guiLabel: "ROSETTA 9-fragment database" }}
                  />
                  <CCP4i2TaskElement
                    itemName="AMPLE_CONTACT_FILE"
                    {...props}
                    qualifiers={{ guiLabel: "Inter-residue restraints" }}
                  />
                  <CCP4i2TaskElement
                    itemName="AMPLE_CONTACT_FORMAT"
                    {...props}
                    qualifiers={{ guiLabel: "Contact file format" }}
                  />
                </>
              )}
            </CCP4i2ContainerElement>
          )}

          {/* Number of processors */}
          <CCP4i2TaskElement
            itemName="AMPLE_NPROC"
            {...props}
            qualifiers={{ guiLabel: "Number of processors" }}
          />
        </CCP4i2Tab>

        {/* ===== Tab 2: Advanced Options ===== */}
        <CCP4i2Tab label="Advanced Options">
          {/* Rebuilding options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Rebuilding options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="AMPLE_REFINE_REBUILD"
              {...props}
              qualifiers={{ guiLabel: "Rebuild REFMAC-refined MR result" }}
            />
            <CCP4i2TaskElement
              itemName="AMPLE_USE_SHELXE"
              {...props}
              qualifiers={{ guiLabel: "Run SHELXE after MR" }}
            />
            <CCP4i2TaskElement
              itemName="AMPLE_SHELXE_REBUILD"
              {...props}
              qualifiers={{ guiLabel: "Rebuild SHELXE traces" }}
            />
          </CCP4i2ContainerElement>

          {/* Ensembling options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Ensembling options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="AMPLE_ENSEMBLING_TM"
              {...props}
              qualifiers={{ guiLabel: "Better but slower ensembling" }}
            />
          </CCP4i2ContainerElement>

          {/* Extra command-line arguments */}
          <CCP4i2TaskElement
            itemName="AMPLE_EXTRA_FLAGS"
            {...props}
            qualifiers={{ guiLabel: "Extra command-line arguments" }}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
