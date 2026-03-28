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
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

/**
 * Task interface for shelxeMR — Model building from Molecular Replacement
 * solution using SHELXE.
 *
 * Parameters from shelxeMR.def.xml:
 *   F_SIGF, FREERFLAG, XYZIN, FSOLVENT,
 *   NTCYCLES, NMCYCLES, TIMEFAC,
 *   SALPHELICE, SBETA, SANTIBETA, PRUNRES, USENFOLD
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Input data
  useTaskItem("F_SIGF");
  useTaskItem("FREERFLAG");
  useTaskItem("XYZIN");
  useTaskItem("FSOLVENT");

  // Run options
  useTaskItem("NTCYCLES");
  useTaskItem("NMCYCLES");
  useTaskItem("TIMEFAC");
  useTaskItem("SALPHELICE");
  useTaskItem("SBETA");
  useTaskItem("SANTIBETA");
  useTaskItem("PRUNRES");
  useTaskItem("USENFOLD");

  return (
    <Paper>
      <CCP4i2Tabs>
        {/* Tab 1: Input Data and Run Parameters */}
        <CCP4i2Tab label="Input Data and Run Parameters" key="input">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Select input data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="FREERFLAG"
              qualifiers={{ guiLabel: "Free R set" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Model used for Molecular Replacement",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN"
              qualifiers={{ guiLabel: "Atomic model" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Solvent Content" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="FSOLVENT"
              qualifiers={{
                guiLabel:
                  "Fraction of solvent in the asymmetric unit:",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* Tab 2: Run Options */}
        <CCP4i2Tab label="Run Options" key="options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Shelxe Run Options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="NTCYCLES"
              qualifiers={{
                guiLabel: "Number of shelxe tracing cycles:",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="NMCYCLES"
              qualifiers={{
                guiLabel:
                  "Number of density modification cycles (per trace cycle):",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="TIMEFAC"
              qualifiers={{
                guiLabel:
                  "Time factor for helix and peptide search (Recommended setting: 1-4",
              }}
            />
            <Typography
              variant="body2"
              sx={{ pl: 2, mb: 1, color: "text.secondary" }}
            >
              , where 1 is the quickest, up to 4 with increasingly
              thorough, but slower, searches)
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="SALPHELICE"
              qualifiers={{ guiLabel: "Search for Alpha Helices" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="SBETA"
              qualifiers={{
                guiLabel: "Search for Beta Sheets (parallel)",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="SANTIBETA"
              qualifiers={{
                guiLabel: "Search for Beta Sheets (anti-parallel)",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="PRUNRES"
              qualifiers={{
                guiLabel:
                  "Optimize correlation coefficient for input model:",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="USENFOLD"
              qualifiers={{
                guiLabel:
                  "Apply n-fold symmetry (NC) to C-alpha traces:",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
