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
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { InlineField } from "../task-elements/inline-field";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Mode and run location
  const { value: ARCIMBOLDO_OPTIONS } = useTaskItem("ARCIMBOLDO_OPTIONS");
  const { value: ARCIMBOLDO_RUN } = useTaskItem("ARCIMBOLDO_RUN");
  const { value: RUN_MODE } = useTaskItem("RUN_MODE");

  // LITE mode
  const { value: LITE_MODELS } = useTaskItem("LITE_MODELS");
  const litePartial = useBoolToggle(useTaskItem, "LITE_PARTIAL");

  // BORGES mode
  const { value: BORGES_LIBRARY } = useTaskItem("BORGES_LIBRARY");

  // SHREDDER mode
  const { value: SHREDDER_OPTIONS } = useTaskItem("SHREDDER_OPTIONS");

  // Developer options
  const { value: DEVELOPER_MODE } = useTaskItem("DEVELOPER_MODE");

  const isLite = ARCIMBOLDO_OPTIONS === "LITE";
  const isBorges = ARCIMBOLDO_OPTIONS === "BORGES";
  const isShredder = ARCIMBOLDO_OPTIONS === "SHREDDER";
  const isGrid =
    ARCIMBOLDO_RUN === "local_grid" || ARCIMBOLDO_RUN === "remote_grid";

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input Data ===== */}
        <CCP4i2Tab label="Input data">
          {/* Mode and run location */}
          <InlineField
            label="Run ARCIMBOLDO"
            width="12rem"
            after={
              <InlineField label="on" width="12rem">
                <CCP4i2TaskElement
                  itemName="ARCIMBOLDO_RUN"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            }
          >
            <CCP4i2TaskElement
              itemName="ARCIMBOLDO_OPTIONS"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>

          <CCP4i2TaskElement
            itemName="COIL_COILED"
            {...props}
            qualifiers={{ guiLabel: "Run in coil coiled mode" }}
          />

          {isShredder && (
            <CCP4i2TaskElement
              itemName="SHREDDER_PREDICTED"
              {...props}
              qualifiers={{ guiLabel: "Run in predicted model mode" }}
            />
          )}

          {/* Grid configuration - visible in grid modes */}
          {isGrid && (
            <Box
              sx={{ pl: 3, display: "flex", flexDirection: "column", gap: 1 }}
            >
              <CCP4i2TaskElement
                itemName="RUN_MODE"
                {...props}
                qualifiers={{ guiLabel: "Grid configuration" }}
              />
              {RUN_MODE === "CUSTOM" && (
                <CCP4i2TaskElement
                  itemName="CONFIG_FILE"
                  {...props}
                  qualifiers={{ guiLabel: "Configuration file" }}
                />
              )}
            </Box>
          )}

          {/* Input data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <InlineField
              label="Asymmetric unit contains"
              width="6rem"
              after={
                <InlineField
                  label="components of molecular weight"
                  width="10rem"
                  hint="Daltons"
                >
                  <CCP4i2TaskElement
                    itemName="MOLECULAR_WEIGHT"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              }
            >
              <CCP4i2TaskElement
                itemName="N_COMPONENTS"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* ---- LITE Model section ---- */}
          {isLite && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Model" }}
              containerHint="FolderLevel"
            >
              <CCP4i2TaskElement
                itemName="LITE_MODELS"
                {...props}
                qualifiers={{ guiLabel: "Search model type" }}
              />

              <InlineField
                label="Expected r.m.s.d. from target"
                hint="Å"
              >
                <CCP4i2TaskElement
                  itemName="LITE_RMSD"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>

              {LITE_MODELS === "HELIX" && (
                <InlineField
                  label="Search for"
                  width="6rem"
                  after={
                    <InlineField
                      label="copies of helix with length"
                      width="6rem"
                      hint="residues"
                    >
                      <CCP4i2TaskElement
                        itemName="HELIX_LENGTH"
                        {...props}
                        qualifiers={{ guiLabel: " " }}
                      />
                    </InlineField>
                  }
                >
                  <CCP4i2TaskElement
                    itemName="N_FRAGMENTS"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              )}

              {LITE_MODELS === "CUSTOM" && (
                <>
                  <InlineField
                    label="Search for"
                    width="6rem"
                    hint="copies of custom model"
                  >
                    <CCP4i2TaskElement
                      itemName="N_FRAGMENTS"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                  <CCP4i2TaskElement
                    itemName="PDB_LITE"
                    {...props}
                    qualifiers={{ guiLabel: "Custom model" }}
                  />
                </>
              )}

              {LITE_MODELS === "HELICES" && (
                <CCP4i2TaskElement
                  itemName="LITE_HELICES_LIST"
                  {...props}
                  qualifiers={{ guiLabel: "Helix lengths" }}
                />
              )}

              {LITE_MODELS === "CUSTOMS" && (
                <CCP4i2TaskElement
                  itemName="LITE_CUSTOMS_LIST"
                  {...props}
                  qualifiers={{ guiLabel: "Custom models" }}
                />
              )}

              <CCP4i2TaskElement
                itemName="LITE_PARTIAL"
                {...props}
                qualifiers={{ guiLabel: "Start from known partial structure" }}
                onChange={litePartial.onChange}
              />
              {litePartial.value && (
                <Box sx={{ pl: 3 }}>
                  <CCP4i2TaskElement
                    itemName="LITE_FIXED"
                    {...props}
                    qualifiers={{
                      guiLabel: "Fixed partial structure directory",
                    }}
                  />
                </Box>
              )}
            </CCP4i2ContainerElement>
          )}

          {/* ---- BORGES Library section ---- */}
          {isBorges && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Borges libraries and model handling" }}
              containerHint="FolderLevel"
            >
              <CCP4i2TaskElement
                itemName="BORGES_LIBRARY"
                {...props}
                qualifiers={{
                  guiLabel: "Library (topology: u - up, d - down)",
                }}
              />
              {BORGES_LIBRARY === "CUSTOM" && (
                <Box sx={{ pl: 3 }}>
                  <CCP4i2TaskElement
                    itemName="BORGES_CUSTOM"
                    {...props}
                    qualifiers={{ guiLabel: "Custom library directory" }}
                  />
                </Box>
              )}
              <CCP4i2TaskElement
                itemName="BORGES_GYRE_T"
                {...props}
                qualifiers={{ guiLabel: "Phaser GYRE option" }}
              />
              <CCP4i2TaskElement
                itemName="BORGES_GIMBLE_T"
                {...props}
                qualifiers={{ guiLabel: "Phaser GIMBLE option" }}
              />
              <CCP4i2TaskElement
                itemName="BORGES_MULTICOPY_T"
                {...props}
                qualifiers={{ guiLabel: "MULTICOPY option" }}
              />
            </CCP4i2ContainerElement>
          )}

          {/* ---- SHREDDER Model section ---- */}
          {isShredder && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Shredder models" }}
              containerHint="FolderLevel"
            >
              <CCP4i2TaskElement
                itemName="PDB_SHREDDER"
                {...props}
                qualifiers={{ guiLabel: "Model to shred" }}
              />

              <InlineField label="Expected r.m.s.d." hint="Å">
                <CCP4i2TaskElement
                  itemName="SHREDDER_RMSD_T"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>

              <CCP4i2TaskElement
                itemName="SHREDDER_CONVERT_T"
                {...props}
                qualifiers={{ guiLabel: "Convert to polyalanine" }}
              />
              <CCP4i2TaskElement
                itemName="SHREDDER_MAKE_T"
                {...props}
                qualifiers={{ guiLabel: "Make all B-factors equal" }}
              />
              <CCP4i2TaskElement
                itemName="SHREDDER_OPTIONS"
                {...props}
                qualifiers={{ guiLabel: "Shredder mode" }}
              />

              {SHREDDER_OPTIONS === "spherical" && (
                <>
                  <CCP4i2TaskElement
                    itemName="SHREDDER_COIL_T"
                    {...props}
                    qualifiers={{ guiLabel: "Maintain coil in the model" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SHREDDER_GYRE_T"
                    {...props}
                    qualifiers={{ guiLabel: "Perform gyre refinement" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SHREDDER_GIMBLE_T"
                    {...props}
                    qualifiers={{ guiLabel: "Perform gimble refinement" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SHREDDER_LLG_T"
                    {...props}
                    qualifiers={{ guiLabel: "Perform LLG-guided pruning" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SHREDDER_COMBINE_T"
                    {...props}
                    qualifiers={{ guiLabel: "Combine phases with alixe" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SHREDDER_MULTICOPY_T"
                    {...props}
                    qualifiers={{ guiLabel: "MULTICOPY option" }}
                  />
                </>
              )}
            </CCP4i2ContainerElement>
          )}
        </CCP4i2Tab>

        {/* ===== Tab 2: Advanced Data ===== */}
        <CCP4i2Tab label="Advanced data">
          {/* TNCS option: checkbox + dropdown inline */}
          <InlineField
            width="auto"
            after={
              <InlineField label="Switch Phaser TNCS option" width="8rem">
                <CCP4i2TaskElement
                  itemName="TNCS_T"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            }
          >
            <CCP4i2TaskElement
              itemName="TNCS"
              {...props}
              qualifiers={{ guiLabel: " " }}
              sx={{ width: "auto" }}
            />
          </InlineField>

          {/* shelxe line */}
          <InlineField label="shelxe_line =" width="auto">
            <Box sx={{ flex: 1, minWidth: "20rem" }}>
              <CCP4i2TaskElement
                itemName="SHELXE_LINE"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </Box>
          </InlineField>

          {/* Add lines to bor-file */}
          <CCP4i2TaskElement
            itemName="KEYWORDS"
            {...props}
            qualifiers={{ guiLabel: "Add lines to bor-file" }}
          />
        </CCP4i2Tab>

        {/* ===== Tab 3: Developer Options ===== */}
        <CCP4i2Tab label="Developer options">
          <CCP4i2TaskElement
            itemName="DEVELOPER_MODE"
            {...props}
            qualifiers={{ guiLabel: "Select run mode" }}
          />
          {DEVELOPER_MODE === "EXISTING" && (
            <CCP4i2TaskElement
              itemName="EXISTING_FOLDER"
              {...props}
              qualifiers={{ guiLabel: "Existing run directory" }}
            />
          )}
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
