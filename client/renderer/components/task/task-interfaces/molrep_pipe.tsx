import React, { useCallback } from "react";
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useFreeRWarning } from "../../../providers/run-check-provider";
import { InlineField } from "../task-elements/inline-field";

/**
 * Task interface for the MOLREP molecular replacement pipeline.
 *
 * Pipeline: MOLREP → optional Sheetbend → REFMAC refinement
 *
 * Three tabs: Input Data and Protocol, Basic Options, Advanced Options
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest, createPeerTask, validation } = useJob(
    job.id
  );

  // --- Input data items ---
  const { item: F_SIGFItem } = useTaskItem("F_SIGF");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  // Free R flag validation warning
  useFreeRWarning({
    job,
    taskName: "molrep_pipe",
    freeRFlag,
    validation,
    createPeerTask,
  });

  // --- Values for conditional visibility ---
  const { value: sgOptions } = useTaskItem("controlParameters.SG_OPTIONS");
  const { value: highPathVar } = useTaskItem("controlParameters.HIGH_PATH_VAR");
  const { value: lowPathVar } = useTaskItem("controlParameters.LOW_PATH_VAR");

  // --- F_SIGF onChange: extract wavelength from digest ---
  const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");

  const handleF_SIGFChange = useCallback(async () => {
    if (!F_SIGFItem?._objectPath) return;
    const digestData = await fetchDigest(F_SIGFItem._objectPath);
    const wavelength = digestData?.wavelengths?.at(-1);
    if (wavelength && wavelength > 0 && wavelength < 9) {
      await forceUpdateWAVELENGTH(wavelength);
    }
  }, [F_SIGFItem?._objectPath, fetchDigest, forceUpdateWAVELENGTH]);

  return (
    <Paper>
      <CCP4i2Tabs>
        {/* ============================================================
            TAB 1: INPUT DATA AND PROTOCOL
            ============================================================ */}
        <CCP4i2Tab label="Input Data and Protocol" key="input">
          {/* Experimental data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Experimental data", initiallyOpen: true }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{ guiLabel: "Reflections" }}
              onChange={handleF_SIGFChange}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="FREERFLAG"
              qualifiers={{ guiLabel: "Free R set" }}
            />
          </CCP4i2ContainerElement>

          {/* Search model */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Search model", initiallyOpen: true }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN"
              qualifiers={{ guiLabel: "Atomic model" }}
            />
          </CCP4i2ContainerElement>

          {/* Sequence of target model */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Sequence of target model",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="ASUIN"
              qualifiers={{ guiLabel: "AU contents" }}
            />
            <Typography variant="body2" sx={{ mt: 1, fontStyle: "italic" }}>
              Select one sequence
            </Typography>
            <InlineField
              label="The number of monomers to search for"
              sx={{ mt: 1 }}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="controlParameters.NMON"
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Fixed Model */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Fixed Model", initiallyOpen: true }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_FIX"
              qualifiers={{ guiLabel: "Atomic model" }}
            />
          </CCP4i2ContainerElement>

          {/* Extra Steps */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Extra Steps", initiallyOpen: true }}
            containerHint="FolderLevel"
          >
            <InlineField
              label="Run"
              width="auto"
              after={
                <InlineField
                  label="shift field refinement followed by"
                  hint="cycles of restrained refinement"
                >
                  <CCP4i2TaskElement
                    {...props}
                    itemName="REFMAC_NCYC"
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              }
            >
              <CCP4i2TaskElement
                {...props}
                itemName="RUNSHEETBEND"
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 2: BASIC OPTIONS
            ============================================================ */}
        <CCP4i2Tab label="Basic Options" key="basic">
          {/* Search parameters */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Search parameters",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <InlineField label="Search in space group(s)" width="20rem">
              <CCP4i2TaskElement
                {...props}
                itemName="controlParameters.SG_OPTIONS"
                qualifiers={{
                  guiLabel: " ",
                  menuText: [
                    "as input reflection file",
                    "specify alternative",
                    "all space groups in Laue group of input reflection file",
                  ],
                  enumerators: ["no", "specify", "laue"],
                }}
              />
            </InlineField>
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.SG"
              qualifiers={{ guiLabel: "Space group" }}
              visibility={() => sgOptions === "specify"}
            />
            <InlineField
              label="Use data in resolution range from low"
              after={
                <InlineField label="to high">
                  <CCP4i2TaskElement
                    {...props}
                    itemName="controlParameters.RESMAX"
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              }
            >
              <CCP4i2TaskElement
                {...props}
                itemName="controlParameters.RESMIN"
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Modify search model */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Modify search model",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.SEQ"
              qualifiers={{
                guiLabel:
                  "Perform alignment and use it to rename residues and trim side chains",
                guiMode: "radio",
                menuText: [
                  "always",
                  "only for sequence identity > 20%",
                  "never",
                ],
                enumerators: ["y", "d", "n"],
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.SURF"
              qualifiers={{
                guiLabel: "B-factors modification options",
                guiMode: "radio",
                menuText: [
                  "Increase B-factor on the molecular surface for all functions",
                  "Increase B-factor on the surface for Packing function only",
                  "Do not do anything",
                  "Set B-factors of all atoms to 20",
                  "Use poly-alanine model with all B-factors 20",
                ],
                enumerators: ["y", "c", "n", "2", "a"],
              }}
            />
          </CCP4i2ContainerElement>

          {/* Customise search procedure */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Customise search procedure",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <Typography
              variant="body2"
              sx={{ fontStyle: "italic", mb: 1 }}
            >
              Number of peaks to analyse
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.NP"
              qualifiers={{ guiLabel: "Number of Rotation Function peaks" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.NPT"
              qualifiers={{ guiLabel: "Number of Translation Function peaks" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 3: ADVANCED OPTIONS
            ============================================================ */}
        <CCP4i2Tab label="Advanced Options" key="advanced">
          {/* Scoring putative solutions */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Scoring putative solutions",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" sx={{ mb: 0.5 }}>
              CC = Correlation Coefficient
            </Typography>
            <Typography variant="body2" sx={{ mb: 1 }}>
              PF = Packing Function
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.SCORE"
              qualifiers={{
                guiLabel: " ",
                guiMode: "radio",
                menuText: [
                  "Use CC times PF as a score and stop translation search if contrast is > 3.0",
                  "Use CC times PF and do not stop",
                  "Use CC and do not stop",
                ],
                enumerators: ["y", "n", "c"],
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.NMON_EXP"
              qualifiers={{
                guiLabel:
                  "Expected number of copies (for contrast calculation only)",
              }}
            />
            <Typography
              variant="body2"
              sx={{ fontStyle: "italic", mt: 2, mb: 1 }}
            >
              Scaling
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.ANISO"
              qualifiers={{
                guiLabel: " ",
                guiMode: "radio",
                menuText: ["anisotropic", "isotropic", "none"],
                enumerators: ["y", "n", "k"],
              }}
            />
          </CCP4i2ContainerElement>

          {/* Filter parameters */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Filter parameters",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            {/* High pass filter */}
            <Typography
              variant="body2"
              sx={{ fontStyle: "italic", mb: 1 }}
            >
              High pass filter parameter (B-add, the B-factor applied to
              input structure amplitudes)
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.HIGH_PATH_VAR"
              qualifiers={{
                guiLabel: " ",
                guiMode: "radio",
                menuText: [
                  "From identity between model and sequence (if sequence given)",
                  "From identity specified manually",
                  "From high resolution limit",
                  "Directly as the value of additional B-factor",
                ],
                enumerators: ["s", "i", "r", "b"],
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.SIM"
              qualifiers={{ guiLabel: "Sequence similarity" }}
              visibility={() => highPathVar === "i"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.BADD"
              qualifiers={{ guiLabel: "Additional B-factor" }}
              visibility={() => highPathVar === "b"}
            />

            {/* Low pass filter */}
            <Typography
              variant="body2"
              sx={{ fontStyle: "italic", mt: 2, mb: 1 }}
            >
              Low pass filter parameter (B-off, the B-factor of the removed
              fraction of structure amplitudes)
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.LOW_PATH_VAR"
              qualifiers={{
                guiLabel: " ",
                guiMode: "radio",
                menuText: [
                  "From completeness of the search model",
                  "From low resolution limit",
                  "Directly as the value of additional B-factor",
                ],
                enumerators: ["c", "r", "b"],
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.BOFF"
              qualifiers={{ guiLabel: "Additional B-factor" }}
              visibility={() => lowPathVar === "b"}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
