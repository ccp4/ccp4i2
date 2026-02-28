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

  // Input Data tab
  const { value: ACORN_PHSIN_TYPE } = useTaskItem("ACORN_PHSIN_TYPE");

  // Advanced tab - boolean toggles
  const extend = useBoolToggle(useTaskItem, "ACORN_EXTEND");
  const bgrid = useBoolToggle(useTaskItem, "ACORN_BGRID");
  const bseed = useBoolToggle(useTaskItem, "ACORN_BSEED");
  const bresol = useBoolToggle(useTaskItem, "ACORN_BRESOL");
  const bexclude = useBoolToggle(useTaskItem, "ACORN_BEXCLUDE");
  const becut = useBoolToggle(useTaskItem, "ACORN_BECUT");
  const customPhase = useBoolToggle(useTaskItem, "ACOPH_CUSTOM");
  const custddm = useBoolToggle(useTaskItem, "ACOPH_CUSTDDM");
  const peaksearch = useBoolToggle(useTaskItem, "ACOMPS_PEAKSEARCH");

  const { value: ACOPH_TRIALS } = useTaskItem("ACOPH_TRIALS");
  const numTrials =
    typeof ACOPH_TRIALS === "number"
      ? ACOPH_TRIALS
      : parseInt(ACOPH_TRIALS as string) || 1;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input Data ===== */}
        <CCP4i2Tab label="Input Data">
          <InlineField label="Run ACORN with" width="auto" sx={{ mb: 1 }}>
            <CCP4i2TaskElement
              itemName="ACORN_PHSIN_TYPE"
              {...props}
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
          </InlineField>

          {/* Reflection Data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflection Data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <CCP4i2TaskElement
              itemName="ABCD"
              {...props}
              qualifiers={{ guiLabel: "Phases" }}
              visibility={() => ACORN_PHSIN_TYPE === "phases"}
            />
          </CCP4i2ContainerElement>

          {/* Model for approximate co-ordinates (model mode only) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model for approximate co-ordinates" }}
            containerHint="FolderLevel"
            visibility={() => ACORN_PHSIN_TYPE === "model"}
          >
            <CCP4i2TaskElement
              itemName="XYZIN"
              {...props}
              qualifiers={{ guiLabel: "Atomic model" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Tab 2: Advanced Acorn Parameters ===== */}
        <CCP4i2Tab label="Advanced Acorn Parameters">
          {/* Data preparation */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Data preparation" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="ACORN_ECALC"
              {...props}
              qualifiers={{
                guiLabel: "Calculate normalised E values using clipper",
              }}
            />
            <CCP4i2TaskElement
              itemName="ACORN_ANISOTROPY"
              {...props}
              qualifiers={{ guiLabel: "Correct for anisotropy" }}
            />
            <CCP4i2TaskElement
              itemName="ACORN_EXTEND"
              {...props}
              qualifiers={{ guiLabel: "Extend data to high resolution" }}
              onChange={extend.onChange}
            />
            {extend.value && (
              <InlineField label="Resolution limit:" hint="Å" sx={{ pl: 3 }}>
                <CCP4i2TaskElement
                  itemName="ACORN_EXTENDRES"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}
          </CCP4i2ContainerElement>

          {/* Refinement parameters */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Refinement parameters" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="ACORN_BGRID"
              {...props}
              qualifiers={{ guiLabel: "Set grid spacing" }}
              onChange={bgrid.onChange}
            />
            {bgrid.value && (
              <InlineField label="Grid spacing:" sx={{ pl: 3 }}>
                <CCP4i2TaskElement
                  itemName="ACOGEN_GRID"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}

            <CCP4i2TaskElement
              itemName="ACORN_BSEED"
              {...props}
              qualifiers={{ guiLabel: "Set random seed" }}
              onChange={bseed.onChange}
            />
            {bseed.value && (
              <InlineField label="Seed value:" sx={{ pl: 3 }}>
                <CCP4i2TaskElement
                  itemName="ACOGEN_SEED"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}

            <CCP4i2TaskElement
              itemName="ACORN_BRESOL"
              {...props}
              qualifiers={{ guiLabel: "Set resolution limits" }}
              onChange={bresol.onChange}
            />
            {bresol.value && (
              <InlineField
                label="Low:"
                sx={{ pl: 3 }}
                after={
                  <InlineField label="High:" hint="Å">
                    <CCP4i2TaskElement
                      itemName="ACOREF_RESOLU"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                }
              >
                <CCP4i2TaskElement
                  itemName="ACOREF_RESOLL"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}

            <CCP4i2TaskElement
              itemName="ACORN_BEXCLUDE"
              {...props}
              qualifiers={{ guiLabel: "Exclude reflections below E value" }}
              onChange={bexclude.onChange}
            />
            {bexclude.value && (
              <InlineField label="Exclude below E =" sx={{ pl: 3 }}>
                <CCP4i2TaskElement
                  itemName="ACOREF_EXCLUDE"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}

            <CCP4i2TaskElement
              itemName="ACORN_BECUT"
              {...props}
              qualifiers={{ guiLabel: "Set E value cutoff for supergrid" }}
              onChange={becut.onChange}
            />
            {becut.value && (
              <InlineField label="E cutoff:" sx={{ pl: 3 }}>
                <CCP4i2TaskElement
                  itemName="ACOREF_ECUT"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}
          </CCP4i2ContainerElement>

          {/* Phase refinement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Phase refinement" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="ACOPH_PATSUP"
              {...props}
              qualifiers={{ guiLabel: "Use Patterson superposition" }}
            />
            <InlineField label="DDM cutoff:">
              <CCP4i2TaskElement
                itemName="ACOPH_CUTDDM"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="PS finish threshold:">
              <CCP4i2TaskElement
                itemName="ACOPH_PSFINISH"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="CC finish threshold:">
              <CCP4i2TaskElement
                itemName="ACOPH_CCFINISH"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>

            {/* Custom phase refinement protocol */}
            <CCP4i2TaskElement
              itemName="ACOPH_CUSTOM"
              {...props}
              qualifiers={{ guiLabel: "Use custom phase refinement protocol" }}
              onChange={customPhase.onChange}
            />
            {customPhase.value && (
              <Box sx={{ pl: 3, display: "flex", flexDirection: "column", gap: 1 }}>
                <CCP4i2TaskElement
                  itemName="ACOPH_CUSTDDM"
                  {...props}
                  qualifiers={{ guiLabel: "Use custom DDM parameters" }}
                  onChange={custddm.onChange}
                />
                {custddm.value && (
                  <Box
                    sx={{
                      pl: 3,
                      display: "flex",
                      flexDirection: "column",
                      gap: 1,
                    }}
                  >
                    <InlineField label="Number of trials:">
                      <CCP4i2TaskElement
                        itemName="ACOPH_TRIALS"
                        {...props}
                        qualifiers={{ guiLabel: " " }}
                      />
                    </InlineField>
                    {Array.from({ length: numTrials }, (_, i) => i + 1).map(
                      (n) => (
                        <Box
                          key={n}
                          sx={{
                            display: "flex",
                            alignItems: "center",
                            gap: 1,
                            flexWrap: "wrap",
                          }}
                        >
                          <Typography
                            variant="body2"
                            sx={{ fontWeight: "bold", minWidth: "4rem" }}
                          >
                            Trial {n}:
                          </Typography>
                          <Box sx={{ width: "8rem" }}>
                            <CCP4i2TaskElement
                              itemName={`ACOPH_NCDDM_${n}`}
                              {...props}
                              qualifiers={{ guiLabel: "Cycles" }}
                            />
                          </Box>
                          <Box sx={{ width: "8rem" }}>
                            <CCP4i2TaskElement
                              itemName={`ACOPH_DDMK_${n}`}
                              {...props}
                              qualifiers={{ guiLabel: "DDM type" }}
                            />
                          </Box>
                          <Box sx={{ width: "14rem" }}>
                            <CCP4i2TaskElement
                              itemName={`ACOPH_REFINE_${n}`}
                              {...props}
                              qualifiers={{ guiLabel: "Refinement" }}
                            />
                          </Box>
                        </Box>
                      )
                    )}
                  </Box>
                )}
              </Box>
            )}
          </CCP4i2ContainerElement>

          {/* Map peak search */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Map peak search" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="ACOMPS_PEAKSEARCH"
              {...props}
              qualifiers={{ guiLabel: "Perform peak search on output map" }}
              onChange={peaksearch.onChange}
            />
            {peaksearch.value && (
              <Box
                sx={{
                  pl: 3,
                  display: "flex",
                  flexDirection: "column",
                  gap: 1,
                }}
              >
                <InlineField label="Maximum peaks:">
                  <CCP4i2TaskElement
                    itemName="ACOMPS_MAXPEAKS"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
                <InlineField label="RMS multiplier:">
                  <CCP4i2TaskElement
                    itemName="ACOMPS_RMSMULT"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              </Box>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
