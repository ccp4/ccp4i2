import React, { useCallback, useEffect, useState } from "react";
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Input Data tab
  const { value: ACORN_PHSIN_TYPE } = useTaskItem("ACORN_PHSIN_TYPE");

  // Advanced tab - boolean toggles
  const { value: ACORN_EXTEND_RAW } = useTaskItem("ACORN_EXTEND");
  const { value: ACORN_BGRID_RAW } = useTaskItem("ACORN_BGRID");
  const { value: ACORN_BSEED_RAW } = useTaskItem("ACORN_BSEED");
  const { value: ACORN_BRESOL_RAW } = useTaskItem("ACORN_BRESOL");
  const { value: ACORN_BEXCLUDE_RAW } = useTaskItem("ACORN_BEXCLUDE");
  const { value: ACORN_BECUT_RAW } = useTaskItem("ACORN_BECUT");
  const { value: ACOPH_CUSTOM_RAW } = useTaskItem("ACOPH_CUSTOM");
  const { value: ACOPH_CUSTDDM_RAW } = useTaskItem("ACOPH_CUSTDDM");
  const { value: ACOMPS_PEAKSEARCH_RAW } = useTaskItem("ACOMPS_PEAKSEARCH");
  const { value: ACOPH_TRIALS } = useTaskItem("ACOPH_TRIALS");

  // Local state for boolean-driven visibility
  const [extend, setExtend] = useState(() => isTruthy(ACORN_EXTEND_RAW));
  const [bgrid, setBgrid] = useState(() => isTruthy(ACORN_BGRID_RAW));
  const [bseed, setBseed] = useState(() => isTruthy(ACORN_BSEED_RAW));
  const [bresol, setBresol] = useState(() => isTruthy(ACORN_BRESOL_RAW));
  const [bexclude, setBexclude] = useState(() => isTruthy(ACORN_BEXCLUDE_RAW));
  const [becut, setBecut] = useState(() => isTruthy(ACORN_BECUT_RAW));
  const [customPhase, setCustomPhase] = useState(() => isTruthy(ACOPH_CUSTOM_RAW));
  const [custddm, setCustddm] = useState(() => isTruthy(ACOPH_CUSTDDM_RAW));
  const [peaksearch, setPeaksearch] = useState(() =>
    isTruthy(ACOMPS_PEAKSEARCH_RAW)
  );

  // Sync from server
  useEffect(() => setExtend(isTruthy(ACORN_EXTEND_RAW)), [ACORN_EXTEND_RAW]);
  useEffect(() => setBgrid(isTruthy(ACORN_BGRID_RAW)), [ACORN_BGRID_RAW]);
  useEffect(() => setBseed(isTruthy(ACORN_BSEED_RAW)), [ACORN_BSEED_RAW]);
  useEffect(() => setBresol(isTruthy(ACORN_BRESOL_RAW)), [ACORN_BRESOL_RAW]);
  useEffect(
    () => setBexclude(isTruthy(ACORN_BEXCLUDE_RAW)),
    [ACORN_BEXCLUDE_RAW]
  );
  useEffect(() => setBecut(isTruthy(ACORN_BECUT_RAW)), [ACORN_BECUT_RAW]);
  useEffect(
    () => setCustomPhase(isTruthy(ACOPH_CUSTOM_RAW)),
    [ACOPH_CUSTOM_RAW]
  );
  useEffect(
    () => setCustddm(isTruthy(ACOPH_CUSTDDM_RAW)),
    [ACOPH_CUSTDDM_RAW]
  );
  useEffect(
    () => setPeaksearch(isTruthy(ACOMPS_PEAKSEARCH_RAW)),
    [ACOMPS_PEAKSEARCH_RAW]
  );

  // onChange handlers
  const handleExtend = useCallback(
    async (item: any) => setExtend(isTruthy(item._value)),
    []
  );
  const handleBgrid = useCallback(
    async (item: any) => setBgrid(isTruthy(item._value)),
    []
  );
  const handleBseed = useCallback(
    async (item: any) => setBseed(isTruthy(item._value)),
    []
  );
  const handleBresol = useCallback(
    async (item: any) => setBresol(isTruthy(item._value)),
    []
  );
  const handleBexclude = useCallback(
    async (item: any) => setBexclude(isTruthy(item._value)),
    []
  );
  const handleBecut = useCallback(
    async (item: any) => setBecut(isTruthy(item._value)),
    []
  );
  const handleCustomPhase = useCallback(
    async (item: any) => setCustomPhase(isTruthy(item._value)),
    []
  );
  const handleCustddm = useCallback(
    async (item: any) => setCustddm(isTruthy(item._value)),
    []
  );
  const handlePeaksearch = useCallback(
    async (item: any) => setPeaksearch(isTruthy(item._value)),
    []
  );

  const numTrials =
    typeof ACOPH_TRIALS === "number"
      ? ACOPH_TRIALS
      : parseInt(ACOPH_TRIALS as string) || 1;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input Data ===== */}
        <CCP4i2Tab label="Input Data">
          {/* Phase input type radio */}
          <Box sx={{ display: "flex", alignItems: "center", gap: 1, mb: 1 }}>
            <Typography variant="body1">Run ACORN with</Typography>
            <CCP4i2TaskElement
              itemName="ACORN_PHSIN_TYPE"
              {...props}
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
          </Box>

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
              onChange={handleExtend}
            />
            {extend && (
              <Box
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                  flexWrap: "wrap",
                  pl: 3,
                }}
              >
                <Typography variant="body1">Resolution limit:</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="ACORN_EXTENDRES"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
                <Typography variant="body1">Å</Typography>
              </Box>
            )}
          </CCP4i2ContainerElement>

          {/* Refinement parameters */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Refinement parameters" }}
            containerHint="FolderLevel"
          >
            {/* Grid */}
            <CCP4i2TaskElement
              itemName="ACORN_BGRID"
              {...props}
              qualifiers={{ guiLabel: "Set grid spacing" }}
              onChange={handleBgrid}
            />
            {bgrid && (
              <Box
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                  flexWrap: "wrap",
                  pl: 3,
                }}
              >
                <Typography variant="body1">Grid spacing:</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="ACOGEN_GRID"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              </Box>
            )}

            {/* Seed */}
            <CCP4i2TaskElement
              itemName="ACORN_BSEED"
              {...props}
              qualifiers={{ guiLabel: "Set random seed" }}
              onChange={handleBseed}
            />
            {bseed && (
              <Box
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                  flexWrap: "wrap",
                  pl: 3,
                }}
              >
                <Typography variant="body1">Seed value:</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="ACOGEN_SEED"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              </Box>
            )}

            {/* Resolution limits */}
            <CCP4i2TaskElement
              itemName="ACORN_BRESOL"
              {...props}
              qualifiers={{ guiLabel: "Set resolution limits" }}
              onChange={handleBresol}
            />
            {bresol && (
              <Box
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                  flexWrap: "wrap",
                  pl: 3,
                }}
              >
                <Typography variant="body1">Low:</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="ACOREF_RESOLL"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
                <Typography variant="body1">High:</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="ACOREF_RESOLU"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
                <Typography variant="body1">Å</Typography>
              </Box>
            )}

            {/* Exclude reflections */}
            <CCP4i2TaskElement
              itemName="ACORN_BEXCLUDE"
              {...props}
              qualifiers={{ guiLabel: "Exclude reflections below E value" }}
              onChange={handleBexclude}
            />
            {bexclude && (
              <Box
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                  flexWrap: "wrap",
                  pl: 3,
                }}
              >
                <Typography variant="body1">Exclude below E =</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="ACOREF_EXCLUDE"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              </Box>
            )}

            {/* E cutoff */}
            <CCP4i2TaskElement
              itemName="ACORN_BECUT"
              {...props}
              qualifiers={{ guiLabel: "Set E value cutoff for supergrid" }}
              onChange={handleBecut}
            />
            {becut && (
              <Box
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                  flexWrap: "wrap",
                  pl: 3,
                }}
              >
                <Typography variant="body1">E cutoff:</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="ACOREF_ECUT"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              </Box>
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
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">DDM cutoff:</Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="ACOPH_CUTDDM"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">PS finish threshold:</Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="ACOPH_PSFINISH"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">CC finish threshold:</Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="ACOPH_CCFINISH"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>

            {/* Custom phase refinement protocol */}
            <CCP4i2TaskElement
              itemName="ACOPH_CUSTOM"
              {...props}
              qualifiers={{ guiLabel: "Use custom phase refinement protocol" }}
              onChange={handleCustomPhase}
            />
            {customPhase && (
              <Box sx={{ pl: 3, display: "flex", flexDirection: "column", gap: 1 }}>
                <CCP4i2TaskElement
                  itemName="ACOPH_CUSTDDM"
                  {...props}
                  qualifiers={{ guiLabel: "Use custom DDM parameters" }}
                  onChange={handleCustddm}
                />
                {custddm && (
                  <Box
                    sx={{
                      pl: 3,
                      display: "flex",
                      flexDirection: "column",
                      gap: 1,
                    }}
                  >
                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        gap: 1,
                      }}
                    >
                      <Typography variant="body1">
                        Number of trials:
                      </Typography>
                      <Box sx={{ width: "8rem" }}>
                        <CCP4i2TaskElement
                          itemName="ACOPH_TRIALS"
                          {...props}
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                    </Box>
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
              onChange={handlePeaksearch}
            />
            {peaksearch && (
              <Box
                sx={{
                  pl: 3,
                  display: "flex",
                  flexDirection: "column",
                  gap: 1,
                }}
              >
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">Maximum peaks:</Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="ACOMPS_MAXPEAKS"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">RMS multiplier:</Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="ACOMPS_RMSMULT"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </Box>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
