import React, { useCallback, useMemo } from "react";
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { InlineField } from "../task-elements/inline-field";

// ---------------------------------------------------------------------------
// Pipeline step definitions (mirrors crank2_basepipe.py)
// ---------------------------------------------------------------------------
const CRANK2_STEPS = [
  "substrdet",
  "phas",
  "handdet",
  "dmfull",
  "building",
  "ref",
];
const SHELX_STEPS = ["substrdet", "phdmmb", "building", "ref"];
const REBUILD_STEPS = [
  "refatompick",
  "handdet",
  "dmfull",
  "building",
  "ref",
];

function getBaseSteps(
  inputPartial: boolean,
  shelxcde: boolean,
  exptype: string | undefined
): string[] {
  let steps: string[];
  if (inputPartial) {
    steps = [...REBUILD_STEPS];
  } else if (shelxcde) {
    steps = [...SHELX_STEPS];
  } else {
    steps = [...CRANK2_STEPS];
  }
  // For SAD, replace 'phas' with 'refatompick'
  if (exptype === "SAD" && steps.includes("phas")) {
    steps[steps.indexOf("phas")] = "refatompick";
  }
  return steps;
}

function checkStartEnd(
  step: string,
  baseSteps: string[],
  startPipeline: string | undefined,
  endPipeline: string | undefined
): boolean {
  const idx = baseSteps.indexOf(step);
  if (idx < 0) return false;
  const startIdx = startPipeline ? baseSteps.indexOf(startPipeline) : 0;
  const endIdx = endPipeline
    ? baseSteps.indexOf(endPipeline)
    : baseSteps.length - 1;
  if (startIdx < 0 || endIdx < 0) return false;
  return idx >= startIdx && idx <= endIdx;
}

// ---------------------------------------------------------------------------
// Anomalous scattering coefficient row
// ---------------------------------------------------------------------------
const ScatteringRow: React.FC<{
  props: CCP4i2TaskInterfaceProps;
  suffix: string;
  onWavelengthChange?: (item: any) => void;
}> = ({ props, suffix, onWavelengthChange }) => (
  <Box
    sx={{
      display: "flex",
      alignItems: "center",
      gap: 1,
      flexWrap: "wrap",
    }}
  >
    <Typography variant="body2">f&apos;:</Typography>
    <Box sx={{ width: "8rem" }}>
      <CCP4i2TaskElement
        {...props}
        itemName={`FPRIME${suffix}`}
        qualifiers={{ guiLabel: " " }}
      />
    </Box>
    <Typography variant="body2">f&apos;&apos;:</Typography>
    <Box sx={{ width: "8rem" }}>
      <CCP4i2TaskElement
        {...props}
        itemName={`FDPRIME${suffix}`}
        qualifiers={{ guiLabel: " " }}
      />
    </Box>
    <Typography variant="body2">wavelength:</Typography>
    <Box sx={{ width: "8rem" }}>
      <CCP4i2TaskElement
        {...props}
        itemName={`WAVELENGTH${suffix}`}
        qualifiers={{ guiLabel: " " }}
        onChange={onWavelengthChange}
      />
    </Box>
    <Box sx={{ width: "8rem" }}>
      <CCP4i2TaskElement
        {...props}
        itemName={`DNAME${suffix}`}
        qualifiers={{ guiLabel: " " }}
      />
    </Box>
  </Box>
);

// ---------------------------------------------------------------------------
// Main interface
// ---------------------------------------------------------------------------
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, callPluginMethod, fetchDigest } = useJob(job.id);

  // =========================================================================
  // useTaskItem hooks — all at top level (React rules)
  // =========================================================================

  // --- Input Data ---
  const { value: startPipeline } = useTaskItem("START_PIPELINE");
  const { value: endPipeline } = useTaskItem("END_PIPELINE");

  const { value: atomType } = useTaskItem("ATOM_TYPE");
  const { forceUpdate: forceUpdateNUM_SUBSTR } = useTaskItem("NUMBER_SUBSTRUCTURE");
  useTaskItem("SUBSTRDET_NUM_DSUL");
  useTaskItem("XYZIN");
  useTaskItem("XYZIN_SUB");
  useTaskItem("XYZIN_SUB_RES");
  useTaskItem("PARTIAL_AS_SUBSTR");
  useTaskItem("SEQIN");

  // Cell params (visible when NON_MTZ)
  useTaskItem("CELL_A");
  useTaskItem("CELL_B");
  useTaskItem("CELL_C");
  useTaskItem("CELL_D");
  useTaskItem("CELL_E");
  useTaskItem("CELL_F");
  useTaskItem("SPACEGROUP");

  // Anomalous data files
  const { item: F_SIGFanomItem } = useTaskItem("F_SIGFanom");
  useTaskItem("F_SIGFanom2");
  useTaskItem("F_SIGFanom3");
  useTaskItem("F_SIGFanom4");
  useTaskItem("F_SIGFanom_nonmtz");
  useTaskItem("F_SIGFanom2_nonmtz");
  useTaskItem("F_SIGFanom3_nonmtz");
  useTaskItem("F_SIGFanom4_nonmtz");

  // Wavelength / f' / f'' (4 datasets)
  const { value: wavelength, forceUpdate: forceUpdateWAVELENGTH } =
    useTaskItem("WAVELENGTH");
  const { forceUpdate: forceUpdateFPRIME } = useTaskItem("FPRIME");
  const { forceUpdate: forceUpdateFDPRIME } = useTaskItem("FDPRIME");
  useTaskItem("WAVELENGTH2");
  useTaskItem("WAVELENGTH3");
  useTaskItem("WAVELENGTH4");
  useTaskItem("FPRIME2");
  useTaskItem("FDPRIME2");
  useTaskItem("FPRIME3");
  useTaskItem("FDPRIME3");
  useTaskItem("FPRIME4");
  useTaskItem("FDPRIME4");
  useTaskItem("DNAME");
  useTaskItem("DNAME2");
  useTaskItem("DNAME3");
  useTaskItem("DNAME4");

  // Native
  useTaskItem("F_SIGFnative");
  useTaskItem("F_SIGFnative_nonmtz");
  useTaskItem("SUBSTR_ATOMS_NATIVE");

  // Free
  const { value: freeVal } = useTaskItem("FREE");
  useTaskItem("FREERFLAG");
  useTaskItem("FREE_RATIO");

  // Phases
  useTaskItem("FPHIN_HL");

  // --- Important Options ---
  useTaskItem("RESIDUES_MON");
  const { value: monomersAsym } = useTaskItem("MONOMERS_ASYM");
  useTaskItem("SOLVENT_CONTENT");
  const { value: exptype } = useTaskItem("EXPTYPE");

  // --- Advanced Options ---
  useTaskItem("FAEST_PROGRAM");
  useTaskItem("SUBSTRDET_PROGRAM");
  useTaskItem("SUBSTRDET_HIGH_RES_CUTOFF");
  useTaskItem("SUBSTRDET_NUM_TRIALS");
  useTaskItem("SUBSTRDET_THRESHOLD_STOP");
  useTaskItem("SUBSTRDET_THRESHOLD_WEAK");
  useTaskItem("SUBSTRDET_MIN_DIST_ATOMS");
  useTaskItem("SUBSTRDET_MIN_DIST_SYMM_ATOMS");
  useTaskItem("SUBSTRDET_NUM_THREADS");
  useTaskItem("KEYWORDS_SUBSTRDET");

  useTaskItem("REFATOMPICK_NUM_ITER");
  useTaskItem("REFATOMPICK_REFCYC");
  useTaskItem("REFATOMPICK_RMS_THRESHOLD");
  useTaskItem("REFATOMPICK_OCC_CUT");

  useTaskItem("HANDDET_THRESHOLD_DISCRIM");
  useTaskItem("HANDDET_DMFULL_DM_PROGRAM");
  useTaskItem("HANDDET_DMFULL_PHCOMB_PROGRAM");

  useTaskItem("DMFULL_DM_PROGRAM");
  useTaskItem("DMFULL_PHCOMB_PROGRAM");
  useTaskItem("DMFULL_DMCYC");
  useTaskItem("DMFULL_THRESHOLD_STOP");
  useTaskItem("KEYWORDS_DMFULL_DM");

  useTaskItem("MB_PROGRAM");
  useTaskItem("COMB_PHDMMB_DMFULL_DM_PROGRAM");
  useTaskItem("COMB_PHDMMB_MINBIGCYC");
  useTaskItem("COMB_PHDMMB_MAXBIGCYC");
  useTaskItem("COMB_PHDMMB_NUM_PARALLEL");
  useTaskItem("COMB_PHDMMB_START_SHELXE");
  useTaskItem("KEYWORDS_COMB_SHELXE");
  useTaskItem("COMB_PHDMMB_EXCLUDE_FREE");
  useTaskItem("COMB_PHDMMB_NCS_DET");
  useTaskItem("COMB_PHDMMB_NCS_DET_MR");
  useTaskItem("COMB_PHDMMB_SKIP_INITIAL_BUILD");
  useTaskItem("COMB_PHDMMB_REBUILD_ONLY");
  useTaskItem("KEYWORDS_COMB_DM");
  useTaskItem("KEYWORDS_MB");

  useTaskItem("MBREF_BIGCYC");
  useTaskItem("MBREF_REF_PROGRAM");
  useTaskItem("MBREF_EXCLUDE_FREE");

  useTaskItem("PHDMMB_DMCYC");
  useTaskItem("PHDMMB_BIGCYC");
  useTaskItem("PHDMMB_THRESHOLD_STOP");
  useTaskItem("PHDMMB_THRESHOLD_HAND_STOP");
  useTaskItem("PHDMMB_THOROUGH_BUILD");
  useTaskItem("ARGUMENTS_SHELXE");

  useTaskItem("PHAS_PROGRAM");
  useTaskItem("PHAS_CYCLES");

  useTaskItem("REF_PROGRAM");
  useTaskItem("REF_CYCLES");
  useTaskItem("REF_EXCLUDE_FREE");
  useTaskItem("KEYWORDS_REF");

  useTaskItem("CLEANUP");
  useTaskItem("RESIDUES_MON_COPY");

  // =========================================================================
  // Boolean toggles (CBoolean pattern via useBoolToggle)
  // =========================================================================
  const inputPartial = useBoolToggle(useTaskItem, "INPUT_PARTIAL");
  const inputSequence = useBoolToggle(useTaskItem, "INPUT_SEQUENCE");
  const inputPhases = useBoolToggle(useTaskItem, "INPUT_PHASES");
  const nonMtz = useBoolToggle(useTaskItem, "NON_MTZ");
  const mad2 = useBoolToggle(useTaskItem, "MAD2");
  const mad3 = useBoolToggle(useTaskItem, "MAD3");
  const mad4 = useBoolToggle(useTaskItem, "MAD4");
  const native = useBoolToggle(useTaskItem, "NATIVE");
  const shelxcde = useBoolToggle(useTaskItem, "SHELXCDE");
  const doHanddet = useBoolToggle(useTaskItem, "DO_HANDDET");
  const useComb = useBoolToggle(useTaskItem, "USE_COMB");

  // =========================================================================
  // Pipeline visibility (mirrors crank2_basepipe.py)
  // =========================================================================
  const baseSteps = useMemo(
    () =>
      getBaseSteps(inputPartial.value, shelxcde.value, exptype as string | undefined),
    [inputPartial.value, shelxcde.value, exptype]
  );

  const check = useCallback(
    (step: string) =>
      checkStartEnd(
        step,
        baseSteps,
        startPipeline as string | undefined,
        endPipeline as string | undefined
      ),
    [baseSteps, startPipeline, endPipeline]
  );

  const showDetection = useMemo(() => check("substrdet"), [check]);
  const showPeakSearch = useMemo(() => check("refatompick"), [check]);
  const showShelxCDE = useMemo(() => check("phdmmb"), [check]);
  const showPhasing = useMemo(() => {
    if (inputPartial.value) return false;
    return check("phas") && !shelxcde.value;
  }, [check, inputPartial.value, shelxcde.value]);
  const showHandDet = useMemo(() => {
    if (inputPartial.value) return false;
    return check("handdet") && !shelxcde.value;
  }, [check, inputPartial.value, shelxcde.value]);
  const showDensityMod = useMemo(() => {
    if (inputPartial.value) return false;
    return check("dmfull") && !shelxcde.value;
  }, [check, inputPartial.value, shelxcde.value]);
  const showModelBuilding = useMemo(() => check("building"), [check]);
  const showRefine = useMemo(() => check("ref"), [check]);

  // =========================================================================
  // onChange handlers
  // =========================================================================

  /** Compute f'/f'' for dataset 1 */
  const computeScatteringFactors = useCallback(
    async (at?: string, wl?: number) => {
      const a = at ?? (atomType as string);
      const w = wl ?? (wavelength as number);
      if (!a || !w || !job || job.status !== 1) return;
      const result = await callPluginMethod("compute_anomalous_scattering", {
        atom_type: String(a),
        wavelength: Number(w),
      });
      if (result && result.fp !== undefined && result.fpp !== undefined) {
        await forceUpdateFPRIME(result.fp);
        await forceUpdateFDPRIME(result.fpp);
      }
    },
    [
      atomType,
      wavelength,
      callPluginMethod,
      forceUpdateFPRIME,
      forceUpdateFDPRIME,
      job?.id,
      job?.status,
    ]
  );

  /** Estimate NUMBER_SUBSTRUCTURE (and NUM_DSUL for S) from sequence + atom type */
  const estimateSubstructure = useCallback(
    async (at?: string) => {
      const a = at ?? (atomType as string);
      if (!a || !job || job.status !== 1) return;
      const result = await callPluginMethod("estimate_num_substructure", {
        atom_type: String(a),
      });
      if (result && result.estimate !== undefined) {
        await forceUpdateNUM_SUBSTR(result.estimate);
      }
    },
    [atomType, callPluginMethod, forceUpdateNUM_SUBSTR, job?.id, job?.status]
  );

  const handleAtomTypeChange = useCallback(
    (item: any) => {
      const at = item?._value ? String(item._value) : undefined;
      computeScatteringFactors(at);
      estimateSubstructure(at);
    },
    [computeScatteringFactors, estimateSubstructure]
  );

  const handleSeqinChange = useCallback(
    () => estimateSubstructure(),
    [estimateSubstructure]
  );

  const handleMonomersAsymChange = useCallback(
    () => estimateSubstructure(),
    [estimateSubstructure]
  );

  const handleWavelengthChange = useCallback(
    (item: any) => {
      const wl = item?._value ? Number(item._value) : undefined;
      computeScatteringFactors(undefined, wl);
    },
    [computeScatteringFactors]
  );

  const handleF_SIGFanomChange = useCallback(async () => {
    if (!F_SIGFanomItem?._objectPath || !job || job.status !== 1) return;
    const digest = await fetchDigest(F_SIGFanomItem._objectPath);
    if (digest?.wavelengths?.length > 0) {
      const wl = digest.wavelengths[digest.wavelengths.length - 1];
      if (wl && wl > 0 && wl < 9) {
        await forceUpdateWAVELENGTH(wl);
        await computeScatteringFactors(undefined, wl);
      }
    }
  }, [
    F_SIGFanomItem?._objectPath,
    fetchDigest,
    forceUpdateWAVELENGTH,
    computeScatteringFactors,
    job?.id,
    job?.status,
  ]);

  // =========================================================================
  // Render
  // =========================================================================
  return (
    <Paper>
      <CCP4i2Tabs>
        {/* ================================================================
            TAB 1: INPUT DATA
            ================================================================ */}
        <CCP4i2Tab label="Input Data" key="input">
          {/* Pipeline selection */}
          <FieldRow>
            <InlineField
              label="Start pipeline with"
              width="14rem"
              after={
                <InlineField label="and end with" width="14rem">
                  <CCP4i2TaskElement
                    {...props}
                    itemName="END_PIPELINE"
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              }
            >
              <CCP4i2TaskElement
                {...props}
                itemName="START_PIPELINE"
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </FieldRow>

          {/* Input partial model */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Input partial model (MR-SAD, model rebuilding)",
            }}
            containerHint="FolderLevel"
            initiallyOpen={false}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="INPUT_PARTIAL"
              qualifiers={{
                guiLabel:
                  "If a partial protein model is available from molecular replacement",
              }}
              onChange={inputPartial.onChange}
            />
            {inputPartial.value && (
              <>
                <CCP4i2TaskElement {...props} itemName="XYZIN" />
                <CCP4i2TaskElement
                  {...props}
                  itemName="PARTIAL_AS_SUBSTR"
                  qualifiers={{
                    guiLabel:
                      "Start from anomalous substructure only",
                  }}
                />
              </>
            )}
          </CCP4i2ContainerElement>

          {/* Input protein sequence */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input protein sequence" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="INPUT_SEQUENCE"
              onChange={inputSequence.onChange}
            />
            {inputSequence.value && (
              <CCP4i2TaskElement {...props} itemName="SEQIN" onChange={handleSeqinChange} />
            )}
            {!inputSequence.value && (
              <CCP4i2TaskElement
                {...props}
                itemName="RESIDUES_MON_COPY"
                qualifiers={{
                  guiLabel: "Number of residues per monomer",
                }}
              />
            )}
          </CCP4i2ContainerElement>

          {/* Input starting phases */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input starting phases" }}
            containerHint="FolderLevel"
            initiallyOpen={false}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="INPUT_PHASES"
              onChange={inputPhases.onChange}
            />
            {inputPhases.value && (
              <CCP4i2TaskElement {...props} itemName="FPHIN_HL" />
            )}
          </CCP4i2ContainerElement>

          {/* Crystal #1 composition and anomalous datasets */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Crystal #1 composition and collected anomalous dataset(s)",
            }}
            containerHint="FolderLevel"
          >
            {/* Atom type + number of atoms */}
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="ATOM_TYPE"
                qualifiers={{ guiLabel: "Substructure atom" }}
                onChange={handleAtomTypeChange}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="NUMBER_SUBSTRUCTURE"
                qualifiers={{
                  guiLabel:
                    "Number of substr. atoms in asymmetric unit",
                }}
              />
            </FieldRow>

            {/* Disulfide count — visible for Sulphur */}
            <CCP4i2TaskElement
              {...props}
              itemName="SUBSTRDET_NUM_DSUL"
              qualifiers={{
                guiLabel:
                  "Number of S-S pairs searched for as 1 supersulfur",
              }}
              visibility={() =>
                String(atomType).toUpperCase() === "S"
              }
            />

            {/* Substructure model (visible when detection not first step) */}
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_SUB"
              qualifiers={{ guiLabel: "Substructure" }}
              visibility={() => !showDetection}
            />

            {/* Cell params when NON_MTZ */}
            {nonMtz.value && (
              <FieldRow equalWidth={false} size="xs">
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_A"
                  qualifiers={{ guiLabel: "a" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_B"
                  qualifiers={{ guiLabel: "b" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_C"
                  qualifiers={{ guiLabel: "c" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_D"
                  qualifiers={{ guiLabel: "\u03B1" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_E"
                  qualifiers={{ guiLabel: "\u03B2" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_F"
                  qualifiers={{ guiLabel: "\u03B3" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="SPACEGROUP"
                  qualifiers={{ guiLabel: "Spacegroup" }}
                />
              </FieldRow>
            )}

            {/* Anomalous data (Friedel pairs) — Dataset 1 */}
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Anomalous data (Friedel pairs)" }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                {...props}
                itemName="NON_MTZ"
                qualifiers={{
                  guiLabel:
                    "Input unmerged/merged SCA/XDS/SHELX format",
                }}
                onChange={nonMtz.onChange}
              />
              {nonMtz.value ? (
                <CCP4i2TaskElement
                  {...props}
                  itemName="F_SIGFanom_nonmtz"
                  qualifiers={{ guiLabel: "Reflections" }}
                />
              ) : (
                <CCP4i2TaskElement
                  {...props}
                  itemName="F_SIGFanom"
                  qualifiers={{ guiLabel: "Reflections" }}
                  onChange={handleF_SIGFanomChange}
                />
              )}
              <ScatteringRow
                props={props}
                suffix=""
                onWavelengthChange={handleWavelengthChange}
              />
            </CCP4i2ContainerElement>

            {/* MAD Dataset 2 */}
            {!inputPartial.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{
                  guiLabel: "Input anomalous data #2 (MAD)",
                }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement
                  {...props}
                  itemName="MAD2"
                  onChange={mad2.onChange}
                />
                {mad2.value && (
                  <>
                    {nonMtz.value ? (
                      <CCP4i2TaskElement
                        {...props}
                        itemName="F_SIGFanom2_nonmtz"
                        qualifiers={{ guiLabel: "Reflections" }}
                      />
                    ) : (
                      <CCP4i2TaskElement
                        {...props}
                        itemName="F_SIGFanom2"
                        qualifiers={{ guiLabel: "Reflections" }}
                      />
                    )}
                    <ScatteringRow props={props} suffix="2" />
                  </>
                )}
              </CCP4i2ContainerElement>
            )}

            {/* MAD Dataset 3 */}
            {!inputPartial.value && mad2.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{
                  guiLabel: "Input anomalous data #3 (MAD)",
                }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement
                  {...props}
                  itemName="MAD3"
                  onChange={mad3.onChange}
                />
                {mad3.value && (
                  <>
                    {nonMtz.value ? (
                      <CCP4i2TaskElement
                        {...props}
                        itemName="F_SIGFanom3_nonmtz"
                        qualifiers={{ guiLabel: "Reflections" }}
                      />
                    ) : (
                      <CCP4i2TaskElement
                        {...props}
                        itemName="F_SIGFanom3"
                        qualifiers={{ guiLabel: "Reflections" }}
                      />
                    )}
                    <ScatteringRow props={props} suffix="3" />
                  </>
                )}
              </CCP4i2ContainerElement>
            )}

            {/* MAD Dataset 4 */}
            {!inputPartial.value && mad3.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{
                  guiLabel: "Input anomalous data #4 (MAD)",
                }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement
                  {...props}
                  itemName="MAD4"
                  onChange={mad4.onChange}
                />
                {mad4.value && (
                  <>
                    {nonMtz.value ? (
                      <CCP4i2TaskElement
                        {...props}
                        itemName="F_SIGFanom4_nonmtz"
                        qualifiers={{ guiLabel: "Reflections" }}
                      />
                    ) : (
                      <CCP4i2TaskElement
                        {...props}
                        itemName="F_SIGFanom4"
                        qualifiers={{ guiLabel: "Reflections" }}
                      />
                    )}
                    <ScatteringRow props={props} suffix="4" />
                  </>
                )}
              </CCP4i2ContainerElement>
            )}
          </CCP4i2ContainerElement>

          {/* Native observations */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Input native observations (Crystal #2)",
            }}
            containerHint="FolderLevel"
            initiallyOpen={false}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="NATIVE"
              onChange={native.onChange}
            />
            {native.value && (
              <>
                {nonMtz.value ? (
                  <CCP4i2TaskElement
                    {...props}
                    itemName="F_SIGFnative_nonmtz"
                    qualifiers={{ guiLabel: "Reflections" }}
                  />
                ) : (
                  <CCP4i2TaskElement
                    {...props}
                    itemName="F_SIGFnative"
                    qualifiers={{ guiLabel: "Reflections" }}
                  />
                )}
                <CCP4i2TaskElement
                  {...props}
                  itemName="SUBSTR_ATOMS_NATIVE"
                  qualifiers={{
                    guiLabel:
                      "Substructure atoms present in the native crystal",
                  }}
                />
              </>
            )}
          </CCP4i2ContainerElement>

          {/* Exclude free */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Exclude free" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="FREE" />
            {freeVal === "existing" && (
              <CCP4i2TaskElement {...props} itemName="FREERFLAG" />
            )}
            {freeVal === "new" && (
              <InlineField label="consisting of" width="6rem" hint="% of reflections">
                <CCP4i2TaskElement
                  {...props}
                  itemName="FREE_RATIO"
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ================================================================
            TAB 2: IMPORTANT OPTIONS
            ================================================================ */}
        <CCP4i2Tab label="Important Options" key="important">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Crystal composition" }}
            containerHint="FolderLevel"
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="RESIDUES_MON"
                qualifiers={{ guiLabel: "Residues/monomer" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="MONOMERS_ASYM"
                qualifiers={{ guiLabel: "NCS copies" }}
                onChange={handleMonomersAsymChange}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_CONTENT"
                qualifiers={{ guiLabel: "Solvent Content" }}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Experiment type" }}
            containerHint="FolderLevel"
            visibility={() => native.value || mad2.value}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="EXPTYPE"
              qualifiers={{ guiLabel: "Phasing method" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ================================================================
            TAB 3: ADVANCED OPTIONS
            ================================================================ */}
        <CCP4i2Tab label="Advanced Options" key="advanced">
          <Typography
            variant="body2"
            sx={{ fontStyle: "italic", mb: 1, pl: 1 }}
          >
            Note: empty input fields mean that internal program defaults
            will be used.
          </Typography>

          {/* SHELXC/D/E toggle */}
          <CCP4i2TaskElement
            {...props}
            itemName="SHELXCDE"
            qualifiers={{ guiLabel: "Use SHELXC/D/E" }}
            onChange={shelxcde.onChange}
            visibility={() => !inputPartial.value}
          />

          {/* Substructure detection */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Substructure detection" }}
            containerHint="FolderLevel"
            visibility={() => showDetection}
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="FAEST_PROGRAM"
                qualifiers={{ guiLabel: "FA Estimation program" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="SUBSTRDET_PROGRAM"
                qualifiers={{ guiLabel: "Detection program" }}
              />
            </FieldRow>
            <CCP4i2TaskElement
              {...props}
              itemName="SUBSTRDET_HIGH_RES_CUTOFF"
              qualifiers={{ guiLabel: "High resolution cutoff" }}
            />
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="SUBSTRDET_NUM_TRIALS"
                qualifiers={{ guiLabel: "Num. trials" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="SUBSTRDET_THRESHOLD_STOP"
                qualifiers={{ guiLabel: "CFOM threshold" }}
              />
            </FieldRow>
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="SUBSTRDET_MIN_DIST_ATOMS"
                qualifiers={{
                  guiLabel: "Minimum distance between atoms",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="SUBSTRDET_MIN_DIST_SYMM_ATOMS"
                qualifiers={{
                  guiLabel: "Atoms in special positions allowed",
                }}
              />
            </FieldRow>
            <CCP4i2TaskElement
              {...props}
              itemName="SUBSTRDET_NUM_THREADS"
              qualifiers={{ guiLabel: "Number of CPU threads" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="KEYWORDS_SUBSTRDET"
              qualifiers={{
                guiLabel: "Custom program keywords (comma separated)",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Substructure improvement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Substructure improvement" }}
            containerHint="FolderLevel"
            visibility={() => showPeakSearch}
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="REFATOMPICK_NUM_ITER"
                qualifiers={{
                  guiLabel: "Max. num. of iterations",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="REFATOMPICK_REFCYC"
                qualifiers={{
                  guiLabel: "Number of ref. cycles per iteration",
                }}
              />
            </FieldRow>
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="REFATOMPICK_RMS_THRESHOLD"
                qualifiers={{
                  guiLabel:
                    "Pick new atoms from anom. maps at peaks above RMS",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="REFATOMPICK_OCC_CUT"
                qualifiers={{
                  guiLabel: "Remove atoms with occupancy below",
                }}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          {/* Substructure phasing */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Substructure phasing" }}
            containerHint="FolderLevel"
            visibility={() => showPhasing}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="PHAS_PROGRAM"
              qualifiers={{ guiLabel: "Phasing program" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="PHAS_CYCLES"
              qualifiers={{ guiLabel: "Refinement cycles" }}
            />
          </CCP4i2ContainerElement>

          {/* Poly-Ala tracing (SHELXE) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Den.mod. & poly-Ala tracing (SHELXE)",
            }}
            containerHint="FolderLevel"
            visibility={() => showShelxCDE}
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="PHDMMB_DMCYC"
                qualifiers={{ guiLabel: "DM cycles" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="PHDMMB_BIGCYC"
                qualifiers={{ guiLabel: "Model building cycles" }}
              />
            </FieldRow>
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="PHDMMB_THRESHOLD_STOP"
                qualifiers={{ guiLabel: "CC threshold" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="PHDMMB_THRESHOLD_HAND_STOP"
                qualifiers={{
                  guiLabel: "Hand determination CC threshold",
                }}
              />
            </FieldRow>
            <CCP4i2TaskElement
              {...props}
              itemName="PHDMMB_THOROUGH_BUILD"
              qualifiers={{ guiLabel: "Use thorough SHELXE building" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="ARGUMENTS_SHELXE"
              qualifiers={{ guiLabel: "SHELXE arguments" }}
            />
          </CCP4i2ContainerElement>

          {/* Hand determination */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Hand determination" }}
            containerHint="FolderLevel"
            visibility={() => showHandDet}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="DO_HANDDET"
              onChange={doHanddet.onChange}
            />
            {doHanddet.value && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="HANDDET_THRESHOLD_DISCRIM"
                  qualifiers={{
                    guiLabel:
                      "Requested hand determination discrimination",
                  }}
                />
                <FieldRow>
                  <CCP4i2TaskElement
                    {...props}
                    itemName="HANDDET_DMFULL_DM_PROGRAM"
                    qualifiers={{
                      guiLabel: "Density modif. program",
                    }}
                  />
                  <CCP4i2TaskElement
                    {...props}
                    itemName="HANDDET_DMFULL_PHCOMB_PROGRAM"
                    qualifiers={{
                      guiLabel: "Phase combination program",
                    }}
                  />
                </FieldRow>
              </>
            )}
          </CCP4i2ContainerElement>

          {/* Density modification */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Density modification" }}
            containerHint="FolderLevel"
            visibility={() => showDensityMod}
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="DMFULL_DM_PROGRAM"
                qualifiers={{ guiLabel: "Density modif. program" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="DMFULL_PHCOMB_PROGRAM"
                qualifiers={{
                  guiLabel: "Phase combination program",
                }}
              />
            </FieldRow>
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="DMFULL_DMCYC"
                qualifiers={{ guiLabel: "Number of iterations" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="DMFULL_THRESHOLD_STOP"
                qualifiers={{ guiLabel: "FOM threshold" }}
              />
            </FieldRow>
            <CCP4i2TaskElement
              {...props}
              itemName="KEYWORDS_DMFULL_DM"
              qualifiers={{
                guiLabel: "Custom options for dens.mod. program",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Model building */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model building" }}
            containerHint="FolderLevel"
            visibility={() => showModelBuilding}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="USE_COMB"
              qualifiers={{
                guiLabel:
                  "Combine phase, model and density modif. information",
              }}
              onChange={useComb.onChange}
            />

            {/* Combined model building */}
            {useComb.value && (
              <>
                <FieldRow>
                  <CCP4i2TaskElement
                    {...props}
                    itemName="MB_PROGRAM"
                    qualifiers={{
                      guiLabel: "Model building program",
                    }}
                  />
                  <CCP4i2TaskElement
                    {...props}
                    itemName="COMB_PHDMMB_DMFULL_DM_PROGRAM"
                    qualifiers={{
                      guiLabel: "Density modif. program",
                    }}
                  />
                </FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="KEYWORDS_MB"
                  qualifiers={{
                    guiLabel: "Custom keywords for building program",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="KEYWORDS_COMB_DM"
                  qualifiers={{
                    guiLabel:
                      "Custom keywords for dens.mod. program",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="COMB_PHDMMB_START_SHELXE"
                  qualifiers={{
                    guiLabel:
                      "Start with a few SHELXE tracing cycles",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="KEYWORDS_COMB_SHELXE"
                  qualifiers={{
                    guiLabel: "with custom keywords",
                  }}
                />
              </>
            )}

            {/* Non-combined (mbref) */}
            {!useComb.value && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="MB_PROGRAM"
                  qualifiers={{
                    guiLabel: "Model building program",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="MBREF_REF_PROGRAM"
                  qualifiers={{
                    guiLabel: "Refinement program",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="MBREF_BIGCYC"
                  qualifiers={{
                    guiLabel: "Building cycles",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="MBREF_EXCLUDE_FREE"
                  qualifiers={{
                    guiLabel: "Exclude free reflections",
                  }}
                />
              </>
            )}

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="COMB_PHDMMB_MINBIGCYC"
                qualifiers={{
                  guiLabel: "Minimum number of cycles",
                }}
                visibility={() => useComb.value}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="COMB_PHDMMB_MAXBIGCYC"
                qualifiers={{
                  guiLabel: "Maximum number of cycles",
                }}
                visibility={() => useComb.value}
              />
            </FieldRow>

            <CCP4i2TaskElement
              {...props}
              itemName="COMB_PHDMMB_NUM_PARALLEL"
              qualifiers={{
                guiLabel:
                  "Parallel model building: simultaneous building and refinement processes",
              }}
              visibility={() => useComb.value}
            />

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="COMB_PHDMMB_SKIP_INITIAL_BUILD"
                qualifiers={{
                  guiLabel: "Skip the first model building cycle",
                }}
                visibility={() => useComb.value}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="COMB_PHDMMB_REBUILD_ONLY"
                qualifiers={{ guiLabel: "Soft rebuilding" }}
                visibility={() => useComb.value}
              />
            </FieldRow>

            <CCP4i2TaskElement
              {...props}
              itemName="COMB_PHDMMB_EXCLUDE_FREE"
              qualifiers={{
                guiLabel:
                  "Exclude the free reflections in model building",
              }}
              visibility={() => useComb.value}
            />
          </CCP4i2ContainerElement>

          {/* Final refinement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Final refinement" }}
            containerHint="FolderLevel"
            visibility={() => showRefine}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="REF_PROGRAM"
              qualifiers={{ guiLabel: "Refinement program" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="REF_CYCLES"
              qualifiers={{ guiLabel: "Refinement cycles" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="REF_EXCLUDE_FREE"
              qualifiers={{
                guiLabel: "Exclude free reflections",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Cleanup */}
          <CCP4i2TaskElement
            {...props}
            itemName="CLEANUP"
            qualifiers={{
              guiLabel: "Remove all intermediate mtz files at the end?",
            }}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
