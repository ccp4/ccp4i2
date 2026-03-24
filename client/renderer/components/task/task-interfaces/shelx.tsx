import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";

const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

// ---------------------------------------------------------------------------
// Pipeline step definitions for SHELX (always SHELXC/D/E mode)
// ---------------------------------------------------------------------------
const SHELX_STEPS = ["substrdet", "phdmmb", "building", "ref"];

function checkStartEnd(
  step: string,
  startPipeline: string | undefined,
  endPipeline: string | undefined
): boolean {
  const idx = SHELX_STEPS.indexOf(step);
  if (idx < 0) return false;
  const startIdx = startPipeline ? SHELX_STEPS.indexOf(startPipeline) : 0;
  const endIdx = endPipeline
    ? SHELX_STEPS.indexOf(endPipeline)
    : SHELX_STEPS.length - 1;
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
      pl: 2,
    }}
  >
    <Typography variant="body2">f&apos;:</Typography>
    <Box sx={{ width: "6rem" }}>
      <CCP4i2TaskElement
        {...props}
        itemName={`FPRIME${suffix}`}
        qualifiers={{ guiLabel: " " }}
      />
    </Box>
    <Typography variant="body2">f&quot;:</Typography>
    <Box sx={{ width: "6rem" }}>
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

/**
 * Task interface for SHELX — SHELXC/D/E phasing and building.
 *
 * SHELX is a simplified variant of the Crank2 pipeline that always uses
 * SHELXC/D/E mode. Compared to the full Crank2 interface it omits:
 * - Input partial model (MR-SAD)
 * - Input starting phases
 * - The SHELXC/D/E toggle (always on)
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, callPluginMethod, fetchDigest } =
    useJob(job.id);

  // =========================================================================
  // useTaskItem hooks — all at top level (React rules)
  // =========================================================================

  // --- Input Data ---
  const { value: startPipeline } = useTaskItem("START_PIPELINE");
  const { value: endPipeline } = useTaskItem("END_PIPELINE");

  const { value: inputSequenceRaw } = useTaskItem("INPUT_SEQUENCE");
  const { value: nonMtzRaw } = useTaskItem("NON_MTZ");
  const { value: mad2Raw } = useTaskItem("MAD2");
  const { value: mad3Raw } = useTaskItem("MAD3");
  const { value: mad4Raw } = useTaskItem("MAD4");
  const { value: nativeRaw } = useTaskItem("NATIVE");

  const { item: ATOM_TYPEItem, value: atomType } = useTaskItem("ATOM_TYPE");
  const { forceUpdate: forceUpdateNUM_SUBSTR } = useTaskItem("NUMBER_SUBSTRUCTURE");
  useTaskItem("SUBSTRDET_NUM_DSUL");
  useTaskItem("SEQIN");
  useTaskItem("RESIDUES_MON_COPY");

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

  // --- Important Options ---
  useTaskItem("RESIDUES_MON");
  const { value: monomersAsym } = useTaskItem("MONOMERS_ASYM");
  useTaskItem("SOLVENT_CONTENT");
  useTaskItem("EXPTYPE");

  // --- Advanced Options: Substructure detection ---
  useTaskItem("SUBSTRDET_HIGH_RES_CUTOFF");
  useTaskItem("SUBSTRDET_HIGH_RES_CUTOFF_CCHALF");
  useTaskItem("SUBSTRDET_NUM_TRIALS");
  useTaskItem("SUBSTRDET_THRESHOLD_STOP");
  useTaskItem("SUBSTRDET_MIN_DIST_ATOMS");
  useTaskItem("SUBSTRDET_MIN_DIST_SYMM_ATOMS");
  useTaskItem("SUBSTRDET_NUM_THREADS");
  useTaskItem("KEYWORDS_SUBSTRDET");

  // --- Advanced Options: SHELXE (phdmmb) ---
  useTaskItem("PHDMMB_DMCYC");
  useTaskItem("PHDMMB_BIGCYC");
  useTaskItem("PHDMMB_THRESHOLD_STOP");
  useTaskItem("PHDMMB_THRESHOLD_HAND_STOP");
  useTaskItem("SUBSTRDET_THRESHOLD_WEAK");
  useTaskItem("ARGUMENTS_SHELXE");

  // --- Advanced Options: Model building ---
  const { value: useCombRaw } = useTaskItem("USE_COMB");
  useTaskItem("MB_PROGRAM");
  useTaskItem("KEYWORDS_MB");
  useTaskItem("MBREF_BIGCYC");
  useTaskItem("MBREF_EXCLUDE_FREE");

  // --- Advanced Options: Final refinement ---
  useTaskItem("REF_PROGRAM");
  useTaskItem("REF_CYCLES");
  useTaskItem("REF_EXCLUDE_FREE");
  useTaskItem("KEYWORDS_REF");

  useTaskItem("CLEANUP");

  // --- Initialization: force SHELXCDE defaults ---
  const { updateNoMutate: updateSHELXCDE } = useTaskItem("SHELXCDE");
  const { updateNoMutate: updateUSE_COMB } = useTaskItem("USE_COMB");
  const { updateNoMutate: updateSHELX_SEPAR } = useTaskItem("SHELX_SEPAR");
  const { updateNoMutate: updateMB_PROGRAM, value: mbProgram } =
    useTaskItem("MB_PROGRAM");
  const { value: shelxcdeVal } = useTaskItem("SHELXCDE");
  const { value: useCombVal } = useTaskItem("USE_COMB");
  const { value: shelxSeparVal } = useTaskItem("SHELX_SEPAR");

  const initDone = useRef(false);

  useEffect(() => {
    initDone.current = false;
  }, [job?.id]);

  useEffect(() => {
    if (initDone.current || !job || job.status !== 1) return;
    const updates: Promise<any>[] = [];
    if (!isTruthy(shelxcdeVal)) updates.push(updateSHELXCDE(true));
    if (isTruthy(useCombVal)) updates.push(updateUSE_COMB(false));
    if (!isTruthy(shelxSeparVal)) updates.push(updateSHELX_SEPAR(true));
    if (mbProgram !== "buccaneer") updates.push(updateMB_PROGRAM("buccaneer"));
    if (updates.length > 0) Promise.all(updates).catch(console.error);
    initDone.current = true;
  }, [job?.status, shelxcdeVal, useCombVal, shelxSeparVal, mbProgram]);

  // =========================================================================
  // Local toggle state (CBoolean pattern — immediate UI)
  // =========================================================================
  const [inputSequence, setInputSequence] = useState(() =>
    isTruthy(inputSequenceRaw)
  );
  const [nonMtz, setNonMtz] = useState(() => isTruthy(nonMtzRaw));
  const [mad2, setMad2] = useState(() => isTruthy(mad2Raw));
  const [mad3, setMad3] = useState(() => isTruthy(mad3Raw));
  const [mad4, setMad4] = useState(() => isTruthy(mad4Raw));
  const [native, setNative] = useState(() => isTruthy(nativeRaw));
  const [useComb, setUseComb] = useState(() => isTruthy(useCombRaw));

  // Sync from server
  useEffect(() => setInputSequence(isTruthy(inputSequenceRaw)), [inputSequenceRaw]);
  useEffect(() => setNonMtz(isTruthy(nonMtzRaw)), [nonMtzRaw]);
  useEffect(() => setMad2(isTruthy(mad2Raw)), [mad2Raw]);
  useEffect(() => setMad3(isTruthy(mad3Raw)), [mad3Raw]);
  useEffect(() => setMad4(isTruthy(mad4Raw)), [mad4Raw]);
  useEffect(() => setNative(isTruthy(nativeRaw)), [nativeRaw]);
  useEffect(() => setUseComb(isTruthy(useCombRaw)), [useCombRaw]);

  const onToggle =
    (setter: (v: boolean) => void) => (item: any) =>
      setter(isTruthy(item._value));

  // =========================================================================
  // Pipeline visibility
  // =========================================================================
  const check = useCallback(
    (step: string) =>
      checkStartEnd(
        step,
        startPipeline as string | undefined,
        endPipeline as string | undefined
      ),
    [startPipeline, endPipeline]
  );

  const showDetection = useMemo(() => check("substrdet"), [check]);
  const showShelxCDE = useMemo(() => check("phdmmb"), [check]);
  const showModelBuilding = useMemo(() => check("building"), [check]);
  const showRefine = useMemo(() => check("ref"), [check]);

  // =========================================================================
  // onChange handlers — compute f'/f'' from atom type + wavelength
  // =========================================================================
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
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">Start pipeline with</Typography>
              <Box sx={{ width: "14rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="START_PIPELINE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">and end with</Typography>
              <Box sx={{ width: "14rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="END_PIPELINE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
          </FieldRow>

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
              onChange={onToggle(setInputSequence)}
            />
            {inputSequence && (
              <CCP4i2TaskElement {...props} itemName="SEQIN" onChange={handleSeqinChange} />
            )}
            {!inputSequence && (
              <CCP4i2TaskElement
                {...props}
                itemName="RESIDUES_MON_COPY"
                qualifiers={{
                  guiLabel: "Number of residues per monomer",
                }}
              />
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

            {/* Cell params when NON_MTZ */}
            {nonMtz && (
              <FieldRow equalWidth={false} size="xs">
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_A"
                  qualifiers={{ guiLabel: "Cell:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_B"
                  qualifiers={{ guiLabel: " " }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_C"
                  qualifiers={{ guiLabel: " " }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_D"
                  qualifiers={{ guiLabel: " " }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_E"
                  qualifiers={{ guiLabel: " " }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CELL_F"
                  qualifiers={{ guiLabel: " " }}
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
                onChange={onToggle(setNonMtz)}
              />
              {nonMtz ? (
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
                onChange={onToggle(setMad2)}
              />
              {mad2 && (
                <>
                  {nonMtz ? (
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

            {/* MAD Dataset 3 */}
            {mad2 && (
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
                  onChange={onToggle(setMad3)}
                />
                {mad3 && (
                  <>
                    {nonMtz ? (
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
            {mad3 && (
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
                  onChange={onToggle(setMad4)}
                />
                {mad4 && (
                  <>
                    {nonMtz ? (
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
              onChange={onToggle(setNative)}
            />
            {native && (
              <>
                {nonMtz ? (
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
              <Box
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                }}
              >
                <Typography variant="body1">
                  consisting of
                </Typography>
                <Box sx={{ width: "6rem" }}>
                  <CCP4i2TaskElement
                    {...props}
                    itemName="FREE_RATIO"
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
                <Typography variant="body1">
                  % of reflections
                </Typography>
              </Box>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ================================================================
            TAB 2: IMPORTANT OPTIONS
            ================================================================ */}
        <CCP4i2Tab label="Important Options" key="important">
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

          <CCP4i2TaskElement
            {...props}
            itemName="EXPTYPE"
            qualifiers={{ guiLabel: "Phasing method" }}
            visibility={() => native || mad2}
          />
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

          {/* Substructure detection */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Substructure detection" }}
            containerHint="FolderLevel"
            visibility={() => showDetection}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SUBSTRDET_HIGH_RES_CUTOFF"
              qualifiers={{ guiLabel: "High resolution cutoff" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="SUBSTRDET_HIGH_RES_CUTOFF_CCHALF"
              qualifiers={{
                guiLabel:
                  "Use CCanom1/2 based cutoff (if available)",
              }}
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

          {/* Density modification and poly-Ala tracing with SHELXE */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Density modification and poly-Ala tracing with SHELXE",
            }}
            containerHint="FolderLevel"
            visibility={() => showShelxCDE}
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="PHDMMB_DMCYC"
                qualifiers={{
                  guiLabel: "Number of density modif. cycles",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="PHDMMB_BIGCYC"
                qualifiers={{
                  guiLabel: "Number of model building cycles",
                }}
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
                  guiLabel: "Other hand CC threshold",
                }}
              />
            </FieldRow>
            <CCP4i2TaskElement
              {...props}
              itemName="SUBSTRDET_THRESHOLD_WEAK"
              qualifiers={{
                guiLabel:
                  "Use thorough building if CFOM from substr. detection is smaller than",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="ARGUMENTS_SHELXE"
              qualifiers={{
                guiLabel: "Custom program arguments (comma separated)",
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
              onChange={onToggle(setUseComb)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="MB_PROGRAM"
              qualifiers={{ guiLabel: "Model building program" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="KEYWORDS_MB"
              qualifiers={{
                guiLabel: "Custom keywords for building program",
              }}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: " " }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                {...props}
                itemName="MBREF_BIGCYC"
                qualifiers={{
                  guiLabel: "Number of building cycles",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="MBREF_EXCLUDE_FREE"
                qualifiers={{
                  guiLabel:
                    "Exclude the free reflections in model building",
                }}
              />
            </CCP4i2ContainerElement>
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
              qualifiers={{ guiLabel: "Exclude free reflections" }}
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
