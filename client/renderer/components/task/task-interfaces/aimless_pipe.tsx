import React, { useCallback, useEffect, useMemo, useState } from "react";
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { FieldRow } from "../task-elements/field-row";

/** Normalize CBoolean values - server may return boolean or string */
const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

/**
 * Hook to manage a CBoolean toggle with local state for immediate UI response.
 * Standard pattern: local state + useEffect sync from server + onChange callback.
 */
function useBoolToggle(
  useTaskItemFn: (name: string) => { value: any },
  itemName: string
) {
  const { value: raw } = useTaskItemFn(itemName);
  const [value, setValue] = useState(() => isTruthy(raw));
  useEffect(() => setValue(isTruthy(raw)), [raw]);
  const onChange = useCallback(async (item: any) => {
    setValue(isTruthy(item._value));
  }, []);
  return { value, onChange };
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // ─── Auto-sync: REFERENCE_FOR_AIMLESS tracks MODE ──────────────────
  const { syncTo: syncToAimlessRef } = useTaskItem("REFERENCE_FOR_AIMLESS");

  // ─── Enum values for conditional rendering ─────────────────────────
  const { value: mode } = useTaskItem("MODE");
  const { value: chooseMode } = useTaskItem("CHOOSE_MODE");
  const { value: referenceDataset } = useTaskItem("REFERENCE_DATASET");
  const { value: sdcorrOptions } = useTaskItem("SDCORRECTION_OPTIONS");
  const { value: scalingProtocol } = useTaskItem("SCALING_PROTOCOL");
  const { value: scalesRotType } = useTaskItem("SCALES_ROTATION_TYPE");
  const { value: scalesBrotType } = useTaskItem("SCALES_BROTATION_TYPE");
  const { value: scalesTileType } = useTaskItem("SCALES_TILETYPE");
  const { value: runMode } = useTaskItem("RUN_MODE");

  // ─── Boolean toggles for conditional visibility ────────────────────
  // Tab 2: Important Options
  const analysisOverride = useBoolToggle(useTaskItem, "ANALYSIS_OVERRIDE");
  const intensitiesOverride = useBoolToggle(
    useTaskItem,
    "INTENSITIES_OVERRIDE"
  );
  const sdcorrOverride = useBoolToggle(useTaskItem, "SDCORRECTION_OVERRIDE");
  const sdcorrRefine = useBoolToggle(useTaskItem, "SDCORRECTION_REFINE");
  const scalingDetails = useBoolToggle(useTaskItem, "SCALING_DETAILS");
  const outputUnmerged = useBoolToggle(useTaskItem, "OUTPUT_UNMERGED");

  // Tab 3: Additional Options
  const removeLattice = useBoolToggle(
    useTaskItem,
    "REMOVE_LATTICE_CENTERING"
  );
  const keepLattice = useBoolToggle(useTaskItem, "KEEP_LATTICE_CENTERING");
  const outlierOverride = useBoolToggle(useTaskItem, "OUTLIER_OVERRIDE");
  const expertOptions = useBoolToggle(useTaskItem, "EXPERT_OPTIONS");

  // ─── Auto-sync REFERENCE_FOR_AIMLESS with MODE ────────────────────
  useEffect(() => {
    syncToAimlessRef(mode === "MATCH");
  }, [mode, syncToAimlessRef]);

  // ─── Visibility helpers ────────────────────────────────────────────
  const vis = useMemo(
    () => ({
      // Tab 1: Symmetry modes
      isChooseMode: () => mode === "CHOOSE",
      isChooseSolution: () =>
        mode === "CHOOSE" && chooseMode === "SOLUTION_NO",
      isChooseSpacegroup: () =>
        mode === "CHOOSE" &&
        ["SPACEGROUP", "REINDEX_SPACE"].includes(chooseMode),
      isReindexSpace: () =>
        mode === "CHOOSE" && chooseMode === "REINDEX_SPACE",
      isChooseLauegroup: () =>
        mode === "CHOOSE" && chooseMode === "LAUEGROUP",
      // Tab 2: SD correction
      showSimilarity: () => sdcorrOptions === "SIMILAR",
      // Tab 3: Scaling protocol
      hasRotationScales: () =>
        ["ROTATION", "SECONDARY"].includes(scalingProtocol),
      hasSecondaryBeam: () => scalingProtocol === "SECONDARY",
      showTileXY: () =>
        scalingProtocol === "SECONDARY" &&
        !["DEFAULT", "NONE"].includes(scalesTileType),
      rotSpacing: () => scalesRotType === "SPACING",
      rotNbins: () => scalesRotType === "NBINS",
      brotSpacing: () => scalesBrotType === "SPACING",
      brotNbins: () => scalesBrotType === "NBINS",
      // Tab 3: Lattice centering
      showLatticeThreshold: () =>
        !removeLattice.value && !keepLattice.value,
      // Tab 3: Run mode
      showBatchList: () => runMode === "BYRANGE",
    }),
    [
      mode,
      chooseMode,
      sdcorrOptions,
      scalingProtocol,
      scalesRotType,
      scalesBrotType,
      scalesTileType,
      removeLattice.value,
      keepLattice.value,
      runMode,
    ]
  );

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ═══════════════ Tab 1: Input Data ═══════════════ */}
        <CCP4i2Tab label="Input Data">
          {/* Select unmerged data files */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Select unmerged data files" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="UNMERGEDFILES" {...props} />
          </CCP4i2ContainerElement>

          {/* Resolution */}
          <CCP4i2TaskElement
            itemName="RESOLUTION_RANGE"
            {...props}
            qualifiers={{ guiLabel: "Resolution range (Å)" }}
          />
          <CCP4i2TaskElement
            itemName="POINTLESS_USE_RESOLUTION_RANGE"
            {...props}
            qualifiers={{
              guiLabel:
                "use explicit resolution range in symmetry determination as well as in scaling",
            }}
          />
          <CCP4i2TaskElement
            itemName="AUTOCUTOFF"
            {...props}
            qualifiers={{
              guiLabel:
                "automatically cut resolution range based on a first Aimless run",
            }}
          />

          {/* Symmetry determination */}
          <CCP4i2TaskElement
            itemName="MODE"
            {...props}
            qualifiers={{ guiLabel: "Options for symmetry determination" }}
          />

          {/* Choose mode options (visible when MODE === "CHOOSE") */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Options for choice of space group or Laue group",
            }}
            containerHint="BlockLevel"
            visibility={vis.isChooseMode}
          >
            <FieldRow>
              <CCP4i2TaskElement
                itemName="CHOOSE_MODE"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
              <CCP4i2TaskElement
                itemName="CHOOSE_SPACEGROUP"
                {...props}
                visibility={vis.isChooseSpacegroup}
              />
              <CCP4i2TaskElement
                itemName="CHOOSE_SOLUTION_NO"
                {...props}
                visibility={vis.isChooseSolution}
              />
              <CCP4i2TaskElement
                itemName="CHOOSE_LAUEGROUP"
                {...props}
                visibility={vis.isChooseLauegroup}
              />
            </FieldRow>
            <CCP4i2TaskElement
              itemName="REINDEX_OPERATOR"
              {...props}
              qualifiers={{ guiLabel: "Reindex operator" }}
              visibility={vis.isReindexSpace}
            />
          </CCP4i2ContainerElement>

          {/* Optional reference data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Optional reference data to resolve indexing ambiguity and space group (also later option to scale against this reference set)",
            }}
            containerHint="FolderLevel"
          >
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">Reference data are</Typography>
              <Box sx={{ width: "12rem" }}>
                <CCP4i2TaskElement
                  itemName="REFERENCE_DATASET"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                {mode === "MATCH"
                  ? "and MUST be defined in next line"
                  : "and is optionally defined in next line"}
              </Typography>
            </Box>

            {referenceDataset === "XYZ" ? (
              <CCP4i2TaskElement
                itemName="XYZIN_REF"
                {...props}
                qualifiers={{ guiLabel: "Atomic model" }}
              />
            ) : (
              <CCP4i2TaskElement
                itemName="HKLIN_REF"
                {...props}
                qualifiers={{ guiLabel: "Reflections" }}
              />
            )}
          </CCP4i2ContainerElement>

          {/* Optional existing FreeR set */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Optional existing FreeR set, define to copy or extend if necessary",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="FREERFLAG"
              {...props}
              qualifiers={{ guiLabel: "Free R set" }}
            />
            <CCP4i2TaskElement
              itemName="CUTRESOLUTION"
              {...props}
              qualifiers={{
                guiLabel:
                  "Cut resolution of FreeR set if necessary to match the data",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ═══════════════ Tab 2: Important Options ═══════════════ */}
        <CCP4i2Tab label="Important Options">
          {/* Pointless symmetry options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Options for symmetry determination in Pointless",
            }}
            containerHint="FolderLevel"
          >
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">
                Maximum resolution for scoring set by CC(1/2) in P1 &gt;
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="CCHALFLIMIT"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                [usual method, default 0.6]
              </Typography>
            </Box>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">
                Maximum resolution for scoring set by I/sigma(I) &gt;
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="ISIGLIMIT"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                [fall-back method, default 6]
              </Typography>
            </Box>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">
                Tolerance for comparing lattices (degrees or equivalent on
                lengths)
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="TOLERANCE"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                [default 2.0]
              </Typography>
            </Box>
          </CCP4i2ContainerElement>

          {/* Aimless scaling and merging options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Options for scaling and merging in Aimless",
            }}
            containerHint="FolderLevel"
          >
            {/* ── Analysis override ── */}
            <CCP4i2TaskElement
              itemName="ANALYSIS_OVERRIDE"
              {...props}
              qualifiers={{
                guiLabel:
                  "override default parameters for estimation of maximum resolution",
              }}
              onChange={analysisOverride.onChange}
            />
            {analysisOverride.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: " " }}
                containerHint="BlockLevel"
              >
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    Set minimum information content (average bits/reflection)
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="INFOCONTENTTHRESHOLD"
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
                  <Typography variant="body1">Set minimum CC(1/2)</Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="CCMINIMUM"
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
                  <Typography variant="body1">
                    Set minimum anomalous CC(1/2)
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="CCANOMMINIMUM"
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
                  <Typography variant="body1">Set minimum MnI/sigI</Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="ISIGMINIMUM"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </CCP4i2ContainerElement>
            )}

            {/* ── Intensities and partials override ── */}
            <CCP4i2TaskElement
              itemName="INTENSITIES_OVERRIDE"
              {...props}
              qualifiers={{
                guiLabel:
                  "override default parameters for selection of intensities and treatment of partials",
              }}
              onChange={intensitiesOverride.onChange}
            />
            {intensitiesOverride.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: " " }}
                containerHint="BlockLevel"
              >
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">Use</Typography>
                  <Box sx={{ width: "16rem" }}>
                    <CCP4i2TaskElement
                      itemName="INTENSITIES_OPTIONS"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                    (profile-fitted for weak intensities, summation for strong)
                  </Typography>
                </Box>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    itemName="PARTIALS_TEST"
                    {...props}
                    qualifiers={{
                      guiLabel:
                        "Only accept partials with total fraction between",
                    }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="PARTIALS_FRACLOW"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">and</Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="PARTIALS_FRACHIGH"
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
                  <CCP4i2TaskElement
                    itemName="PARTIALS_SCALE"
                    {...props}
                    qualifiers={{ guiLabel: "Scale partials in range" }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="PARTIALS_SCALE_MIN"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    to lower acceptance limit
                  </Typography>
                </Box>
                <CCP4i2TaskElement
                  itemName="ACCEPT_OVERLOADS"
                  {...props}
                  qualifiers={{ guiLabel: "accept overloaded observations" }}
                />
                <CCP4i2TaskElement
                  itemName="ACCEPT_EDGES"
                  {...props}
                  qualifiers={{
                    guiLabel:
                      "accept observations on edge of tile or detector",
                  }}
                />
                <CCP4i2TaskElement
                  itemName="ACCEPT_XDS_MISFITS"
                  {...props}
                  qualifiers={{
                    guiLabel:
                      "accept observations flagged by XDS as outliers (MISFITS)",
                  }}
                />
              </CCP4i2ContainerElement>
            )}

            {/* ── SD correction override ── */}
            <CCP4i2TaskElement
              itemName="SDCORRECTION_OVERRIDE"
              {...props}
              qualifiers={{
                guiLabel: "override default parameters for SD correction",
              }}
              onChange={sdcorrOverride.onChange}
            />
            {sdcorrOverride.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: " " }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement
                  itemName="SDCORRECTION_REFINE"
                  {...props}
                  qualifiers={{ guiLabel: "refine SDcorrection parameters" }}
                  onChange={sdcorrRefine.onChange}
                />
                {sdcorrRefine.value && (
                  <>
                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        gap: 1,
                        flexWrap: "wrap",
                      }}
                    >
                      <Box sx={{ width: "12rem" }}>
                        <CCP4i2TaskElement
                          itemName="SDCORRECTION_OPTIONS"
                          {...props}
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        {sdcorrOptions === "SAME"
                          ? "same parameters for each run"
                          : sdcorrOptions === "SIMILAR"
                            ? "similar parameters for each run"
                            : "different parameters for each run"}
                      </Typography>
                    </Box>
                    {sdcorrOptions === "SIMILAR" && (
                      <Box sx={{ pl: 3 }}>
                        <Typography
                          variant="body2"
                          sx={{ fontStyle: "italic", mb: 0.5 }}
                        >
                          Similarity restraint SDs:
                        </Typography>
                        <FieldRow equalWidth={false} size="xs">
                          <CCP4i2TaskElement
                            itemName="SDCORRECTION_SIMILARITY_SDFAC"
                            {...props}
                            qualifiers={{ guiLabel: "sdFac" }}
                          />
                          <CCP4i2TaskElement
                            itemName="SDCORRECTION_SIMILARITY_SDB"
                            {...props}
                            qualifiers={{ guiLabel: "sdB" }}
                          />
                          <CCP4i2TaskElement
                            itemName="SDCORRECTION_SIMILARITY_SDADD"
                            {...props}
                            qualifiers={{ guiLabel: "sdAdd" }}
                          />
                        </FieldRow>
                      </Box>
                    )}
                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        gap: 1,
                        flexWrap: "wrap",
                      }}
                    >
                      <CCP4i2TaskElement
                        itemName="SDCORRECTION_FIXSDB"
                        {...props}
                        qualifiers={{ guiLabel: "fix sdB parameters" }}
                      />
                      <Typography variant="body1">
                        SD to tie sdB to zero
                      </Typography>
                      <Box sx={{ width: "8rem" }}>
                        <CCP4i2TaskElement
                          itemName="SDCORRECTION_TIESDB_SD"
                          {...props}
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                    </Box>
                  </>
                )}
                {!sdcorrRefine.value && (
                  <>
                    <CCP4i2TaskElement
                      itemName="SDCORRECTION_SET"
                      {...props}
                      qualifiers={{
                        guiLabel: "Set SD correction parameters",
                      }}
                    />
                    <FieldRow equalWidth={false} size="xs">
                      <CCP4i2TaskElement
                        itemName="SDCORRECTION_SDFAC"
                        {...props}
                        qualifiers={{ guiLabel: "sdFac" }}
                      />
                      <CCP4i2TaskElement
                        itemName="SDCORRECTION_SDB"
                        {...props}
                        qualifiers={{ guiLabel: "sdB" }}
                      />
                      <CCP4i2TaskElement
                        itemName="SDCORRECTION_SDADD"
                        {...props}
                        qualifiers={{ guiLabel: "sdAdd" }}
                      />
                    </FieldRow>
                  </>
                )}
              </CCP4i2ContainerElement>
            )}

            {/* ── Scaling details override ── */}
            <CCP4i2TaskElement
              itemName="SCALING_DETAILS"
              {...props}
              qualifiers={{
                guiLabel: "override default parameters for scaling details",
              }}
              onChange={scalingDetails.onChange}
            />
            {scalingDetails.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: " " }}
                containerHint="BlockLevel"
              >
                {/* Cycles */}
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    itemName="CYCLES_FLAG"
                    {...props}
                    qualifiers={{ guiLabel: "Refine scale factors for" }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="CYCLES_N"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">cycles</Typography>
                </Box>

                {/* Parallel */}
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    itemName="PARALLEL"
                    {...props}
                    qualifiers={{ guiLabel: "use multiple processors" }}
                  />
                  <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                    determined from number of reflections
                  </Typography>
                </Box>

                {/* Selection thresholds */}
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    itemName="SELECT1"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="SELECT_IOVSDMIN"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    minimum I/sd for 1st round scaling
                  </Typography>
                  <CCP4i2TaskElement
                    itemName="SELECT2"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="SELECT_EMIN"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    minimum E for 2nd round scaling
                  </Typography>
                </Box>

                {/* Tie restraints */}
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    itemName="TIE_ROTATION"
                    {...props}
                    qualifiers={{
                      guiLabel:
                        "Restrain neighbouring scale factors on rotation axis with SD",
                    }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="TIE_ROTATION_SD"
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
                  <CCP4i2TaskElement
                    itemName="TIE_SURFACE"
                    {...props}
                    qualifiers={{
                      guiLabel:
                        "Restrain surface parameters to a sphere with SD",
                    }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="TIE_SURFACE_SD"
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
                  <CCP4i2TaskElement
                    itemName="TIE_BFACTOR"
                    {...props}
                    qualifiers={{
                      guiLabel:
                        "Restrain neighbouring B-factors on rotation axis with SD",
                    }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="TIE_BFACTOR_SD"
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
                  <CCP4i2TaskElement
                    itemName="TIE_BZERO"
                    {...props}
                    qualifiers={{
                      guiLabel: "Restrain B-factors to zero with SD",
                    }}
                  />
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="TIE_BZERO_SD"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </CCP4i2ContainerElement>
            )}

            {/* ── Output options ── */}
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 2,
                flexWrap: "wrap",
                mt: 1,
              }}
            >
              <CCP4i2TaskElement
                itemName="OUTPUT_UNMERGED"
                {...props}
                qualifiers={{
                  guiLabel: "output unmerged data as well as merged",
                }}
                onChange={outputUnmerged.onChange}
              />
              {outputUnmerged.value && (
                <CCP4i2TaskElement
                  itemName="ORIGINAL_HKL"
                  {...props}
                  qualifiers={{
                    guiLabel: "output original hkl instead of reduced hkl",
                  }}
                />
              )}
            </Box>
          </CCP4i2ContainerElement>

          {/* FreeR extension options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Options for FreeR set extension" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" sx={{ fontStyle: "italic" }}>
              If you are extending an existing FreeR set, it must match the
              observed data in unit cell and Laue group
            </Typography>
            <Typography variant="body2" sx={{ fontStyle: "italic" }}>
              The unit cells should match to the lower resolution of the two
              datasets
            </Typography>
            <CCP4i2TaskElement
              itemName="OVERRIDE_CELL_DIFFERENCE"
              {...props}
              qualifiers={{
                guiLabel:
                  "allow existing freeR set to have different unit cells",
              }}
            />
            <Typography
              variant="body2"
              sx={{ fontStyle: "italic", color: "error.main" }}
            >
              Be sure you know what you are doing: the cells must be very
              similar even if outside the test limits
            </Typography>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ═══════════════ Tab 3: Additional Options ═══════════════ */}
        <CCP4i2Tab label="Additional Options">
          {/* Options for symmetry determination in Pointless */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Options for symmetry determination in Pointless",
            }}
            containerHint="FolderLevel"
          >
            <Typography
              variant="body2"
              sx={{ fontStyle: "italic", color: "primary.main", mb: 0.5 }}
            >
              Choice of cell setting conventions:
            </Typography>
            <CCP4i2TaskElement
              itemName="SET_SETTING"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </CCP4i2ContainerElement>

          {/* Option to remove centred lattice absences in Pointless */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Option to remove centred lattice absences in Pointless",
            }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" sx={{ fontStyle: "italic" }}>
              If the input data belong to a primitive lattice (P), the data
              are checked for additional centred lattice symmetry
            </Typography>
            <Typography variant="body2" sx={{ fontStyle: "italic" }}>
              By default, extra centred lattice reflections will be removed
              if the estimated &quot;probability&quot; of a centred lattice
            </Typography>
            <Typography variant="body2" sx={{ fontStyle: "italic", mb: 1 }}>
              is greater than the THRESHOLD set here or by default. KEEP and
              REMOVE options are unconditional
            </Typography>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="LATTICE_CENTERING_THRESHOLD"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">
                Probability threshold for removing centred lattice
                reflections
              </Typography>
            </Box>
            <CCP4i2TaskElement
              itemName="KEEP_LATTICE_CENTERING"
              {...props}
              qualifiers={{
                guiLabel: "Always keep centred lattice reflections (KEEP)",
              }}
              onChange={keepLattice.onChange}
            />
            <CCP4i2TaskElement
              itemName="REMOVE_LATTICE_CENTERING"
              {...props}
              qualifiers={{
                guiLabel:
                  "Unconditionally remove centred lattice absences (REMOVE)",
              }}
              onChange={removeLattice.onChange}
            />
            {removeLattice.value && (
              <Box sx={{ pl: 3 }}>
                <CCP4i2TaskElement
                  itemName="LATTICE_CENTERING"
                  {...props}
                  qualifiers={{ guiLabel: "Lattice centering type" }}
                />
              </Box>
            )}
          </CCP4i2ContainerElement>

          {/* Options for scaling in Aimless */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Options for scaling in Aimless" }}
            containerHint="FolderLevel"
          >
            {/* Scale [dropdown] with relative B-factor [checkbox] */}
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">Scale</Typography>
              <Box sx={{ width: "20rem" }}>
                <CCP4i2TaskElement
                  itemName="SCALING_PROTOCOL"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              {scalingProtocol !== "ONLYMERGE" && (
                <CCP4i2TaskElement
                  itemName="BFACTOR_SCALE"
                  {...props}
                  qualifiers={{ guiLabel: "with relative B-factor" }}
                />
              )}
            </Box>

            {/* Rotation scales (ROTATION or SECONDARY protocols) */}
            {vis.hasRotationScales() && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: " " }}
                containerHint="BlockLevel"
              >
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    Define scale ranges along rotation axis by
                  </Typography>
                  <Box sx={{ width: "12rem" }}>
                    <CCP4i2TaskElement
                      itemName="SCALES_ROTATION_TYPE"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <CCP4i2TaskElement
                    itemName="SCALES_ROTATION_SPACING"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                    visibility={vis.rotSpacing}
                  />
                  <CCP4i2TaskElement
                    itemName="SCALES_ROTATION_NBINS"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                    visibility={vis.rotNbins}
                  />
                  {scalesRotType === "SPACING" && (
                    <Typography variant="body1">degrees</Typography>
                  )}
                </Box>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    Define B-factor ranges along rotation axis by
                  </Typography>
                  <Box sx={{ width: "12rem" }}>
                    <CCP4i2TaskElement
                      itemName="SCALES_BROTATION_TYPE"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <CCP4i2TaskElement
                    itemName="SCALES_BROTATION_SPACING"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                    visibility={vis.brotSpacing}
                  />
                  <CCP4i2TaskElement
                    itemName="SCALES_BROTATION_NBINS"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                    visibility={vis.brotNbins}
                  />
                  {scalesBrotType === "SPACING" && (
                    <Typography variant="body1">degrees</Typography>
                  )}
                </Box>
              </CCP4i2ContainerElement>
            )}

            {/* Secondary beam correction (SECONDARY protocol only) */}
            {vis.hasSecondaryBeam() && (
              <>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    Maximum order of spherical harmonics for secondary beam
                    correction (eg 4 or 6)
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="SCALES_SECONDARY_NSPHHARMONICS"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
                <CCP4i2ContainerElement
                  {...props}
                  itemName=""
                  qualifiers={{ guiLabel: " " }}
                  containerHint="BlockLevel"
                >
                  <Box
                    sx={{
                      display: "flex",
                      alignItems: "center",
                      gap: 1,
                      flexWrap: "wrap",
                    }}
                  >
                    <Typography variant="body1">
                      Tile scaling for CCD detectors
                    </Typography>
                    <Box sx={{ width: "12rem" }}>
                      <CCP4i2TaskElement
                        itemName="SCALES_TILETYPE"
                        {...props}
                        qualifiers={{ guiLabel: " " }}
                      />
                    </Box>
                  </Box>
                  {vis.showTileXY() && (
                    <FieldRow equalWidth={false} size="xs">
                      <CCP4i2TaskElement
                        itemName="SCALES_NTILEX"
                        {...props}
                        qualifiers={{ guiLabel: "Tiles on X" }}
                      />
                      <CCP4i2TaskElement
                        itemName="SCALES_NTILEY"
                        {...props}
                        qualifiers={{ guiLabel: "Tiles on Y" }}
                      />
                    </FieldRow>
                  )}
                </CCP4i2ContainerElement>
              </>
            )}

            {/* Outlier rejection (within scaling section) */}
            <CCP4i2TaskElement
              itemName="OUTLIER_OVERRIDE"
              {...props}
              qualifiers={{
                guiLabel:
                  "override default parameters for outlier rejection",
              }}
              onChange={outlierOverride.onChange}
            />
            {outlierOverride.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: " " }}
                containerHint="BlockLevel"
              >
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    Reject outliers if &gt;
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="OUTLIER_SDMAX"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    from mean, or &gt;
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="OUTLIER_SDMAX2"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    if 2 observations
                  </Typography>
                </Box>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    Reject outliers between I+ and I- if &gt;
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="OUTLIER_SDMAXALL"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">,</Typography>
                  <CCP4i2TaskElement
                    itemName="OUTLIER_SDMAXALL_ADJUST"
                    {...props}
                    qualifiers={{
                      guiLabel:
                        "increase for large anomalous differences",
                    }}
                  />
                </Box>
                <CCP4i2TaskElement
                  itemName="OUTLIER_COMBINE"
                  {...props}
                  qualifiers={{
                    guiLabel: "compare outliers within each datasets",
                  }}
                />
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    Set maximum E to reject unreasonably large intensities
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="OUTLIER_EMAX"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </CCP4i2ContainerElement>
            )}
          </CCP4i2ContainerElement>

          {/* Options for both Pointless and Aimless */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Options for both Pointless and Aimless",
            }}
            containerHint="FolderLevel"
          >
            <Typography
              variant="body2"
              sx={{ fontStyle: "italic", color: "primary.main", mb: 0.5 }}
            >
              Override automatic definition of runs to mark discontinuities
              in data
            </Typography>
            <CCP4i2TaskElement
              itemName="RUN_MODE"
              {...props}
              qualifiers={{
                guiLabel: "Run selection options",
                guiMode: "radio",
              }}
            />
            <CCP4i2TaskElement
              itemName="RUN_BATCHLIST"
              {...props}
              qualifiers={{
                guiLabel:
                  "Batch ranges to define runs (after any renumbering in Pointless), with optional high resolution limit",
              }}
              visibility={vis.showBatchList}
            />
          </CCP4i2ContainerElement>

          {/* Expert options, not for normal use */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Expert options, not for normal use",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="EXPERT_OPTIONS"
              {...props}
              qualifiers={{ guiLabel: "show expert options" }}
              onChange={expertOptions.onChange}
            />
            {expertOptions.value && (
              <Box sx={{ pl: 3 }}>
                <CCP4i2TaskElement
                  itemName="ALLOW_NONCHIRAL"
                  {...props}
                  qualifiers={{
                    guiLabel: "Allow non-chiral space groups",
                  }}
                />
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">
                    CC(1/2) inner shell disaster limit
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="CCHALFDISASTER"
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
                  <Typography variant="body1">
                    Minimum multiplicity, inner disaster limit
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="MULTIPLICITYDISASTER"
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
