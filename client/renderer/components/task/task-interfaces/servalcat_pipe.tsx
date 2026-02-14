import React, { useMemo } from "react";
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";
import { useFreeRWarning } from "../../../providers/run-check-provider";
import { CChainSelectElement } from "../task-elements/cchainselect";

/** Normalize CBoolean values - server may return boolean or string */
const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const {
    useTaskItem,
    createPeerTask,
    validation,
    useFileDigest,
  } = useJob(job.id);

  // --- Reflection data ---
  const { item: HKLINItem, value: HKLINValue } = useTaskItem("HKLIN");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  useFreeRWarning({
    job,
    taskName: "servalcat_pipe",
    freeRFlag,
    validation,
    createPeerTask,
  });

  // --- Values for visibility conditions ---
  const { value: dataMethod } = useTaskItem("DATA_METHOD");
  const { value: mergedOrUnmerged } = useTaskItem("MERGED_OR_UNMERGED");
  const { value: hydrUse } = useTaskItem("HYDR_USE");
  const { value: addWaters } = useTaskItem("ADD_WATERS");
  const { value: useAnomalous } = useTaskItem("USEANOMALOUS");
  const { value: weightOpt } = useTaskItem("WEIGHT_OPT");
  const { value: bfacSetUse } = useTaskItem("BFACSETUSE");
  const { value: randomizeUse } = useTaskItem("RANDOMIZEUSE");
  const { value: occupancyRefinement } = useTaskItem("OCCUPANCY_REFINEMENT");
  const { value: occupancyGroups } = useTaskItem("OCCUPANCY_GROUPS");
  const { value: occupancyComplete } = useTaskItem("OCCUPANCY_COMPLETE");
  const { value: occupancyIncomplete } = useTaskItem("OCCUPANCY_INCOMPLETE");
  const { value: useJelly } = useTaskItem("USE_JELLY");
  const { value: prosmartProteinToggle } = useTaskItem(
    "prosmartProtein.TOGGLE"
  );
  const { value: prosmartProteinAdvanced } = useTaskItem(
    "prosmartProtein.ADVANCED"
  );
  const { value: prosmartNucleicAcidToggle } = useTaskItem(
    "prosmartNucleicAcid.TOGGLE"
  );
  const { value: prosmartNucleicAcidAdvanced } = useTaskItem(
    "prosmartNucleicAcid.ADVANCED"
  );
  const { value: libgToggle } = useTaskItem("libg.TOGGLE");
  const { value: libgOption } = useTaskItem("libg.OPTION");
  const { value: libgAdvanced } = useTaskItem("libg.ADVANCED");
  const { value: platonyzerToggle } = useTaskItem("platonyzer.TOGGLE");
  const { value: metalCoordRun } = useTaskItem(
    "metalCoordPipeline.RUN_METALCOORD"
  );
  const { value: metalCoordGenOrUse } = useTaskItem(
    "metalCoordPipeline.GENERATE_OR_USE"
  );
  const { value: metalCoordAdvanced } = useTaskItem(
    "metalCoordPipeline.TOGGLE_ADVANCED"
  );
  const { value: scatteringFactors } = useTaskItem("SCATTERING_FACTORS");
  const { value: resCustom } = useTaskItem("RES_CUSTOM");
  const { value: blurUse } = useTaskItem("BLURUSE");
  const { value: runAdpAnalysis } = useTaskItem("RUN_ADP_ANALYSIS");
  const { value: runCoordAdpDev } = useTaskItem(
    "monitor.RUN_COORDADPDEV_ANALYSIS"
  );

  // Chain composition from XYZIN digest (for ProSMART/libg section visibility)
  const { item: XYZINItem, value: XYZINValue } = useTaskItem("XYZIN");
  const xyzinDigestPath =
    XYZINValue?.dbFileId && XYZINItem?._objectPath
      ? XYZINItem._objectPath
      : "";
  const { data: xyzinDigest } = useFileDigest(xyzinDigestPath);
  const xyzinComposition = xyzinDigest?.composition;
  const hasProteinChains =
    Array.isArray(xyzinComposition?.peptides) &&
    xyzinComposition.peptides.length > 0;
  const hasNucleotideChains =
    Array.isArray(xyzinComposition?.nucleics) &&
    xyzinComposition.nucleics.length > 0;
  const hasMetalSites =
    Array.isArray(xyzinComposition?.ligands) &&
    xyzinComposition.ligands.some(
      (l: any) => l.isMetal || l.type === "metal"
    );

  // Derived state
  const isXtal = dataMethod === "xtal";
  const isSpa = dataMethod === "spa";
  const isMerged = mergedOrUnmerged === "merged";

  const intensitiesAvailable = useMemo(() => {
    return [1, 3].includes(HKLINValue?.contentFlag);
  }, [HKLINValue?.contentFlag]);

  return (
    <Paper>
      <CCP4i2Tabs>
        {/* ============================================================
            TAB 1: INPUT DATA
            ============================================================ */}
        <CCP4i2Tab label="Input Data" key="input">
          {/* Data source */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Data source" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="DATA_METHOD"
              qualifiers={{ guiLabel: "Data type:" }}
            />
            {isXtal && (
              <CCP4i2TaskElement
                {...props}
                itemName="MERGED_OR_UNMERGED"
                qualifiers={{ guiLabel: "Diffraction data are" }}
              />
            )}
          </CCP4i2ContainerElement>

          {/* Main inputs */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Main inputs" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="XYZIN" />

            {/* X-ray merged inputs */}
            {isXtal && isMerged && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="HKLIN"
                  qualifiers={{ guiLabel: "Reflections" }}
                />
                {intensitiesAvailable ? (
                  <CCP4i2TaskElement
                    {...props}
                    itemName="F_SIGF_OR_I_SIGI"
                    qualifiers={{ guiLabel: "Refinement against" }}
                  />
                ) : (
                  <Typography
                    variant="body1"
                    sx={{ fontWeight: "medium", mt: 1 }}
                  >
                    Using <strong>amplitudes</strong>
                  </Typography>
                )}
                <CCP4i2TaskElement {...props} itemName="FREERFLAG" />
              </>
            )}

            {/* X-ray unmerged inputs */}
            {isXtal && !isMerged && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="HKLIN_UNMERGED"
                  qualifiers={{ guiLabel: "Unmerged reflection data" }}
                />
                <CCP4i2TaskElement {...props} itemName="FREERFLAG" />
              </>
            )}

            {/* SPA map inputs */}
            {isSpa && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="MAPIN1"
                  qualifiers={{ guiLabel: "Half map 1" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="MAPIN2"
                  qualifiers={{ guiLabel: "Half map 2" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="MAPMASK"
                  qualifiers={{ guiLabel: "Mask" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="RES_MIN"
                  qualifiers={{ guiLabel: "Resolution:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="MASK_RADIUS"
                  qualifiers={{ guiLabel: "Mask radius:" }}
                />
              </>
            )}

            {isXtal && (
              <CCP4i2TaskElement
                {...props}
                itemName="USE_TWIN"
                qualifiers={{ guiLabel: "Twin refinement" }}
              />
            )}
          </CCP4i2ContainerElement>

          {/* Additional geometry dictionaries */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Additional geometry dictionaries" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="DICT_LIST" />
          </CCP4i2ContainerElement>

          {/* Options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="NCYCLES"
              qualifiers={{ guiLabel: "Number of refinement cycles:" }}
            />

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="HYDR_USE"
                qualifiers={{
                  guiLabel: "Use riding hydrogens during refinement",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="HYDR_ALL"
                visibility={() => isTruthy(hydrUse)}
              />
            </FieldRow>

            {isXtal && (
              <>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    {...props}
                    itemName="ADD_WATERS"
                    qualifiers={{
                      guiLabel: "Add waters and then perform",
                    }}
                    sx={{ width: "auto" }}
                  />
                  {isTruthy(addWaters) && (
                    <>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="NCYCLES_AFTER_ADD_WATERS"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        further refinement cycles.
                      </Typography>
                    </>
                  )}
                </Box>
              </>
            )}
          </CCP4i2ContainerElement>

          {/* Anomalous signal (xtal only, when anomalous data available) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Anomalous signal" }}
            containerHint="FolderLevel"
            visibility={() =>
              isXtal &&
              HKLINValue?.contentFlag !== undefined &&
              [1, 2].includes(HKLINValue.contentFlag)
            }
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="USEANOMALOUS"
                qualifiers={{ guiLabel: "Use anomalous" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="USEANOMALOUSFOR"
                qualifiers={{ guiLabel: "Use for:" }}
                visibility={() => isTruthy(useAnomalous)}
              />
            </FieldRow>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 2: PARAMETERISATION
            ============================================================ */}
        <CCP4i2Tab label="Parameterisation" key="parameterisation">
          {/* Atomic displacement parameters (ADPs) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Atomic displacement parameters (ADPs)",
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
              <Box sx={{ width: "10rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="B_REFINEMENT_MODE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">ADPs</Typography>
            </Box>
          </CCP4i2ContainerElement>

          {/* Scaling */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Scaling" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="NO_SOLVENT"
              qualifiers={{
                guiLabel: "Do not consider bulk solvent contribution",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Conformer groups and occupancy refinement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Conformer groups and occupancy refinement",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_GROUPS"
              qualifiers={{
                guiLabel:
                  "Specify partial occupancy groups (alternative conformers)",
              }}
            />
            {isTruthy(occupancyGroups) && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="OCCUPANCY_SELECTION"
                />

                <CCP4i2TaskElement
                  {...props}
                  itemName="OCCUPANCY_COMPLETE"
                  qualifiers={{
                    guiLabel:
                      "Specify overlapping alternative conformer groups (constrain occupancies to sum to one)",
                  }}
                />
                {isTruthy(occupancyComplete) && (
                  <CCP4i2TaskElement
                    {...props}
                    itemName="OCCUPANCY_COMPLETE_TABLE"
                  />
                )}

                <CCP4i2TaskElement
                  {...props}
                  itemName="OCCUPANCY_INCOMPLETE"
                  qualifiers={{
                    guiLabel:
                      "Specify overlapping alternative conformer groups (occupancies sum to less than one)",
                  }}
                />
                {isTruthy(occupancyIncomplete) && (
                  <CCP4i2TaskElement
                    {...props}
                    itemName="OCCUPANCY_INCOMPLETE_TABLE"
                  />
                )}
              </>
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
                {...props}
                itemName="OCCUPANCY_REFINEMENT"
                qualifiers={{
                  guiLabel:
                    "Perform refinement of atomic occupancies every",
                }}
                sx={{ width: "auto" }}
              />
              {isTruthy(occupancyRefinement) && (
                <>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="OCCUPANCY_NCYCLE"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">cycle.</Typography>
                </>
              )}
            </Box>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 3: RESTRAINTS
            ============================================================ */}
        <CCP4i2Tab label="Restraints" key="restraints">
          {/* Weights */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Weights" }}
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
                Weigh the experimental data using
              </Typography>
              <Box sx={{ width: "10rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="WEIGHT_OPT"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">
                weight versus the restraints
              </Typography>
            </Box>
            <CCP4i2TaskElement
              {...props}
              itemName="WEIGHT"
              qualifiers={{ guiLabel: "Weight:" }}
              visibility={() => weightOpt === "MANUAL"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="WEIGHT_NO_ADJUST"
              qualifiers={{
                guiLabel: "Do not adjust weight during refinement",
              }}
            />
            <FieldRow>
              <Typography variant="body1">
                Bond RMSZ range for weight adjustment:
              </Typography>
              <Box sx={{ width: "6rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Box sx={{ width: "6rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </FieldRow>
          </CCP4i2ContainerElement>

          {/* Non-Crystallographic Symmetry (NCS) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Non-Crystallographic Symmetry (NCS)",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="USE_NCS"
              qualifiers={{
                guiLabel:
                  "Use local non-crystallographic symmetry (NCS) restraints",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Covalent links */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Covalent links" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="FIND_LINKS"
              qualifiers={{
                guiLabel:
                  "Detect and apply covalent linkages based on the current atomic coordinates",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Jelly-body */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Jelly-body" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="USE_JELLY"
              qualifiers={{ guiLabel: "Use jelly-body restraints" }}
            />
            {isTruthy(useJelly) && (
              <>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1">with sigma:</Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="JELLY_SIGMA"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    and max distance:
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="JELLY_DIST"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
                <CCP4i2TaskElement
                  {...props}
                  itemName="JELLY_ONLY"
                  qualifiers={{
                    guiLabel: "Jelly-body refinement only",
                  }}
                />
              </>
            )}
          </CCP4i2ContainerElement>

          {/* MetalCoord External Restraints for Metals */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "MetalCoord External Restraints for Metals",
            }}
            containerHint="FolderLevel"
          >
            {!hasMetalSites && (
              <Typography
                variant="body2"
                sx={{ fontStyle: "italic", fontWeight: "bold" }}
              >
                No monomers including metal sites were found in the
                input atomic model.
              </Typography>
            )}
            <CCP4i2TaskElement
              {...props}
              itemName="metalCoordPipeline.RUN_METALCOORD"
              qualifiers={{
                guiLabel:
                  "Apply MetalCoord restraints for metal sites:",
              }}
            />
            {isTruthy(metalCoordRun) && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="metalCoordPipeline.GENERATE_OR_USE"
                  qualifiers={{ guiLabel: "Restraints:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="metalCoordPipeline.LINKS"
                  qualifiers={{ guiLabel: "Link records:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="metalCoordPipeline.METALCOORD_RESTRAINTS"
                  qualifiers={{
                    guiLabel: "MetalCoord restraints file:",
                  }}
                  visibility={() => metalCoordGenOrUse === "USE"}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="metalCoordPipeline.TOGGLE_ADVANCED"
                  qualifiers={{
                    guiLabel: "Show advanced MetalCoord options",
                  }}
                />
                {isTruthy(metalCoordAdvanced) && (
                  <CCP4i2ContainerElement
                    {...props}
                    itemName=""
                    qualifiers={{
                      guiLabel: "Advanced MetalCoord options",
                    }}
                    containerHint="BlockLevel"
                  >
                    <CCP4i2TaskElement
                      {...props}
                      itemName="metalCoordWrapper.controlParameters.MAXIMUM_COORDINATION_NUMBER"
                      qualifiers={{
                        guiLabel:
                          "Maximum coordination number:",
                      }}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      itemName="metalCoordWrapper.controlParameters.MINIMUM_SAMPLE_SIZE"
                      qualifiers={{
                        guiLabel: "Minimum sample size:",
                      }}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      itemName="metalCoordWrapper.controlParameters.DISTANCE_THRESHOLD"
                      qualifiers={{
                        guiLabel: "Distance threshold:",
                      }}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      itemName="metalCoordWrapper.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD"
                      qualifiers={{
                        guiLabel:
                          "Procrustes distance threshold:",
                      }}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      itemName="metalCoordWrapper.controlParameters.IDEAL_ANGLES"
                      qualifiers={{
                        guiLabel: "Use ideal angles",
                      }}
                    />
                  </CCP4i2ContainerElement>
                )}
              </>
            )}
          </CCP4i2ContainerElement>

          {/* ProSMART External Restraints for Protein Chains */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "ProSMART External Restraints for Protein Chains",
            }}
            containerHint="FolderLevel"
          >
            {!hasProteinChains && (
              <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                Input atomic model contains no protein chains
              </Typography>
            )}
            {hasProteinChains && (
              <>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    {...props}
                    itemName="prosmartProtein.TOGGLE"
                    qualifiers={{
                      guiLabel:
                        "Generate and apply restraints for protein chains using homologous models",
                    }}
                    sx={{ width: "auto" }}
                  />
                  <Box sx={{ minWidth: "8rem", flex: "0 1 auto" }}>
                    <CChainSelectElement
                      job={job}
                      itemName="prosmartProtein.CHAINLIST_1"
                      options={xyzinComposition?.peptides || []}
                      label=" "
                      visibility={() =>
                        isTruthy(prosmartProteinToggle)
                      }
                    />
                  </Box>
                </Box>

                {isTruthy(prosmartProteinToggle) && (
                  <>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="prosmartProtein.REFERENCE_MODELS"
                    />

                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        gap: 1,
                        flexWrap: "wrap",
                      }}
                    >
                      <Typography variant="body1">Use</Typography>
                      <Box sx={{ width: "14rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.ALL_BEST"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        chain(s) from reference model(s). Minimum
                        sequence identity:
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.SEQID"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">%</Typography>
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
                        Generate restraints between
                      </Typography>
                      <Box sx={{ width: "14rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.SIDE_MAIN"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        atom-pairs. Interatomic distance range:
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.RMIN"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">to</Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.RMAX"
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
                        Apply restraints up to
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.DMAX"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        with sigma range:
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.SGMN"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">to</Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.SGMX"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        and robustness parameter (alpha)
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.ALPHA"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                    </Box>

                    <CCP4i2TaskElement
                      {...props}
                      itemName="prosmartProtein.ADVANCED"
                      qualifiers={{
                        guiLabel: "Show advanced options",
                      }}
                    />

                    {isTruthy(prosmartProteinAdvanced) && (
                      <CCP4i2ContainerElement
                        {...props}
                        itemName=""
                        qualifiers={{
                          guiLabel:
                            "Advanced ProSMART protein options",
                        }}
                        containerHint="BlockLevel"
                      >
                        <FieldRow>
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartProtein.TOGGLE_BFAC"
                            qualifiers={{
                              guiLabel:
                                "Remove restraints for atoms with B-factor above",
                            }}
                          />
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartProtein.BFAC"
                            qualifiers={{ guiLabel: " " }}
                          />
                        </FieldRow>
                        <FieldRow>
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartProtein.TOGGLE_ALT"
                            qualifiers={{
                              guiLabel:
                                "Don't generate restraints for atoms with occupancy below",
                            }}
                          />
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartProtein.OCCUPANCY"
                            qualifiers={{ guiLabel: " " }}
                          />
                        </FieldRow>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.KEYWORDS"
                          qualifiers={{
                            guiLabel:
                              "Additional ProSMART keywords:",
                            guiMode: "multiLine",
                          }}
                        />
                      </CCP4i2ContainerElement>
                    )}
                  </>
                )}
              </>
            )}
          </CCP4i2ContainerElement>

          {/* ProSMART External Restraints for Nucleic Acids */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "ProSMART External Restraints for Nucleic Acids",
            }}
            containerHint="FolderLevel"
          >
            {!hasNucleotideChains && (
              <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                Input atomic model contains no nucleotide chains
              </Typography>
            )}
            {hasNucleotideChains && (
              <>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    {...props}
                    itemName="prosmartNucleicAcid.TOGGLE"
                    qualifiers={{
                      guiLabel:
                        "Generate restraints for nucleic acid chain(s):",
                    }}
                    sx={{ width: "auto" }}
                  />
                  <Box sx={{ minWidth: "8rem", flex: "0 1 auto" }}>
                    <CChainSelectElement
                      job={job}
                      itemName="prosmartNucleicAcid.CHAINLIST_1"
                      options={xyzinComposition?.nucleics || []}
                      label=" "
                      visibility={() =>
                        isTruthy(prosmartNucleicAcidToggle)
                      }
                    />
                  </Box>
                  {isTruthy(prosmartNucleicAcidToggle) && (
                    <Typography variant="body1">
                      using homologous model(s):
                    </Typography>
                  )}
                </Box>

                {isTruthy(prosmartNucleicAcidToggle) && (
                  <>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="prosmartNucleicAcid.REFERENCE_MODELS"
                    />

                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        gap: 1,
                        flexWrap: "wrap",
                      }}
                    >
                      <Typography variant="body1">Use</Typography>
                      <Box sx={{ width: "14rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.ALL_BEST"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        chain(s) from reference model(s). Minimum
                        sequence identity:
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.SEQID"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">%</Typography>
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
                        Generate restraints between
                      </Typography>
                      <Box sx={{ width: "14rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.SIDE_MAIN"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        atom-pairs. Interatomic distance range:
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.RMIN"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">to</Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.RMAX"
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
                        Apply restraints up to
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.DMAX"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        with weight
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.WEIGHT"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                      <Typography variant="body1">
                        and robustness parameter (alpha)
                      </Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.ALPHA"
                          qualifiers={{ guiLabel: " " }}
                        />
                      </Box>
                    </Box>

                    <CCP4i2TaskElement
                      {...props}
                      itemName="prosmartNucleicAcid.ADVANCED"
                      qualifiers={{
                        guiLabel: "Show advanced options",
                      }}
                    />

                    {isTruthy(prosmartNucleicAcidAdvanced) && (
                      <CCP4i2ContainerElement
                        {...props}
                        itemName=""
                        qualifiers={{
                          guiLabel:
                            "Advanced ProSMART nucleic acid options",
                        }}
                        containerHint="BlockLevel"
                      >
                        <FieldRow>
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartNucleicAcid.TOGGLE_BFAC"
                            qualifiers={{
                              guiLabel:
                                "Remove restraints for atoms with B-factor above",
                            }}
                          />
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartNucleicAcid.BFAC"
                            qualifiers={{ guiLabel: " " }}
                          />
                        </FieldRow>
                        <FieldRow>
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartNucleicAcid.TOGGLE_ALT"
                            qualifiers={{
                              guiLabel:
                                "Don't generate restraints for atoms with occupancy below",
                            }}
                          />
                          <CCP4i2TaskElement
                            {...props}
                            itemName="prosmartNucleicAcid.OCCUPANCY"
                            qualifiers={{ guiLabel: " " }}
                          />
                        </FieldRow>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartNucleicAcid.KEYWORDS"
                          qualifiers={{
                            guiLabel:
                              "Additional ProSMART keywords:",
                            guiMode: "multiLine",
                          }}
                        />
                      </CCP4i2ContainerElement>
                    )}
                  </>
                )}
              </>
            )}
          </CCP4i2ContainerElement>

          {/* ADP Restraints */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "ADP Restraints" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="ADPR_WEIGHT"
              qualifiers={{ guiLabel: "ADP restraint weight:" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="MAX_DIST_FOR_ADP_RESTRAINT"
              qualifiers={{
                guiLabel: "Maximum distance for ADP restraint:",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="ADP_RESTRAINT_NO_LONG_RANGE"
              qualifiers={{
                guiLabel: "No long range for ADP restraint",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Platonyzer */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Platonyzer" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="platonyzer.TOGGLE"
              qualifiers={{
                guiLabel: "Use Platonyzer metal restraints",
              }}
            />
            {isTruthy(platonyzerToggle) && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="platonyzer.MODE"
                  qualifiers={{ guiLabel: "Metal type:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="platonyzer.RM_VDW"
                  qualifiers={{
                    guiLabel:
                      "Remove VDW restraints for metal sites",
                  }}
                />
              </>
            )}
          </CCP4i2ContainerElement>

          {/* libg - nucleic acid base restraints */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Libg - nucleic acid base restraints",
            }}
            containerHint="FolderLevel"
          >
            {!hasNucleotideChains && (
              <Typography variant="body2" sx={{ fontStyle: "italic" }}>
                Input atomic model contains no nucleotide chains
              </Typography>
            )}
            {hasNucleotideChains && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="libg.TOGGLE"
                  qualifiers={{
                    guiLabel:
                      "Generate nucleic acid base restraints using libg",
                  }}
                />
                {isTruthy(libgToggle) && (
                  <>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="libg.OPTION"
                      qualifiers={{
                        guiLabel: "Restraint types:",
                      }}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      itemName="libg.BP"
                      qualifiers={{
                        guiLabel: "Include base pair restraints",
                      }}
                      visibility={() => libgOption === "MANUAL"}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      itemName="libg.ADVANCED"
                      qualifiers={{
                        guiLabel: "Show advanced options",
                      }}
                    />
                    {isTruthy(libgAdvanced) && (
                      <CCP4i2TaskElement
                        {...props}
                        itemName="libg.KEYWORDS"
                        qualifiers={{
                          guiLabel: "Additional libg keywords:",
                          guiMode: "multiLine",
                        }}
                      />
                    )}
                  </>
                )}
              </>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 4: ADVANCED
            ============================================================ */}
        <CCP4i2Tab label="Advanced" key="advanced">
          {/* Resolution and experiment (xtal only) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Resolution and experiment" }}
            containerHint="FolderLevel"
            visibility={() => isXtal}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="RES_CUSTOM"
              qualifiers={{
                guiLabel: "Use custom resolution limits",
              }}
            />
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="RES_MIN"
                qualifiers={{ guiLabel: "High resolution (dmin):" }}
                visibility={() => isTruthy(resCustom)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="RES_MAX"
                qualifiers={{ guiLabel: "Low resolution (dmax):" }}
                visibility={() => isTruthy(resCustom)}
              />
            </FieldRow>

            <CCP4i2TaskElement
              {...props}
              itemName="FREERFLAG_NUMBER"
              qualifiers={{
                guiLabel: "FreeR flag number for test set:",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="SCATTERING_FACTORS"
              qualifiers={{ guiLabel: "Diffraction experiment type:" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="SCATTERING_ELECTRON"
              qualifiers={{ guiLabel: "Form factor calculation:" }}
              visibility={() => scatteringFactors === "electron"}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="USE_WORK_IN_EST"
              qualifiers={{
                guiLabel:
                  "Use work reflections in maximum likelihood parameter estimates",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Hydrogen and charge options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Hydrogen and charge options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="H_OUT"
              qualifiers={{
                guiLabel:
                  "Write hydrogen atoms in the output model",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="H_REFINE"
              qualifiers={{
                guiLabel: "Refine hydrogen positions",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="KEEP_CHARGES"
              qualifiers={{
                guiLabel:
                  "Keep charges, i.e. use scattering factor for charged atoms where relevant",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Structure model modification before refinement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Structure model modification before refinement",
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
              <CCP4i2TaskElement
                {...props}
                itemName="BFACSETUSE"
                qualifiers={{
                  guiLabel: "Reset all ADPs at start",
                }}
                sx={{ width: "auto" }}
              />
              {isTruthy(bfacSetUse) && (
                <>
                  <Typography variant="body1">
                    to fixed value:
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="BFACSET"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </>
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
              <CCP4i2TaskElement
                {...props}
                itemName="RANDOMIZEUSE"
                qualifiers={{
                  guiLabel: "Shake coordinates at start",
                }}
                sx={{ width: "auto" }}
              />
              {isTruthy(randomizeUse) && (
                <>
                  <Typography variant="body1">
                    with RMSD:
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="RANDOMIZE"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </>
              )}
            </Box>
          </CCP4i2ContainerElement>

          {/* Additional keywords */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Additional keywords" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SERVALCAT_KEYWORD_FILE"
            />
            <CCP4i2TaskElement
              {...props}
              itemName="EXTRA_SERVALCAT_OPTIONS"
              qualifiers={{
                guiLabel: "Extra servalcat command line options:",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Validation and Analysis */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Validation and Analysis" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_IRIS"
              qualifiers={{
                guiLabel: "Generate Iris validation report",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_RAMACHANDRAN"
              qualifiers={{
                guiLabel: "Generate Ramachandran plots",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_MOLPROBITY"
              qualifiers={{
                guiLabel: "Run MolProbity to analyse geometry",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="RUN_ADP_ANALYSIS"
              qualifiers={{ guiLabel: "Run ADP analysis" }}
            />
            {isTruthy(runAdpAnalysis) && (
              <>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body2">
                    Atoms with a B-value lower than{" "}
                    <em>the first quartile - factor * interquartile_range</em>
                    {" "}or higher than{" "}
                    <em>the third quartile + factor * interquartile_range</em>
                    {" "}to be reported. Factor:
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="ADP_IQR_FACTOR"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </>
            )}

            <CCP4i2TaskElement
              {...props}
              itemName="monitor.RUN_COORDADPDEV_ANALYSIS"
              qualifiers={{
                guiLabel:
                  "Run analysis of changes in coordinates and ADPs",
              }}
            />
            {isTruthy(runCoordAdpDev) && (
              <>
                <CCP4i2TaskElement
                  {...props}
                  itemName="monitor.MIN_COORDDEV"
                  qualifiers={{
                    guiLabel:
                      "Minimum shift of atom coordinates to be reported:",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="monitor.MIN_ADPDEV"
                  qualifiers={{
                    guiLabel:
                      "Minimum shift of B-values to be reported:",
                  }}
                />
              </>
            )}
          </CCP4i2ContainerElement>

          {/* SPA-specific options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "SPA-specific options" }}
            containerHint="FolderLevel"
            visibility={() => isSpa}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="PIXEL_SIZE"
              qualifiers={{
                guiLabel: "Pixel size (\u00C5/pixel):",
              }}
            />
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="POINTGROUP"
                qualifiers={{ guiLabel: "Point group:" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="IGNORE_SYMMETRY"
                qualifiers={{
                  guiLabel: "Ignore symmetry in model file",
                }}
              />
            </FieldRow>

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Helical parameters" }}
              containerHint="BlockLevel"
            >
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="TWIST"
                  qualifiers={{
                    guiLabel: "Helical twist (degrees):",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="RISE"
                  qualifiers={{
                    guiLabel: "Helical rise (\u00C5):",
                  }}
                />
              </FieldRow>
            </CCP4i2ContainerElement>

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Symmetry origin" }}
              containerHint="BlockLevel"
            >
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="CENTER_X"
                  qualifiers={{ guiLabel: "X:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CENTER_Y"
                  qualifiers={{ guiLabel: "Y:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="CENTER_Z"
                  qualifiers={{ guiLabel: "Z:" }}
                />
              </FieldRow>
            </CCP4i2ContainerElement>

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Symmetry axes" }}
              containerHint="BlockLevel"
            >
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="AXIS1_X"
                  qualifiers={{ guiLabel: "Axis 1 X:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="AXIS1_Y"
                  qualifiers={{ guiLabel: "Axis 1 Y:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="AXIS1_Z"
                  qualifiers={{ guiLabel: "Axis 1 Z:" }}
                />
              </FieldRow>
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="AXIS2_X"
                  qualifiers={{ guiLabel: "Axis 2 X:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="AXIS2_Y"
                  qualifiers={{ guiLabel: "Axis 2 Y:" }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="AXIS2_Z"
                  qualifiers={{ guiLabel: "Axis 2 Z:" }}
                />
              </FieldRow>
            </CCP4i2ContainerElement>

            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="BLURUSE"
                qualifiers={{
                  guiLabel: "Blur map with B-value:",
                }}
                sx={{ width: "auto" }}
              />
              {isTruthy(blurUse) && (
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    {...props}
                    itemName="BLUR"
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              )}
            </Box>

            <CCP4i2TaskElement
              {...props}
              itemName="CROSS_VALIDATION"
              qualifiers={{
                guiLabel: "Run cross validation with half maps",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
