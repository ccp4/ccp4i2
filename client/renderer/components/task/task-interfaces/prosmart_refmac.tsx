import React, { useCallback } from "react";
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
  const { useTaskItem, fetchDigest, createPeerTask, validation, useFileDigest } =
    useJob(job.id);

  // --- F_SIGF onChange: extract wavelength, set anomalous/twinning flags ---
  const { item: F_SIGFItem, value: F_SIGFValue } = useTaskItem("F_SIGF");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  useFreeRWarning({
    job,
    taskName: "prosmart_refmac",
    freeRFlag,
    validation,
    createPeerTask,
  });

  const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");
  const { forceUpdate: forceUpdateUSEANOMALOUS } = useTaskItem("USEANOMALOUS");
  const { forceUpdate: forceUpdateUSE_TWIN } = useTaskItem("USE_TWIN");

  // --- Values for visibility conditions ---
  const { value: refinementMode } = useTaskItem("REFINEMENT_MODE");
  const { value: hydrUse } = useTaskItem("HYDR_USE");
  const { value: addWaters } = useTaskItem("ADD_WATERS");
  const { value: useAnomalous } = useTaskItem("USEANOMALOUS");
  const { value: bRefinementMode } = useTaskItem("B_REFINEMENT_MODE");
  const { value: scaleType } = useTaskItem("SCALE_TYPE");
  const { value: solventMaskType } = useTaskItem("SOLVENT_MASK_TYPE");
  const { value: solventAdvanced } = useTaskItem("SOLVENT_ADVANCED");
  const { value: tlsMode } = useTaskItem("TLSMODE");
  const { value: bfacSetUse } = useTaskItem("BFACSETUSE");
  const { value: weightOpt } = useTaskItem("WEIGHT_OPT");
  const { value: occupancyRefinement } = useTaskItem("OCCUPANCY_REFINEMENT");
  const { value: occupancyGroups } = useTaskItem("OCCUPANCY_GROUPS");
  const { value: occupancyComplete } = useTaskItem("OCCUPANCY_COMPLETE");
  const { value: occupancyIncomplete } = useTaskItem("OCCUPANCY_INCOMPLETE");
  const { value: useNcs } = useTaskItem("USE_NCS");
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
  const { value: platonyzerToggle } = useTaskItem("platonyzer.TOGGLE");
  const { value: mapSharp } = useTaskItem("MAP_SHARP");
  const { value: mapSharpCustom } = useTaskItem("MAP_SHARP_CUSTOM");
  const { value: scatteringFactors } = useTaskItem("SCATTERING_FACTORS");
  const { value: resCustom } = useTaskItem("RES_CUSTOM");
  const { value: hdInitToggle } = useTaskItem("HD_INIT_TOGGLE");

  // Chain composition from XYZIN digest (for ProSMART section visibility)
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

  // --- F_SIGF onChange handler ---
  const handleF_SIGFChange = useCallback(async () => {
    if (!job || job.status !== 1) return;

    if (F_SIGFItem?._objectPath) {
      const digestData = await fetchDigest(F_SIGFItem._objectPath);
      const wavelength = digestData?.wavelengths?.at(-1);
      if (
        wavelength &&
        wavelength > 0 &&
        wavelength < 9 &&
        forceUpdateWAVELENGTH
      ) {
        await forceUpdateWAVELENGTH(wavelength);
      }
    }

    if (F_SIGFValue) {
      try {
        const contentFlag = F_SIGFValue.contentFlag;
        if (![1, 2].includes(contentFlag) && forceUpdateUSEANOMALOUS) {
          await forceUpdateUSEANOMALOUS(false);
        }
        if (![3].includes(contentFlag) && forceUpdateUSE_TWIN) {
          await forceUpdateUSE_TWIN(false);
        }
      } catch (error) {
        console.error("Error processing F_SIGF change:", error);
      }
    }
  }, [
    F_SIGFItem?._objectPath,
    F_SIGFValue,
    job,
    fetchDigest,
    forceUpdateWAVELENGTH,
    forceUpdateUSEANOMALOUS,
    forceUpdateUSE_TWIN,
  ]);

  return (
    <Paper>
      <CCP4i2Tabs>
        {/* ============================================================
            TAB 1: INPUT DATA
            ============================================================ */}
        <CCP4i2Tab label="Input data" key="input">
          {/* Main inputs */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Main inputs" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="XYZIN" />
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              onChange={handleF_SIGFChange}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="WAVELENGTH"
              qualifiers={{ guiLabel: "Wavelength" }}
            />
            <CCP4i2TaskElement {...props} itemName="FREERFLAG" />
          </CCP4i2ContainerElement>

          {/* Experimental phase information */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Experimental phase information" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="ABCD" />
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
              itemName="REFINEMENT_MODE"
              qualifiers={{
                guiLabel: "Refinement mode (restrained or rigid body):",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="NCYCLES"
              qualifiers={{
                guiLabel: "Number of restrained refinement cycles:",
              }}
              visibility={() => refinementMode !== "RIGID"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="NCYCRIGID"
              qualifiers={{
                guiLabel: "Number of rigid body refinement cycles:",
              }}
              visibility={() => refinementMode === "RIGID"}
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

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="ADD_WATERS"
                qualifiers={{ guiLabel: "Add waters" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="REFPRO_RSR_RWORK_LIMIT"
                qualifiers={{ guiLabel: "or lower" }}
                visibility={() => isTruthy(addWaters)}
              />
            </FieldRow>

            <CCP4i2TaskElement
              {...props}
              itemName="USE_TWIN"
              qualifiers={{ guiLabel: "Crystal is twinned" }}
            />
          </CCP4i2ContainerElement>

          {/* Anomalous signal */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Anomalous signal" }}
            containerHint="FolderLevel"
            visibility={() =>
              F_SIGFItem?.contentFlag && [1, 2].includes(F_SIGFItem.contentFlag)
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
                qualifiers={{ guiLabel: "Use for" }}
                visibility={() => isTruthy(useAnomalous)}
              />
            </FieldRow>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 2: PARAMETERISATION
            ============================================================ */}
        <CCP4i2Tab label="Parameterisation" key="parameterisation">
          {/* B-factors */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "B-factors" }}
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
              <Typography variant="body1">Refine</Typography>
              <Box sx={{ width: "16rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="B_REFINEMENT_MODE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">B-factors</Typography>
            </Box>
          </CCP4i2ContainerElement>

          {/* Scaling */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Scaling" }}
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
              <Typography variant="body1">Use</Typography>
              <Box sx={{ width: "10rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="SCALE_TYPE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">scaling, with</Typography>
              <Box sx={{ width: "10rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="SOLVENT_MASK_TYPE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">solvent mask</Typography>
            </Box>

            <CCP4i2TaskElement
              {...props}
              itemName="SOLVENT_ADVANCED"
              qualifiers={{ guiLabel: "Use custom solvent mask parameters" }}
              visibility={() => solventMaskType === "EXPLICIT"}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Custom solvent mask parameters" }}
              containerHint="BlockLevel"
              visibility={() =>
                solventMaskType === "EXPLICIT" && isTruthy(solventAdvanced)
              }
            >
              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_VDW_RADIUS"
                qualifiers={{
                  guiLabel:
                    "Increase VDW Radius of non-ion atoms by",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_IONIC_RADIUS"
                qualifiers={{
                  guiLabel:
                    "Increase VDW Radius of potential ion atoms by",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_SHRINK"
                qualifiers={{
                  guiLabel: "Shrink the mask area by a factor of",
                }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* TLS */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Translation libration screw (TLS)" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="TLSMODE"
              qualifiers={{ guiLabel: "TLS parameters:" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="NTLSCYCLES"
              qualifiers={{ guiLabel: "Number of TLS cycles:" }}
              visibility={() => tlsMode !== "NONE"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="TLSIN"
              qualifiers={{ guiLabel: "TLS group definitions:" }}
              visibility={() => tlsMode === "FILE"}
            />

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
                qualifiers={{ guiLabel: "Reset all B-factors at start" }}
                sx={{ width: "auto" }}
              />
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="BFACSET"
                  qualifiers={{ guiLabel: "to fixed value:" }}
                  visibility={() => isTruthy(bfacSetUse)}
                />
              </Box>
            </Box>

            <CCP4i2TaskElement
              {...props}
              itemName="TLSOUT_ADDU"
              qualifiers={{
                guiLabel:
                  "Add TLS contribution to output B-factors (only for analysis and deposition)",
              }}
              visibility={() => tlsMode !== "NONE"}
            />
          </CCP4i2ContainerElement>

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
                Weight restraints versus experimental data using
              </Typography>
              <Box sx={{ width: "10rem" }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="WEIGHT_OPT"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">weighting</Typography>
            </Box>
            <CCP4i2TaskElement
              {...props}
              itemName="controlParameters.WEIGHT"
              qualifiers={{ guiLabel: "Weight:" }}
              visibility={() => weightOpt === "MANUAL"}
            />
          </CCP4i2ContainerElement>

          {/* Conformer and occupancy refinement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Conformer and occupancy refinement" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_REFINEMENT"
              qualifiers={{ guiLabel: "Refine conformer occupancies" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_GROUPS"
              qualifiers={{ guiLabel: "Specify partial occupancy groups" }}
              visibility={() => isTruthy(occupancyRefinement)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_SELECTION"
              visibility={() =>
                isTruthy(occupancyRefinement) && isTruthy(occupancyGroups)
              }
            />
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_COMPLETE"
              qualifiers={{ guiLabel: "Complete groups (sum to 1.0)" }}
              visibility={() =>
                isTruthy(occupancyRefinement) && isTruthy(occupancyGroups)
              }
            />
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_COMPLETE_TABLE"
              visibility={() =>
                isTruthy(occupancyRefinement) &&
                isTruthy(occupancyGroups) &&
                isTruthy(occupancyComplete)
              }
            />
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_INCOMPLETE"
              qualifiers={{ guiLabel: "Incomplete groups (sum to < 1.0)" }}
              visibility={() =>
                isTruthy(occupancyRefinement) && isTruthy(occupancyGroups)
              }
            />
            <CCP4i2TaskElement
              {...props}
              itemName="OCCUPANCY_INCOMPLETE_TABLE"
              visibility={() =>
                isTruthy(occupancyRefinement) &&
                isTruthy(occupancyGroups) &&
                isTruthy(occupancyIncomplete)
              }
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 3: RESTRAINTS
            ============================================================ */}
        <CCP4i2Tab label="Restraints" key="restraints">
          {/* NCS */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Non-crystallographic symmetry (NCS)",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="USE_NCS"
              qualifiers={{
                guiLabel:
                  "Use non-crystallographic symmetry (NCS) restraints",
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
              {isTruthy(useNcs) && (
                <>
                  <Typography variant="body1">Use automatic</Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="NCS_TYPE"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">NCS restraints</Typography>
                </>
              )}
            </Box>
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
                <Typography variant="body1">and max distance:</Typography>
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    {...props}
                    itemName="JELLY_DIST"
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              </Box>
            )}
          </CCP4i2ContainerElement>

          {/* ProSMART - protein */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "ProSMART - protein" }}
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
                      guiLabel: "Generate restraints for protein chain(s):",
                    }}
                    sx={{ width: "auto" }}
                  />
                  <Box sx={{ minWidth: "8rem", flex: "0 1 auto" }}>
                    <CChainSelectElement
                      job={job}
                      itemName="prosmartProtein.CHAINLIST_1"
                      options={xyzinComposition?.peptides || []}
                      label=" "
                      visibility={() => isTruthy(prosmartProteinToggle)}
                    />
                  </Box>
                  {isTruthy(prosmartProteinToggle) && (
                    <Typography variant="body1">
                      using homologous model(s):
                    </Typography>
                  )}
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
                        chain(s) from reference model(s). Minimum sequence
                        identity:
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
                      <Typography variant="body1">with weight</Typography>
                      <Box sx={{ width: "6rem" }}>
                        <CCP4i2TaskElement
                          {...props}
                          itemName="prosmartProtein.WEIGHT"
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
                      qualifiers={{ guiLabel: "Show advanced options" }}
                    />

                    {isTruthy(prosmartProteinAdvanced) && (
                      <CCP4i2ContainerElement
                        {...props}
                        itemName=""
                        qualifiers={{
                          guiLabel: "Advanced ProSMART protein options",
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
                            guiLabel: "Additional ProSMART keywords:",
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

          {/* ProSMART - nucleic acids */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "ProSMART - nucleic acids" }}
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
                        chain(s) from reference model(s). Minimum sequence
                        identity:
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
                      <Typography variant="body1">with weight</Typography>
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
                      qualifiers={{ guiLabel: "Show advanced options" }}
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

          {/* External restraints */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "External restraints" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="EXTERNAL_RESTRAINTS_FILE"
              qualifiers={{ guiLabel: "External restraints file:" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="MAKE_LINK"
              qualifiers={{
                guiLabel: "Detect covalent linkages based on coordinates",
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
              qualifiers={{ guiLabel: "Use Platonyzer metal restraints" }}
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
                    guiLabel: "Remove VDW restraints for metal sites",
                  }}
                />
              </>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 4: OUTPUT
            ============================================================ */}
        <CCP4i2Tab label="Output" key="output">
          {/* Output options */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Output options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="OUTPUT_HYDROGENS"
              qualifiers={{
                guiLabel: "Output calculated riding hydrogens to file",
              }}
            />
          </CCP4i2ContainerElement>

          {/* Map calculation */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Map calculation" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="MAP_SHARP"
              qualifiers={{
                guiLabel: "Perform map sharpening when calculating maps",
              }}
            />
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="MAP_SHARP_CUSTOM"
                qualifiers={{
                  guiLabel: "Use custom sharpening parameter (B-factor)",
                }}
                visibility={() => isTruthy(mapSharp)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="BSHARP"
                qualifiers={{ guiLabel: "B-factor" }}
                visibility={() =>
                  isTruthy(mapSharp) && isTruthy(mapSharpCustom)
                }
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          {/* Validation and analysis */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Validation and analysis" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_IRIS"
              qualifiers={{ guiLabel: "Run IRIS for validation" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_BAVERAGE"
              qualifiers={{ guiLabel: "Analyse B-factor distributions" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_RAMACHANDRAN"
              qualifiers={{ guiLabel: "Calculate Ramachandran plots" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_MOLPROBITY"
              qualifiers={{
                guiLabel: "Run MolProbity to analyse geometry",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ============================================================
            TAB 5: ADVANCED
            ============================================================ */}
        <CCP4i2Tab label="Advanced" key="advanced">
          {/* Experiment */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Experiment" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SCATTERING_FACTORS"
              qualifiers={{ guiLabel: "Diffraction experiment type:" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="SCATTERING_ELECTRON"
              qualifiers={{ guiLabel: "Form factor calculation:" }}
              visibility={() => scatteringFactors === "ELECTRON"}
            />

            {/* Neutron refinement options */}
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Neutron refinement options" }}
              containerHint="BlockLevel"
              visibility={() => scatteringFactors === "NEUTRON"}
            >
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="HYDR_USE"
                  qualifiers={{
                    guiLabel: "Use hydrogens during refinement",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="HYDR_ALL"
                  visibility={() => isTruthy(hydrUse)}
                />
              </FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="HD_INIT_TOGGLE"
                qualifiers={{
                  guiLabel: "Initialise hydrogen/deuterium fractions",
                }}
                visibility={() => isTruthy(hydrUse)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="HD_INIT"
                visibility={() =>
                  isTruthy(hydrUse) && isTruthy(hdInitToggle)
                }
              />
              <CCP4i2TaskElement
                {...props}
                itemName="HD_FRACTION"
                qualifiers={{
                  guiLabel: "Refine hydrogen/deuterium fractions",
                }}
                visibility={() => isTruthy(hydrUse)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="HD_FRACTION_TYPE"
                qualifiers={{ guiLabel: "for" }}
                visibility={() => isTruthy(hydrUse)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="H_REFINE"
                qualifiers={{ guiLabel: "Refine hydrogen positions" }}
                visibility={() => isTruthy(hydrUse)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="H_REFINE_SELECT"
                qualifiers={{ guiLabel: "for" }}
                visibility={() => isTruthy(hydrUse)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="H_TORSION"
                qualifiers={{
                  guiLabel: "Use hydrogen torsion angle restraints",
                }}
                visibility={() => isTruthy(hydrUse)}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* Resolution */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Resolution" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="RES_CUSTOM"
              qualifiers={{ guiLabel: "Use custom resolution limits" }}
            />
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="RES_MIN"
                qualifiers={{ guiLabel: "Low (dmax):" }}
                visibility={() => isTruthy(resCustom)}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="RES_MAX"
                qualifiers={{ guiLabel: "High (dmin):" }}
                visibility={() => isTruthy(resCustom)}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          {/* Reset B-factors */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reset B-factors" }}
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
                qualifiers={{ guiLabel: "Reset all B-factors at start" }}
                sx={{ width: "auto" }}
              />
              {isTruthy(bfacSetUse) && (
                <>
                  <Typography variant="body1">to fixed value:</Typography>
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
          </CCP4i2ContainerElement>

          {/* Sequence and structure */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Sequence and structure" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="ASUIN" />
          </CCP4i2ContainerElement>

          {/* Extra keywords */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Extra keywords" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="EXTRAREFMACKEYWORDS"
              qualifiers={{ guiLabel: " ", guiMode: "multiLine" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="REFMAC_KEYWORD_FILE"
            />
          </CCP4i2ContainerElement>

          {/* Cleanup */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Cleanup" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="REFMAC_CLEANUP"
              qualifiers={{
                guiLabel: "Clean up intermediate files at end of job",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
