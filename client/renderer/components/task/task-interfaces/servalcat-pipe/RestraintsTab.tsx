import React from "react";
import { Box, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { CChainSelectElement } from "../../task-elements/cchainselect";
import { isTruthy } from "../../task-elements/shared-hooks";

interface RestraintsTabProps extends CCP4i2TaskInterfaceProps {
  weightOpt: any;
  useJelly: any;
  prosmartProteinToggle: any;
  prosmartProteinAdvanced: any;
  prosmartNucleicAcidToggle: any;
  prosmartNucleicAcidAdvanced: any;
  libgToggle: any;
  libgOption: any;
  libgAdvanced: any;
  platonyzerToggle: any;
  metalCoordRun: any;
  metalCoordGenOrUse: any;
  metalCoordAdvanced: any;
  hasProteinChains: boolean;
  hasNucleotideChains: boolean;
  hasMetalSites: boolean;
  xyzinComposition: any;
}

export const RestraintsTab: React.FC<RestraintsTabProps> = (props) => {
  const {
    weightOpt,
    useJelly,
    prosmartProteinToggle,
    prosmartProteinAdvanced,
    prosmartNucleicAcidToggle,
    prosmartNucleicAcidAdvanced,
    libgToggle,
    libgOption,
    libgAdvanced,
    platonyzerToggle,
    metalCoordRun,
    metalCoordGenOrUse,
    metalCoordAdvanced,
    hasProteinChains,
    hasNucleotideChains,
    hasMetalSites,
    xyzinComposition,
    ...taskProps
  } = props;

  return (
    <>
      {/* Weights */}
      <CCP4i2ContainerElement
        {...taskProps}
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
              {...taskProps}
              itemName="WEIGHT_OPT"
              qualifiers={{ guiLabel: " " }}
            />
          </Box>
          <Typography variant="body1">
            weight versus the restraints
          </Typography>
        </Box>
        <CCP4i2TaskElement
          {...taskProps}
          itemName="WEIGHT"
          qualifiers={{ guiLabel: "Weight:" }}
          visibility={() => weightOpt === "MANUAL"}
        />
        <CCP4i2TaskElement
          {...taskProps}
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
              {...taskProps}
              itemName="WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN"
              qualifiers={{ guiLabel: " " }}
            />
          </Box>
          <Box sx={{ width: "6rem" }}>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX"
              qualifiers={{ guiLabel: " " }}
            />
          </Box>
        </FieldRow>
      </CCP4i2ContainerElement>

      {/* Non-Crystallographic Symmetry (NCS) */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Non-Crystallographic Symmetry (NCS)",
        }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="USE_NCS"
          qualifiers={{
            guiLabel:
              "Use local non-crystallographic symmetry (NCS) restraints",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Covalent links */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Covalent links" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="FIND_LINKS"
          qualifiers={{
            guiLabel:
              "Detect and apply covalent linkages based on the current atomic coordinates",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Jelly-body */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Jelly-body" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
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
                  {...taskProps}
                  itemName="JELLY_SIGMA"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">
                and max distance:
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="JELLY_DIST"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
            <CCP4i2TaskElement
              {...taskProps}
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
        {...taskProps}
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
          {...taskProps}
          itemName="metalCoordPipeline.RUN_METALCOORD"
          qualifiers={{
            guiLabel:
              "Apply MetalCoord restraints for metal sites:",
          }}
        />
        {isTruthy(metalCoordRun) && (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="metalCoordPipeline.GENERATE_OR_USE"
              qualifiers={{ guiLabel: "Restraints:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="metalCoordPipeline.LINKS"
              qualifiers={{ guiLabel: "Link records:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="metalCoordPipeline.METALCOORD_RESTRAINTS"
              qualifiers={{
                guiLabel: "MetalCoord restraints file:",
              }}
              visibility={() => metalCoordGenOrUse === "USE"}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="metalCoordPipeline.TOGGLE_ADVANCED"
              qualifiers={{
                guiLabel: "Show advanced MetalCoord options",
              }}
            />
            {isTruthy(metalCoordAdvanced) && (
              <CCP4i2ContainerElement
                {...taskProps}
                itemName=""
                qualifiers={{
                  guiLabel: "Advanced MetalCoord options",
                }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="metalCoordWrapper.controlParameters.MAXIMUM_COORDINATION_NUMBER"
                  qualifiers={{
                    guiLabel:
                      "Maximum coordination number:",
                  }}
                />
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="metalCoordWrapper.controlParameters.MINIMUM_SAMPLE_SIZE"
                  qualifiers={{
                    guiLabel: "Minimum sample size:",
                  }}
                />
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="metalCoordWrapper.controlParameters.DISTANCE_THRESHOLD"
                  qualifiers={{
                    guiLabel: "Distance threshold:",
                  }}
                />
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="metalCoordWrapper.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD"
                  qualifiers={{
                    guiLabel:
                      "Procrustes distance threshold:",
                  }}
                />
                <CCP4i2TaskElement
                  {...taskProps}
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
        {...taskProps}
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
                {...taskProps}
                itemName="prosmartProtein.TOGGLE"
                qualifiers={{
                  guiLabel:
                    "Generate and apply restraints for protein chains using homologous models",
                }}
                sx={{ width: "auto" }}
              />
              <Box sx={{ minWidth: "8rem", flex: "0 1 auto" }}>
                <CChainSelectElement
                  job={taskProps.job}
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
                  {...taskProps}
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
                      {...taskProps}
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
                      {...taskProps}
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
                      {...taskProps}
                      itemName="prosmartProtein.SIDE_MAIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    atom-pairs. Interatomic distance range:
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="prosmartProtein.RMIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">to</Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
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
                      {...taskProps}
                      itemName="prosmartProtein.DMAX"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    with sigma range:
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="prosmartProtein.SGMN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">to</Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="prosmartProtein.SGMX"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    and robustness parameter (alpha)
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="prosmartProtein.ALPHA"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>

                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="prosmartProtein.ADVANCED"
                  qualifiers={{
                    guiLabel: "Show advanced options",
                  }}
                />

                {isTruthy(prosmartProteinAdvanced) && (
                  <CCP4i2ContainerElement
                    {...taskProps}
                    itemName=""
                    qualifiers={{
                      guiLabel:
                        "Advanced ProSMART protein options",
                    }}
                    containerHint="BlockLevel"
                  >
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartProtein.TOGGLE_BFAC"
                        qualifiers={{
                          guiLabel:
                            "Remove restraints for atoms with B-factor above",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartProtein.BFAC"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartProtein.TOGGLE_ALT"
                        qualifiers={{
                          guiLabel:
                            "Don't generate restraints for atoms with occupancy below",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartProtein.OCCUPANCY"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <CCP4i2TaskElement
                      {...taskProps}
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
        {...taskProps}
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
                {...taskProps}
                itemName="prosmartNucleicAcid.TOGGLE"
                qualifiers={{
                  guiLabel:
                    "Generate restraints for nucleic acid chain(s):",
                }}
                sx={{ width: "auto" }}
              />
              <Box sx={{ minWidth: "8rem", flex: "0 1 auto" }}>
                <CChainSelectElement
                  job={taskProps.job}
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
                  {...taskProps}
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
                      {...taskProps}
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
                      {...taskProps}
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
                      {...taskProps}
                      itemName="prosmartNucleicAcid.SIDE_MAIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    atom-pairs. Interatomic distance range:
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="prosmartNucleicAcid.RMIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">to</Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
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
                      {...taskProps}
                      itemName="prosmartNucleicAcid.DMAX"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    with weight
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="prosmartNucleicAcid.WEIGHT"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    and robustness parameter (alpha)
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...taskProps}
                      itemName="prosmartNucleicAcid.ALPHA"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>

                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="prosmartNucleicAcid.ADVANCED"
                  qualifiers={{
                    guiLabel: "Show advanced options",
                  }}
                />

                {isTruthy(prosmartNucleicAcidAdvanced) && (
                  <CCP4i2ContainerElement
                    {...taskProps}
                    itemName=""
                    qualifiers={{
                      guiLabel:
                        "Advanced ProSMART nucleic acid options",
                    }}
                    containerHint="BlockLevel"
                  >
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartNucleicAcid.TOGGLE_BFAC"
                        qualifiers={{
                          guiLabel:
                            "Remove restraints for atoms with B-factor above",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartNucleicAcid.BFAC"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartNucleicAcid.TOGGLE_ALT"
                        qualifiers={{
                          guiLabel:
                            "Don't generate restraints for atoms with occupancy below",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...taskProps}
                        itemName="prosmartNucleicAcid.OCCUPANCY"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <CCP4i2TaskElement
                      {...taskProps}
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
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "ADP Restraints" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="ADPR_WEIGHT"
          qualifiers={{ guiLabel: "ADP restraint weight:" }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="MAX_DIST_FOR_ADP_RESTRAINT"
          qualifiers={{
            guiLabel: "Maximum distance for ADP restraint:",
          }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="ADP_RESTRAINT_NO_LONG_RANGE"
          qualifiers={{
            guiLabel: "No long range for ADP restraint",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Platonyzer */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Platonyzer" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="platonyzer.TOGGLE"
          qualifiers={{
            guiLabel: "Use Platonyzer metal restraints",
          }}
        />
        {isTruthy(platonyzerToggle) && (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="platonyzer.MODE"
              qualifiers={{ guiLabel: "Metal type:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
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
        {...taskProps}
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
              {...taskProps}
              itemName="libg.TOGGLE"
              qualifiers={{
                guiLabel:
                  "Generate nucleic acid base restraints using libg",
              }}
            />
            {isTruthy(libgToggle) && (
              <>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="libg.OPTION"
                  qualifiers={{
                    guiLabel: "Restraint types:",
                  }}
                />
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="libg.BP"
                  qualifiers={{
                    guiLabel: "Include base pair restraints",
                  }}
                  visibility={() => libgOption === "MANUAL"}
                />
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="libg.ADVANCED"
                  qualifiers={{
                    guiLabel: "Show advanced options",
                  }}
                />
                {isTruthy(libgAdvanced) && (
                  <CCP4i2TaskElement
                    {...taskProps}
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
    </>
  );
};
