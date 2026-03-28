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
import { Box, Typography } from "@mui/material";
import { Job } from "../../../../types/models";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { InlineField } from "../../task-elements/inline-field";
import { isTruthy } from "../../task-elements/shared-hooks";
import { CChainSelectElement } from "../../task-elements/cchainselect";

interface RestraintsTabProps extends CCP4i2TaskInterfaceProps {
  useNcs: any;
  useJelly: any;
  prosmartProteinToggle: any;
  prosmartProteinAdvanced: any;
  prosmartNucleicAcidToggle: any;
  prosmartNucleicAcidAdvanced: any;
  platonyzerToggle: any;
  hasProteinChains: boolean;
  hasNucleotideChains: boolean;
  xyzinComposition: any;
}

export const RestraintsTab: React.FC<RestraintsTabProps> = (props) => {
  const {
    useNcs,
    useJelly,
    prosmartProteinToggle,
    prosmartProteinAdvanced,
    prosmartNucleicAcidToggle,
    prosmartNucleicAcidAdvanced,
    platonyzerToggle,
    hasProteinChains,
    hasNucleotideChains,
    xyzinComposition,
    job,
    ...taskProps
  } = props;

  // Re-assemble taskProps with job for passing to elements
  const allProps = { ...taskProps, job };

  return (
    <>
      {/* NCS */}
      <CCP4i2ContainerElement
        {...allProps}
        itemName=""
        qualifiers={{
          guiLabel: "Non-crystallographic symmetry (NCS)",
        }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...allProps}
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
                  {...allProps}
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
        {...allProps}
        itemName=""
        qualifiers={{ guiLabel: "Jelly-body" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...allProps}
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
                {...allProps}
                itemName="JELLY_SIGMA"
                qualifiers={{ guiLabel: " " }}
              />
            </Box>
            <Typography variant="body1">and max distance:</Typography>
            <Box sx={{ width: "8rem" }}>
              <CCP4i2TaskElement
                {...allProps}
                itemName="JELLY_DIST"
                qualifiers={{ guiLabel: " " }}
              />
            </Box>
          </Box>
        )}
      </CCP4i2ContainerElement>

      {/* ProSMART - protein */}
      <CCP4i2ContainerElement
        {...allProps}
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
                {...allProps}
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
                  {...allProps}
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
                      {...allProps}
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
                      {...allProps}
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
                      {...allProps}
                      itemName="prosmartProtein.SIDE_MAIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    atom-pairs. Interatomic distance range:
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
                      itemName="prosmartProtein.RMIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">to</Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
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
                      {...allProps}
                      itemName="prosmartProtein.DMAX"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">with weight</Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
                      itemName="prosmartProtein.WEIGHT"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    and robustness parameter (alpha)
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
                      itemName="prosmartProtein.ALPHA"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>

                <CCP4i2TaskElement
                  {...allProps}
                  itemName="prosmartProtein.ADVANCED"
                  qualifiers={{ guiLabel: "Show advanced options" }}
                />

                {isTruthy(prosmartProteinAdvanced) && (
                  <CCP4i2ContainerElement
                    {...allProps}
                    itemName=""
                    qualifiers={{
                      guiLabel: "Advanced ProSMART protein options",
                    }}
                    containerHint="BlockLevel"
                  >
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartProtein.TOGGLE_BFAC"
                        qualifiers={{
                          guiLabel:
                            "Remove restraints for atoms with B-factor above",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartProtein.BFAC"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartProtein.TOGGLE_ALT"
                        qualifiers={{
                          guiLabel:
                            "Don't generate restraints for atoms with occupancy below",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartProtein.OCCUPANCY"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <CCP4i2TaskElement
                      {...allProps}
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
        {...allProps}
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
                {...allProps}
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
                  {...allProps}
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
                      {...allProps}
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
                      {...allProps}
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
                      {...allProps}
                      itemName="prosmartNucleicAcid.SIDE_MAIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    atom-pairs. Interatomic distance range:
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
                      itemName="prosmartNucleicAcid.RMIN"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">to</Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
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
                      {...allProps}
                      itemName="prosmartNucleicAcid.DMAX"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">with weight</Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
                      itemName="prosmartNucleicAcid.WEIGHT"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                  <Typography variant="body1">
                    and robustness parameter (alpha)
                  </Typography>
                  <Box sx={{ width: "6rem" }}>
                    <CCP4i2TaskElement
                      {...allProps}
                      itemName="prosmartNucleicAcid.ALPHA"
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>

                <CCP4i2TaskElement
                  {...allProps}
                  itemName="prosmartNucleicAcid.ADVANCED"
                  qualifiers={{ guiLabel: "Show advanced options" }}
                />

                {isTruthy(prosmartNucleicAcidAdvanced) && (
                  <CCP4i2ContainerElement
                    {...allProps}
                    itemName=""
                    qualifiers={{
                      guiLabel:
                        "Advanced ProSMART nucleic acid options",
                    }}
                    containerHint="BlockLevel"
                  >
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartNucleicAcid.TOGGLE_BFAC"
                        qualifiers={{
                          guiLabel:
                            "Remove restraints for atoms with B-factor above",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartNucleicAcid.BFAC"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <FieldRow>
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartNucleicAcid.TOGGLE_ALT"
                        qualifiers={{
                          guiLabel:
                            "Don't generate restraints for atoms with occupancy below",
                        }}
                      />
                      <CCP4i2TaskElement
                        {...allProps}
                        itemName="prosmartNucleicAcid.OCCUPANCY"
                        qualifiers={{ guiLabel: " " }}
                      />
                    </FieldRow>
                    <CCP4i2TaskElement
                      {...allProps}
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
        {...allProps}
        itemName=""
        qualifiers={{ guiLabel: "External restraints" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...allProps}
          itemName="EXTERNAL_RESTRAINTS_FILE"
          qualifiers={{ guiLabel: "External restraints file:" }}
        />
        <CCP4i2TaskElement
          {...allProps}
          itemName="MAKE_LINK"
          qualifiers={{
            guiLabel: "Detect covalent linkages based on coordinates",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Platonyzer */}
      <CCP4i2ContainerElement
        {...allProps}
        itemName=""
        qualifiers={{ guiLabel: "Platonyzer" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...allProps}
          itemName="platonyzer.TOGGLE"
          qualifiers={{ guiLabel: "Use Platonyzer metal restraints" }}
        />
        {isTruthy(platonyzerToggle) && (
          <>
            <CCP4i2TaskElement
              {...allProps}
              itemName="platonyzer.MODE"
              qualifiers={{ guiLabel: "Metal type:" }}
            />
            <CCP4i2TaskElement
              {...allProps}
              itemName="platonyzer.RM_VDW"
              qualifiers={{
                guiLabel: "Remove VDW restraints for metal sites",
              }}
            />
          </>
        )}
      </CCP4i2ContainerElement>
    </>
  );
};
