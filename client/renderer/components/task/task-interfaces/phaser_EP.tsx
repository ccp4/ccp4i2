/*
 * Copyright (C) 2025-2026 Newcastle University
 * Copyright (C) 2025-2026 University of York
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
import React, { useCallback } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";

/**
 * Task interface component for Phaser Experimental Phasing (EP).
 *
 * Provides functionality for:
 * - Heavy atom location using various methods (search, manual input, partial model)
 * - Phase calculation and density modification
 * - Integration with Parrot and ModelCraft for automated building
 * - Automatic wavelength detection from reflection files
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest } = useJob(job.id);

  // Get task items at top level (Rules of Hooks)
  const { item: F_SIGFItem } = useTaskItem("F_SIGF");
  const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");
  const { value: PARTIALMODELORMAP_value } = useTaskItem("PARTIALMODELORMAP");
  const { value: COMP_BY_value } = useTaskItem("COMP_BY");
  const { value: RUNMODELCRAFT_value } = useTaskItem("RUNMODELCRAFT");

  // Handle F_SIGF file change - extract wavelength from digest
  const handleF_SIGFChange = useCallback(async () => {
    if (!F_SIGFItem?._objectPath) return;

    const digestData = await fetchDigest(F_SIGFItem._objectPath);
    const wavelength = digestData?.wavelengths?.at(-1);
    if (wavelength && wavelength > 0 && wavelength < 9) {
      await forceUpdateWAVELENGTH(wavelength);
    }
  }, [F_SIGFItem?._objectPath, fetchDigest, forceUpdateWAVELENGTH]);

  // Visibility conditions
  const visibility = {
    isPartialModel: () => PARTIALMODELORMAP_value === "MODEL",
    isPartialMap: () => PARTIALMODELORMAP_value === "MAP",
    isNoPartial: () => PARTIALMODELORMAP_value === "NONE",
    isSearch: () => PARTIALMODELORMAP_value === "SEARCH",
    showASUFile: () => COMP_BY_value === "ASU",
    showMolecularWeights: () => COMP_BY_value === "MW",
    hasModelCraft: () => RUNMODELCRAFT_value === true,
  };

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Input data" key="input">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Observed data",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{
                guiLabel: "Reflections",
                toolTip: "Anomalous reflection data",
              }}
              onChange={handleF_SIGFChange}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="FREERFLAG"
              qualifiers={{
                guiLabel: "Free R flags",
                toolTip: "Test set flags for cross-validation",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Heavy atom input",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="PARTIALMODELORMAP"
              qualifiers={{
                guiLabel: "How to provide heavy atoms",
                toolTip:
                  "Choose whether to search for heavy atoms, provide coordinates, or use a partial model",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_HA"
              qualifiers={{
                guiLabel: "Heavy atom coordinates",
                toolTip: "PDB file containing heavy atom positions",
              }}
              visibility={visibility.isNoPartial}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_PARTIAL"
              qualifiers={{
                guiLabel: "Partial model coordinates",
                toolTip:
                  "Coordinates of a partial protein/nucleic acid model",
              }}
              visibility={visibility.isPartialModel}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="MAPCOEFF_PARTIAL"
              qualifiers={{
                guiLabel: "Map coefficients",
                toolTip:
                  "Map coefficients for partial protein/nucleic acid model",
              }}
              visibility={visibility.isPartialMap}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Heavy atom search",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
            visibility={visibility.isSearch}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SFAC"
              qualifiers={{
                guiLabel: "Atom type to find",
                toolTip: "Element symbol for the heavy atom (e.g. SE, S, BR)",
              }}
            />

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="FIND"
                qualifiers={{
                  guiLabel: "Number to find",
                  toolTip: "Number of heavy atom sites to search for",
                }}
              />

              <CCP4i2TaskElement
                {...props}
                itemName="NTRY"
                qualifiers={{
                  guiLabel: "SHELXD trials",
                  toolTip: "Number of SHELXD trials for heavy atom search",
                }}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Scattering in the crystal",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="COMP_BY"
              qualifiers={{
                guiLabel: "How to specify scattering content",
                toolTip:
                  "Method for defining the scattering content of the asymmetric unit",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="ASUFILE"
              qualifiers={{
                guiLabel: "CCP4i2 ASU file",
                toolTip: "Asymmetric unit content file",
              }}
              visibility={visibility.showASUFile}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="ASU_PROTEIN_MW"
              qualifiers={{
                guiLabel: "Protein (Da)",
                toolTip:
                  "Molecular weight of protein in the asymmetric unit",
              }}
              visibility={visibility.showMolecularWeights}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="ASU_NUCLEICACID_MW"
              qualifiers={{
                guiLabel: "Nucleic acid (Da)",
                toolTip:
                  "Molecular weight of nucleic acid in the asymmetric unit",
              }}
              visibility={visibility.showMolecularWeights}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Parameters" key="params">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Basic controls",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="WAVELENGTH"
              qualifiers={{
                guiLabel: "Wavelength",
                toolTip: "X-ray wavelength in Angstroms",
              }}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Resolution",
                initiallyOpen: true,
              }}
              containerHint="RowLevel"
            >
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="RESOLUTION_LOW"
                  qualifiers={{
                    guiLabel: "Low resolution limit",
                    toolTip: "Low resolution cutoff in Angstroms",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="RESOLUTION_HIGH"
                  qualifiers={{
                    guiLabel: "High resolution limit",
                    toolTip: "High resolution cutoff in Angstroms",
                  }}
                />
              </FieldRow>
            </CCP4i2ContainerElement>

            <CCP4i2TaskElement
              {...props}
              itemName="ELEMENTS"
              qualifiers={{
                guiLabel: "Elements for HA completion",
                toolTip:
                  "Additional element types to search for during heavy atom completion",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Extra steps",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="RUNPARROT"
                qualifiers={{
                  guiLabel: "Run Parrot",
                  toolTip: "Run Parrot density modification",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="RUNMODELCRAFT"
                qualifiers={{
                  guiLabel: "Run ModelCraft",
                  toolTip: "Run automated model building with ModelCraft",
                }}
              />
            </FieldRow>

            <CCP4i2TaskElement
              {...props}
              itemName="MODELCRAFT_ITERATIONS"
              qualifiers={{
                guiLabel: "ModelCraft cycles",
                toolTip: "Number of ModelCraft iterations",
              }}
              visibility={visibility.hasModelCraft}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Similarity of partial model",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
            visibility={visibility.isPartialModel}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="PART_VARI"
              qualifiers={{
                guiLabel: "How to specify similarity",
                toolTip: "Choose between sequence identity or coordinate RMSD",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="PART_DEVI"
              qualifiers={{
                guiLabel: "Sequence identity (0.0-1.0) or RMSD (Angstroms)",
                toolTip:
                  "Expected similarity between partial model and target",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Keywords" key="keywords">
          <CCP4i2TaskElement
            {...props}
            itemName="keywords"
            qualifiers={{
              guiLabel: "Additional keywords",
              guiMode: "multiLine",
              toolTip: "Additional Phaser keywords for advanced control",
            }}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
