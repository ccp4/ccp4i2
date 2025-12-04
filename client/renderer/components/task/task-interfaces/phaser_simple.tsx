import React, { useCallback, useRef, useEffect } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";

/**
 * Task interface component for Phaser Simple - Molecular Replacement.
 *
 * Phaser Simple is used for:
 * - Automated molecular replacement using known search models
 * - Finding multiple copies of molecules in the asymmetric unit
 * - Integration with refinement and model building pipelines
 * - Handling both structure factors and intensities
 * - Spacegroup testing and validation
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Use refs to track processed values and prevent cycles
  const lastProcessedF_SIGFValue = useRef<any>(null);
  const initializationDone = useRef(false);
  const currentJobId = useRef<number | null>(null);

  // Get task items
  const { item: F_SIGFItem } = useTaskItem("F_SIGF");
  const { update: updateF_OR_I } = useTaskItem("F_OR_I");
  const { value: INPUT_FIXED_value } = useTaskItem("INPUT_FIXED");
  const { value: COMP_BY_value } = useTaskItem("COMP_BY");
  const { value: SGALT_SELECT_value } = useTaskItem("SGALT_SELECT");
  const { value: ID_RMS_value } = useTaskItem("ID_RMS");

  // Stable initialization function (runs once per job)
  const handleInitialization = useCallback(async () => {
    if (
      initializationDone.current ||
      !F_SIGFItem ||
      !updateF_OR_I ||
      job?.status !== 1
    ) {
      return;
    }

    const currentFlag = F_SIGFItem.contentFlag;
    if ([2, 4].includes(currentFlag)) {
      try {
        await updateF_OR_I("F");
        initializationDone.current = true;
      } catch (error) {
        console.error("Error during initialization:", error);
      }
    }
  }, [F_SIGFItem, updateF_OR_I, job?.status]);

  // Reset initialization when job changes
  useEffect(() => {
    if (currentJobId.current !== job?.id) {
      initializationDone.current = false;
      lastProcessedF_SIGFValue.current = null;
      currentJobId.current = job?.id || null;
    }
  }, [job?.id]);

  // Run initialization once when component mounts or job changes
  useEffect(() => {
    if (!initializationDone.current && job?.id) {
      handleInitialization();
    }
  }, [handleInitialization, job?.id]);

  // Stable change handler with cycle prevention
  const handleF_SIGFChange = useCallback(async () => {
    if (!F_SIGFItem || !updateF_OR_I || job?.status !== 1) return;

    // Prevent processing the same value multiple times
    const currentValue = F_SIGFItem.contentFlag;
    if (lastProcessedF_SIGFValue.current === currentValue) return;

    if ([2, 4].includes(currentValue)) {
      try {
        lastProcessedF_SIGFValue.current = currentValue;
        await updateF_OR_I("F");
      } catch (error) {
        console.error("Error updating F_OR_I:", error);
      }
    }
  }, [F_SIGFItem, updateF_OR_I, job?.status]);

  // Visibility conditions (stable references)
  const visibility = {
    showF_OR_I: () =>
      F_SIGFItem?.contentFlag && [1, 3].includes(F_SIGFItem.contentFlag),
    showXYZIN_FIXED: () => INPUT_FIXED_value === true,
    showASUFile: () => COMP_BY_value === "ASU",
    showMolecularWeights: () => COMP_BY_value === "MW",
    showSpacegroupList: () => SGALT_SELECT_value === "LIST",
    showSequenceIdentity: () => ID_RMS_value === "ID",
    showRMSD: () => ID_RMS_value === "RMS",
  };

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Main inputs" key="main">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Key files",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="F_OR_I"
              qualifiers={{
                guiLabel: "Use Fs or Is",
                toolTip:
                  "Choose between structure factors (F) or intensities (I)",
              }}
              visibility={visibility.showF_OR_I}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{
                guiLabel: "Reflections",
                toolTip: "Reflection data file",
              }}
              onChange={handleF_SIGFChange} // Event-driven updates
            />

            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN"
              qualifiers={{
                guiLabel: "Search coordinates",
                toolTip: "Search model coordinates for molecular replacement",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="INPUT_FIXED"
              qualifiers={{
                guiLabel: "Have known partial model",
                toolTip:
                  "Include a known partial model that will remain fixed during the search",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_FIXED"
              qualifiers={{
                guiLabel: "Known partial model",
                toolTip: "Coordinates of the fixed partial model",
              }}
              visibility={visibility.showXYZIN_FIXED}
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
              guiLabel: "Basic parameters",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="NCOPIES"
              qualifiers={{
                guiLabel: "Copies to find",
                toolTip:
                  "Number of copies of the search model to find in the asymmetric unit",
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

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Extra steps",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
            >
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="RUNSHEETBEND"
                  qualifiers={{
                    guiLabel: "Sheet bend",
                    toolTip: "Run sheet bend correction",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="RUNREFMAC"
                  qualifiers={{
                    guiLabel: "Refmac",
                    toolTip:
                      "Run Refmac refinement after molecular replacement",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="RUNCOOT"
                  qualifiers={{
                    guiLabel: "Coot add-waters",
                    toolTip: "Run Coot to add water molecules",
                  }}
                />
              </FieldRow>
            </CCP4i2ContainerElement>
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
              itemName="ASU_NUCLEICACID_MW"
              qualifiers={{
                guiLabel: "Nucleic acid (Da)",
                toolTip:
                  "Molecular weight of nucleic acid in the asymmetric unit",
              }}
              visibility={visibility.showMolecularWeights}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="ASU_PROTEIN_MW"
              qualifiers={{
                guiLabel: "Protein (Da)",
                toolTip: "Molecular weight of protein in the asymmetric unit",
              }}
              visibility={visibility.showMolecularWeights}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Spacegroups",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SGALT_SELECT"
              qualifiers={{
                guiLabel: "How spacegroups are chosen",
                toolTip: "Method for selecting spacegroups to test",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="SGALT_TEST"
              qualifiers={{
                guiLabel: "List of spacegroups to try",
                toolTip:
                  "Specific spacegroups to test during molecular replacement",
              }}
              visibility={visibility.showSpacegroupList}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Similarity of search model",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="ID_RMS"
              qualifiers={{
                guiLabel: "How to specify similarity",
                toolTip: "Choose between sequence identity or coordinate RMSD",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="SEARCHSEQUENCEIDENTITY"
              qualifiers={{
                guiLabel: "Sequence identity (0.0-1.0)",
                toolTip: "Sequence identity between search model and target",
              }}
              visibility={visibility.showSequenceIdentity}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="SEARCHRMS"
              qualifiers={{
                guiLabel: "Expected coordinate RMSD (Angstroms)",
                toolTip:
                  "Expected RMSD between search model and target coordinates",
              }}
              visibility={visibility.showRMSD}
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
