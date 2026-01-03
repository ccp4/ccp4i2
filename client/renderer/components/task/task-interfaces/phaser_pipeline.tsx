import React, { useCallback, useEffect, useRef } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";

/**
 * Task interface component for Phaser Pipeline - Automated Molecular Replacement Pipeline.
 *
 * Phaser Pipeline is used for:
 * - Automated molecular replacement with ensemble search models
 * - Multi-step phasing and model building workflow
 * - Integration with refinement and density modification
 * - Automated space group testing and validation
 * - Post-MR model building with sheet bend correction
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Use refs to track processed states and prevent cycles
  const initializationDone = useRef(false);
  const currentJobId = useRef<number | null>(null);
  const lastProcessedF_SIGFValue = useRef<any>(null);

  // Get task items
  const { item: F_SIGFItem, value: F_SIGFValue } = useTaskItem("F_SIGF");
  const { update: updateF_OR_I } = useTaskItem("F_OR_I");
  const { value: COMP_BY_value } = useTaskItem("COMP_BY");
  const { value: SGALT_SELECT_value } = useTaskItem("SGALT_SELECT");

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

  // Stable change handler with cycle prevention
  const handleF_SIGFChange = useCallback(async () => {
    if (!F_SIGFItem || !updateF_OR_I || job?.status !== 1) return;

    // Prevent processing the same value multiple times
    const currentFlag = F_SIGFItem.contentFlag;
    if (lastProcessedF_SIGFValue.current === currentFlag) return;

    if ([2, 4].includes(currentFlag)) {
      try {
        lastProcessedF_SIGFValue.current = currentFlag;
        await updateF_OR_I("F");
      } catch (error) {
        console.error("Error updating F_OR_I:", error);
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

  // Visibility conditions (stable references)
  const visibility = {
    showF_OR_I: () =>
      F_SIGFItem?.contentFlag && [1, 3].includes(F_SIGFItem.contentFlag),
    showASUFile: () => COMP_BY_value === "ASU",
    showMolecularWeights: () => COMP_BY_value === "MW",
    showSpacegroupList: () => SGALT_SELECT_value === "LIST",
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
                toolTip:
                  "Reflection data file containing observed structure factors or intensities",
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

            <CCP4i2TaskElement
              {...props}
              itemName="ENSEMBLES"
              qualifiers={{
                guiLabel: "Ensembles",
                toolTip: "Ensemble search models for molecular replacement",
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
                    toolTip:
                      "Run sheet bend correction for beta-strand alignment",
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
                    toolTip: "Run Coot to automatically add water molecules",
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
