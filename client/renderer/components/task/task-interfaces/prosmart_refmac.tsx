import React, { useCallback } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";
import { useFreeRWarning } from "../../../providers/run-check-provider";

/**
 * Task interface component for Prosmart-Refmac - Prosmart-guided Refinement.
 *
 * Prosmart-Refmac is used for:
 * - Macromolecular structure refinement with external restraints
 * - Integration with Prosmart for protein restraint generation
 * - Advanced B-factor and TLS parameter refinement
 * - Comprehensive validation and geometry analysis
 * - Anomalous signal and twinning handling
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest, createPeerTask, validation } = useJob(job.id);

  // Get task items with update functions and/or values
  const { item: F_SIGFItem, value: F_SIGFValue } = useTaskItem("F_SIGF");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  // Use centralized FreeR warning hook
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

  // Get current values for visibility conditions
  const { value: hydrUseValue } = useTaskItem("HYDR_USE");
  const { value: solventMaskTypeValue } = useTaskItem("SOLVENT_MASK_TYPE");
  const { value: solventAdvancedValue } = useTaskItem("SOLVENT_ADVANCED");
  const { value: tlsModeValue } = useTaskItem("TLSMODE");
  const { value: bfacSetUseValue } = useTaskItem("BFACSETUSE");
  const { value: weightOptValue } = useTaskItem("WEIGHT_OPT");
  const { value: useNcsValue } = useTaskItem("USE_NCS");
  const { value: useJellyValue } = useTaskItem("USE_JELLY");
  const { value: mapSharpValue } = useTaskItem("MAP_SHARP");
  const { value: mapSharpCustomValue } = useTaskItem("MAP_SHARP_CUSTOM");
  const { value: scatteringFactorsValue } = useTaskItem("SCATTERING_FACTORS");
  const { value: resCustomValue } = useTaskItem("RES_CUSTOM");
  const { value: useAnomalousValue } = useTaskItem("USEANOMALOUS");

  // Handle F_SIGF file changes - extract wavelength and update anomalous/twinning flags
  const handleF_SIGFChange = useCallback(async () => {
    if (!job || job.status !== 1) return;

    // Fetch digest and extract wavelength
    if (F_SIGFItem?._objectPath) {
      const digestData = await fetchDigest(F_SIGFItem._objectPath);
      const wavelength = digestData?.wavelengths?.at(-1);
      if (wavelength && wavelength > 0 && wavelength < 9 && forceUpdateWAVELENGTH) {
        await forceUpdateWAVELENGTH(wavelength);
      }
    }

    // Handle anomalous signal and twinning based on content flag
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
  }, [F_SIGFItem?._objectPath, F_SIGFValue, job, fetchDigest, forceUpdateWAVELENGTH, forceUpdateUSEANOMALOUS, forceUpdateUSE_TWIN]);

  // Visibility conditions
  const visibility = {
    showAnomalousSignal: () =>
      F_SIGFItem?.contentFlag && [1, 2].includes(F_SIGFItem.contentFlag),
    showUseAnomalousFor: () =>
      F_SIGFItem?.contentFlag &&
      [1, 2].includes(F_SIGFItem.contentFlag) &&
      useAnomalousValue,
    showTwinRefinement: () =>
      F_SIGFItem?.contentFlag && [3].includes(F_SIGFItem.contentFlag),
    showHydrogenOptions: () => hydrUseValue,
    showSolventAdvanced: () => solventMaskTypeValue === "EXPLICIT",
    showCustomSolventParams: () =>
      solventMaskTypeValue === "EXPLICIT" && solventAdvancedValue,
    showTLSOptions: () => tlsModeValue !== "NONE",
    showTLSFile: () => tlsModeValue === "FILE",
    showBFactorReset: () => bfacSetUseValue,
    showManualWeight: () => weightOptValue === "MANUAL",
    showNCSOptions: () => useNcsValue,
    showJellyOptions: () => useJellyValue,
    showMapSharpening: () => mapSharpValue,
    showCustomBFactor: () => mapSharpValue && mapSharpCustomValue,
    showElectronFormFactor: () => scatteringFactorsValue === "ELECTRON",
    showCustomResolution: () => resCustomValue,
  };

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Input data" key="input">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Main inputs",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN"
              qualifiers={{
                guiLabel: "Coordinates",
                toolTip: "Input macromolecular coordinates for refinement",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{
                guiLabel: "Reflections",
                toolTip: "Observed reflection data for refinement",
              }}
              onChange={handleF_SIGFChange}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="WAVELENGTH"
              qualifiers={{
                guiLabel: "Wavelength",
                toolTip: "Wavelength",
              }}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Anomalous signal",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
              visibility={visibility.showAnomalousSignal}
            >
              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="USEANOMALOUS"
                  qualifiers={{
                    guiLabel: "Use anomalous",
                    toolTip: "Use anomalous differences in refinement",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="USEANOMALOUSFOR"
                  qualifiers={{
                    guiLabel: "Use for",
                    toolTip: "How to use anomalous differences",
                  }}
                  visibility={visibility.showUseAnomalousFor}
                />
              </FieldRow>
            </CCP4i2ContainerElement>

            <CCP4i2TaskElement
              {...props}
              itemName="USE_TWIN"
              qualifiers={{
                guiLabel: "Twin refinement",
                toolTip: "Handle crystal twinning during refinement",
              }}
              visibility={visibility.showTwinRefinement}
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
              itemName="DICT_LIST"
              qualifiers={{
                guiLabel: "Additional dictionaries",
                toolTip:
                  "Additional geometry dictionaries for non-standard residues",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Options",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="NCYCLES"
              qualifiers={{
                guiLabel: "Cycles",
                toolTip: "Number of refinement cycles",
              }}
            />

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="HYDR_USE"
                qualifiers={{
                  guiLabel: "Use hydrogens during refinement",
                  toolTip: "Include hydrogen atoms in refinement calculations",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="HYDR_ALL"
                qualifiers={{
                  guiLabel: "Hydrogen type",
                  toolTip: "Use all hydrogen atoms",
                }}
                visibility={visibility.showHydrogenOptions}
              />
            </FieldRow>

            <CCP4i2TaskElement
              {...props}
              itemName="ADD_WATERS"
              qualifiers={{
                guiLabel: "Add waters",
                toolTip: "Automatically add water molecules during refinement",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="USE_TWIN"
              qualifiers={{
                guiLabel: "Crystal is twinned",
                toolTip: "Handle crystal twinning",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Parameterisation" key="parameterisation">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "B-factors",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="B_REFINEMENT_MODE"
              qualifiers={{
                guiLabel: "B-factors",
                toolTip: "B-factor refinement mode",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Scaling",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="SCALE_TYPE"
                qualifiers={{
                  guiLabel: "Use ",
                  toolTip: "Scaling type for structure factors",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_MASK_TYPE"
                qualifiers={{
                  guiLabel: "solvent scaling, with mask type",
                  toolTip: "Type of solvent mask for bulk solvent correction",
                }}
              />
            </FieldRow>

            <CCP4i2TaskElement
              {...props}
              itemName="SOLVENT_ADVANCED"
              qualifiers={{
                guiLabel: "Use custom solvent mask parameters",
                toolTip: "Enable custom parameters for solvent masking",
              }}
              visibility={visibility.showSolventAdvanced}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Custom parameters",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
              visibility={visibility.showCustomSolventParams}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_VDW_RADIUS"
                qualifiers={{
                  guiLabel: "Increase VDW Radius of non-ion atoms by ",
                  toolTip:
                    "Additional radius for non-ionic atoms in solvent mask",
                }}
              />

              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_IONIC_RADIUS"
                qualifiers={{
                  guiLabel: "Increase VDW Radius of potential ion atoms by ",
                  toolTip: "Additional radius for ionic atoms in solvent mask",
                }}
              />

              <CCP4i2TaskElement
                {...props}
                itemName="SOLVENT_SHRINK"
                qualifiers={{
                  guiLabel: "Shrink the mask area by a factor of",
                  toolTip: "Shrinkage factor for solvent mask",
                }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Translation libration screw (TLS)",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="TLSMODE"
              qualifiers={{
                guiLabel: "TLS parameters",
                toolTip: "TLS parameter refinement mode",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="NTLSCYCLES"
              qualifiers={{
                guiLabel: "Number of TLS cycles",
                toolTip: "Number of TLS refinement cycles",
              }}
              visibility={visibility.showTLSOptions}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Custom parameters",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
              visibility={visibility.showTLSOptions}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="TLSIN"
                qualifiers={{
                  guiLabel: "TLS coefficients",
                  toolTip: "Input file containing TLS parameters",
                }}
                visibility={visibility.showTLSFile}
              />

              <FieldRow>
                <CCP4i2TaskElement
                  {...props}
                  itemName="BFACSETUSE"
                  qualifiers={{
                    guiLabel: "Reset all B-factors at start ",
                    toolTip:
                      "Reset all B-factors to a fixed value before TLS refinement",
                  }}
                />
                <CCP4i2TaskElement
                  {...props}
                  itemName="BFACSET"
                  qualifiers={{
                    guiLabel: "...to a value of",
                    toolTip: "B-factor value for reset",
                  }}
                  visibility={visibility.showBFactorReset}
                />
              </FieldRow>

              <CCP4i2TaskElement
                {...props}
                itemName="TLSOUT_ADDU"
                qualifiers={{
                  guiLabel:
                    "Add TLS contribution to output B-factors (only for analysis and deposition)",
                  toolTip: "Include TLS contribution in final B-factors",
                }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Restraints" key="restraints">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Weights",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="WEIGHT_OPT"
                qualifiers={{
                  guiLabel: "Weight restraints versus experimental data using",
                  toolTip:
                    "Method for weighting geometric restraints against experimental data",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="controlParameters.WEIGHT"
                qualifiers={{
                  guiLabel: "Weight",
                  toolTip: "Manual weight value",
                }}
                visibility={visibility.showManualWeight}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Non-crystallographic symmetry (NCS)",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="USE_NCS"
              qualifiers={{
                guiLabel: "Use NCS",
                toolTip: "Apply non-crystallographic symmetry restraints",
              }}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "NCS",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
              visibility={visibility.showNCSOptions}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="NCS_TYPE"
                qualifiers={{
                  guiLabel: "Type",
                  toolTip: "Type of NCS restraints",
                }}
              />

              <CCP4i2TaskElement
                {...props}
                itemName="NCS_AUTO"
                qualifiers={{
                  guiLabel: "Auto",
                  toolTip: "Automatically detect NCS relationships",
                }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Jelly-body",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="USE_JELLY"
              qualifiers={{
                guiLabel: "Use jelly body",
                toolTip: "Apply jelly-body restraints for flexible refinement",
              }}
            />

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Jelly body",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
              visibility={visibility.showJellyOptions}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="JELLY_SIGMA"
                qualifiers={{
                  guiLabel: "Sigma",
                  toolTip: "Sigma value for jelly-body restraints",
                }}
              />

              <CCP4i2TaskElement
                {...props}
                itemName="JELLY_DIST"
                qualifiers={{
                  guiLabel: "Dist",
                  toolTip: "Distance cutoff for jelly-body restraints",
                }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName="prosmartProtein"
            qualifiers={{
              guiLabel: "Prosmart - protein",
              toolTip: "Prosmart restraints for protein chains",
            }}
            containerHint="FolderLevel"
          />
        </CCP4i2Tab>

        <CCP4i2Tab label="Output" key="output">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Output options",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="OUTPUT_HYDROGENS"
              qualifiers={{
                guiLabel: "Output calculated riding hydrogens to file",
                toolTip:
                  "Include calculated hydrogen atoms in output coordinates",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Map calculation",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="MAP_SHARP"
              qualifiers={{
                guiLabel: "Perform map sharpening when calculating maps",
                toolTip: "Apply B-factor sharpening to improve map quality",
              }}
            />

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="MAP_SHARP_CUSTOM"
                qualifiers={{
                  guiLabel: "Use custom sharpening parameter (B-factor)",
                  toolTip: "Specify custom B-factor for map sharpening",
                }}
                visibility={visibility.showMapSharpening}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="BSHARP"
                qualifiers={{
                  guiLabel: "B-factor",
                  toolTip: "B-factor value for custom map sharpening",
                }}
                visibility={visibility.showCustomBFactor}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Validation and analysis",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_BAVERAGE"
              qualifiers={{
                guiLabel: "Analyse B-factor distributions",
                toolTip: "Generate B-factor distribution analysis",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_RAMACHANDRAN"
              qualifiers={{
                guiLabel: "Calculate Ramachandran plots",
                toolTip: "Generate Ramachandran plot validation",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="VALIDATE_MOLPROBITY"
              qualifiers={{
                guiLabel: "Run MolProbity to analyse geometry",
                toolTip:
                  "Perform comprehensive geometry validation with MolProbity",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Advanced" key="advanced">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Experiment",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SCATTERING_FACTORS"
              qualifiers={{
                guiLabel: "Diffraction experiment type",
                toolTip: "Type of diffraction experiment (X-ray or electron)",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="SCATTERING_ELECTRON"
              qualifiers={{
                guiLabel: "Form factor calculation",
                toolTip:
                  "Method for electron scattering form factor calculation",
              }}
              visibility={visibility.showElectronFormFactor}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Resolution",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="RES_CUSTOM"
              qualifiers={{
                guiLabel: "Use custom resolution",
                toolTip: "Override automatic resolution limits",
              }}
            />

            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="RES_MIN"
                qualifiers={{
                  guiLabel: "min",
                  toolTip: "Minimum resolution limit",
                }}
                visibility={visibility.showCustomResolution}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="RES_MAX"
                qualifiers={{
                  guiLabel: "max",
                  toolTip: "Maximum resolution limit",
                }}
                visibility={visibility.showCustomResolution}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Reset B-factors",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="BFACSETUSE"
              qualifiers={{
                guiLabel: "Reset all B-factors at start ",
                toolTip:
                  "Reset all B-factors to a fixed value before refinement",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="BFACSET"
              qualifiers={{
                guiLabel: "...to a value of",
                toolTip: "B-factor value for reset",
              }}
              visibility={visibility.showBFactorReset}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Extra keywords",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="EXTRAREFMACKEYWORDS"
              qualifiers={{
                guiLabel: " ",
                guiMode: "multiLine",
                toolTip: "Additional Refmac keywords for advanced control",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="REFMAC_KEYWORD_FILE"
              qualifiers={{
                guiLabel: "",
                toolTip: "File containing additional Refmac keywords",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
