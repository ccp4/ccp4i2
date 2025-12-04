import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: ADD_WATERS } = useTaskItem("ADD_WATERS");
  const { value: BFACSETUSE } = useTaskItem("BFACSETUSE");
  const { value: DATA_METHOD } = useTaskItem("DATA_METHOD");
  const { value: HKLIN_IS_I_SIGI } = useTaskItem("HKLIN_IS_I_SIGI");
  const { value: HYDR_USE } = useTaskItem("HYDR_USE");
  const { value: MERGED_OR_UNMERGED } = useTaskItem("MERGED_OR_UNMERGED");
  const { value: RANDOMIZEUSE } = useTaskItem("RANDOMIZEUSE");
  const { value: RES_CUSTOM } = useTaskItem("RES_CUSTOM");
  const { value: RUN_ADP_ANALYSIS } = useTaskItem("RUN_ADP_ANALYSIS");
  const { value: RUN_COORDADPDEV_ANALYSIS } = useTaskItem("RUN_COORDADPDEV_ANALYSIS");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Main inputs
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          </CCP4i2ContainerElement>
          {(DATA_METHOD === "xtal") && (
            <>
              {(MERGED_OR_UNMERGED === "merged") && (
                <CCP4i2TaskElement itemName="HKLIN" {...props} />
              )}
              {(MERGED_OR_UNMERGED === "unmerged") && (
                <CCP4i2TaskElement itemName="HKLIN_UNMERGED" {...props} />
              )}
              <CCP4i2TaskElement itemName="MERGED_OR_UNMERGED" {...props} qualifiers={{ guiLabel: "Diffraction data are" }} />
              <CCP4i2TaskElement itemName="F_SIGF_OR_I_SIGI" {...props} qualifiers={{ guiLabel: "Refinement against" }} />
              <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
              <CCP4i2TaskElement itemName="USE_TWIN" {...props} qualifiers={{ guiLabel: "Twin refinement" }} />
            </>
          )}
          {(DATA_METHOD === "spa") && (
            <>
              <CCP4i2TaskElement itemName="MAPIN1" {...props} qualifiers={{ guiLabel: "Half map 1" }} />
              <CCP4i2TaskElement itemName="MAPIN2" {...props} qualifiers={{ guiLabel: "Half map 2" }} />
              <CCP4i2TaskElement itemName="MAPMASK" {...props} qualifiers={{ guiLabel: "Mask" }} />
              <CCP4i2TaskElement itemName="RES_MIN" {...props} qualifiers={{ guiLabel: "Resolution:" }} />
              <CCP4i2TaskElement itemName="MASK_RADIUS" {...props} qualifiers={{ guiLabel: "Mask radius:" }} />
            </>
          )}
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Additional geometry dictionaries
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="DICT_LIST" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="NCYCLES" {...props} qualifiers={{ guiLabel: "Number of refinement cycles:" }} />
            {(HYDR_USE === false) && (
              <CCP4i2TaskElement itemName="HYDR_USE" {...props} qualifiers={{ guiLabel: "Use riding hydrogens during refinement" }} />
            )}
            {(HYDR_USE === true) && (
              <CCP4i2TaskElement itemName="HYDR_ALL" {...props} qualifiers={{ guiLabel: "Use riding hydrogens during refinement" }} />
            )}
            {(DATA_METHOD === "xtal") && (
              <CCP4i2TaskElement itemName="ADD_WATERS" {...props} qualifiers={{ guiLabel: "Add waters" }} />
            )}
            {(ADD_WATERS === true) && (
              <CCP4i2TaskElement itemName="NCYCLES_AFTER_ADD_WATERS" {...props} qualifiers={{ guiLabel: "refinement cycles" }} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Parameterisation">
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Restraints">
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced">
          {(DATA_METHOD === "xtal") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="RES_CUSTOM" {...props} qualifiers={{ guiLabel: "Use custom resolution limits" }} />
              {(RES_CUSTOM === true) && (
                <CCP4i2TaskElement itemName="RES_MAX" {...props} qualifiers={{ guiLabel: "lowest (dmax):" }} />
              )}
              <CCP4i2TaskElement itemName="FREERFLAG_NUMBER" {...props} qualifiers={{ guiLabel: "FreeR flag number for test set:" }} />
              <CCP4i2TaskElement itemName="SCATTERING_FACTORS" {...props} qualifiers={{ guiLabel: "Diffraction experiment type:" }} />
              <CCP4i2TaskElement itemName="USE_WORK_IN_EST" {...props} qualifiers={{ guiLabel: "Use work reflections in maximum likelihood parameter estimates" }} />
            </CCP4i2ContainerElement>
          )}
          {(DATA_METHOD === "spa") && (
            <CCP4i2TaskElement itemName="CROSS_VALIDATION" {...props} qualifiers={{ guiLabel: "Cross validation with half maps" }} />
          )}
          <CCP4i2TaskElement itemName="H_OUT" {...props} qualifiers={{ guiLabel: "Write hydrogen atoms in the output model" }} />
          {(HYDR_USE === true) && (
            <CCP4i2TaskElement itemName="H_REFINE" {...props} qualifiers={{ guiLabel: "Refine hydrogen positions" }} />
          )}
          <CCP4i2TaskElement itemName="KEEP_CHARGES" {...props} qualifiers={{ guiLabel: "Keep charges, i.e. use scattering factor for charged atoms where relevant" }} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Structure model modification before refinement
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="BFACSETUSE" {...props} qualifiers={{ guiLabel: "Reset all ADPs at start" }} />
            {(BFACSETUSE === true) && (
              <CCP4i2TaskElement itemName="BFACSET" {...props} qualifiers={{ guiLabel: "&nbsp;to fixed value:" }} />
            )}
            <CCP4i2TaskElement itemName="RANDOMIZEUSE" {...props} qualifiers={{ guiLabel: "Shake coordinates at start" }} />
            {(RANDOMIZEUSE === true) && (
              <CCP4i2TaskElement itemName="RANDOMIZE" {...props} qualifiers={{ guiLabel: "&nbsp;with specified RMSD:" }} />
            )}
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Additional keywords
          </Typography>
          <CCP4i2TaskElement itemName="SERVALCAT_KEYWORD_FILE" {...props} />
          <CCP4i2TaskElement itemName="EXTRA_SERVALCAT_OPTIONS" {...props} qualifiers={{ guiLabel: "Extra servalcat command line options:" }} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Validation and Analysis
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="VALIDATE_IRIS" {...props} qualifiers={{ guiLabel: "Generate Iris validation report" }} />
            <CCP4i2TaskElement itemName="VALIDATE_RAMACHANDRAN" {...props} qualifiers={{ guiLabel: "Generate Ramachandran plots" }} />
            <CCP4i2TaskElement itemName="VALIDATE_MOLPROBITY" {...props} qualifiers={{ guiLabel: "Run MolProbity to analyse geometry" }} />
            <CCP4i2TaskElement itemName="RUN_ADP_ANALYSIS" {...props} qualifiers={{ guiLabel: "Run ADP analysis" }} />
            {(RUN_ADP_ANALYSIS === true) && (
              <CCP4i2TaskElement itemName="ADP_IQR_FACTOR" {...props} qualifiers={{ guiLabel: "Atoms with a B-value lower than the first quartile - factor * interquartile_rangeor higher than the third quartile + factor * interquartile_range to be reported. Factor:" }} />
            )}
            <CCP4i2TaskElement itemName="RUN_COORDADPDEV_ANALYSIS" {...props} qualifiers={{ guiLabel: "Run analysis of changes in coordinates and ADPs" }} />
            {(RUN_COORDADPDEV_ANALYSIS === true) && (
              <CCP4i2TaskElement itemName="monitor.MIN_COORDDEV" {...props} qualifiers={{ guiLabel: "Minimum shift of atom coordinates to be reported:" }} />
            )}
            {(RUN_COORDADPDEV_ANALYSIS === true) && (
              <CCP4i2TaskElement itemName="monitor.MIN_ADPDEV" {...props} qualifiers={{ guiLabel: "Minimum shift of B-values to be reported:" }} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;