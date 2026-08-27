/**
 * ProSMART — structural comparison and restraint generation.
 *
 * The standalone wrapper, not the `prosmartProtein` / `prosmartNucleicAcid`
 * blocks that prosmart_refmac and servalcat_pipe embed. Those are compact
 * containers defined inline in each pipeline's def.xml (RMAX, SEQID, DMAX,
 * WEIGHT...); this task exposes ProSMART's own option set, which is both larger
 * and differently named (RESTRAIN_RMAX, RESTRAIN_SEQID, ...) and adds the
 * alignment, clustering and H-bond options the pipelines never surface.
 *
 * The Qt GUI drew almost all of this with `autoGenerate`, naming only the four
 * input files explicitly, and grouped the rest into five tabs by include-list.
 * Those groupings are reproduced here.
 *
 * One departure from Qt: `autoGenerate` had no way to hide a parameter that
 * does not apply, so it showed all sixty at once. The conditional visibility
 * below is new, and each condition is taken from the def.xml's own semantics
 * (CLUSTER_PERFORM gates the clustering options, H_BOND gates the H-bond ones,
 * and PROGRAM_MODE says whether alignment or restraint generation is running at
 * all) rather than invented.
 */
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { isTruthy } from "../task-elements/shared-hooks";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  const { value: programMode } = useTaskItem("PROGRAM_MODE");
  const { value: libMode } = useTaskItem("LIB_MODE");
  const { value: clusterPerform } = useTaskItem("CLUSTER_PERFORM");
  const { value: hBond } = useTaskItem("H_BOND");
  const { value: bfacFilter } = useTaskItem("RESTRAIN_BFAC_FILTER");

  // PROGRAM_MODE is one of ALIGN_RESTRAIN / ALIGN / RESTRAIN.
  const aligning = programMode !== "RESTRAIN";
  const restraining = programMode !== "ALIGN";

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Input data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Structures to compare" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="TARGET_MODEL"
              qualifiers={{ guiLabel: "Target" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="CHAINLIST_1"
              qualifiers={{ guiLabel: "Restrict to target chains" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="REFERENCE_MODELS"
              qualifiers={{ guiLabel: "Reference model(s)" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Additional input" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="EXT_FILE"
              qualifiers={{
                guiLabel: "External keywords file",
                toolTip: "File containing additional ProSMART keywords",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="FRAGLIB"
              qualifiers={{ guiLabel: "Fragment library" }}
              visibility={() => libMode !== "NONE"}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Simple options" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement {...props} itemName="PROGRAM_MODE" />
            <CCP4i2TaskElement
              {...props}
              itemName="NUCLEIC_ACID"
              qualifiers={{ guiLabel: "Structures are nucleic acid" }}
            />
            <CCP4i2TaskElement {...props} itemName="LIB_MODE" />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Additional keywords" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="KEYWORDS"
              qualifiers={{
                guiLabel: "Extra ProSMART keywords",
                guiMode: "multiLine",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Alignment">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Alignment" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
            visibility={() => aligning}
          >
            <CCP4i2TaskElement {...props} itemName="ALIGN_MODE" />
            <CCP4i2TaskElement {...props} itemName="FRAGLEN" />
            <CCP4i2TaskElement {...props} itemName="ALIGN_THRESHOLD" />
            <CCP4i2TaskElement {...props} itemName="HELIX_CUTOFF" />
            <CCP4i2TaskElement {...props} itemName="HELIX_PENALTY" />
            <CCP4i2TaskElement {...props} itemName="REWARD_SEQ" />
            <CCP4i2TaskElement {...props} itemName="ALIGN_REFINE" />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Rigid substructure">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Rigid substructure identification" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement {...props} itemName="CLUSTER_PERFORM" />
            <CCP4i2TaskElement
              {...props}
              itemName="CLUSTER_SCORE"
              visibility={() => isTruthy(clusterPerform)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="CLUSTER_ANGLE"
              visibility={() => isTruthy(clusterPerform)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="CLUSTER_MIN"
              visibility={() => isTruthy(clusterPerform)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="CLUSTER_LINK"
              visibility={() => isTruthy(clusterPerform)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="CLUSTER_RIGID"
              visibility={() => isTruthy(clusterPerform)}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Restraint generation">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Which restraints to generate" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
            visibility={() => restraining}
          >
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_MAIN_VS_SIDE" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_ALL_VS_BEST" />
            <CCP4i2TaskElement
              {...props}
              itemName="RESTRAIN_SEQID"
              qualifiers={{ guiLabel: "Minimum sequence identity (%)" }}
            />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_TO_SELF" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_SELF" />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Distance range and weighting" }}
            containerHint="FolderLevel"
            visibility={() => restraining}
          >
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_RMIN" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_RMAX" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_SIGMA" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_MIN_SIGMA" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_SIGMATYPE" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_SCALE_SIGMAS" />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Filtering" }}
            containerHint="FolderLevel"
            visibility={() => restraining}
          >
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_MAIN_CUTOFF" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_SIDE_CUTOFF" />
            <CCP4i2TaskElement
              {...props}
              itemName="RESTRAIN_OUTLIER_THRESHOLD"
            />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_BFAC_FILTER" />
            <CCP4i2TaskElement
              {...props}
              itemName="RESTRAIN_BFAC_ALPHA"
              qualifiers={{ guiLabel: "B-factor threshold (sigma)" }}
              visibility={() => isTruthy(bfacFilter)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="RESTRAIN_OCCUP"
              qualifiers={{ guiLabel: "Minimum atom occupancy" }}
            />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_ALT" />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Refmac output" }}
            containerHint="FolderLevel"
            visibility={() => restraining}
          >
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_REFMAC_TYPE" />
            <CCP4i2TaskElement {...props} itemName="RESTRAIN_RM_BONDS" />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="H-bond restraints">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "H-bond and generic restraints" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement {...props} itemName="H_BOND" />
            <CCP4i2TaskElement
              {...props}
              itemName="H_BOND_OPT"
              visibility={() => isTruthy(hBond)}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="H_BOND_OVERRIDE"
              visibility={() => isTruthy(hBond)}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Secondary structure" }}
            containerHint="FolderLevel"
            visibility={() => isTruthy(hBond)}
          >
            <CCP4i2TaskElement {...props} itemName="H_HELIX" />
            <CCP4i2TaskElement {...props} itemName="H_ALPHA" />
            <CCP4i2TaskElement {...props} itemName="H_PI" />
            <CCP4i2TaskElement {...props} itemName="H_310" />
            <CCP4i2TaskElement {...props} itemName="H_SHEET" />
            <CCP4i2TaskElement {...props} itemName="H_STRICT" />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Geometry" }}
            containerHint="FolderLevel"
            visibility={() => isTruthy(hBond)}
          >
            <CCP4i2TaskElement {...props} itemName="H_DIST" />
            <CCP4i2TaskElement {...props} itemName="H_MIN" />
            <CCP4i2TaskElement {...props} itemName="H_MAX" />
            <CCP4i2TaskElement {...props} itemName="H_MIN_SEP" />
            <CCP4i2TaskElement {...props} itemName="H_MAX_SEP" />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Scoring and output">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Scoring" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement {...props} itemName="SUPERPOSE_THRESHOLD" />
            <CCP4i2TaskElement {...props} itemName="INCLUDE_MAIN" />
            <CCP4i2TaskElement {...props} itemName="PERFORM_FLIPS" />
            <CCP4i2TaskElement {...props} itemName="DISPLAY_AS_DEGREES" />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Output" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="OUTPUT_DM" />
            <CCP4i2TaskElement
              {...props}
              itemName="OUTPUT_PDB_CHAIN_RESTRAINTS"
            />
            <CCP4i2TaskElement {...props} itemName="MERGE_CHAINS" />
            <CCP4i2TaskElement {...props} itemName="RENAME_CHAIN" />
            <CCP4i2TaskElement {...props} itemName="IS_NMR_MD_ENSEMBLE" />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
