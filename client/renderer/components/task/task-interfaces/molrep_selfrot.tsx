/**
 * Molrep — self-rotation function.
 *
 * Rewritten against the task's own def.xml. What was here before named
 * parameters that do not exist in it and never have: RESO_LOW, RESO_HIGH,
 * ANGLE_STEP, RADIUS, NPEAKS, THETA_MIN, THETA_MAX, PHI_STEP, USE_ADVANCED,
 * WAVELENGTH and a "keywords" container. `useTaskItem` returns a null item for
 * a name the container does not hold, so eleven of the interface's thirteen
 * fields rendered as nothing while looking like a finished form. It also
 * fetched an MTZ digest on every F_SIGF change to write a WAVELENGTH that has
 * no home here.
 *
 * The real vocabulary is Molrep's: NP and NPT peak counts, SEQ/SURF/ANISO/SCORE
 * mode selectors, and the B-add and B-off filter parameters whose sub-options
 * depend on HIGH_PATH_VAR. Groupings follow the Qt GUI.
 *
 * PERFORM lives in guiParameters and the plugin fixes it to 'srf' for this
 * task; Qt had its selector commented out for the same reason. The options its
 * subframes gated (PRF, SCORE, NMON_EXP, ANISO) are shown unconditionally here
 * rather than behind a control the user cannot reach.
 */
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // B-add is specified one of four ways, and each way needs a different field.
  const { value: highPathVar } = useTaskItem("HIGH_PATH_VAR");
  const { value: lowPathVar } = useTaskItem("LOW_PATH_VAR");

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Input data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Experimental data" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{ guiLabel: "Reflections" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model and sequence (optional)" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN"
              qualifiers={{ guiLabel: "Search model" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_FIX"
              qualifiers={{ guiLabel: "Fixed model" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="F_PHI_MAP"
              qualifiers={{ guiLabel: "Map coefficients" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="SEQIN"
              qualifiers={{ guiLabel: "Sequence" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="ASUIN"
              qualifiers={{ guiLabel: "Asymmetric unit contents" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="NMON"
              qualifiers={{ guiLabel: "Number of monomers to search for" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Basic options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Perform alignment and use it to rename residues and trim side chains",
            }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SEQ"
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "B-factor modification" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SURF"
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Number of peaks to analyse" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="NP"
              qualifiers={{ guiLabel: "Rotation function peaks" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="NPT"
              qualifiers={{ guiLabel: "Translation function peaks" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Advanced options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Search method" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <Typography variant="body2" sx={{ fontStyle: "italic" }}>
              SAPTF = spherically averaged phased translation function; PRF =
              phased rotation function; RF(M) uses f-obs from the density
              outside the fixed model, RF(S) from inside a sphere.
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="PRF"
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Scoring putative solutions" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" sx={{ fontStyle: "italic" }}>
              CC = correlation coefficient; PF = packing function.
            </Typography>
            <CCP4i2TaskElement
              {...props}
              itemName="SCORE"
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="NMON_EXP"
              qualifiers={{
                guiLabel:
                  "Expected number of copies (for contrast calculation only)",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="ANISO"
              qualifiers={{ guiLabel: "Scaling", guiMode: "radio" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "High pass filter (B-add, applied to input structure amplitudes)",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="HIGH_PATH_VAR"
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="SIM"
              qualifiers={{
                guiLabel:
                  "Identity between search and target sequences (0 to 1)",
              }}
              visibility={() => highPathVar === "i"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="RESMAX"
              qualifiers={{ guiLabel: "High resolution limit (Å)" }}
              visibility={() => highPathVar === "r"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="BADD"
              qualifiers={{ guiLabel: "B-add" }}
              visibility={() => highPathVar === "b"}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Low pass filter (B-off, for the removed fraction of amplitudes)",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="LOW_PATH_VAR"
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="RESMIN"
              qualifiers={{ guiLabel: "Low resolution limit (Å)" }}
              visibility={() => lowPathVar === "r"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="BOFF"
              qualifiers={{ guiLabel: "B-off" }}
              visibility={() => lowPathVar === "b"}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Space group" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement {...props} itemName="SG_OPTIONS" />
            <CCP4i2TaskElement {...props} itemName="SG" />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
