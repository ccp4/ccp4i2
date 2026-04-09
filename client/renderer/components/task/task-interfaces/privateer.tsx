import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const newSugar = useBoolToggle(useTaskItem, "NEW_SUGAR");
  const glyTouCan = useBoolToggle(useTaskItem, "GLYTOUCAN");
  const { value: RING_TYPE } = useTaskItem("RING_TYPE");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          {/* --- Model and experimental data --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model and experimental data" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" color="text.secondary">
              Your model must contain monosaccharides. Please check the documentation for a list of known monosaccharide three-letter codes. Reflection data are required for calculating correlation between model and unbiased density.
            </Typography>
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="XYZIN" {...props} />
              <CCP4i2TaskElement itemName="F_SIGF" {...props} />
              <InlineField label="Change the mask radius around the sugar atoms to" hint="Angstroems">
                <CCP4i2TaskElement itemName="RADIUSIN" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>

              <CCP4i2TaskElement itemName="NEW_SUGAR" {...props} qualifiers={{ guiLabel: "A sugar I want to validate is not yet part of the Chemical Component Dictionary" }} />
            </CCP4i2ContainerElement>

            {/* Define a new sugar — shown only when NEW_SUGAR is true */}
            {newSugar.value && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                  Please define your sugar using the fields below. If your sugar is a polysaccharide, you will have to enter information corresponding to just one ring.
                </Typography>
                <CCP4i2ContainerElement
                  {...props}
                  itemName=""
                  qualifiers={{ initiallyOpen: true }}
                  containerHint="BlockLevel"
                >
                  <InlineField label="Analyse sugar with code">
                    <CCP4i2TaskElement itemName="CODEIN" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                  <InlineField label="and type">
                    <CCP4i2TaskElement itemName="ANOMER" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                  <CCP4i2TaskElement itemName="HAND" {...props} />
                  <CCP4i2TaskElement itemName="RING_TYPE" {...props} />

                  {/* Conformation — pyranose */}
                  {RING_TYPE === "pyranose" && (
                    <InlineField label="Expected minimal energy ring conformation:">
                      <CCP4i2TaskElement itemName="CONFORMATION_PYRANOSE" {...props} qualifiers={{ guiLabel: " " }} />
                    </InlineField>
                  )}

                  {/* Conformation — furanose */}
                  {RING_TYPE === "furanose" && (
                    <InlineField label="Expected minimal energy ring conformation:">
                      <CCP4i2TaskElement itemName="CONFORMATION_FURANOSE" {...props} qualifiers={{ guiLabel: " " }} />
                    </InlineField>
                  )}

                  {/* Ring atoms — pyranose (6 atoms) */}
                  {RING_TYPE === "pyranose" && (
                    <InlineField label="Ring atoms:">
                      <CCP4i2TaskElement itemName="RING_OXYGEN" {...props} qualifiers={{ guiLabel: " " }} />
                    </InlineField>
                  )}
                  {RING_TYPE === "pyranose" && (
                    <>
                      <CCP4i2TaskElement itemName="RING_C1" {...props} />
                      <CCP4i2TaskElement itemName="RING_C2" {...props} />
                      <CCP4i2TaskElement itemName="RING_C3" {...props} />
                      <CCP4i2TaskElement itemName="RING_C4" {...props} />
                      <CCP4i2TaskElement itemName="RING_C5" {...props} />
                    </>
                  )}

                  {/* Ring atoms — furanose (5 atoms) */}
                  {RING_TYPE === "furanose" && (
                    <InlineField label="Ring atoms:">
                      <CCP4i2TaskElement itemName="RING_OXYGEN" {...props} qualifiers={{ guiLabel: " " }} />
                    </InlineField>
                  )}
                  {RING_TYPE === "furanose" && (
                    <>
                      <CCP4i2TaskElement itemName="RING_C1" {...props} />
                      <CCP4i2TaskElement itemName="RING_C2" {...props} />
                      <CCP4i2TaskElement itemName="RING_C3" {...props} />
                      <CCP4i2TaskElement itemName="RING_C4" {...props} />
                    </>
                  )}
                </CCP4i2ContainerElement>
              </>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="settings" label="Settings">
          {/* --- Glycosylation analysis --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Glycosylation analysis" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" color="text.secondary">
              Here you can set up how the diagrams will look like. The schemes will follow the Essentials of glycobiology 3rd edition notation with a choice of colours. They are vector-based and can be used in publications straight away.
            </Typography>
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <InlineField label="Create glycan diagrams in">
                <CCP4i2TaskElement itemName="OLDSTYLE" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              <InlineField label="style, in">
                <CCP4i2TaskElement itemName="VERTICAL" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              <InlineField label="Colour them using the">
                <CCP4i2TaskElement itemName="ESSENTIALS" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              <InlineField label="colour scheme with">
                <CCP4i2TaskElement itemName="INVERT" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              <InlineField label="Validate glycan structures assuming a" hint="expression system">
                <CCP4i2TaskElement itemName="EXPRESSION" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>

              <CCP4i2TaskElement itemName="GLYTOUCAN" {...props} qualifiers={{ guiLabel: "Validate Glycans with GlyTouCan and GlyConnect databases" }} />
              {glyTouCan.value && (
                <>
                  <CCP4i2TaskElement itemName="CLOSESTMATCH" {...props} qualifiers={{ guiLabel: "GlyConnect: Don't look for the closest match, if modelled glycan is not found in the database." }} />
                  <CCP4i2TaskElement itemName="ALLPERMUTATIONS" {...props} qualifiers={{ guiLabel: "Conduct all possible permutation combinations (WARNING: Will take a long time to finish, should only be really used for O-Glycans)." }} />
                </>
              )}
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Performance settings --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Performance settings" }}
            containerHint="FolderLevel"
          >
            <Typography variant="body2" color="text.secondary">
              Here you can tweak performance related settings, such as whether to run Privateer with a single CPU thread or customize the number of CPU threads to run Privateer with. Privateer is tweaked by default to try to obtain best performance possible that is offered by the computer, i.e. use maximum number of threads available.
            </Typography>
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <InlineField label="Run Privateer with" hint="threads.">
                <CCP4i2TaskElement itemName="NUMTHREADS" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
