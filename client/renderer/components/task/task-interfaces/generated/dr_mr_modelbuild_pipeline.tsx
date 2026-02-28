import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";
import { useBoolToggle } from "../../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MERGED_OR_UNMERGED } = useTaskItem("MERGED_OR_UNMERGED");
  const { value: XYZINORMRBUMP } = useTaskItem("XYZINORMRBUMP");
  const { value: LIGANDAS } = useTaskItem("LIGANDAS");
  const searchPdb = useBoolToggle(useTaskItem, "SEARCH_PDB");
  const searchAfdb = useBoolToggle(useTaskItem, "SEARCH_AFDB");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data and Protocol">
          {/* --- Input data --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input data" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Input data type (run import merged, use merged from 'import merged' job or unmerged:">
              <CCP4i2TaskElement itemName="MERGED_OR_UNMERGED" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            {MERGED_OR_UNMERGED === "UNMERGED" && (
              <CCP4i2TaskElement itemName="UNMERGEDFILES" {...props} />
            )}
            {MERGED_OR_UNMERGED === "MERGED" && (
              <>
                <CCP4i2TaskElement itemName="F_SIGF_IN" {...props} />
                <CCP4i2TaskElement itemName="FREER_IN" {...props} />
              </>
            )}
            {MERGED_OR_UNMERGED === "MERGED_F" && (
              <CCP4i2TaskElement itemName="HKLIN" {...props} />
            )}
          </CCP4i2ContainerElement>

          {/* --- Search model --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Search model" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZINORMRBUMP" {...props} />
            {XYZINORMRBUMP === "XYZINPUT" && (
              <CCP4i2TaskElement itemName="XYZIN" {...props} />
            )}
            {XYZINORMRBUMP === "MRBUMP" && (
              <>
                <CCP4i2TaskElement itemName="SEARCH_PDB" {...props} qualifiers={{ guiLabel: "Search PDB for possible MR search models" }} />
                {searchPdb.value && (
                  <>
                    <Typography variant="body2" color="text.secondary" sx={{ pl: 3 }}>
                      Non-redundancy level for homologue search:
                    </Typography>
                    <InlineField label="" sx={{ pl: 3 }}>
                      <CCP4i2TaskElement itemName="REDUNDANCYLEVEL" {...props} qualifiers={{ guiLabel: " " }} />
                    </InlineField>
                  </>
                )}
                <CCP4i2TaskElement itemName="SEARCH_AFDB" {...props} qualifiers={{ guiLabel: "Search EBI-AFDB for possible MR search models" }} />
                {searchAfdb.value && (
                  <>
                    <Typography variant="body2" color="text.secondary" sx={{ pl: 3 }}>
                      EBI-AFDB pLDDT residue score cut-off:
                    </Typography>
                    <InlineField label="" sx={{ pl: 3 }}>
                      <CCP4i2TaskElement itemName="AFDBLEVEL" {...props} qualifiers={{ guiLabel: " " }} />
                    </InlineField>
                  </>
                )}
                <Typography variant="body2" color="text.secondary">
                  Maximum no. of search models to create:
                </Typography>
                <CCP4i2TaskElement itemName="MRMAX" {...props} />
              </>
            )}
          </CCP4i2ContainerElement>

          {/* --- Sequence of target model --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Sequence of target model" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="ASUIN" {...props} />
            <InlineField label="The number of monomers to search for">
              <CCP4i2TaskElement itemName="NMON" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* --- Options --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="AUTOCUTOFF" {...props} qualifiers={{ guiLabel: "Run Aimless twice, first to find resolution limit" }} />
            <CCP4i2TaskElement itemName="RUNACORN" {...props} qualifiers={{ guiLabel: "Run phase refinement with acorn before model building" }} />
            <InlineField label="Run" hint="cycles of restrained refinement after MR">
              <CCP4i2TaskElement itemName="REFMAC_NCYC" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Run" hint="model building pipeline iterations">
              <CCP4i2TaskElement itemName="BUCC_NCYC" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* --- Ligand geometry --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Ligand geometry" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Format in which geometry will be specified:">
              <CCP4i2TaskElement itemName="LIGANDAS" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            {LIGANDAS === "MOL" && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ initiallyOpen: true }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement itemName="MOLIN" {...props} />
              </CCP4i2ContainerElement>
            )}
            {LIGANDAS === "DICT" && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ initiallyOpen: true }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement itemName="DICTIN" {...props} />
              </CCP4i2ContainerElement>
            )}
            {LIGANDAS === "SMILES" && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ initiallyOpen: true }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement itemName="SMILESIN" {...props} />
              </CCP4i2ContainerElement>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="advanced" label="Advanced options">
          {/* --- Advanced model building pipeline options --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced model building pipeline options" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Use" hint="model building pipeline (Buccaneer or Modelcraft)">
              <CCP4i2TaskElement itemName="BUCCANEER_OR_MODELCRAFT" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <Typography variant="body2" color="text.secondary">
              The modelcraft pipeline is now the default option and recommended in most cases
            </Typography>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
