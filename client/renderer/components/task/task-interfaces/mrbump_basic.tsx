import React from "react";
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { InlineField } from "../task-elements/inline-field";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  const searchPdb = useBoolToggle(useTaskItem, "SEARCH_PDB");
  const searchAfdb = useBoolToggle(useTaskItem, "SEARCH_AFDB");
  const includeLocal = useBoolToggle(useTaskItem, "INCLUDE");

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input Data ===== */}
        <CCP4i2Tab label="Input Data">
          {/* Target Sequence */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Target Sequence" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="ASUIN"
              {...props}
              qualifiers={{ guiLabel: "AU contents" }}
            />
            <Typography variant="body2" sx={{ fontStyle: "italic", my: 0.5 }}>
              If a suitable ASU is not available above, you can press the cross
              &amp; then button to quickly create one.
            </Typography>
            <InlineField label="The number of monomers to search for" hint=" ">
              <CCP4i2TaskElement
                itemName="NMON"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Experimental Data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Experimental Data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <CCP4i2TaskElement
              itemName="FREERFLAG"
              {...props}
              qualifiers={{ guiLabel: "Free R set" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Tab 2: Search Models ===== */}
        <CCP4i2Tab label="Search Models">
          {/* Model databases */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model databases" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="SEARCH_PDB"
              {...props}
              qualifiers={{
                guiLabel: "Search PDB for possible MR search models",
              }}
              onChange={searchPdb.onChange}
            />
            {searchPdb.value && (
              <InlineField
                label="Non-redundancy level for homologue search:"
                sx={{ pl: 3 }}
              >
                <CCP4i2TaskElement
                  itemName="REDUNDANCYLEVEL"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}

            <CCP4i2TaskElement
              itemName="SEARCH_AFDB"
              {...props}
              qualifiers={{
                guiLabel: "Search EBI-AFDB for possible MR search models",
              }}
              onChange={searchAfdb.onChange}
            />
            {searchAfdb.value && (
              <InlineField
                label="EBI-AFDB pLDDT residue score cut-off:"
                sx={{ pl: 3 }}
              >
                <CCP4i2TaskElement
                  itemName="AFDBLEVEL"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
            )}

            <InlineField
              label="Maximum no. of search models to create:"
              sx={{ mt: 1 }}
            >
              <CCP4i2TaskElement
                itemName="MRMAX"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Optional Settings */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Optional Settings" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="HHPREDIN"
              {...props}
              qualifiers={{ guiLabel: "HHPred hhr file" }}
            />
            <CCP4i2TaskElement
              itemName="PDBLOCAL"
              {...props}
              qualifiers={{ guiLabel: "Path to local PDB mirror" }}
            />
          </CCP4i2ContainerElement>

          {/* Local coordinate files */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Local coordinate files to be used as search models",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="INCLUDE"
              {...props}
              qualifiers={{ guiLabel: "Include local files" }}
              onChange={includeLocal.onChange}
            />
            {includeLocal.value && (
              <CCP4i2TaskElement itemName="XYZIN_LIST" {...props} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Tab 3: Options ===== */}
        <CCP4i2Tab label="Options">
          {/* Molecular Replacement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Molecular Replacement" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Number of cores for Phaser (maximum=10)">
              <CCP4i2TaskElement
                itemName="PJOBS"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Refinement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Refinement" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Number of refinement cycles in Refmac">
              <CCP4i2TaskElement
                itemName="NCYC"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Model Building */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model Building" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="BUCC"
              {...props}
              qualifiers={{
                guiLabel: "Run Buccaneer after refinement",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
