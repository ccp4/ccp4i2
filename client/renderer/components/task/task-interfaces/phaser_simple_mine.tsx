import { LinearProgress, Paper, Stack, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useApi } from "../../../api";
import { useJob, usePrevious } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useCallback, useEffect, useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { value: F_SIGFValue } = useTaskItem("F_SIGF");
  const { value: F_OR_IValue, update: updateF_or_I } = useTaskItem("F_OR_I");
  const { value: INPUT_FIXedValue } = useTaskItem("INPUT_FIXED");
  const { value: COMP_BYValue } = useTaskItem("COMP_BY");
  const { value: ID_RMSValue } = useTaskItem("ID_RMS");
  const { value: SGALT_SELECTValue } = useTaskItem("SGALT_SELECT");

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Key files" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="F_OR_I"
            qualifiers={{ guiLabel: "Use Fs or Is" }}
            visibility={() => {
              return [1, 3].includes(F_SIGFValue.contentFlag);
            }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="F_SIGF"
            qualifiers={{ guiLabel: "Reflections" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN"
            qualifiers={{ guiLabel: "Search coordinates" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="INPUT_FIXED"
            qualifiers={{ guiLabel: "Have known partial model" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN_FIXED"
            qualifiers={{ guiLabel: "Known partial model" }}
            visibility={() => {
              return INPUT_FIXedValue;
            }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="FREERFLAG"
            qualifiers={{ guiLabel: "Free R flags" }}
          />
        </CCP4i2ContainerElement>
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Basic parameters" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="NCOPIES"
            qualifiers={{ guiLabel: "Copies to find" }}
          />
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Resolution" }}
            containerHint="BlockLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="RESOLUTION_LOW"
              qualifiers={{ guiLabel: "Low" }}
            />
            <CCP4i2TaskElement
              itemName="RESOLUTION_HIGH"
              {...props}
              qualifiers={{ guiLabel: "High" }}
            />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Extra steps" }}
            containerHint="BlockLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="RUNSHEETBEND"
              qualifiers={{ guiLabel: "Sheet bend" }}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="RUNREFMAC"
              qualifiers={{ guiLabel: "Refmac" }}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2ContainerElement>
        <CCP4i2ContainerElement
          itemName=""
          {...props}
          qualifiers={{ guiLabel: "Scattering in the crystal" }}
          key="Scattering"
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="COMP_BY"
            sx={{ minWidth: "100%" }}
            qualifiers={{ guiLabel: "How to specify scattering content" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="ASUFILE"
            qualifiers={{ guiLabel: "CCP4i2 ASU file" }}
            visibility={() => {
              return COMP_BYValue === "ASU";
            }}
          />
          {/* FieldRow distributes children equally */}
          <FieldRow>
            <CCP4i2TaskElement
              {...props}
              itemName="ASU_NUCLEICACID_MW"
              qualifiers={{ guiLabel: "nucleic acid (Da)" }}
              visibility={() => {
                return COMP_BYValue === "MW";
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="ASU_PROTEIN_MW"
              qualifiers={{ guiLabel: "protein (Da)" }}
              visibility={() => {
                return COMP_BYValue === "MW";
              }}
            />
          </FieldRow>
        </CCP4i2ContainerElement>
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Spacegroups" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="SGALT_SELECT"
            qualifiers={{ guiLabel: "How spacegroups are chosen" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="SGALT_TEST"
            qualifiers={{ guiLabel: "List of spacegroups to try" }}
            visibility={() => {
              return SGALT_SELECTValue === "LIST";
            }}
          />
        </CCP4i2ContainerElement>
        <CCP4i2ContainerElement
          itemName=""
          {...props}
          qualifiers={{ guiLabel: "Similarity of search model" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="ID_RMS"
            qualifiers={{
              guiLabel: "How to specify similarity (i.e. sequence or coords)",
            }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="SEARCHSEQUENCEIDENTITY"
            qualifiers={{ guiLabel: "Sequence identity (0.0-1.0)" }}
            visibility={() => {
              return ID_RMSValue == "ID";
            }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="SEARCHRMS"
            qualifiers={{ guiLabel: "Expected coordinate RMSD (angstroms)" }}
            visibility={() => {
              return ID_RMSValue == "RMS";
            }}
          />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
      <CCP4i2Tab label="Keywords" key="2">
        <CCP4i2ContainerElement
          {...props}
          itemName="keywords"
          qualifiers={{ guiLabel: "" }}
          containerHint="BlockLevel"
          initiallyOpen={true}
        />
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};

export default TaskInterface;
