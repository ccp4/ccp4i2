import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useApi } from "../../../api";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { Grid2 } from "@mui/material";
import { RDKitView } from "../../rdkit-view";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { job } = props;
  const { useTaskItem, useFileContent } = useJob(job.id);
  const { value: MOLSMILESORSKETCH } = useTaskItem("MOLSMILESORSKETCH");
  const { value: ATOMMATCHOPTION } = useTaskItem("ATOMMATCHOPTION");
  const { data: MOLINContent, mutate: mutateMOLINContent } =
    useFileContent("MOLIN");
  const { data: SMILESFILEContent } = useFileContent("SMILESFILE");

  ATOMMATCHOPTION;
  return (
    <>
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Main inputs" key="1">
          <CCP4i2TaskElement
            {...props}
            itemName="MOLSMILESORSKETCH"
            sx={{ width: "100%" }}
            qualifiers={{ guiLabel: "Ligand geometry provided as" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="SMILESIN"
            qualifiers={{ guiLabel: "Smiles" }}
            sx={{ width: "100%" }}
            visibility={() => {
              return MOLSMILESORSKETCH === "SMILES";
            }}
          />
          <Grid2 container>
            <Grid2 size={{ xs: 12, sm: 6 }}>
              <CCP4i2TaskElement
                {...props}
                itemName="MOLIN"
                qualifiers={{ guiLabel: "MDL Mol File" }}
                visibility={() => {
                  return MOLSMILESORSKETCH === "MOL";
                }}
              />
            </Grid2>
            <Grid2 size={{ xs: 12, sm: 6 }}>
              {MOLSMILESORSKETCH === "MOL" && MOLINContent && (
                <RDKitView smiles={MOLINContent} />
              )}
            </Grid2>
          </Grid2>
          <CCP4i2TaskElement
            {...props}
            itemName="DICTIN2"
            qualifiers={{ guiLabel: "CIF Dictionary" }}
            visibility={() => {
              return MOLSMILESORSKETCH === "DICT";
            }}
          />
          <Grid2 container>
            <Grid2 size={{ xs: 12, sm: 6 }}>
              <CCP4i2TaskElement
                {...props}
                itemName="SMILESFILEIN"
                qualifiers={{ guiLabel: "File containg SMILES string" }}
                visibility={() => {
                  return MOLSMILESORSKETCH === "SMILESFILE";
                }}
              />
            </Grid2>
            <Grid2 size={{ xs: 12, sm: 6 }}>
              {MOLSMILESORSKETCH === "SMILESFILE" && SMILESFILEContent && (
                <RDKitView smiles={SMILESFILEContent.trim()} />
              )}
            </Grid2>
          </Grid2>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Atom/residue naming",
              initiallyOpen: true,
            }}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="TLC"
              qualifiers={{ guiLabel: "Three letter code for result" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="ATOMMATCHOPTION"
              qualifiers={{ guiLabel: "Atom name matching" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="MATCHTLC"
              qualifiers={{ guiLabel: "Three letter code to match names to" }}
              visibility={() => {
                return ATOMMATCHOPTION === "MONLIBCODE";
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="DICTIN"
              qualifiers={{ guiLabel: "Local dictionary" }}
              visibility={() => {
                return ATOMMATCHOPTION === "LOCALDICT";
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>{" "}
    </>
  );
};

export default TaskInterface;
