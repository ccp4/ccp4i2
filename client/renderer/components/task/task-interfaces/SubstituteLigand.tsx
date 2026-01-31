import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  const { value: ligandAs } = useTaskItem("LIGANDAS");
  const { value: obsAs } = useTaskItem("OBSAS");

  return (
    <>
      <CCP4i2ContainerElement
        itemName=""
        key="Ligand geometry"
        {...props}
        qualifiers={{ guiLabel: "Ligand geometry" }}
        containerHint="BlockLevel"
      >
        <CCP4i2TaskElement {...props} itemName="LIGANDAS" />
        <CCP4i2TaskElement
          {...props}
          itemName="SMILESIN"
          visibility={() => ligandAs === "SMILES"}
          qualifiers={{ guiMode: "multiLine" }}
        />
        <CCP4i2TaskElement
          {...props}
          itemName="MOLIN"
          visibility={() => ligandAs === "MOL"}
        />
        <CCP4i2TaskElement
          {...props}
          itemName="DICTIN"
          visibility={() => ligandAs === "DICT"}
        />
      </CCP4i2ContainerElement>
      <CCP4i2ContainerElement
        itemName=""
        key="Refinement type"
        {...props}
        qualifiers={{ guiLabel: "Refinement type" }}
        containerHint="BlockLevel"
      >
        <CCP4i2TaskElement {...props} itemName="PIPELINE" />
      </CCP4i2ContainerElement>
      <CCP4i2ContainerElement
        itemName=""
        key="Starting coordinates"
        {...props}
        qualifiers={{ guiLabel: "Starting coordinates" }}
        containerHint="BlockLevel"
      >
        <CCP4i2TaskElement {...props} itemName="XYZIN" />
      </CCP4i2ContainerElement>
      <CCP4i2ContainerElement
        itemName=""
        key="Reflection data"
        {...props}
        qualifiers={{ guiLabel: "Reflection data" }}
        containerHint="BlockLevel"
      >
        <CCP4i2TaskElement key="OBSAS" {...props} itemName="OBSAS" />
        <CCP4i2TaskElement
          key="UNMERGED"
          {...props}
          itemName="UNMERGEDFILES"
          visibility={() => obsAs === "UNMERGED"}
        />
        <CCP4i2TaskElement
          key="MERGED"
          {...props}
          itemName="F_SIGF_IN"
          visibility={() => obsAs === "MERGED"}
        />
      </CCP4i2ContainerElement>
      <CCP4i2ContainerElement
        itemName=""
        key="Free R"
        {...props}
        qualifiers={{ guiLabel: "Free R flag" }}
        containerHint="BlockLevel"
      >
        <CCP4i2TaskElement {...props} itemName="FREERFLAG_IN" />
      </CCP4i2ContainerElement>
    </>
  );
};

export default TaskInterface;
