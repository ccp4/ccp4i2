import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { InlineField } from "../../task-elements/inline-field";
import { useBoolToggle } from "../../task-elements/shared-hooks";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MAPTYPE } = useTaskItem("MAPTYPE");
  const { value: scaleOn, onChange: onScaleChange } = useBoolToggle(
    useTaskItem,
    "SCALE"
  );
  const { value: mapOutputOn, onChange: onMapOutputChange } = useBoolToggle(
    useTaskItem,
    "MAP_OUTPUT"
  );

  if (!container) return <LinearProgress />;

  const isFobsFobs = MAPTYPE === "fobsfobs";

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Type of map */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Type of map to create" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="MAPTYPE" {...props} />
      </CCP4i2ContainerElement>

      {/* First dataset */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "First dataset" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="F_SIGF1" {...props} />
        <CCP4i2TaskElement itemName="ABCD1" {...props} />
      </CCP4i2ContainerElement>

      {/* Second dataset - only for Fobs-Fobs */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Second dataset" }}
        containerHint="FolderLevel"
        visibility={isFobsFobs}
      >
        <CCP4i2TaskElement itemName="F_SIGF2" {...props} />
        <CCP4i2TaskElement itemName="ABCD2" {...props} />
      </CCP4i2ContainerElement>

      {/* Basic operations */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Basic operations on map coefficients" }}
        containerHint="FolderLevel"
      >
        <InlineField
          label="Sharpen or blur the output map using an isotropic B-factor:"
          hint="Use a negative B-factor to sharpen the map, or a positive one to blur it"
        >
          <CCP4i2TaskElement
            itemName="B_VALUE"
            {...props}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
        <CCP4i2TaskElement
          itemName="SCALE"
          {...props}
          qualifiers={{ guiLabel: "Scale the observations" }}
          visibility={isFobsFobs}
          onChange={onScaleChange}
        />
        <CCP4i2TaskElement
          itemName="F1_TO_F2"
          {...props}
          visibility={scaleOn && isFobsFobs}
        />
        <InlineField label="Choose a different resolution limit">
          <CCP4i2TaskElement
            itemName="RESOLUTION"
            {...props}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
      </CCP4i2ContainerElement>

      {/* Output options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Output options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="MAP_OUTPUT"
          {...props}
          qualifiers={{ guiLabel: "Produce a map file too" }}
          onChange={onMapOutputChange}
        />
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          Grid sampling should not be changed unless you know what you are
          doing.
        </Typography>
        <CCP4i2TaskElement
          itemName="INDEX_U"
          {...props}
          qualifiers={{ guiLabel: "u" }}
          visibility={mapOutputOn}
        />
        <CCP4i2TaskElement
          itemName="INDEX_V"
          {...props}
          qualifiers={{ guiLabel: "v" }}
          visibility={mapOutputOn}
        />
        <CCP4i2TaskElement
          itemName="INDEX_W"
          {...props}
          qualifiers={{ guiLabel: "w" }}
          visibility={mapOutputOn}
        />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
