import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { CCP4i2ContainerElement } from "./ccontainer";
import { Typography } from "@mui/material";

export const CRangeElement: React.FC<CCP4i2TaskElementProps> = (props) => (
  <CCP4i2ContainerElement
    {...props}
    itemName=""
    qualifiers={props.qualifiers}
    containerHint="RowLevel"
  >
    <Typography variant="body2">Resolution</Typography>
    <CCP4i2TaskElement job={props.job} itemName={`${props.itemName}.start`} />
    <CCP4i2TaskElement job={props.job} itemName={`${props.itemName}.end`} />
  </CCP4i2ContainerElement>
);
