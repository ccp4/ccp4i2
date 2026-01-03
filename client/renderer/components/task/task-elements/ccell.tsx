import { CCP4i2TaskElementProps } from "./task-element";
import { CCP4i2ContainerElement } from "./ccontainer";
import { FIELD_SIZES } from "./field-sizes";

export const CCellElement: React.FC<CCP4i2TaskElementProps> = (props) => (
  <CCP4i2ContainerElement
    {...props}
    qualifiers={props.qualifiers}
    containerHint="RowLevel"
    elementSx={{ width: FIELD_SIZES.xs }}
  />
);
