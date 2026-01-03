import { CCP4i2TaskElementProps } from "./task-element";
import { CSimpleElement } from "./csimple";
import { useInferredVisibility } from "./hooks/useInferredVisibility";

export const CIntElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const isVisible = useInferredVisibility(props.visibility);

  return isVisible ? <CSimpleElement {...props} type="int" /> : null;
};
