import { CCP4i2TaskElementProps } from "./task-element";
import { CSimpleElement } from "./csimple";
import { useMemo } from "react";

export const CIntElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const inferredVisibility = useMemo(() => {
    if (!props.visibility) return true;
    if (typeof props.visibility === "function") {
      return props.visibility();
    }
    return props.visibility;
  }, [props.visibility]);

  return inferredVisibility ? <CSimpleElement {...props} type="int" /> : null;
};
