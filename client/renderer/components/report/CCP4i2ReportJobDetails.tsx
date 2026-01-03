import { useEffect, useState } from "react";
import { CCP4i2ReportElementProps, cssToDict } from "./CCP4i2ReportElement";
import $ from "jquery";
import { SimpleObjectTable } from "../simple-object-table";
/*
 *Handles  <CCP4i2ReportJobDetails key="JobDetails_249" class="" style="" creationtime="19:23 20-Mar-2025" finishtime="19:25 20-Mar-2025" status="Finished"/>
 *
 * */
export const CCP4i2ReportJobDetails: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const [object, setObject] = useState<any | null>(null);

  useEffect(() => {
    const object: any = {};
    $.each(props.item.attributes, function (index, attr) {
      if (!["key", "class", "style"].includes(attr.name))
        object[attr.name] = attr.value;
      setObject(object);
    });
  }, [props.item, props.job]);

  return props.item && <SimpleObjectTable sx={{ mb: 4 }} object={object} />;
};
