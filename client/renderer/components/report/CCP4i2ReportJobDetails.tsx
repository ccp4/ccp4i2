/*
 * Copyright (C) 2025 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
