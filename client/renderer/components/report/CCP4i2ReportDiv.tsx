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
import { Fragment, useEffect, useState } from "react";
import $ from "jquery";
import { Grid2 } from "@mui/material";
import {
  CCP4i2ReportElement,
  CCP4i2ReportElementProps,
  cssToDict,
} from "./CCP4i2ReportElement";

export const CCP4i2ReportDiv: React.FC<CCP4i2ReportElementProps> = (props) => {
  const [nFloatingChildren, setNFloatingChildren] = useState(1);

  useEffect(() => {
    if (props.item) {
      let nFloatingChildren = 0;
      for (var child of $(props.item).children()) {
        try {
          const attrValue = $(child).attr("style");
          if (attrValue === undefined) {
            continue;
          }
          var childCssDict = {}; //cssToDict(attrValue);
          if (Object.keys(childCssDict).includes("float")) {
            const oldStyle = attrValue;
            const fixedStyle = oldStyle
              .replace("float:left;", "")
              .replace("float:right;", "");
            //console.log({ oldStyle, fixedStyle });
            $(child).attr("style", fixedStyle);
            nFloatingChildren += 1;
          }
        } catch (err) {}
      }
      setNFloatingChildren(nFloatingChildren);
    }
  }, [props.item]);

  return nFloatingChildren > 0 ? (
    <Grid2 container>
      {$(props.item)
        .children()
        .map((iChild, child) => (
          <Grid2 key={iChild} size={{ xs: 12 / nFloatingChildren }}>
            <CCP4i2ReportElement iItem={iChild} item={child} job={props.job} />
          </Grid2>
        ))}
    </Grid2>
  ) : (
    <Fragment>
      {$(props.item)
        .children()
        .map((iChild, child) => (
          <div key={iChild} className="CCP4i2ReportDiv">
            <CCP4i2ReportElement iItem={iChild} item={child} job={props.job} />
          </div>
        ))}
    </Fragment>
  );
};
