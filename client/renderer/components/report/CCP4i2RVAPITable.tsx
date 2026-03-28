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
import $ from "jquery";
import { CCP4i2RVAPIRow } from "./CCP4i2RVAPIRow";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";

export const CCP4i2RVAPITable: React.FC<CCP4i2ReportElementProps> = (props) => {
  const [rows, setRows] = useState<JQuery<HTMLTableRowElement>>();

  useEffect(() => {
    setRows($(props.item).find("tr"));
  }, [props.item, props.job]);

  return (
    <div>
      {rows &&
        rows.map((iRow: number, row) => (
          <CCP4i2RVAPIRow
            key={`${iRow}`}
            iItem={iRow}
            item={row}
            job={props.job}
          />
        ))}
    </div>
  );
};
