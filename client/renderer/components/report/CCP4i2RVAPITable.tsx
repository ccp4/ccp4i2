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
