import { useEffect, useState } from "react";
import $ from "jquery";
import { Grid, Grid2 } from "@mui/material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { CCP4i2ApplicationOutputView } from "./CCP4i2ApplicationOutputView";

export const CCP4i2RVAPIRow: React.FC<CCP4i2ReportElementProps> = (props) => {
  const [content, setContent] = useState<
    JQuery<React.ReactElement> | undefined
  >();

  useEffect(() => {
    let newContent = $(props.item)
      .find("td")
      .map((iCol, col) => {
        if ($(col).find("div[data-renderer]").length > 0) {
          // For RVAPI graphs, extract the ID from the ccp4_data element
          const dataElement = $(col).find("ccp4\\:ccp4_data, ccp4_data, ns0\\:ccp4_data").first();
          const graphKey = dataElement.attr("id");

          return (
            <Grid2 size={{ xs: 6 }} key={iCol}>
              <CCP4i2ApplicationOutputView
                output={col}
                jobId={props.job?.id}
                graphId={graphKey}
              />
            </Grid2>
          );
        } else {
          return (
            <Grid2
              size={{ xs: 6 }}
              key={iCol}
              dangerouslySetInnerHTML={{ __html: col.innerHTML }}
            />
          );
        }
      });
    setContent(newContent);
  }, [props.job, props.item]);

  return <Grid2 container>{content}</Grid2>;
};
