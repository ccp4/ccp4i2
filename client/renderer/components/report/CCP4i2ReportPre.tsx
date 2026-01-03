import { useEffect, useState } from "react";
import $ from "jquery";
import { CCP4i2ReportElementProps, cssToDict } from "./CCP4i2ReportElement";

export const CCP4i2ReportPre: React.FC<CCP4i2ReportElementProps> = (props) => {
  const [style, setStyle] = useState({});
  const [innerHTML, setInnerHTML] = useState("");

  useEffect(() => {
    const style = $(props.item).attr("style");
    if (style) {
      //setStyle(cssToDict(style));
    }
    setInnerHTML(props.item.innerHTML);
  }, [props.item, props.job]);

  return <pre style={style} dangerouslySetInnerHTML={{ __html: innerHTML }} />;
};
