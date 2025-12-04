import { useState, useMemo } from "react";
import $ from "jquery";
import { Grid2, ListItem, Paper } from "@mui/material";

import {
  CCP4i2ReportElement,
  CCP4i2ReportElementProps,
} from "./CCP4i2ReportElement";

export const CCP4i2ReportObjectGallery: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const [selected, setSelected] = useState<number | null>(null);
  const childItems = useMemo(() => {
    if (props.item && props.job) {
      const childGraphs = $(props.item).children().toArray();
      if (childGraphs.length > 0) {
        setSelected(0);
      }
      return childGraphs;
    }
    return [];
  }, [props.item, props.job]);
  return (
    <Paper
      sx={{
        borderColor: "primary.main",
        borderTopLeftRadius: "2px",
        borderBottomLeftRadius: "2px",
      }}
    >
      <Grid2 container>
        <Grid2 size={{ xs: 12, sm: 6 }}>
          <Paper
            sx={{
              maxHeight: "20rem",
              overflowY: "auto",
              borderRadius: "2px",
              height: "100%",
            }}
          >
            {childItems &&
              childItems.map((childItem: any, iItem: number) => (
                <ListItem
                  key={`${iItem}`}
                  {...props}
                  sx={iItem === selected ? { border: "2px solid black" } : {}}
                  onClick={() => {
                    setSelected(iItem);
                  }}
                >
                  {$(childItem).attr("title")
                    ? $(childItem).attr("title")
                    : $(childItem).attr("key")
                    ? $(childItem).attr("key")
                    : `Object ${iItem}`}
                </ListItem>
              ))}
          </Paper>
        </Grid2>
        <Grid2 size={{ xs: 12, sm: 6 }}>
          <Paper
            sx={{
              maxHeight: "35rem",
              overflowY: "auto",
              borderRadius: "2px",
              height: "100%",
            }}
          >
            {childItems &&
              childItems.map(
                (childItem: any, iItem: number) =>
                  iItem == selected && (
                    <CCP4i2ReportElement
                      key={`${iItem}`}
                      {...props}
                      iItem={iItem}
                      item={childItem}
                    />
                  )
              )}
          </Paper>
        </Grid2>
      </Grid2>
    </Paper>
  );
};
