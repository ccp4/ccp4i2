import { useState, useMemo } from "react";
import $ from "jquery";
import { galleryItemLabel } from "../../lib/report-gallery";
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
  // One object is not a gallery: show it, without a list to choose from
  // (Make Ligand's single 2D drawing sat beside a list of one, titled with
  // its internal key, "Div_3").
  if (childItems.length === 1) {
    return (
      <Paper sx={{ borderRadius: "2px", display: "inline-block" }}>
        <CCP4i2ReportElement {...props} iItem={0} item={childItems[0]} />
      </Paper>
    );
  }
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
                  sx={iItem === selected ? { border: "2px solid black" } : {}}
                  onClick={() => {
                    setSelected(iItem);
                  }}
                >
                  {galleryItemLabel(childItem, iItem)}
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
