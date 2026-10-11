import { useState, useMemo } from "react";
import $ from "jquery";
import { galleryItemLabel } from "../../lib/report-gallery";
import {
  Grid2,
  List,
  ListItemButton,
  ListItemText,
  Paper,
} from "@mui/material";

import {
  CCP4i2ReportElement,
  CCP4i2ReportElementProps,
} from "./CCP4i2ReportElement";

/**
 * A gallery of report objects (graphs, mostly): a dense list of their
 * titles down the left, the chosen one drawn at full width on the right.
 *
 * The list takes a third of the width (a quarter on wide screens) and the
 * object the rest; both columns share one height, so a long list scrolls
 * beside the plot rather than stretching the panel. The choice is shown
 * the way MUI shows a selected list item, not by a black border.
 */
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
    <Paper variant="outlined" sx={{ borderRadius: 1, my: 1, overflow: "hidden" }}>
      <Grid2 container sx={{ minHeight: "22rem" }}>
        <Grid2
          size={{ xs: 12, sm: 4, md: 3 }}
          sx={{
            borderRight: { sm: 1 },
            borderBottom: { xs: 1, sm: 0 },
            borderColor: { xs: "divider", sm: "divider" },
            maxHeight: "32rem",
            overflowY: "auto",
          }}
        >
          <List dense disablePadding>
            {childItems.map((childItem: any, iItem: number) => (
              <ListItemButton
                key={`${iItem}`}
                selected={iItem === selected}
                onClick={() => setSelected(iItem)}
                sx={{ py: 0.5 }}
              >
                <ListItemText
                  primary={galleryItemLabel(childItem, iItem)}
                  primaryTypographyProps={{ variant: "body2", noWrap: false }}
                />
              </ListItemButton>
            ))}
          </List>
        </Grid2>
        <Grid2
          size={{ xs: 12, sm: 8, md: 9 }}
          sx={{ maxHeight: "32rem", overflowY: "auto", p: 1 }}
        >
          {childItems.map(
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
        </Grid2>
      </Grid2>
    </Paper>
  );
};
