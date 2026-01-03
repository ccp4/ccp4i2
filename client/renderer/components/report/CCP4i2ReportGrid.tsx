/**
 * Grid layout components for CCP4i2 reports.
 *
 * These components render the GridContainer and GridItem elements
 * from the Python report module as MUI Grid2 components.
 *
 * Usage in Python report:
 *   grid = self.addGrid(spacing=2)
 *   left = grid.addItem(xs=12, md=8)
 *   left.append(self.addTable(...))
 */

import { useMemo } from "react";
import $ from "jquery";
import { Grid2, Box } from "@mui/material";
import { Job } from "../../types/models";
import { CCP4i2ReportElement } from "./CCP4i2ReportElement";

export interface CCP4i2ReportGridProps {
  iItem: number;
  item: HTMLElement;
  job: Job;
}

/**
 * Grid container component - maps to <Grid2 container>
 *
 * Attributes from Python:
 *   - spacing: Gap between items (0-10)
 *   - direction: row, column, row-reverse, column-reverse
 *   - justifyContent: flex-start, center, flex-end, space-between, etc.
 *   - alignItems: flex-start, center, flex-end, stretch, baseline
 *   - wrap: wrap or nowrap
 */
export const CCP4i2ReportGridContainer: React.FC<CCP4i2ReportGridProps> = ({
  iItem,
  item,
  job,
}) => {
  const gridProps = useMemo(() => {
    const $item = $(item);
    const props = {
      spacing: parseInt($item.attr("spacing") || "2", 10),
      direction: ($item.attr("direction") || "row") as
        | "row"
        | "column"
        | "row-reverse"
        | "column-reverse",
      justifyContent: ($item.attr("justifyContent") || "flex-start") as
        | "flex-start"
        | "center"
        | "flex-end"
        | "space-between"
        | "space-around"
        | "space-evenly",
      alignItems: ($item.attr("alignItems") || "stretch") as
        | "flex-start"
        | "center"
        | "flex-end"
        | "stretch"
        | "baseline",
      wrap: $item.attr("wrap") !== "nowrap" ? "wrap" : "nowrap",
    };
    console.log("[Grid] Container props:", props, "children:", $item.children().length);
    return props;
  }, [item]);

  const children = useMemo(() => {
    return $(item)
      .children()
      .map((i: number, child: any) => (
        <CCP4i2ReportElement key={`${iItem}_${i}`} iItem={i} item={child} job={job} />
      ))
      .toArray();
  }, [item, iItem, job]);

  return (
    <Grid2
      container
      spacing={gridProps.spacing}
      direction={gridProps.direction}
      justifyContent={gridProps.justifyContent}
      alignItems={gridProps.alignItems}
      sx={{ flexWrap: gridProps.wrap }}
    >
      {children}
    </Grid2>
  );
};

/**
 * Grid item component - maps to <Grid2 size={{...}}>
 *
 * Attributes from Python:
 *   - xs: Columns at extra-small (0-600px)
 *   - sm: Columns at small (600-900px)
 *   - md: Columns at medium (900-1200px)
 *   - lg: Columns at large (1200-1536px)
 *   - xl: Columns at extra-large (1536px+)
 *
 * Special values:
 *   - 'auto': Size based on content
 *   - 'grow': Grow to fill available space
 */
export const CCP4i2ReportGridItem: React.FC<CCP4i2ReportGridProps> = ({
  iItem,
  item,
  job,
}) => {
  const sizeProps = useMemo(() => {
    const $item = $(item);
    const parseSize = (value: string | undefined): number | "auto" | "grow" | undefined => {
      if (!value) return undefined;
      if (value === "auto") return "auto";
      if (value === "grow" || value === "true") return "grow";
      const num = parseInt(value, 10);
      return isNaN(num) ? undefined : num;
    };

    const xs = parseSize($item.attr("xs"));
    const sm = parseSize($item.attr("sm"));
    const md = parseSize($item.attr("md"));
    const lg = parseSize($item.attr("lg"));
    const xl = parseSize($item.attr("xl"));

    // Build size object, only including defined values
    const size: Record<string, number | "auto" | "grow"> = {};
    if (xs !== undefined) size.xs = xs;
    if (sm !== undefined) size.sm = sm;
    if (md !== undefined) size.md = md;
    if (lg !== undefined) size.lg = lg;
    if (xl !== undefined) size.xl = xl;

    // Default to full width if nothing specified
    if (Object.keys(size).length === 0) {
      size.xs = 12;
    }

    console.log("[Grid] Item size:", size, "children:", $item.children().length);
    return size;
  }, [item]);

  const children = useMemo(() => {
    return $(item)
      .children()
      .map((i: number, child: any) => (
        <CCP4i2ReportElement key={`${iItem}_${i}`} iItem={i} item={child} job={job} />
      ))
      .toArray();
  }, [item, iItem, job]);

  return <Grid2 size={sizeProps}>{children}</Grid2>;
};

/**
 * Grid row component - semantic wrapper for a row in a grid
 *
 * This is essentially the same as GridContainer but provides
 * semantic clarity when building complex multi-row layouts.
 */
export const CCP4i2ReportGridRow: React.FC<CCP4i2ReportGridProps> = (props) => {
  return <CCP4i2ReportGridContainer {...props} />;
};

/**
 * Diagnostic display component for report errors/warnings
 *
 * Renders diagnostic messages collected during report generation.
 */
export const CCP4i2ReportDiagnostics: React.FC<CCP4i2ReportGridProps> = ({
  item,
}) => {
  const diagnostics = useMemo(() => {
    const $item = $(item);
    const hasErrors = $item.attr("hasErrors") === "true";
    const hasWarnings = $item.attr("hasWarnings") === "true";
    const count = parseInt($item.attr("count") || "0", 10);

    const items = $item
      .children("CCP4i2ReportDiagnostic")
      .map((_: number, child: any) => {
        const $child = $(child);
        return {
          level: $child.attr("level") || "info",
          code: $child.attr("code") || "",
          message: $child.attr("message") || "",
          location: $child.attr("location") || "",
        };
      })
      .toArray();

    return { hasErrors, hasWarnings, count, items };
  }, [item]);

  if (diagnostics.count === 0) return null;

  return (
    <Box
      sx={{
        mt: 2,
        p: 2,
        borderRadius: 1,
        backgroundColor: diagnostics.hasErrors
          ? "error.light"
          : diagnostics.hasWarnings
          ? "warning.light"
          : "info.light",
      }}
    >
      {diagnostics.items.map((diag, i) => (
        <Box
          key={i}
          sx={{
            display: "flex",
            alignItems: "flex-start",
            gap: 1,
            mb: 1,
            "&:last-child": { mb: 0 },
          }}
        >
          <Box
            component="span"
            sx={{
              fontWeight: "bold",
              color:
                diag.level === "error" || diag.level === "critical"
                  ? "error.main"
                  : diag.level === "warning"
                  ? "warning.main"
                  : "info.main",
            }}
          >
            [{diag.code}]
          </Box>
          <Box component="span">{diag.message}</Box>
          {diag.location && (
            <Box component="span" sx={{ color: "text.secondary", fontSize: "0.85em" }}>
              ({diag.location})
            </Box>
          )}
        </Box>
      ))}
    </Box>
  );
};
