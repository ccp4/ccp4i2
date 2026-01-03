import React, { Children, PropsWithChildren } from "react";
import { Box, SxProps, Theme } from "@mui/material";
import { FIELD_SPACING } from "./field-sizes";

interface FieldRowProps {
  /** Additional sx props */
  sx?: SxProps<Theme>;
  /**
   * When true (default), children share equal width like Grid2 columns.
   * When false, children size to their content/maxWidth.
   */
  equalWidth?: boolean;
}

/**
 * Simple wrapper for placing multiple fields in a horizontal row.
 *
 * By default, fields share equal width like Grid2 columns. On narrow
 * screens, fields wrap to the next line. Use `equalWidth={false}` for
 * content-based sizing when fields form an inline sentence.
 *
 * @example
 * // Equal-width columns (default, like Grid2 xs={6})
 * <FieldRow>
 *   <CCP4i2TaskElement itemName="HYDR_USE" />
 *   <CCP4i2TaskElement itemName="HYDR_ALL" />
 * </FieldRow>
 *
 * @example
 * // Content-based sizing for inline sentences
 * <FieldRow equalWidth={false}>
 *   <CCP4i2TaskElement itemName="SCALE_TYPE" />
 *   <CCP4i2TaskElement itemName="SOLVENT_MASK_TYPE" />
 * </FieldRow>
 */
export const FieldRow: React.FC<PropsWithChildren<FieldRowProps>> = ({
  children,
  sx,
  equalWidth = true,
}) => (
  <Box
    sx={{
      display: "flex",
      flexDirection: "row",
      flexWrap: "wrap",
      gap: FIELD_SPACING.columnGap,
      alignItems: "flex-start",
      ...sx,
    }}
  >
    {equalWidth
      ? Children.map(children, (child) => (
          <Box
            sx={{
              flex: "1 1 0",
              minWidth: "12rem", // Allow wrapping on narrow screens
            }}
          >
            {child}
          </Box>
        ))
      : children}
  </Box>
);

export default FieldRow;
