import React, { Children, PropsWithChildren } from "react";
import { Box, SxProps, Theme } from "@mui/material";
import {
  FIELD_SPACING,
  FIELD_SIZES,
  FieldSize,
  getContainerSizeStyles,
} from "./field-sizes";

interface FieldRowProps {
  /** Additional sx props */
  sx?: SxProps<Theme>;
  /**
   * When true (default), children share equal width like Grid2 columns.
   * When false, children use the `size` prop or their natural width.
   */
  equalWidth?: boolean;
  /**
   * Constrain each child to this size. Only applies when equalWidth={false}.
   * Children are full-width by default; use this to limit their width.
   */
  size?: FieldSize;
}

/**
 * Simple wrapper for placing multiple fields in a horizontal row.
 *
 * DESIGN: Base widgets are full-width by default. FieldRow controls sizing.
 *
 * @example
 * // Equal-width columns (default) - children fill available space equally
 * <FieldRow>
 *   <CCP4i2TaskElement itemName="HYDR_USE" />
 *   <CCP4i2TaskElement itemName="HYDR_ALL" />
 * </FieldRow>
 *
 * @example
 * // Constrained children - each child is limited to 'xs' width (8rem)
 * <FieldRow equalWidth={false} size="xs">
 *   <CCP4i2TaskElement itemName="CELL_A" />
 *   <CCP4i2TaskElement itemName="CELL_B" />
 *   <CCP4i2TaskElement itemName="CELL_C" />
 * </FieldRow>
 *
 * @example
 * // Natural sizing - children use their content width (not recommended)
 * <FieldRow equalWidth={false}>
 *   <CCP4i2TaskElement itemName="SCALE_TYPE" />
 * </FieldRow>
 */
export const FieldRow: React.FC<PropsWithChildren<FieldRowProps>> = ({
  children,
  sx,
  equalWidth = true,
  size,
}) => {
  // When equalWidth is true, children share space equally
  if (equalWidth) {
    return (
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
        {Children.map(children, (child) => (
          <Box
            sx={{
              flex: "1 1 0",
              minWidth: "12rem", // Allow wrapping on narrow screens
            }}
          >
            {child}
          </Box>
        ))}
      </Box>
    );
  }

  // When equalWidth is false, use size prop to constrain children
  // If no size is given, children use their natural (full) width
  const childSx = size ? getContainerSizeStyles(size) : undefined;

  return (
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
      {childSx
        ? Children.map(children, (child) => <Box sx={childSx}>{child}</Box>)
        : children}
    </Box>
  );
};

export default FieldRow;
