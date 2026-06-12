import React, { PropsWithChildren, ReactNode } from "react";
import { Box, Card, Stack, SxProps, Theme, Typography } from "@mui/material";

/**
 * Shared visual chrome for "boxed" task widgets (lists, cells, column groups,
 * space groups, row-level containers).
 *
 * Before FieldShell, every boxed widget invented its own frame — border 1 vs 2
 * vs "3px solid", Card vs Box, radius 1 vs 2, ad-hoc hover states, and only some
 * coloured the border by validation state. FieldShell is the single source of
 * truth for that frame so all boxed widgets look identical and respond to
 * validation/hover the same way.
 *
 * Scope is deliberately narrow: FieldShell owns the FRAME and the optional
 * HEADER (title + action slot + error trigger). Each widget keeps full control
 * of its own content.
 */

/** Standard chrome tokens — keep boxed widgets visually consistent. */
export const FIELD_SHELL = {
  /** Border radius (MUI spacing units → 8px). */
  borderRadius: 2,
  /** Border width in px. */
  borderWidth: 1,
  /** Header padding. */
  headerPx: 2,
  headerPy: 1,
  /** Content padding. */
  contentPx: 2,
  contentPy: 1.5,
} as const;

export interface FieldShellProps extends PropsWithChildren {
  /** Heading shown in the header row. Omit for a frame with no header. */
  title?: ReactNode;
  /**
   * Validation-driven border colour token (e.g. "error.light",
   * "warning.light", or "divider"). Defaults to "divider".
   */
  borderColor?: string;
  /** Right-aligned header content — buttons, toggles, etc. */
  action?: ReactNode;
  /** Validation error trigger element, rendered after the action slot. */
  errorTrigger?: ReactNode;
  /** Tints the frame as selected (e.g. CColumnGroup selection). */
  selected?: boolean;
  /** Lift border to primary.light + shadow on hover. Default true. */
  hoverable?: boolean;
  /** Override sx on the outer Card. */
  sx?: SxProps<Theme>;
  /** Override sx on the content Box. */
  contentSx?: SxProps<Theme>;
}

export const FieldShell: React.FC<FieldShellProps> = ({
  title,
  borderColor = "divider",
  action,
  errorTrigger,
  selected = false,
  hoverable = true,
  sx,
  contentSx,
  children,
}) => {
  const hasHeader = title != null || action != null || errorTrigger != null;

  return (
    <Card
      variant="outlined"
      sx={{
        mx: 1,
        borderWidth: FIELD_SHELL.borderWidth,
        borderColor,
        borderRadius: FIELD_SHELL.borderRadius,
        boxShadow: "none",
        bgcolor: selected ? "action.selected" : "background.paper",
        transition: "border-color 0.2s ease, box-shadow 0.2s ease",
        ...(hoverable && {
          "&:hover": {
            borderColor:
              borderColor === "divider" ? "primary.light" : borderColor,
            boxShadow: 1,
          },
        }),
        ...sx,
      }}
    >
      {hasHeader && (
        <Stack
          direction="row"
          alignItems="center"
          justifyContent="space-between"
          spacing={1}
          sx={{
            px: FIELD_SHELL.headerPx,
            py: FIELD_SHELL.headerPy,
            minHeight: "2.5rem",
          }}
        >
          {typeof title === "string" ? (
            <Typography variant="subtitle2" sx={{ fontWeight: 500 }} noWrap>
              {title}
            </Typography>
          ) : (
            title ?? <Box />
          )}
          <Stack direction="row" alignItems="center" spacing={0.5}>
            {action}
            {errorTrigger}
          </Stack>
        </Stack>
      )}
      <Box
        sx={{
          px: FIELD_SHELL.contentPx,
          py: FIELD_SHELL.contentPy,
          pt: hasHeader ? 0 : FIELD_SHELL.contentPy,
          ...contentSx,
        }}
      >
        {children}
      </Box>
    </Card>
  );
};

export default FieldShell;
