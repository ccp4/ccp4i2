import React, { PropsWithChildren, useCallback, useState } from "react";
import {
  Collapse,
  IconButton,
  Stack,
  Typography,
  SxProps,
  Theme,
} from "@mui/material";
import {
  ExpandMore as ExpandMoreIcon,
  ChevronRight as ChevronRightIcon,
} from "@mui/icons-material";

export interface ExpandableSectionProps extends PropsWithChildren {
  /**
   * Whether the section is currently expanded.
   * If not provided, the component manages its own state.
   */
  expanded?: boolean;

  /**
   * Callback when expand/collapse state changes.
   */
  onToggle?: (expanded: boolean) => void;

  /**
   * Whether the section is disabled (cannot be toggled).
   */
  disabled?: boolean;

  /**
   * Whether the section should be expanded initially.
   * Only used when `expanded` prop is not provided (uncontrolled mode).
   * @default false
   */
  defaultExpanded?: boolean;

  /**
   * Force the section to be expanded (e.g., due to validation errors).
   * When true, the user cannot collapse the section.
   */
  forceExpanded?: boolean;

  /**
   * Title shown in the collapsed header.
   */
  title?: string;

  /**
   * Title shown when force expanded (e.g., "Required Options (Error)").
   */
  forceExpandedTitle?: string;

  /**
   * Whether the section is in an error state (affects styling).
   */
  hasError?: boolean;

  /**
   * Additional styles for the container.
   */
  sx?: SxProps<Theme>;

  /**
   * Additional styles for the content area.
   */
  contentSx?: SxProps<Theme>;
}

/**
 * A reusable expandable/collapsible section component.
 *
 * Supports both controlled and uncontrolled modes:
 * - Controlled: Pass `expanded` and `onToggle` props
 * - Uncontrolled: Component manages its own state, use `defaultExpanded`
 *
 * The `forceExpanded` prop can be used to force the section open
 * (e.g., when there are validation errors that need user attention).
 */
export const ExpandableSection: React.FC<ExpandableSectionProps> = ({
  children,
  expanded: controlledExpanded,
  onToggle,
  disabled = false,
  defaultExpanded = false,
  forceExpanded = false,
  title = "Additional Options",
  forceExpandedTitle = "Required Options (Error)",
  hasError = false,
  sx,
  contentSx,
}) => {
  // Internal state for uncontrolled mode
  const [internalExpanded, setInternalExpanded] = useState(defaultExpanded);

  // Determine if we're in controlled mode
  const isControlled = controlledExpanded !== undefined;

  // Calculate effective expanded state
  const isExpanded =
    forceExpanded || (isControlled ? controlledExpanded : internalExpanded);

  // Can the user toggle this section?
  const canToggle = !disabled && !forceExpanded;

  const handleToggle = useCallback(() => {
    if (!canToggle) return;

    const newExpanded = !isExpanded;

    if (isControlled) {
      onToggle?.(newExpanded);
    } else {
      setInternalExpanded(newExpanded);
      onToggle?.(newExpanded);
    }
  }, [canToggle, isExpanded, isControlled, onToggle]);

  // Don't render anything if there are no children
  if (!children) return null;

  const displayTitle = hasError ? forceExpandedTitle : title;

  return (
    <Collapse in={isExpanded} timeout={200}>
      <Stack
        sx={{
          px: 2,
          pb: 1,
          pt: 0,
          backgroundColor: hasError ? "error.lighter" : "background.paper",
          borderTop: "1px solid",
          borderTopColor: hasError ? "error.light" : "divider",
          borderBottomLeftRadius: "0.4rem",
          borderBottomRightRadius: "0.4rem",
          ...sx,
        }}
        spacing={0.5}
      >
        <Stack
          direction="row"
          alignItems="center"
          justifyContent="space-between"
          sx={{ cursor: canToggle ? "pointer" : "default" }}
          onClick={handleToggle}
        >
          <Typography
            variant="caption"
            color={hasError ? "error.main" : "text.secondary"}
            sx={{
              fontWeight: 500,
              textTransform: "uppercase",
              letterSpacing: 0.5,
              mb: 0.5,
            }}
          >
            {displayTitle}
          </Typography>

          {canToggle && (
            <IconButton
              size="small"
              sx={{
                p: 0,
                transition: "transform 0.2s ease-in-out",
                transform: isExpanded ? "rotate(90deg)" : "rotate(0deg)",
              }}
              aria-label={isExpanded ? "Collapse section" : "Expand section"}
            >
              <ChevronRightIcon fontSize="small" />
            </IconButton>
          )}
        </Stack>

        <Stack sx={contentSx}>{children}</Stack>
      </Stack>
    </Collapse>
  );
};

ExpandableSection.displayName = "ExpandableSection";

export default ExpandableSection;
