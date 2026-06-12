/**
 * Centralized sizing system for task interface widgets.
 *
 * DESIGN PRINCIPLE: Base widgets (TextField, Autocomplete) are full-width by default.
 * Container components (FieldRow, FieldContainer, etc.) control the width of their children.
 *
 * This separation of concerns means:
 * - Base widgets just render content, no width logic
 * - Containers/layouts are the single source of truth for sizing
 * - No duplication of width constraints across widgets
 */

/**
 * Container width constraints.
 * Use these when wrapping fields in containers to control their width.
 */
export const FIELD_SIZES = {
  /** Auto-size based on content (no fixed width) - for checkboxes, toggles */
  auto: 'auto',
  /** Very short - single digits, cell parameters (8rem / 128px) */
  xs: '8rem',
  /** Short - small numbers, short enums (12rem / 192px) */
  sm: '12rem',
  /** Medium - default for standalone fields (20rem / 320px) */
  md: '20rem',
  /** Long - file paths, longer text (32rem / 512px) */
  lg: '32rem',
  /** Full container width */
  full: '100%',
} as const;

export type FieldSize = keyof typeof FIELD_SIZES;

/**
 * Standard spacing values used across field components.
 * Ensures consistent margins between fields.
 */
export const FIELD_SPACING = {
  /** Vertical gap between rows of fields (MUI spacing units: 1 = 8px) */
  rowGap: 1,
  /** Horizontal gap between fields in a row */
  columnGap: 2,
  /** Left margin for fields (indentation from container edge) */
  marginLeft: 2,
} as const;

/**
 * Default styles for full-width base widgets.
 * Base widgets use these by default - containers override as needed.
 */
export const FULL_WIDTH_FIELD_STYLES = {
  width: '100%',
  minWidth: 0, // Allow shrinking in flex containers
} as const;

/**
 * Get CSS properties for constraining a container to a specific size.
 * Use this in container components (FieldRow, FieldContainer) to limit child width.
 *
 * @param size - The FieldSize to constrain to
 * @returns CSS properties object for the container
 */
export function getContainerSizeStyles(size: FieldSize) {
  if (size === 'auto') {
    return {
      width: 'auto',
    };
  }

  if (size === 'full') {
    return {
      width: '100%',
      flex: '1 1 0',
      minWidth: 0,
    };
  }

  // Constrained sizes: set maxWidth on the container
  return {
    width: '100%',
    maxWidth: FIELD_SIZES[size],
    minWidth: 0,
    flexShrink: 1,
  };
}
