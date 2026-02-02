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

// =============================================================================
// LEGACY COMPATIBILITY - Remove after migration
// =============================================================================

/**
 * @deprecated Use FULL_WIDTH_FIELD_STYLES for base widgets,
 * or getContainerSizeStyles for containers.
 *
 * This function is kept for backwards compatibility during migration.
 */
export function inferFieldSize(item: any, qualifiers: any): FieldSize {
  // File selectors always get full width
  if (item?._baseClass === 'CDataFile') {
    return 'full';
  }

  // Lists get full width
  if (item?._class === 'CList' || item?._baseClass === 'CList') {
    return 'full';
  }

  // Multiline text gets full width
  if (qualifiers?.guiMode === 'multiLine') {
    return 'full';
  }

  // Sequence strings are typically long
  if (item?._class === 'CSequenceString') {
    return 'lg';
  }

  // File paths need more space
  if (item?._class === 'CFilePath') {
    return 'lg';
  }

  // Enumerators - size based on option count
  const enumCount = qualifiers?.enumerators?.length || 0;
  if (enumCount > 0) {
    // Radio buttons with few options can be smaller
    if (qualifiers?.guiMode === 'radio' || qualifiers?.guiMode === 'multiLineRadio') {
      return enumCount <= 4 ? 'md' : 'lg';
    }
    // Autocomplete/dropdown
    if (enumCount <= 3) return 'sm';
    if (enumCount <= 8) return 'md';
    return 'lg';
  }

  // Numeric types are typically short
  if (item?._class === 'CInt' || item?._class === 'CFloat') {
    return 'md';
  }

  // Boolean checkboxes use medium size to match other fields
  if (item?._class === 'CBoolean') {
    return 'md';
  }

  // Cell parameters
  if (item?._class === 'CCellLength' || item?._class === 'CCellAngle') {
    return 'xs';
  }

  // Default to medium
  return 'md';
}

/**
 * @deprecated Use FULL_WIDTH_FIELD_STYLES for base widgets,
 * or getContainerSizeStyles for containers.
 */
export function getFieldSizeStyles(size: FieldSize) {
  if (size === 'auto') {
    return {
      width: 'auto',
      minWidth: 0,
    };
  }

  if (size === 'full') {
    return {
      width: '100%',
      flexGrow: 1,
      minWidth: 0,
    };
  }

  const maxWidth = FIELD_SIZES[size];
  return {
    width: '100%',
    maxWidth,
    minWidth: 0,
    flexShrink: 1,
  };
}
