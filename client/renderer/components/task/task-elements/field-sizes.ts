/**
 * Centralized field sizing system for task interface widgets.
 *
 * This provides consistent, responsive widths for form fields across
 * all task interfaces. Fields use maxWidth constraints so they can
 * shrink on narrow screens while maintaining alignment on wider screens.
 */

export const FIELD_SIZES = {
  /** Auto-size based on content (no fixed width) - for checkboxes, toggles */
  auto: 'auto',
  /** Very short - single digits (8rem / 128px) */
  xs: '8rem',
  /** Short - small numbers, short enums (12rem / 192px) */
  sm: '12rem',
  /** Medium - default for most fields (20rem / 320px) */
  md: '20rem',
  /** Long - file paths, longer text, multi-option dropdowns (32rem / 512px) */
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
 * Infer the appropriate field size based on item class and qualifiers.
 *
 * This provides sensible defaults so most fields don't need explicit sizing.
 * Can be overridden by passing a `size` prop to the component.
 *
 * @param item - The task item object (contains _class, _baseClass)
 * @param qualifiers - Qualifiers object (contains enumerators, guiMode, etc.)
 * @returns The inferred FieldSize
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
  // (MUI TextField wraps checkbox+label in a box that needs consistent sizing)
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
 * Get the CSS properties for a given field size.
 * Uses maxWidth constraints so fields can shrink on narrow screens
 * while maintaining consistent sizing on wider screens.
 *
 * @param size - The FieldSize to get styles for
 * @returns CSS properties object
 */
export function getFieldSizeStyles(size: FieldSize) {
  // Auto-sizing: no width constraints, let content determine size
  if (size === 'auto') {
    return {
      width: 'auto',
      minWidth: 0, // Allow shrinking in flex containers
    };
  }

  // Full width: fill container
  if (size === 'full') {
    return {
      width: '100%',
      flexGrow: 1,
      minWidth: 0, // Allow shrinking in flex containers
    };
  }

  // Responsive sizes: fill available space up to maxWidth
  // This allows fields to shrink on narrow screens while
  // maintaining consistent max widths on wider screens
  const maxWidth = FIELD_SIZES[size];
  return {
    width: '100%',
    maxWidth,
    minWidth: 0, // Allow shrinking below content width in flex containers
    flexShrink: 1, // Allow shrinking when container is narrow
  };
}
