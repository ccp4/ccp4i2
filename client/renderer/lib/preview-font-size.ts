/**
 * The file preview's text size (#790). Ctrl/Cmd and + or - zoom the whole
 * window, but testers did not know that and asked for buttons. The chosen
 * size is kept per browser so it outlives the dialog.
 */
const FONT_SIZE_KEY = "ccp4i2.ui.previewFontSize";
const FONT_SIZES = [10, 11, 12, 13, 14, 16, 18, 20, 24, 28];
export const DEFAULT_FONT_SIZE = 14;

export function readFontSize(): number {
  if (typeof window === "undefined") return DEFAULT_FONT_SIZE;
  try {
    const stored = Number(window.localStorage.getItem(FONT_SIZE_KEY));
    return FONT_SIZES.includes(stored) ? stored : DEFAULT_FONT_SIZE;
  } catch {
    return DEFAULT_FONT_SIZE;
  }
}

export function writeFontSize(size: number): void {
  try {
    window.localStorage.setItem(FONT_SIZE_KEY, String(size));
  } catch {
    /* private mode or quota: the size still applies for this session */
  }
}

/** The next size up or down the scale, or the same one at either end. */
export function stepFontSize(current: number, direction: 1 | -1): number {
  const index = FONT_SIZES.indexOf(current);
  const from = index === -1 ? FONT_SIZES.indexOf(DEFAULT_FONT_SIZE) : index;
  return FONT_SIZES[from + direction] ?? current;
}
