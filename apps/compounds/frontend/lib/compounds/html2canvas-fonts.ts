/**
 * onclone hook for html2canvas captures.
 *
 * Production wraps Roboto via next/font/google, which resolves to a
 * CSS-variable family name (`__Roboto_<hash>__`, plus a metric-tuned
 * `__Roboto_<hash>_fallback__`). html2canvas reads that family name via
 * `getComputedStyle`, then can't resolve it to a real font when it
 * measures glyph widths on its own canvas. Result: character widths
 * diverge between DOM and capture, and adjacent words overlap in the
 * captured PNG ("Lipinski compliance" → "Lipinskicompliance",
 * "OVCAR8 IC50" → "OVCAR8IC50" with a strikethrough-looking glyph
 * collision).
 *
 * Fix: in the `onclone` hook (which runs against the off-screen clone
 * html2canvas builds before measuring), pin the body's font-family to a
 * pure system stack. CSS inheritance carries it to descendants; child
 * elements that set their own fontFamily explicitly (e.g. the monospace
 * value columns via inline `style={{ fontFamily: MONOSPACE_FONT_STACK }}`)
 * keep their override because inline declarations beat parent inheritance.
 */

const CAPTURE_SANS_SERIF =
  'system-ui, -apple-system, "Segoe UI", "Helvetica Neue", Arial, sans-serif';

export function pinCaptureFonts(clonedDoc: Document): void {
  clonedDoc.body.style.setProperty('font-family', CAPTURE_SANS_SERIF, 'important');
}
