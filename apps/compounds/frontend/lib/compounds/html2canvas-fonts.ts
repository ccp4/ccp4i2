/**
 * onclone hook for html2canvas captures.
 *
 * Production wraps Roboto via next/font/google, which resolves to a
 * CSS-variable family name (`__Roboto_<hash>__`, plus a metric-tuned
 * `__Roboto_<hash>_fallback__`). html2canvas reads that family name via
 * `getComputedStyle`, then can't resolve it to a real font when it
 * measures glyph widths on its own canvas. Result: character widths
 * diverge between DOM and capture, and adjacent words overlap in the
 * captured PNG ("OVCAR8 IC50" → "OVCARGI50" as the 8/space/IC glyphs
 * crash into each other; "1.28×10⁷" → "1.28e7"-shaped blob).
 *
 * A body-only !important override isn't enough because MUI emits
 * class-level font-family rules on Typography variants that beat
 * inheritance. The reliable fix is to walk every element in the clone
 * and rewrite its inline `font-family` to a system-resolvable stack —
 * keeping monospace where the element's *original* computed style was
 * monospace, sans-serif everywhere else.
 *
 * Cost: a single sync pass over the clone subtree. Capture path is
 * already user-initiated and visibly synchronous (html2canvas itself
 * dominates the timeline), so the extra walk is invisible.
 */

const CAPTURE_SANS_SERIF =
  'system-ui, -apple-system, "Segoe UI", "Helvetica Neue", Arial, sans-serif';

const CAPTURE_MONOSPACE =
  '"Roboto Mono", Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace';

const MONOSPACE_TOKENS = /Mono|Menlo|Monaco|Consolas|Courier|monospace/i;

export function pinCaptureFonts(clonedDoc: Document): void {
  const win = clonedDoc.defaultView;
  if (!win) return;
  // Body first, so any element that genuinely inherits picks up the
  // sans stack rather than the next/font wrapper.
  clonedDoc.body.style.setProperty('font-family', CAPTURE_SANS_SERIF, 'important');
  clonedDoc.body.querySelectorAll<HTMLElement>('*').forEach((el) => {
    const computed = win.getComputedStyle(el).fontFamily;
    const target = MONOSPACE_TOKENS.test(computed) ? CAPTURE_MONOSPACE : CAPTURE_SANS_SERIF;
    el.style.setProperty('font-family', target, 'important');
  });
}
