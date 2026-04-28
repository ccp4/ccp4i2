/**
 * Shared categorisation palette for scatter view. One source of truth so
 * the chip-strip swatch matches the point colour for the same selection.
 *
 * The palette is colour-blind-friendly (distinct hues at multiple
 * lightness levels) and 8 entries deep — fewer than 8 simultaneous
 * categories means we never wrap.
 */

export const SCATTER_CATEGORY_COLOURS: ReadonlyArray<{ fill: string; border: string }> = [
  { fill: 'rgba(231, 76,  60,  0.7)',  border: 'rgba(231, 76,  60,  1)' },
  { fill: 'rgba( 46, 134, 193, 0.7)',  border: 'rgba( 46, 134, 193, 1)' },
  { fill: 'rgba(241, 196,  15, 0.85)', border: 'rgba(212, 172,  13, 1)' },
  { fill: 'rgba( 39, 174,  96, 0.7)',  border: 'rgba( 39, 174,  96, 1)' },
  { fill: 'rgba(155,  89, 182, 0.7)',  border: 'rgba(155,  89, 182, 1)' },
  { fill: 'rgba(230, 126,  34, 0.7)',  border: 'rgba(230, 126,  34, 1)' },
  { fill: 'rgba( 26, 188, 156, 0.7)',  border: 'rgba( 26, 188, 156, 1)' },
  { fill: 'rgba(231,  76, 153, 0.7)',  border: 'rgba(231,  76, 153, 1)' },
];
