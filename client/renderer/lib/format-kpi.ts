/**
 * A KPI value for a chip. Every numeric KPI reaches the client through the
 * float table, counts included (the gleaner stores a CInt as a float), so an
 * integral value is printed as the integer it is; anything else to three
 * significant figures.
 */
export function formatKpiValue(value: number): string {
  if (!Number.isFinite(value)) return String(value);
  if (Number.isInteger(value)) return String(value);
  return value.toPrecision(3);
}
