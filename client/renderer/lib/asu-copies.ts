/**
 * The copy numbers in the Define AU contents table, and the Matthews table
 * beside it (#687).
 *
 * Matthews is run on the table's smallest whole set: the copies divided by
 * their common factor. A table holding 2 A + 2 B is analysed as multiples of
 * 1 A + 1 B, so its rows read 1x, 2x, 3x ... and stay put when the copies are
 * changed; the row for the copies as they stand is marked "current", and
 * clicking another row rescales every entry. Earlier the analysis was of the
 * whole contents, so after choosing "2x" the table recomputed and the row
 * the user had clicked became "1x" (Paul's objection in #687).
 */

export interface MatthewsRow {
  nmol_in_asu: number;
  percent_solvent: number;
  prob_matth: number;
}

function gcd(a: number, b: number): number {
  a = Math.abs(Math.trunc(a));
  b = Math.abs(Math.trunc(b));
  while (b) [a, b] = [b, a % b];
  return a;
}

/** Copies as the table rows hold them; an unset or non-positive value is 1. */
export function copiesOfRows(rows: any[]): number[] {
  return rows.map((row) => {
    const n = Number(row?._value?.nCopies?._value);
    return Number.isFinite(n) && n > 0 ? Math.trunc(n) : 1;
  });
}

/** The common factor of the copies: 1 for an empty table. */
export function commonFactor(copies: number[]): number {
  if (copies.length === 0) return 1;
  return copies.reduce((acc, n) => gcd(acc, n), 0) || 1;
}

/** The copies after choosing `multiple` of the smallest whole set. */
export function scaledCopies(copies: number[], multiple: number): number[] {
  const factor = commonFactor(copies);
  return copies.map((n) => (n / factor) * multiple);
}

/** The most probable row; null for no rows. */
export function bestMatthewsRow(rows: MatthewsRow[] | null | undefined): MatthewsRow | null {
  if (!rows || rows.length === 0) return null;
  return rows.reduce((best, row) => (row.prob_matth > best.prob_matth ? row : best), rows[0]);
}

/** The entry names as the table rows hold them; an unnamed row is "Sequence N". */
export function namesOfRows(rows: any[]): string[] {
  return rows.map((row, index) => String(row?._value?.name?._value || `Sequence ${index + 1}`));
}

/**
 * What choosing `multiple` of the smallest whole set means for each entry,
 * spelled out ("Chain_A_C × 2, Chain_B_D × 2"), so a Matthews row says what
 * the table would hold and not only a multiplier (#687).
 */
export function describeCopies(names: string[], copies: number[], multiple: number): string {
  const scaled = scaledCopies(copies, multiple);
  return names.map((name, i) => `${name} × ${scaled[i] ?? multiple}`).join(", ");
}
