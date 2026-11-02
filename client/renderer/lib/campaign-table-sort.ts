/**
 * Sorting the campaign table's rows by a column.
 *
 * Name sorts alphabetically, numbers within names by value, so "x0010" comes
 * after "x0009" and "x10" after "x9". Resolution and the R factors sort by
 * value, best (smallest) first. Eval. sorts by how many sites have a verdict.
 *
 * A site column sorts by what the dataset holds at that site, most
 * interesting first:
 *
 *   1. confirmed hit, PanDDA event       5. PanDDA event, not yet judged
 *   2. confirmed hit, no event           6. rejected, PanDDA event
 *   3. possible hit, PanDDA event        7. rejected, no event
 *   4. possible hit, no event            8. nothing here (no event, no verdict)
 *                                        9. model in another frame: no box
 *
 * "Possible" is the cell's effective verdict "unclear", which also covers
 * reviewers who disagree (lib/site-verdicts).
 *
 * Whichever way a column sorts, a dataset with no value (no refinement yet,
 * so no resolution) goes last, and ties keep name order. The rows are sorted
 * before the virtualiser sees them, so sorting and windowing do not interact;
 * the table keys the virtualiser by project id so measured row heights follow
 * their rows.
 */
import type { CampaignSite, MemberProjectWithSummary, SiteCell } from "../types/campaigns";
import { cellFor } from "./site-matrix";

export type SortColumn =
  | { kind: "name" }
  | { kind: "kpi"; key: "highResLimit" | "RFactor" | "RFree" }
  | { kind: "evaluated" }
  | { kind: "site"; siteUuid: string };

export type SortDirection = "asc" | "desc";

export interface SortState {
  column: SortColumn;
  direction: SortDirection;
}

export const sameColumn = (a: SortColumn, b: SortColumn): boolean =>
  a.kind === b.kind &&
  (a.kind !== "kpi" || a.key === (b as typeof a).key) &&
  (a.kind !== "site" || a.siteUuid === (b as typeof a).siteUuid);

/** A column's first direction: best first, which for every column here is "asc". */
export const firstDirection = (_column: SortColumn): SortDirection => "asc";

/** Clicking a heading: a new column sorts its first way; the same column
 *  turns round; a third click goes back to the table's own order (null). */
export function nextSort(current: SortState | null, column: SortColumn): SortState | null {
  if (!current || !sameColumn(current.column, column)) {
    return { column, direction: firstDirection(column) };
  }
  if (current.direction === firstDirection(column)) {
    return { column, direction: current.direction === "asc" ? "desc" : "asc" };
  }
  return null;
}

/** Where a dataset's cell at a site ranks: 1 most interesting, 9 least. */
export function siteRank(cell: SiteCell | null, frameMismatch?: string | null): number {
  if (frameMismatch) return 9;
  const event = Boolean(cell?.event);
  switch (cell?.verdict ?? null) {
    case "hit":
      return event ? 1 : 2;
    case "unclear":
      return event ? 3 : 4;
    case null:
      return event ? 5 : 8;
    case "empty":
      return event ? 6 : 7;
    default:
      return 8;
  }
}

const nameCollator = new Intl.Collator(undefined, { numeric: true, sensitivity: "base" });

const kpiValue = (project: MemberProjectWithSummary, key: string): number | null => {
  const raw = (project.kpis as Record<string, unknown> | undefined)?.[key];
  const value = typeof raw === "number" ? raw : parseFloat(String(raw ?? ""));
  return Number.isFinite(value) ? value : null;
};

const evaluatedValue = (project: MemberProjectWithSummary, sites: CampaignSite[]): number =>
  sites.filter((site) => Boolean(cellFor(project, site)?.verdict)).length;

/** The rows in the order `sort` asks for; `projects` itself when null. */
export function sortProjects(
  projects: MemberProjectWithSummary[],
  sort: SortState | null,
  sites: CampaignSite[]
): MemberProjectWithSummary[] {
  if (!sort) return projects;
  const sign = sort.direction === "asc" ? 1 : -1;
  const byName = (a: MemberProjectWithSummary, b: MemberProjectWithSummary) =>
    nameCollator.compare(a.name ?? "", b.name ?? "");

  // A sort key per row, or null for "no value": those go last either way.
  const { column } = sort;
  let key: (project: MemberProjectWithSummary) => number | null;
  if (column.kind === "name") {
    return [...projects].sort((a, b) => sign * byName(a, b));
  } else if (column.kind === "kpi") {
    key = (project) => kpiValue(project, column.key);
  } else if (column.kind === "evaluated") {
    // More sites judged is further on: that is "best first"
    key = (project) => -evaluatedValue(project, sites);
  } else {
    const site = sites.find((s) => s.uuid === column.siteUuid);
    key = (project) => (site ? siteRank(cellFor(project, site), project.frame_mismatch) : null);
  }
  return [...projects].sort((a, b) => {
    const ka = key(a);
    const kb = key(b);
    if (ka === null && kb === null) return byName(a, b);
    if (ka === null) return 1;
    if (kb === null) return -1;
    return sign * (ka - kb) || byName(a, b);
  });
}
