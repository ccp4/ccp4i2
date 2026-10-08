/**
 * The campaign overview's site matrix: one column per campaign site, one box
 * per dataset at that site.
 *
 * Two independent things are shown in one box, because they answer different
 * questions and both matter when scanning a campaign:
 *
 * * **Fill** says whether PanDDA found an event at the site in this dataset.
 *   A filled box is an event; an outlined box is a verdict with no event (a
 *   hit PanDDA missed, or a site somebody checked by hand); nothing at all is
 *   neither.
 * * **Colour** says what a reviewer decided: grey not yet looked at, red
 *   rejected (empty), amber uncertain, green accepted (hit).
 *
 * So a grey filled box is work outstanding (an event nobody has judged), and
 * a green outlined box is a finding PanDDA did not make. A row whose model is
 * in a different frame from the campaign's sites cannot be placed against
 * them at all, and says so rather than showing boxes that would be wrong.
 *
 * Everything here is pure so it can be tested without rendering the table.
 */

import {
  CampaignSite,
  CurrentModelJob,
  MemberProjectWithSummary,
  SiteCell,
  SiteVerdict,
} from "../types/campaigns";
import { siteViewUrl } from "./site-verdicts";

/** How a dataset's box at one site is drawn. */
export type SiteCellKind =
  /** PanDDA found an event here. */
  | "filled"
  /** No event, but somebody gave a verdict here. */
  | "outlined"
  /** No event and no verdict. */
  | "none"
  /** The dataset's model is in another frame: no box can be placed. */
  | "mismatch";

/** MUI palette keys for the four verdict states. */
export type VerdictColourKey = "grey" | "error" | "warning" | "success";

/** Legend order: the decision states in the order a reviewer moves through. */
export const VERDICT_LEGEND: { verdict: SiteVerdict | null; label: string }[] =
  [
    { verdict: null, label: "Not evaluated" },
    { verdict: "hit", label: "Hit" },
    { verdict: "unclear", label: "Uncertain" },
    { verdict: "empty", label: "Rejected" },
  ];

/**
 * The sites that get a column, in the campaign's order.
 *
 * A site without a uuid cannot be matched to `site_cells` (which are keyed by
 * it), so it gets no column; that only happens against an older server, when
 * the overview falls back to the chips anyway.
 */
export function matrixSites(sites: CampaignSite[] | undefined): CampaignSite[] {
  return [...(sites ?? [])]
    .filter((site) => Boolean(site.uuid))
    .sort((a, b) => (a.order ?? 0) - (b.order ?? 0) || a.id - b.id);
}

/**
 * Whether the overview can draw the matrix at all: there are sites with
 * uuids, and the server sends per-site cells. Otherwise the verdict chips
 * column is used instead.
 */
export function canShowMatrix(
  sites: CampaignSite[] | undefined,
  projects: MemberProjectWithSummary[]
): boolean {
  return (
    matrixSites(sites).length > 0 &&
    projects.some((project) => project.site_cells !== undefined)
  );
}

/** One dataset's cell at one site, or null when the server sent none. */
export function cellFor(
  project: MemberProjectWithSummary,
  site: CampaignSite
): SiteCell | null {
  if (!site.uuid) return null;
  return project.site_cells?.[site.uuid] ?? null;
}

export function siteCellKind(
  cell: SiteCell | null,
  frameMismatch?: string | null
): SiteCellKind {
  if (frameMismatch) return "mismatch";
  if (!cell) return "none";
  if (cell.event) return "filled";
  if (cell.verdict) return "outlined";
  return "none";
}

export function verdictColourKey(
  verdict: SiteVerdict | null | undefined
): VerdictColourKey {
  switch (verdict) {
    case "hit":
      return "success";
    case "unclear":
      return "warning";
    case "empty":
      return "error";
    default:
      return "grey";
  }
}

export function verdictWords(verdict: SiteVerdict | null | undefined): string {
  switch (verdict) {
    case "hit":
      return "hit (accepted)";
    case "unclear":
      return "uncertain";
    case "empty":
      return "rejected (empty)";
    default:
      return "not evaluated";
  }
}

/** "event 2, hit probability 0.87, d = 1.4 Å" */
export function eventDescription(event: NonNullable<SiteCell["event"]>): string {
  const parts = [`event ${event.event_idx}`];
  if (event.hit_probability !== null && event.hit_probability !== undefined) {
    parts.push(`hit probability ${event.hit_probability.toFixed(2)}`);
  }
  parts.push(`d = ${event.distance.toFixed(1)} Å`);
  return parts.join(", ");
}

export const NO_MODEL_JOB_NOTE = "No refinement or dimple run yet";

/**
 * The tooltip for a box, as lines. For a mismatch it is the reason only; for
 * a blank cell it is still useful (the site name, "not evaluated") because a
 * hover over an empty column is how somebody checks which site it is.
 */
export function siteCellTooltip(
  siteName: string,
  cell: SiteCell | null,
  frameMismatch: string | null | undefined,
  hasModelJob: boolean
): string[] {
  if (frameMismatch) return [siteName, frameMismatch];
  const lines = [siteName, `Verdict: ${verdictWords(cell?.verdict)}`];
  if (cell?.event) {
    lines.push(`PanDDA ${eventDescription(cell.event)}`);
  } else {
    lines.push("No PanDDA event");
  }
  lines.push(hasModelJob ? "Click to view" : NO_MODEL_JOB_NOTE);
  return lines;
}

/**
 * Where clicking a box goes: the campaign's Moorhen view on this dataset's
 * current model, at this site (by the site's numeric id, which is what the
 * Moorhen page's `site` parameter takes). Null when there is no model to
 * open, or no campaign.
 */
export function siteCellUrl(
  campaignId: number | undefined,
  modelJob: CurrentModelJob | null | undefined,
  site: CampaignSite
): string | null {
  return siteViewUrl(campaignId, modelJob?.id, site.id);
}

/**
 * How many of the campaign's sites have a verdict in this dataset, counted
 * over the sites shown, so the count and the boxes beside it always agree.
 */
export function evaluatedCount(
  project: MemberProjectWithSummary,
  sites: CampaignSite[]
): { done: number; total: number; findings: number } {
  const verdicts = sites.map((site) => cellFor(project, site)?.verdict);
  return {
    done: verdicts.filter(Boolean).length,
    total: sites.length,
    /** Hits and uncertains: the verdicts that are something found. */
    findings: verdicts.filter((v) => v === "hit" || v === "unclear").length,
  };
}
