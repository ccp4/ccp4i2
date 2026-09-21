/**
 * What the campaign overview shows in a dataset's Sites cell.
 *
 * A row shows chips for the sites where something was found, and how much of
 * the dataset is still unlooked at. Both are deliberately partial views:
 *
 * * Only `hit` and `unclear` get a chip. A rich campaign has 30-40 sites and
 *   most are empty for most datasets, so chipping every verdict would fill the
 *   column with the least interesting one. `empty` and "nobody has looked"
 *   both render as nothing here — they are distinguished in the per-dataset
 *   view, and in the count below.
 * * The count appears only while the dataset is incomplete. "40/40" on every
 *   finished row is noise; what the eye is looking for down the column is the
 *   rows still carrying work.
 */

import { CampaignJobInfo, SiteEvaluationSummary } from "../types/campaigns";

/** Hits first, then unclears; alphabetical within each. */
export function sortVerdicts(
  evaluations: SiteEvaluationSummary[]
): SiteEvaluationSummary[] {
  return [...evaluations].sort((a, b) => {
    if ((a.verdict === "hit") !== (b.verdict === "hit")) {
      return a.verdict === "hit" ? -1 : 1;
    }
    return a.site_name.localeCompare(b.site_name);
  });
}

export interface VerdictChips {
  /** The chips to render, hits first. */
  shown: SiteEvaluationSummary[];
  /** How many were held back by the cap. */
  overflow: number;
}

/**
 * The chips a row shows, capped so one busy dataset cannot stretch the column.
 *
 * Hits are ordered first so that a dataset with many unclears never has its
 * hits pushed out of the cap — the hit is the thing somebody is scanning for.
 * The server sends them in this order already; re-sorting here means the
 * guarantee does not depend on that.
 */
export function verdictChips(
  evaluations: SiteEvaluationSummary[] | undefined,
  maxChips = 3
): VerdictChips {
  const ordered = sortVerdicts(evaluations ?? []);
  return {
    shown: ordered.slice(0, maxChips),
    overflow: Math.max(0, ordered.length - maxChips),
  };
}

/**
 * "12/40" while a dataset is part-evaluated, or null when there is nothing
 * useful to say — all sites looked at, none looked at, or no sites at all.
 *
 * The denominator is the campaign's site count rather than anything
 * per-dataset, so the fractions are comparable down the column.
 */
export function evaluationProgress(
  evaluated: number | undefined,
  total: number | undefined
): string | null {
  if (!total || total <= 0) return null;
  const done = evaluated ?? 0;
  if (done <= 0) return null;
  if (done >= total) return null;
  return `${done}/${total}`;
}

/**
 * Where a verdict chip navigates: the campaign's Moorhen view, opened on this
 * dataset's job at this site.
 *
 * Returns null when there is no job to open — a dataset with no job has
 * nothing to show at a site, so its chips are not links.
 */
export function siteViewUrl(
  campaignId: number | undefined,
  jobId: number | undefined,
  siteId: number
): string | null {
  if (!campaignId || !jobId) return null;
  return `/ccp4i2/moorhen-page/campaign/${campaignId}?job=${jobId}&site=${siteId}`;
}

/** Job status 6 is "Finished" (see STATUS_COLORS in the overview table). */
const FINISHED = 6;

/**
 * The job a site chip should open.
 *
 * The latest finished top-level job, because that is the refined result
 * somebody clicking a verdict wants to look at; the latest top-level job of
 * any status if none has finished, so a running or failed dataset still opens
 * on something rather than nothing. Sub-jobs are skipped — they are steps of
 * a pipeline, not the dataset's result.
 */
export function preferredJobId(
  jobs: CampaignJobInfo[] | undefined
): number | undefined {
  const topLevel = (jobs ?? []).filter((job) => !job.number.includes("."));
  if (topLevel.length === 0) return undefined;

  const byNumber = [...topLevel].sort(
    (a, b) => Number(a.number) - Number(b.number)
  );
  const finished = byNumber.filter((job) => job.status === FINISHED);
  const chosen = finished.length ? finished : byNumber;
  return chosen[chosen.length - 1].id;
}
