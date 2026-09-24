/**
 * What the campaign overview shows in a dataset's Sites cell.
 *
 * A row shows chips for the sites where something was found, and a one-line
 * note saying how far through the dataset the evaluation has got. Both are
 * deliberately partial views:
 *
 * * Only `hit` and `unclear` get a chip. A rich campaign has 30-40 sites and
 *   most are empty for most datasets, so chipping every verdict would fill the
 *   column with the least interesting one. `empty` renders as no chip.
 * * The note carries the rest. Chips alone cannot separate a dataset examined
 *   throughout and found empty from one nobody has opened — both have no hits
 *   and no unclears — and those are opposite states: one is a finished result,
 *   the other is work outstanding. So the note speaks whenever the chips do
 *   not already answer the question, and stays quiet only when they do (a
 *   fully evaluated dataset that has chips) or when the campaign has no sites
 *   to evaluate against.
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

/** How loudly a dataset's note should be rendered. */
export type EvaluationTone =
  /** Part-way through: work outstanding. */
  | "progress"
  /** Finished, and nothing was found. */
  | "clear"
  /** Nobody has looked yet. */
  | "untouched";

export interface EvaluationNote {
  /** The caption under the chips. */
  text: string;
  tone: EvaluationTone;
  /** The longer form, for the tooltip. */
  detail: string;
}

/**
 * The note under a dataset's chips, or null when there is nothing to say.
 *
 * The denominator is the campaign's site count rather than anything
 * per-dataset, so the fractions are comparable down the column.
 *
 * Silent in two cases. A campaign with no sites has nothing to be part-way
 * through. And a fully evaluated dataset that has chips is already legible
 * from the chips alone — "40/40" on those rows is the noise that would make
 * the notes on the other rows hard to pick out.
 */
export function evaluationNote(
  evaluated: number | undefined,
  total: number | undefined,
  hasVerdicts: boolean
): EvaluationNote | null {
  if (!total || total <= 0) return null;

  // Clamped, because a row is better described as finished than as silent if
  // the two counts ever disagree.
  const done = Math.min(Math.max(evaluated ?? 0, 0), total);

  if (done <= 0) {
    return {
      text: "not evaluated",
      tone: "untouched",
      detail: `None of the campaign's ${total} sites has been looked at in this dataset.`,
    };
  }

  if (done < total) {
    return {
      text: `${done}/${total}`,
      tone: "progress",
      detail: `${done} of ${total} sites evaluated \u2014 ${
        total - done
      } still to look at.`,
    };
  }

  if (hasVerdicts) return null;

  return {
    text: `all ${total} empty`,
    tone: "clear",
    detail: `All ${total} sites evaluated in this dataset, and nothing was found at any of them.`,
  };
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
