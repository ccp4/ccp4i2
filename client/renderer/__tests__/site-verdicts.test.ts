/**
 * The rules behind a dataset's Sites cell.
 *
 * These pin the two judgements that make the column readable: which verdicts
 * earn a chip (and in what order, so hits are never the ones dropped), and
 * when the completeness count is worth showing at all.
 */

import { describe, it, expect } from "vitest";
import {
  evaluationProgress,
  preferredJobId,
  siteViewUrl,
  sortVerdicts,
  verdictChips,
} from "../lib/site-verdicts";
import { CampaignJobInfo, SiteEvaluationSummary } from "../types/campaigns";

const hit = (id: number, name: string): SiteEvaluationSummary => ({
  site_id: id,
  site_name: name,
  verdict: "hit",
});
const unclear = (id: number, name: string): SiteEvaluationSummary => ({
  site_id: id,
  site_name: name,
  verdict: "unclear",
});

describe("ordering", () => {
  it("puts hits before unclears", () => {
    const ordered = sortVerdicts([unclear(1, "A"), hit(2, "Z")]);
    expect(ordered.map((e) => e.site_id)).toEqual([2, 1]);
  });

  it("is alphabetical within a verdict", () => {
    const ordered = sortVerdicts([hit(1, "Pocket B"), hit(2, "Pocket A")]);
    expect(ordered.map((e) => e.site_name)).toEqual(["Pocket A", "Pocket B"]);
  });

  it("does not mutate what it is given", () => {
    const input = [unclear(1, "A"), hit(2, "Z")];
    sortVerdicts(input);
    expect(input.map((e) => e.site_id)).toEqual([1, 2]);
  });
});

describe("which chips a row shows", () => {
  it("caps the list", () => {
    const { shown, overflow } = verdictChips(
      [hit(1, "A"), hit(2, "B"), hit(3, "C"), hit(4, "D")],
      3
    );
    expect(shown).toHaveLength(3);
    expect(overflow).toBe(1);
  });

  it("never drops a hit in favour of an unclear", () => {
    // The point of ordering before capping: a dataset with many unclears
    // would otherwise hide the one result somebody is scanning for.
    const { shown } = verdictChips(
      [
        unclear(1, "A"),
        unclear(2, "B"),
        unclear(3, "C"),
        unclear(4, "D"),
        hit(5, "Z"),
      ],
      3
    );
    expect(shown.some((e) => e.verdict === "hit")).toBe(true);
  });

  it("reports no overflow when everything fits", () => {
    expect(verdictChips([hit(1, "A")], 3)).toEqual({
      shown: [hit(1, "A")],
      overflow: 0,
    });
  });

  it("handles a dataset with no verdicts at all", () => {
    expect(verdictChips(undefined)).toEqual({ shown: [], overflow: 0 });
    expect(verdictChips([])).toEqual({ shown: [], overflow: 0 });
  });
});

describe("the completeness count", () => {
  it("shows a fraction while a dataset is part-evaluated", () => {
    expect(evaluationProgress(12, 40)).toBe("12/40");
  });

  it("is silent once every site has been looked at", () => {
    // "40/40" on every finished row is noise; the eye is looking for the
    // rows that still carry work.
    expect(evaluationProgress(40, 40)).toBeNull();
  });

  it("is silent when nobody has looked yet", () => {
    expect(evaluationProgress(0, 40)).toBeNull();
    expect(evaluationProgress(undefined, 40)).toBeNull();
  });

  it("is silent when the campaign has no sites", () => {
    expect(evaluationProgress(0, 0)).toBeNull();
    expect(evaluationProgress(3, undefined)).toBeNull();
  });

  it("does not exceed the denominator when data disagrees", () => {
    expect(evaluationProgress(41, 40)).toBeNull();
  });
});

describe("where a chip navigates", () => {
  it("opens the campaign view on that job and site", () => {
    expect(siteViewUrl(7, 42, 3)).toBe(
      "/ccp4i2/moorhen-page/campaign/7?job=42&site=3"
    );
  });

  it("is not a link when the dataset has no job", () => {
    // Nothing to show at a site without a job, so the chip does not pretend.
    expect(siteViewUrl(7, undefined, 3)).toBeNull();
  });

  it("is not a link outside a campaign", () => {
    expect(siteViewUrl(undefined, 42, 3)).toBeNull();
  });
});

describe("which job a chip opens", () => {
  const job = (
    id: number,
    number: string,
    status: number
  ): CampaignJobInfo =>
    ({ id, number, status, task_name: "servalcat_pipe" } as CampaignJobInfo);

  it("prefers the latest finished job", () => {
    expect(
      preferredJobId([job(1, "1", 6), job(2, "2", 6), job(3, "3", 5)])
    ).toBe(2);
  });

  it("falls back to the latest job when none has finished", () => {
    // A running or failed dataset still opens on something.
    expect(preferredJobId([job(1, "1", 5), job(2, "2", 3)])).toBe(2);
  });

  it("ignores sub-jobs", () => {
    // Steps of a pipeline, not the dataset's result.
    expect(preferredJobId([job(1, "1", 6), job(9, "1.2", 6)])).toBe(1);
  });

  it("orders numerically, not as text", () => {
    expect(preferredJobId([job(1, "9", 6), job(2, "10", 6)])).toBe(2);
  });

  it("has nothing to open when the dataset has no jobs", () => {
    expect(preferredJobId([])).toBeUndefined();
    expect(preferredJobId(undefined)).toBeUndefined();
  });
});
