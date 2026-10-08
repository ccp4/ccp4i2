// @vitest-environment node
/**
 * The rules behind the campaign overview's site matrix: which box a dataset
 * gets at a site, its colour, what hovering says, and where a click goes.
 */

import { describe, it, expect } from "vitest";
import {
  NO_MODEL_JOB_NOTE,
  canShowMatrix,
  cellFor,
  evaluatedCount,
  eventDescription,
  matrixSites,
  siteCellKind,
  siteCellTooltip,
  siteCellUrl,
  verdictColourKey,
  verdictWords,
} from "../lib/site-matrix";
import {
  CampaignSite,
  MemberProjectWithSummary,
  SiteCell,
} from "../types/campaigns";

const site = (id: number, uuid?: string, order = id): CampaignSite => ({
  id,
  uuid,
  name: `Site ${id}`,
  origin: [0, 0, 0],
  order,
});

const event = {
  event_idx: 2,
  hit_probability: 0.874,
  distance: 1.42,
  has_pose: true,
};

const project = (
  overrides: Partial<MemberProjectWithSummary> = {}
): MemberProjectWithSummary =>
  ({ id: 1, name: "ds", jobs: [], kpis: {}, ...overrides } as MemberProjectWithSummary);

describe("siteCellKind", () => {
  it("fills a box for an event, whatever the verdict", () => {
    expect(siteCellKind({ event, verdict: null })).toBe("filled");
    expect(siteCellKind({ event, verdict: "empty" })).toBe("filled");
  });

  it("outlines a verdict with no event (a hit PanDDA missed)", () => {
    expect(siteCellKind({ event: null, verdict: "hit" })).toBe("outlined");
    expect(siteCellKind({ event: null, verdict: "empty" })).toBe("outlined");
  });

  it("draws nothing for neither, or for a missing cell", () => {
    expect(siteCellKind({ event: null, verdict: null })).toBe("none");
    expect(siteCellKind(null)).toBe("none");
  });

  it("gives way to a frame mismatch", () => {
    expect(siteCellKind({ event, verdict: "hit" }, "different cell")).toBe(
      "mismatch"
    );
    expect(siteCellKind(null, "different cell")).toBe("mismatch");
  });
});

describe("verdict colours and words", () => {
  it("maps each verdict to a palette key", () => {
    expect(verdictColourKey(null)).toBe("grey");
    expect(verdictColourKey(undefined)).toBe("grey");
    expect(verdictColourKey("empty")).toBe("error");
    expect(verdictColourKey("unclear")).toBe("warning");
    expect(verdictColourKey("hit")).toBe("success");
  });

  it("says each verdict in words", () => {
    expect(verdictWords(null)).toBe("not evaluated");
    expect(verdictWords("empty")).toBe("rejected (empty)");
    expect(verdictWords("unclear")).toBe("uncertain");
    expect(verdictWords("hit")).toBe("hit (accepted)");
  });
});

describe("siteCellTooltip", () => {
  it("describes an event", () => {
    expect(eventDescription(event)).toBe(
      "event 2, hit probability 0.87, d = 1.4 Å"
    );
    const lines = siteCellTooltip(
      "Pocket A",
      { event, verdict: "hit" },
      null,
      true
    );
    expect(lines).toEqual([
      "Pocket A",
      "Verdict: hit (accepted)",
      "PanDDA event 2, hit probability 0.87, d = 1.4 Å",
      "Click to view",
    ]);
  });

  it("leaves out an unknown hit probability", () => {
    expect(eventDescription({ ...event, hit_probability: null })).toBe(
      "event 2, d = 1.4 Å"
    );
  });

  it("says when there is no model to open", () => {
    const lines = siteCellTooltip(
      "Pocket A",
      { event: null, verdict: "unclear" },
      null,
      false
    );
    expect(lines).toContain("No PanDDA event");
    expect(lines).toContain("Verdict: uncertain");
    expect(lines[lines.length - 1]).toBe(NO_MODEL_JOB_NOTE);
  });

  it("gives only the reason for a frame mismatch", () => {
    expect(
      siteCellTooltip("Pocket A", null, "Cell differs by 4%", true)
    ).toEqual(["Pocket A", "Cell differs by 4%"]);
  });
});

describe("siteCellUrl", () => {
  const model = { id: 42, uuid: "j", number: "3", task_name: "i2Dimple" };

  it("opens the current model at the site, by the site's numeric id", () => {
    expect(siteCellUrl(7, model, site(5, "u5"))).toBe(
      "/ccp4i2/moorhen-page/campaign/7?job=42&site=5"
    );
  });

  it("is not a link without a model job or a campaign", () => {
    expect(siteCellUrl(7, null, site(5, "u5"))).toBeNull();
    expect(siteCellUrl(undefined, model, site(5, "u5"))).toBeNull();
  });
});

describe("matrix columns", () => {
  it("orders sites and drops those without a uuid", () => {
    const sites = [site(1, "a", 2), site(2, undefined, 0), site(3, "c", 1)];
    expect(matrixSites(sites).map((s) => s.id)).toEqual([3, 1]);
  });

  it("falls back to chips when the server sends no cells", () => {
    const sites = [site(1, "a")];
    expect(canShowMatrix(sites, [project()])).toBe(false);
    expect(canShowMatrix(sites, [project({ site_cells: {} })])).toBe(true);
    expect(canShowMatrix([], [project({ site_cells: {} })])).toBe(false);
    expect(canShowMatrix([site(1)], [project({ site_cells: {} })])).toBe(
      false
    );
  });

  it("looks cells up by uuid and counts verdicts", () => {
    const cells: Record<string, SiteCell> = {
      a: { event, verdict: null },
      b: { event: null, verdict: "empty" },
      c: { event: null, verdict: "hit" },
    };
    const p = project({ site_cells: cells });
    const sites = [site(1, "a"), site(2, "b"), site(3, "c"), site(4, "d")];
    expect(cellFor(p, sites[0])).toBe(cells.a);
    expect(cellFor(p, sites[3])).toBeNull();
    expect(evaluatedCount(p, sites)).toEqual({
      done: 2,
      total: 4,
      findings: 1,
    });
  });
});
