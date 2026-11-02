/**
 * The campaign table sorts by a clicked column (lib/campaign-table-sort).
 * A site column ranks what each dataset holds there, most interesting first:
 * confirmed PanDDA hit, confirmed hit without an event, possible hit with and
 * without an event, an event nobody has judged, rejected with and without an
 * event, then nothing there, then a model in another frame.
 */
import { describe, expect, it } from "vitest";
import {
  nextSort,
  siteRank,
  sortProjects,
} from "../lib/campaign-table-sort";
import type {
  CampaignSite,
  MemberProjectWithSummary,
  SiteCell,
  SiteVerdict,
} from "../types/campaigns";

const SITE = { id: 1, uuid: "s1", name: "Site 1", order: 1 } as CampaignSite;
const EVENT = { event_idx: 1, hit_probability: 0.9, distance: 1.2, has_pose: true };

const cell = (verdict: SiteVerdict | null, event: boolean): SiteCell => ({
  verdict,
  event: event ? EVENT : null,
});

const project = (
  name: string,
  extra: Partial<MemberProjectWithSummary> = {}
): MemberProjectWithSummary =>
  ({ id: name.length * 1000 + name.charCodeAt(name.length - 1), name, kpis: {}, ...extra }) as any;

const names = (rows: MemberProjectWithSummary[]) => rows.map((p) => p.name);

describe("siteRank", () => {
  it("orders every combination of verdict and event", () => {
    const order = [
      cell("hit", true),
      cell("hit", false),
      cell("unclear", true),
      cell("unclear", false),
      cell(null, true),
      cell("empty", true),
      cell("empty", false),
      cell(null, false),
    ].map((c) => siteRank(c));
    expect(order).toEqual([1, 2, 3, 4, 5, 6, 7, 8]);
  });

  it("puts no cell with an empty one, and another frame after both", () => {
    expect(siteRank(null)).toBe(8);
    expect(siteRank(cell("hit", true), "different cell")).toBe(9);
  });
});

describe("sortProjects", () => {
  it("sorts names alphabetically, numbers by value", () => {
    const rows = ["x0010", "x0002", "x0009", "X0001"].map((n) => project(n));
    expect(names(sortProjects(rows, { column: { kind: "name" }, direction: "asc" }, []))).toEqual([
      "X0001",
      "x0002",
      "x0009",
      "x0010",
    ]);
  });

  it("sorts a site column most interesting first, ties by name", () => {
    const rows = [
      project("x05", { site_cells: { s1: cell("empty", false) } }),
      project("x02", { site_cells: { s1: cell("hit", false) } }),
      project("x04", { site_cells: { s1: cell(null, true) } }),
      project("x01", { site_cells: { s1: cell("hit", true) } }),
      project("x06", {}),
      project("x03", { site_cells: { s1: cell("unclear", true) } }),
      project("x00", { site_cells: { s1: cell("hit", true) } }),
    ];
    const sorted = sortProjects(
      rows,
      { column: { kind: "site", siteUuid: "s1" }, direction: "asc" },
      [SITE]
    );
    expect(names(sorted)).toEqual(["x00", "x01", "x02", "x03", "x04", "x05", "x06"]);
  });

  it("sorts resolution best first, with no value last whichever way", () => {
    const rows = [
      project("a", { kpis: { highResLimit: 2.1 } }),
      project("b", { kpis: {} }),
      project("c", { kpis: { highResLimit: 1.4 } }),
    ];
    const column = { kind: "kpi", key: "highResLimit" } as const;
    expect(names(sortProjects(rows, { column, direction: "asc" }, []))).toEqual(["c", "a", "b"]);
    expect(names(sortProjects(rows, { column, direction: "desc" }, []))).toEqual(["a", "c", "b"]);
  });

  it("keeps the given order when not sorted", () => {
    const rows = [project("b"), project("a")];
    expect(sortProjects(rows, null, [])).toBe(rows);
  });
});

describe("nextSort", () => {
  it("goes best first, then the other way, then back to the table's own order", () => {
    const column = { kind: "name" } as const;
    const first = nextSort(null, column);
    expect(first).toEqual({ column, direction: "asc" });
    const second = nextSort(first, column);
    expect(second).toEqual({ column, direction: "desc" });
    expect(nextSort(second, column)).toBeNull();
  });

  it("starts a new column best first", () => {
    const sorted = nextSort({ column: { kind: "name" }, direction: "desc" }, {
      kind: "site",
      siteUuid: "s1",
    });
    expect(sorted?.direction).toBe("asc");
  });
});
