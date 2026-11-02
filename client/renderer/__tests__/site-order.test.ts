/**
 * The order every list of a campaign's sites is shown in (lib/site-matrix,
 * compareSites): by name, alphabetically, numbers as numbers.
 *
 * NAMES/EXPECTED are shared with the server's
 * ccp4i2/tests/unit/lib/test_campaign_site_order.py, whose key must agree
 * with this collator so that the server's "first site" is the one shown
 * first here. Change one, change both.
 */
import { describe, it, expect } from "vitest";
import { compareSites, matrixSites, sortSites } from "../lib/site-matrix";
import type { CampaignSite } from "../types/campaigns";

const NAMES = [
  "Site 10", "Pocket_A_Tyr127", "site 1", "Pocket-B", "Interface_AB_top",
  "x0010", "Pocket", "Hinge_B", "Pocket_A-Met", "Site 2", "Pocket B",
  "pocket_a", "Back_A", "Interface_AB_low", "Pocket_A_Trp", "x0009",
  "Pocket_A_Tyr12", "Pocket_A",
];

const EXPECTED = [
  "Back_A", "Hinge_B", "Interface_AB_low", "Interface_AB_top", "Pocket",
  "Pocket B", "pocket_a", "Pocket_A", "Pocket_A_Trp", "Pocket_A_Tyr12",
  "Pocket_A_Tyr127", "Pocket_A-Met", "Pocket-B", "site 1", "Site 2",
  "Site 10", "x0009", "x0010",
];

const site = (name: string, id: number, extra: Partial<CampaignSite> = {}): CampaignSite =>
  ({ id, name, uuid: `u${id}`, origin: [0, 0, 0], order: 100 - id, ...extra }) as CampaignSite;

describe("sortSites", () => {
  it("is alphabetical, numbers as numbers, the same order as the server", () => {
    const sites = NAMES.map((name, i) => site(name, i + 1));
    expect(sortSites(sites).map((s) => s.name)).toEqual(EXPECTED);
  });

  it("ignores the stored order: names group a site with its subsites", () => {
    const sites = [site("Pocket_B", 1), site("Other", 2), site("Pocket_A_Trp", 3), site("Pocket_A", 4)];
    expect(sortSites(sites).map((s) => s.name)).toEqual(["Other", "Pocket_A", "Pocket_A_Trp", "Pocket_B"]);
  });

  it("does not touch its input, and copes with none", () => {
    const sites = [site("b", 1), site("a", 2)];
    sortSites(sites);
    expect(sites.map((s) => s.name)).toEqual(["b", "a"]);
    expect(sortSites(undefined)).toEqual([]);
  });

  it("breaks a tie of names equal but for case by creation", () => {
    expect(compareSites(site("POCKET", 3), site("pocket", 7))).toBeLessThan(0);
  });
});

describe("matrixSites", () => {
  it("gives the table's site columns in the same order, still dropping uuid-less sites", () => {
    const sites = [site("Zeta", 1), site("alpha", 2), site("Mid", 3, { uuid: undefined })];
    expect(matrixSites(sites).map((s) => s.name)).toEqual(["alpha", "Zeta"]);
  });
});
