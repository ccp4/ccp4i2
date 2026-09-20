/**
 * The Sites cell, rendered.
 *
 * The rules live in lib/site-verdicts and are tested there; these cover what
 * only rendering shows — that an empty verdict and an unevaluated site are
 * equally invisible, that the cell disappears entirely rather than leaving a
 * stray fragment, and that a chip is a link only when there is a job to open.
 */

import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import { SiteVerdictChips } from "../components/campaigns/site-verdict-chips";
import {
  CampaignJobInfo,
  MemberProjectWithSummary,
} from "../types/campaigns";

const job = (id: number, number: string, status = 6) =>
  ({ id, number, status, task_name: "servalcat_pipe" } as CampaignJobInfo);

const project = (
  overrides: Partial<MemberProjectWithSummary> = {}
): MemberProjectWithSummary =>
  ({
    id: 1,
    name: "ds_a",
    jobs: [job(42, "1")],
    kpis: {},
    ...overrides,
  } as MemberProjectWithSummary);

describe("SiteVerdictChips", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("chips a hit and an unclear", () => {
    render(
      <SiteVerdictChips
        campaignId={7}
        project={project({
          site_evaluations: [
            { site_id: 1, site_name: "Pocket A", verdict: "hit" },
            { site_id: 2, site_name: "Pocket B", verdict: "unclear" },
          ],
          sites_evaluated: 2,
          sites_total: 2,
        })}
      />
    );

    expect(screen.getByText("Pocket A")).toBeTruthy();
    expect(screen.getByText("Pocket B")).toBeTruthy();
  });

  it("shows nothing for a dataset examined and found empty throughout", () => {
    // Empties are not sent, so this looks the same as unevaluated here -- the
    // count is what tells them apart, and it is silent once complete.
    const { container } = render(
      <SiteVerdictChips
        campaignId={7}
        project={project({
          site_evaluations: [],
          sites_evaluated: 40,
          sites_total: 40,
        })}
      />
    );

    expect(container.firstChild).toBeNull();
  });

  it("shows the count while a dataset is part-evaluated", () => {
    render(
      <SiteVerdictChips
        campaignId={7}
        project={project({
          site_evaluations: [],
          sites_evaluated: 12,
          sites_total: 40,
        })}
      />
    );

    expect(screen.getByText("12/40")).toBeTruthy();
  });

  it("collapses the tail into a +N", () => {
    render(
      <SiteVerdictChips
        campaignId={7}
        maxChips={2}
        project={project({
          site_evaluations: [
            { site_id: 1, site_name: "A", verdict: "hit" },
            { site_id: 2, site_name: "B", verdict: "hit" },
            { site_id: 3, site_name: "C", verdict: "hit" },
            { site_id: 4, site_name: "D", verdict: "hit" },
          ],
          sites_evaluated: 4,
          sites_total: 4,
        })}
      />
    );

    expect(screen.getByText("+2")).toBeTruthy();
  });

  it("opens the campaign view on that job and site", () => {
    const open = vi.spyOn(window, "open").mockImplementation(() => null as never);

    render(
      <SiteVerdictChips
        campaignId={7}
        project={project({
          site_evaluations: [
            { site_id: 3, site_name: "Pocket A", verdict: "hit" },
          ],
          sites_evaluated: 1,
          sites_total: 4,
        })}
      />
    );
    fireEvent.click(screen.getByText("Pocket A"));

    expect(open).toHaveBeenCalledWith(
      "/ccp4i2/moorhen-page/campaign/7?job=42&site=3",
      "_blank"
    );
  });

  it("does not navigate when the dataset has no job", () => {
    const open = vi.spyOn(window, "open").mockImplementation(() => null as never);

    render(
      <SiteVerdictChips
        campaignId={7}
        project={project({
          jobs: [],
          site_evaluations: [
            { site_id: 3, site_name: "Pocket A", verdict: "hit" },
          ],
          sites_evaluated: 1,
          sites_total: 4,
        })}
      />
    );
    fireEvent.click(screen.getByText("Pocket A"));

    expect(open).not.toHaveBeenCalled();
  });
});
