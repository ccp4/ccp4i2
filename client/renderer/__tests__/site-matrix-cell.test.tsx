/**
 * One box of the campaign overview's site matrix, rendered. The rules are
 * tested in site-matrix.test.ts; these cover what only rendering shows: the
 * box drawn, that it is a link only with a model to open, and the dash for a
 * dataset whose model is in another frame.
 */

import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import {
  SiteMatrixCellContent,
  SiteMatrixLegend,
} from "../components/campaigns/site-matrix";
import {
  CampaignSite,
  MemberProjectWithSummary,
} from "../types/campaigns";

const site: CampaignSite = {
  id: 5,
  uuid: "u5",
  name: "Pocket A",
  origin: [0, 0, 0],
};

const model = { id: 42, uuid: "j", number: "3", task_name: "i2Dimple" };

const project = (
  overrides: Partial<MemberProjectWithSummary> = {}
): MemberProjectWithSummary =>
  ({
    id: 1,
    name: "ds",
    jobs: [],
    kpis: {},
    current_model_job: model,
    frame_mismatch: null,
    site_cells: {
      u5: {
        event: {
          event_idx: 1,
          hit_probability: 0.9,
          distance: 2,
          has_pose: false,
        },
        verdict: "hit",
      },
    },
    ...overrides,
  } as MemberProjectWithSummary);

describe("SiteMatrixCellContent", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("fills a green box for an accepted event and opens the site", () => {
    const open = vi.spyOn(window, "open").mockImplementation(() => null as never);
    render(<SiteMatrixCellContent project={project()} site={site} campaignId={7} />);

    const box = screen.getByTestId("site-box");
    expect(box.getAttribute("data-kind")).toBe("filled");
    expect(box.getAttribute("data-colour")).toBe("success");

    fireEvent.click(screen.getByRole("link"));
    expect(open).toHaveBeenCalledWith(
      "/ccp4i2/moorhen-page/campaign/7?job=42&site=5",
      "_blank"
    );
  });

  it("outlines a verdict with no event", () => {
    render(
      <SiteMatrixCellContent
        project={project({
          site_cells: { u5: { event: null, verdict: "empty" } },
        })}
        site={site}
        campaignId={7}
      />
    );
    const box = screen.getByTestId("site-box");
    expect(box.getAttribute("data-kind")).toBe("outlined");
    expect(box.getAttribute("data-colour")).toBe("error");
  });

  it("is not a link without a model job", () => {
    render(
      <SiteMatrixCellContent
        project={project({ current_model_job: null })}
        site={site}
        campaignId={7}
      />
    );
    expect(screen.getByTestId("site-box")).toBeTruthy();
    expect(screen.queryByRole("link")).toBeNull();
    expect(
      screen.getByLabelText(/No refinement or dimple run yet/)
    ).toBeTruthy();
  });

  it("draws nothing for an empty cell", () => {
    render(
      <SiteMatrixCellContent
        project={project({ site_cells: { u5: { event: null, verdict: null } } })}
        site={site}
        campaignId={7}
      />
    );
    expect(screen.queryByTestId("site-box")).toBeNull();
    expect(screen.queryByRole("link")).toBeNull();
  });

  it("shows a dash with the reason for a frame mismatch", () => {
    render(
      <SiteMatrixCellContent
        project={project({ frame_mismatch: "Different space group" })}
        site={site}
        campaignId={7}
      />
    );
    expect(screen.queryByTestId("site-box")).toBeNull();
    expect(screen.getByText("–")).toBeTruthy();
    expect(screen.getByLabelText(/Different space group/)).toBeTruthy();
  });
});

describe("SiteMatrixLegend", () => {
  it("names every verdict", () => {
    render(<SiteMatrixLegend />);
    for (const label of ["Not evaluated", "Hit", "Uncertain", "Rejected"]) {
      expect(screen.getByText(label)).toBeTruthy();
    }
  });
});
