/**
 * The chip that follows the cursor while a job is dragged out of the job list.
 *
 * It used to show a generic CCP4i2 diamond squeezed into a 20px circular
 * avatar, whatever job was being dragged and whatever the "show job icons"
 * preference said — so the icon did not match the row it came from, and the
 * circular mask clipped it top and bottom (Paul Bond, issue #432). The chip
 * now reuses the job list's own avatar, and obeys the same preference as the
 * rows do.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";

vi.mock("../theme/theme-provider", () => ({
  useTheme: () => ({ customColors: { ui: { lightBlue: "#9cf" } } }),
}));

import { JobDragChip } from "../components/job-drag-chip";
import { writeUiPreference } from "../lib/ui-preferences";
import { Job } from "../types/models";

const JOB = {
  id: 12,
  uuid: "u12",
  number: "1",
  title: "Import merged",
  task_name: "import_merged",
  status: 6,
} as unknown as Job;

beforeEach(() => {
  window.localStorage.clear();
});

describe("the dragged-job chip", () => {
  it("names the job being dragged", () => {
    render(<JobDragChip job={JOB} />);
    expect(screen.getByText(/Import merged/)).toBeTruthy();
  });

  it("shows that job's own icon, not a generic CCP4i2 one", () => {
    render(<JobDragChip job={JOB} />);
    const icon = screen.getByAltText("import_merged");
    expect(icon.getAttribute("src")).toBe("/svgicons/import_merged.svg");
  });

  it("drops the icon when job icons are turned off", () => {
    writeUiPreference("showJobIcons", false);
    render(<JobDragChip job={JOB} />);
    expect(screen.queryByAltText("import_merged")).toBeNull();
    // The chip itself stays — it is what tells you what you are dragging.
    expect(screen.getByText(/Import merged/)).toBeTruthy();
  });
});
