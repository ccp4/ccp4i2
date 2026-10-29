// @vitest-environment jsdom
/**
 * A render error stays inside its panel (#614): the boundary shows what went
 * wrong with a way to try again, and clears itself when the panel is pointed
 * at something else.
 */
import React, { useState } from "react";
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";

import {
  PanelErrorBoundary,
  describePanelError,
} from "../components/panel-error-boundary";

function Explodes({ when }: { when: boolean }) {
  if (when) throw new Error("boom in the panel");
  return <div>panel content</div>;
}

describe("PanelErrorBoundary", () => {
  beforeEach(() => vi.spyOn(console, "error").mockImplementation(() => {}));
  afterEach(() => {
    cleanup();
    vi.restoreAllMocks();
  });

  it("confines a render error to the panel and names it", () => {
    render(
      <div>
        <div>elsewhere on the page</div>
        <PanelErrorBoundary name="job list">
          <Explodes when />
        </PanelErrorBoundary>
      </div>,
    );
    expect(screen.getByText("elsewhere on the page")).toBeTruthy();
    expect(
      screen.getByText("Something went wrong showing the job list"),
    ).toBeTruthy();
    expect(screen.getByText("boom in the panel")).toBeTruthy();
  });

  it("redraws the panel on Try again", () => {
    // Read at render time, so the same element can fail once and then draw:
    // Try again re-renders the children as they were given.
    const fault = { on: true };
    function Flaky() {
      if (fault.on) throw new Error("transient");
      return <div>panel content</div>;
    }
    render(
      <PanelErrorBoundary name="job list">
        <Flaky />
      </PanelErrorBoundary>,
    );
    expect(screen.queryByText("panel content")).toBeNull();
    fault.on = false;
    fireEvent.click(screen.getByRole("button", { name: "Try again" }));
    expect(screen.getByText("panel content")).toBeTruthy();
  });

  it("clears the error when a reset key changes", () => {
    function Harness() {
      const [jobId, setJobId] = useState(1);
      return (
        <div>
          <button onClick={() => setJobId(2)}>next job</button>
          <PanelErrorBoundary name="job's current tab" resetKeys={[jobId]}>
            <Explodes when={jobId === 1} />
          </PanelErrorBoundary>
        </div>
      );
    }
    render(<Harness />);
    expect(screen.queryByText("panel content")).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "next job" }));
    expect(screen.getByText("panel content")).toBeTruthy();
  });

  it("formats a report with the panel, the page and both stacks", () => {
    const error = new Error("boom");
    error.stack = "Error: boom\n    at render";
    const text = describePanelError(
      "job list",
      error,
      "\n    at Explodes",
      "Agent/1",
      "http://x/job/3",
    );
    expect(text).toContain("CCP4i2 panel error: job list");
    expect(text).toContain("Page: http://x/job/3");
    expect(text).toContain("Error: boom");
    expect(text).toContain("Component stack:\nat Explodes");
  });
});
