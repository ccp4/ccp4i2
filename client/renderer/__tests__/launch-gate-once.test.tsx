/**
 * #732: "Starting CCP4i2..." showed on every visit to the projects page,
 * because the gate wrapping it probed /health from scratch on each mount.
 * Once the backend has answered in this session the gate stays open.
 */
import React from "react";
import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { LaunchGate } from "../components/launch-gate";

afterEach(() => {
  vi.unstubAllGlobals();
});

describe("LaunchGate", () => {
  it("shows the starting state only until the server first answers", async () => {
    vi.stubGlobal("fetch", vi.fn(async () => ({ ok: true })));
    const first = render(
      <LaunchGate>
        <div>projects</div>
      </LaunchGate>
    );
    expect(screen.getByText(/Starting CCP4i2/)).toBeTruthy();
    expect(await screen.findByText("projects")).toBeTruthy();
    first.unmount();

    // Coming back to the page: no probe, no starting screen, even if the
    // health endpoint would now be slow to answer.
    const fetchMock = vi.fn(() => new Promise(() => {}));
    vi.stubGlobal("fetch", fetchMock);
    render(
      <LaunchGate>
        <div>projects</div>
      </LaunchGate>
    );
    expect(screen.queryByText(/Starting CCP4i2/)).toBeNull();
    expect(screen.getByText("projects")).toBeTruthy();
    expect(fetchMock).not.toHaveBeenCalled();
  });
});
