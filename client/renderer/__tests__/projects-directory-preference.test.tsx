/**
 * The Preferences panel for the default projects directory.
 *
 * Its two buttons write through the config API and then re-read it, so the
 * panel shows what the server actually did rather than what was asked for.
 * That matters because the thing being changed is resolved on the server: a
 * reset has no value to echo back, and the server can refuse outright.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";

const apiGet = vi.fn();
const apiPatch = vi.fn();
vi.mock("../api-fetch", () => ({
  apiGet: (...args: any[]) => apiGet(...args),
  apiPatch: (...args: any[]) => apiPatch(...args),
}));

import { ProjectsDirectory } from "@/components/projects-directory";

const DEFAULT_DIR = "/home/someone/.ccp4i2x/projects";
const CHOSEN_DIR = "/data/crystallography";

/** The GET body, with whatever the server currently reports as the default. */
const state = (directory: string, editable = true) => ({
  data: { directory, default: DEFAULT_DIR, editable },
});

function mockElectronPicker(picked: string | null) {
  (window as any).electronAPI = { invoke: vi.fn(async () => picked) };
}

beforeEach(() => {
  apiGet.mockReset();
  apiPatch.mockReset();
  delete (window as any).electronAPI;
});

describe("ProjectsDirectory", () => {
  it("shows the current default and says so when it is the built-in one", async () => {
    apiGet.mockResolvedValue(state(DEFAULT_DIR));
    render(<ProjectsDirectory />);

    await screen.findByDisplayValue(DEFAULT_DIR);
    expect(screen.getByText("This is the default.")).toBeTruthy();
    // Nothing to reset to: the button is there but inert.
    expect(screen.getByRole("button", { name: "Reset" }).hasAttribute("disabled")).toBe(
      true
    );
  });

  it("writes a browsed directory and shows what came back", async () => {
    apiGet
      .mockResolvedValueOnce(state(DEFAULT_DIR))
      .mockResolvedValue(state(CHOSEN_DIR));
    apiPatch.mockResolvedValue({ data: { directory: CHOSEN_DIR } });
    mockElectronPicker(CHOSEN_DIR);

    render(<ProjectsDirectory />);
    await screen.findByDisplayValue(DEFAULT_DIR);
    fireEvent.click(screen.getByRole("button", { name: "Change" }));

    await waitFor(() =>
      expect(apiPatch).toHaveBeenCalledWith("config/default-project-parent/set/", {
        directory: CHOSEN_DIR,
      })
    );
    await screen.findByDisplayValue(CHOSEN_DIR);
    expect(screen.getByText(`Reset restores ${DEFAULT_DIR}`)).toBeTruthy();
  });

  it("resets by asking for no directory at all", async () => {
    apiGet
      .mockResolvedValueOnce(state(CHOSEN_DIR))
      .mockResolvedValue(state(DEFAULT_DIR));
    apiPatch.mockResolvedValue({ data: { directory: DEFAULT_DIR } });

    render(<ProjectsDirectory />);
    await screen.findByDisplayValue(CHOSEN_DIR);
    fireEvent.click(screen.getByRole("button", { name: "Reset" }));

    await waitFor(() =>
      expect(apiPatch).toHaveBeenCalledWith("config/default-project-parent/set/", {
        directory: null,
      })
    );
    await screen.findByDisplayValue(DEFAULT_DIR);
  });

  it("accepts a typed path on blur", async () => {
    apiGet
      .mockResolvedValueOnce(state(DEFAULT_DIR))
      .mockResolvedValue(state(CHOSEN_DIR));
    apiPatch.mockResolvedValue({ data: { directory: CHOSEN_DIR } });

    render(<ProjectsDirectory />);
    const field = await screen.findByDisplayValue(DEFAULT_DIR);
    fireEvent.change(field, { target: { value: CHOSEN_DIR } });
    fireEvent.blur(field);

    await waitFor(() =>
      expect(apiPatch).toHaveBeenCalledWith("config/default-project-parent/set/", {
        directory: CHOSEN_DIR,
      })
    );
  });

  it("reports a refusal instead of appearing to have done nothing", async () => {
    apiGet.mockResolvedValue(state(DEFAULT_DIR));
    apiPatch.mockRejectedValue(new Error("Cannot use [/nope]: Permission denied"));
    mockElectronPicker("/nope");

    render(<ProjectsDirectory />);
    await screen.findByDisplayValue(DEFAULT_DIR);
    fireEvent.click(screen.getByRole("button", { name: "Change" }));

    await screen.findByText("Cannot use [/nope]: Permission denied");
    // The panel still shows the setting as it really is.
    expect(screen.getByDisplayValue(DEFAULT_DIR)).toBeTruthy();
  });

  it("is read-only where the server says the setting is not editable", async () => {
    apiGet.mockResolvedValue(state("/srv/projects", false));
    render(<ProjectsDirectory />);

    await screen.findByDisplayValue("/srv/projects");
    expect(screen.queryByRole("button", { name: "Change" })).toBeNull();
    expect(screen.queryByRole("button", { name: "Reset" })).toBeNull();
    expect(screen.getByText(/CCP4I2_PROJECTS_DIR/)).toBeTruthy();
  });
});
