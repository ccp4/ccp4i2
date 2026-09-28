/**
 * The i2run dialog displays the recipe the server rendered: the CCP4 setup to
 * source, the variables the app's own server is running with, where to run
 * from, and the command.
 *
 * It used to assemble a command itself from Electron's cwd and config, as
 * `${cwd}/ccp4i2/i2run/i2run.sh <args> --dbFile <...>/db.sqlite3` -- a script
 * that does not exist and a flag i2run does not accept -- and bailed out to ""
 * whenever cwd was missing, which is always in a packaged build. So the dialog
 * rendered blank and offered an empty clipboard. It also said nothing about the
 * environment, so a command copied out of a GUI running on a non-default
 * projects directory addressed a different database, silently.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";

import { I2RunDialog } from "../components/i2run-dialog";

const COMMAND =
  'ccp4-python -m ccp4i2.cli.i2run freerflag --project_name gamma --FRAC "0.07"';

const writeText = vi.fn(async (_text: string) => {});

beforeEach(() => {
  writeText.mockClear();
  Object.assign(navigator, { clipboard: { writeText } });
});

const copied = () => writeText.mock.calls[0][0];

describe("I2RunDialog", () => {
  it("shows the command with no Electron main process to ask", () => {
    // window.electronAPI is undefined here, as it is in the web build.
    render(<I2RunDialog open command={COMMAND} onClose={() => {}} />);
    expect(screen.getByText(/ccp4i2\.cli\.i2run/)).toBeTruthy();
    expect(screen.getByText(/--FRAC/)).toBeTruthy();
  });

  it("prefixes a cd when the command has to be run from somewhere", () => {
    render(
      <I2RunDialog
        open
        command={COMMAND}
        workingDirectory="/opt/ccp4i2/server"
        onClose={() => {}}
      />
    );
    expect(screen.getByText(/cd '\/opt\/ccp4i2\/server'/)).toBeTruthy();
  });

  it("exports the variables the server runs with, sourcing CCP4 first", async () => {
    render(
      <I2RunDialog
        open
        command={COMMAND}
        workingDirectory="/opt/ccp4i2/server"
        ccp4Setup="/opt/ccp4/bin/ccp4.setup-sh"
        environment={{ CCP4I2_PROJECTS_DIR: "/data/My Projects" }}
        platform="darwin"
        onClose={() => {}}
      />
    );
    fireEvent.click(screen.getByRole("button", { name: /copy/i }));
    await waitFor(() => expect(writeText).toHaveBeenCalledTimes(1));
    expect(copied()).toBe(
      [
        "source '/opt/ccp4/bin/ccp4.setup-sh'",
        "export CCP4I2_PROJECTS_DIR='/data/My Projects'",
        "cd '/opt/ccp4i2/server'",
        COMMAND,
      ].join("\n")
    );
  });

  it("uses PowerShell syntax on Windows, which has no setup script", async () => {
    render(
      <I2RunDialog
        open
        command={COMMAND}
        workingDirectory={"C:\\ccp4i2\\server"}
        ccp4Setup={null}
        environment={{ CCP4I2_HOME: "C:\\Users\\me\\.ccp4i2x" }}
        platform="win32"
        onClose={() => {}}
      />
    );
    fireEvent.click(screen.getByRole("button", { name: /copy/i }));
    await waitFor(() => expect(writeText).toHaveBeenCalledTimes(1));
    expect(copied()).toBe(
      [
        "# Run this from the CCP4 command prompt.",
        "$env:CCP4I2_HOME = 'C:\\Users\\me\\.ccp4i2x'",
        "cd 'C:\\ccp4i2\\server'",
        COMMAND,
      ].join("\n")
    );
  });

  it("quotes a value containing a quote rather than breaking out of it", async () => {
    render(
      <I2RunDialog
        open
        command={COMMAND}
        environment={{ CCP4I2_PROJECTS_DIR: "/data/it's here" }}
        platform="linux"
        onClose={() => {}}
      />
    );
    fireEvent.click(screen.getByRole("button", { name: /copy/i }));
    await waitFor(() => expect(writeText).toHaveBeenCalledTimes(1));
    expect(copied()).toBe(
      `export CCP4I2_PROJECTS_DIR='/data/it'\\''s here'\n${COMMAND}`
    );
  });

  it("says nothing about the environment when there is nothing to say", () => {
    render(
      <I2RunDialog
        open
        command={COMMAND}
        environment={{}}
        platform="darwin"
        onClose={() => {}}
      />
    );
    expect(screen.getByText(COMMAND).textContent).toBe(COMMAND);
  });

  it("disables Copy rather than offering an empty clipboard", () => {
    render(<I2RunDialog open command="" onClose={() => {}} />);
    const copy = screen.getByRole("button", { name: /copy/i });
    expect(copy.hasAttribute("disabled")).toBe(true);
  });
});
