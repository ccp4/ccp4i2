/**
 * The About dialog's "This session" section.
 *
 * Both ports are picked afresh at launch (detect-port from 3000 up) and the
 * local-session token dies with the process, so a user who wants to reach
 * the backend from a terminal has no way to discover either. The dialog is
 * the one place that knows both, and it must hand them over verbatim --
 * including a curl line whose port and bearer token match this launch.
 */
import { describe, expect, it, vi, beforeEach, afterEach } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import HelpMenu from "@/components/help-menu";

const TOKEN = "a".repeat(64);

/** Stand in for the preload's contextBridge surface. */
function installLocalSession(token?: string) {
  if (token === undefined) {
    delete (window as any).ccp4i2LocalSession;
    return;
  }
  (window as any).ccp4i2LocalSession = {
    token,
    userEmail: "martin@ccp4i2.invalid",
  };
}

/**
 * Stand in for the Electron main process: "get-config" is answered with the
 * ports it chose, delivered on the "message-from-main" channel the renderer
 * already listens to elsewhere.
 */
function installElectronApi(config: Record<string, unknown> | null) {
  const listeners: Record<string, Function[]> = {};
  (window as any).electronAPI = {
    sendMessage: vi.fn((channel: string) => {
      if (channel !== "get-config" || !config) return;
      for (const fn of listeners["message-from-main"] ?? []) {
        fn({}, { message: "get-config", config });
      }
    }),
    onMessage: (channel: string, fn: Function) => {
      (listeners[channel] ??= []).push(fn);
    },
    removeMessageListener: (channel: string, fn: Function) => {
      listeners[channel] = (listeners[channel] ?? []).filter((f) => f !== fn);
    },
    sendSync: vi.fn(),
    invoke: vi.fn(),
  };
}

const openAbout = () => {
  fireEvent.click(screen.getByRole("button", { name: "Help" }));
  fireEvent.click(screen.getByText("About CCP4i2"));
};

/** The value of the read-only field carrying the given label. */
const fieldValue = (label: RegExp | string) =>
  (screen.getByLabelText(label) as HTMLInputElement).value;

describe("About dialog: this session", () => {
  beforeEach(() => {
    // Build info is best-effort; a rejecting fetch exercises that path too.
    vi.stubGlobal("fetch", vi.fn(() => Promise.reject(new Error("no build info"))));
  });

  afterEach(() => {
    delete (window as any).electronAPI;
    delete (window as any).ccp4i2LocalSession;
    vi.unstubAllGlobals();
  });

  it("reports the ports the main process chose, not the defaults", async () => {
    installLocalSession(TOKEN);
    installElectronApi({ UVICORN_PORT: 3101, NEXT_PORT: 3100 });
    render(<HelpMenu />);
    openAbout();

    await waitFor(() => expect(fieldValue(/Backend/)).toBe("3101"));
    expect(fieldValue(/Front-end/)).toBe("3100");
  });

  it("offers the token and a curl line carrying that exact token and port", async () => {
    installLocalSession(TOKEN);
    installElectronApi({ UVICORN_PORT: 3101, NEXT_PORT: 3100 });
    render(<HelpMenu />);
    openAbout();

    await waitFor(() => expect(fieldValue(/Session token/)).toBe(TOKEN));
    const command = fieldValue(/Example request/);
    expect(command).toContain(`Authorization: Bearer ${TOKEN}`);
    expect(command).toContain("http://localhost:3101/api/ccp4i2/projects/");
  });

  it("copies a value to the clipboard when its button is pressed", async () => {
    const writeText = vi.fn(() => Promise.resolve());
    Object.assign(navigator, { clipboard: { writeText } });
    installLocalSession(TOKEN);
    installElectronApi({ UVICORN_PORT: 3101, NEXT_PORT: 3100 });
    render(<HelpMenu />);
    openAbout();

    await waitFor(() => expect(fieldValue(/Session token/)).toBe(TOKEN));
    fireEvent.click(screen.getByRole("button", { name: /Copy session token/i }));
    expect(writeText).toHaveBeenCalledWith(TOKEN);
    // The button flips to a tick once the clipboard promise settles; wait for
    // it so the state update happens inside the test, not after it.
    await screen.findByRole("button", { name: /Copied/i });
  });

  it("says nothing about a session in the web build, where there is none", async () => {
    installLocalSession(undefined);
    render(<HelpMenu />);
    openAbout();

    await waitFor(() => expect(screen.getByText(/graphical environment/)).toBeTruthy());
    expect(screen.queryByText("This session")).toBeNull();
    expect(screen.queryByLabelText(/Session token/)).toBeNull();
  });
});
