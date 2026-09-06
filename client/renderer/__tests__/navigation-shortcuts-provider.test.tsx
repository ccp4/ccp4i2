/**
 * NavigationShortcutsProvider is mounted once in app/ccp4i2/layout.tsx so
 * back/forward keyboard shortcuts work on every route. It must translate the
 * platform chords into history.back()/forward(), ignore unmodified keys, and
 * register exactly one keydown listener even if a second copy is nested
 * inside the first (otherwise one keystroke would navigate twice).
 */
import React from "react";
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { NavigationShortcutsProvider } from "../providers/navigation-shortcuts-provider";

describe("NavigationShortcutsProvider", () => {
  let back: ReturnType<typeof vi.spyOn>;
  let forward: ReturnType<typeof vi.spyOn>;

  beforeEach(() => {
    back = vi.spyOn(window.history, "back").mockImplementation(() => {});
    forward = vi.spyOn(window.history, "forward").mockImplementation(() => {});
  });

  afterEach(() => {
    back.mockRestore();
    forward.mockRestore();
  });

  it("renders its children", () => {
    const { getByText } = render(
      <NavigationShortcutsProvider>
        <span>child</span>
      </NavigationShortcutsProvider>
    );
    expect(getByText("child")).toBeInTheDocument();
  });

  it("Alt+Left / Alt+Right go back / forward (Windows, Linux)", () => {
    render(<NavigationShortcutsProvider />);
    fireEvent.keyDown(window, { key: "ArrowLeft", altKey: true });
    expect(back).toHaveBeenCalledTimes(1);
    expect(forward).not.toHaveBeenCalled();
    fireEvent.keyDown(window, { key: "ArrowRight", altKey: true });
    expect(forward).toHaveBeenCalledTimes(1);
  });

  it("Cmd+[ / Cmd+] go back / forward (Mac)", () => {
    render(<NavigationShortcutsProvider />);
    fireEvent.keyDown(window, { key: "[", metaKey: true });
    expect(back).toHaveBeenCalledTimes(1);
    fireEvent.keyDown(window, { key: "]", metaKey: true });
    expect(forward).toHaveBeenCalledTimes(1);
  });

  it("prevents the browser default so Electron does not also navigate", () => {
    render(<NavigationShortcutsProvider />);
    const event = new KeyboardEvent("keydown", {
      key: "ArrowLeft",
      altKey: true,
      cancelable: true,
    });
    window.dispatchEvent(event);
    expect(event.defaultPrevented).toBe(true);
  });

  it("ignores the same keys without a modifier, and other modifier chords", () => {
    render(<NavigationShortcutsProvider />);
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "[" });
    fireEvent.keyDown(window, { key: "ArrowLeft", ctrlKey: true });
    fireEvent.keyDown(window, { key: "[", altKey: true });
    expect(back).not.toHaveBeenCalled();
    expect(forward).not.toHaveBeenCalled();
  });

  it("navigates once per keystroke even when nested", () => {
    render(
      <NavigationShortcutsProvider>
        <NavigationShortcutsProvider>
          <span>inner</span>
        </NavigationShortcutsProvider>
      </NavigationShortcutsProvider>
    );
    fireEvent.keyDown(window, { key: "ArrowLeft", altKey: true });
    expect(back).toHaveBeenCalledTimes(1);
  });

  it("removes its listener on unmount", () => {
    const { unmount } = render(<NavigationShortcutsProvider />);
    unmount();
    fireEvent.keyDown(window, { key: "ArrowLeft", altKey: true });
    expect(back).not.toHaveBeenCalled();
  });
});
