/**
 * Zoom keys, read per window (main/ccp4i2-zoom.ts).
 *
 * Zoom was bound with globalShortcut "CommandOrControl+Plus/-/0": refused
 * under Wayland, so Ctrl+/Ctrl- did nothing on Ubuntu 24.04 (alpha tracker,
 * #790); "Plus" needed Shift, so Ctrl+= never zoomed; and the app took those
 * keys from every other application while it ran.
 */
import { describe, expect, it, vi } from "vitest";

vi.mock("electron", () => ({ BrowserWindow: {}, Menu: {} }));

import { zoomStepForKey } from "../../main/ccp4i2-zoom";

const key = (k: string, mods: Partial<{ control: boolean; meta: boolean; alt: boolean }> = {}) => ({
  type: "keyDown" as const,
  key: k,
  control: false,
  meta: false,
  alt: false,
  ...mods,
});

describe("zoomStepForKey", () => {
  it("zooms in on Ctrl with +, =, or numpad +, Shift or not", () => {
    for (const k of ["+", "="]) expect(zoomStepForKey(key(k, { control: true }), false)).toBe("in");
  });

  it("zooms out on Ctrl with - (or _ when Shift is held), and resets on 0", () => {
    expect(zoomStepForKey(key("-", { control: true }), false)).toBe("out");
    expect(zoomStepForKey(key("_", { control: true }), false)).toBe("out");
    expect(zoomStepForKey(key("0", { control: true }), false)).toBe("reset");
  });

  it("uses Cmd, not Ctrl, on macOS", () => {
    expect(zoomStepForKey(key("=", { meta: true }), true)).toBe("in");
    expect(zoomStepForKey(key("=", { control: true }), true)).toBeNull();
  });

  it("leaves other keys, Alt chords, key-ups and bare keys alone", () => {
    expect(zoomStepForKey(key("r", { control: true }), false)).toBeNull();
    expect(zoomStepForKey(key("=", { control: true, alt: true }), false)).toBeNull();
    expect(zoomStepForKey({ ...key("=", { control: true }), type: "keyUp" } as any, false)).toBeNull();
    expect(zoomStepForKey(key("="), false)).toBeNull();
  });
});
