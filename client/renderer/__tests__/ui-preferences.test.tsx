/**
 * Per-browser interface preferences (View menu: show job icons) survive a
 * reload and reach every component showing them at once.
 */
import { describe, expect, it, beforeEach } from "vitest";
import { act, renderHook } from "@testing-library/react";
import {
  readUiPreference,
  useUiPreference,
  writeUiPreference,
} from "../lib/ui-preferences";

beforeEach(() => {
  window.localStorage.clear();
});

describe("ui preferences", () => {
  it("defaults to showing job icons", () => {
    expect(readUiPreference("showJobIcons")).toBe(true);
  });

  it("remembers a choice across reads", () => {
    writeUiPreference("showJobIcons", false);
    expect(readUiPreference("showJobIcons")).toBe(false);
    expect(window.localStorage.getItem("ccp4i2.ui.showJobIcons")).toBe("false");
  });

  it("every hook instance in the window follows a change from any of them", () => {
    const menu = renderHook(() => useUiPreference("showJobIcons"));
    const jobList = renderHook(() => useUiPreference("showJobIcons"));
    expect(jobList.result.current[0]).toBe(true);

    act(() => menu.result.current[1](false));

    expect(menu.result.current[0]).toBe(false);
    expect(jobList.result.current[0]).toBe(false);
  });

  it("a new hook starts from what was stored", () => {
    writeUiPreference("showJobIcons", false);
    const { result } = renderHook(() => useUiPreference("showJobIcons"));
    expect(result.current[0]).toBe(false);
  });
});
