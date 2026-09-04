/**
 * Small, per-browser interface preferences: things like "show job icons"
 * that belong to how one person likes the screen, not to the project or the
 * server. Kept in localStorage so they work in the web build as well as the
 * desktop app, and broadcast within the window so every component showing
 * the same preference updates at once.
 */
import { useCallback, useEffect, useState } from "react";

const PREFIX = "ccp4i2.ui.";
const EVENT = "ccp4i2-ui-preference";

export type UiPreferenceKey = "showJobIcons";

const DEFAULTS: Record<UiPreferenceKey, boolean> = {
  showJobIcons: true,
};

export function readUiPreference(key: UiPreferenceKey): boolean {
  if (typeof window === "undefined") return DEFAULTS[key];
  try {
    const raw = window.localStorage.getItem(PREFIX + key);
    return raw === null ? DEFAULTS[key] : raw === "true";
  } catch {
    return DEFAULTS[key];
  }
}

export function writeUiPreference(key: UiPreferenceKey, value: boolean): void {
  if (typeof window === "undefined") return;
  try {
    window.localStorage.setItem(PREFIX + key, String(value));
  } catch {
    /* private mode or quota: the in-memory state still updates below */
  }
  window.dispatchEvent(new CustomEvent(EVENT, { detail: { key, value } }));
}

/** Read a preference and re-render when it changes anywhere in this window. */
export function useUiPreference(key: UiPreferenceKey): [boolean, (value: boolean) => void] {
  const [value, setValue] = useState<boolean>(() => readUiPreference(key));
  useEffect(() => {
    const onChange = (ev: Event) => {
      const detail = (ev as CustomEvent).detail;
      if (detail?.key === key) setValue(!!detail.value);
    };
    window.addEventListener(EVENT, onChange);
    return () => window.removeEventListener(EVENT, onChange);
  }, [key]);
  const set = useCallback((next: boolean) => writeUiPreference(key, next), [key]);
  return [value, set];
}
