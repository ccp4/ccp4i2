/**
 * An item whose _class the client does not know must still render as its
 * _baseClass: the server generates a CContainer subclass per repeatable
 * PHIL scope, and before this fallback such an item rendered as its class
 * name in plain text.
 */
import { describe, it, expect } from "vitest";

import { resolveComponentEntry } from "../components/task/task-elements/component-dispatch";

const registry = {
  CContainer: { component: "container" },
  CFloat: { component: "float" },
  CEnsemble: { component: "ensemble" },
};

describe("resolveComponentEntry", () => {
  it("dispatches a known class directly", () => {
    expect(resolveComponentEntry(registry, { _class: "CEnsemble", _baseClass: "CData" }))
      .toBe(registry.CEnsemble);
  });

  it("falls back to the base class for an unknown subclass", () => {
    expect(resolveComponentEntry(registry, {
      _class: "phaser__composition__chain",
      _baseClass: "CContainer",
    })).toBe(registry.CContainer);
  });

  it("yields nothing for an item with no usable class", () => {
    expect(resolveComponentEntry(registry, { _class: "Mystery", _baseClass: "Mystery" }))
      .toBeUndefined();
    expect(resolveComponentEntry(registry, undefined)).toBeUndefined();
  });
});
