/**
 * A container draws a child unless a qualifier says not to: `hidden` (from a
 * PHIL .style such as Phaser's keywords.general.mute) or an expertLevel above
 * the chosen one -- and a container with nothing drawable beneath it is not
 * drawn either.
 */
import { describe, it, expect } from "vitest";

import {
  shouldRenderChild,
  passesExpertLevel,
} from "../components/task/task-elements/container-filters";

const leaf = (qualifiers: Record<string, any> = {}) => ({
  _baseClass: "CFloat",
  _qualifiers: qualifiers,
});
const container = (children: Record<string, any>, qualifiers = {}) => ({
  _baseClass: "CContainer",
  _qualifiers: qualifiers,
  _value: children,
});

describe("shouldRenderChild", () => {
  it("draws an ordinary leaf", () => {
    expect(shouldRenderChild(leaf(), undefined)).toBe(true);
    expect(shouldRenderChild(leaf(), 0)).toBe(true);
  });

  it("never draws a hidden item, whatever the expert level", () => {
    expect(shouldRenderChild(leaf({ hidden: true }), undefined)).toBe(false);
    expect(shouldRenderChild(leaf({ hidden: true }), 9)).toBe(false);
  });

  it("filters by expert level only when one is chosen", () => {
    expect(shouldRenderChild(leaf({ expertLevel: 2 }), undefined)).toBe(true);
    expect(shouldRenderChild(leaf({ expertLevel: 2 }), 1)).toBe(false);
    expect(shouldRenderChild(leaf({ expertLevel: 2 }), 2)).toBe(true);
  });

  it("does not draw a container whose only children are hidden or expert", () => {
    const c = container({ a: leaf({ hidden: true }), b: leaf({ expertLevel: 3 }) });
    expect(passesExpertLevel(c, 0)).toBe(false);
    expect(shouldRenderChild(container({ a: leaf({ hidden: true }) }), 0)).toBe(false);
    expect(shouldRenderChild(container({ a: leaf() }), 0)).toBe(true);
  });

  it("reads a def.xml expert level from guiDefinition", () => {
    expect(shouldRenderChild(leaf({ guiDefinition: { expertLevel: 1 } }), 0)).toBe(false);
  });
});
