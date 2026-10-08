// @vitest-environment node
/**
 * The `event` URL parameter a filled site-matrix box carries, and what the
 * campaign viewer loads from a receipt's answer for that event.
 */

import { describe, expect, it } from "vitest";
import {
  EventEvidence,
  eventParam,
  evidenceLoadPlan,
  hexToRgb01,
  parseEventParam,
} from "../lib/event-evidence";

const file = (id: number, name: string) => ({
  id,
  uuid: `u${id}`,
  name,
  type: null,
  annotation: "",
});

const evidence = (overrides: Partial<EventEvidence> = {}): EventEvidence => ({
  receipt_job_id: 424,
  dtag: "xtal-0001",
  event_idx: 1,
  position: 0,
  site_idx: 1,
  centroid: [1, 2, 3],
  ligand_id: "5KX",
  display_contour: 0.4328,
  optimal_contour: 1.34,
  contour: 0.4328,
  event_map: file(1306, "event_1_map.map"),
  pose: file(1307, "event_1_pose.pdb"),
  dictionary: file(1309, "DICT.cif"),
  has_map: true,
  has_pose: true,
  colour: "#3f51b5",
  radius: 12,
  ...overrides,
});

describe("event parameter", () => {
  it("round-trips a receipt and event number", () => {
    expect(eventParam(424, 1)).toBe("424:1");
    expect(parseEventParam("424:1")).toEqual({ receiptJobId: 424, eventIdx: 1 });
  });

  it("is absent without both halves", () => {
    expect(eventParam(null, 1)).toBeNull();
    expect(eventParam(424, undefined)).toBeNull();
  });

  it("rejects anything malformed", () => {
    for (const bad of [null, "", "424", "424:", ":1", "a:1", "424:1:2", "4.2:1"]) {
      expect(parseEventParam(bad)).toBeNull();
    }
  });
});

describe("hexToRgb01", () => {
  it("reads the receipt's event map colour", () => {
    expect(hexToRgb01("#3f51b5")).toEqual({ r: 63 / 255, g: 81 / 255, b: 181 / 255 });
    expect(hexToRgb01("blue")).toBeNull();
  });
});

describe("evidenceLoadPlan", () => {
  it("loads the map at the receipt's contour and the pose with its dictionary", () => {
    expect(evidenceLoadPlan(evidence())).toEqual({
      map: {
        fileId: 1306,
        name: "Event 1 map",
        contourLevel: 0.4328,
        colour: "#3f51b5",
        radius: 12,
      },
      pose: { fileId: 1307, name: "Event 1 autobuild", dictionaryFileId: 1309 },
    });
  });

  it("loads what exists and nothing else", () => {
    const plan = evidenceLoadPlan(
      evidence({ pose: null, has_pose: false, dictionary: null, contour: null })
    );
    expect(plan.pose).toBeNull();
    expect(plan.map?.contourLevel).toBeNull();
    expect(evidenceLoadPlan(evidence({ event_map: null, has_map: false })).map).toBeNull();
  });

  it("loads a pose without a dictionary when the receipt has none", () => {
    expect(evidenceLoadPlan(evidence({ dictionary: null })).pose?.dictionaryFileId).toBeNull();
  });
});
