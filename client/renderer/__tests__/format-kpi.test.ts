// @vitest-environment node
/**
 * KPI chips show labels, never field names (#596). The server sends a label
 * for each key; a key it did not label is turned into words here, the same
 * way server/ccp4i2/lib/kpi_labels.py:humanise_key does.
 */
import { describe, expect, it } from "vitest";
import {
  formatKpiValue,
  humaniseKpiKey,
  kpiLabel,
} from "../lib/format-kpi";

describe("kpiLabel", () => {
  const labels = {
    spaceGroup: "Space group",
    highResLimit: "High resolution (Å)",
    rMeas: "Rmeas",
  };

  it("uses the server's label for a key", () => {
    expect(kpiLabel("spaceGroup", labels)).toBe("Space group");
    expect(kpiLabel("highResLimit", labels)).toBe("High resolution (Å)");
    expect(kpiLabel("rMeas", labels)).toBe("Rmeas");
  });

  it("turns an unlabelled key into words rather than showing it", () => {
    expect(kpiLabel("highResLimit")).toBe("High res limit");
    expect(kpiLabel("someNewMetric", labels)).toBe("Some new metric");
    expect(kpiLabel("RFree", null)).toBe("R free");
  });
});

describe("humaniseKpiKey agrees with the server's humanise_key", () => {
  // The same cases as server/ccp4i2/tests/unit/lib/test_kpi_labels.py
  it.each([
    ["highResLimit", "High res limit"],
    ["nEvents", "N events"],
    ["someNewMetric", "Some new metric"],
    ["mean_phase_error", "Mean phase error"],
    ["XMLFile", "XML file"],
    ["Hand1Score", "Hand1 score"],
    ["rmsd", "Rmsd"],
    ["TFZ", "TFZ"],
    ["", ""],
  ])("%s -> %s", (key, words) => {
    expect(humaniseKpiKey(key)).toBe(words);
  });
});

describe("formatKpiValue", () => {
  it("prints a count as an integer and a measurement to 3 s.f.", () => {
    expect(formatKpiValue(42)).toBe("42");
    expect(formatKpiValue(0.21345)).toBe("0.213");
  });
});
