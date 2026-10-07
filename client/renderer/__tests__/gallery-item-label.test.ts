// @vitest-environment jsdom
import { describe, expect, it } from "vitest";
import $ from "jquery";
import { galleryItemLabel } from "../lib/report-gallery";

// As the server serialises a Phaser report's gallery (#677): the graph's
// title is on its ccp4_data child, the graph itself has only a key.
const REPORT = `<CCP4i2ReportObjectGallery xmlns:ns0="http://www.ccp4.ac.uk/ccp4ns">
  <CCP4i2ReportFlotGraph key="PhaserGraph0" class="">
    <ns0:ccp4_data title="Figures of Merit" id="data_PhaserGraph0"/>
  </CCP4i2ReportFlotGraph>
  <CCP4i2ReportDiv key="Div_3" title="Ligand drawing"/>
  <CCP4i2ReportDiv key="Div_4"/>
</CCP4i2ReportObjectGallery>`;

describe("galleryItemLabel", () => {
  const items = $($.parseXML(REPORT)).find("CCP4i2ReportObjectGallery").children().toArray();

  it("lists a graph by its data's title, not its key", () => {
    expect(galleryItemLabel(items[0], 0)).toBe("Figures of Merit");
  });

  it("prefers an item's own title", () => {
    expect(galleryItemLabel(items[1], 1)).toBe("Ligand drawing");
  });

  it("falls back to a number, never the key", () => {
    expect(galleryItemLabel(items[2], 2)).toBe("Object 3");
  });
});
