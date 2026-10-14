import { describe, expect, it } from "vitest";
import { fileKind, formatModified, formatSize } from "../lib/file-kinds";

describe("fileKind", () => {
  it("opens PanDDA maps in Moorhen and offers no text preview of them", () => {
    const kind = fileKind("xtal-0001-event_1_1-BDC_0.27_map.native.ccp4");
    expect(kind.moorhen).toBe("map");
    expect(kind.preview).toBeNull();
  });

  it("treats .cif as ambiguous for Moorhen and previews it as CIF", () => {
    const kind = fileKind("LIG.cif");
    expect(kind.moorhen).toBe("cif");
    expect(kind.preview).toBe("cif");
  });

  it("gives an MTZ its header, ViewHKL, and no Moorhen", () => {
    const kind = fileKind("hklout.mtz");
    expect(kind).toMatchObject({ preview: "mtz", previewLabel: "MTZ header", viewhkl: true, moorhen: null });
  });

  it("previews images, tables and molfiles in their own views", () => {
    expect(fileKind("plot.PNG").preview).toBe("image");
    expect(fileKind("pandda_analyse_events.csv").previewLabel).toBe("Table");
    expect(fileKind("ligand.mol").preview).toBe("molblock-raw");
  });

  it("opens HTML and PDF in a new tab", () => {
    expect(fileKind("report.html").newTab).toBe("html");
    expect(fileKind("manual.pdf").newTab).toBe("pdf");
    expect(fileKind("manual.pdf").preview).toBeNull();
  });

  it("offers a text preview of unknown extensions and none of archives", () => {
    expect(fileKind("events.yaml").preview).toBe("text");
    expect(fileKind("stdout").preview).toBe("text");
    expect(fileKind("something.weird").preview).toBe("text");
    expect(fileKind("bundle.zip").preview).toBeNull();
  });
});

describe("formatting", () => {
  it("sizes in binary units", () => {
    expect(formatSize(512)).toBe("512 B");
    expect(formatSize(1536)).toBe("1.5 KB");
    expect(formatSize(3 * 1024 * 1024)).toBe("3.0 MB");
    expect(formatSize(undefined)).toBe("");
  });

  it("says how long ago a file was written", () => {
    const now = Date.UTC(2026, 9, 14, 12, 0, 0);
    const s = now / 1000;
    expect(formatModified(s - 10, now)).toBe("just now");
    expect(formatModified(s - 300, now)).toBe("5 min ago");
    expect(formatModified(s - 7200, now)).toBe("2 h ago");
    expect(formatModified(undefined, now)).toBe("");
  });
});
