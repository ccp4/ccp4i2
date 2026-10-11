// @vitest-environment jsdom
import React from "react";
import { describe, it, expect, vi } from "vitest";

// The gallery draws its chosen object through the recursive element, which
// imports every report component; a stub keeps the test on the gallery.
vi.mock("../components/report/CCP4i2ReportElement", () => ({
  CCP4i2ReportElement: (props: any) => <div data-testid="object">{props.item?.textContent}</div>,
  cssToDict: (css: string) => Object.fromEntries(
    css.split(";").filter(Boolean).map((kv) => kv.split(":").map((x) => x.trim()))),
}));
import { render, screen } from "@testing-library/react";
import { CCP4i2ReportText } from "../components/report/CCP4i2ReportText";
import { CCP4i2ReportProgressBar } from "../components/report/CCP4i2ReportProgressBar";
import { CCP4i2ReportObjectGallery } from "../components/report/CCP4i2ReportObjectGallery";

/**
 * The running Phaser report looked crushed: every addText was a bare inline
 * span that dropped its style, the progress bar read "0/0 (0.0%)" in red
 * between stages, and the gallery gave half its width to a list of titles.
 */
const element = (xml: string) =>
  new DOMParser().parseFromString(xml, "text/xml").documentElement as unknown as HTMLElement;

const job: any = { id: 1, status: 3 };

describe("report text", () => {
  it("is a block that keeps the style the report asked for", () => {
    const { container } = render(
      <CCP4i2ReportText iItem={0} job={job}
        item={element('<CCP4i2ReportText style="font-weight:bold;color:green;">Now: ROTATION FUNCTION</CCP4i2ReportText>')} />
    );
    const node = container.firstElementChild as HTMLElement;
    expect(node.tagName).toBe("DIV");
    expect(node.style.fontWeight).toBe("bold");
    expect(node.style.color).toBe("green");
    expect(node.textContent).toBe("Now: ROTATION FUNCTION");
  });

  it("renders nothing for an empty line", () => {
    const { container } = render(
      <CCP4i2ReportText iItem={0} job={job} item={element("<CCP4i2ReportText> </CCP4i2ReportText>")} />
    );
    expect(container.firstElementChild).toBeNull();
  });
});

describe("report progress bar", () => {
  it("is busy, with no counts, when there is no maximum yet", () => {
    render(
      <CCP4i2ReportProgressBar iItem={0} job={job}
        item={element('<CCP4i2ReportProgress label="Generating Statistics" value="0" max="0"/>')} />
    );
    expect(screen.getByText("Generating Statistics")).toBeTruthy();
    expect(screen.queryByText(/0\/0/)).toBeNull();
    expect(document.querySelector(".MuiLinearProgress-indeterminate")).not.toBeNull();
  });

  it("shows the fraction when it has one", () => {
    render(
      <CCP4i2ReportProgressBar iItem={0} job={job}
        item={element('<CCP4i2ReportProgress label="Cycle" value="3" max="10"/>')} />
    );
    expect(screen.getByText("3/10 (30.0%)")).toBeTruthy();
    expect(document.querySelector(".MuiLinearProgress-determinate")).not.toBeNull();
  });
});

describe("report object gallery", () => {
  it("lists the titles in a narrow selectable column with the first chosen", () => {
    const gallery = element(
      '<CCP4i2ReportObjectGallery>' +
      '<CCP4i2ReportText title="Cell Content Analysis">one</CCP4i2ReportText>' +
      '<CCP4i2ReportText title="Intensity distribution for Data">two</CCP4i2ReportText>' +
      '</CCP4i2ReportObjectGallery>');
    render(<CCP4i2ReportObjectGallery iItem={0} job={job} item={gallery} />);
    const buttons = document.querySelectorAll(".MuiListItemButton-root");
    expect(buttons.length).toBe(2);
    expect(buttons[0].classList.contains("Mui-selected")).toBe(true);
    expect(screen.getByText("Cell Content Analysis")).toBeTruthy();
    expect(document.querySelector(".MuiGrid2-grid-md-3")).not.toBeNull();
  });
});
