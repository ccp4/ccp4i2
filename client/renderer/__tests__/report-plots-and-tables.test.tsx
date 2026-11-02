/**
 * Report plots and tables, as the Servalcat author found them (alpha tracker,
 * 2026-10-16):
 * - a plot's legend ignored the corner the report asked for;
 * - "Resolution (&Aring;)" was drawn literally on the axis;
 * - axis and legend text kept Chart.js's #666, unreadable in the dark theme;
 * - a grouped column heading ("RMSD" over Bond/Angle/Chiral) was read as one
 *   column, so every heading after it sat over the wrong data: "Weight" over
 *   the chiral RMSD, which looked like weights shifted by a cycle.
 */
import React from "react";
import { describe, it, expect } from "vitest";
import { render, within } from "@testing-library/react";
import { CCP4Table, parseXML } from "../components/report/ChartLib";
import { CCP4i2ReportTable } from "../components/report/CCP4i2ReportTable";

const GRAPH = `<ccp4_data title="Correlations">
  <headers separator=" ">1/d^2 CCIwork</headers>
  <data>0.0 0.95
0.1 0.9
0.2 0.8</data>
  <plot>
    <title>Correlations</title>
    <xlabel>Resolution (&amp;Aring;)</xlabel>
    <ylabel>Correlation</ylabel>
    <xscale>oneoversqrt</xscale>
    <legendposition x="0" y="0"/>
    <plotline xcol="1" ycol="2"><colour>orange</colour></plotline>
  </plot>
</ccp4_data>`;

const plotOptionsOf = async (xml: string) => {
  const parsed: any = await parseXML(xml);
  const table = new CCP4Table(parsed.CCP4Table);
  const plot = Array.isArray(table.plot) ? table.plot[0] : table.plot!;
  return table.plotOptions(plot, { text: "rgb(1, 2, 3)", grid: "rgb(4, 5, 6)" }) as any;
};

describe("report plots", () => {
  it("puts the legend where the report asked", async () => {
    const options = await plotOptionsOf(GRAPH);
    expect(options.plugins.legendPlacement).toEqual({ x: 0, y: 0 });
  });

  it("asks for no placement when the report gives none", async () => {
    const options = await plotOptionsOf(GRAPH.replace('<legendposition x="0" y="0"/>', ""));
    expect(options.plugins.legendPlacement).toBeUndefined();
  });

  it("draws the entity as the character it names", async () => {
    const options = await plotOptionsOf(GRAPH);
    expect(options.scales.x.title.text).toBe("Resolution (Å)");
  });

  it("labels infinite resolution, not 'Infinity'", async () => {
    const options = await plotOptionsOf(GRAPH);
    expect(options.scales.x.ticks.callback(0, 0)).toBe("∞");
    expect(options.scales.x.ticks.callback(0.25, 1)).toBe("2.00");
  });

  it("takes text and grid colours from the theme", async () => {
    const options = await plotOptionsOf(GRAPH);
    expect(options.plugins.title.color).toBe("rgb(1, 2, 3)");
    expect(options.plugins.legend.labels.color).toBe("rgb(1, 2, 3)");
    for (const scale of [options.scales.x, options.scales.yAxisLeft]) {
      expect(scale.ticks.color).toBe("rgb(1, 2, 3)");
      expect(scale.title.color).toBe("rgb(1, 2, 3)");
      expect(scale.grid.color).toBe("rgb(4, 5, 6)");
    }
  });
});

const tableItem = (thead: string) => {
  const wrapper = document.createElement("div");
  wrapper.setAttribute("transpose", "False");
  wrapper.innerHTML = `<table>${thead}<tbody>
    <tr><td>1</td><td>0.016</td><td>1.9</td><td>0.108</td><td>0.10</td></tr>
  </tbody></table>`;
  return wrapper;
};

const GROUPED = `<thead>
  <tr><th>Cycle</th><th colspan="3">RMSD</th><th>Weight</th></tr>
  <tr><th></th><th>Bond</th><th>Angle</th><th>Chiral</th><th></th></tr>
</thead>`;

describe("report tables with grouped headings", () => {
  it("spans the group and keeps every data cell", () => {
    const { container } = render(
      <CCP4i2ReportTable iItem={0} item={tableItem(GROUPED)} job={{} as any} />
    );
    const [groups, subtitles] = Array.from(container.querySelectorAll("thead tr"));
    const rmsd = within(groups as HTMLElement).getByText("RMSD").closest("th")!;
    expect(rmsd.getAttribute("colspan")).toBe("3");
    expect(within(subtitles as HTMLElement).getByText("Chiral")).toBeTruthy();
    const cells = Array.from(container.querySelectorAll("tbody tr td")).map((td) => td.textContent);
    expect(cells).toEqual(["1", "0.016", "1.9", "0.108", "0.10"]);
  });

  it("does not drop data cells when only the grouped row is sent", () => {
    const { container } = render(
      <CCP4i2ReportTable
        iItem={0}
        item={tableItem(`<thead><tr><th>Cycle</th><th colspan="3">RMSD</th><th>Weight</th></tr></thead>`)}
        job={{} as any}
      />
    );
    const cells = Array.from(container.querySelectorAll("tbody tr td")).map((td) => td.textContent);
    expect(cells).toEqual(["1", "0.016", "1.9", "0.108", "0.10"]);
  });
});
