import { assert } from "console";
import $ from "jquery";

const lineColors = ["red", "green", "blue", "cyan", "magenta", "yellow"];

export interface CCP4i2Chart {
  plots: CCP4i2Plot[];
  datas: { [key: string]: (number | string)[][] };
  headers: string[];
}

export interface CCP4i2Axis {
  plotLine: CCP4i2PlotLine;
  xAxisId?: number;
  yAxisId?: number;
  dataKey: number;
  position: "insideBottom" | "insideLeft" | "insideRight";
  orientation?: "left" | "right";
  angle?: number;
  domain: [number | "auto", number | "auto"];
  tickFormatter?: (value: number) => string;
}

interface CCP4i2Plot {
  plottype: string;
  title: string;
  xLabel: string;
  yLabel: string;
  ryLabel: string;
  plotLines: CCP4i2PlotLine[];
  barCharts: CCP4i2PlotLine[];
  yAxes: { [key: string]: CCP4i2Axis };
  xAxes: { [key: string]: CCP4i2Axis };
  leftYAxes: { [key: string]: CCP4i2Axis };
  rightYAxes: { [key: string]: CCP4i2Axis };
  xintegral: string;
  yintegral: string;
  xscale: string;
  polygons: any[];
  circles: any[];
  lines: any[];
  fixaspectratio: boolean;
  xrange: { min: number | "auto"; max: number | "auto" };
  yrange: { min: number | "auto"; max: number | "auto" };
}

interface CCP4i2PlotLine {
  type: string;
  xcol: number;
  ycol: number;
  dataid: string;
  color?: string;
  linestyle?: string;
  fillcolour?: string;
  label: string | false;
  symbolsize: number;
  xAxisId?: number;
  yAxisId?: number;
}

interface CCP4i2Circle {
  xpos: number;
  ypos: number;
  radius: number;
  linecolour: string | undefined;
  alpha: number | undefined;
}

interface CCP4i2Line {
  x1: number;
  x2: number;
  y1: number;
  y2: number;
  linestyle: string | undefined;
  linecolour: string | undefined;
}

const parseHeaders = (item: Element): string[] => {
  const headersNode = $(item).find("headers");
  let separator: string | RegExp | undefined = headersNode.attr("separator");
  if (!separator || separator === " ") separator = /\s+|\t+/;

  return headersNode.text().trim().split(separator).map(decodeHTMLEntities);
};

const parseData = (item: Element): { [key: string]: (number | string)[][] } => {
  const datas: { [key: string]: (number | string)[][] } = {};
  $(item)
    .find("data")
    .toArray()
    .forEach((data, iData) => {
      let dataId = $(data).attr("id") || `${iData}`;
      datas[dataId] = $(data)
        .text()
        .trim()
        .split(/\n|\r/)
        .map((rowString) => {
          return rowString
            .trim()
            .split(/\s+|\t+/)
            .map((item) => {
              const value = parseFloat(item);
              return Number.isNaN(value) ? item : value;
            });
        });
    });
  return datas;
};

const parsePlotLines = (
  plot: Element,
  newChart: CCP4i2Chart
): CCP4i2PlotLine[] => {
  const plotLines: CCP4i2PlotLine[] = [];
  const allTypes = $(plot).find("plotline, barchart").toArray();

  allTypes.forEach((plotLine, iPlotLine) => {
    let newPlotLine: CCP4i2PlotLine | null = null;

    if (plotLine.nodeName === "plotline") {
      newPlotLine = parsePlotLine(plotLine, newChart, iPlotLine);
    } else if (plotLine.nodeName === "barchart") {
      newPlotLine = parseBarChart(plotLine, newChart, iPlotLine);
    }

    if (newPlotLine) {
      plotLines.push(newPlotLine);
    }
  });

  return plotLines;
};

const parsePlotLine = (
  plotLine: Element,
  newChart: CCP4i2Chart,
  iPlotLine: number
): CCP4i2PlotLine => {
  const xcol = $(plotLine).attr("xcol");
  const ycol = $(plotLine).attr("ycol");

  let newPlotLine: CCP4i2PlotLine = {
    type: plotLine.nodeName,
    xcol: parseInt(xcol || "0"),
    ycol: parseInt(ycol || "0"),
    dataid: $(plotLine).attr("ycol") as string,
    color: $(plotLine).attr("colour")
      ? $(plotLine).attr("colour")
      : lineColors[iPlotLine % 6],
    linestyle: $(plotLine).find("linestyle").text(),
    label: $(plotLine).find("label").text() || false,
    symbolsize: parseInt($(plotLine).find("symbolsize").text() || "1"),
  };

  finalizePlotLine(newPlotLine, newChart, plotLine, iPlotLine);
  return newPlotLine;
};

const parseBarChart = (
  plotLine: Element,
  newChart: CCP4i2Chart,
  iPlotLine: number
): CCP4i2PlotLine => {
  const col = $(plotLine).attr("col");
  const tcol = $(plotLine).attr("tcol");

  let newPlotLine: CCP4i2PlotLine = {
    type: plotLine.nodeName,
    xcol: parseInt(col || "0"),
    ycol: parseInt(tcol || "0"),
    dataid:
      $(plotLine).attr("tcol") || ""
        ? $(plotLine).attr("tcol") || ""
        : Object.keys(newChart.datas)[0] || "",
    color: $(plotLine).attr("colour")
      ? $(plotLine).attr("colour")
      : lineColors[iPlotLine % 6],
    fillcolour:
      $(plotLine).attr("fillcolour") || ""
        ? $(plotLine).attr("fillcolour") || ""
        : lineColors[iPlotLine % 6],
    linestyle: $(plotLine).find("linestyle").text(),
    label: $(plotLine).find("label").text() || false,
    symbolsize: parseInt($(plotLine).find("symbolsize").text() || "1"),
  };

  finalizePlotLine(newPlotLine, newChart, plotLine, iPlotLine);
  return newPlotLine;
};

const finalizePlotLine = (
  newPlotLine: CCP4i2PlotLine,
  newChart: CCP4i2Chart,
  plotLine: Element,
  iPlotLine: number
) => {
  if (
    !newPlotLine.dataid ||
    !Object.keys(newChart.datas).includes(newPlotLine.dataid)
  ) {
    newPlotLine.dataid = Object.keys(newChart.datas)[0];
  }
};

const parseAxes = (newPlot: CCP4i2Plot, newChart: CCP4i2Chart) => {
  newPlot.plotLines.forEach((newPlotLine) => {
    const xAxisId = newPlotLine.xcol;
    if (!newPlot.xAxes[`${xAxisId}`]) {
      newPlot.xAxes[`${xAxisId}`] = {
        plotLine: newPlotLine,
        xAxisId: xAxisId,
        dataKey: newPlotLine.xcol,
        position: "insideBottom",
        domain: [newPlot.xrange.min, newPlot.xrange.max],
      };
    }
    if (newPlot.xscale === "oneoversqrt") {
      newPlot.xAxes[`${xAxisId}`].tickFormatter = (value: number) => {
        return Math.pow(value, -0.5).toFixed(2);
      };
    }
    newPlotLine.xAxisId = xAxisId;

    const yAxisId = newPlotLine.ycol;
    if (!newPlot.yAxes[`${yAxisId}`]) {
      newPlot.yAxes[`${yAxisId}`] = {
        plotLine: newPlotLine,
        yAxisId: yAxisId,
        dataKey: newPlotLine.ycol,
        position:
          $(newPlotLine).attr("rightaxis") === "true"
            ? "insideRight"
            : "insideLeft",
        orientation:
          $(newPlotLine).attr("rightaxis") === "true" ? "right" : "left",
        angle: -90,
        domain: [newPlot.yrange.min, newPlot.yrange.max],
      };
    }
    newPlotLine.yAxisId = yAxisId;
  });

  Object.keys(newPlot.yAxes).forEach((axisKey) => {
    const axis = newPlot.yAxes[axisKey];
    if (axis.orientation === "right") {
      newPlot.rightYAxes[axisKey] = axis;
    } else {
      newPlot.leftYAxes[axisKey] = axis;
    }
  });

  [newPlot.xAxes, newPlot.leftYAxes, newPlot.rightYAxes].forEach((axisSet) => {
    let globalMin = 1e30;
    let globalMax = -1e30;
    Object.keys(axisSet).forEach((axisKey) => {
      const axis = axisSet[axisKey];
      const relevantData = newChart.datas[axis.plotLine.dataid].map(
        (item: any) => parseFloat(item[axis.dataKey])
      );
      const filteredData = relevantData.filter(
        (entry: number) => !Number.isNaN(entry)
      );
      let axisMin = axis.domain[0];
      if (axisMin === "auto") {
        axisMin = Math.min(1e30, Math.min(...filteredData));
      }
      let axisMax = axis.domain[1];
      if (axisMax === "auto") {
        axisMax = Math.max(-1e30, Math.max(...filteredData));
      }
      globalMin = Math.min(globalMin, axisMin);
      globalMax = Math.max(globalMax, axisMax);
    });
    Object.keys(axisSet).forEach((axisKey) => {
      const axis = axisSet[axisKey];
      axis.domain = [
        globalMin !== 1e30 ? globalMin : "auto",
        globalMax !== -1e30 ? globalMax : "auto",
      ];
    });
  });
};

const parsePlot = (plot: Element, newChart: CCP4i2Chart): CCP4i2Plot => {
  const newPlot: CCP4i2Plot = {
    plottype: $(plot).find("plottype").text().trim(),
    title: $(plot).find("title").text().trim(),
    xLabel: decodeHTMLEntities($(plot).find("xlabel").text().trim()),
    yLabel: decodeHTMLEntities($(plot).find("ylabel").text().trim()),
    ryLabel: decodeHTMLEntities($(plot).find("rylabel").text().trim()),
    plotLines: [],
    barCharts: [],
    yAxes: {},
    xAxes: {},
    leftYAxes: {},
    rightYAxes: {},
    xintegral: $(plot).find("xintegral").text(),
    yintegral: $(plot).find("yintegral").text(),
    xscale: $(plot).find("xscale").text(),
    polygons: $(plot)
      .find("polygon")
      .toArray()
      .map((polygon) => ({
        points: $(polygon)
          .text()
          .trim()
          .split(/\s+|\t+/),
        fill: $(polygon).attr("fillcolour"),
        stroke: $(polygon).attr("linecolour"),
        fillOpacity: $(polygon).attr("alpha"),
      })),
    circles: $(plot)
      .find("circle")
      .toArray()
      .map((circle) => ({
        xpos: parseFloat($(circle).attr("xpos") || "0"),
        ypos: parseFloat($(circle).attr("ypos") || "0"),
        radius: parseFloat($(circle).attr("radius") || "0"),
        linecolour: $(circle).attr("linecolour"),
        alpha: parseFloat($(circle).attr("alpha") || "1"),
      })),
    lines: $(plot)
      .find("line")
      .toArray()
      .map((line) => ({
        x1: parseFloat($(line).attr("x1") || "0"),
        x2: parseFloat($(line).attr("x2") || "0"),
        y1: parseFloat($(line).attr("y1") || "0"),
        y2: parseFloat($(line).attr("y2") || "0"),
        linestyle: $(line).attr("linestyle"),
        linecolour: $(line).attr("linecolour"),
      })),
    fixaspectratio: false,
    xrange: {
      min: $(plot).find("xrange").attr("min")
        ? parseFloat($(plot).find("xrange").attr("min") || "auto")
        : "auto",
      max: $(plot).find("xrange").attr("max")
        ? parseFloat($(plot).find("xrange").attr("max") || "auto")
        : "auto",
    },
    yrange: {
      min: $(plot).find("yrange").attr("min")
        ? parseFloat($(plot).find("yrange").attr("min") || "")
        : "auto",
      max: $(plot).find("yrange").attr("max")
        ? parseFloat($(plot).find("yrange").attr("max") || "")
        : "auto",
    },
  };

  newPlot.plotLines = parsePlotLines(plot, newChart);
  parseAxes(newPlot, newChart);

  return newPlot;
};

export const parseTables = (newTables: Element[]): CCP4i2Chart[] => {
  return newTables.map((item) => {
    const newChart: CCP4i2Chart = {
      plots: [],
      datas: parseData(item),
      headers: parseHeaders(item),
    };

    $(item)
      .find("plot")
      .toArray()
      .forEach((plot) => {
        newChart.plots.push(parsePlot(plot, newChart));
      });

    return newChart;
  });
};

export const decodeHTMLEntities = (str: string) => {
  if (str && typeof str === "string") {
    str = str.replace(/&nbsp;/g, " ");
  }
  return str;
};
