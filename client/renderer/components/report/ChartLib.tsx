import { Chart, ChartData, ChartOptions, ChartTypeRegistry } from "chart.js";
import { parseStringPromise } from "xml2js";

// Extend Chart.js types to include the custom backgroundImage plugin
declare module "chart.js" {
  interface PluginOptionsByType<TType extends keyof ChartTypeRegistry> {
    backgroundImage?: {
      image: HTMLImageElement;
    };
  }
}

const colours = [
  "rgba(255, 99, 132, 1)",
  "rgba(54, 162, 235, 1)",
  "rgba(255, 206, 86, 1)",
  "rgba(75, 192, 192, 1)",
  "rgba(153, 102, 255, 1)",
  "rgba(255, 159, 64, 1)",
];

export interface CCP4ApplicationOutput {
  CCP4Table?: CCP4Table | CCP4Table[];
  Fonts?: Fonts;
  CCP4Surface?: Surface;
  title?: string;
}

interface CCP4TableInterface {
  headers?: Header;
  data?: Data | Data[];
  plot?: Plot | Plot[];
  title?: string;
}

export class CCP4Table implements CCP4TableInterface {
  headers?: Header;
  data?: Data | Data[];
  plot?: Plot | Plot[];
  title?: string;
  _parsedDataBlocks?: any;
  constructor(xmlFromJson: CCP4Table) {
    this.headers = xmlFromJson.headers;
    this.data = xmlFromJson.data;
    this.plot = xmlFromJson.plot;
    this.title = xmlFromJson.title;
    this._parsedDataBlocks = parsedDataBlocks(this);
  }
  get parsedDataBlocks(): any {
    return this._parsedDataBlocks;
  }
  get allPlots(): null | Plot[] {
    if (!this.plot) return null;
    if (Array.isArray(this.plot)) return this.plot;
    else return [this.plot];
  }
  get allHeaders(): null | string[] {
    if (!this.headers) return null;
    let headers: string[];
    if (typeof this.headers === "string" || this.headers instanceof String) {
      headers = this.headers.split(" ");
    } else {
      //console.log("Non string headers: ", this.headers);
      let separator: string | RegExp | undefined = this.headers.separator;
      if (!separator) separator = /\s+/;
      headers = this.headers["_"]
        .split(separator)
        .map((header: string) => header.trim())
        .filter((header: string) => header.length > 0);
    }
    //console.log("Headers is", { headers });
    return headers;
  }

  plotData(selectedPlot: Plot): ChartData<"bar" | "scatter"> | null {
    if (!selectedPlot || !this.parsedDataBlocks || !this.allHeaders)
      return null;
    const datasets = extractDatasets(
      selectedPlot,
      this.parsedDataBlocks,
      this.allHeaders
    );
    if (!datasets) return null;
    const result: any = {
      datasets: datasets,
    };
    return result;
  }

  plotOptions(selectedPlot: Plot): null | ChartOptions {
    if (!selectedPlot) return null;
    const result: ChartOptions = {
      animation: false,
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        title: {
          display: true,
          text: this.title,
        },
      },
    };

    installDefaultScales(selectedPlot, result);

    checkForRightAxis(selectedPlot, result);

    handleRangeSpecifiers(selectedPlot, result);

    handleLegend(selectedPlot, result);

    if (selectedPlot.polygon) {
      addPolygons(selectedPlot, result);
    }

    if (selectedPlot.circle) {
      addCircles(selectedPlot, result);
    }

    if (selectedPlot.line) {
      addLines(selectedPlot, result);
    }

    if (selectedPlot?.text) {
      addTexts(selectedPlot, result);
    }

    if (selectedPlot?.background) {
      result.plugins = result.plugins || {};

      // Extract path from the full URL
      let imagePath: string;
      try {
        const url = new URL(selectedPlot.background);
        imagePath = url.pathname; // This gets just the path part
      } catch (error) {
        // Fallback if URL parsing fails - assume it's already a path
        imagePath = selectedPlot.background;
      }
      const img = new Image();
      img.src = imagePath;
      result.plugins.backgroundImage = {
        image: img,
      };
    }

    if (selectedPlot?.xscale === "oneoversqrt") {
      handleOneOverSqrt(selectedPlot, result);
    }

    //Custom tick scenarios
    if (selectedPlot?.xintegral) {
      handleXIntegral(selectedPlot, result);
    }

    if (selectedPlot?.customXLabels) {
      handleCustomXLabels(selectedPlot, result);
    }

    if (
      (selectedPlot?.xlabel === "Phi" && selectedPlot?.ylabel === "Psi") ||
      selectedPlot?.fixaspectratio === "true"
    ) {
      result.maintainAspectRatio = true;
      result.aspectRatio = 1;
    }

    //console.log({ options: result });
    return result;
  }
}

export interface Header {
  _: string;
  separator?: string;
}

export interface Data {
  _: string;
  separator?: string;
  id?: string;
}

export class Plot {
  title?: string;
  plotline?: PlotLine[] | PlotLine;
  histogram?: Histogram[];
  barchart?: BarChart[];
  xlabel?: string;
  ylabel?: string;
  rylabel?: string;
  description?: string;
  xintegral?: string;
  xrange?: { min?: number; max?: number };
  yrange?: { min?: number; max?: number };
  polygon?: any[] | any;
  circle?: any[] | any;
  line?: any[] | any;
  text?: any[] | any;
  xscale?: string;
  customXTicks?: string;
  customXLabels?: string;
  showlegend?: "true" | "false";
  fixaspectratio?: "true" | "false";
  background?: string;
}

export interface PlotLine {
  xcol: number;
  ycol: number;
  dataid?: string;
  rightaxis?: string;
  colour?: string;
  linestyle?: string;
  showlegend?: "true" | "false";
}

export interface Histogram {
  col: number;
  colour?: string;
  nbins?: number;
  binwidth?: number;
}

export interface BarChart {
  col?: number;
  tcol?: number;
  colour?: string;
  width?: string;
  rightaxis?: string;
  dataid?: string;
}

export interface Fonts {
  titleFont?: FontLine;
  legendFont?: FontLine;
  axesTickerFont?: FontLine;
}

export interface FontLine {
  _: string;
  family?: string;
  size?: number;
  weight?: string;
  slant?: string;
}

export interface Surface {
  _: string;
  rows?: number;
  columns?: number;
  title?: string;
}

export const parseXML = async (xml: string): Promise<CCP4ApplicationOutput> => {
  const nsStripped = stripNamespaces(xml);
  const tablised = changeTagName(nsStripped, "ccp4_data", "CCP4Table");
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(tablised, "application/xml");
  const oldElements = xmlDoc.getElementsByTagName("CCP4Table");
  const firstTable = oldElements[0];
  // Serialize back to string
  const serializer = new XMLSerializer();
  const firstTableString = serializer.serializeToString(firstTable);

  const result = await parseStringPromise(firstTableString, {
    explicitArray: false,
    mergeAttrs: true,
  });
  return result as CCP4ApplicationOutput;
};

function stripNamespaces(xmlString: string): string {
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(xmlString, "application/xml");
  function processNode(node: Node): void {
    if (node.nodeType === Node.ELEMENT_NODE) {
      const element = node as Element;

      // Create a new element without namespace
      const newElement = xmlDoc.createElement(element.localName);

      // Copy attributes (without namespace)
      Array.from(element.attributes).forEach((attr) => {
        newElement.setAttribute(attr.name, attr.value);
      });

      // Move child nodes
      while (element.firstChild) {
        newElement.appendChild(element.firstChild);
      }

      // Replace the old element with the new one
      const parent = element.parentNode;
      if (parent) {
        parent.replaceChild(newElement, element);
        // Process child nodes
        newElement.childNodes.forEach(processNode);
      }
    }
  }

  processNode(xmlDoc);

  // Serialize back to string
  const serializer = new XMLSerializer();
  return serializer.serializeToString(xmlDoc);
}

/**
 * Changes the tag name of elements in an XML string within a specified namespace.
 *
 * @param xmlString - The XML string to be modified.
 * @param namespaceURI - The namespace URI of the elements to be changed.
 * @param oldTagName - The current tag name of the elements to be changed.
 * @param newTagName - The new tag name to replace the old tag name.
 * @returns The modified XML string with the tag names changed.
 */
function changeTagNameNS(
  xmlString: string,
  namespaceURI: string,
  oldTagName: string,
  newTagName: string
): string {
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(xmlString, "application/xml");
  //console.log({ xmlDoc });
  // Get all elements with the old tag name in the namespace
  const oldElements = xmlDoc.getElementsByTagNameNS(namespaceURI, oldTagName);
  //console.log({ oldElements });

  // Convert NodeList to array (since NodeList is live)
  const elementsArray: Element[] = Array.from(oldElements);

  elementsArray.forEach((oldElement) => {
    // Create a new element with the same namespace
    const newElement = xmlDoc.createElementNS(namespaceURI, newTagName);

    // Copy attributes
    Array.from(oldElement.attributes).forEach((attr) => {
      newElement.setAttributeNS(attr.namespaceURI, attr.name, attr.value);
    });

    // Move child nodes to the new element
    while (oldElement.firstChild) {
      newElement.appendChild(oldElement.firstChild);
    }

    // Replace the old element with the new one
    const parent = oldElement.parentNode;
    if (parent) {
      parent.replaceChild(newElement, oldElement);
    }
  });

  // Serialize XML back to string
  const serializer = new XMLSerializer();
  return serializer.serializeToString(xmlDoc);
}

/**
 * Changes the tag name of elements in an XML string.
 *
 * @param {string} xmlString - The XML string to process.
 * @param {string} oldTagName - The old tag name to be replaced.
 * @param {string} newTagName - The new tag name to replace with.
 * @param {string} [namespaceURI] - Optional namespace URI to match and preserve.
 * @returns {string} - The modified XML string with updated tag names.
 *
 * @example
 * ```typescript
 * const xml = `<root><oldTag>Content</oldTag></root>`;
 * const result = changeTagName(xml, "oldTag", "newTag");
 * console.log(result); // <root><newTag>Content</newTag></root>
 * ```
 */
function changeTagName(
  xmlString: string,
  oldTagName: string,
  newTagName: string,
  namespaceURI?: string
): string {
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(xmlString, "application/xml");
  //console.log("In changeTagName", { xmlDoc });
  function processNode(node: Node): void {
    if (node.nodeType === Node.ELEMENT_NODE) {
      const element = node as Element;

      // Check if the element matches the old tag name and optional namespace
      if (
        element.localName === oldTagName &&
        (!namespaceURI || element.namespaceURI === namespaceURI)
      ) {
        const newElement = namespaceURI
          ? xmlDoc.createElementNS(namespaceURI, newTagName) // Preserve namespace if provided
          : xmlDoc.createElement(newTagName); // Otherwise, create without namespace

        // Copy attributes
        Array.from(element.attributes).forEach((attr) => {
          newElement.setAttribute(attr.name, attr.value);
        });

        // Move child nodes
        while (element.firstChild) {
          newElement.appendChild(element.firstChild);
        }

        // Replace the old element with the new one
        const parent = element.parentNode;
        if (parent) {
          parent.replaceChild(newElement, element);
          processNode(newElement); // Recursive call to process children
        }
      }
    }

    // Process child nodes
    node.childNodes.forEach(processNode);
  }

  processNode(xmlDoc.documentElement);

  // Serialize back to string
  const serializer = new XMLSerializer();
  return serializer.serializeToString(xmlDoc);
}

const handleLegend = (selectedPlot: Plot, result: ChartOptions) => {
  const plotlines = Array.isArray(selectedPlot.plotline)
    ? selectedPlot.plotline
    : selectedPlot.plotline
    ? [selectedPlot.plotline]
    : [];
  if (
    plotlines.filter((plotline: PlotLine) => plotline.showlegend === "false")
      .length > 0
  ) {
    if (!result.plugins) result.plugins = {};
    result.plugins.legend = {
      display: false,
    };
  } else {
    if (!result.plugins) result.plugins = {};
    result.plugins.legend = {
      display: true,
      position: "chartArea",
    };
  }
};

/**
 * Configures the default scales for a given chart.
 *
 * @param {Plot} selectedPlot - The plot object containing the labels for the axes.
 * @param {ChartOptions} result - The chart options object where the scales will be set.
 *
 * The function sets up the x-axis and y-axis scales with the following properties:
 * - `x` axis:
 *   - `axis`: "x"
 *   - `type`: "linear"
 *   - `position`: "bottom"
 *   - `title`: Displays the x-axis label from `selectedPlot.xlabel`
 *   - `ticks`: Formats the tick values, displaying integers as they are and non-integers with a precision of 3
 * - `yAxisLeft` axis:
 *   - `axis`: "y"
 *   - `type`: "linear"
 *   - `position`: "left"
 *   - `title`: Displays the y-axis label from `selectedPlot.ylabel`
 *   - `ticks`: Formats the tick values, displaying integers as they are and non-integers with a precision of 3
 */
const installDefaultScales = (selectedPlot: Plot, result: ChartOptions) => {
  result.scales = {
    x: {
      axis: "x",
      type: "linear",
      position: "bottom",
      title: {
        display: true,
        text: selectedPlot?.xlabel ? selectedPlot.xlabel : "",
      },
      ticks: {
        callback: (value: string | number, index: number) => {
          if (typeof value === "string") {
            return value;
          } else if (Number.isInteger(value)) {
            return value;
          }
          return value.toPrecision(3);
        },
      },
    },
    yAxisLeft: {
      axis: "y",
      type: "linear",
      position: "left",
      title: {
        display: true,
        text: selectedPlot?.ylabel ? selectedPlot.ylabel : "",
      },
      ticks: {
        callback: (value: string | number, index: number) => {
          if (typeof value === "string") {
            return value;
          } else if (Number.isInteger(value)) {
            return value;
          }
          return value.toPrecision(3);
        },
      },
    },
  };
};

/**
 * Checks for occurrence of plotLines which are specified to have a right axis.
 *
 * @param {Plot} selectedPlot - The plot object containing plotline information.
 * @param {ChartOptions} result - The chart options object to which plot lines will be added.
 *
 * @remarks
 * - If the selected plot does not have any plotline, the function returns early.
 * - The function processes the plotline(s) and filters those that should be displayed on the right axis.
 * - If there are plotlines for the right axis, it ensures the `scales` property exists in the result and adds the right axis configuration.
 * - The right axis configuration includes axis type, position, grid display, and title.
 * - The ticks for the right axis are configured to display integer values as they are and non-integer values with a precision of 3 significant digits.
 */
export const checkForRightAxis = (selectedPlot: Plot, result: ChartOptions) => {
  if (!selectedPlot.plotline) return;
  const plotlineList: PlotLine[] = Array.isArray(selectedPlot?.plotline)
    ? (selectedPlot.plotline as PlotLine[])
    : [selectedPlot.plotline as PlotLine];
  const rightAxisPlotlines = plotlineList.filter(
    (plotline: PlotLine) => plotline.rightaxis === "true"
  );
  if (rightAxisPlotlines.length > 0) {
    if (!result.scales) result.scales = {};
    result.scales.yAxisRight = {
      axis: "y",
      type: "linear",
      position: "right",
      grid: { display: false },
      title: {
        display: true,
        text: selectedPlot?.rylabel ? selectedPlot.rylabel : "",
      },
      ticks: {
        callback: (value: string | number, index: number) => {
          if (typeof value === "string") {
            return value;
          } else if (Number.isInteger(value)) {
            return value;
          }
          return value.toPrecision(3);
        },
      },
    };
  }
};

/**
 * Decodes HTML entities in a string (e.g., &nbsp; -> space, &lt; -> <)
 * @param text - The text containing HTML entities
 * @returns The decoded text
 */
const decodeHTMLEntities = (text: string): string => {
  const textarea = document.createElement('textarea');
  textarea.innerHTML = text;
  return textarea.value;
};

/**
 * Extracts the datasets for a chart from the provided data blocks.
 *
 * @param {Plot} selectedPlot - The plot object containing the chart configuration.
 * @param {object} parsedDataBlocks - The parsed data blocks from the CCP4 XML output.
 * @param {string[]} allHeaders - An array of all header strings.
 * @returns {any[]} An array of objects representing the datasets for the chart.
 */
export const extractDatasets = (
  selectedPlot: Plot,
  parsedDataBlocks: any,
  allHeaders: string[]
) => {
  if (!selectedPlot.plotline && !selectedPlot.barchart) return null;
  const nPlotlines: number = Array.isArray(selectedPlot.plotline)
    ? selectedPlot.plotline.length
    : selectedPlot.plotline
    ? 1
    : 0;
  const nBarcharts: number = Array.isArray(selectedPlot.barchart)
    ? selectedPlot.barchart.length
    : selectedPlot.barchart
    ? 1
    : 0;
  if (nPlotlines + nBarcharts === 0) return null;

  const plotlineList = Array.isArray(selectedPlot.plotline)
    ? selectedPlot.plotline
    : selectedPlot.plotline
    ? [selectedPlot.plotline]
    : [];

  const plotlineDatasets = plotlineList.map(
    (plotline: PlotLine, iPlotline: number) => {
      let dataAsGrid: any[][] = [];
      if (
        plotline.dataid &&
        Object.keys(parsedDataBlocks).includes(plotline.dataid)
      ) {
        dataAsGrid = parsedDataBlocks[plotline.dataid];
      } else if (Object.keys(parsedDataBlocks).includes("_")) {
        dataAsGrid = parsedDataBlocks["_"];
      } else return null;

      const dataset: any = extractPlotLineDataset(
        dataAsGrid,
        allHeaders,
        plotline,
        iPlotline
      );
      return dataset;
    }
  );

  const barChartList = Array.isArray(selectedPlot.barchart)
    ? selectedPlot.barchart
    : selectedPlot.barchart
    ? [selectedPlot.barchart]
    : [];
  const barChartDatasets = barChartList.map(
    (barChart: BarChart, iBarChart: number) => {
      let dataAsGrid: any[][] = [];
      if (
        barChart.dataid &&
        Object.keys(parsedDataBlocks).includes(barChart.dataid)
      ) {
        dataAsGrid = parsedDataBlocks[barChart.dataid];
      } else if (Object.keys(parsedDataBlocks).includes("_")) {
        dataAsGrid = parsedDataBlocks["_"];
      } else return null;

      const dataset: any = extractBarChartDataset(
        dataAsGrid,
        allHeaders,
        barChart,
        iBarChart
      );
      return dataset;
    }
  );

  return plotlineDatasets.concat(barChartDatasets);
};

/**
 * Extracts a dataset for a plot line from the provided grid data.
 *
 * @param dataAsGrid - A 2D array representing the grid data.
 * @param allHeaders - An array of all header strings.
 * @param plotline - An object representing the plot line configuration.
 * @param iPlotline - The index of the current plot line.
 * @returns An object representing the dataset for the plot line.
 */
export const extractPlotLineDataset = (
  dataAsGrid: any[][],
  allHeaders: string[],
  plotline: PlotLine,
  iPlotline: number
) => {
  // Decode HTML entities in label (e.g., &nbsp; -> space)
  const rawLabel = allHeaders[plotline.ycol - 1];
  const decodedLabel = rawLabel ? decodeHTMLEntities(rawLabel) : rawLabel;

  const result = {
    label: decodedLabel,
    labels: dataAsGrid.map((row: any) => row[parseInt(`${plotline.xcol}`) - 1]),
    yAxisID: plotline.rightaxis
      ? plotline.rightaxis === "true"
        ? "yAxisRight"
        : "yAxisLeft"
      : "yAxisLeft",
    data: dataAsGrid.map((row) => ({
      x: row[parseInt(`${plotline.xcol}`) - 1],
      y: row[parseInt(`${plotline.ycol}`) - 1],
    })),
    backgroundColor: plotline?.colour?.startsWith
      ? plotline.colour.startsWith("#")
        ? hexToRGBA(plotline.colour, 0.5)
        : colorNameToRGBA(plotline.colour, 0.5)
      : colorNameToRGBA(colours[iPlotline % colours.length], 0.5),
    borderColor: plotline.colour
      ? plotline.colour
      : colours[iPlotline % colours.length],
    showLine: plotline.linestyle !== ".",
  };
  return result;
};

/**
 * Extracts a dataset for a bar chart from the provided grid data.
 *
 * @param dataAsGrid - A 2D array representing the grid data.
 * @param allHeaders - An array of all header strings.
 * @param barChart - An object representing the bar chart configuration.
 * @param iBarChart - The index of the current bar chart.
 * @returns An object representing the dataset for the bar chart, or null if the required columns are not specified.
 */
export const extractBarChartDataset = (
  dataAsGrid: any[][],
  allHeaders: string[],
  barChart: BarChart,
  iBarChart: number
) => {
  if (!barChart?.col) return null;
  if (!barChart?.tcol) return null;

  // Decode HTML entities in label (e.g., &nbsp; -> space)
  const rawLabel = allHeaders[barChart.tcol - 1];
  const decodedLabel = rawLabel ? decodeHTMLEntities(rawLabel) : rawLabel;

  const result = {
    label: decodedLabel,
    type: "bar",
    labels: dataAsGrid.map((row: any) => row[parseInt(`${barChart.col}`) - 1]),
    yAxisID: barChart.rightaxis
      ? barChart.rightaxis === "true"
        ? "yAxisRight"
        : "yAxisLeft"
      : "yAxisLeft",
    data: dataAsGrid.map((row) => ({
      x: row[parseInt(`${barChart.col}`) - 1],
      y: row[parseInt(`${barChart.tcol}`) - 1],
    })),
    backgroundColor: barChart.colour
      ? barChart.colour.startsWith("#")
        ? hexToRGBA(barChart.colour, 0.5)
        : colorNameToRGBA(barChart.colour, 0.5)
      : colorNameToRGBA(colours[iBarChart % colours.length], 0.5),
    borderColor: barChart.colour
      ? barChart.colour
      : colours[iBarChart % colours.length],
    showLine: true,
  };
  return result;
};

/**
 * Adds lines to the chart options for annotation purposes.
 *
 * @param selectedPlot - The plot object containing line data.
 * @param result - The chart options object to which the lines will be added.
 *
 * The function checks if the selected plot has line annotations. If it does, it processes
 * each line annotation and adds it to the chart options under the `plugins.annotation.annotations` property.
 * Each line annotation is configured with properties such as start and end points, draw time, and color.
 */
export const addLines = (selectedPlot: Plot, result: ChartOptions): void => {
  if (!selectedPlot.line) return;
  let lines: any[] = [];
  if (Array.isArray(selectedPlot.line)) {
    lines = selectedPlot.line;
  } else {
    lines = [selectedPlot.line];
  }
  lines.forEach((line: any, iLine: number) => {
    if (!result.plugins) result.plugins = {};
    if (!result.plugins.annotation) result.plugins.annotation = {};
    if (!result.plugins.annotation.annotations)
      result.plugins.annotation.annotations = {};
    const lineObject = {
      type: "line",
      xMin: parseFloat(line.x1),
      yMin: parseFloat(line.y1),
      xMax: parseFloat(line.x2),
      yMax: parseFloat(line.y2),
      drawTime: "beforeDatasetsDraw",
      borderColor: line.linecolour
        ? line.alpha
          ? line.fillcolour.startsWith("#")
            ? hexToRGBA(line.fillcolour, parseFloat(line.alpha))
            : colorNameToRGBA(line.fillcolour, parseFloat(line.alpha))
          : line.fillcolour
        : "rgba(0,0,0,0)",
    };
    (result.plugins.annotation.annotations as Record<string, any>)[
      `line-${iLine}`
    ] = lineObject;
  });
};

/**
 * Adds polygons to the chart options for annotation purposes.
 *
 * @param {Plot} selectedPlot - The plot object containing polygon data.
 * @param {ChartOptions} result - The chart options object to which the polygons will be added.
 */
export const addPolygons = (selectedPlot: Plot, result: ChartOptions) => {
  let polygons: any[] = [];
  if (Array.isArray(selectedPlot.polygon)) {
    polygons = selectedPlot.polygon;
  } else {
    polygons = [selectedPlot.polygon];
  }
  polygons.forEach((polygon: any, iPolygon: number) => {
    if (!result.plugins) result.plugins = {};
    if (!result.plugins.annotation) result.plugins.annotation = {};
    if (!result.plugins.annotation.annotations)
      result.plugins.annotation.annotations = {};
    const boxObject = {
      type: "box",
      xMax: Math.max(
        ...polygon["_"]
          .split(" ")
          .filter((item: string, iItem: number) => !(iItem % 2))
          .map((item: string) => parseFloat(item))
      ),
      xMin: Math.min(
        ...polygon["_"]
          .split(" ")
          .filter((item: string, iItem: number) => !(iItem % 2))
          .map((item: string) => parseFloat(item))
      ),
      yMax: Math.max(
        ...polygon["_"]
          .split(" ")
          .filter((item: string, iItem: number) => iItem % 2)
          .map((item: string) => parseFloat(item))
      ),
      yMin: Math.min(
        ...polygon["_"]
          .split(" ")
          .filter((item: string, iItem: number) => iItem % 2)
          .map((item: string) => parseFloat(item))
      ),
      drawTime: "beforeDatasetsDraw",
      backgroundColor: polygon.fillcolour
        ? polygon.alpha
          ? polygon.fillcolour.startsWith("#")
            ? hexToRGBA(polygon.fillcolour, parseFloat(polygon.alpha))
            : colorNameToRGBA(polygon.fillcolour, parseFloat(polygon.alpha))
          : polygon.fillcolour
        : "rgba(0,0,0,0)",
      strokeColor: polygon.linecolour || "rgba(0,0,0,1)",
    };
    (result.plugins.annotation.annotations as Record<string, any>)[
      `polygon-${iPolygon}`
    ] = boxObject;
  });
};

/**
 * Adds circle annotations to the given chart options based on the selected plot.
 *
 * @param {Plot} selectedPlot - The plot containing circle data to be added as annotations.
 * @param {ChartOptions} result - The chart options object where the circle annotations will be added.
 */
export const addCircles = (selectedPlot: Plot, result: ChartOptions) => {
  let circles: any[] = [];
  if (Array.isArray(selectedPlot.circle)) {
    circles = selectedPlot.circle;
  } else {
    circles = [selectedPlot.circle];
  }
  circles.forEach((circle: any, iCircle: number) => {
    if (!result.plugins) result.plugins = {};
    if (!result.plugins.annotation) result.plugins.annotation = {};
    if (!result.plugins.annotation.annotations)
      result.plugins.annotation.annotations = {};
    const circleObject = {
      type: "ellipse",
      xMax: parseFloat(circle.xpos) + parseFloat(circle.radius),
      xMin: parseFloat(circle.xpos) - parseFloat(circle.radius),
      yMax: parseFloat(circle.ypos) + parseFloat(circle.radius),
      yMin: parseFloat(circle.ypos) - parseFloat(circle.radius),
      drawTime: "beforeDatasetsDraw",
      backgroundColor: circle.fillcolour
        ? circle.alpha
          ? circle.fillcolour.startsWith("#")
            ? hexToRGBA(circle.fillcolour, parseFloat(circle.alpha))
            : colorNameToRGBA(circle.fillcolour, parseFloat(circle.alpha))
          : circle.fillcolour
        : "rgba(0,0,0,0)",
      strokeColor: circle.linecolour,
    };
    (result.plugins.annotation.annotations as Record<string, any>)[
      `circle-${iCircle}`
    ] = circleObject;
  });
};

/**
 * Adds text annotations to the given chart options based on the selected plot.
 *
 * @param selectedPlot - The plot object containing text annotations to be added.
 * @param result - The chart options object where the text annotations will be added.
 *
 * The function checks if the selected plot has text annotations. If it does, it processes
 * each text annotation and adds it to the chart options under the `plugins.annotation.annotations` property.
 * Each text annotation is configured with properties such as position, background color, content, color, font, padding, and border radius.
 */
export const addTexts = (selectedPlot: Plot, result: ChartOptions) => {
  if (selectedPlot?.text) {
    let texts: any[] = [];
    if (Array.isArray(selectedPlot.text)) {
      texts = selectedPlot.text;
    } else {
      texts = [selectedPlot.text];
    }
    texts.forEach((text, iText) => {
      const textObject = {
        type: "label",
        xValue: parseFloat(text.xpos), // X position (match a data point)
        yValue: parseFloat(text.ypos), // Y position
        backgroundColor: "rgba(0,0,0, 0.0)",
        content: text["_"],
        color: text.colour,
        font: {
          size: 14,
          weight: "bold",
        },
        padding: 6,
        borderRadius: 4,
      };
      if (!result.plugins) result.plugins = {};
      if (!result.plugins.annotation) result.plugins.annotation = {};
      if (!result.plugins.annotation.annotations)
        result.plugins.annotation.annotations = {};
      (result.plugins.annotation.annotations as Record<string, any>)[
        `text-${iText}`
      ] = textObject;
    });
  }
};

/**
 * Adjusts the chart options for a scatter plot to display the x-axis values as the inverse of their square roots.
 *
 * @param selectedPlot - The plot object that is selected.
 * @param result - The chart options for a scatter plot that will be modified.
 *
 * This function modifies the provided `result` object to:
 * - Set the x-axis tick labels to the inverse of the square root of the original values, formatted to 3 significant digits.
 * - Customize the tooltip to display the x-axis value as the inverse of its square root, along with the dataset label and y-axis value.
 */
export const handleOneOverSqrt = (selectedPlot: Plot, result: ChartOptions) => {
  if (!result.scales) result.scales = {};
  if (!result.scales.x) result.scales.x = {};
  result.scales.x.ticks = {
    callback: (value: any, index: number) =>
      (1 / Math.sqrt(value)).toPrecision(3), // Hide non-integer labels
  };
  if (!result.plugins) result.plugins = {};
  if (!result.plugins) result.plugins = {};
  result.plugins.tooltip = {
    callbacks: {
      label: (tooltipItem: any) =>
        `Res: ${(1 / Math.sqrt(tooltipItem.raw.x)).toPrecision(3)}, ${
          tooltipItem.dataset.label
        }: ${tooltipItem.raw.y}`,
    },
  };
};

/**
 * Adjusts the x-axis scale of a scatter chart to display integer steps if the selected plot requires it.
 *
 * @param selectedPlot - The plot object containing configuration details.
 * @param result - The chart options object to be modified.
 *
 * @remarks
 * This function modifies the `result` object to ensure that the x-axis displays integer steps
 * if the `xintegral` property of the `selectedPlot` is set to "true". It sets the `stepSize` to 1
 * and hides non-integer labels on the x-axis.
 *
 * @example
 * ```typescript
 * const selectedPlot = { xintegral: "true" };
 * const result = { scales: { x: {} } };
 * handleXIntegral(selectedPlot, result);
 * // result.scales.x.ticks will be set to force integer steps and hide non-integer labels
 * ```
 */
export const handleXIntegral = (selectedPlot: Plot, result: ChartOptions) => {
  if (!result.scales) result.scales = {};
  if (!result.scales.x) result.scales.x = {};
  if (selectedPlot.xintegral && selectedPlot.xintegral === "true") {
    result.scales.x.ticks = {
      stepSize: 1, // Force integer steps
      callback: (value: any) => (Number.isInteger(value) ? value : null), // Hide non-integer labels
    };
  }
};

/**
 * Updates the x-axis tick labels of a scatter chart with custom labels.
 *
 * @param selectedPlot - The plot object containing custom X labels.
 * @param result - The chart options object for a scatter chart.
 *
 * @remarks
 * If `selectedPlot.customXLabels` is defined, this function will split the labels by commas
 * and assign them to the x-axis ticks based on their index. If not defined, it will return an empty string.
 *
 * @example
 * ```typescript
 * const plot = { customXLabels: "Label1,Label2,Label3" };
 * const chartOptions = { scales: { x: { ticks: {} } } };
 * handleCustomXLabels(plot, chartOptions);
 * // chartOptions.scales.x.ticks.callback will now use the custom labels
 * ```
 */
export const handleCustomXLabels = (
  selectedPlot: Plot,
  result: ChartOptions
) => {
  //@ts-ignore
  result.scales.x.ticks = {
    callback: (value: string | number, index: number) => {
      if (selectedPlot.customXLabels) {
        return selectedPlot.customXLabels.split(",")[index]; // Hide non-integer labels
      }
      return "";
    },
  };
};

/**
 * Updates the `result` ChartOptions object with the range specifiers from the `selectedPlot`.
 *
 * @param selectedPlot - The plot object containing the range specifiers for x and y axes.
 * @param result - The ChartOptions object to be updated with the range specifiers.
 *
 * The function checks for the presence of `xrange` and `yrange` properties in the `selectedPlot` object.
 * If these properties are present, it updates the corresponding `min` and `max` values in the `result` object.
 *
 * - If `selectedPlot.xrange.min` is defined, it sets `result.scales.x.min` to this value.
 * - If `selectedPlot.xrange.max` is defined, it sets `result.scales.x.max` to this value.
 * - If `selectedPlot.yrange.min` is defined, it sets `result.scales.yAxisLeft.min` to this value.
 * - If `selectedPlot.yrange.max` is defined, it sets `result.scales.yAxisLeft.max` to this value.
 */
export const handleRangeSpecifiers = (
  selectedPlot: Plot,
  result: ChartOptions
) => {
  if (selectedPlot?.xrange?.min) {
    if (!result?.scales) result.scales = {};
    if (!result?.scales.x) result.scales.x = {};
    result.scales.x.min = selectedPlot.xrange.min;
  }
  if (selectedPlot?.xrange?.max) {
    if (!result?.scales) result.scales = {};
    if (!result?.scales.x) result.scales.x = {};
    result.scales.x.max = selectedPlot.xrange.max;
  }
  if (selectedPlot?.yrange?.min) {
    if (!result?.scales) result.scales = {};
    if (!result?.scales.yAxisLeft) result.scales.yAxisLeft = {};
    result.scales.yAxisLeft.min = selectedPlot.yrange.min;
  }
  if (selectedPlot?.yrange?.max) {
    if (!result?.scales) result.scales = {};
    if (!result?.scales.yAxisLeft) result.scales.yAxisLeft = {};
    result.scales.yAxisLeft.max = selectedPlot.yrange.max;
  }
};

/**
 * Parses the data blocks from a CCP4Table graph object.
 *
 * @param {CCP4Table} graph - The graph object containing data blocks to be parsed.
 * @returns {object | null} - An object where each key is a data block ID and the value is a 2D array of parsed rows, or null if the graph is not provided.
 *
 * The function performs the following steps:
 * 1. Checks if the graph object is provided.
 * 2. Initializes an empty object to store data blocks.
 * 3. Checks if the graph data is an array, if not, wraps it in an array.
 * 4. Iterates over each data block and processes it:
 *    - If the data block is an object, extracts the actual rows and the data block ID.
 *    - Removes leading newline characters from the actual rows.
 *    - Splits the rows by newline characters and then splits each row by whitespace, filtering out empty items.
 *    - Stores the parsed rows in the parsedDataBlocks object with the data block ID as the key.
 * 5. Returns the parsedDataBlocks object.
 * 6. Returns null if the graph object is not provided.
 */
const parsedDataBlocks = (graph: CCP4Table) => {
  if (graph) {
    let dataBlocks: any = {};
    if (Array.isArray(graph.data)) {
      dataBlocks = graph.data;
    } else {
      dataBlocks = [graph.data];
    }
    let parsedDataBlocks: any = {};

    dataBlocks.forEach((dataBlock: any) => {
      let actualRows: string = dataBlock;
      let dataBlockId: string = "_";
      if (typeof dataBlock === "object") {
        actualRows = dataBlock["_"];
        if (dataBlock.id) {
          dataBlockId = dataBlock.id;
        }
      }
      const data = actualRows.replace(/^\n/, ""); // Remove leading newline;
      const blockRows = data.split("\n");
      const parsedBlockRows: any[][] = blockRows.map((row: string) =>
        row.split(/\s+/).filter((item: string) => item.trim().length > 0)
      );
      parsedDataBlocks[dataBlockId] = parsedBlockRows;
    });
    return parsedDataBlocks;
  }
  return null;
};

/**
 * Converts a color name to an RGBA string with the specified alpha value.
 *
 * @param colorName - The name of the color to convert (e.g., "red", "blue").
 * @param alpha - The alpha value for the RGBA color (default is 0.5). Must be between 0 and 1.
 * @returns The RGBA string representation of the color.
 * @throws Will throw an error if the alpha value is not between 0 and 1 or if the color name is invalid.
 */
export function colorNameToRGBA(
  colorName: string,
  alpha: number = 0.5
): string {
  // Ensure alpha value is between 0 and 1
  if (alpha < 0 || alpha > 1) {
    throw new Error("Alpha value must be between 0 and 1");
  }

  // Create a temporary element to access the computed color value
  const tempElement = document.createElement("div");
  tempElement.style.color = colorName;
  document.body.appendChild(tempElement);

  // Get the computed RGB color
  if (typeof window !== "undefined") {
    const rgb = window.getComputedStyle(tempElement).color;

    // Remove the temporary element from the DOM
    document.body.removeChild(tempElement);

    // rgb will be in format "rgb(r, g, b)"
    const rgbValues = rgb.match(/\d+/g); // Extracts the numeric values

    if (!rgbValues) {
      throw new Error("Invalid color name");
    }

    // Convert the extracted values to an RGBA string
    return `rgba(${rgbValues[0]}, ${rgbValues[1]}, ${rgbValues[2]}, ${alpha})`;
  }
  return "";
}

/**
 * Converts a hex color code to an RGBA color string.
 *
 * @param {string} hex - The hex color code in the format #RRGGBB or #RGB.
 * @param {number} [alpha=0.5] - The alpha value for the RGBA color, between 0 and 1.
 * @returns {string} The RGBA color string.
 * @throws {Error} If the alpha value is not between 0 and 1.
 * @throws {Error} If the hex color code is invalid.
 */
export function hexToRGBA(hex: string, alpha: number = 0.5): string {
  // Ensure alpha value is between 0 and 1
  if (alpha < 0 || alpha > 1) {
    throw new Error("Alpha value must be between 0 and 1");
  }

  // Check if the hex is in the valid format (#RRGGBB or #RGB)
  const hexRegex = /^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/;
  if (!hexRegex.test(hex)) {
    throw new Error("Invalid hex color");
  }

  // If the hex is 3 characters, convert it to 6 characters (e.g. #FFF -> #FFFFFF)
  if (hex.length === 4) {
    hex = `#${hex[1]}${hex[1]}${hex[2]}${hex[2]}${hex[3]}${hex[3]}`;
  }

  // Extract RGB values from the hex color
  const r = parseInt(hex.slice(1, 3), 16);
  const g = parseInt(hex.slice(3, 5), 16);
  const b = parseInt(hex.slice(5, 7), 16);

  // Return the RGBA string
  return `rgba(${r}, ${g}, ${b}, ${alpha})`;
}

// Improved custom plugin to add background image
const backgroundImagePlugin = {
  id: "backgroundImage",
  beforeDraw: (chart: any, args: any, options: any) => {
    if (options.image) {
      const ctx = chart.ctx;
      const canvas = chart.canvas;
      const chartArea = chart.chartArea;

      // Check if image is loaded before drawing
      if (options.image.complete && options.image.naturalHeight !== 0) {
        // Save the context
        ctx.save();

        try {
          // Set global alpha for transparency if needed
          ctx.globalAlpha = options.opacity || 1.0;

          // Draw the image within the chart area
          ctx.drawImage(
            options.image,
            chartArea.left,
            chartArea.top,
            chartArea.right - chartArea.left,
            chartArea.bottom - chartArea.top
          );
        } catch (error) {
          console.error("Error drawing background image:", error);
        } finally {
          // Always restore the context
          ctx.restore();
        }
      } else {
        // Image not loaded yet, schedule a redraw
        console.log("Background image not ready, scheduling redraw");
        setTimeout(() => {
          if (chart && chart.update) {
            chart.update("none");
          }
        }, 100);
      }
    }
  },
};

// Register the plugin
Chart.register(backgroundImagePlugin);
