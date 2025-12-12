import { Editor } from "@monaco-editor/react";
import { useEffect, useMemo, useRef, useState } from "react";
import { CCP4Table, parseXML, Plot } from "./ChartLib";

import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  BarElement,
  BarController,
  ScatterController,
  LineElement,
  Title,
  Tooltip,
  Legend,
  ChartOptions,
  ChartData,
} from "chart.js";
import { Scatter } from "react-chartjs-2";

import annotationPlugin from "chartjs-plugin-annotation";
import {
  Autocomplete,
  Box,
  Button,
  Card,
  CardContent,
  CardHeader,
  Dialog,
  DialogContent,
  DialogTitle,
  IconButton,
  Menu,
  MenuItem,
  TextField,
  Typography,
} from "@mui/material";
import {
  OpenInNew,
  Download,
  Edit as EditIcon,
  Close,
} from "@mui/icons-material";
import { useTheme } from "../../theme/theme-provider";
import * as XLSX from "xlsx";
import jsPDF from "jspdf";
import autoTable from "jspdf-autotable";

// Simple client-side download utility - no external dependencies needed
const downloadBlob = (blob: Blob, filename: string) => {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

// Decode HTML entities (e.g., &nbsp; -> space, &amp; -> &)
const decodeHTMLEntities = (text: string): string => {
  const textarea = document.createElement("textarea");
  textarea.innerHTML = text;
  return textarea.value;
};

// Simple CSV generation - handles arrays of objects
const jsonToCSV = (data: any[]): string => {
  if (!data || data.length === 0) return "";

  const headers = Object.keys(data[0]);
  const escapeCSV = (val: any): string => {
    if (val === null || val === undefined) return "";
    const str = String(val);
    // Escape quotes and wrap in quotes if contains comma, quote, or newline
    if (str.includes(",") || str.includes('"') || str.includes("\n")) {
      return `"${str.replace(/"/g, '""')}"`;
    }
    return str;
  };

  const headerRow = headers.map(escapeCSV).join(",");
  const dataRows = data.map(row =>
    headers.map(header => escapeCSV(row[header])).join(",")
  );

  return [headerRow, ...dataRows].join("\n");
};

// Register components for both Scatter and Bar charts
ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  BarElement,
  BarController,
  ScatterController,
  LineElement,
  Title,
  Tooltip,
  Legend,
  annotationPlugin
);

interface CCP4i2ApplicationOutputViewProps {
  output: Element;
  height?: string | number;
  jobId?: number;
  graphId?: string;
}

const colours = [
  "rgba(255, 99, 132, 1)",
  "rgba(54, 162, 235, 1)",
  "rgba(255, 206, 86, 1)",
  "rgba(75, 192, 192, 1)",
  "rgba(153, 102, 255, 1)",
  "rgba(255, 159, 64, 1)",
];
interface ChartArgs {
  options: ChartOptions;
  data: any;
}
export const CCP4i2ApplicationOutputView: React.FC<
  CCP4i2ApplicationOutputViewProps
> = ({ output, height = "400px", jobId, graphId }) => {
  const [table, setTable] = useState<CCP4Table | null>(null);
  const [selectedPlot, setSelectedPlot] = useState<Plot | null>(null);
  const [exportMenuAnchor, setExportMenuAnchor] = useState<null | HTMLElement>(null);
  const [editDialogOpen, setEditDialogOpen] = useState(false);
  const [editedOptions, setEditedOptions] = useState<ChartOptions | null>(null);
  const { mode } = useTheme();

  const chartRef = useRef<any>(null);

  useEffect(() => {
    const asyncEffect = async () => {
      const serializedOutput = new XMLSerializer().serializeToString(output);
      try {
        const parsedOutput = await parseXML(serializedOutput);
        if (
          parsedOutput &&
          parsedOutput.CCP4Table &&
          !Array.isArray(parsedOutput.CCP4Table)
        ) {
          const graph = new CCP4Table(parsedOutput.CCP4Table);
          if (graph) {
            setTable(graph);
            setSelectedPlot(
              Array.isArray(graph.plot)
                ? graph.plot[0]
                : graph.plot
                  ? graph.plot
                  : null
            );
          }
        }
      } catch (err) {
        console.log(err);
      }
    };
    asyncEffect();
  }, [output]);

  const allPlots = useMemo(() => {
    return table ? table.allPlots : [];
  }, [table]);

  // Export handlers
  const handleExportPNG = () => {
    if (chartRef.current) {
      const url = chartRef.current.toBase64Image();
      const link = document.createElement("a");
      link.download = `${selectedPlot?.title || "chart"}.png`;
      link.href = url;
      link.click();
    }
    setExportMenuAnchor(null);
  };

  const handleExportCSV = () => {
    if (!table || !selectedPlot) return;

    const parsedBlocks = table.parsedDataBlocks;
    if (!parsedBlocks) return;

    const headers = table.allHeaders;
    if (!headers || headers.length === 0) return;

    // Decode HTML entities in headers
    const decodedHeaders = headers.map(h => decodeHTMLEntities(h));

    // Convert parsed data blocks to array of objects
    const dataArray: any[] = [];
    Object.values(parsedBlocks).forEach((block: any) => {
      if (Array.isArray(block)) {
        block.forEach((row: any[]) => {
          if (Array.isArray(row) && row.length > 0) {
            const rowObj: any = {};
            decodedHeaders.forEach((header, i) => {
              rowObj[header] = row[i] ?? "";
            });
            dataArray.push(rowObj);
          }
        });
      }
    });

    if (dataArray.length === 0) return;

    const csv = jsonToCSV(dataArray);
    const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
    downloadBlob(blob, `${selectedPlot.title || "data"}.csv`);
    setExportMenuAnchor(null);
  };

  const handleExportExcel = () => {
    if (!table || !selectedPlot) return;

    const parsedBlocks = table.parsedDataBlocks;
    if (!parsedBlocks) return;

    const headers = table.allHeaders;
    if (!headers || headers.length === 0) return;

    // Decode HTML entities in headers
    const decodedHeaders = headers.map(h => decodeHTMLEntities(h));

    // Convert parsed data blocks to array of objects
    const dataArray: any[] = [];
    Object.values(parsedBlocks).forEach((block: any) => {
      if (Array.isArray(block)) {
        block.forEach((row: any[]) => {
          if (Array.isArray(row) && row.length > 0) {
            const rowObj: any = {};
            decodedHeaders.forEach((header, i) => {
              rowObj[header] = row[i] ?? "";
            });
            dataArray.push(rowObj);
          }
        });
      }
    });

    if (dataArray.length === 0) return;

    const worksheet = XLSX.utils.json_to_sheet(dataArray);
    const workbook = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(workbook, worksheet, "Data");

    const excelBuffer = XLSX.write(workbook, { bookType: "xlsx", type: "array" });
    const blob = new Blob([excelBuffer], { type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet" });
    downloadBlob(blob, `${selectedPlot.title || "data"}.xlsx`);
    setExportMenuAnchor(null);
  };

  const handleExportPDF = () => {
    if (!table || !selectedPlot || !chartRef.current) return;

    const doc = new jsPDF();

    // Add title
    doc.setFontSize(16);
    doc.text(selectedPlot.title || "Chart", 14, 20);

    // Add chart image with maintained aspect ratio
    const canvas = chartRef.current.canvas;
    if (canvas) {
      const imgData = chartRef.current.toBase64Image();

      // Calculate dimensions to fit within PDF while maintaining aspect ratio
      const canvasAspectRatio = canvas.width / canvas.height;
      const maxWidth = 180; // Max width in PDF units
      const maxHeight = 100; // Max height in PDF units

      let imgWidth = maxWidth;
      let imgHeight = maxWidth / canvasAspectRatio;

      // If height exceeds max, scale down based on height instead
      if (imgHeight > maxHeight) {
        imgHeight = maxHeight;
        imgWidth = maxHeight * canvasAspectRatio;
      }

      // Center the image horizontally
      const xOffset = 14 + (maxWidth - imgWidth) / 2;

      doc.addImage(imgData, "PNG", xOffset, 30, imgWidth, imgHeight);
    }

    // Add data table
    const parsedBlocks = table.parsedDataBlocks;
    const headers = table.allHeaders;

    if (parsedBlocks && headers && headers.length > 0) {
      // Decode HTML entities in headers
      const decodedHeaders = headers.map(h => decodeHTMLEntities(h));

      // Convert parsed data blocks to array of arrays for the table
      const rows: any[][] = [];
      Object.values(parsedBlocks).forEach((block: any) => {
        if (Array.isArray(block)) {
          block.forEach((row: any[]) => {
            if (Array.isArray(row) && row.length > 0) {
              rows.push(row);
            }
          });
        }
      });

      if (rows.length > 0) {
        autoTable(doc, {
          head: [decodedHeaders],
          body: rows,
          startY: 140,
          theme: "grid",
          styles: { fontSize: 8 },
        });
      }
    }

    doc.save(`${selectedPlot.title || "chart"}.pdf`);
    setExportMenuAnchor(null);
  };

  const handleExportSVG = () => {
    if (!chartRef.current) return;

    // Get the canvas element
    const canvas = chartRef.current.canvas;
    if (!canvas) return;

    // Create SVG with embedded canvas image
    const svgString = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${canvas.width}" height="${canvas.height}">
  <image href="${canvas.toDataURL()}" width="${canvas.width}" height="${canvas.height}"/>
</svg>`;

    const blob = new Blob([svgString], { type: "image/svg+xml" });
    downloadBlob(blob, `${selectedPlot?.title || "chart"}.svg`);
    setExportMenuAnchor(null);
  };

  const chartArgs: ChartArgs | null = useMemo(() => {
    if (selectedPlot && table) {
      const options: ChartOptions | null = table.plotOptions(selectedPlot);
      const data: ChartData<"bar" | "scatter"> | null =
        table.plotData(selectedPlot);
      if (!options || !data) {
        return null;
      }
      return { data, options };
    }
    return null;
  }, [selectedPlot, table]);

  return (
    <>
      <Card variant="outlined" sx={{ height, display: 'flex', flexDirection: 'column' }}>
        <CardHeader
          title={
            allPlots &&
            selectedPlot &&
            allPlots.length > 0 && (
              <Autocomplete
                sx={{
                  px: 0,
                  py: 0,
                }}
                options={allPlots}
                getOptionKey={(option) => allPlots.indexOf(option)}
                getOptionLabel={(option) => option.title || "Unnamed plot"}
                onChange={(event, newValue) => {
                  setSelectedPlot(newValue);
                }}
                disabled={allPlots.length === 1}
                value={selectedPlot}
                renderInput={(params) => (
                  <TextField {...params} size="small" label="Plot" />
                )}
              />
            )
          }
          action={
            <>
              {jobId && graphId && (
                <Button
                  onClick={() => {
                    window.open(
                      `/graph-viewer/${jobId}/${graphId}`,
                      "_blank",
                      "width=1000,height=700"
                    );
                  }}
                  title="Open in standalone viewer"
                >
                  <OpenInNew />
                </Button>
              )}
              <Button
                onClick={(e) => setExportMenuAnchor(e.currentTarget)}
                title="Export data or image"
              >
                <Download />
              </Button>
              <Button
                onClick={() => setEditDialogOpen(true)}
                title="Edit chart properties"
              >
                <EditIcon />
              </Button>
            </>
          }
        />
        <CardContent sx={{ flex: 1, display: 'flex', flexDirection: 'column', padding: 2, minHeight: 0 }}>
          {chartArgs && (
            <Box sx={{ flex: 1, position: 'relative', minHeight: 0, width: '100%', height: '100%' }}>
              <Scatter
                ref={chartRef}
                key={
                  allPlots && selectedPlot
                    ? `plot-${allPlots.indexOf(selectedPlot)}`
                    : "scatter-plot"
                }
                options={editedOptions || chartArgs.options as any}
                data={chartArgs.data}
              />
            </Box>
          )}
        </CardContent>
      </Card>

      {/* Export Menu */}
      <Menu
        anchorEl={exportMenuAnchor}
        open={Boolean(exportMenuAnchor)}
        onClose={() => setExportMenuAnchor(null)}
      >
        <MenuItem onClick={handleExportPNG}>Export Chart as PNG</MenuItem>
        <MenuItem onClick={handleExportSVG}>Export Chart as SVG (raster embedded)</MenuItem>
        <MenuItem onClick={handleExportPDF}>Export as PDF (Chart + Data Table)</MenuItem>
        <MenuItem onClick={handleExportCSV}>Export Data as CSV</MenuItem>
        <MenuItem onClick={handleExportExcel}>Export Data as Excel</MenuItem>
      </Menu>

      {/* Edit Dialog */}
      <Dialog
        fullWidth
        maxWidth="md"
        open={editDialogOpen}
        onClose={() => setEditDialogOpen(false)}
      >
        <DialogTitle>
          Edit Chart Properties
          <IconButton
            onClick={() => setEditDialogOpen(false)}
            sx={{ position: 'absolute', right: 8, top: 8 }}
          >
            <Close />
          </IconButton>
        </DialogTitle>
        <DialogContent>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
            Edit the chart configuration JSON below. Changes will be applied when you close this dialog.
          </Typography>
          <Editor
            height="400px"
            defaultLanguage="json"
            value={JSON.stringify(editedOptions || chartArgs?.options || {}, null, 2)}
            onChange={(value) => {
              try {
                if (value) {
                  const parsed = JSON.parse(value);
                  setEditedOptions(parsed);
                }
              } catch (e) {
                // Invalid JSON, ignore
              }
            }}
            theme={mode === "dark" ? "vs-dark" : "light"}
            options={{
              minimap: { enabled: false },
              scrollBeyondLastLine: false,
            }}
          />
          <Box sx={{ mt: 2, display: 'flex', gap: 1 }}>
            <Button
              variant="outlined"
              onClick={() => {
                setEditedOptions(null);
                setEditDialogOpen(false);
              }}
            >
              Reset to Default
            </Button>
            <Button
              variant="contained"
              onClick={() => setEditDialogOpen(false)}
            >
              Apply Changes
            </Button>
          </Box>
        </DialogContent>
      </Dialog>
    </>
  );
};
