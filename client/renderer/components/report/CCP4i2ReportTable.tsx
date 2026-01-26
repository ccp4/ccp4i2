import { useEffect, useState, useCallback } from "react";
import {
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
  Box,
  IconButton,
  Tooltip,
} from "@mui/material";
import { Download } from "@mui/icons-material";
import { useTheme, alpha } from "@mui/material/styles";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import $ from "jquery";

interface TableColumn {
  key: string;
  title: string;
  titleText: string; // Plain text version for CSV export
}

interface CellData {
  html: string;
  text: string;
  isHeader: boolean;
}

export const CCP4i2ReportTable: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const theme = useTheme();
  const [data, setData] = useState<Record<string, CellData>[]>([]);
  const [columns, setColumns] = useState<TableColumn[]>([]);

  useEffect(() => {
    let newData: Record<string, CellData>[] = [];
    let newColumns: TableColumn[] = [];
    const transposeAttr = $(props.item).attr("transpose");
    const isHeaderRow = transposeAttr === "False";

    $(props.item)
      .find("tr")
      .each((_iRow, tableRow) => {
        const dataItem: Record<string, CellData> = {};
        let colIndex = 0;

        // Process all cells (th and td) in order
        $(tableRow)
          .find("th, td")
          .each((_i, cell) => {
            const $cell = $(cell);
            const isHeader = cell.tagName.toLowerCase() === "th";

            // For transpose="False" mode, th elements in header rows define columns
            if (isHeaderRow && isHeader) {
              while (colIndex >= newColumns.length) {
                newColumns.push({ key: "", title: "", titleText: "" });
              }
              newColumns[colIndex] = {
                key: "col_" + colIndex,
                title: cell.innerHTML,
                titleText: $cell.text(),
              };
            }

            dataItem["col_" + colIndex] = {
              html: $cell.html() || "",
              text: $cell.text(),
              isHeader: isHeader,
            };
            colIndex++;
          });

        // Only add rows that have data cells (td)
        if ($(tableRow).find("td").length > 0) {
          newData.push(dataItem);
        }
      });

    // If no explicit columns were defined, create them from the data
    if (newColumns.length === 0 && newData.length > 0) {
      const colCount = Object.keys(newData[0]).length;
      for (let iCol = 0; iCol < colCount; iCol++) {
        newColumns.push({
          key: "col_" + iCol,
          title: "",
          titleText: "",
        });
      }
    }

    setData(newData);
    setColumns(newColumns);
  }, [props.job, props.item]);

  const downloadCSV = useCallback(() => {
    if (columns.length === 0 || data.length === 0) return;

    // Build CSV content
    const escapeCSV = (value: string) => {
      if (value.includes(",") || value.includes('"') || value.includes("\n")) {
        return `"${value.replace(/"/g, '""')}"`;
      }
      return value;
    };

    const headerRow = columns.map((col) => escapeCSV(col.titleText)).join(",");
    const dataRows = data
      .map((row) =>
        columns
          .map((col) => {
            const cell = row[col.key];
            return escapeCSV(cell?.text || "");
          })
          .join(",")
      )
      .join("\n");

    const csvContent = headerRow + "\n" + dataRows;

    // Create and trigger download
    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.setAttribute(
      "download",
      `${props.job?.task_name || "table"}_report.csv`
    );
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  }, [columns, data, props.job]);

  if (columns.length === 0 || data.length === 0) {
    return null;
  }

  // Check if we have any actual column titles
  const hasColumnTitles = columns.some((col) => col.title.trim() !== "");

  // Theme-aware colors for zebra striping that work in both light and dark modes
  const isDark = theme.palette.mode === "dark";
  const zebraColor = isDark
    ? alpha(theme.palette.common.white, 0.03)
    : alpha(theme.palette.common.black, 0.02);
  const hoverColor = isDark
    ? alpha(theme.palette.primary.main, 0.15)
    : alpha(theme.palette.primary.main, 0.08);
  const rowHeaderBg = isDark
    ? alpha(theme.palette.common.white, 0.06)
    : alpha(theme.palette.common.black, 0.04);
  const rowHeaderHoverBg = isDark
    ? alpha(theme.palette.primary.main, 0.2)
    : alpha(theme.palette.primary.main, 0.12);

  return (
    <Box sx={{ mx: 2, my: 0.5 }}>
      {/* Toolbar with download button */}
      <Box
        sx={{
          display: "flex",
          justifyContent: "flex-end",
          mb: 0.5,
        }}
      >
        <Tooltip title="Download as CSV">
          <IconButton
            size="small"
            onClick={downloadCSV}
            sx={{
              color: "text.secondary",
              "&:hover": {
                color: "primary.main",
                bgcolor: alpha(theme.palette.primary.main, 0.08),
              },
            }}
          >
            <Download fontSize="small" />
          </IconButton>
        </Tooltip>
      </Box>

      <TableContainer
        component={Paper}
        elevation={0}
        sx={{
          borderRadius: 1,
          border: 1,
          borderColor: "divider",
          overflow: "hidden",
        }}
      >
        <Table size="small">
          {hasColumnTitles && (
            <TableHead>
              <TableRow
                sx={{
                  bgcolor: alpha(theme.palette.primary.main, 0.08),
                  "& .MuiTableCell-head": {
                    fontWeight: 600,
                    color: "text.primary",
                    borderBottom: 2,
                    borderColor: "primary.main",
                    py: 1.5,
                    px: 2,
                  },
                }}
              >
                {columns.map((col) => (
                  <TableCell
                    key={col.key}
                    dangerouslySetInnerHTML={{ __html: col.title }}
                  />
                ))}
              </TableRow>
            </TableHead>
          )}
          <TableBody>
            {data.map((row, rowIndex) => {
              const isOddRow = rowIndex % 2 === 1;
              return (
                <TableRow
                  key={rowIndex}
                  sx={{
                    bgcolor: isOddRow ? zebraColor : "transparent",
                    transition: "background-color 0.15s ease-in-out",
                    "&:hover": {
                      bgcolor: hoverColor,
                    },
                    "&:last-child td, &:last-child th": { borderBottom: 0 },
                  }}
                >
                  {columns.map((col) => {
                    const cell = row[col.key];
                    const isHeader = cell?.isHeader;
                    return (
                      <TableCell
                        key={col.key}
                        component={isHeader ? "th" : "td"}
                        scope={isHeader ? "row" : undefined}
                        sx={{
                          py: 1,
                          px: 2,
                          color: isHeader ? "text.primary" : "text.secondary",
                          fontWeight: isHeader ? 600 : 400,
                          bgcolor: isHeader ? rowHeaderBg : "transparent",
                          transition: "background-color 0.15s ease-in-out",
                          "tr:hover &": isHeader
                            ? { bgcolor: rowHeaderHoverBg }
                            : {},
                        }}
                        dangerouslySetInnerHTML={{ __html: cell?.html || "" }}
                      />
                    );
                  })}
                </TableRow>
              );
            })}
          </TableBody>
        </Table>
      </TableContainer>
    </Box>
  );
};
