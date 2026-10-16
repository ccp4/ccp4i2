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

interface HeaderCell {
  html: string;
  text: string;
  span: number;
}

export const CCP4i2ReportTable: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const theme = useTheme();
  const [data, setData] = useState<Record<string, CellData>[]>([]);
  const [columns, setColumns] = useState<TableColumn[]>([]);
  // Header rows as written, spans and all: a grouped heading ("RMSD
  // (bond/angle/chiral)" over Bond, Angle, Chiral) spans its columns, and
  // reading it as one column put every heading after it over the wrong data.
  const [headerRows, setHeaderRows] = useState<HeaderCell[][]>([]);

  useEffect(() => {
    let newData: Record<string, CellData>[] = [];
    let newColumns: TableColumn[] = [];
    const newHeaderRows: HeaderCell[][] = [];
    const transposeAttr = $(props.item).attr("transpose");
    const isHeaderRow = transposeAttr === "False";

    $(props.item)
      .find("tr")
      .each((_iRow, tableRow) => {
        // For transpose="False" mode, rows of th elements are the column headings
        if (isHeaderRow && $(tableRow).find("td").length === 0) {
          const headerRow: HeaderCell[] = [];
          $(tableRow)
            .find("th")
            .each((_i, cell) => {
              const span = parseInt($(cell).attr("colspan") || "1", 10);
              headerRow.push({
                html: cell.innerHTML,
                text: $(cell).text(),
                span: Number.isFinite(span) && span > 0 ? span : 1,
              });
            });
          if (headerRow.length > 0) newHeaderRows.push(headerRow);
          return;
        }

        const dataItem: Record<string, CellData> = {};
        let colIndex = 0;

        // Process all cells (th and td) in order
        $(tableRow)
          .find("th, td")
          .each((_i, cell) => {
            const $cell = $(cell);
            const isHeader = cell.tagName.toLowerCase() === "th";
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

    // One column per spanned position. Its title, for the CSV, joins every
    // heading above it: "RMSD (bond/angle/chiral) Bond".
    if (newHeaderRows.length > 0) {
      const nColumns = Math.max(
        ...newHeaderRows.map((row) => row.reduce((n, cell) => n + cell.span, 0))
      );
      const titles: string[][] = Array.from({ length: nColumns }, () => []);
      const htmls: string[] = Array(nColumns).fill("");
      newHeaderRows.forEach((row) => {
        let position = 0;
        row.forEach((cell) => {
          for (let k = 0; k < cell.span && position + k < nColumns; k++) {
            if (cell.text.trim()) titles[position + k].push(cell.text.trim());
            if (cell.html.trim()) htmls[position + k] = cell.html;
          }
          position += cell.span;
        });
      });
      newColumns = titles.map((parts, iCol) => ({
        key: "col_" + iCol,
        title: htmls[iCol],
        titleText: parts.join(" "),
      }));
    }

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
    setHeaderRows(newHeaderRows);
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
              {headerRows.map((headerRow, iHeaderRow) => (
                <TableRow
                  key={iHeaderRow}
                  sx={{
                    bgcolor: alpha(theme.palette.primary.main, 0.08),
                    "& .MuiTableCell-head": {
                      fontWeight: 600,
                      color: "text.primary",
                      // The rule goes under the last heading row only
                      borderBottom: iHeaderRow === headerRows.length - 1 ? 2 : 0,
                      borderColor: "primary.main",
                      whiteSpace: "pre-line",
                      py: iHeaderRow === headerRows.length - 1 ? 1.5 : 0.75,
                      px: 2,
                    },
                  }}
                >
                  {headerRow.map((cell, iCell) => (
                    <TableCell
                      key={iCell}
                      colSpan={cell.span > 1 ? cell.span : undefined}
                      align={cell.span > 1 ? "center" : undefined}
                      dangerouslySetInnerHTML={{ __html: cell.html }}
                    />
                  ))}
                </TableRow>
              ))}
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
