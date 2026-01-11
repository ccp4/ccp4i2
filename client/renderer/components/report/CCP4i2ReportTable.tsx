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

export const CCP4i2ReportTable: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const theme = useTheme();
  const [data, setData] = useState<any[]>([]);
  const [columns, setColumns] = useState<TableColumn[]>([]);

  useEffect(() => {
    let newData: any[] = [];
    let newColumns: TableColumn[] = [];

    $(props.item)
      .find("tr")
      .each((iRow, tableRow) => {
        var dataItem: any = { key: iRow };
        if ($(props.item).attr("transpose") === "True") {
          $(tableRow)
            .find("th")
            .each((iColumn, tableData) => {
              dataItem["col_" + iColumn] = $(tableData).html();
              dataItem["col_" + iColumn + "_text"] = $(tableData).text();
            });
        }
        const thCount = Object.keys(dataItem).filter(
          (k) => !k.endsWith("_text") && k !== "key"
        ).length;
        $(tableRow)
          .find("td")
          .each((iColumn, tableData) => {
            dataItem["col_" + (thCount + iColumn - 1)] = $(tableData).html();
            dataItem["col_" + (thCount + iColumn - 1) + "_text"] =
              $(tableData).text();
          });
        if ($(tableRow).find("td").length > 0) {
          newData.push(dataItem);
        }
        if ($(props.item).attr("transpose") === "False") {
          $(tableRow)
            .find("th")
            .each((iColumn, tableData) => {
              while (iColumn >= newColumns.length) {
                newColumns.push({ key: "", title: "", titleText: "" });
              }
              newColumns[iColumn] = {
                key: "col_" + iColumn,
                title: tableData.innerHTML,
                titleText: $(tableData).text(),
              };
            });
        }
      });

    if (newColumns.length === 0 && newData.length > 0) {
      const colCount = Object.keys(newData[0]).filter(
        (k) => !k.endsWith("_text") && k !== "key"
      ).length;
      for (var iCol = 0; iCol < colCount; iCol++) {
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
          .map((col) => escapeCSV(row[col.key + "_text"] || row[col.key] || ""))
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

  return (
    <Box sx={{ mx: 2, my: 0.5 }}>
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
              <TableCell sx={{ width: 48, p: 0.5 }}>
                <Tooltip title="Download as CSV">
                  <IconButton
                    size="small"
                    onClick={downloadCSV}
                    sx={{
                      color: "primary.main",
                      "&:hover": {
                        bgcolor: alpha(theme.palette.primary.main, 0.12),
                      },
                    }}
                  >
                    <Download fontSize="small" />
                  </IconButton>
                </Tooltip>
              </TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {data.map((row) => (
              <TableRow
                key={row.key}
                sx={{
                  "&:hover": {
                    bgcolor: alpha(theme.palette.primary.main, 0.04),
                  },
                  "&:last-child td": { borderBottom: 0 },
                }}
              >
                {columns.map((col) => (
                  <TableCell
                    key={col.key}
                    sx={{ py: 1, px: 2, color: "text.secondary" }}
                    dangerouslySetInnerHTML={{ __html: row[col.key] || "" }}
                  />
                ))}
                <TableCell sx={{ width: 48 }} />
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>
    </Box>
  );
};
