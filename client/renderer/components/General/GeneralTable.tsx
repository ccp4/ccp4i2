/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import { useEffect, useState, ReactNode, useMemo, useCallback } from "react";
import {
  Box,
  Button,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Typography,
  MenuItem,
  InputLabel,
  Select,
  Toolbar,
  AppBar,
  Autocomplete,
  TextField,
} from "@mui/material";
import {
  ArrowBack,
  ArrowForward,
  Delete,
  Download,
  Edit,
  Visibility,
} from "@mui/icons-material";
import { ElaborateSearch } from "./SearchObjects";
import {
  SxProps,
  TableProps,
  TableRowProps,
  TableCellProps,
} from "@mui/material";
import { useTheme } from "@mui/material/styles";
import ExcelJS from "exceljs";

export interface GeneralTableColumn {
  title: string;
  key: string;
  dataIndex: string;
  searchable?: boolean | null;
  render?: (item?: any, row?: any, index?: number) => ReactNode;
}

export interface GeneralTableProps {
  dataSource: any[] | null;
  columns?: GeneralTableColumn[];
  pagination?: boolean | undefined;
  extra?: ReactNode[] | undefined;
  sx?: SxProps | undefined;
  tableProps?: TableProps | undefined;
  tableRowProps?: TableRowProps | undefined;
  tableHeadProps?: TableRowProps | undefined;
  tableCellProps?: TableCellProps | undefined;
  title?: ReactNode | string | null;
  locked?: boolean | undefined;
  onAdd?: Function | undefined;
  onEdit?: Function | undefined;
  onDelete?: Function | undefined;
  modelName?: string | undefined | null;
  onRow?: (item: any, iItem?: number) => void | undefined;
  rowClassName?: (item: any, iItem?: number) => string;
  dataIndex?: string | Function | undefined;
  selectColumns?: boolean;
  size?: any;
  showToolbar?: boolean;
}

export const GeneralTable = (props: GeneralTableProps) => {
  const theme = useTheme();
  const [filteredItems, setFilteredItems] = useState<any[] | null>(null);
  const [paginatedItems, setPaginatedItems] = useState<any[] | null>(null);
  const [paginationLength, setPaginationLength] = useState(30);
  const [filterPredicate, setFilterPredicate] = useState<string | null>(null);
  const [paginationStart, setPaginationStart] = useState<number | null>(null);
  const [selectedColumnKeys, setSelectedColumnKeys] = useState<string[] | null>(
    null
  );

  const dataAsHTMLTable = useCallback(() => {
    const tableAsHTML = `\
<table>
                <thead>
                    <tr>
                        ${props.columns
                          ?.map((column) => {
                            return `<th>${column.title}</th>`;
                          })
                          .join("\n")}
                    </tr>
                </thead>
                <tbody>
                    ${props.dataSource
                      ?.map((row, iRow) => {
                        let rowOfCells = props.columns
                          ?.map((column) => {
                            let finalValue:
                              | string
                              | number
                              | String
                              | Number
                              | null = null;
                            let renderedValue = null;
                            if (Object.keys(column).includes("render")) {
                              renderedValue = (column.render as Function)(
                                row[column.dataIndex],
                                row,
                                iRow
                              );
                              if (
                                renderedValue === null ||
                                typeof renderedValue === "string"
                              ) {
                                finalValue = renderedValue;
                              }
                            }
                            if (finalValue === null) {
                              const maybeStringOrNumber = row[column.dataIndex];
                              finalValue =
                                typeof maybeStringOrNumber === "string" ||
                                typeof maybeStringOrNumber === "number" ||
                                maybeStringOrNumber instanceof String ||
                                maybeStringOrNumber instanceof Number
                                  ? maybeStringOrNumber
                                  : null;
                            }
                            if (finalValue === null) {
                              try {
                                const maybeStringOrNumber =
                                  row[column.dataIndex].props.value;
                                finalValue =
                                  typeof maybeStringOrNumber === "string" ||
                                  typeof maybeStringOrNumber === "number" ||
                                  maybeStringOrNumber instanceof String ||
                                  maybeStringOrNumber instanceof Number
                                    ? maybeStringOrNumber
                                    : null;
                              } catch (err) {
                                finalValue = null;
                              }
                            }
                            if (finalValue === null) {
                              try {
                                finalValue = JSON.stringify(
                                  row[column.dataIndex]
                                );
                              } catch (err) {
                                finalValue = "";
                              }
                            }
                            return `<td>${finalValue}</td>`;
                          })
                          .join("");
                        return `<tr>${rowOfCells}</tr>`;
                      })
                      .join("\n")}
                </tbody >
            </table>
`;
    return tableAsHTML;
  }, [props.dataSource, props.columns]);

  useEffect(() => {
    if (props.columns && props.columns.length > 0) {
      setSelectedColumnKeys(props.columns.map((item) => item.key));
    }
  }, []);

  const selectedColumns = useMemo(() => {
    if (!selectedColumnKeys) return props.columns;
    return (props.columns as GeneralTableColumn[]).filter((item) =>
      selectedColumnKeys.includes(item.key)
    );
  }, [props.columns, selectedColumnKeys]);

  const allColumnKeys = useMemo(() => {
    if (!props.columns) return [];
    return props.columns.map((item) => item.key);
  }, [props.columns]);

  useEffect(() => {
    if (props.dataSource) {
      if (filterPredicate && filterPredicate.length > 0) {
        const newFilteredItems: any[] = [];
        const lcFilterPredicate = filterPredicate.toLowerCase();
        props.dataSource.forEach((djangoItem, iDjangoItem) => {
          for (let column of props.columns as GeneralTableColumn[]) {
            if (
              Object.keys(column).includes("searchable") &&
              column.searchable
            ) {
              let finalValue: string | number | String | Number | null = null;
              let renderedValue = null;
              if (Object.keys(column).includes("render")) {
                renderedValue = (column.render as Function)(
                  djangoItem[column.dataIndex],
                  djangoItem,
                  iDjangoItem
                );
                if (
                  renderedValue === null ||
                  typeof renderedValue === "string"
                ) {
                  finalValue = renderedValue;
                }
              }
              if (finalValue === null) {
                finalValue =
                  typeof djangoItem[column.dataIndex] === "string" ||
                  typeof djangoItem[column.dataIndex] === "number" ||
                  djangoItem[column.dataIndex] instanceof String ||
                  djangoItem[column.dataIndex] instanceof Number
                    ? djangoItem[column.dataIndex]
                    : null;
              }
              if (finalValue === null) {
                try {
                  finalValue = JSON.stringify(djangoItem[column.dataIndex]);
                } catch (err) {
                  finalValue = "";
                }
              }
              if (`${finalValue}`.toLowerCase().includes(lcFilterPredicate)) {
                newFilteredItems.push(djangoItem);
                break;
              }
            }
          }
        });
        setFilteredItems(newFilteredItems);
      } else {
        setFilteredItems(props.dataSource);
      }
    }
  }, [props.dataSource, filterPredicate]);

  useEffect(() => {
    setPaginationStart(0);
  }, [filteredItems]);

  useEffect(() => {
    if (filteredItems) {
      if (
        props.pagination &&
        paginationLength > 0 &&
        paginationStart !== null
      ) {
        setPaginatedItems(
          filteredItems.slice(
            paginationStart,
            paginationStart + paginationLength
          )
        );
      } else {
        setPaginatedItems(filteredItems);
      }
    }
  }, [filteredItems, paginationLength, paginationStart]);

  return (
    <Paper sx={props.sx}>
      {props.selectColumns && selectedColumnKeys && (
        <Autocomplete
          multiple={true}
          value={selectedColumnKeys}
          options={allColumnKeys}
          onChange={(changeEvent, newValue) => {
            setSelectedColumnKeys(newValue);
          }}
          renderInput={(params) => (
            <TextField {...params} label="Columns to include" />
          )}
        />
      )}
      {props.showToolbar && (
        <AppBar position="static">
          <Toolbar>
            {props.pagination &&
              filteredItems &&
              paginationLength < filteredItems.length && (
                <InputLabel
                  key="Showing"
                  sx={{
                    color: theme.palette.common.black,
                    marginLeft: "0.5rem",
                    display: { sm: "none", lg: "block" },
                  }}
                >
                  {props.modelName}:Showing {paginationStart} to{" "}
                  {(paginationStart as number) + paginationLength} of{" "}
                </InputLabel>
              )}
            {filterPredicate && filterPredicate.length > 0 && (
              <InputLabel
                key="from"
                sx={{
                  color: theme.palette.common.black,
                  marginLeft: "0.5rem",
                  display: { sm: "none", lg: "block" },
                }}
              >
                filtered from {props.dataSource && props.dataSource.length}{" "}
                total
              </InputLabel>
            )}
            <Typography
              key="spacer"
              variant="h6"
              noWrap
              component="div"
              sx={{ flexGrow: 1, display: { xs: "none", sm: "block" } }}
            >
              {props.title && <span>{props.title}</span>}
            </Typography>
            {props.onAdd && (
              <Button
                variant="contained"
                color="primary"
                sx={{ width: "15rem", marginX: "0.5rem" }}
                disabled={props.locked}
                onClick={() => {
                  (props.onAdd as Function)();
                }}
              >
                {props.locked
                  ? `Unlock to add ${props.modelName} `
                  : `Add new ${props.modelName} `}
              </Button>
            )}
            <Typography
              key="secondSpacer"
              variant="h6"
              noWrap
              component="div"
              sx={{ flexGrow: 1, display: { xs: "none", sm: "block" } }}
            />
            {props.extra &&
              (props.extra as ReactNode[]).map((extraItem) => extraItem)}
            {false && (
              <Button
                key="AsHTML"
                variant="contained"
                size="small"
                sx={{
                  paddingX: "0px",
                  minWidth: "2rem",
                  maxWidth: "4rem",
                  marginRight: "0.5rem",
                }}
                disabled={!filteredItems || filteredItems?.length == 0}
                onClick={(ev) => {
                  const htmlText = dataAsHTMLTable();
                  // Create blob link to download
                  const url = URL.createObjectURL(new Blob([htmlText]));
                  const link = document.createElement("a");
                  link.href = url;
                  link.setAttribute("download", "DataSeries.html");

                  // Append to html link element page
                  document.body.appendChild(link);

                  // Start download
                  link.click();

                  // Clean up and remove the link
                  link.parentNode?.removeChild(link);
                }}
              >
                <Download />
                HTML
              </Button>
            )}
            <Button
              key="xlsx"
              variant="contained"
              size="small"
              sx={{
                paddingX: "0px",
                minWidth: "2rem",
                maxWidth: "4rem",
                marginRight: "0.5rem",
              }}
              disabled={!filteredItems || filteredItems.length == 0}
              onClick={async (ev) => {
                const htmlText = dataAsHTMLTable();

                // Parse HTML table and build ExcelJS workbook
                const parser = new DOMParser();
                const doc = parser.parseFromString(htmlText, "text/html");
                const table = doc.querySelector("table");
                if (!table) return;

                const wb = new ExcelJS.Workbook();
                const ws = wb.addWorksheet("Sheet1");

                table.querySelectorAll("tr").forEach((tr) => {
                  const rowData: (string | null)[] = [];
                  tr.querySelectorAll("th, td").forEach((cell) => {
                    rowData.push(cell.textContent);
                  });
                  ws.addRow(rowData);
                });

                const buffer = await wb.xlsx.writeBuffer();
                const blob = new Blob([buffer], {
                  type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                });
                const url = URL.createObjectURL(blob);
                const link = document.createElement("a");
                link.href = url;
                link.setAttribute("download", "Download.xlsx");
                document.body.appendChild(link);
                link.click();
                link.parentNode?.removeChild(link);
                URL.revokeObjectURL(url);
              }}
            >
              <Download />
              XLSX
            </Button>
            <Button
              key="csv"
              variant="contained"
              size="small"
              sx={{
                paddingX: "0px",
                minWidth: "2rem",
                maxWidth: "4rem",
                marginRight: "0.5rem",
              }}
              disabled={!filteredItems || filteredItems.length == 0}
              onClick={async (ev) => {
                const htmlText = dataAsHTMLTable();

                // Parse HTML table and build ExcelJS workbook for CSV export
                const parser = new DOMParser();
                const doc = parser.parseFromString(htmlText, "text/html");
                const table = doc.querySelector("table");
                if (!table) return;

                const wb = new ExcelJS.Workbook();
                const ws = wb.addWorksheet("Sheet1");

                table.querySelectorAll("tr").forEach((tr) => {
                  const rowData: (string | null)[] = [];
                  tr.querySelectorAll("th, td").forEach((cell) => {
                    rowData.push(cell.textContent);
                  });
                  ws.addRow(rowData);
                });

                const buffer = await wb.csv.writeBuffer();
                const blob = new Blob([buffer], { type: "text/csv;charset=utf-8;" });
                const url = URL.createObjectURL(blob);
                const link = document.createElement("a");
                link.href = url;
                link.setAttribute("download", "Download.csv");
                document.body.appendChild(link);
                link.click();
                link.parentNode?.removeChild(link);
                URL.revokeObjectURL(url);
              }}
            >
              <Download />
              CSV
            </Button>
            {props.pagination && (
              <>
                <InputLabel
                  sx={{
                    color: theme.palette.common.black,
                    display: { sm: "none", lg: "block" },
                  }}
                >
                  Show
                </InputLabel>
                <Select
                  size="small"
                  sx={{ minWidth: "10rem", color: theme.palette.common.black }}
                  value={paginationLength}
                  onChange={(e) => {
                    setPaginationLength(e.target.value as number);
                  }}
                >
                  <MenuItem key={0} sx={{ minWidth: "10rem" }} value={0}>
                    All
                  </MenuItem>
                  <MenuItem key={10} sx={{ minWidth: "10rem" }} value={10}>
                    10
                  </MenuItem>
                  <MenuItem key={30} sx={{ minWidth: "10rem" }} value={30}>
                    30
                  </MenuItem>
                  <MenuItem key={100} sx={{ minWidth: "10rem" }} value={100}>
                    100
                  </MenuItem>
                  <MenuItem key={300} sx={{ minWidth: "10rem" }} value={300}>
                    300
                  </MenuItem>
                  <MenuItem key={1000} sx={{ minWidth: "10rem" }} value={1000}>
                    1000
                  </MenuItem>
                </Select>
                <InputLabel
                  sx={{
                    color: theme.palette.common.black,
                    marginLeft: "0.5rem",
                    display: { sm: "none", lg: "block" },
                  }}
                >
                  {" "}
                  of {filteredItems && filteredItems.length}{" "}
                </InputLabel>
                <Box
                  sx={{
                    display: { sm: "none", lg: "block", marginLeft: "0.5rem" },
                  }}
                >
                  <Button
                    key="back"
                    variant="contained"
                    size="small"
                    sx={{ paddingX: "0px", minWidth: "2rem", maxWidth: "4rem" }}
                    disabled={
                      !filteredItems ||
                      (paginationStart as number) + paginationLength >=
                        filteredItems.length
                    }
                    onClick={() => {
                      setPaginationStart(
                        (paginationStart as number) + paginationLength
                      );
                    }}
                  >
                    <ArrowBack
                      sx={{
                        color: theme.palette.common.black,
                        paddingX: "1px",
                      }}
                    />
                  </Button>
                  <Button
                    key="forward"
                    variant="contained"
                    size="small"
                    sx={{ paddingX: "0px", minWidth: "2rem", maxWidth: "4rem" }}
                    disabled={
                      (paginationStart as number) - paginationLength < 0
                    }
                    onClick={() => {
                      setPaginationStart(
                        (paginationStart as number) - paginationLength
                      );
                    }}
                  >
                    <ArrowForward
                      sx={{
                        color: theme.palette.common.black,
                        paddingX: "1px",
                      }}
                    />
                  </Button>
                </Box>
              </>
            )}
            <ElaborateSearch
              searchValue={filterPredicate}
              setSearchValue={setFilterPredicate}
            />
          </Toolbar>
        </AppBar>
      )}
      <Table {...props.tableProps}>
        <TableHead>
          <TableRow {...(props.tableHeadProps || props.tableRowProps)}>
            {selectedColumns?.map((column, iColumn) => (
              <TableCell
                variant="head"
                padding="none"
                {...props.tableCellProps}
                size={props.size}
                key={`${iColumn}`}
              >
                {column.title}
              </TableCell>
            ))}
            {props.onEdit && (
              <TableCell
                padding="none"
                {...props.tableCellProps}
                key="View"
                sx={{ maxWidth: "10rem" }}
              >
                View/Edit
              </TableCell>
            )}
            {!props.locked && props.onDelete && (
              <TableCell
                padding="none"
                {...props.tableCellProps}
                key="Delete"
                sx={{ maxWidth: "10rem" }}
              >
                Delete
              </TableCell>
            )}
          </TableRow>
        </TableHead>
        <TableBody>
          {paginatedItems &&
            paginatedItems.length > 0 &&
            paginatedItems.map((djangoItem, iDjangoItem) => {
              const tableRowProps: any = { ...props.tableRowProps };
              if (props.rowClassName) {
                tableRowProps.className = props.rowClassName(
                  djangoItem,
                  iDjangoItem
                );
              }
              return (
                <TableRow
                  {...tableRowProps}
                  key={
                    typeof props.dataIndex === "function"
                      ? props.dataIndex(djangoItem)
                      : props.dataIndex
                        ? djangoItem[props.dataIndex]
                        : iDjangoItem
                  }
                  hover
                  onClick={(e: any) => {
                    if (props.onRow) {
                      props.onRow(djangoItem, iDjangoItem);
                    }
                  }}
                >
                  {selectedColumns?.map((column, iColumn) => (
                    <TableCell {...props.tableCellProps} key={iColumn}>
                      {Object.keys(column).includes("render")
                        ? (
                            column.render as (
                              item: any,
                              row?: any,
                              index?: number
                            ) => ReactNode
                          )(
                            djangoItem[column.dataIndex],
                            djangoItem,
                            iDjangoItem
                          )
                        : djangoItem[column.dataIndex]}
                    </TableCell>
                  ))}
                  {props.onEdit && (
                    <TableCell
                      padding="none"
                      {...props.tableCellProps}
                      key="View"
                      sx={{ maxWidth: "10rem" }}
                    >
                      {
                        <Button
                          onClick={(e) => {
                            if (props.onEdit) {
                              e.stopPropagation();
                              props.onEdit(djangoItem);
                            }
                          }}
                        >
                          {props.locked ? <Visibility /> : <Edit />}
                        </Button>
                      }
                    </TableCell>
                  )}
                  {!props.locked && props.onDelete && (
                    <TableCell
                      padding="none"
                      {...props.tableCellProps}
                      key="Delete"
                      sx={{ maxWidth: "10rem" }}
                    >
                      {
                        <Button
                          onClick={(e: any) => {
                            e.stopPropagation();
                            if (props.onDelete) {
                              //e.stopPropagation()
                              (props.onDelete as Function)(djangoItem);
                            }
                          }}
                        >
                          {<Delete />}
                        </Button>
                      }
                    </TableCell>
                  )}
                </TableRow>
              );
            })}
        </TableBody>
      </Table>
    </Paper>
  );
};
GeneralTable.defaultProps = {
  size: "small",
  pagination: true,
  selectColumns: false,
};
