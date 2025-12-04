import React from "react";
import { SxProps, Table, TableBody, TableCell, TableRow } from "@mui/material";

interface SimpleObjectTableProps {
  object: any | null;
  sx?: SxProps;
}

// Utility function to check if a value is a simple type
const isSimpleValue = (value: any): boolean => {
  return (
    value === null ||
    value === undefined ||
    typeof value === "string" ||
    typeof value === "number" ||
    typeof value === "boolean"
  );
};

// Utility function to format value for display
const formatValue = (value: any): string => {
  if (value === null) return "null";
  if (value === undefined) return "undefined";
  if (isSimpleValue(value)) return String(value);
  return JSON.stringify(value, null, 2);
};

export const SimpleObjectTable: React.FC<SimpleObjectTableProps> = ({
  object,
  sx,
}) => {
  if (!object) {
    return null;
  }

  return (
    <Table sx={sx || { mb: 2 }} size="small">
      <TableBody>
        {Object.keys(object).map((key: string) => (
          <TableRow key={key}>
            <TableCell
              variant="head"
              sx={{ textAlign: "left", fontWeight: "medium" }}
            >
              {key}
            </TableCell>
            <TableCell
              variant="body"
              sx={{
                textAlign: "left", // Changed from center for better readability
                fontFamily: isSimpleValue(object[key])
                  ? "inherit"
                  : "monospace",
                whiteSpace: isSimpleValue(object[key]) ? "normal" : "pre-wrap",
                wordBreak: "break-word",
              }}
            >
              {formatValue(object[key])}
            </TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
};
