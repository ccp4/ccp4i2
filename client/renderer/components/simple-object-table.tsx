import React from "react";
import {
  Box,
  Paper,
  SxProps,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableRow,
} from "@mui/material";
import { useTheme, alpha } from "@mui/material/styles";

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
  const theme = useTheme();

  if (!object) {
    return null;
  }

  return (
    <Box sx={sx}>
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
          <TableBody>
            {Object.keys(object).map((key: string) => (
              <TableRow
                key={key}
                sx={{
                  "&:hover": {
                    bgcolor: alpha(theme.palette.primary.main, 0.04),
                  },
                  "&:last-child td": { borderBottom: 0 },
                }}
              >
                <TableCell
                  sx={{
                    fontWeight: 600,
                    color: "text.primary",
                    bgcolor: alpha(theme.palette.primary.main, 0.08),
                    borderRight: 1,
                    borderColor: "divider",
                    py: 1,
                    px: 2,
                    width: "30%",
                  }}
                >
                  {key}
                </TableCell>
                <TableCell
                  sx={{
                    color: "text.secondary",
                    py: 1,
                    px: 2,
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
      </TableContainer>
    </Box>
  );
};
