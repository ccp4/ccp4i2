/*
 * Copyright (C) 2025 Newcastle University
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
import React from "react";
import {
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
} from "@mui/material";

export interface CsvTableProps {
  csvText: string;
}

function parseCSV(csvText: string): string[][] {
  return csvText
    .trim()
    .split(/\r?\n/)
    .map((line) =>
      line.split(",").map((cell) => cell.trim().replace(/^"|"$/g, ""))
    );
}

export const CsvTable: React.FC<CsvTableProps> = ({ csvText }) => {
  const rows = parseCSV(csvText);
  if (rows.length === 0) return null;
  const header = rows[0];
  const dataRows = rows.slice(1);

  return (
    <TableContainer
      component={Paper}
      sx={{ width: "100%", maxHeight: 500, overflow: "auto" }}
    >
      <Table size="small" sx={{ width: "100%", tableLayout: "fixed" }}>
        <TableHead>
          <TableRow>
            {header.map((cell, idx) => (
              <TableCell
                key={idx}
                sx={{
                  fontWeight: "bold",
                  wordBreak: "break-word",
                  maxWidth: 150,
                }}
              >
                {cell}
              </TableCell>
            ))}
          </TableRow>
        </TableHead>
        <TableBody>
          {dataRows.map((row, ridx) => (
            <TableRow key={ridx}>
              {row.map((cell, cidx) => (
                <TableCell
                  key={cidx}
                  sx={{ wordBreak: "break-word", maxWidth: 150 }}
                >
                  {cell}
                </TableCell>
              ))}
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  );
};
