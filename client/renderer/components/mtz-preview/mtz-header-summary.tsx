/*
 * Copyright (C) 2026 Newcastle University
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
"use client";
import {
  Box,
  Chip,
  Divider,
  FormControl,
  InputLabel,
  MenuItem,
  Select,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableRow,
  Typography,
  Accordion,
  AccordionSummary,
  AccordionDetails,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import type { MtzHeader } from "../../lib/mtz-parser";

interface MtzHeaderSummaryProps {
  header: MtzHeader;
  onColorColumnChange?: (columnLabel: string | null) => void;
  selectedColorColumn?: string | null;
  /** Column label auto-detected by the parser (shown as default in dropdown) */
  autoDetectedColumn?: string | null;
}

const columnTypeNames: Record<string, string> = {
  H: "Index",
  J: "Intensity",
  F: "Amplitude",
  D: "Anomalous diff",
  Q: "Sigma",
  P: "Phase",
  W: "Weight",
  A: "HL coeff",
  B: "Batch",
  I: "Integer",
  R: "Real",
  K: "Anom I",
  G: "Anom F",
  L: "Anom D",
  M: "Anom sigma",
  E: "Normalised F",
  Y: "M/ISYM",
};

/** Compact display of a cell parameter: italic label, monospace value, small unit */
function CellParam({ label, value, unit }: { label: string; value: string; unit: string }) {
  return (
    <Box sx={{ display: "inline-flex", alignItems: "baseline", mr: 2, mb: 0.25 }}>
      <Typography
        component="span"
        sx={{ fontStyle: "italic", fontWeight: 600, fontSize: "0.85rem", mr: 0.5, color: "text.secondary" }}
      >
        {label}
      </Typography>
      <Typography component="span" sx={{ fontFamily: "monospace", fontSize: "0.9rem" }}>
        {value}
      </Typography>
      <Typography component="span" sx={{ color: "text.disabled", fontSize: "0.7rem", ml: 0.25 }}>
        {unit}
      </Typography>
    </Box>
  );
}

export const MtzHeaderSummary: React.FC<MtzHeaderSummaryProps> = ({
  header,
  onColorColumnChange,
  selectedColorColumn,
  autoDetectedColumn,
}) => {
  // Columns eligible for coloring: any numeric column except indices and batch
  const colorableTypes = new Set(["J", "F", "E", "D", "Q", "R", "K", "G", "W"]);
  const colorColumns = header.columns.filter((c) => colorableTypes.has(c.type));

  // Format resolution string
  const resolutionStr = header.resolution
    ? `${header.resolution[0].toFixed(1)}\u2009\u2013\u2009${header.resolution[1].toFixed(2)}\u2009\u00C5`
    : null;

  // Datasets with wavelength
  const datasetsWithWavelength = header.datasets.filter(
    (ds) => ds.wavelength !== undefined && ds.wavelength > 0
  );

  // Distinct wavelengths
  const wavelengths = [...new Set(datasetsWithWavelength.map((ds) => ds.wavelength!))];

  return (
    <Stack spacing={1.5} sx={{ height: "100%", overflowY: "auto" }}>
      {/* Provenance / title */}
      {header.title && (
        <Typography variant="body2" color="text.secondary" sx={{ fontStyle: "italic", lineHeight: 1.3 }}>
          {header.title}
        </Typography>
      )}

      {/* Space group — prominent */}
      {header.spaceGroup && (
        <Box>
          <Typography
            sx={{
              fontFamily: "monospace",
              fontSize: "1.1rem",
              fontWeight: 600,
              letterSpacing: "0.02em",
            }}
          >
            {header.spaceGroup}
            {header.spaceGroupNumber ? (
              <Typography
                component="span"
                sx={{ color: "text.secondary", fontWeight: 400, fontSize: "0.85rem", ml: 0.75 }}
              >
                #{header.spaceGroupNumber}
              </Typography>
            ) : null}
          </Typography>
        </Box>
      )}

      {/* Key stats in a compact row */}
      <Stack direction="row" spacing={0.75} flexWrap="wrap" useFlexGap alignItems="center">
        <Chip
          label={`${header.nReflections.toLocaleString()} reflections`}
          size="small"
          variant="outlined"
        />
        <Chip
          label={header.isMerged ? "Merged" : "Unmerged"}
          size="small"
          variant="outlined"
          color={header.isMerged ? "default" : "info"}
        />
        {!header.isMerged && header.nBatches > 0 && (
          <Chip
            label={`${header.nBatches} batches`}
            size="small"
            variant="outlined"
          />
        )}
        {resolutionStr && (
          <Chip label={resolutionStr} size="small" variant="outlined" />
        )}
        {wavelengths.length === 1 && (
          <Chip
            label={`\u03BB\u2009${wavelengths[0].toFixed(4)}\u2009\u00C5`}
            size="small"
            variant="outlined"
          />
        )}
      </Stack>

      <Divider />

      {/* Unit cell — compact inline layout */}
      {header.cell && (
        <Box>
          <Typography variant="caption" color="text.secondary" sx={{ mb: 0.5, display: "block" }}>
            Unit Cell
          </Typography>
          <Box sx={{ display: "flex", flexWrap: "wrap" }}>
            <CellParam label="a" value={header.cell.a.toFixed(2)} unit={"\u00C5"} />
            <CellParam label="b" value={header.cell.b.toFixed(2)} unit={"\u00C5"} />
            <CellParam label="c" value={header.cell.c.toFixed(2)} unit={"\u00C5"} />
          </Box>
          <Box sx={{ display: "flex", flexWrap: "wrap" }}>
            <CellParam label={"\u03B1"} value={header.cell.alpha.toFixed(2)} unit={"\u00B0"} />
            <CellParam label={"\u03B2"} value={header.cell.beta.toFixed(2)} unit={"\u00B0"} />
            <CellParam label={"\u03B3"} value={header.cell.gamma.toFixed(2)} unit={"\u00B0"} />
          </Box>
        </Box>
      )}

      {/* Datasets — shown if more than one or if multi-wavelength */}
      {header.datasets.length > 0 && (wavelengths.length > 1 || header.datasets.length > 1) && (
        <>
          <Divider />
          <Accordion defaultExpanded={header.datasets.length <= 3} disableGutters elevation={0}>
            <AccordionSummary expandIcon={<ExpandMoreIcon />} sx={{ minHeight: 0, px: 0, "& .MuiAccordionSummary-content": { my: 0.5 } }}>
              <Typography variant="caption" color="text.secondary">
                Datasets ({header.datasets.length})
              </Typography>
            </AccordionSummary>
            <AccordionDetails sx={{ p: 0 }}>
              <Table size="small" padding="none">
                <TableBody>
                  {header.datasets.map((ds) => {
                    const name = ds.datasetName || ds.crystalName || `Dataset ${ds.id}`;
                    const parts: string[] = [];
                    if (ds.crystalName && ds.crystalName !== name) parts.push(ds.crystalName);
                    if (ds.wavelength !== undefined && ds.wavelength > 0)
                      parts.push(`\u03BB\u2009${ds.wavelength.toFixed(4)}\u2009\u00C5`);
                    return (
                      <TableRow key={ds.id}>
                        <TableCell sx={{ py: 0.25, px: 0.5, fontWeight: 500, fontSize: "0.8rem" }}>
                          {name}
                        </TableCell>
                        <TableCell sx={{ py: 0.25, px: 0.5, color: "text.secondary", fontSize: "0.75rem" }}>
                          {parts.join(" \u00B7 ")}
                        </TableCell>
                      </TableRow>
                    );
                  })}
                </TableBody>
              </Table>
            </AccordionDetails>
          </Accordion>
        </>
      )}

      {/* Columns */}
      <Divider />
      <Accordion disableGutters elevation={0}>
        <AccordionSummary expandIcon={<ExpandMoreIcon />} sx={{ minHeight: 0, px: 0, "& .MuiAccordionSummary-content": { my: 0.5 } }}>
          <Typography variant="caption" color="text.secondary">
            Columns ({header.columns.length})
          </Typography>
        </AccordionSummary>
        <AccordionDetails sx={{ p: 0, maxHeight: 300, overflowY: "auto" }}>
          <Table size="small" stickyHeader>
            <TableBody>
              {header.columns.map((col) => (
                <TableRow key={col.label} hover>
                  <TableCell sx={{ fontFamily: "monospace", fontWeight: 500, py: 0.25, fontSize: "0.8rem" }}>
                    {col.label}
                  </TableCell>
                  <TableCell sx={{ py: 0.25 }}>
                    <Chip
                      label={col.type}
                      size="small"
                      variant="outlined"
                      sx={{ height: 18, fontSize: "0.65rem" }}
                      title={columnTypeNames[col.type] || col.type}
                    />
                  </TableCell>
                  <TableCell
                    sx={{ py: 0.25, color: "text.secondary", fontSize: "0.7rem" }}
                  >
                    {columnTypeNames[col.type] || ""}
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </AccordionDetails>
      </Accordion>

      {/* Color column selector */}
      {colorColumns.length > 0 && onColorColumnChange && (
        <>
          <Divider />
          <FormControl size="small" fullWidth>
            <InputLabel>Color by</InputLabel>
            <Select
              value={selectedColorColumn === null ? "__auto__" : selectedColorColumn}
              label="Color by"
              onChange={(e) => {
                const v = e.target.value;
                if (v === "__auto__") onColorColumnChange(null);
                else if (v === "") onColorColumnChange("");
                else onColorColumnChange(v);
              }}
            >
              <MenuItem value="__auto__">
                <em>Auto{autoDetectedColumn ? ` (${autoDetectedColumn})` : ""}</em>
              </MenuItem>
              <MenuItem value="">
                <em>None (uniform color)</em>
              </MenuItem>
              {colorColumns.map((col) => (
                <MenuItem key={col.label} value={col.label}>
                  {col.label} ({columnTypeNames[col.type] || col.type})
                </MenuItem>
              ))}
            </Select>
          </FormControl>
        </>
      )}
    </Stack>
  );
};
