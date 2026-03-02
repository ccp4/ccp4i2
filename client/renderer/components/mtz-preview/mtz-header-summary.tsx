"use client";
import {
  Box,
  Chip,
  FormControl,
  InputLabel,
  MenuItem,
  Paper,
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

const cellParamLabels = ["a", "b", "c", "\u03B1", "\u03B2", "\u03B3"] as const;
const cellParamUnits = ["\u00C5", "\u00C5", "\u00C5", "\u00B0", "\u00B0", "\u00B0"] as const;

function formatCell(cell: NonNullable<MtzHeader["cell"]>) {
  const values = [cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma];
  return values.map((v, i) => ({
    label: cellParamLabels[i],
    value: v.toFixed(2),
    unit: cellParamUnits[i],
  }));
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

export const MtzHeaderSummary: React.FC<MtzHeaderSummaryProps> = ({
  header,
  onColorColumnChange,
  selectedColorColumn,
  autoDetectedColumn,
}) => {
  // Columns eligible for coloring: any numeric column except indices and batch
  const colorableTypes = new Set(["J", "F", "E", "D", "Q", "R", "K", "G", "W"]);
  const colorColumns = header.columns.filter((c) => colorableTypes.has(c.type));

  return (
    <Stack spacing={2} sx={{ height: "100%", overflowY: "auto" }}>
      {/* Title */}
      {header.title && (
        <Typography variant="body2" color="text.secondary" sx={{ fontStyle: "italic" }}>
          {header.title}
        </Typography>
      )}

      {/* Key stats */}
      <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
        {header.spaceGroup && (
          <Chip
            label={`${header.spaceGroup}${header.spaceGroupNumber ? ` (#${header.spaceGroupNumber})` : ""}`}
            size="small"
            color="primary"
            variant="outlined"
          />
        )}
        <Chip
          label={`${header.nReflections.toLocaleString()} reflections`}
          size="small"
          variant="outlined"
        />
        <Chip
          label={header.isMerged ? "Merged" : "Unmerged"}
          size="small"
          variant="outlined"
        />
        {header.resolution && (
          <Chip
            label={`${header.resolution[0].toFixed(1)} - ${header.resolution[1].toFixed(2)} \u00C5`}
            size="small"
            variant="outlined"
          />
        )}
      </Stack>

      {/* Unit Cell - compact 3-column layout */}
      {header.cell && (
        <Paper variant="outlined" sx={{ p: 1.5 }}>
          <Typography variant="subtitle2" gutterBottom>
            Unit Cell
          </Typography>
          <Table size="small" padding="none">
            <TableBody>
              <TableRow sx={{ "& td": { borderBottom: "none", py: 0.25 } }}>
                {formatCell(header.cell).slice(0, 3).map(({ label, value, unit }) => (
                  <TableCell key={label} sx={{ px: 0.5 }}>
                    <Typography component="span" sx={{ fontWeight: 500, fontStyle: "italic", mr: 0.5 }}>
                      {label}
                    </Typography>
                    <Typography component="span" sx={{ fontFamily: "monospace" }}>
                      {value}
                    </Typography>
                    <Typography component="span" sx={{ color: "text.secondary", ml: 0.25, fontSize: "0.8rem" }}>
                      {unit}
                    </Typography>
                  </TableCell>
                ))}
              </TableRow>
              <TableRow sx={{ "& td": { borderBottom: "none", py: 0.25 } }}>
                {formatCell(header.cell).slice(3, 6).map(({ label, value, unit }) => (
                  <TableCell key={label} sx={{ px: 0.5 }}>
                    <Typography component="span" sx={{ fontWeight: 500, fontStyle: "italic", mr: 0.5 }}>
                      {label}
                    </Typography>
                    <Typography component="span" sx={{ fontFamily: "monospace" }}>
                      {value}
                    </Typography>
                    <Typography component="span" sx={{ color: "text.secondary", ml: 0.25, fontSize: "0.8rem" }}>
                      {unit}
                    </Typography>
                  </TableCell>
                ))}
              </TableRow>
            </TableBody>
          </Table>
        </Paper>
      )}

      {/* Datasets */}
      {header.datasets.length > 0 && (
        <Accordion defaultExpanded={header.datasets.length <= 2} disableGutters>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography variant="subtitle2">
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
                    parts.push(`\u03BB ${ds.wavelength.toFixed(4)}\u00C5`);
                  return (
                    <TableRow key={ds.id}>
                      <TableCell sx={{ py: 0.25, px: 1, fontWeight: 500, fontSize: "0.8rem" }}>
                        {name}
                      </TableCell>
                      <TableCell sx={{ py: 0.25, px: 1, color: "text.secondary", fontSize: "0.75rem" }}>
                        {parts.join(" \u00B7 ")}
                      </TableCell>
                    </TableRow>
                  );
                })}
              </TableBody>
            </Table>
          </AccordionDetails>
        </Accordion>
      )}

      {/* Columns */}
      <Accordion disableGutters>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography variant="subtitle2">
            Columns ({header.columns.length})
          </Typography>
        </AccordionSummary>
        <AccordionDetails sx={{ p: 0, maxHeight: 300, overflowY: "auto" }}>
          <Table size="small" stickyHeader>
            <TableBody>
              {header.columns.map((col) => (
                <TableRow key={col.label} hover>
                  <TableCell sx={{ fontFamily: "monospace", fontWeight: 500, py: 0.5 }}>
                    {col.label}
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <Chip
                      label={col.type}
                      size="small"
                      variant="outlined"
                      sx={{ height: 20, fontSize: "0.7rem" }}
                      title={columnTypeNames[col.type] || col.type}
                    />
                  </TableCell>
                  <TableCell
                    sx={{ py: 0.5, color: "text.secondary", fontSize: "0.75rem" }}
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
        <Box sx={{ pt: 1 }}>
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
        </Box>
      )}
    </Stack>
  );
};
