import { useEffect, useMemo, useState } from "react";
import {
  Alert,
  Box,
  Chip,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Typography,
  useTheme,
} from "@mui/material";

export interface SegmentRow {
  index: number;
  BEG: number;
  END: number;
  KIND: string;
  length: number;
  first: string;
  last: string;
  reference: string;
  problems: string[];
}

export interface SequenceView {
  ok: boolean;
  error?: string;
  name?: string;
  sequence?: string;
  length?: number;
  segments?: SegmentRow[];
  uncovered?: [number, number][];
  overlaps?: [number, number, number[]][];
  problems?: string[];
  modelled?: number[] | null;
}

const LINE = 60;

/** Theme colour for each segment kind; unknown kinds fall back to grey. */
const useKindColours = (): Record<string, string> => {
  const theme = useTheme();
  return useMemo(
    () => ({
      tag: theme.palette.grey[400],
      linker: theme.palette.grey[300],
      natural: theme.palette.primary.light,
      engineered: theme.palette.secondary.light,
      designed: theme.palette.success.light,
      scaffold: theme.palette.info.light,
      synthetic_peptide: theme.palette.warning.light,
    }),
    [theme]
  );
};

/** The segment a 1-based position falls in (the first, if several). */
const kindAt = (segments: SegmentRow[], pos: number): SegmentRow | undefined =>
  segments.find((s) => s.BEG <= pos && pos <= s.END);

const fmtRange = ([a, b]: [number, number]) => (a === b ? `${a}` : `${a}-${b}`);

/**
 * Read-only picture of one molecule record: its AU-contents sequence in
 * lines of 60 with residue numbers and a tick every 10, each residue tinted
 * by its segment (uncovered residues hatched), a thin track of what the
 * model covers, a legend, a table of the segments, and the problems the
 * parsing shows. Fetched from the server (sequence_view) whenever the
 * caller's refresh key changes; it never writes anything.
 */
export const DepositionSequenceView: React.FC<{
  sequenceName: string;
  title: string;
  refreshKey: string;
  fetchView: (name: string) => Promise<SequenceView | null>;
}> = ({ sequenceName, title, refreshKey, fetchView }) => {
  const theme = useTheme();
  const colours = useKindColours();
  const [view, setView] = useState<SequenceView | null>(null);

  useEffect(() => {
    let live = true;
    fetchView(sequenceName).then((v) => {
      if (live) setView(v);
    });
    return () => {
      live = false;
    };
  }, [sequenceName, refreshKey, fetchView]);

  if (!view) return null;
  if (!view.ok || !view.sequence) {
    return (
      <Alert severity="warning" sx={{ mb: 1 }}>
        {title}: {view.error ?? `no sequence found for ${sequenceName}`}
      </Alert>
    );
  }

  const seq = view.sequence;
  const segments = view.segments ?? [];
  const modelled = view.modelled ? new Set(view.modelled) : null;
  const overlapPositions = new Set<number>();
  for (const [a, b] of view.overlaps ?? []) for (let p = a; p <= b; p++) overlapPositions.add(p);
  const lines: number[] = [];
  for (let start = 1; start <= seq.length; start += LINE) lines.push(start);
  const hatch = `repeating-linear-gradient(45deg, ${theme.palette.grey[300]} 0 3px, transparent 3px 6px)`;
  const kinds = Array.from(new Set(segments.map((s) => s.KIND)));

  return (
    <Paper variant="outlined" sx={{ p: 1, mb: 1 }}>
      <Typography variant="subtitle2">
        {title} ({seq.length} residues)
      </Typography>
      <Box sx={{ fontFamily: "monospace", fontSize: "0.8rem", overflowX: "auto", my: 0.5 }}>
        {lines.map((start) => {
          const end = Math.min(start + LINE - 1, seq.length);
          return (
            <Box key={start} sx={{ display: "flex", alignItems: "flex-start", mb: 0.5 }}>
              <Box sx={{ width: "3.5em", textAlign: "right", pr: 1, color: "text.secondary" }}>{start}</Box>
              <Box>
                <Box sx={{ display: "flex" }}>
                  {Array.from(seq.slice(start - 1, end)).map((res, k) => {
                    const pos = start + k;
                    const seg = kindAt(segments, pos);
                    const tenth = pos % 10 === 0;
                    return (
                      <Box
                        key={pos}
                        title={`${res}${pos}${seg ? " " + seg.KIND : " (no segment)"}`}
                        sx={{
                          width: "1ch",
                          textAlign: "center",
                          background: seg ? colours[seg.KIND] ?? theme.palette.grey[200] : hatch,
                          outline: overlapPositions.has(pos) ? `1px solid ${theme.palette.error.main}` : "none",
                          borderRight: tenth ? `1px solid ${theme.palette.text.disabled}` : "none",
                          color: theme.palette.getContrastText(
                            seg ? colours[seg.KIND] ?? theme.palette.grey[200] : theme.palette.background.paper
                          ),
                        }}
                      >
                        {res}
                      </Box>
                    );
                  })}
                </Box>
                {modelled && (
                  <Box sx={{ display: "flex", height: 3, mt: "1px" }}>
                    {Array.from({ length: end - start + 1 }, (_, k) => (
                      <Box
                        key={k}
                        sx={{
                          width: "1ch",
                          background: modelled.has(start + k) ? theme.palette.success.main : "transparent",
                        }}
                      />
                    ))}
                  </Box>
                )}
              </Box>
            </Box>
          );
        })}
      </Box>
      <Box sx={{ display: "flex", flexWrap: "wrap", gap: 0.5, mb: 0.5 }}>
        {kinds.map((kind) => (
          <Chip key={kind} size="small" label={kind} sx={{ background: colours[kind] ?? theme.palette.grey[200] }} />
        ))}
        {(view.uncovered ?? []).length > 0 && <Chip size="small" label="no segment" sx={{ background: hatch }} />}
        {modelled && (
          <Chip
            size="small"
            variant="outlined"
            label="bar underneath: modelled in at least one chain"
            sx={{ borderColor: theme.palette.success.main }}
          />
        )}
      </Box>
      <Table size="small">
        <TableHead>
          <TableRow>
            <TableCell>Segment</TableCell>
            <TableCell>Kind</TableCell>
            <TableCell>Range</TableCell>
            <TableCell>Length</TableCell>
            <TableCell>Residues</TableCell>
            <TableCell>Reference</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {segments.map((s) => (
            <TableRow key={s.index}>
              <TableCell>{s.index + 1}</TableCell>
              <TableCell>{s.KIND}</TableCell>
              <TableCell>
                {s.BEG}-{s.END}
              </TableCell>
              <TableCell>{s.length}</TableCell>
              <TableCell sx={{ fontFamily: "monospace" }}>
                {s.length > 10 ? `${s.first}...${s.last}` : s.first}
              </TableCell>
              <TableCell>{s.reference}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
      {segments.length === 0 && (
        <Typography variant="body2" color="text.secondary">
          No segments yet: use Fill from the AU contents, or add them below.
        </Typography>
      )}
      {(view.problems ?? []).map((p) => (
        <Alert key={p} severity="warning" sx={{ mt: 0.5, py: 0 }}>
          {p}
        </Alert>
      ))}
      {(view.uncovered ?? []).length > 0 && (
        <Typography variant="caption" color="text.secondary">
          Uncovered: {(view.uncovered ?? []).map(fmtRange).join(", ")}
        </Typography>
      )}
    </Paper>
  );
};
