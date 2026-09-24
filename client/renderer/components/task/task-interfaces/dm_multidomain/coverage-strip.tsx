import React, { useMemo } from "react";
import { Box, Typography, useTheme } from "@mui/material";

import { Assembly, ChainInfo, ParsedBody, referenceBounds, roleLabel } from "./dm-spec";

/**
 * The rigid bodies drawn on the reference copy's residues.
 *
 * Not a diagram of the data model — that would be the explanatory paragraph
 * with boxes round it, and would still have to be read and held in the head.
 * This is a picture of the user's own structure, and it can be wrong in a
 * visible way: a gap is residues no body claims (they get no averaging), an
 * overlap is two bodies claiming the same residues (their masks get split and
 * each operator is fitted to atoms the other also owns), and a body that
 * crosses entities shows up as one colour on two tracks.
 */

interface CoverageStripProps {
  assembly: Assembly;
  bodies: ParsedBody[];
  chains: ChainInfo[];
  /** Index of a body to pick out, e.g. the row the pointer is on. */
  highlight?: number | null;
}

const TRACK_HEIGHT = 18;
const ROW_GAP = 14;
const LABEL_WIDTH = 96;
const AXIS_PAD = 34;

export const bodyPalette = (theme: any): string[] => [
  theme.palette.primary.main,
  theme.palette.secondary.main,
  theme.palette.success.main,
  theme.palette.warning.main,
  theme.palette.info.main,
  theme.palette.error.main,
];

export const bodyColour = (theme: any, index: number): string => {
  const palette = bodyPalette(theme);
  return palette[index % palette.length];
};

export const CoverageStrip: React.FC<CoverageStripProps> = ({
  assembly,
  bodies,
  chains,
  highlight,
}) => {
  const theme = useTheme();

  const tracks = useMemo(
    () =>
      (assembly.roles.length ? assembly.roles : [])
        .map((role) => ({ role, bounds: referenceBounds(assembly, role, chains) }))
        .filter((track) => track.bounds !== null) as Array<{
        role: string;
        bounds: { lo: number; hi: number };
      }>,
    [assembly, chains]
  );

  if (!tracks.length) return null;

  const width = 640;
  const plotWidth = width - LABEL_WIDTH - AXIS_PAD;
  const height = tracks.length * (TRACK_HEIGHT + ROW_GAP);

  return (
    <Box sx={{ mt: 1 }}>
      <Box
        component="svg"
        viewBox={`0 0 ${width} ${height + 8}`}
        sx={{ width: "100%", maxWidth: width, height: "auto", display: "block" }}
        role="img"
        aria-label="Rigid bodies drawn on the residues of the reference copy"
      >
        <defs>
          <pattern
            id="dm-clash"
            width="6"
            height="6"
            patternUnits="userSpaceOnUse"
            patternTransform="rotate(45)"
          >
            <rect width="6" height="6" fill={theme.palette.error.main} fillOpacity={0.18} />
            <line
              x1="0"
              y1="0"
              x2="0"
              y2="6"
              stroke={theme.palette.error.main}
              strokeWidth="2"
            />
          </pattern>
        </defs>

        {tracks.map((track, row) => {
          const { lo, hi } = track.bounds;
          const span = Math.max(1, hi - lo);
          const y = row * (TRACK_HEIGHT + ROW_GAP);
          const x = (residue: number) =>
            LABEL_WIDTH + ((Math.min(Math.max(residue, lo), hi) - lo) / span) * plotWidth;

          // Residues claimed more than once, found by sweeping the segments of
          // this track: the overlap is the thing worth drawing differently.
          const spans = bodies.flatMap((body, index) =>
            body.segments
              .filter((s) => s.role === track.role)
              .map((s) => ({ lo: s.lo, hi: s.hi, index }))
          );
          const clashes: Array<{ lo: number; hi: number }> = [];
          for (let a = 0; a < spans.length; a += 1) {
            for (let b = a + 1; b < spans.length; b += 1) {
              if (spans[a].index === spans[b].index) continue;
              const start = Math.max(spans[a].lo, spans[b].lo);
              const end = Math.min(spans[a].hi, spans[b].hi);
              if (start <= end) clashes.push({ lo: start, hi: end });
            }
          }

          return (
            <g key={track.role}>
              <text
                x={0}
                y={y + TRACK_HEIGHT - 4}
                fill={theme.palette.text.primary}
                fontSize="12"
                fontFamily={theme.typography.fontFamily}
              >
                {roleLabel(track.role, assembly)}
              </text>

              {/* the chain itself: anything not covered is a gap you can see */}
              <rect
                x={LABEL_WIDTH}
                y={y}
                width={plotWidth}
                height={TRACK_HEIGHT}
                rx={3}
                fill={theme.palette.action.hover}
                stroke={theme.palette.divider}
              />

              {spans.map((s, i) => (
                <rect
                  key={`${s.index}-${i}`}
                  x={x(s.lo)}
                  y={y}
                  width={Math.max(2, x(s.hi) - x(s.lo))}
                  height={TRACK_HEIGHT}
                  rx={3}
                  fill={bodyColour(theme, s.index)}
                  fillOpacity={highlight == null || highlight === s.index ? 0.85 : 0.3}
                >
                  <title>
                    {`body ${s.index + 1}: ${s.lo}–${s.hi}`}
                  </title>
                </rect>
              ))}

              {clashes.map((clash, i) => (
                <rect
                  key={`clash-${i}`}
                  x={x(clash.lo)}
                  y={y}
                  width={Math.max(2, x(clash.hi) - x(clash.lo))}
                  height={TRACK_HEIGHT}
                  rx={3}
                  fill="url(#dm-clash)"
                  stroke={theme.palette.error.main}
                >
                  <title>{`residues ${clash.lo}–${clash.hi} are claimed by two bodies`}</title>
                </rect>
              ))}

              <text
                x={LABEL_WIDTH}
                y={y + TRACK_HEIGHT + 11}
                fill={theme.palette.text.secondary}
                fontSize="10"
                fontFamily={theme.typography.fontFamily}
              >
                {lo}
              </text>
              <text
                x={LABEL_WIDTH + plotWidth}
                y={y + TRACK_HEIGHT + 11}
                textAnchor="end"
                fill={theme.palette.text.secondary}
                fontSize="10"
                fontFamily={theme.typography.fontFamily}
              >
                {hi}
              </text>
            </g>
          );
        })}
      </Box>

      <Box sx={{ display: "flex", gap: 1.5, flexWrap: "wrap", mt: 0.5 }}>
        {bodies.map((body, index) => (
          <Box key={index} sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
            <Box
              sx={{
                width: 10,
                height: 10,
                borderRadius: "2px",
                bgcolor: bodyColour(theme, index),
                opacity: body.mode === "exclude" ? 0.35 : 1,
              }}
            />
            <Typography variant="caption" color="text.secondary">
              body {index + 1}
              {body.mode !== "average" ? ` (${body.mode})` : ""}
            </Typography>
          </Box>
        ))}
      </Box>
    </Box>
  );
};

export default CoverageStrip;
