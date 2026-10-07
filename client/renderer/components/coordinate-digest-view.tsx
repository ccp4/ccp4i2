/**
 * The Digest of a coordinate file (file menu -> DIGEST), as tables rather
 * than the raw JSON it used to be shown as: one row per chain with its type,
 * residue count and range, then the ligands (#681).
 *
 * A chain's range is the range of its polymer residues where it has any.
 * The whole-chain range runs over waters and ligands numbered after (or
 * before) the polymer, so "1-38" could describe a 298-residue chain.
 */
import React from "react";
import {
  Box,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Typography,
} from "@mui/material";

export interface CoordinateChainDetail {
  id: string;
  type: string;
  nResidues: number;
  nAtoms?: number;
  firstRes: string;
  lastRes: string;
  ligandCount?: number;
  hasAltConf?: boolean;
  // Present from servers that split polymer from water and other residues
  nPolymer?: number;
  polymerFirstRes?: string;
  polymerLastRes?: string;
  nWater?: number;
  nOther?: number;
}

export interface CoordinateDigest {
  composition?: {
    nModels?: number;
    nChains?: number;
    nResidues?: number;
    nAtoms?: number;
    chainDetails?: CoordinateChainDetail[];
    ligands?: Array<{
      chain: string;
      name: string;
      seqNum: number;
      atomCount: number;
    }>;
    residueNameCounts?: Record<string, number>;
  };
  cell?: Record<string, number>;
  spaceGroup?: string;
}

export interface ChainDigestRow {
  chain: string;
  type: string;
  residues: number;
  range: string;
  waters: number | null;
  other: number | null;
}

const POLYMER_TYPE_LABELS: Record<string, string> = {
  protein: "Protein",
  nucleic: "Nucleic acid",
  saccharide: "Saccharide",
};

function rangeOf(first?: string, last?: string): string {
  if (!first && !last) return "";
  if (!last || first === last) return first || last || "";
  return `${first}–${last}`;
}

/** One table row per chain, from the digest's chainDetails. */
export function chainDigestRows(digest: CoordinateDigest): ChainDigestRow[] {
  return (digest.composition?.chainDetails ?? []).map((d) => {
    const split = d.nPolymer !== undefined;
    const hasPolymer = split && (d.nPolymer ?? 0) > 0;
    let type = POLYMER_TYPE_LABELS[d.type] ?? "Ligand/other";
    if (split && !hasPolymer) {
      const w = d.nWater ?? 0;
      const o = d.nOther ?? 0;
      type = o === 0 ? "Water" : w === 0 ? "Ligand/other" : "Ligand/water";
    } else if (d.type === "solvent") {
      type = "Ligand/water";
    }
    return {
      chain: d.id === "" ? "(blank)" : d.id,
      type,
      residues: hasPolymer ? (d.nPolymer as number) : d.nResidues,
      range: hasPolymer
        ? rangeOf(d.polymerFirstRes, d.polymerLastRes)
        : rangeOf(d.firstRes, d.lastRes),
      waters: hasPolymer ? d.nWater ?? 0 : null,
      other: hasPolymer ? d.nOther ?? 0 : null,
    };
  });
}

/** True when a digest is a coordinate file's (has per-chain details). */
export function isCoordinateDigest(digest: any): digest is CoordinateDigest {
  return Array.isArray(digest?.composition?.chainDetails);
}

const cellText = (cell?: Record<string, number>) =>
  cell && ["a", "b", "c", "alpha", "beta", "gamma"].every((k) => k in cell)
    ? `${cell.a} ${cell.b} ${cell.c}  ${cell.alpha} ${cell.beta} ${cell.gamma}`
    : null;

export const CoordinateDigestView: React.FC<{ digest: CoordinateDigest }> = ({
  digest,
}) => {
  const rows = chainDigestRows(digest);
  const composition = digest.composition ?? {};
  const ligands = composition.ligands ?? [];
  const showSplit = rows.some((r) => r.waters !== null);
  const cell = cellText(digest.cell);
  const residueNames = Object.entries(composition.residueNameCounts ?? {});

  return (
    <Stack spacing={2}>
      <Typography variant="body2" color="text.secondary">
        {[
          digest.spaceGroup && `Space group ${digest.spaceGroup}`,
          cell && `Cell ${cell}`,
          composition.nAtoms !== undefined && `${composition.nAtoms} atoms`,
          composition.nModels !== undefined &&
            composition.nModels > 1 &&
            `${composition.nModels} models (first shown)`,
        ]
          .filter(Boolean)
          .join(" · ")}
      </Typography>

      <TableContainer>
        <Table size="small" aria-label="Chains">
          <TableHead>
            <TableRow>
              <TableCell>Chain</TableCell>
              <TableCell>Type</TableCell>
              <TableCell align="right">Residues</TableCell>
              <TableCell>Range</TableCell>
              {showSplit && <TableCell align="right">Waters</TableCell>}
              {showSplit && (
                <TableCell align="right">Ligands/ions</TableCell>
              )}
            </TableRow>
          </TableHead>
          <TableBody>
            {rows.map((r, i) => (
              <TableRow key={`${r.chain}-${i}`}>
                <TableCell sx={{ fontFamily: "monospace" }}>{r.chain}</TableCell>
                <TableCell>{r.type}</TableCell>
                <TableCell align="right">{r.residues}</TableCell>
                <TableCell sx={{ fontFamily: "monospace" }}>{r.range}</TableCell>
                {showSplit && (
                  <TableCell align="right">{r.waters ?? ""}</TableCell>
                )}
                {showSplit && (
                  <TableCell align="right">{r.other ?? ""}</TableCell>
                )}
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>

      {ligands.length > 0 && (
        <Box>
          <Typography variant="subtitle2">Ligands</Typography>
          <Table size="small" aria-label="Ligands">
            <TableHead>
              <TableRow>
                <TableCell>Name</TableCell>
                <TableCell>Chain</TableCell>
                <TableCell align="right">Residue</TableCell>
                <TableCell align="right">Atoms</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {ligands.map((l, i) => (
                <TableRow key={`${l.chain}-${l.seqNum}-${l.name}-${i}`}>
                  <TableCell sx={{ fontFamily: "monospace" }}>{l.name}</TableCell>
                  <TableCell sx={{ fontFamily: "monospace" }}>
                    {l.chain === "" ? "(blank)" : l.chain}
                  </TableCell>
                  <TableCell align="right">{l.seqNum}</TableCell>
                  <TableCell align="right">{l.atomCount}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </Box>
      )}

      {residueNames.length > 0 && (
        <Box>
          <Typography variant="subtitle2">Residue names</Typography>
          <Typography
            variant="body2"
            color="text.secondary"
            sx={{ fontFamily: "monospace" }}
          >
            {residueNames
              .sort(([a], [b]) => a.localeCompare(b))
              .map(([name, n]) => `${name} ×${n}`)
              .join("  ")}
          </Typography>
        </Box>
      )}
    </Stack>
  );
};
