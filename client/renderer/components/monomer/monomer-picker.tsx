"use client";

/**
 * A clickable 2D depiction of one monomer, for choosing an atom or a bond.
 *
 * Every drawn element carries its dictionary atom name, so a click yields
 * "NZ" or "C4A" directly -- which is exactly what the AceDRG task parameters
 * hold. Nothing here knows about jobs or parameters: it takes a monomer and
 * reports what was clicked, so the AceDRG tasks can share it.
 *
 * RDKit is used only to lay the molecule out. Its own get_svg() would key
 * the drawing by RDKit atom index, forcing a name -> index -> name round
 * trip, and gives no way to make bonds clickable.
 */

import { useMemo } from "react";
import { useTheme } from "@mui/material/styles";
import Box from "@mui/material/Box";
import Typography from "@mui/material/Typography";

import { useRDKit } from "../../providers/rdkit-provider";
import {
  buildMolblock,
  readCoordinates,
  type MonomerAtomDetail,
  type MonomerBond,
} from "../../lib/monomer-molblock";
import { bondKey, bondLines, layoutMonomer, type MonomerLayout } from "./monomer-layout";

export type PickMode = "atom" | "bond" | "none";

export interface MonomerPickerProps {
  atomDetails: MonomerAtomDetail[];
  bonds: MonomerBond[];
  /** What a click means. "none" renders a plain, non-interactive depiction. */
  mode?: PickMode;
  /** The atom currently chosen, drawn as selected. */
  selectedAtom?: string | null;
  /** The bond currently chosen, drawn as selected; order does not matter. */
  selectedBond?: { atom1: string; atom2: string } | null;
  /** Atoms to mark without selecting them, e.g. ones queued for deletion. */
  markedAtoms?: string[];
  onPickAtom?: (name: string) => void;
  onPickBond?: (atom1: string, atom2: string) => void;
  width?: number;
  height?: number;
  /** Shown above the drawing, e.g. the monomer code. */
  title?: string;
  emptyMessage?: string;
}

// Enough to tell the common elements apart at a glance. Anything else takes
// the default text colour, which follows the theme into dark mode.
const ELEMENT_COLOURS: Record<string, string> = {
  N: "#2563eb",
  O: "#dc2626",
  S: "#ca8a04",
  P: "#ea580c",
  F: "#16a34a",
  CL: "#16a34a",
  BR: "#a16207",
  I: "#7c3aed",
};

const PADDING = 26;
const HIT_RADIUS = 13;
const LABEL_RADIUS = 11;
const DOUBLE_BOND_GAP = 5;

function superscriptCharge(charge: number): string {
  if (!charge) return "";
  const sign = charge > 0 ? "+" : "-";
  const size = Math.abs(charge);
  return size === 1 ? sign : `${size}${sign}`;
}

export function MonomerPicker({
  atomDetails,
  bonds,
  mode = "none",
  selectedAtom = null,
  selectedBond = null,
  markedAtoms,
  onPickAtom,
  onPickBond,
  width = 320,
  height = 260,
  title,
  emptyMessage = "No monomer selected",
}: MonomerPickerProps) {
  const theme = useTheme();
  const { rdkitModule, isLoading } = useRDKit();

  const layout: MonomerLayout | null = useMemo(() => {
    if (!rdkitModule || atomDetails.length === 0) return null;
    const molblock = buildMolblock(atomDetails, bonds);
    // MinimalLib hands back raw WASM objects that must be released by hand;
    // a component that re-lays-out on every keystroke would otherwise leak.
    const mol = (rdkitModule as any).get_mol(molblock);
    if (!mol) return null;
    try {
      mol.set_new_coords(true);
      const coordinates = readCoordinates(mol.get_molblock());
      return layoutMonomer(atomDetails, bonds, coordinates, {
        width,
        height,
        padding: PADDING,
      });
    } catch {
      return null;
    } finally {
      mol.delete();
    }
  }, [rdkitModule, atomDetails, bonds, width, height]);

  const marked = useMemo(() => new Set(markedAtoms ?? []), [markedAtoms]);
  const selectedBondKey = selectedBond
    ? bondKey(selectedBond.atom1, selectedBond.atom2)
    : null;

  const lineColour = theme.palette.text.primary;
  const selectionColour = theme.palette.primary.main;
  const markColour = theme.palette.warning.main;

  if (atomDetails.length === 0) {
    return <Placeholder height={height} title={title} message={emptyMessage} />;
  }
  if (!layout) {
    return (
      <Placeholder
        height={height}
        title={title}
        message={
          isLoading
            ? "Preparing the chemistry toolkit..."
            : // Not a silent blank: say which of the two things went wrong.
              rdkitModule
              ? "This monomer could not be drawn from its dictionary"
              : "The chemistry toolkit could not be loaded"
        }
      />
    );
  }

  const atomsPickable = mode === "atom" && Boolean(onPickAtom);
  const bondsPickable = mode === "bond" && Boolean(onPickBond);

  return (
    <Box>
      {title && (
        <Typography variant="subtitle2" sx={{ mb: 0.5 }}>
          {title}
        </Typography>
      )}
      <Box
        component="svg"
        viewBox={`0 0 ${layout.width} ${layout.height}`}
        width="100%"
        role="group"
        aria-label={title ? `${title} structure` : "Monomer structure"}
        sx={{
          maxWidth: layout.width,
          border: 1,
          borderColor: "divider",
          borderRadius: 1,
          bgcolor: "background.paper",
          display: "block",
        }}
      >
        {layout.bonds.map((bond) => {
          const key = bondKey(bond.atom1, bond.atom2);
          const isSelected = key === selectedBondKey;
          const label = `${bond.atom1}-${bond.atom2}`;
          return (
            <g
              key={key}
              data-bond={label}
              role={bondsPickable ? "button" : undefined}
              tabIndex={bondsPickable ? 0 : undefined}
              aria-label={bondsPickable ? `Bond ${label}` : undefined}
              aria-pressed={bondsPickable ? isSelected : undefined}
              onClick={bondsPickable ? () => onPickBond!(bond.atom1, bond.atom2) : undefined}
              onKeyDown={
                bondsPickable
                  ? (event) => {
                      if (event.key === "Enter" || event.key === " ") {
                        event.preventDefault();
                        onPickBond!(bond.atom1, bond.atom2);
                      }
                    }
                  : undefined
              }
              style={{ cursor: bondsPickable ? "pointer" : "default" }}
            >
              {/* An invisible fat line under the bond: a 2px stroke is far
                  too thin to hit with a mouse, let alone a trackpad. */}
              {bondsPickable && (
                <line
                  x1={bond.x1}
                  y1={bond.y1}
                  x2={bond.x2}
                  y2={bond.y2}
                  stroke="transparent"
                  strokeWidth={12}
                />
              )}
              {bondLines(bond, DOUBLE_BOND_GAP).map((line, index) => (
                <line
                  key={index}
                  x1={line.x1}
                  y1={line.y1}
                  x2={line.x2}
                  y2={line.y2}
                  stroke={isSelected ? selectionColour : lineColour}
                  strokeWidth={isSelected ? 3 : 1.6}
                  strokeLinecap="round"
                />
              ))}
            </g>
          );
        })}

        {layout.atoms.map((atom) => {
          const isSelected = atom.name === selectedAtom;
          const isMarked = marked.has(atom.name);
          const colour =
            ELEMENT_COLOURS[atom.element?.toUpperCase() ?? ""] ?? lineColour;
          return (
            <g
              key={atom.name}
              data-atom={atom.name}
              role={atomsPickable ? "button" : undefined}
              tabIndex={atomsPickable ? 0 : undefined}
              aria-label={atomsPickable ? `Atom ${atom.name}` : undefined}
              aria-pressed={atomsPickable ? isSelected : undefined}
              onClick={atomsPickable ? () => onPickAtom!(atom.name) : undefined}
              onKeyDown={
                atomsPickable
                  ? (event) => {
                      if (event.key === "Enter" || event.key === " ") {
                        event.preventDefault();
                        onPickAtom!(atom.name);
                      }
                    }
                  : undefined
              }
              style={{ cursor: atomsPickable ? "pointer" : "default" }}
            >
              {atomsPickable && (
                <circle cx={atom.x} cy={atom.y} r={HIT_RADIUS} fill="transparent" />
              )}
              {/* Opaque disc so the label is legible where a bond runs under it. */}
              <circle
                cx={atom.x}
                cy={atom.y}
                r={LABEL_RADIUS}
                fill={theme.palette.background.paper}
                stroke={
                  isSelected ? selectionColour : isMarked ? markColour : "transparent"
                }
                strokeWidth={isSelected || isMarked ? 2.5 : 0}
              />
              {/* The dictionary name, not the element symbol: the name is
                  what the user is choosing and what the task stores. */}
              <text
                x={atom.x}
                y={atom.y}
                textAnchor="middle"
                dominantBaseline="central"
                fontSize={atom.name.length > 3 ? 8.5 : 10}
                fontFamily={theme.typography.fontFamily}
                fill={isSelected ? selectionColour : colour}
                fontWeight={isSelected ? 700 : 500}
                style={{ userSelect: "none", pointerEvents: "none" }}
              >
                {atom.name}
                {atom.charge !== 0 && (
                  <tspan fontSize={7} dy={-4}>
                    {superscriptCharge(atom.charge)}
                  </tspan>
                )}
              </text>
            </g>
          );
        })}
      </Box>
    </Box>
  );
}

function Placeholder({
  height,
  title,
  message,
}: {
  height: number;
  title?: string;
  message: string;
}) {
  return (
    <Box>
      {title && (
        <Typography variant="subtitle2" sx={{ mb: 0.5 }}>
          {title}
        </Typography>
      )}
      <Box
        sx={{
          height,
          border: 1,
          borderColor: "divider",
          borderRadius: 1,
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          px: 2,
          textAlign: "center",
        }}
      >
        <Typography variant="body2" color="text.secondary">
          {message}
        </Typography>
      </Box>
    </Box>
  );
}

export default MonomerPicker;
