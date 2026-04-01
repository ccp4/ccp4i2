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

/**
 * Structure Info Modal for Moorhen viewer.
 *
 * Displays key metadata extracted from a loaded coordinate set (PDB/mmCIF):
 * - Identity: PDB ID, title, keywords, experimental method
 * - Crystal: space group, unit cell, resolution
 * - Refinement: R-work, R-free (mmCIF only)
 * - Composition: chains, ligands, atom/residue counts
 * - Citation: title, DOI, PubMed (mmCIF only)
 *
 * All parsing is done client-side using Moorhen's WASM gemmi and coot bindings.
 */

import React, { useEffect, useState } from "react";
import { moorhen } from "moorhen/types/moorhen";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  IconButton,
  Typography,
  Box,
  Table,
  TableBody,
  TableRow,
  TableCell,
  TableHead,
  Chip,
  Stack,
  Divider,
  CircularProgress,
  Link,
} from "@mui/material";
import { Close as CloseIcon } from "@mui/icons-material";
import { useTheme } from "../../theme/theme-provider";

// ── Types ──────────────────────────────────────────────────────────────────

interface StructureInfo {
  // Identity
  pdbId: string;
  title: string;
  keywords: string;
  experimentalMethod: string;
  depositionDate: string;

  // Crystal
  spacegroup: string;
  cell: { a: number; b: number; c: number; alpha: number; beta: number; gamma: number } | null;
  resolution: number;

  // Refinement (mmCIF only)
  rWork: string;
  rFree: string;
  resolutionLow: string;
  numReflections: string;

  // Composition
  chains: ChainSummary[];
  ligands: LigandSummary[];
  nModels: number;
  nAtoms: number;
  nResidues: number;
  nWater: number;

  // Citation (mmCIF only)
  citationTitle: string;
  doi: string;
  pubmedId: string;
  authorJournal: string;
}

interface ChainSummary {
  id: string;
  type: "protein" | "nucleic" | "saccharide" | "water" | "other";
  nResidues: number;
  nAtoms: number;
  sequence: string;
}

interface LigandSummary {
  chainId: string;
  name: string;
  seqNum: number;
  nAtoms: number;
}

interface StructureInfoModalProps {
  open: boolean;
  onClose: () => void;
  molecule: moorhen.Molecule | null;
}

// ── mmCIF text parsing ─────────────────────────────────────────────────────

/**
 * Extract a single-value mmCIF tag from raw text.
 * Handles both inline values and multi-line ;...;\n values.
 */
function extractCifValue(text: string, tag: string): string {
  // Try inline: _tag value
  const inlineRe = new RegExp(`^${tag.replace(/\./g, "\\.")}\\s+(.+)$`, "m");
  const inlineMatch = text.match(inlineRe);
  if (inlineMatch) {
    let val = inlineMatch[1].trim();
    // Strip quotes
    if ((val.startsWith("'") && val.endsWith("'")) ||
        (val.startsWith('"') && val.endsWith('"'))) {
      val = val.slice(1, -1);
    }
    return val === "?" || val === "." ? "" : val;
  }

  // Try multi-line: _tag\n;value\n;
  const multiRe = new RegExp(`^${tag.replace(/\./g, "\\.")}\\s*\\n;([\\s\\S]*?)\\n;`, "m");
  const multiMatch = text.match(multiRe);
  if (multiMatch) {
    return multiMatch[1].trim();
  }

  return "";
}

/**
 * Parse PDB REMARK records for R-factors (fallback for PDB format).
 */
function parseRemarksForRFactors(text: string): { rWork: string; rFree: string } {
  let rWork = "";
  let rFree = "";

  // REMARK   3   R VALUE            (WORKING SET) : 0.1735
  const rWorkMatch = text.match(/R VALUE\s+\(WORKING SET\)\s*:\s*([\d.]+)/);
  if (rWorkMatch) rWork = rWorkMatch[1];

  // REMARK   3   FREE R VALUE                     : 0.2147
  const rFreeMatch = text.match(/FREE R VALUE\s+:\s*([\d.]+)/);
  if (rFreeMatch) rFree = rFreeMatch[1];

  return { rWork, rFree };
}

// ── Amino acid / nucleic acid classification ───────────────────────────────

const AMINO_ACIDS = new Set([
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  "MSE", "SEC", "PYL", "UNK",
]);

const NUCLEIC_ACIDS = new Set([
  "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU",
  "ADE", "CYT", "GUA", "THY", "URA",
]);

const WATER = new Set(["HOH", "WAT", "H2O", "DOD"]);

const ONE_LETTER: Record<string, string> = {
  ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E",
  GLY: "G", HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F",
  PRO: "P", SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V", MSE: "M",
  SEC: "U", PYL: "O",
};

// ── Extraction logic ───────────────────────────────────────────────────────

async function extractStructureInfo(molecule: moorhen.Molecule): Promise<StructureInfo> {
  const mol = molecule as any; // access internal properties

  const info: StructureInfo = {
    pdbId: "", title: "", keywords: "", experimentalMethod: "", depositionDate: "",
    spacegroup: "", cell: null, resolution: -1,
    rWork: "", rFree: "", resolutionLow: "", numReflections: "",
    chains: [], ligands: [], nModels: 0, nAtoms: 0, nResidues: 0, nWater: 0,
    citationTitle: "", doi: "", pubmedId: "", authorJournal: "",
  };

  // 1. fetchHeaderInfo() — coot-side extraction (title, compound)
  try {
    const header = await molecule.fetchHeaderInfo(false);
    if (header) {
      info.title = header.title || "";
      if (header.compound_lines?.length && !info.title) {
        info.title = header.compound_lines.join(" ").trim();
      }
    }
  } catch {
    // fetchHeaderInfo may fail for some molecules — continue
  }

  // 2. gemmiStructure.get_info() — mmCIF-style metadata tags
  const gemmiStruct = mol.gemmiStructure;
  if (gemmiStruct && !gemmiStruct.isDeleted()) {
    try {
      const getTag = (tag: string): string => {
        try { return gemmiStruct.get_info(tag) || ""; }
        catch { return ""; }
      };

      info.pdbId = info.pdbId || getTag("_entry.id");
      if (!info.title) info.title = getTag("_struct.title");
      info.experimentalMethod = getTag("_exptl.method");
      info.depositionDate = getTag("_pdbx_database_status.recvd_initial_deposition_date");
      info.keywords = getTag("_struct_keywords.text");

      // Cell from gemmi if not already set
      if (!info.cell) {
        const cell = gemmiStruct.cell;
        if (cell && cell.a > 0) {
          info.cell = {
            a: cell.a, b: cell.b, c: cell.c,
            alpha: cell.alpha, beta: cell.beta, gamma: cell.gamma,
          };
        }
      }

      // PDB ID fallback from structure name
      if (!info.pdbId && gemmiStruct.name) {
        info.pdbId = gemmiStruct.name;
      }

      // Walk the model to extract composition
      const models = gemmiStruct.models;
      info.nModels = models.size();

      if (info.nModels > 0) {
        const model = models.get(0);
        const chainCount = model.chains.size();

        for (let ci = 0; ci < chainCount; ci++) {
          const chain = model.chains.get(ci);
          const chainId = chain.name;
          const residues = chain.residues;
          const nRes = residues.size();

          let chainType: ChainSummary["type"] = "other";
          let nAtoms = 0;
          let nWaterRes = 0;
          let seqParts: string[] = [];
          let ligands: LigandSummary[] = [];
          let hasProtein = false;
          let hasNucleic = false;

          for (let ri = 0; ri < nRes; ri++) {
            const residue = residues.get(ri);
            const resName = residue.name;
            const atomCount = residue.atoms.size();
            nAtoms += atomCount;

            if (AMINO_ACIDS.has(resName)) {
              hasProtein = true;
              seqParts.push(ONE_LETTER[resName] || "X");
            } else if (NUCLEIC_ACIDS.has(resName)) {
              hasNucleic = true;
              seqParts.push(resName.length <= 2 ? resName : resName.charAt(0));
            } else if (WATER.has(resName)) {
              nWaterRes++;
            } else if (atomCount > 4) {
              // Likely a ligand
              ligands.push({
                chainId, name: resName,
                seqNum: residue.seqid?.num?.value ?? 0,
                nAtoms: atomCount,
              });
            }

            residue.delete();
          }

          if (hasProtein) chainType = "protein";
          else if (hasNucleic) chainType = "nucleic";
          else if (nWaterRes > 0 && ligands.length === 0) chainType = "water";

          info.nAtoms += nAtoms;
          info.nResidues += nRes;
          info.nWater += nWaterRes;
          info.ligands.push(...ligands);

          // Only add non-water chains (or chains with ligands alongside water)
          if (chainType !== "water") {
            info.chains.push({
              id: chainId, type: chainType,
              nResidues: nRes - nWaterRes,
              nAtoms, sequence: seqParts.join(""),
            });
          }

          chain.delete();
          residues.delete();
        }
        model.delete();
      }
    } catch (e) {
      console.warn("Error extracting gemmi structure info:", e);
    }
  }

  // 3. Parse raw coordinate text for R-factors, citation (mmCIF or PDB remarks)
  try {
    const coordText = await molecule.getAtoms();
    if (coordText) {
      const ismmCIF = mol.coordsFormat === "mmcif";

      if (ismmCIF) {
        // mmCIF: extract structured fields
        if (!info.rWork) info.rWork = extractCifValue(coordText, "_refine.ls_R_factor_R_work");
        if (!info.rFree) info.rFree = extractCifValue(coordText, "_refine.ls_R_factor_R_free");
        if (!info.resolutionLow) info.resolutionLow = extractCifValue(coordText, "_refine.ls_d_res_low");
        if (!info.numReflections) info.numReflections = extractCifValue(coordText, "_reflns.number_obs");
        if (!info.citationTitle) info.citationTitle = extractCifValue(coordText, "_citation.title");
        if (!info.doi) info.doi = extractCifValue(coordText, "_citation.pdbx_database_id_DOI");
        if (!info.pubmedId) info.pubmedId = extractCifValue(coordText, "_citation.pdbx_database_id_PubMed");
      } else {
        // PDB format: parse REMARK records
        const { rWork, rFree } = parseRemarksForRFactors(coordText);
        if (!info.rWork) info.rWork = rWork;
        if (!info.rFree) info.rFree = rFree;
      }
    }
  } catch {
    // getAtoms() may fail — continue with what we have
  }

  return info;
}

// ── Formatting helpers ─────────────────────────────────────────────────────

function formatCell(cell: StructureInfo["cell"]): string {
  if (!cell || cell.a <= 0) return "";
  const fmt = (n: number) => n.toFixed(2);
  return `${fmt(cell.a)}  ${fmt(cell.b)}  ${fmt(cell.c)}    ${fmt(cell.alpha)}  ${fmt(cell.beta)}  ${fmt(cell.gamma)}`;
}

function formatResolution(res: number): string {
  if (res <= 0) return "";
  return `${res.toFixed(2)} \u00C5`;
}

function formatRFactor(val: string): string {
  if (!val) return "";
  const n = parseFloat(val);
  if (isNaN(n)) return val;
  // Display as percentage if < 1, otherwise as-is
  return n < 1 ? (n * 100).toFixed(1) + "%" : val;
}

function chainTypeLabel(type: ChainSummary["type"]): string {
  switch (type) {
    case "protein": return "Protein";
    case "nucleic": return "Nucleic acid";
    case "saccharide": return "Saccharide";
    default: return "Other";
  }
}

function chainTypeColor(type: ChainSummary["type"]): "primary" | "secondary" | "success" | "default" {
  switch (type) {
    case "protein": return "primary";
    case "nucleic": return "secondary";
    case "saccharide": return "success";
    default: return "default";
  }
}

// ── Sub-components ─────────────────────────────────────────────────────────

const InfoRow: React.FC<{ label: string; value: React.ReactNode; hide?: boolean }> = ({
  label, value, hide,
}) => {
  if (hide || !value || value === "") return null;
  return (
    <TableRow sx={{ "& td": { borderBottom: "none", py: 0.5 } }}>
      <TableCell sx={{ fontWeight: 500, color: "text.secondary", width: 160, whiteSpace: "nowrap", verticalAlign: "top" }}>
        {label}
      </TableCell>
      <TableCell>{value}</TableCell>
    </TableRow>
  );
};

const SectionHeader: React.FC<{ title: string }> = ({ title }) => (
  <Box sx={{ mt: 2, mb: 0.5 }}>
    <Typography variant="subtitle2" color="primary" sx={{ fontWeight: 600, textTransform: "uppercase", letterSpacing: 0.5, fontSize: "0.75rem" }}>
      {title}
    </Typography>
    <Divider />
  </Box>
);

const SequenceBlock: React.FC<{ sequence: string }> = ({ sequence }) => {
  if (!sequence) return null;
  // Format in blocks of 10 with position markers
  const blocks: string[] = [];
  for (let i = 0; i < sequence.length; i += 10) {
    blocks.push(sequence.substring(i, i + 10));
  }
  return (
    <Typography
      component="pre"
      sx={{
        fontFamily: "monospace", fontSize: "0.75rem", lineHeight: 1.6,
        whiteSpace: "pre-wrap", wordBreak: "break-all",
        m: 0, p: 0.5, bgcolor: "action.hover", borderRadius: 0.5,
        maxHeight: 80, overflow: "auto",
      }}
    >
      {blocks.join(" ")}
    </Typography>
  );
};

// ── Main modal ─────────────────────────────────────────────────────────────

export const StructureInfoModal: React.FC<StructureInfoModalProps> = ({
  open, onClose, molecule,
}) => {
  const { customColors } = useTheme();
  const [info, setInfo] = useState<StructureInfo | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!open || !molecule) {
      setInfo(null);
      setError(null);
      return;
    }

    let cancelled = false;
    setLoading(true);
    setError(null);

    extractStructureInfo(molecule).then((result) => {
      if (!cancelled) {
        setInfo(result);
        setLoading(false);
      }
    }).catch((err) => {
      if (!cancelled) {
        setError(err instanceof Error ? err.message : "Failed to extract structure info");
        setLoading(false);
      }
    });

    return () => { cancelled = true; };
  }, [open, molecule]);

  const displayTitle = info?.title || info?.pdbId || molecule?.name || "Structure";

  return (
    <Dialog
      open={open}
      onClose={onClose}
      maxWidth="md"
      fullWidth
      PaperProps={{
        sx: { maxHeight: "85vh" },
      }}
    >
      <DialogTitle sx={{ display: "flex", alignItems: "center", justifyContent: "space-between", pb: 0 }}>
        <Box sx={{ flex: 1, minWidth: 0 }}>
          <Stack direction="row" alignItems="center" spacing={1}>
            <Typography variant="h6" sx={{ fontWeight: 600 }}>
              Structure Info
            </Typography>
            {info?.pdbId && (
              <Chip label={info.pdbId} size="small" color="primary" variant="outlined"
                sx={{ fontFamily: "monospace", fontWeight: 600 }} />
            )}
          </Stack>
        </Box>
        <IconButton onClick={onClose} size="small">
          <CloseIcon />
        </IconButton>
      </DialogTitle>

      <DialogContent sx={{ pt: 1 }}>
        {loading && (
          <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
            <CircularProgress size={32} />
          </Box>
        )}

        {error && (
          <Typography color="error" sx={{ py: 2 }}>{error}</Typography>
        )}

        {info && !loading && (
          <>
            {/* ── Identity ── */}
            <SectionHeader title="Identity" />
            <Table size="small">
              <TableBody>
                <InfoRow label="Title" value={info.title} />
                <InfoRow label="Method" value={info.experimentalMethod} />
                <InfoRow label="Keywords" value={info.keywords} />
                <InfoRow label="Deposited" value={info.depositionDate} />
              </TableBody>
            </Table>

            {/* ── Crystal ── */}
            {(info.spacegroup || info.cell || info.resolution > 0) && (
              <>
                <SectionHeader title="Crystal" />
                <Table size="small">
                  <TableBody>
                    <InfoRow label="Space group" value={
                      info.spacegroup ? (
                        <Typography sx={{ fontFamily: "monospace" }}>{info.spacegroup}</Typography>
                      ) : null
                    } />
                    <InfoRow label="Unit cell" value={
                      info.cell && info.cell.a > 0 ? (
                        <Box>
                          <Typography sx={{ fontFamily: "monospace", fontSize: "0.875rem" }}>
                            {formatCell(info.cell)}
                          </Typography>
                          <Typography variant="caption" color="text.secondary">
                            a  b  c (\u00C5)    alpha  beta  gamma (\u00B0)
                          </Typography>
                        </Box>
                      ) : null
                    } />
                    <InfoRow label="Resolution" value={formatResolution(info.resolution)} />
                  </TableBody>
                </Table>
              </>
            )}

            {/* ── Refinement ── */}
            {(info.rWork || info.rFree) && (
              <>
                <SectionHeader title="Refinement" />
                <Table size="small">
                  <TableBody>
                    <InfoRow label="R-work" value={formatRFactor(info.rWork)} />
                    <InfoRow label="R-free" value={formatRFactor(info.rFree)} />
                    <InfoRow label="Resolution range" value={
                      info.resolutionLow && info.resolution > 0
                        ? `${info.resolutionLow} - ${info.resolution.toFixed(2)} \u00C5`
                        : null
                    } />
                    <InfoRow label="Reflections" value={
                      info.numReflections ? parseInt(info.numReflections).toLocaleString() : null
                    } />
                  </TableBody>
                </Table>
              </>
            )}

            {/* ── Composition ── */}
            <SectionHeader title="Composition" />
            <Table size="small">
              <TableBody>
                <InfoRow label="Models" value={info.nModels > 1 ? String(info.nModels) : null} />
                <InfoRow label="Total atoms" value={info.nAtoms > 0 ? info.nAtoms.toLocaleString() : null} />
                <InfoRow label="Water" value={info.nWater > 0 ? `${info.nWater.toLocaleString()} molecules` : null} />
              </TableBody>
            </Table>

            {/* Chains table */}
            {info.chains.length > 0 && (
              <Table size="small" sx={{ mt: 1 }}>
                <TableHead>
                  <TableRow>
                    <TableCell sx={{ fontWeight: 600, py: 0.5 }}>Chain</TableCell>
                    <TableCell sx={{ fontWeight: 600, py: 0.5 }}>Type</TableCell>
                    <TableCell sx={{ fontWeight: 600, py: 0.5 }} align="right">Residues</TableCell>
                    <TableCell sx={{ fontWeight: 600, py: 0.5 }} align="right">Atoms</TableCell>
                    <TableCell sx={{ fontWeight: 600, py: 0.5 }}>Sequence</TableCell>
                  </TableRow>
                </TableHead>
                <TableBody>
                  {info.chains.map((chain) => (
                    <TableRow key={chain.id} sx={{ "& td": { py: 0.5 } }}>
                      <TableCell>
                        <Typography sx={{ fontFamily: "monospace", fontWeight: 600 }}>
                          {chain.id}
                        </Typography>
                      </TableCell>
                      <TableCell>
                        <Chip label={chainTypeLabel(chain.type)} size="small"
                          color={chainTypeColor(chain.type)} variant="outlined"
                          sx={{ fontSize: "0.7rem", height: 20 }} />
                      </TableCell>
                      <TableCell align="right">{chain.nResidues}</TableCell>
                      <TableCell align="right">{chain.nAtoms.toLocaleString()}</TableCell>
                      <TableCell sx={{ maxWidth: 300 }}>
                        <SequenceBlock sequence={chain.sequence} />
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            )}

            {/* Ligands */}
            {info.ligands.length > 0 && (
              <>
                <Typography variant="subtitle2" sx={{ mt: 1.5, mb: 0.5, fontWeight: 600, fontSize: "0.8rem" }}>
                  Ligands
                </Typography>
                <Stack direction="row" spacing={0.5} sx={{ flexWrap: "wrap", gap: 0.5 }}>
                  {info.ligands.map((lig, i) => (
                    <Chip
                      key={`${lig.chainId}-${lig.name}-${lig.seqNum}-${i}`}
                      label={`${lig.name} (${lig.chainId}:${lig.seqNum}, ${lig.nAtoms} atoms)`}
                      size="small"
                      variant="outlined"
                      color="warning"
                      sx={{ fontSize: "0.75rem" }}
                    />
                  ))}
                </Stack>
              </>
            )}

            {/* ── Citation ── */}
            {(info.citationTitle || info.doi || info.pubmedId || info.authorJournal) && (
              <>
                <SectionHeader title="Citation" />
                <Table size="small">
                  <TableBody>
                    <InfoRow label="Citation" value={info.citationTitle} />
                    <InfoRow label="Authors" value={info.authorJournal} />
                    <InfoRow label="DOI" value={
                      info.doi ? (
                        <Link href={`https://doi.org/${info.doi}`} target="_blank" rel="noopener">
                          {info.doi}
                        </Link>
                      ) : null
                    } />
                    <InfoRow label="PubMed" value={
                      info.pubmedId ? (
                        <Link href={`https://pubmed.ncbi.nlm.nih.gov/${info.pubmedId}/`} target="_blank" rel="noopener">
                          {info.pubmedId}
                        </Link>
                      ) : null
                    } />
                  </TableBody>
                </Table>
              </>
            )}

            {/* ── Fallback notice ── */}
            {!info.pdbId && !info.experimentalMethod && !info.rWork && (
              <Box sx={{ mt: 2, p: 1.5, bgcolor: "action.hover", borderRadius: 1 }}>
                <Typography variant="body2" color="text.secondary">
                  This coordinate set does not contain PDB/mmCIF header records.
                  It may be an intermediate from a local refinement run.
                  Only composition information extracted from the atomic model is shown above.
                </Typography>
              </Box>
            )}
          </>
        )}
      </DialogContent>
    </Dialog>
  );
};
