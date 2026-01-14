"use client";

import { useEffect, useState } from "react";
import { Box, Skeleton, Typography } from "@mui/material";
import { useCCP4i2Window } from "../../app-context";

interface SmilesViewProps {
  smiles: string;
  width?: number;
  height?: number;
}

/**
 * Small inline SMILES renderer using RDKit WASM.
 * Renders a molecule as SVG from a SMILES string.
 */
export function SmilesView({ smiles, width = 100, height = 75 }: SmilesViewProps) {
  const { rdkitModule } = useCCP4i2Window();
  const [dataURI, setDataURI] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!rdkitModule || !smiles) {
      setDataURI(null);
      return;
    }

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (!mol) {
        setError("Invalid SMILES");
        setDataURI(null);
        return;
      }

      // Get SVG with custom dimensions
      const svg = mol.get_svg(width, height);
      mol.delete();

      const blob = new Blob([svg], { type: "image/svg+xml" });
      const url = URL.createObjectURL(blob);
      setDataURI(url);
      setError(null);

      // Cleanup blob URL when component unmounts or smiles changes
      return () => URL.revokeObjectURL(url);
    } catch (err) {
      console.error("RDKit error:", err);
      setError("Failed to render");
      setDataURI(null);
    }
  }, [smiles, rdkitModule, width, height]);

  if (!rdkitModule) {
    return (
      <Skeleton
        variant="rectangular"
        width={width}
        height={height}
        animation="wave"
      />
    );
  }

  if (error) {
    return (
      <Box
        sx={{
          width,
          height,
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          bgcolor: "background.paper",
          border: "1px solid",
          borderColor: "divider",
          borderRadius: 1,
        }}
      >
        <Typography variant="caption" color="text.secondary">
          {error}
        </Typography>
      </Box>
    );
  }

  if (!dataURI) {
    return (
      <Skeleton
        variant="rectangular"
        width={width}
        height={height}
        animation="wave"
      />
    );
  }

  return (
    <Box
      sx={{
        width,
        height,
        bgcolor: "background.paper",
        borderRadius: 1,
        overflow: "hidden",
      }}
    >
      <img
        src={dataURI}
        alt="Molecule structure"
        style={{
          width: "100%",
          height: "100%",
          objectFit: "contain",
        }}
      />
    </Box>
  );
}
