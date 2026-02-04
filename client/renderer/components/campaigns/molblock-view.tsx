"use client";

import { useEffect, useState } from "react";
import { Box, Skeleton, Tooltip, Typography } from "@mui/material";
import { useRDKit } from "../../providers/rdkit-provider";

interface MolBlockViewProps {
  molblock: string;
  name?: string;
  width?: number;
  height?: number;
}

/**
 * Renders a 2D molecule structure from MolBlock format using RDKit WASM.
 *
 * MolBlock is the MDL Molfile format which includes 2D coordinates,
 * making it ideal for rendering ligand structures from CIF dictionaries.
 */
export function MolBlockView({ molblock, name, width = 200, height = 150 }: MolBlockViewProps) {
  const { rdkitModule } = useRDKit();
  const [dataURI, setDataURI] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!rdkitModule || !molblock) {
      setDataURI(null);
      return;
    }

    let currentUrl: string | null = null;

    try {
      // RDKit's get_mol can accept MolBlock format directly
      const mol = rdkitModule.get_mol(molblock);
      if (!mol) {
        setError("Invalid molecule");
        setDataURI(null);
        return;
      }

      // Get SVG with custom dimensions
      const svg = mol.get_svg(width, height);
      mol.delete();

      const blob = new Blob([svg], { type: "image/svg+xml" });
      currentUrl = URL.createObjectURL(blob);
      setDataURI(currentUrl);
      setError(null);
    } catch (err) {
      console.error("RDKit error rendering MolBlock:", err);
      setError("Failed to render");
      setDataURI(null);
    }

    // Cleanup blob URL when component unmounts or molblock changes
    return () => {
      if (currentUrl) {
        URL.revokeObjectURL(currentUrl);
      }
    };
  }, [molblock, rdkitModule, width, height]);

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
    <Tooltip title={name || "Ligand structure"} placement="top">
      <Box
        sx={{
          width,
          height,
          bgcolor: "background.paper",
          borderRadius: 1,
          overflow: "hidden",
          border: "1px solid",
          borderColor: "divider",
        }}
      >
        <img
          src={dataURI}
          alt={name || "Ligand structure"}
          style={{
            width: "100%",
            height: "100%",
            objectFit: "contain",
          }}
        />
      </Box>
    </Tooltip>
  );
}
