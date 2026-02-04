"use client";

import { useEffect, useState } from "react";
import { Box, Skeleton, Typography } from "@mui/material";
import { MolBlockView } from "./molblock-view";
import { apiGet } from "../../api-fetch";

interface Ligand2DViewProps {
  /** File ID of the ligand dictionary CIF file */
  fileId: number;
  /** Display name for the ligand (fallback if ligand_code not returned by API) */
  name?: string;
  width?: number;
  height?: number;
  /** Callback when ligand code is loaded from the CIF file */
  onLigandCodeLoaded?: (code: string) => void;
}

/**
 * Fetches and renders a 2D structure of a ligand from its dictionary CIF file.
 *
 * This component:
 * 1. Calls the /files/{id}/molblock/ API to convert CIF to MolBlock
 * 2. Renders the MolBlock using RDKit WASM
 * 3. Displays the ligand code extracted from the CIF file
 */
export function Ligand2DView({
  fileId,
  name,
  width = 200,
  height = 150,
  onLigandCodeLoaded,
}: Ligand2DViewProps) {
  const [molblock, setMolblock] = useState<string | null>(null);
  const [ligandCode, setLigandCode] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    async function fetchMolblock() {
      setLoading(true);
      setError(null);

      try {
        const response = await apiGet(`files/${fileId}/molblock/`);

        if (cancelled) return;

        // API returns {success: true, data: {molblock, ligand_code}} or {success: false, error: "..."}
        if (response?.success && response.data?.molblock) {
          setMolblock(response.data.molblock);
          if (response.data.ligand_code) {
            setLigandCode(response.data.ligand_code);
            onLigandCodeLoaded?.(response.data.ligand_code);
          }
        } else if (response?.error) {
          setError(response.error);
        } else if (!response?.success) {
          setError(response?.message || "Failed to convert ligand");
        } else {
          setError("No molblock returned");
        }
      } catch (err) {
        if (cancelled) return;
        console.error("Failed to fetch molblock:", err);
        setError("Failed to load structure");
      } finally {
        if (!cancelled) {
          setLoading(false);
        }
      }
    }

    fetchMolblock();

    return () => {
      cancelled = true;
    };
  }, [fileId, onLigandCodeLoaded]);

  // Use ligand code from CIF if available, otherwise fall back to provided name
  const displayName = ligandCode || name;

  if (loading) {
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
        <Typography variant="caption" color="text.secondary" sx={{ textAlign: "center", p: 1 }}>
          {error}
        </Typography>
      </Box>
    );
  }

  if (!molblock) {
    return null;
  }

  return <MolBlockView molblock={molblock} name={displayName} width={width} height={height} />;
}
