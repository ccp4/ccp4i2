import { useEffect, useState } from "react";
import { useRDKit } from "../providers/rdkit-provider";
import { Typography } from "@mui/material";

interface RDKitViewProps {
  smiles: string;
  width?: number;
  height?: number;
}

/**
 * RDKitView - Renders a molecule structure from a SMILES string using RDKit WASM.
 */
export const RDKitView: React.FC<RDKitViewProps> = ({
  smiles,
  width = 250,
  height = 250
}) => {
  const { rdkitModule, isLoading, error: rdkitError } = useRDKit();
  const [dataURI, setDataURI] = useState<string | null>(null);

  useEffect(() => {
    if (!rdkitModule || !smiles) {
      setDataURI(null);
      return;
    }

    let currentUrl: string | null = null;

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (mol) {
        const svg = mol.get_svg(width, height);
        mol.delete();
        const blob = new Blob([svg], { type: "image/svg+xml" });
        currentUrl = URL.createObjectURL(blob);
        setDataURI(currentUrl);
      } else {
        setDataURI(null);
      }
    } catch (error) {
      console.error("RDKit error:", error);
      setDataURI(null);
    }

    // Cleanup blob URL when component unmounts or smiles changes
    return () => {
      if (currentUrl) {
        URL.revokeObjectURL(currentUrl);
      }
    };
  }, [smiles, rdkitModule, width, height]);

  if (rdkitError) {
    return (
      <Typography sx={{ color: "error.main" }}>
        Failed to load molecule viewer
      </Typography>
    );
  }

  if (isLoading) {
    return (
      <Typography sx={{ color: "text.secondary" }}>
        Loading molecule viewer...
      </Typography>
    );
  }

  return dataURI ? (
    <img src={dataURI} alt="Molecule structure" style={{ height, width }} />
  ) : (
    <Typography sx={{ color: "error.main" }}>No valid SMILES</Typography>
  );
};
