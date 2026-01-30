import React, { useCallback, useEffect, useState } from "react";
import { Box, Typography } from "@mui/material";

import { CCP4i2TaskElementProps } from "./task-element";
import { CSimpleTextFieldElement } from "./csimple-textfield";
import { useInferredVisibility } from "./hooks/useInferredVisibility";
import { useJob } from "../../../utils";
import { useRDKit } from "../../../providers/rdkit-provider";

/**
 * Molecule preview component that renders SMILES using RDKit.
 * Handles loading states and invalid SMILES gracefully.
 */
interface MoleculePreviewProps {
  smiles: string;
  size?: number;
}

const MoleculePreview: React.FC<MoleculePreviewProps> = ({
  smiles,
  size = 200,
}) => {
  const { rdkitModule, isLoading, error: rdkitError } = useRDKit();
  const [svgUrl, setSvgUrl] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    // Cleanup previous URL
    let currentUrl: string | null = null;

    if (!rdkitModule) {
      setError(null);
      setSvgUrl(null);
      return;
    }

    if (!smiles || smiles.trim() === "") {
      setError(null);
      setSvgUrl(null);
      return;
    }

    try {
      const mol = rdkitModule.get_mol(smiles.trim());
      if (mol) {
        const svg = mol.get_svg(size, size);
        mol.delete();
        const blob = new Blob([svg], { type: "image/svg+xml" });
        currentUrl = URL.createObjectURL(blob);
        setSvgUrl(currentUrl);
        setError(null);
      } else {
        setError("Invalid SMILES");
        setSvgUrl(null);
      }
    } catch (err) {
      console.error("RDKit error:", err);
      setError("Invalid SMILES");
      setSvgUrl(null);
    }

    // Cleanup blob URL on unmount or when smiles changes
    return () => {
      if (currentUrl) {
        URL.revokeObjectURL(currentUrl);
      }
    };
  }, [smiles, rdkitModule, size]);

  if (rdkitError) {
    return (
      <Typography variant="caption" color="error">
        Failed to load molecule viewer
      </Typography>
    );
  }

  if (isLoading || !rdkitModule) {
    return (
      <Typography variant="caption" color="text.secondary">
        Loading molecule viewer...
      </Typography>
    );
  }

  if (!smiles || smiles.trim() === "") {
    return (
      <Typography variant="caption" color="text.secondary">
        Enter a SMILES string to see structure preview
      </Typography>
    );
  }

  if (error) {
    return (
      <Typography variant="caption" color="error">
        {error}
      </Typography>
    );
  }

  if (svgUrl) {
    return (
      <img
        src={svgUrl}
        alt="Molecule structure"
        style={{ width: size, height: size }}
      />
    );
  }

  return null;
};

/**
 * CSMILESStringElement - A specialized input component for SMILES chemical notation strings.
 *
 * Features:
 * - Text input field for entering/editing SMILES strings
 * - Real-time molecule structure preview using RDKit
 * - Supports multiLine guiMode for longer SMILES strings
 * - Visual feedback for invalid SMILES
 * - onChange callback to update preview when value changes
 */
export const CSMILESStringElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, onChange: parentOnChange } = props;
  const isVisible = useInferredVisibility(props.visibility);

  // Get the current SMILES value from the task item
  const { useTaskItem } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  // Local state for the preview - starts with item value and updates on change
  const [previewSmiles, setPreviewSmiles] = useState<string>(
    item?._value?.toString() || ""
  );

  // Sync with item value when it changes externally
  useEffect(() => {
    setPreviewSmiles(item?._value?.toString() || "");
  }, [item?._value]);

  // Handle onChange - update local preview and propagate to parent
  const handleChange = useCallback(
    (updatedItem: any) => {
      // Update local preview immediately with the new value
      if (updatedItem?._value !== undefined) {
        setPreviewSmiles(updatedItem._value?.toString() || "");
      }
      // Propagate to parent onChange if provided
      if (parentOnChange) {
        parentOnChange(updatedItem);
      }
    },
    [parentOnChange]
  );

  if (!isVisible) {
    return null;
  }

  return (
    <Box sx={{ display: "flex", flexDirection: "column", gap: 1 }}>
      {/* SMILES text input - using CSimpleTextFieldElement directly */}
      <CSimpleTextFieldElement
        {...props}
        type="text"
        onChange={handleChange}
      />

      {/* Molecule structure preview */}
      <Box
        sx={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          minHeight: 100,
          border: "1px solid",
          borderColor: "divider",
          borderRadius: 1,
          backgroundColor: "background.paper",
          p: 1,
        }}
      >
        <MoleculePreview smiles={previewSmiles} size={200} />
      </Box>
    </Box>
  );
};

CSMILESStringElement.displayName = "CSMILESStringElement";

export default CSMILESStringElement;
