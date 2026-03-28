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
/**
 * Copy View Link Button component for Moorhen viewer.
 *
 * Captures the current view state (camera position, orientation, zoom,
 * clip/fog planes, and molecule/map visibility) and copies a shareable
 * URL to the clipboard.
 */

import React, { useCallback } from "react";
import { Button, Tooltip } from "@mui/material";
import { Link as LinkIcon } from "@mui/icons-material";
import { usePopcorn } from "../../providers/popcorn-provider";

interface CopyViewLinkButtonProps {
  getViewUrl: () => string;
  disabled?: boolean;
}

export const CopyViewLinkButton: React.FC<CopyViewLinkButtonProps> = ({
  getViewUrl,
  disabled = false,
}) => {
  const { setMessage } = usePopcorn();

  const handleCopyViewLink = useCallback(async () => {
    try {
      const url = getViewUrl();
      await navigator.clipboard.writeText(url);
      setMessage("View link copied to clipboard", "success");
    } catch (error) {
      console.error("Failed to copy view link:", error);
      setMessage("Failed to copy view link", "error");
    }
  }, [getViewUrl, setMessage]);

  return (
    <Tooltip title="Copy a link that restores this exact view">
      <span>
        <Button
          variant="outlined"
          size="small"
          startIcon={<LinkIcon />}
          onClick={handleCopyViewLink}
          disabled={disabled}
          sx={{
            fontSize: "0.75rem",
            textTransform: "none",
          }}
        >
          Copy View Link
        </Button>
      </span>
    </Tooltip>
  );
};
