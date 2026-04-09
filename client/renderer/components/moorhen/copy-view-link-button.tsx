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
