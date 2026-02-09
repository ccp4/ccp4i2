/**
 * Paste View Link field for syncing camera between Moorhen sessions.
 *
 * Accepts a pasted URL with a ?view= parameter and applies only the
 * camera parameters (origin, quaternion, zoom, clip, fog) without
 * changing representations or visibility, which may not apply across
 * different scenes.
 */

import React, { useState, useCallback } from "react";
import { useDispatch } from "react-redux";
import { TextField, InputAdornment, IconButton, Tooltip } from "@mui/material";
import { ContentPaste as PasteIcon } from "@mui/icons-material";
import {
  setOrigin,
  setQuat,
  setZoom,
  setClipStart,
  setClipEnd,
  setFogStart,
  setFogEnd,
  setRequestDrawScene,
} from "moorhen";
import { decodeViewState } from "../../lib/moorhen-view-state";
import { MOORHEN_DEFAULTS } from "../../types/moorhen-view-state";
import { usePopcorn } from "../../providers/popcorn-provider";

export const PasteViewLinkField: React.FC = () => {
  const dispatch = useDispatch();
  const { setMessage } = usePopcorn();
  const [value, setValue] = useState("");

  const applyViewFromUrl = useCallback(
    (url: string) => {
      try {
        const parsed = new URL(url);
        const viewParam = parsed.searchParams.get("view");
        if (!viewParam) {
          setMessage("No view parameter found in URL");
          return;
        }

        const viewState = decodeViewState(viewParam);
        if (!viewState) {
          setMessage("Failed to decode view state");
          return;
        }

        // Apply camera params only - skip representations (r), molecule (m)
        // and map (p) visibility which are scene-specific
        dispatch(setOrigin(viewState.o));
        dispatch(setQuat(viewState.q));
        dispatch(setZoom(viewState.z));
        dispatch(setClipStart(viewState.cs ?? MOORHEN_DEFAULTS.clipStart));
        dispatch(setClipEnd(viewState.ce ?? MOORHEN_DEFAULTS.clipEnd));
        dispatch(setFogStart(viewState.fs ?? MOORHEN_DEFAULTS.fogStart));
        dispatch(setFogEnd(viewState.fe ?? MOORHEN_DEFAULTS.fogEnd));
        dispatch(setRequestDrawScene(true));

        setValue("");
        setMessage("View applied");
      } catch {
        setMessage("Invalid URL");
      }
    },
    [dispatch, setMessage]
  );

  const handlePaste = useCallback(
    (e: React.ClipboardEvent) => {
      const text = e.clipboardData.getData("text").trim();
      if (text) {
        e.preventDefault();
        setValue(text);
        applyViewFromUrl(text);
      }
    },
    [applyViewFromUrl]
  );

  const handleReadClipboard = useCallback(async () => {
    try {
      const text = (await navigator.clipboard.readText()).trim();
      if (text) {
        setValue(text);
        applyViewFromUrl(text);
      }
    } catch {
      setMessage("Cannot read clipboard");
    }
  }, [applyViewFromUrl, setMessage]);

  return (
    <TextField
      size="small"
      placeholder="Paste view link..."
      value={value}
      onChange={(e) => setValue(e.target.value)}
      onPaste={handlePaste}
      onKeyDown={(e) => {
        if (e.key === "Enter" && value.trim()) {
          applyViewFromUrl(value.trim());
        }
      }}
      fullWidth
      InputProps={{
        endAdornment: (
          <InputAdornment position="end">
            <Tooltip title="Paste view link from clipboard">
              <IconButton size="small" onClick={handleReadClipboard} edge="end">
                <PasteIcon fontSize="small" />
              </IconButton>
            </Tooltip>
          </InputAdornment>
        ),
        sx: { fontSize: "0.75rem" },
      }}
    />
  );
};
