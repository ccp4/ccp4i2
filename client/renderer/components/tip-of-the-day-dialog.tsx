"use client";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
  CircularProgress,
  Stack,
  Typography,
} from "@mui/material";
import { EmojiObjects, Refresh } from "@mui/icons-material";
import React, { useCallback, useEffect, useState } from "react";
import { apiGet } from "../api-fetch";

interface TipData {
  id: number;
  html: string;
  count: number;
}

/**
 * Tip of the Day dialog. Fetches a random tip's body HTML from the tips API and
 * renders it; "Next tip" fetches another. Content is bundled, first-party HTML
 * (server/ccp4i2/tipsOfTheDay/*.html) with images served via the tips route.
 */
export function TipOfTheDayDialog({
  open,
  onClose,
}: {
  open: boolean;
  onClose: () => void;
}) {
  const [tip, setTip] = useState<TipData | null>(null);
  const [loading, setLoading] = useState(false);

  const fetchTip = useCallback(async () => {
    setLoading(true);
    try {
      const resp = await apiGet<any>("tips/");
      const data = resp?.data ?? resp;
      if (data?.html != null) setTip(data as TipData);
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    if (open && !tip) fetchTip();
  }, [open, tip, fetchTip]);

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle>
        <Stack direction="row" spacing={1} alignItems="center">
          <EmojiObjects fontSize="small" color="warning" />
          <span>Tip of the day</span>
        </Stack>
      </DialogTitle>
      <DialogContent dividers>
        {loading || !tip ? (
          <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
            <CircularProgress />
          </Box>
        ) : (
          <Box
            sx={{
              "& img": { maxWidth: "100%", height: "auto" },
              "& div": { width: "auto !important" },
            }}
            // Bundled first-party tip HTML.
            dangerouslySetInnerHTML={{ __html: tip.html }}
          />
        )}
      </DialogContent>
      <DialogActions sx={{ justifyContent: "space-between", px: 3 }}>
        <Typography variant="caption" color="text.secondary">
          {tip ? `Tip ${tip.id} of ${tip.count}` : ""}
        </Typography>
        <Stack direction="row" spacing={1}>
          <Button startIcon={<Refresh />} onClick={fetchTip} disabled={loading}>
            Next tip
          </Button>
          <Button onClick={onClose}>Close</Button>
        </Stack>
      </DialogActions>
    </Dialog>
  );
}
