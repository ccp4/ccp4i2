"use client";
import { ImportProjectContent } from "@/components/import-project-content";
import { Box } from "@mui/material";
import { useTopBar } from "@/providers/top-bar-context";

export default function ImportProjectPage() {
  useTopBar({ title: "Import Project" });
  // The content area scrolls; the card centres itself with `margin: auto`,
  // which — unlike alignItems: center — collapses to 0 when the card is
  // taller than the area, so its top is never clipped.
  return (
    <Box sx={{ flex: 1, minHeight: 0, overflow: "auto", display: "flex" }}>
      <ImportProjectContent />
    </Box>
  );
}
