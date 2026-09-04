"use client";

import { ReactNode, useState } from "react";
import { Paper, SxProps, Theme } from "@mui/material";

interface DropZoneProps {
  /** Called with every file dropped or picked, whether one or many. */
  onFilesSelected: (files: FileList) => void;
  accept?: string;
  /** Whether the underlying file input allows selecting more than one file. */
  multiple?: boolean;
  disabled?: boolean;
  /** Border tint hint independent of drag-over state, e.g. once a file is staged. */
  accent?: "default" | "success";
  children: ReactNode;
  sx?: SxProps<Theme>;
}

/**
 * Chrome shared by every drag-and-drop file target in the app: a dashed
 * Paper that opens a file picker on click and accepts a native drop. Content
 * (placeholder text, staged file list, chips, etc.) is left to the caller via
 * children, since that varies per use. Uses theme palette tokens throughout
 * (never hardcoded hex) so hover/drag states read correctly in dark mode too.
 */
export function DropZone({
  onFilesSelected,
  accept,
  multiple = false,
  disabled = false,
  accent = "default",
  children,
  sx,
}: DropZoneProps) {
  const [isDragOver, setIsDragOver] = useState(false);

  const openFilePicker = () => {
    if (disabled) return;
    const input = document.createElement("input");
    input.type = "file";
    input.multiple = multiple;
    if (accept) input.accept = accept;
    input.onchange = (e) => {
      const files = (e.target as HTMLInputElement).files;
      if (files && files.length > 0) onFilesSelected(files);
    };
    input.click();
  };

  return (
    <Paper
      variant="outlined"
      sx={{
        border: "2px dashed",
        borderColor: isDragOver
          ? "primary.main"
          : accent === "success"
            ? "success.main"
            : "divider",
        backgroundColor: isDragOver ? "action.hover" : "transparent",
        textAlign: "center",
        cursor: disabled ? "default" : "pointer",
        opacity: disabled ? 0.5 : 1,
        transition: "all 0.2s ease",
        "&:hover": disabled ? undefined : { borderColor: "primary.main" },
        ...sx,
      }}
      onClick={openFilePicker}
      onDragOver={(e) => {
        e.preventDefault();
        if (!disabled) setIsDragOver(true);
      }}
      onDragLeave={() => setIsDragOver(false)}
      onDrop={(e) => {
        e.preventDefault();
        setIsDragOver(false);
        if (!disabled && e.dataTransfer.files.length > 0) {
          onFilesSelected(e.dataTransfer.files);
        }
      }}
    >
      {children}
    </Paper>
  );
}
