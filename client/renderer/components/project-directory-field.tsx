"use client";
import { Button, Stack, TextField, ToggleButton, ToggleButtonGroup } from "@mui/material";
import { Folder } from "@mui/icons-material";
import type { UseProjectDirectoryReturn } from "../hooks/use-project-directory";

interface ProjectDirectoryFieldProps extends UseProjectDirectoryReturn {
  mode: "create" | "move";
  /** Required in move mode: the project's current directory path. */
  currentDirectory?: string;
}

/**
 * Directory picker UI shared between the new-project and edit-project forms.
 *
 * Pair with `useProjectDirectory` in the parent component — pass the hook's
 * return value (spread) along with `mode` and, in move mode, `currentDirectory`.
 *
 * Create mode: Default/Custom toggle → optional parent-directory picker →
 *   computed directory preview.
 * Move mode: current directory (read-only) → optional new-location picker
 *   with Browse / Change / Clear buttons.
 */
export function ProjectDirectoryField({
  mode,
  currentDirectory,
  electronAPIAvailable,
  parentDirectory,
  computedDirectory,
  defaultDirectory,
  customMode,
  setCustomMode,
  moveMode,
  setMoveMode,
  directoryError,
  handleBrowse,
  clearParent,
}: ProjectDirectoryFieldProps) {
  // True when the project already lives at its computed default location —
  // used to suppress the "Move to Default" option in that case.
  const isAtDefault =
    !!defaultDirectory && !!currentDirectory && currentDirectory === defaultDirectory;

  return (
    <Stack spacing={2}>
      {/* ── Create mode: Default / Custom toggle ── */}
      {mode === "create" && (
        <ToggleButtonGroup
          exclusive
          value={customMode}
          onChange={(_e, value) => {
            if (value !== null) setCustomMode(value);
          }}
          fullWidth
        >
          <ToggleButton value={false}>Default Directory</ToggleButton>
          <ToggleButton value={true}>Custom Directory</ToggleButton>
        </ToggleButtonGroup>
      )}

      {/* ── Move mode: adaptive toggle ──
            Already at default → 2 options (Keep Directory / Move to Custom)
            At a custom path   → 3 options (Keep Directory / Move to Default / Move to Custom) */}
      {mode === "move" && isAtDefault && (
        <ToggleButtonGroup
          exclusive
          value={moveMode === "custom" ? "custom" : "keep"}
          onChange={(_e, value) => {
            if (value !== null) setMoveMode(value === "custom" ? "custom" : "keep");
          }}
          fullWidth
        >
          <ToggleButton value="keep">Keep Directory</ToggleButton>
          <ToggleButton value="custom">Move to Custom</ToggleButton>
        </ToggleButtonGroup>
      )}
      {mode === "move" && !isAtDefault && (
        <ToggleButtonGroup
          exclusive
          value={moveMode}
          onChange={(_e, value) => {
            if (value !== null) setMoveMode(value);
          }}
          fullWidth
        >
          <ToggleButton value="keep">Keep Directory</ToggleButton>
          <ToggleButton value="default">Move to Default</ToggleButton>
          <ToggleButton value="custom">Move to Custom</ToggleButton>
        </ToggleButtonGroup>
      )}

      {/* ── Move mode: current directory (read-only) ── */}
      {mode === "move" && currentDirectory && (
        <TextField label="Current directory" value={currentDirectory} disabled />
      )}

      {/* ── Move mode: "Move to Default" — show target path and any error ── */}
      {mode === "move" && moveMode === "default" && (
        <TextField
          label="New project directory"
          value={defaultDirectory}
          disabled
          error={directoryError.length > 0}
          helperText={directoryError || " "}
        />
      )}

      {/* ── Parent-directory picker row ──
            Both modes show the parent directory (not the full computed path) so the
            layout is consistent: pick a parent here, see the full path below. */}
      {electronAPIAvailable && (mode === "create" ? customMode : moveMode === "custom") && (
        <Stack direction="row" spacing={2} sx={{ alignItems: "center" }}>
          <TextField
            label="Parent directory"
            value={parentDirectory}
            disabled
            sx={{ flexGrow: 1 }}
            required
          />
          <Button variant="outlined" startIcon={<Folder />} onClick={handleBrowse}>
            {parentDirectory ? "Change" : "Select"}
          </Button>
          {mode === "move" && parentDirectory && (
            <Button variant="text" onClick={clearParent}>
              Clear
            </Button>
          )}
        </Stack>
      )}

      {/* ── Resulting directory preview ──
            Create mode: always shown (acts as preview even in default mode).
            Move+custom mode: always shown once the section is active so the
            "Directory is required" / "already exists" errors have a place to appear. */}
      {(mode === "create" || (mode === "move" && moveMode === "custom")) && (
        <Stack direction="row">
          <TextField
            label="Resulting project directory"
            value={mode === "create" ? computedDirectory : (parentDirectory ? computedDirectory : "")}
            disabled
            error={directoryError.length > 0}
            helperText={directoryError}
            sx={{ flexGrow: 1 }}
          />
        </Stack>
      )}

    </Stack>
  );
}
