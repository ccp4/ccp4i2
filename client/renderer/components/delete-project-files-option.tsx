"use client";
import { useState } from "react";
import { Box, Checkbox, FormControlLabel, Typography } from "@mui/material";

interface DeleteProjectFilesOptionProps {
  /** Called with the current choice whenever it changes. Starts false. */
  onChange: (deleteFiles: boolean) => void;
  /** How many projects the dialog is about; only affects the wording. */
  count?: number;
}

/**
 * The "also delete the files on disk" choice in the delete-project dialog.
 *
 * Qt-era CCP4i2 offered this checkbox and defaulted to leaving the files;
 * so does this. It keeps its own state so it can live inside the shared
 * delete dialog, whose content is handed over as ready-made elements; the
 * caller reads the choice back through onChange when Delete is pressed.
 */
export function DeleteProjectFilesOption({
  onChange,
  count = 1,
}: DeleteProjectFilesOptionProps) {
  const [checked, setChecked] = useState(false);
  const plural = count !== 1;
  return (
    <Box sx={{ mt: 2 }}>
      <FormControlLabel
        control={
          <Checkbox
            checked={checked}
            onChange={(event) => {
              setChecked(event.target.checked);
              onChange(event.target.checked);
            }}
          />
        }
        label={
          plural
            ? "Also delete the projects' files on disk"
            : "Also delete the project's files on disk"
        }
      />
      <Typography variant="body2" color="text.secondary">
        {checked
          ? plural
            ? "The project folders and everything in them will be removed."
            : "The project folder and everything in it will be removed."
          : plural
          ? "The project folders are left where they are; a project can be brought back with Import → a project folder."
          : "The project folder is left where it is; the project can be brought back with Import → a project folder."}
      </Typography>
    </Box>
  );
}
