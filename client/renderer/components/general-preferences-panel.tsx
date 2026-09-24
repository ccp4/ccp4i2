"use client";
/**
 * Preferences → General: running-app preferences that are neither program
 * locations nor credentials.
 *
 * Currently just the import-provenance prompt toggle. This preference used to
 * live in the View menu, but the View menu is for per-view display toggles;
 * an "ask me / don't ask me" preference belongs with the other settings in
 * Edit → Preferences (Paul, issue #416).
 */
import { FormControlLabel, Stack, Switch, Typography } from "@mui/material";
import { useUiPreference } from "../lib/ui-preferences";

export function GeneralPreferencesPanel() {
  const [captureImportProvenance, setCaptureImportProvenance] = useUiPreference(
    "captureImportProvenance"
  );

  return (
    <Stack spacing={3} sx={{ maxWidth: 720, mx: "auto", p: 2 }}>
      <Typography variant="h6">General</Typography>
      <Stack spacing={0.5}>
        <FormControlLabel
          control={
            <Switch
              checked={captureImportProvenance}
              onChange={(e) => setCaptureImportProvenance(e.target.checked)}
            />
          }
          label="Ask for a note when importing a file"
        />
        <Typography variant="body2" color="text.secondary" sx={{ pl: 6 }}>
          When on, importing a file offers a short prompt to record where it
          came from and why. The note is stored with the import for provenance.
        </Typography>
      </Stack>
    </Stack>
  );
}
