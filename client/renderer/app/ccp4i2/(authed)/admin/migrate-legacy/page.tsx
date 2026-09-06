"use client";
import { Alert, AlertTitle, Button, Paper, Stack, Typography } from "@mui/material";
import { Upload } from "@mui/icons-material";
import { useRouter } from "next/navigation";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";

/**
 * In-place migration is withdrawn for the alpha.
 *
 * The page is kept rather than deleted so that a bookmark, or a link from the
 * issue tracker, lands on an explanation instead of a 404 — and so that
 * re-enabling it later is a matter of restoring one import. The component it
 * used, `migrate-legacy-content.tsx`, is untouched; the server endpoints it
 * calls now return 403 unless CCP4I2_ALLOW_INPLACE_MIGRATION is set.
 */
export default function MigrateLegacyPage() {
  const router = useRouter();

  return (
    <Stack
      sx={{
        height: "100vh",
        "@supports (height: 100dvh)": { height: "100dvh" },
        overflow: "hidden",
      }}
    >
      <CCP4i2TopBar
        title="Migrate legacy CCP4i2"
        showBackButton
        backPath="/ccp4i2"
      />
      <Paper sx={{ flex: 1, minHeight: 0, overflow: "auto", p: 3 }}>
        <Stack spacing={3} sx={{ maxWidth: "48rem" }}>
          <Alert severity="info">
            <AlertTitle>Not available in this alpha</AlertTitle>
            Migrating a legacy CCP4i2 installation in place would adopt your
            existing projects where they sit, alongside the CCP4i2 you are still
            using day to day. Until we are confident that cannot damage them, the
            option is withdrawn.
          </Alert>

          <Typography variant="body1">
            Importing a <strong>copy</strong> does the same job safely and is
            fully supported: point it at a project directory or a project zip and
            it copies what it needs, leaving the original untouched. Import one
            project at a time, and your existing CCP4i2 keeps working exactly as
            it does now.
          </Typography>

          <Typography variant="body2" color="text.secondary">
            If you have a legacy database holding many projects and copying each
            one is impractical, please say so on the alpha feedback document —
            that is exactly the case we need to hear about before deciding how
            in-place migration should work.
          </Typography>

          <Stack direction="row" spacing={2}>
            <Button
              variant="contained"
              startIcon={<Upload />}
              onClick={() => router.push("/ccp4i2/import-project")}
            >
              Import a project
            </Button>
            <Button variant="outlined" onClick={() => router.push("/ccp4i2")}>
              Back to projects
            </Button>
          </Stack>
        </Stack>
      </Paper>
    </Stack>
  );
}
