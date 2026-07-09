"use client";
import { useRouter } from "next/navigation";
import { Button, Toolbar, Tooltip } from "@mui/material";
import { Add, Upload, SettingsBackupRestore } from "@mui/icons-material";
import { useAuth } from "@/lib/compounds/auth-context";

export default function ProjectsToolbar() {
  const router = useRouter();
  const { canAdminister } = useAuth();

  function importProjects() {
    router.push("/ccp4i2/import-project");
  }

  return (
    <Toolbar sx={{ gap: 2 }}>
      <Tooltip title="Start a new project">
        <Button
          variant="contained"
          startIcon={<Add />}
          onClick={() => router.push("/ccp4i2/new-project")}
        >
          New
        </Button>
      </Tooltip>
      <Tooltip title="Import existing projects">
        <Button
          variant="outlined"
          startIcon={<Upload />}
          onClick={importProjects}
        >
          Import
        </Button>
      </Tooltip>
      {canAdminister && (
        <Tooltip title="Migrate a legacy CCP4i2 database into this installation">
          <Button
            variant="outlined"
            startIcon={<SettingsBackupRestore />}
            onClick={() => router.push("/ccp4i2/admin/migrate-legacy")}
          >
            Migrate legacy
          </Button>
        </Tooltip>
      )}
    </Toolbar>
  );
}
