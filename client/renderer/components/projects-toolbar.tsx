"use client";
import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { Button, Toolbar, Tooltip } from "@mui/material";
import { Add, Upload, SettingsBackupRestore, Restore } from "@mui/icons-material";
import { useAuth } from "@/lib/compounds/auth-context";
import { isElectron } from "../utils/platform";
import { RecoverProjectsDialog } from "./recover-projects-dialog";

export default function ProjectsToolbar() {
  const router = useRouter();
  const { canAdminister } = useAuth();
  const [recoverOpen, setRecoverOpen] = useState(false);
  // Recovery reads project directories on the machine running the server and
  // needs a native folder picker when the registry is gone, so it is
  // desktop-only. Resolved in an effect so the first client render matches the
  // server-rendered markup.
  const [canRecover, setCanRecover] = useState(false);
  useEffect(() => setCanRecover(isElectron()), []);

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
      {canRecover && (
        <Tooltip title="Rebuild the database from the records kept in your project folders">
          <Button
            variant="outlined"
            startIcon={<Restore />}
            onClick={() => setRecoverOpen(true)}
          >
            Recover
          </Button>
        </Tooltip>
      )}
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
      {canRecover && (
        <RecoverProjectsDialog
          open={recoverOpen}
          onClose={() => setRecoverOpen(false)}
          onRecovered={() => router.refresh()}
        />
      )}
    </Toolbar>
  );
}
