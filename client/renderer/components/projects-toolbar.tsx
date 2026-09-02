"use client";
import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { Button, Toolbar, Tooltip } from "@mui/material";
import { Add, Upload, Restore } from "@mui/icons-material";
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
      {/* "Migrate legacy" is withdrawn for the alpha - see welcome-chooser.tsx.
          Import (which copies) covers the same need without touching the
          user's existing CCP4i2 installation. */}
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
