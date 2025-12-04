"use client";
import { useRouter } from "next/navigation";
import { Button, Toolbar, Tooltip } from "@mui/material";
import { Add, Upload } from "@mui/icons-material";

export default function ProjectsToolbar() {
  const router = useRouter();

  function importProjects() {
    router.push("/import-project");
  }

  return (
    <Toolbar sx={{ gap: 2 }}>
      <Tooltip title="Start a new project">
        <Button
          variant="contained"
          startIcon={<Add />}
          onClick={() => router.push("/new-project")}
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
    </Toolbar>
  );
}
