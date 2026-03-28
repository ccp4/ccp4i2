/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
"use client";
import { useCallback, useState } from "react";
import { Button, Menu, MenuItem } from "@mui/material";
import { useRouter } from "next/navigation";
import { useRunningProcesses } from "../providers/running-processes";

export default function UtilMenu() {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const { jobsAndProcessesDialogOpen, setJobsAndProcessesDialogOpen } =
    useRunningProcesses();
  const router = useRouter();
  const open = Boolean(anchorEl);
  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
  };
  const handleSystemAdministratorToolsClick = () => {
    router.push("/ccp4i2/config");
    handleClose();
  };
  const handleRunningJobsClick = useCallback(async () => {
    setJobsAndProcessesDialogOpen(true);
    handleClose();
  }, [setJobsAndProcessesDialogOpen]);

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        Utilities
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <MenuItem onClick={handleClose}>Copy demo data to project</MenuItem>
        <MenuItem onClick={handleRunningJobsClick}>
          Running jobs and processes {`${jobsAndProcessesDialogOpen}`}
        </MenuItem>
        <MenuItem onClick={handleClose}>Manage imported files</MenuItem>
        <MenuItem onClick={handleClose}>Send error report</MenuItem>
        <MenuItem onClick={handleSystemAdministratorToolsClick}>
          System administrator tools
        </MenuItem>
        <MenuItem onClick={handleClose}>Developer tools</MenuItem>
      </Menu>
    </>
  );
}
