"use client";
import { useCallback, useState } from "react";
import { Button, Menu, MenuItem } from "@mui/material";
import { useRunningProcesses } from "../providers/running-processes";

export default function UtilMenu() {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const { setJobsAndProcessesDialogOpen } = useRunningProcesses();
  const open = Boolean(anchorEl);
  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
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
        <MenuItem onClick={handleRunningJobsClick}>
          Running jobs and processes
        </MenuItem>
      </Menu>
    </>
  );
}
