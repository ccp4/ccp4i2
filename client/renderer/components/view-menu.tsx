"use client";
import { useState, useRef } from "react";
import {
  Button,
  ListItemText,
  Menu,
  MenuItem,
} from "@mui/material";
import {
  YoutubeSearchedFor,
  ZoomIn,
  ZoomOut,
  Refresh,
  SyncAlt,
} from "@mui/icons-material";
import { ThemeToggle } from "./theme-toggle";
import { DevModeToggle } from "./dev-mode-toggle";
import { useCCP4i2Window } from "../app-context";
import { useApi } from "../api";
import { Project } from "../types/models";
import { useRouter } from "next/dist/client/components/navigation";
import LanIcon from "@mui/icons-material/Lan";
import { CCP4i2MenuItem } from "./menu-item";

export default function ViewMenu() {
  const { projectId } = useCCP4i2Window();
  const api = useApi();
  const { data: project } = api.get<Project>(`projects/${projectId}`);
  const router = useRouter();

  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const themeToggleRef = useRef<HTMLButtonElement>(null);
  const open = Boolean(anchorEl);
  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
  };

  const handleReload = () => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("reload-window", {});
    } else {
      window.location.reload();
    }
    handleClose();
  };

  const handleForceReload = () => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("force-reload-window", {});
    } else {
      window.location.reload();
    }
    handleClose();
  };

  const handleZoomIn = () => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("zoom-in", {});
    }
    handleClose();
  };

  const handleZoomOut = () => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("zoom-out", {});
    }
    handleClose();
  };

  const handleZoomReset = () => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("zoom-reset", {});
    }
    handleClose();
  };

  const handleThemeToggle = () => {
    if (themeToggleRef.current) {
      themeToggleRef.current.click();
    }
    handleClose();
  };

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        View
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <CCP4i2MenuItem
          text="Reload"
          icon={Refresh}
          onClick={handleReload}
          shortcut="Ctrl+R"
        />
        <CCP4i2MenuItem
          text="Force Reload"
          icon={SyncAlt}
          onClick={handleForceReload}
          shortcut="Ctrl+Shift+R"
        />
        <CCP4i2MenuItem
          text="Zoom In"
          icon={ZoomIn}
          onClick={handleZoomIn}
          shortcut="Ctrl++"
        />
        <CCP4i2MenuItem
          text="Zoom Out"
          icon={ZoomOut}
          onClick={handleZoomOut}
          shortcut="Ctrl+-"
        />
        <CCP4i2MenuItem
          text="Reset Zoom"
          icon={YoutubeSearchedFor}
          onClick={handleZoomReset}
          shortcut="Ctrl+0"
        />
        <MenuItem onClick={handleThemeToggle}>
          <ThemeToggle ref={themeToggleRef} />
          <ListItemText>Toggle Theme</ListItemText>
        </MenuItem>
        <MenuItem>
          <DevModeToggle />
        </MenuItem>
        {project && (
          <CCP4i2MenuItem
            text="Project Network"
            icon={LanIcon}
            onClick={() => {
              router.push(`/ccp4i2/project/${projectId}/network`);
              handleClose();
            }}
          />
        )}
      </Menu>
    </>
  );
}
