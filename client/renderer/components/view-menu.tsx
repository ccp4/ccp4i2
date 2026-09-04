"use client";
import { useState } from "react";
import { Button, Menu } from "@mui/material";
import {
  Brightness4,
  Brightness7,
  DeveloperMode,
  Refresh,
  SyncAlt,
  Visibility,
  VisibilityOff,
  YoutubeSearchedFor,
  ZoomIn,
  ZoomOut,
} from "@mui/icons-material";
import { useTheme } from "../theme/theme-provider";
import { useUiPreference } from "../lib/ui-preferences";
import { useCCP4i2Window } from "../app-context";
import { useApi } from "../api";
import { Project } from "../types/models";
import { useRouter } from "next/dist/client/components/navigation";
import LanIcon from "@mui/icons-material/Lan";
import { CCP4i2MenuItem } from "./menu-item";

export default function ViewMenu() {
  const { projectId, devMode, setDevMode } = useCCP4i2Window();
  const { mode, setTheme } = useTheme();
  const [showJobIcons, setShowJobIcons] = useUiPreference("showJobIcons");
  const api = useApi();
  const { data: project } = api.get<Project>(`projects/${projectId}`);
  const router = useRouter();

  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
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

  const handleSwitchTheme = () => {
    const next = mode === "light" ? "dark" : "light";
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("set-theme-mode", { theme: next });
    }
    setTheme(next);
    handleClose();
  };
  const handleSwitchDevMode = () => {
    const enabled = !devMode;
    setDevMode(enabled);
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("set-dev-mode", { enabled });
    }
    handleClose();
  };
  const handleSwitchJobIcons = () => {
    setShowJobIcons(!showJobIcons);
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
          secondary="Ctrl+R"
        />
        <CCP4i2MenuItem
          text="Force Reload"
          icon={SyncAlt}
          onClick={handleForceReload}
          secondary="Ctrl+Shift+R"
        />
        <CCP4i2MenuItem
          text="Zoom In"
          icon={ZoomIn}
          onClick={handleZoomIn}
          secondary="Ctrl++"
        />
        <CCP4i2MenuItem
          text="Zoom Out"
          icon={ZoomOut}
          onClick={handleZoomOut}
          secondary="Ctrl+-"
        />
        <CCP4i2MenuItem
          text="Reset Zoom"
          icon={YoutubeSearchedFor}
          onClick={handleZoomReset}
          secondary="Ctrl+0"
        />
        {/* Plain items whose text says what a click will do, rather than
            controls embedded in a menu row (where the hover of the inner
            switch made it unclear what to click -- Paul). */}
        <CCP4i2MenuItem
          text={mode === "light" ? "Switch to Dark Mode" : "Switch to Light Mode"}
          icon={mode === "light" ? Brightness4 : Brightness7}
          onClick={handleSwitchTheme}
        />
        <CCP4i2MenuItem
          text={devMode ? "Turn Dev Mode Off" : "Turn Dev Mode On"}
          icon={DeveloperMode}
          onClick={handleSwitchDevMode}
        />
        <CCP4i2MenuItem
          text={showJobIcons ? "Hide Job Icons" : "Show Job Icons"}
          icon={showJobIcons ? VisibilityOff : Visibility}
          onClick={handleSwitchJobIcons}
        />
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
