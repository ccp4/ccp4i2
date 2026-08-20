"use client";
import { useState } from "react";
import { Button, Menu } from "@mui/material";
import {
  Search,
  Settings,
  ContentCut,
  ContentCopy,
  ContentPaste,
  SelectAll,
} from "@mui/icons-material";
import { useRouter } from "next/navigation";
import { useFindInPage } from "../providers/find-in-page-provider";
import { CCP4i2MenuItem } from "./menu-item";

export default function EditMenu() {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);
  const findInPage = useFindInPage();
  const router = useRouter();
  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
  };

  const handleFind = () => {
    handleClose();
    findInPage?.open();
  };

  // In Electron use webContents native methods via IPC; fall back to execCommand in web mode.
  const execEdit = (ipcChannel: string, legacyCommand: string) => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage(ipcChannel, {});
    } else {
      // execCommand is deprecated but has no viable alternative for cut/paste in plain web
      (document as any).execCommand(legacyCommand);
    }
    handleClose();
  };

  const handleCut = () => execEdit("edit-cut", "cut");
  const handleCopy = () => execEdit("edit-copy", "copy");
  const handlePaste = () => execEdit("edit-paste", "paste");
  const handleSelectAll = () => execEdit("edit-select-all", "selectAll");

  const handlePreferences = () => {
    handleClose();
    router.push("/ccp4i2/preferences");
  };

  let ctrlOrCmd = "Ctrl+";
  if (
    typeof navigator !== "undefined" &&
    navigator?.platform?.includes("Mac")
  ) {
    ctrlOrCmd = "\u2318"; // Command key symbol for Mac
  }

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        Edit
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <CCP4i2MenuItem
          text="Cut"
          icon={ContentCut}
          onClick={handleCut}
          shortcut={`${ctrlOrCmd}X`}
        />
        <CCP4i2MenuItem
          text="Copy"
          icon={ContentCopy}
          onClick={handleCopy}
          shortcut={`${ctrlOrCmd}C`}
        />
        <CCP4i2MenuItem
          text="Paste"
          icon={ContentPaste}
          onClick={handlePaste}
          shortcut={`${ctrlOrCmd}V`}
        />
        <CCP4i2MenuItem
          text="Select All"
          icon={SelectAll}
          onClick={handleSelectAll}
          shortcut={`${ctrlOrCmd}A`}
        />
        <CCP4i2MenuItem
          text="Find"
          icon={Search}
          onClick={handleFind}
          shortcut={`${ctrlOrCmd}F`}
        />
        <CCP4i2MenuItem
          text="Preferences"
          icon={Settings}
          onClick={handlePreferences}
        />
      </Menu>
    </>
  );
}
