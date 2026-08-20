"use client";
import { useState } from "react";
import {
  Button,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Typography,
} from "@mui/material";
import { Search, Settings, ContentCut, ContentCopy, ContentPaste, SelectAll } from "@mui/icons-material";
import { useRouter } from "next/navigation";
import { useFindInPage } from "../providers/find-in-page-provider";

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

  const handleCut = () => {
    document.execCommand("cut");
    handleClose();
  };

  const handleCopy = () => {
    document.execCommand("copy");
    handleClose();
  };

  const handlePaste = () => {
    document.execCommand("paste");
    handleClose();
  };

  const handleSelectAll = () => {
    document.execCommand("selectAll");
    handleClose();
  };

  const handlePreferences = () => {
    handleClose();
    router.push("/ccp4i2/preferences");
  };

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        Edit
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <MenuItem onClick={handleCut}>
          <ListItemIcon>
            <ContentCut fontSize="small" />
          </ListItemIcon>
          <ListItemText>Cut</ListItemText>
          <Typography variant="body2" color="textSecondary">
            (Ctrl+X)
          </Typography>
        </MenuItem>
        <MenuItem onClick={handleCopy}>
          <ListItemIcon>
            <ContentCopy fontSize="small" />
          </ListItemIcon>
          <ListItemText>Copy</ListItemText>
          <Typography variant="body2" color="textSecondary">
            (Ctrl+C)
          </Typography>
        </MenuItem>
        <MenuItem onClick={handlePaste}>
          <ListItemIcon>
            <ContentPaste fontSize="small" />
          </ListItemIcon>
          <ListItemText>Paste</ListItemText>
          <Typography variant="body2" color="textSecondary">
            (Ctrl+V)
          </Typography>
        </MenuItem>
        <MenuItem onClick={handleSelectAll}>
          <ListItemIcon>
            <SelectAll fontSize="small" />
          </ListItemIcon>
          <ListItemText>Select All</ListItemText>
          <Typography variant="body2" color="textSecondary">
            (Ctrl+A)
          </Typography>
        </MenuItem>
        <MenuItem onClick={handleFind}>
          <ListItemIcon>
            <Search fontSize="small" />
          </ListItemIcon>
          <ListItemText>Find</ListItemText>
          <Typography variant="body2" color="textSecondary">
            ({typeof navigator !== "undefined" && navigator?.platform?.includes("Mac") ? "\u2318" : "Ctrl+"}F)
          </Typography>
        </MenuItem>
        <MenuItem onClick={handlePreferences}>
          <ListItemIcon>
            <Settings fontSize="small" />
          </ListItemIcon>
          <ListItemText>Preferences</ListItemText>
        </MenuItem>
      </Menu>
    </>
  );
}
