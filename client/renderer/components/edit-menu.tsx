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
import { Search, Settings } from "@mui/icons-material";
import { useFindInPage } from "../providers/find-in-page-provider";

export default function EditMenu() {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);
  const findInPage = useFindInPage();
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

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        Edit
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <MenuItem onClick={handleFind}>
          <ListItemIcon>
            <Search fontSize="small" />
          </ListItemIcon>
          <ListItemText>Find</ListItemText>
          <Typography variant="body2" color="textSecondary">
            ({typeof navigator !== "undefined" && navigator?.platform?.includes("Mac") ? "\u2318" : "Ctrl+"}F)
          </Typography>
        </MenuItem>
        <MenuItem onClick={handleClose}>
          <ListItemIcon>
            <Settings fontSize="small" />
          </ListItemIcon>
          <ListItemText>Preferences</ListItemText>
        </MenuItem>
      </Menu>
    </>
  );
}
