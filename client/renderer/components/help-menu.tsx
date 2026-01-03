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
import { Help, PlayCircle } from "@mui/icons-material";

export default function HelpMenu() {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);
  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
  };

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        Help/Tutorials
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <MenuItem onClick={handleClose}>About</MenuItem>
        <MenuItem onClick={handleClose}>Tutorials</MenuItem>
        <MenuItem onClick={handleClose}>Quickstart - 10 minute intro</MenuItem>
        <MenuItem onClick={handleClose}>
          Quick expert - more quick hints
        </MenuItem>
        <MenuItem onClick={handleClose}>
          <ListItemIcon>
            <PlayCircle fontSize="small" />
          </ListItemIcon>
          <ListItemText>View YouTube video</ListItemText>
        </MenuItem>
        <MenuItem>Task documentation</MenuItem>
        <MenuItem>Task info from CCP4 Cloud</MenuItem>
        <MenuItem>
          <ListItemIcon>
            <Help fontSize="small" />
          </ListItemIcon>
          <ListItemText>CCP4i2 Help</ListItemText>
          <Typography variant="body2" color="textSecondary">
            (F1)
          </Typography>
        </MenuItem>
        <MenuItem>Tips of the day</MenuItem>
      </Menu>
    </>
  );
}
