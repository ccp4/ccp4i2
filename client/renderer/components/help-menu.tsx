"use client";
import { useEffect, useState } from "react";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Link,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Stack,
  Typography,
} from "@mui/material";
import { Help, Info, MenuBook, EmojiObjects } from "@mui/icons-material";
import { TipOfTheDayDialog } from "./tip-of-the-day-dialog";

const CCP4I2_HELP_URL = "https://www.ccp4.ac.uk/html/ccp4i2.html";

interface VersionInfo {
  web?: { buildTimestamp?: string; gitCommit?: string };
}

export default function HelpMenu() {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);
  const [aboutOpen, setAboutOpen] = useState(false);
  const [tipOpen, setTipOpen] = useState(false);
  const [version, setVersion] = useState<VersionInfo | null>(null);

  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
  };

  // Fetch build info once the About dialog is first opened.
  useEffect(() => {
    if (aboutOpen && !version) {
      fetch("/api/version")
        .then((r) => (r.ok ? r.json() : null))
        .then((data) => data && setVersion(data))
        .catch(() => {
          /* build info is best-effort */
        });
    }
  }, [aboutOpen, version]);

  const handleAbout = () => {
    handleClose();
    setAboutOpen(true);
  };

  const handleHelp = () => {
    handleClose();
    window.open(CCP4I2_HELP_URL, "_blank", "noopener,noreferrer");
  };

  const handleTip = () => {
    handleClose();
    setTipOpen(true);
  };

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        Help/Tutorials
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <MenuItem onClick={handleHelp}>
          <ListItemIcon>
            <MenuBook fontSize="small" />
          </ListItemIcon>
          <ListItemText>CCP4i2 documentation</ListItemText>
          <Typography variant="body2" color="textSecondary">
            (F1)
          </Typography>
        </MenuItem>
        <MenuItem onClick={handleTip}>
          <ListItemIcon>
            <EmojiObjects fontSize="small" />
          </ListItemIcon>
          <ListItemText>Tip of the day</ListItemText>
        </MenuItem>
        <MenuItem onClick={handleAbout}>
          <ListItemIcon>
            <Info fontSize="small" />
          </ListItemIcon>
          <ListItemText>About CCP4i2</ListItemText>
        </MenuItem>
      </Menu>

      <TipOfTheDayDialog open={tipOpen} onClose={() => setTipOpen(false)} />

      <Dialog
        open={aboutOpen}
        onClose={() => setAboutOpen(false)}
        maxWidth="xs"
        fullWidth
      >
        <DialogTitle sx={{ display: "flex", alignItems: "center", gap: 1 }}>
          <Help color="primary" />
          About CCP4i2
        </DialogTitle>
        <DialogContent>
          <Stack spacing={1.5} sx={{ pt: 1 }}>
            <Typography variant="body1">
              CCP4i2 &mdash; a graphical environment for crystallographic
              computing.
            </Typography>
            {version?.web?.buildTimestamp &&
              version.web.buildTimestamp !== "dev" && (
                <Typography variant="body2" color="text.secondary">
                  <strong>Build:</strong> {version.web.buildTimestamp}
                </Typography>
              )}
            {version?.web?.gitCommit &&
              version.web.gitCommit !== "unknown" && (
                <Typography variant="body2" color="text.secondary">
                  <strong>Commit:</strong> {version.web.gitCommit}
                </Typography>
              )}
            <Typography variant="body2" color="text.secondary">
              <Link
                href="https://www.ccp4.ac.uk"
                target="_blank"
                rel="noopener noreferrer"
              >
                www.ccp4.ac.uk
              </Link>
            </Typography>
          </Stack>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setAboutOpen(false)}>Close</Button>
        </DialogActions>
      </Dialog>
    </>
  );
}
