"use client";
import { useCallback, useEffect, useState } from "react";
import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Divider,
  IconButton,
  Link,
  Menu,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  Help,
  Info,
  MenuBook,
  EmojiObjects,
  ContentCopy,
  Check,
  SystemUpdateAlt,
} from "@mui/icons-material";
import { hasLocalSessionToken } from "@ccp4/ccp4i2-api";
import { TipOfTheDayDialog } from "./tip-of-the-day-dialog";
import { CCP4i2MenuItem } from "./menu-item";

const CCP4I2_HELP_URL = "https://www.ccp4.ac.uk/html/ccp4i2.html";

interface VersionInfo {
  web?: { buildTimestamp?: string; gitCommit?: string };
}

/**
 * What the About dialog can tell a user about the running backend.
 *
 * Both ports are chosen at launch (detect-port from 3000 upwards), so they
 * differ between runs and between machines: a tester asked to "curl the
 * API" has no way to guess them. The token is the per-launch local-session
 * secret the Electron preload exposes on window; it is what the Django
 * LocalSessionAuthMiddleware checks on every request, so anything outside
 * the app -- curl, a notebook, a script -- needs it to get past 401.
 */
interface SessionInfo {
  uvicornPort?: number;
  nextPort?: number;
  token?: string;
  userEmail?: string;
}

/**
 * A read-only value with a copy button. The value goes in a TextField
 * rather than a Typography so it is selectable and horizontally scrollable
 * -- a 64-character token would otherwise wrap into an unreadable block.
 */
function CopyableField({
  label,
  value,
  monospace = true,
}: {
  label: string;
  value: string;
  monospace?: boolean;
}) {
  const [copied, setCopied] = useState(false);

  const handleCopy = useCallback(() => {
    navigator.clipboard
      .writeText(value)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 1500);
      })
      .catch(() => {
        /* clipboard can be denied; the field is still selectable */
      });
  }, [value]);

  return (
    <Box sx={{ display: "flex", alignItems: "flex-end", gap: 0.5 }}>
      <TextField
        label={label}
        value={value}
        size="small"
        fullWidth
        variant="outlined"
        InputProps={{
          readOnly: true,
          sx: monospace
            ? { fontFamily: "monospace", fontSize: "0.8rem" }
            : undefined,
        }}
        onFocus={(event) => event.target.select()}
      />
      <Tooltip title={copied ? "Copied" : `Copy ${label.toLowerCase()}`}>
        <IconButton onClick={handleCopy} size="small" sx={{ mb: 0.25 }}>
          {copied ? <Check fontSize="small" color="success" /> : <ContentCopy fontSize="small" />}
        </IconButton>
      </Tooltip>
    </Box>
  );
}

export default function HelpMenu() {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);
  const [aboutOpen, setAboutOpen] = useState(false);
  const [tipOpen, setTipOpen] = useState(false);
  const [version, setVersion] = useState<VersionInfo | null>(null);
  const [session, setSession] = useState<SessionInfo | null>(null);
  // Only the packaged Electron app can self-update; hide the item in the web build.
  const [isElectron, setIsElectron] = useState(false);
  useEffect(() => {
    setIsElectron(typeof window !== "undefined" && !!window.electronAPI);
  }, []);

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

  // Collect the backend session details when the dialog first opens.
  //
  // The token comes straight off the preload-exposed window surface; the
  // ports come from the main process via the existing "get-config" IPC,
  // which already reports UVICORN_PORT and NEXT_PORT. In the web build
  // there is no electronAPI and no local session, so this stays null and
  // the section is simply not rendered.
  useEffect(() => {
    if (!aboutOpen || session) return;
    if (typeof window === "undefined") return;

    const local = hasLocalSessionToken() ? window.ccp4i2LocalSession : undefined;

    if (!window.electronAPI) {
      if (local) setSession({ token: local.token, userEmail: local.userEmail });
      return;
    }

    const handler = (_event: any, data: any) => {
      if (data?.message !== "get-config") return;
      setSession({
        uvicornPort: data.config?.UVICORN_PORT,
        nextPort: data.config?.NEXT_PORT,
        token: local?.token,
        userEmail: local?.userEmail,
      });
    };

    window.electronAPI.onMessage("message-from-main", handler);
    window.electronAPI.sendMessage("get-config");
    return () => {
      window.electronAPI.removeMessageListener("message-from-main", handler);
    };
  }, [aboutOpen, session]);

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

  // Fire-and-forget: the main process (ccp4i2-updater) runs the check and shows
  // the result — "up to date", "downloading", or "not available here" — in a
  // native dialog, so there is nothing to plumb back into the renderer.
  const handleCheckForUpdates = () => {
    handleClose();
    window.electronAPI?.sendMessage("check-for-updates");
  };

  // A command that works as pasted: the backend port and the bearer token
  // this launch actually uses, against an endpoint that needs auth (so a
  // 200 proves the token was accepted, which /health would not).
  const curlCommand =
    session?.uvicornPort && session?.token
      ? `curl -H "Authorization: Bearer ${session.token}" \\\n  http://localhost:${session.uvicornPort}/api/ccp4i2/projects/`
      : null;

  const showSession = Boolean(
    session && (session.uvicornPort || session.token)
  );

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        Help
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <CCP4i2MenuItem
          text="CCP4i2 documentation"
          icon={MenuBook}
          onClick={handleHelp}
          secondary="F1"
        />
        <CCP4i2MenuItem
          text="Tip of the day"
          icon={EmojiObjects}
          onClick={handleTip}
        />
        {isElectron && (
          <CCP4i2MenuItem
            text="Check for Updates…"
            icon={SystemUpdateAlt}
            onClick={handleCheckForUpdates}
          />
        )}
        <CCP4i2MenuItem text="About CCP4i2" icon={Info} onClick={handleAbout} />
      </Menu>

      <TipOfTheDayDialog open={tipOpen} onClose={() => setTipOpen(false)} />

      <Dialog
        open={aboutOpen}
        onClose={() => setAboutOpen(false)}
        maxWidth={showSession ? "sm" : "xs"}
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
            {version?.web?.gitCommit && version.web.gitCommit !== "unknown" && (
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

            {showSession && (
              <>
                <Divider sx={{ mt: 1 }} />
                <Typography variant="subtitle2">This session</Typography>
                <Typography variant="body2" color="text.secondary">
                  Ports and the access token are chosen afresh each time the
                  app starts, and are gone when it quits. Copy them here to
                  reach the backend from a terminal, a script or a notebook.
                </Typography>

                {session?.uvicornPort && (
                  <CopyableField
                    label="Backend (Django) port"
                    value={String(session.uvicornPort)}
                  />
                )}
                {session?.nextPort && (
                  <CopyableField
                    label="Front-end (Next.js) port"
                    value={String(session.nextPort)}
                  />
                )}
                {session?.userEmail && (
                  <CopyableField
                    label="Session user"
                    value={session.userEmail}
                    monospace={false}
                  />
                )}
                {session?.token && (
                  <CopyableField
                    label="Session token"
                    value={session.token}
                  />
                )}
                {curlCommand && (
                  <CopyableField
                    label="Example request"
                    value={curlCommand}
                  />
                )}
                {session?.token && (
                  <Typography variant="caption" color="warning.main">
                    The token grants full access to your projects. Do not paste
                    it into a bug report, a shared document or a chat.
                  </Typography>
                )}
              </>
            )}
          </Stack>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setAboutOpen(false)}>Close</Button>
        </DialogActions>
      </Dialog>
    </>
  );
}
