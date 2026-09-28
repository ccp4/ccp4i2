import React, { useMemo, useState } from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  DialogContentText,
  Button,
  Box,
} from "@mui/material";
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import CheckIcon from "@mui/icons-material/Check";
import CloseIcon from "@mui/icons-material/Close";

interface I2RunDialogProps {
  open: boolean;
  /**
   * The command as the server rendered it: a complete, runnable invocation
   * (``ccp4-python -m ccp4i2.cli.i2run <task> ...``).
   */
  command: string;
  /**
   * Where the command has to be run from, or null when it runs from anywhere.
   * A dev checkout needs its own ``server/`` directory on sys.path ahead of
   * the legacy ``ccp4i2`` that the CCP4 bundle ships.
   */
  workingDirectory?: string | null;
  /**
   * Variables the server is running with that a fresh terminal will not have.
   * These are not decoration: a run without CCP4I2_PROJECTS_DIR or
   * CCP4I2_HOME can address a different database and give no sign of it.
   */
  environment?: Record<string, string> | null;
  /** CCP4 setup script to source first; null on Windows, which has none. */
  ccp4Setup?: string | null;
  /** The server's sys.platform, which decides the shell syntax. */
  platform?: string | null;
  onClose: () => void;
}

/**
 * Single-quote for POSIX shells: inside single quotes everything is literal,
 * so the only thing to handle is a single quote itself.
 */
const posixQuote = (value: string) => `'${value.split("'").join(`'\\''`)}'`;

/** Single-quote for PowerShell, where a literal quote is doubled. */
const powershellQuote = (value: string) =>
  `'${value.split("'").join("''")}'`;

const isWindows = (platform?: string | null) =>
  (platform ?? "").startsWith("win");

export const I2RunDialog: React.FC<I2RunDialogProps> = ({
  open,
  command,
  workingDirectory,
  environment,
  ccp4Setup,
  platform,
  onClose,
}) => {
  const [copied, setCopied] = useState(false);

  // The whole recipe, in the order it has to be run, quoted for the shell the
  // server says it is on. Copy takes all of it, so pasting works rather than
  // working only for the last line.
  //
  // This used to be assembled here from Electron's cwd and config, as
  // `${cwd}/ccp4i2/i2run/i2run.sh <args> --dbFile <...>/db.sqlite3` -- a
  // script that does not exist and a flag i2run does not accept. In a packaged
  // build cwd is "" (there is no server directory on disk), so the dialog
  // rendered an empty string. The server knows which invocation works and
  // which variables it is itself running with, so it now says.
  const script = useMemo(() => {
    if (!command) return "";
    const windows = isWindows(platform);
    const quote = windows ? powershellQuote : posixQuote;
    const lines: string[] = [];

    if (ccp4Setup) {
      lines.push(`source ${posixQuote(ccp4Setup)}`);
    } else if (windows) {
      lines.push("# Run this from the CCP4 command prompt.");
    }

    for (const [name, value] of Object.entries(environment ?? {})) {
      lines.push(
        windows
          ? `$env:${name} = ${quote(value)}`
          : `export ${name}=${quote(value)}`
      );
    }

    if (workingDirectory) lines.push(`cd ${quote(workingDirectory)}`);
    lines.push(command);
    return lines.join("\n");
  }, [command, workingDirectory, environment, ccp4Setup, platform]);

  const handleCopy = async () => {
    if (!script) return;
    try {
      await navigator.clipboard.writeText(script);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch {
      // Clipboard access can be refused; the command is on screen to select.
      setCopied(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle>i2run command</DialogTitle>
      <DialogContent>
        <DialogContentText sx={{ mb: 2 }}>
          This runs the job as configured. The lines above the command put the
          terminal in the same environment as the app — without them it can
          address a different database.
        </DialogContentText>
        <Box
          component="pre"
          sx={{
            m: 0,
            p: 2,
            borderRadius: 1,
            bgcolor: "action.hover",
            fontFamily: "monospace",
            fontSize: "0.85rem",
            whiteSpace: "pre-wrap",
            wordBreak: "break-all",
            userSelect: "text",
          }}
        >
          {script}
        </Box>
      </DialogContent>
      <DialogActions>
        <Button
          onClick={handleCopy}
          variant="outlined"
          disabled={!script}
          startIcon={copied ? <CheckIcon /> : <ContentCopyIcon />}
        >
          {copied ? "Copied" : "Copy"}
        </Button>
        <Button onClick={onClose} variant="contained" startIcon={<CloseIcon />}>
          Dismiss
        </Button>
      </DialogActions>
    </Dialog>
  );
};
