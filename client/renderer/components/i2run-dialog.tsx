import React, { useCallback, useEffect, useMemo, useState } from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
} from "@mui/material";
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import CloseIcon from "@mui/icons-material/Close";

interface I2RunDialogProps {
  open: boolean;
  command: string;
  onClose: () => void;
}

export const I2RunDialog: React.FC<I2RunDialogProps> = ({
  open,
  command,
  onClose,
}) => {
  const [cwd, setCwd] = useState<any | null>(null);
  const [config, setConfig] = useState<any | null>(null);
  const messageHandler = useCallback((event: any, data: any) => {
    if (data.message === "get-cwd") {
      setCwd(data.cwd);
    } else if (data.message === "get-config") {
      setConfig(data.config);
    }
  }, []);

  useEffect(() => {
    if (window.electronAPI) {
      window.electronAPI.sendMessage("get-cwd");
      window.electronAPI.sendMessage("get-config");
      window.electronAPI.onMessage("message-from-main", messageHandler);

      return () => {
        window.electronAPI.removeMessageListener(
          "message-from-main",
          messageHandler
        );
      };
    } else {
      console.log("window.electron is not available");
    }
  }, []);

  const fullCommand = useMemo(() => {
    if (!config || !cwd) return "";
    return `${cwd}/ccp4i2/i2run/i2run.sh ${command} --dbFile ${config.CCP4I2_PROJECTS_DIR}/db.sqlite3`;
  }, [config, cwd, command]);

  const handleCopy = () => {
    if (fullCommand) {
      navigator.clipboard.writeText(fullCommand);
    }
  };

  return (
    <Dialog open={open} onClose={onClose}>
      <DialogTitle>i2run Command</DialogTitle>
      <DialogContent>
        <code style={{ wordBreak: "break-all", whiteSpace: "pre-wrap" }}>
          {fullCommand}
        </code>
      </DialogContent>
      <DialogActions>
        <Button
          onClick={handleCopy}
          variant="outlined"
          startIcon={<ContentCopyIcon />}
        >
          Copy
        </Button>
        <Button onClick={onClose} variant="contained" startIcon={<CloseIcon />}>
          Dismiss
        </Button>
      </DialogActions>
    </Dialog>
  );
};
