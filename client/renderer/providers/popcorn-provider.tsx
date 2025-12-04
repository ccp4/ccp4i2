import React, {
  createContext,
  useContext,
  useState,
  useCallback,
  ReactNode,
} from "react";
import Snackbar from "@mui/material/Snackbar";
import Alert from "@mui/material/Alert";
import Dialog from "@mui/material/Dialog";
import DialogTitle from "@mui/material/DialogTitle";
import DialogContent from "@mui/material/DialogContent";
import DialogContentText from "@mui/material/DialogContentText";
import DialogActions from "@mui/material/DialogActions";
import Button from "@mui/material/Button";

type MessageSeverity = "info" | "success" | "warning" | "error";

type PopcornContextType = {
  setMessage: (msg: string, severity?: MessageSeverity) => void;
  setError: (msg: string) => void; // Modal error - requires user acknowledgment
};

const PopcornContext = createContext<PopcornContextType | undefined>(undefined);

export const usePopcorn = () => {
  const ctx = useContext(PopcornContext);
  if (!ctx) throw new Error("usePopcorn must be used within a PopcornProvider");
  return ctx;
};

export const PopcornProvider: React.FC<{ children: ReactNode }> = ({
  children,
}) => {
  const [snackbarMessage, setSnackbarMessage] = useState<string | null>(null);
  const [snackbarSeverity, setSnackbarSeverity] = useState<MessageSeverity>("info");
  const [errorDialogMessage, setErrorDialogMessage] = useState<string | null>(null);

  // Show a snackbar message (auto-hides after duration based on severity)
  const showMessage = useCallback((msg: string, severity: MessageSeverity = "info") => {
    setSnackbarSeverity(severity);
    setSnackbarMessage(msg);
  }, []);

  // Show a modal error dialog (requires user to click OK)
  const showError = useCallback((msg: string) => {
    setErrorDialogMessage(msg);
  }, []);

  const handleSnackbarClose = useCallback(() => {
    setSnackbarMessage(null);
  }, []);

  const handleErrorDialogClose = useCallback(() => {
    setErrorDialogMessage(null);
  }, []);

  // Determine auto-hide duration based on severity
  const getAutoHideDuration = (severity: MessageSeverity) => {
    switch (severity) {
      case "error":
        return 6000; // Errors stay longer
      case "warning":
        return 5000;
      case "success":
        return 3000;
      default:
        return 4000;
    }
  };

  return (
    <PopcornContext.Provider value={{ setMessage: showMessage, setError: showError }}>
      {/* Snackbar for regular messages */}
      <Snackbar
        open={!!snackbarMessage}
        anchorOrigin={{ vertical: "top", horizontal: "center" }}
        onClose={handleSnackbarClose}
        autoHideDuration={getAutoHideDuration(snackbarSeverity)}
      >
        <Alert
          onClose={handleSnackbarClose}
          severity={snackbarSeverity}
          variant="filled"
          sx={{ width: "100%" }}
        >
          {snackbarMessage}
        </Alert>
      </Snackbar>

      {/* Modal dialog for critical errors */}
      <Dialog
        open={!!errorDialogMessage}
        onClose={handleErrorDialogClose}
        aria-labelledby="error-dialog-title"
        aria-describedby="error-dialog-description"
      >
        <DialogTitle id="error-dialog-title" sx={{ color: "error.main" }}>
          Error
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="error-dialog-description">
            {errorDialogMessage}
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleErrorDialogClose} color="primary" autoFocus>
            OK
          </Button>
        </DialogActions>
      </Dialog>

      {children}
    </PopcornContext.Provider>
  );
};

// Usage in a child component:
// const { setMessage, setError } = usePopcorn();
// setMessage("Operation completed"); // Shows snackbar
// setMessage("Warning", "warning"); // Shows warning snackbar
// setError("Critical error occurred"); // Shows modal dialog
