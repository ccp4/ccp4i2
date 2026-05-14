import React, {
  createContext,
  useContext,
  useState,
  useCallback,
  useMemo,
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

export type SnackbarAction = {
  label: string;
  onClick: () => void;
};

type PopcornContextType = {
  setMessage: (
    msg: string,
    severity?: MessageSeverity,
    action?: SnackbarAction,
  ) => void;
  setError: (msg: string) => void; // Modal error - requires user acknowledgment
};

const PopcornContext = createContext<PopcornContextType | undefined>(undefined);

const noopPopcorn: PopcornContextType = {
  setMessage: () => {},
  setError: () => {},
};

// /ccp4i2/config sits outside the (authed) provider stack on purpose so it
// remains reachable when Django isn't running. Components that work both
// inside and outside that tree (ConfigContent etc.) get a no-op fallback
// rather than throwing.
export const usePopcorn = () => useContext(PopcornContext) ?? noopPopcorn;

export const PopcornProvider: React.FC<{ children: ReactNode }> = ({
  children,
}) => {
  const [snackbarMessage, setSnackbarMessage] = useState<string | null>(null);
  const [snackbarSeverity, setSnackbarSeverity] = useState<MessageSeverity>("info");
  const [snackbarAction, setSnackbarAction] = useState<SnackbarAction | null>(null);
  const [errorDialogMessage, setErrorDialogMessage] = useState<string | null>(null);

  // Show a snackbar message (auto-hides after duration based on severity).
  // Optional `action` renders an inline button (e.g. "Sign in" on auth errors)
  // so users can recover without typing a URL.
  const showMessage = useCallback(
    (msg: string, severity: MessageSeverity = "info", action?: SnackbarAction) => {
      setSnackbarSeverity(severity);
      setSnackbarMessage(msg);
      setSnackbarAction(action ?? null);
    },
    [],
  );

  // Show a modal error dialog (requires user to click OK)
  const showError = useCallback((msg: string) => {
    setErrorDialogMessage(msg);
  }, []);

  const handleSnackbarClose = useCallback(() => {
    setSnackbarMessage(null);
    setSnackbarAction(null);
  }, []);

  const handleErrorDialogClose = useCallback(() => {
    setErrorDialogMessage(null);
  }, []);

  // Determine auto-hide duration based on severity.
  // Snackbars with an action button don't auto-hide -- the user has a
  // decision to make (e.g. clicking Sign in), and dismissing it for them
  // defeats the purpose.
  function getAutoHideDuration(
    severity: MessageSeverity,
    hasAction: boolean,
  ): number | null {
    if (hasAction) return null;
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
  }

  const contextValue = useMemo(
    () => ({ setMessage: showMessage, setError: showError }),
    [showMessage, showError]
  );

  return (
    <PopcornContext.Provider value={contextValue}>
      {/* Snackbar for regular messages */}
      <Snackbar
        open={!!snackbarMessage}
        anchorOrigin={{ vertical: "top", horizontal: "center" }}
        onClose={handleSnackbarClose}
        autoHideDuration={getAutoHideDuration(snackbarSeverity, snackbarAction !== null)}
      >
        <Alert
          onClose={handleSnackbarClose}
          severity={snackbarSeverity}
          variant="filled"
          sx={{ width: "100%" }}
          action={
            snackbarAction ? (
              <Button
                color="inherit"
                size="small"
                onClick={() => {
                  const fn = snackbarAction.onClick;
                  handleSnackbarClose();
                  fn();
                }}
              >
                {snackbarAction.label}
              </Button>
            ) : undefined
          }
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
