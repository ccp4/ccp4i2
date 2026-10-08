"use client";
import {
  createContext,
  Dispatch,
  PropsWithChildren,
  useContext,
  useReducer,
} from "react";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
} from "@mui/material";

interface DeleteDialogState {
  open: boolean;
  what?: string;
  children?: React.ReactNode[];
  deleteDisabled?: boolean;
  onDelete?: () => void;
  onCancel?: () => void;
}

interface DeleteDialogAction {
  /**
   * "update" changes whether Delete is enabled while the dialog is open: the
   * dialog's own content may hold a choice that decides it.
   */
  type: "show" | "hide" | "update";
  what?: string;
  children?: React.ReactNode[];
  deleteDisabled?: boolean;
  onDelete?: () => void;
  onCancel?: () => void;
}

const DeleteDialogContext = createContext<Dispatch<DeleteDialogAction> | null>(
  null
);

interface DeleteDialogProviderProps extends PropsWithChildren {
  deleteDisabled?: boolean;
}

export function DeleteDialogProvider(props: DeleteDialogProviderProps) {
  const [state, dispatch] = useReducer(deleteDialogReducer, { open: false });

  function handleCancel() {
    if (state.onCancel) state.onCancel();
    dispatch({ type: "hide" });
  }

  function handleDelete() {
    if (state.onDelete) state.onDelete?.();
    dispatch({ type: "hide" });
  }

  // The dialog sits inside the provider so that its content can reach it.
  return (
    <DeleteDialogContext.Provider value={dispatch}>
      <Dialog open={state.open}>
        <DialogTitle>{`Delete ${state.what}?`}</DialogTitle>
        <DialogContent>
          <DialogContentText key="DialogContentText">
            This action cannot be undone.
          </DialogContentText>
          {state.children}
        </DialogContent>
        <DialogActions>
          <Button key="Close" autoFocus onClick={handleCancel}>
            Cancel
          </Button>
          <Button
            key="Delete"
            disabled={state.deleteDisabled}
            variant="contained"
            onClick={handleDelete}
          >
            Delete
          </Button>
        </DialogActions>
      </Dialog>
      {props.children}
    </DeleteDialogContext.Provider>
  );
}

export function useDeleteDialog() {
  return useContext(DeleteDialogContext);
}

export function deleteDialogReducer(
  state: DeleteDialogState,
  action: DeleteDialogAction
): DeleteDialogState {
  switch (action.type) {
    case "show":
      return {
        open: true,
        what: action.what,
        children: action.children,
        deleteDisabled: action.deleteDisabled,
        onDelete: action.onDelete,
        onCancel: action.onCancel,
      };
    case "update":
      return state.open
        ? { ...state, deleteDisabled: action.deleteDisabled }
        : state;
    case "hide":
      return { open: false };
  }
}
