/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
  type: "show" | "hide";
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

  return (
    <>
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
      <DeleteDialogContext.Provider value={dispatch}>
        {props.children}
      </DeleteDialogContext.Provider>
    </>
  );
}

export function useDeleteDialog() {
  return useContext(DeleteDialogContext);
}

function deleteDialogReducer(
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
    case "hide":
      return { open: false };
  }
}
