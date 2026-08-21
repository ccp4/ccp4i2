"use client";
/**
 * Inline credential banner for a task interface — the *primary* surface for
 * getting a credential set.
 *
 * The reasoning (see `docs/CREDENTIALS_DESIGN.md`): users do not browse
 * Preferences speculatively. They discover that a token is needed at the moment
 * a task will not run, so the fix belongs where the problem appears. A task
 * that needs a credential renders this at the top of its interface; when the
 * credential is present and verified it shrinks to a single quiet line.
 *
 * Purely status-driven — it asks the credential API what is configured and
 * never handles secret material itself.
 */
import { Alert, AlertTitle, Button, Stack } from "@mui/material";
import React, { useCallback, useEffect, useState } from "react";

import {
  CredentialDescriptor,
  CredentialDialog,
  fetchCredentials,
} from "../../credential-dialog";
import { apiPost } from "../../../api-fetch";

interface CredentialAlertProps {
  /** Credential name, e.g. "pdb_redo". */
  name: string;
  /** Called after the credential changes, so the task can revalidate. */
  onChanged?: () => void;
}

export const CredentialAlert: React.FC<CredentialAlertProps> = ({
  name,
  onChanged,
}) => {
  const [credential, setCredential] = useState<CredentialDescriptor | null>(null);
  const [loaded, setLoaded] = useState(false);
  const [dialogOpen, setDialogOpen] = useState(false);
  const [testing, setTesting] = useState(false);

  const load = useCallback(async () => {
    try {
      const { credentials } = await fetchCredentials();
      setCredential(credentials.find((c) => c.name === name) ?? null);
    } catch {
      setCredential(null);
    } finally {
      setLoaded(true);
    }
  }, [name]);

  useEffect(() => {
    load();
  }, [load]);

  const handleSaved = useCallback(async () => {
    await load();
    onChanged?.();
  }, [load, onChanged]);

  const handleTest = useCallback(async () => {
    setTesting(true);
    try {
      await apiPost<any>(`config/credentials/${name}/validate/`, {});
      await load();
      onChanged?.();
    } catch {
      /* the reloaded descriptor carries the outcome */
    } finally {
      setTesting(false);
    }
  }, [name, load, onChanged]);

  // Say nothing until we know — a banner that flashes "no token" on every task
  // open would be worse than no banner.
  if (!loaded || !credential) return null;

  const dialog = dialogOpen && (
    <CredentialDialog
      open
      name={name}
      onClose={() => setDialogOpen(false)}
      onSaved={handleSaved}
    />
  );

  if (!credential.isSet) {
    return (
      <>
        <Alert
          severity="warning"
          sx={{ mb: 1 }}
          action={
            credential.editable ? (
              <Button
                size="small"
                variant="contained"
                onClick={() => setDialogOpen(true)}
              >
                Set token…
              </Button>
            ) : undefined
          }
        >
          <AlertTitle>{credential.label} needs a token</AlertTitle>
          {credential.description}
        </Alert>
        {dialog}
      </>
    );
  }

  if (credential.lastValidationOk === false) {
    return (
      <>
        <Alert
          severity="error"
          sx={{ mb: 1 }}
          action={
            credential.editable ? (
              <Button size="small" onClick={() => setDialogOpen(true)}>
                Replace
              </Button>
            ) : undefined
          }
        >
          <AlertTitle>{credential.label} rejected your token</AlertTitle>
          {credential.lastValidationMessage}
        </Alert>
        {dialog}
      </>
    );
  }

  // Configured: one quiet line, with the means to check it still works.
  return (
    <>
      <Alert
        severity="success"
        variant="outlined"
        sx={{ mb: 1, py: 0 }}
        action={
          credential.canValidate ? (
            <Stack direction="row" spacing={1}>
              <Button size="small" onClick={handleTest} disabled={testing}>
                {testing ? "Testing…" : "Test"}
              </Button>
              {credential.editable && (
                <Button size="small" onClick={() => setDialogOpen(true)}>
                  Replace
                </Button>
              )}
            </Stack>
          ) : undefined
        }
      >
        {credential.label} token configured
        {credential.hint ? ` (ends ...${credential.hint})` : ""}.
      </Alert>
      {dialog}
    </>
  );
};

export default CredentialAlert;
