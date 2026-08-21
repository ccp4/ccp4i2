"use client";
/**
 * Preferences → Credentials: the *management* surface for API tokens and
 * passwords.
 *
 * Deliberately the secondary surface, not the primary one. Nobody visits
 * Preferences speculatively; users meet the need when a task refuses to run, so
 * the task itself carries the "Set token…" affordance (see
 * `credential-alert.tsx` and the run dialog's quick actions). This panel is
 * where you come back to re-test, replace or remove a credential.
 *
 * Reads/writes the credential API. No secret material is ever fetched — the
 * panel shows only whether something is stored, where it came from, and how the
 * last live check went.
 */
import {
  Alert,
  Box,
  Button,
  Chip,
  CircularProgress,
  Stack,
  Tooltip,
  Typography,
} from "@mui/material";
import { CheckCircle, ErrorOutline, HelpOutline, Lock, LockOpen } from "@mui/icons-material";
import React, { useCallback, useEffect, useState } from "react";

import { apiPost } from "../api-fetch";
import {
  CredentialDescriptor,
  CredentialDialog,
  fetchCredentials,
} from "./credential-dialog";

/** Where a value came from, in words the user can act on. */
const SOURCE_LABEL: Record<string, string> = {
  environment: "from environment",
  session: "this session only",
  keychain: "system keychain",
  file: "local file",
  preferences: "legacy preferences file",
  unset: "not set",
};

function statusChip(credential: CredentialDescriptor) {
  if (!credential.isSet) {
    return <Chip size="small" label="Not set" />;
  }
  if (credential.lastValidationOk === true) {
    return (
      <Chip
        size="small"
        color="success"
        icon={<CheckCircle />}
        label="Verified"
      />
    );
  }
  if (credential.lastValidationOk === false) {
    return (
      <Chip
        size="small"
        color="error"
        icon={<ErrorOutline />}
        label="Rejected"
      />
    );
  }
  return <Chip size="small" color="info" icon={<HelpOutline />} label="Set, untested" />;
}

export function CredentialsPanel() {
  const [credentials, setCredentials] = useState<CredentialDescriptor[]>([]);
  const [editable, setEditable] = useState(false);
  const [secure, setSecure] = useState(false);
  const [storageLabel, setStorageLabel] = useState("");
  const [loading, setLoading] = useState(true);
  const [busy, setBusy] = useState("");
  const [dialogFor, setDialogFor] = useState<string | null>(null);
  const [notice, setNotice] = useState<{ ok: boolean; message: string } | null>(
    null
  );

  const load = useCallback(async () => {
    const data = await fetchCredentials();
    setCredentials(data.credentials);
    setEditable(data.editable);
    setSecure(data.secure);
    setStorageLabel(data.storageLabel);
  }, []);

  useEffect(() => {
    (async () => {
      setLoading(true);
      try {
        await load();
      } finally {
        setLoading(false);
      }
    })();
  }, [load]);

  const handleTest = useCallback(
    async (name: string) => {
      setBusy(name);
      setNotice(null);
      try {
        const resp = await apiPost<any>(`config/credentials/${name}/validate/`, {});
        const data = resp?.data ?? resp;
        setNotice({ ok: Boolean(data?.ok), message: data?.message ?? "" });
        await load();
      } catch {
        setNotice({ ok: false, message: "Could not run the test." });
      } finally {
        setBusy("");
      }
    },
    [load]
  );

  const handleClear = useCallback(
    async (name: string) => {
      setBusy(name);
      setNotice(null);
      try {
        await apiPost<any>(`config/credentials/${name}/clear/`, {});
        await load();
      } catch {
        setNotice({ ok: false, message: "Could not remove the credential." });
      } finally {
        setBusy("");
      }
    },
    [load]
  );

  if (loading) {
    return (
      <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
        <CircularProgress />
      </Box>
    );
  }

  return (
    <Stack spacing={3} sx={{ maxWidth: 720, mx: "auto", p: 2 }}>
      <Box>
        <Typography variant="h6">Credentials</Typography>
        <Typography variant="body2" color="text.secondary">
          Tokens and passwords for services CCP4i2 talks to on your behalf.
          Secrets are stored outside your project data and are never included in
          an exported project.
        </Typography>
        {!editable && (
          <Typography variant="body2" color="warning.main" sx={{ mt: 1 }}>
            Read-only: in a server deployment credentials are supplied as
            environment variables.
          </Typography>
        )}
        {editable && storageLabel && (
          <Stack
            direction="row"
            spacing={0.5}
            alignItems="center"
            sx={{ mt: 1, color: secure ? "text.secondary" : "warning.main" }}
          >
            {secure ? <Lock fontSize="small" /> : <LockOpen fontSize="small" />}
            <Typography variant="body2">
              Stored in {storageLabel}.
            </Typography>
          </Stack>
        )}
      </Box>

      {notice && (
        <Alert severity={notice.ok ? "success" : "warning"}>{notice.message}</Alert>
      )}

      <Stack spacing={2}>
        {credentials.map((credential) => (
          <Stack
            key={credential.name}
            direction="row"
            spacing={1}
            alignItems="center"
            sx={{
              border: 1,
              borderColor: "divider",
              borderRadius: 1,
              p: 1.5,
            }}
          >
            <Box sx={{ flexGrow: 1, minWidth: 0 }}>
              <Stack direction="row" spacing={1} alignItems="center">
                <Typography variant="subtitle2">{credential.label}</Typography>
                {statusChip(credential)}
              </Stack>
              <Typography variant="caption" color="text.secondary">
                {credential.isSet
                  ? `${SOURCE_LABEL[credential.source] ?? credential.source}${
                      credential.hint ? ` · ends ...${credential.hint}` : ""
                    }`
                  : "No credential configured."}
                {credential.lastValidationMessage
                  ? ` · ${credential.lastValidationMessage}`
                  : ""}
              </Typography>
            </Box>

            {credential.canValidate && credential.isSet && (
              <Button
                size="small"
                onClick={() => handleTest(credential.name)}
                disabled={busy !== ""}
              >
                Test
              </Button>
            )}
            <Tooltip
              title={
                credential.editable
                  ? ""
                  : "Supplied by this deployment's configuration"
              }
            >
              <span>
                <Button
                  size="small"
                  variant="outlined"
                  onClick={() => setDialogFor(credential.name)}
                  disabled={!credential.editable || busy !== ""}
                >
                  {credential.isSet ? "Replace" : "Set…"}
                </Button>
              </span>
            </Tooltip>
            {credential.isSet && credential.editable && (
              <Button
                size="small"
                color="error"
                onClick={() => handleClear(credential.name)}
                disabled={busy !== ""}
              >
                Clear
              </Button>
            )}
          </Stack>
        ))}
        {credentials.length === 0 && (
          <Typography variant="body2" color="text.secondary">
            No services requiring credentials are registered.
          </Typography>
        )}
      </Stack>

      {dialogFor && (
        <CredentialDialog
          open
          name={dialogFor}
          onClose={() => setDialogFor(null)}
          onSaved={load}
        />
      )}
    </Stack>
  );
}

export default CredentialsPanel;
