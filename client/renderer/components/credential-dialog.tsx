"use client";
/**
 * Generic credential dialog — the one place a user types a secret.
 *
 * Entirely descriptor-driven: it renders whatever fields the server declares
 * for a credential, so ssh passwords and OAuth tokens will reuse it unchanged
 * (see `docs/CREDENTIALS_DESIGN.md`).
 *
 * Three deliberate design points:
 *
 * 1. **Test before Save.** Without a live check a user cannot tell a mistyped
 *    secret from a service outage from a bad input file. The Test button is
 *    what makes an invisible value trustworthy.
 * 2. **Paste both halves at once.** Copying two long opaque strings separately
 *    from a web page is exactly where people fail, so a single paste
 *    containing all the values is split across the fields.
 * 3. **Say where it will be kept.** The server reports the real storage
 *    backend; if the platform offers no secret store we say so rather than
 *    implying protection we are not providing.
 *
 * Stored values are never pre-filled — the server never returns them. An
 * existing credential shows a placeholder instead.
 */
import {
  Alert,
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton,
  InputAdornment,
  Link,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import {
  CheckCircle,
  ErrorOutline,
  Launch,
  Visibility,
  VisibilityOff,
} from "@mui/icons-material";
import React, { useCallback, useEffect, useMemo, useState } from "react";

import { apiGet, apiPost } from "../api-fetch";

export interface CredentialFieldDescriptor {
  name: string;
  label: string;
  secret: boolean;
  help: string;
  isSet: boolean;
}

export interface CredentialDescriptor {
  name: string;
  label: string;
  description: string;
  kind: string;
  signupUrl: string;
  pasteSplit: boolean;
  canValidate: boolean;
  fields: CredentialFieldDescriptor[];
  isSet: boolean;
  source: string;
  secure: boolean;
  storageLabel: string;
  editable: boolean;
  hint: string;
  lastValidated: number | null;
  lastValidationOk: boolean | null;
  lastValidationMessage: string;
}

/** Fetch every credential descriptor. Never returns secret material. */
export async function fetchCredentials(): Promise<{
  credentials: CredentialDescriptor[];
  editable: boolean;
  secure: boolean;
  storageLabel: string;
}> {
  const resp = await apiGet<any>("config/credentials/");
  const data = resp?.data ?? resp;
  return {
    credentials: data?.credentials ?? [],
    editable: Boolean(data?.editable),
    secure: Boolean(data?.secure),
    storageLabel: data?.storageLabel ?? "",
  };
}

/**
 * Split a single pasted blob into one value per field.
 *
 * PDB-REDO's token page presents an id and a secret; users paste them together
 * far more often than one might hope. Splits on whitespace, colons, commas and
 * semicolons, and only accepts the result if it yields exactly one token per
 * field — otherwise the paste is treated as an ordinary single value.
 */
export function splitPastedValue(text: string, fieldCount: number): string[] | null {
  const parts = text
    .trim()
    .split(/[\s:,;]+/)
    .filter(Boolean);
  return parts.length === fieldCount ? parts : null;
}

interface CredentialDialogProps {
  open: boolean;
  /** Credential name, e.g. "pdb_redo". */
  name: string;
  onClose: () => void;
  /** Called after a successful save, so callers can revalidate their own state. */
  onSaved?: () => void;
}

export const CredentialDialog: React.FC<CredentialDialogProps> = ({
  open,
  name,
  onClose,
  onSaved,
}) => {
  const [descriptor, setDescriptor] = useState<CredentialDescriptor | null>(null);
  const [values, setValues] = useState<Record<string, string>>({});
  const [reveal, setReveal] = useState<Record<string, boolean>>({});
  const [busy, setBusy] = useState<"" | "testing" | "saving">("");
  const [result, setResult] = useState<{ ok: boolean; message: string } | null>(
    null
  );
  const [error, setError] = useState("");

  useEffect(() => {
    if (!open) return;
    setValues({});
    setReveal({});
    setResult(null);
    setError("");
    fetchCredentials()
      .then(({ credentials }) => {
        setDescriptor(credentials.find((c) => c.name === name) ?? null);
      })
      .catch(() => setError("Could not load credential settings."));
  }, [open, name]);

  const complete = useMemo(
    () =>
      Boolean(descriptor) &&
      descriptor!.fields.every((f) => (values[f.name] ?? "").trim().length > 0),
    [descriptor, values]
  );

  const handlePaste = useCallback(
    (event: React.ClipboardEvent<HTMLInputElement>) => {
      if (!descriptor?.pasteSplit) return;
      const text = event.clipboardData.getData("text");
      const parts = splitPastedValue(text, descriptor.fields.length);
      if (!parts) return;
      event.preventDefault();
      const next: Record<string, string> = {};
      descriptor.fields.forEach((f, i) => {
        next[f.name] = parts[i];
      });
      setValues(next);
      setResult(null);
    },
    [descriptor]
  );

  const handleTest = useCallback(async () => {
    if (!descriptor) return;
    setBusy("testing");
    setResult(null);
    setError("");
    try {
      const resp = await apiPost<any>(
        `config/credentials/${descriptor.name}/validate/`,
        { values }
      );
      const data = resp?.data ?? resp;
      setResult({ ok: Boolean(data?.ok), message: data?.message ?? "" });
    } catch {
      setError("Could not run the test.");
    } finally {
      setBusy("");
    }
  }, [descriptor, values]);

  const handleSave = useCallback(async () => {
    if (!descriptor) return;
    setBusy("saving");
    setError("");
    try {
      await apiPost<any>(`config/credentials/${descriptor.name}/set/`, {
        values,
        persistence: "keychain",
      });
      onSaved?.();
      onClose();
    } catch {
      setError("Could not save the credential.");
    } finally {
      setBusy("");
    }
  }, [descriptor, values, onSaved, onClose]);

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle>{descriptor?.label ?? "Credential"}</DialogTitle>
      <DialogContent>
        <Stack spacing={2} sx={{ mt: 1 }}>
          {descriptor?.description && (
            <Typography variant="body2" color="text.secondary">
              {descriptor.description}
            </Typography>
          )}

          {descriptor?.signupUrl && (
            <Link
              href={descriptor.signupUrl}
              target="_blank"
              rel="noopener"
              variant="body2"
              sx={{ display: "inline-flex", alignItems: "center", gap: 0.5 }}
            >
              Get a token at {new URL(descriptor.signupUrl).host}
              <Launch fontSize="inherit" />
            </Link>
          )}

          {descriptor?.isSet && (
            <Alert severity="info" variant="outlined">
              A credential is already stored
              {descriptor.hint ? ` (ending ...${descriptor.hint})` : ""}. Entering
              new values replaces it; existing values are never shown.
            </Alert>
          )}

          {descriptor?.fields.map((field, index) => (
            <TextField
              key={field.name}
              label={field.label}
              helperText={field.help}
              value={values[field.name] ?? ""}
              onChange={(e) => {
                setValues((prev) => ({ ...prev, [field.name]: e.target.value }));
                setResult(null);
              }}
              onPaste={index === 0 ? handlePaste : undefined}
              type={field.secret && !reveal[field.name] ? "password" : "text"}
              placeholder={field.isSet ? "•".repeat(12) : ""}
              autoComplete="off"
              fullWidth
              size="small"
              slotProps={
                field.secret
                  ? {
                      input: {
                        endAdornment: (
                          <InputAdornment position="end">
                            <IconButton
                              size="small"
                              edge="end"
                              aria-label={
                                reveal[field.name]
                                  ? `Hide ${field.label}`
                                  : `Show ${field.label}`
                              }
                              onClick={() =>
                                setReveal((prev) => ({
                                  ...prev,
                                  [field.name]: !prev[field.name],
                                }))
                              }
                            >
                              {reveal[field.name] ? (
                                <VisibilityOff fontSize="small" />
                              ) : (
                                <Visibility fontSize="small" />
                              )}
                            </IconButton>
                          </InputAdornment>
                        ),
                      },
                    }
                  : undefined
              }
            />
          ))}

          {descriptor?.pasteSplit && descriptor.fields.length > 1 && (
            <Typography variant="caption" color="text.secondary">
              You can paste both values at once into the first box.
            </Typography>
          )}

          {result && (
            <Alert
              severity={result.ok ? "success" : "warning"}
              icon={result.ok ? <CheckCircle /> : <ErrorOutline />}
            >
              {result.message}
            </Alert>
          )}

          {error && <Alert severity="error">{error}</Alert>}

          {descriptor?.storageLabel && (
            <Typography variant="caption" color="text.secondary">
              {descriptor.secure ? "Will be stored in " : "Warning: will be stored in "}
              {descriptor.storageLabel}.
            </Typography>
          )}
        </Stack>
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Cancel</Button>
        {descriptor?.canValidate && (
          <Button onClick={handleTest} disabled={!complete || busy !== ""}>
            {busy === "testing" ? "Testing…" : "Test"}
          </Button>
        )}
        <Button
          variant="contained"
          onClick={handleSave}
          disabled={!complete || busy !== ""}
        >
          {busy === "saving" ? "Saving…" : "Save"}
        </Button>
      </DialogActions>
    </Dialog>
  );
};

export default CredentialDialog;
