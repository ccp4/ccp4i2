/**
 * The Session tab of the CCP4i2 side panel, shown when the window is a
 * recorded Moorhen session for one job: what was loaded, "save to this
 * job", what has been saved, and Finish. Only the user ends a session.
 */
import React, { useEffect, useMemo, useState } from "react";
import {
  Alert,
  Box,
  Button,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Divider,
  FormControl,
  InputLabel,
  List,
  ListItem,
  ListItemText,
  MenuItem,
  Select,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { moorhen } from "moorhen/types/moorhen";
import { usePopcorn } from "../../providers/popcorn-provider";
import type { MoorhenSessionApi } from "../../hooks/use-moorhen-session";

const KIND_LABELS: Record<string, string> = {
  dictionary: "Dictionary",
  coordinates: "Coordinates",
  map_2fofc: "Map",
  map_fofc: "Difference map",
  map_anom: "Anomalous map",
  map: "Real-space map",
};

export interface MoorhenSessionPanelProps {
  session: MoorhenSessionApi;
  molecules: moorhen.Molecule[];
}

export const MoorhenSessionPanel: React.FC<MoorhenSessionPanelProps> = ({
  session,
  molecules,
}) => {
  const { setMessage } = usePopcorn();
  const [selectedMolNo, setSelectedMolNo] = useState<number | "">("");
  const [annotation, setAnnotation] = useState("");
  const [saving, setSaving] = useState(false);
  const [confirmFinish, setConfirmFinish] = useState(false);
  const [finishing, setFinishing] = useState(false);

  const { state, job, isOpen, error } = session;
  const outputs = state?.outputs ?? [];
  const loadPlan = state?.load_plan ?? [];
  const finished = !!state?.session?.finished;

  // Moorhen types molNo as nullable; only molecules with one can be saved.
  const saveable = useMemo(
    () => molecules.filter((m) => m.molNo != null) as (moorhen.Molecule & { molNo: number })[],
    [molecules],
  );

  // Keep a sensible molecule selected as molecules come and go.
  useEffect(() => {
    if (saveable.length === 0) {
      setSelectedMolNo("");
      return;
    }
    if (selectedMolNo === "" || !saveable.some((m) => m.molNo === selectedMolNo)) {
      setSelectedMolNo(saveable[0].molNo);
    }
  }, [saveable, selectedMolNo]);

  const selectedMolecule = useMemo(
    () => saveable.find((m) => m.molNo === selectedMolNo) ?? null,
    [saveable, selectedMolNo],
  );

  const handleSave = async () => {
    if (!selectedMolecule) return;
    setSaving(true);
    try {
      const dropped = await session.saveMolecule(
        selectedMolecule,
        annotation.trim() || selectedMolecule.name,
      );
      setMessage(`Saved ${dropped.name} to job ${job?.number ?? session.jobId}`, "success");
      setAnnotation("");
    } catch (err) {
      setMessage(`Save failed: ${err instanceof Error ? err.message : String(err)}`, "error");
    } finally {
      setSaving(false);
    }
  };

  const handleFinish = async () => {
    setFinishing(true);
    try {
      const disposition = await session.finish();
      setConfirmFinish(false);
      if (disposition === "deleted") {
        setMessage("Nothing was saved; the job has been discarded. You can close this window.", "info");
      } else {
        setMessage("Session finished; the job is harvesting what you saved. You can close this window.", "success");
      }
    } catch (err) {
      setMessage(`Finish failed: ${err instanceof Error ? err.message : String(err)}`, "error");
    } finally {
      setFinishing(false);
    }
  };

  const statusChip = finished ? (
    <Chip size="small" color="default" label="Finished" />
  ) : isOpen ? (
    <Chip size="small" color="success" label="Session open" />
  ) : (
    <Chip size="small" color="warning" label={state ? "Not running" : "Connecting"} />
  );

  return (
    <Box sx={{ p: 1.5, height: "100%", overflowY: "auto" }}>
      <Stack spacing={1.5}>
        <Box>
          <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
            Moorhen session
          </Typography>
          <Typography variant="body2" color="text.secondary">
            {job ? `Job ${job.number}: ${job.title}` : `Job ${session.jobId}`}
          </Typography>
          <Box sx={{ mt: 0.5 }}>{statusChip}</Box>
        </Box>

        {error && <Alert severity="error">{error}</Alert>}
        {finished && (
          <Alert severity="info">
            This session is finished. Anything you saved is being filed into the
            job. You can close this window.
          </Alert>
        )}
        {!finished && state && !isOpen && (
          <Alert severity="warning">
            The job is no longer running, so nothing more can be saved here.
          </Alert>
        )}

        <Divider />
        <Typography variant="subtitle2">Loaded at start</Typography>
        {loadPlan.length === 0 ? (
          <Typography variant="body2" color="text.secondary">
            Nothing: this session started empty.
          </Typography>
        ) : (
          <List dense disablePadding>
            {loadPlan.map((entry, index) => (
              <ListItem key={`${entry.param}-${index}`} disableGutters sx={{ py: 0 }}>
                <ListItemText
                  primary={entry.label}
                  secondary={KIND_LABELS[entry.kind] ?? entry.kind}
                  primaryTypographyProps={{ variant: "body2", noWrap: true }}
                  secondaryTypographyProps={{ variant: "caption" }}
                />
              </ListItem>
            ))}
          </List>
        )}

        <Divider />
        <Typography variant="subtitle2">Save to this job</Typography>
        <FormControl size="small" fullWidth disabled={!isOpen || saveable.length === 0}>
          <InputLabel id="session-save-molecule">Molecule</InputLabel>
          <Select
            labelId="session-save-molecule"
            label="Molecule"
            value={selectedMolNo}
            onChange={(ev) => setSelectedMolNo(ev.target.value as number)}
          >
            {saveable.map((m) => (
              <MenuItem key={m.molNo} value={m.molNo}>
                {m.name}
              </MenuItem>
            ))}
          </Select>
        </FormControl>
        <TextField
          size="small"
          fullWidth
          label="Annotation"
          placeholder={selectedMolecule?.name ?? "What this model is"}
          value={annotation}
          onChange={(ev) => setAnnotation(ev.target.value)}
          disabled={!isOpen}
        />
        <Button
          variant="contained"
          size="small"
          onClick={handleSave}
          disabled={!isOpen || !selectedMolecule || saving}
        >
          {saving ? "Saving..." : "Save model to this job"}
        </Button>

        <Divider />
        <Typography variant="subtitle2">Saved so far ({outputs.length})</Typography>
        {outputs.length === 0 ? (
          <Typography variant="body2" color="text.secondary">
            Nothing yet. A session with nothing saved is discarded when finished.
          </Typography>
        ) : (
          <List dense disablePadding>
            {outputs.map((output) => (
              <ListItem key={output.number} disableGutters sx={{ py: 0 }}>
                <ListItemText
                  primary={output.annotation || output.name}
                  secondary={`${output.name} (${output.kind})`}
                  primaryTypographyProps={{ variant: "body2", noWrap: true }}
                  secondaryTypographyProps={{ variant: "caption" }}
                />
              </ListItem>
            ))}
          </List>
        )}

        <Divider />
        <Button
          variant="outlined"
          color={outputs.length > 0 ? "primary" : "warning"}
          size="small"
          onClick={() => setConfirmFinish(true)}
          disabled={!isOpen || finishing}
        >
          Finish session
        </Button>
        <Typography variant="caption" color="text.secondary">
          Closing this window without finishing keeps the session open if you
          have saved anything; finish it later from the job menu.
        </Typography>
      </Stack>

      <Dialog open={confirmFinish} onClose={() => setConfirmFinish(false)}>
        <DialogTitle>Finish this session?</DialogTitle>
        <DialogContent>
          {outputs.length > 0 ? (
            <Typography variant="body2">
              {outputs.length === 1
                ? "The one model you saved"
                : `The ${outputs.length} files you saved`}{" "}
              will be filed as the job&apos;s outputs. Anything not saved is lost.
            </Typography>
          ) : (
            <Typography variant="body2">
              Nothing has been saved, so the job will be discarded.
            </Typography>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setConfirmFinish(false)} disabled={finishing}>
            Keep working
          </Button>
          <Button variant="contained" onClick={handleFinish} disabled={finishing}>
            {finishing ? "Finishing..." : "Finish"}
          </Button>
        </DialogActions>
      </Dialog>
    </Box>
  );
};
