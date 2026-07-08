"use client";
import React, { useState } from "react";
import {
  Alert,
  AlertTitle,
  Box,
  Button,
  Chip,
  CircularProgress,
  FormControl,
  FormControlLabel,
  FormLabel,
  Checkbox,
  Radio,
  RadioGroup,
  Stack,
  Step,
  StepLabel,
  Stepper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TextField,
  Typography,
} from "@mui/material";
import {
  ArrowForward,
  CheckCircle,
  ContentCopy,
  Error as ErrorIcon,
  FolderSpecial,
  PlayArrow,
  Warning,
} from "@mui/icons-material";
import { useApi } from "../api";
import { apiJson } from "../api-fetch";
import { usePopcorn } from "../providers/popcorn-provider";
import { useAuth } from "@/lib/compounds/auth-context";
import type {
  ImportResult,
  ImportStatus,
  StructureIssue,
  ValidateReport,
} from "../types/migration";

const STEPS = ["Source", "Validate", "Dry run", "Import"];

const DEFAULT_DB_PATH = "~/.CCP4I2/db/database.sqlite";

type CopyMode = "copy" | "schlurp";

/** Friendly labels for the machine-readable issue types. */
const ISSUE_TITLES: Record<string, string> = {
  nested_project: "Nested project",
  dest_collision: "Destination name collision",
  case_clash: "Name clash (case-insensitive filesystem)",
  unicode_clash: "Name clash (unicode normalisation)",
  shared_directory: "Shared project directory",
  shared_subtree: "Shared content directory",
  not_project_root: "Not a project root",
  dest_unwritable: "Destination not writable",
  insufficient_space: "Insufficient disk space",
  unreadable_source: "Unreadable source",
  i1_dir_overlap: "i1 directory overlap",
};

function issueColor(severity: string): "error" | "warning" | "info" {
  if (severity === "blocking") return "error";
  if (severity === "warning") return "warning";
  return "info";
}

export const MigrateLegacyContent: React.FC = () => {
  const api = useApi();
  const { setMessage, setError } = usePopcorn();
  const { canAdminister } = useAuth();

  const [activeStep, setActiveStep] = useState(0);

  // Step 1 — source & mode
  const [dbPath, setDbPath] = useState(DEFAULT_DB_PATH);
  const [copyMode, setCopyMode] = useState<CopyMode>("copy");
  const [remapFrom, setRemapFrom] = useState("");
  const [remapTo, setRemapTo] = useState("");
  const [acknowledged, setAcknowledged] = useState(false);

  // Results
  const [report, setReport] = useState<ValidateReport | null>(null);
  const [dryRun, setDryRun] = useState<ImportResult | null>(null);
  const [imported, setImported] = useState<ImportResult | null>(null);
  const [status, setStatus] = useState<ImportStatus | null>(null);
  const [busy, setBusy] = useState(false);

  if (!canAdminister) {
    return (
      <Alert severity="error" sx={{ maxWidth: 640, mx: "auto", mt: 4 }}>
        <AlertTitle>Administrator access required</AlertTitle>
        Migrating a legacy database is an administrative operation. Please sign in
        as an administrator to continue.
      </Alert>
    );
  }

  const copyFiles = copyMode === "copy";

  /** Build the shared multipart body for all three endpoints. */
  const buildForm = (extra: Record<string, string> = {}) => {
    const fd = new FormData();
    fd.append("db_path", dbPath.trim());
    fd.append("copy_files", String(copyFiles));
    if (remapFrom.trim() && remapTo.trim()) {
      fd.append("remap_from", remapFrom.trim());
      fd.append("remap_to", remapTo.trim());
    }
    Object.entries(extra).forEach(([k, v]) => fd.append(k, v));
    return fd;
  };

  const blocking = report?.summary.blocking_issues ?? 0;
  const canAdvanceFromValidate = !!report && (blocking === 0 || acknowledged);

  // ---- Step actions -------------------------------------------------------

  const runValidate = async () => {
    setBusy(true);
    setReport(null);
    setDryRun(null);
    try {
      const r = await api.post<ValidateReport>(
        "admin/validate-sqlite/",
        buildForm()
      );
      setReport(r);
      setAcknowledged(false);
      setActiveStep(1);
    } catch (e: any) {
      setError(`Validation failed: ${e?.message ?? e}`);
    } finally {
      setBusy(false);
    }
  };

  const runDryRun = async () => {
    setBusy(true);
    setDryRun(null);
    try {
      const r = await api.post<ImportResult>(
        "admin/import-sqlite/",
        buildForm({ dry_run: "true", allow_structural_issues: "true" })
      );
      setDryRun(r);
      setActiveStep(2);
    } catch (e: any) {
      setError(`Dry run failed: ${friendlyError(e)}`);
    } finally {
      setBusy(false);
    }
  };

  const runImport = async () => {
    setBusy(true);
    try {
      const r = await api.post<ImportResult>(
        "admin/import-sqlite/",
        buildForm({ dry_run: "false", allow_structural_issues: "true" })
      );
      setImported(r);
      setActiveStep(3);
      setMessage("Migration complete", "success");
      // Confirmation counts (one-shot GET; api.get is an SWR hook, so use apiJson).
      const s = await apiJson<ImportStatus>(
        "/api/proxy/ccp4i2/admin/import-status/"
      ).catch(() => null);
      if (s) setStatus(s);
    } catch (e: any) {
      setError(`Import failed: ${friendlyError(e)}`);
    } finally {
      setBusy(false);
    }
  };

  // ---- Render -------------------------------------------------------------

  return (
    <Box sx={{ maxWidth: 900, mx: "auto" }}>
      <Typography variant="h4" gutterBottom>
        Migrate a legacy CCP4i2 database
      </Typography>
      <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
        Bring projects, jobs and files from an old (Qt-era) CCP4i2 installation
        into this one. Migration moves the database records; project
        directories are either copied in or adopted where they are.
      </Typography>

      <Stepper activeStep={activeStep} sx={{ mb: 4 }}>
        {STEPS.map((label) => (
          <Step key={label}>
            <StepLabel>{label}</StepLabel>
          </Step>
        ))}
      </Stepper>

      {activeStep === 0 && (
        <SourceStep
          dbPath={dbPath}
          setDbPath={setDbPath}
          copyMode={copyMode}
          setCopyMode={setCopyMode}
          remapFrom={remapFrom}
          setRemapFrom={setRemapFrom}
          remapTo={remapTo}
          setRemapTo={setRemapTo}
          busy={busy}
          onValidate={runValidate}
        />
      )}

      {activeStep === 1 && report && (
        <ValidateStep
          report={report}
          copyFiles={copyFiles}
          acknowledged={acknowledged}
          setAcknowledged={setAcknowledged}
          canAdvance={canAdvanceFromValidate}
          busy={busy}
          onBack={() => setActiveStep(0)}
          onNext={runDryRun}
        />
      )}

      {activeStep === 2 && dryRun && (
        <ResultStep
          title="Dry run — nothing was written"
          result={dryRun}
          busy={busy}
          onBack={() => setActiveStep(1)}
          primaryLabel="Run import"
          primaryIcon={<PlayArrow />}
          onPrimary={runImport}
        />
      )}

      {activeStep === 3 && imported && (
        <DoneStep result={imported} status={status} />
      )}
    </Box>
  );
};

// ---------------------------------------------------------------------------
// Step 1 — source & copy mode
// ---------------------------------------------------------------------------

const SourceStep: React.FC<{
  dbPath: string;
  setDbPath: (v: string) => void;
  copyMode: CopyMode;
  setCopyMode: (v: CopyMode) => void;
  remapFrom: string;
  setRemapFrom: (v: string) => void;
  remapTo: string;
  setRemapTo: (v: string) => void;
  busy: boolean;
  onValidate: () => void;
}> = (p) => (
  <Stack spacing={3}>
    <TextField
      label="Legacy database path"
      helperText="Path on this machine to the old CCP4i2 SQLite database"
      value={p.dbPath}
      onChange={(e) => p.setDbPath(e.target.value)}
      fullWidth
    />

    <FormControl>
      <FormLabel>How should project files be handled?</FormLabel>
      <RadioGroup
        value={p.copyMode}
        onChange={(e) => p.setCopyMode(e.target.value as CopyMode)}
      >
        <FormControlLabel
          value="copy"
          control={<Radio />}
          label={
            <Box>
              <Typography>Copy projects into this installation (recommended)</Typography>
              <Typography variant="caption" color="text.secondary">
                Your old installation is left untouched; you can delete it
                afterwards. Uses extra disk space.
              </Typography>
            </Box>
          }
        />
        <FormControlLabel
          value="schlurp"
          control={<Radio />}
          label={
            <Box>
              <Typography>Adopt projects where they are</Typography>
              <Typography variant="caption" color="text.secondary">
                Faster, no extra disk, but this installation reads from your old
                directories — don&apos;t delete them. (Nested projects are still
                copied out.)
              </Typography>
            </Box>
          }
        />
      </RadioGroup>
    </FormControl>

    <Box>
      <Typography variant="subtitle2" gutterBottom>
        Path remapping (optional)
      </Typography>
      <Typography variant="caption" color="text.secondary">
        If the database was created on another machine or under a different home
        directory, remap the old path prefix to where the files live now.
      </Typography>
      <Stack direction="row" spacing={2} sx={{ mt: 1 }}>
        <TextField
          label="From prefix"
          value={p.remapFrom}
          onChange={(e) => p.setRemapFrom(e.target.value)}
          fullWidth
          size="small"
        />
        <TextField
          label="To prefix"
          value={p.remapTo}
          onChange={(e) => p.setRemapTo(e.target.value)}
          fullWidth
          size="small"
        />
      </Stack>
    </Box>

    <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
      <Button
        variant="contained"
        endIcon={p.busy ? <CircularProgress size={16} /> : <ArrowForward />}
        disabled={p.busy || !p.dbPath.trim()}
        onClick={p.onValidate}
      >
        Validate
      </Button>
    </Box>
  </Stack>
);

// ---------------------------------------------------------------------------
// Step 2 — validation report + structure/plan
// ---------------------------------------------------------------------------

const ValidateStep: React.FC<{
  report: ValidateReport;
  copyFiles: boolean;
  acknowledged: boolean;
  setAcknowledged: (v: boolean) => void;
  canAdvance: boolean;
  busy: boolean;
  onBack: () => void;
  onNext: () => void;
}> = ({ report, acknowledged, setAcknowledged, canAdvance, busy, onBack, onNext }) => {
  const s = report.summary;
  const blocking = s.blocking_issues;
  return (
    <Stack spacing={3}>
      <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
        <RatioChip label="Projects on disk" value={s.projects_on_disk} />
        <RatioChip label="Jobs on disk" value={s.jobs_on_disk} />
        <RatioChip label="Files on disk" value={s.files_on_disk} />
        <RatioChip label="Import sources" value={s.import_sources_on_disk} />
      </Stack>

      {(s.integrity_issues > 0 || s.data_quality_issues > 0) && (
        <Alert severity="warning">
          <AlertTitle>Data checks</AlertTitle>
          {report.integrity.issues.concat(report.data_quality.issues).map((m, i) => (
            <Typography key={i} variant="body2">
              • {m}
            </Typography>
          ))}
        </Alert>
      )}

      {(report.projects.dir_missing.length > 0 ||
        report.jobs.dir_missing_count > 0) && (
        <Alert severity="warning">
          <AlertTitle>Some directories are missing on disk</AlertTitle>
          {report.projects.dir_missing.slice(0, 10).map((d, i) => (
            <Typography key={i} variant="body2">
              • {d.project}: {d.directory}
            </Typography>
          ))}
          {report.jobs.dir_missing_count > 0 && (
            <Typography variant="body2">
              • {report.jobs.dir_missing_count} job director
              {report.jobs.dir_missing_count === 1 ? "y" : "ies"} not found
            </Typography>
          )}
        </Alert>
      )}

      <IssueList issues={report.structure.issues} />

      <PlanTable report={report} />

      {blocking > 0 && (
        <Alert severity="error">
          <AlertTitle>
            {blocking} issue{blocking === 1 ? "" : "s"} need
            {blocking === 1 ? "s" : ""} your attention
          </AlertTitle>
          These must be acknowledged before continuing. Review the issues above.
          <FormControlLabel
            sx={{ display: "block", mt: 1 }}
            control={
              <Checkbox
                checked={acknowledged}
                onChange={(e) => setAcknowledged(e.target.checked)}
              />
            }
            label="I understand and want to proceed anyway"
          />
        </Alert>
      )}

      <Box sx={{ display: "flex", justifyContent: "space-between" }}>
        <Button onClick={onBack} disabled={busy}>
          Back
        </Button>
        <Button
          variant="contained"
          endIcon={busy ? <CircularProgress size={16} /> : <ArrowForward />}
          disabled={busy || !canAdvance}
          onClick={onNext}
        >
          Dry run
        </Button>
      </Box>
    </Stack>
  );
};

const RatioChip: React.FC<{ label: string; value: string }> = ({ label, value }) => {
  const [have, total] = value.split("/").map((x) => parseInt(x, 10));
  const complete = Number.isFinite(have) && have === total;
  return (
    <Chip
      icon={complete ? <CheckCircle /> : <Warning />}
      color={complete ? "success" : "warning"}
      variant="outlined"
      label={`${label}: ${value}`}
    />
  );
};

const IssueList: React.FC<{ issues: StructureIssue[] }> = ({ issues }) => {
  if (!issues.length) return null;
  return (
    <Stack spacing={1}>
      <Typography variant="subtitle1">Structural checks</Typography>
      {issues.map((issue, i) => (
        <Alert key={i} severity={issueColor(issue.severity)} icon={<Warning />}>
          <AlertTitle sx={{ mb: 0 }}>
            {ISSUE_TITLES[issue.type] ?? issue.type}
          </AlertTitle>
          <Typography variant="body2">{issue.detail}</Typography>
          {issue.resolution && (
            <Typography variant="body2" sx={{ mt: 0.5, fontStyle: "italic" }}>
              → {issue.resolution}
            </Typography>
          )}
        </Alert>
      ))}
    </Stack>
  );
};

const PlanTable: React.FC<{ report: ValidateReport }> = ({ report }) => {
  const ps = report.summary.plan_summary;
  return (
    <Box>
      <Typography variant="subtitle1" gutterBottom>
        Migration plan — {ps.in_place} in place, {ps.copy} copied
        {ps.copied_due_to_nesting > 0
          ? ` (${ps.copied_due_to_nesting} un-nested)`
          : ""}
      </Typography>
      <TableContainer sx={{ maxHeight: 320, overflow: "auto" }}>
        <Table size="small" stickyHeader>
          <TableHead>
            <TableRow>
              <TableCell>Project</TableCell>
              <TableCell>Action</TableCell>
              <TableCell>Destination</TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {report.plan.map((e) => (
              <TableRow key={e.legacy_project_id}>
                <TableCell>
                  {e.name}
                  {e.renamed_to && (
                    <Typography variant="caption" color="text.secondary">
                      {" "}
                      → {e.renamed_to}
                    </Typography>
                  )}
                  {!e.source_exists && (
                    <Chip
                      size="small"
                      color="warning"
                      label="source missing"
                      sx={{ ml: 1 }}
                    />
                  )}
                </TableCell>
                <TableCell>
                  {e.mode === "copy" ? (
                    <Chip
                      size="small"
                      icon={e.reason === "nested" ? <FolderSpecial /> : <ContentCopy />}
                      label={e.reason === "nested" ? "copy (un-nest)" : "copy"}
                      color="primary"
                      variant="outlined"
                    />
                  ) : (
                    <Chip size="small" label="in place" variant="outlined" />
                  )}
                </TableCell>
                <TableCell>
                  <Typography variant="caption" sx={{ wordBreak: "break-all" }}>
                    {e.dest}
                  </Typography>
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>
    </Box>
  );
};

// ---------------------------------------------------------------------------
// Step 3 — dry-run result
// ---------------------------------------------------------------------------

const ResultStep: React.FC<{
  title: string;
  result: ImportResult;
  busy: boolean;
  onBack: () => void;
  primaryLabel: string;
  primaryIcon: React.ReactNode;
  onPrimary: () => void;
}> = ({ title, result, busy, onBack, primaryLabel, primaryIcon, onPrimary }) => (
  <Stack spacing={3}>
    <Alert severity="info">
      <AlertTitle>{title}</AlertTitle>
      This is what a real import would do. Review the counts, then proceed.
    </Alert>
    <StatsTable stats={result.stats} />
    {result.errors.length > 0 && (
      <Alert severity="error">
        <AlertTitle>{result.errors.length} error(s)</AlertTitle>
        {result.errors.slice(0, 20).map((e, i) => (
          <Typography key={i} variant="body2">
            • {e}
          </Typography>
        ))}
      </Alert>
    )}
    <Box sx={{ display: "flex", justifyContent: "space-between" }}>
      <Button onClick={onBack} disabled={busy}>
        Back
      </Button>
      <Button
        variant="contained"
        color="primary"
        endIcon={busy ? <CircularProgress size={16} /> : primaryIcon}
        disabled={busy}
        onClick={onPrimary}
      >
        {primaryLabel}
      </Button>
    </Box>
  </Stack>
);

// ---------------------------------------------------------------------------
// Step 4 — done
// ---------------------------------------------------------------------------

const DoneStep: React.FC<{ result: ImportResult; status: ImportStatus | null }> = ({
  result,
  status,
}) => (
  <Stack spacing={3}>
    <Alert severity={result.errors.length ? "warning" : "success"} icon={<CheckCircle />}>
      <AlertTitle>
        {result.errors.length ? "Migration completed with warnings" : "Migration complete"}
      </AlertTitle>
      {result.errors.length
        ? `${result.errors.length} record(s) could not be imported — see below.`
        : "Your legacy projects are now available in this installation."}
    </Alert>
    <StatsTable stats={result.stats} />
    {status && (
      <Box>
        <Typography variant="subtitle1" gutterBottom>
          Database now contains
        </Typography>
        <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
          <Chip label={`${status.projects.total} projects`} />
          <Chip label={`${status.jobs.total} jobs`} />
          <Chip label={`${status.files.total} files`} />
        </Stack>
      </Box>
    )}
    {result.errors.length > 0 && (
      <Alert severity="error" icon={<ErrorIcon />}>
        {result.errors.slice(0, 20).map((e, i) => (
          <Typography key={i} variant="body2">
            • {e}
          </Typography>
        ))}
      </Alert>
    )}
  </Stack>
);

// ---------------------------------------------------------------------------
// Shared bits
// ---------------------------------------------------------------------------

const StatsTable: React.FC<{ stats: Record<string, number> }> = ({ stats }) => {
  const rows = Object.entries(stats).filter(([, v]) => v);
  if (!rows.length)
    return (
      <Typography variant="body2" color="text.secondary">
        No records to import.
      </Typography>
    );
  return (
    <TableContainer sx={{ maxHeight: 320, overflow: "auto" }}>
      <Table size="small" stickyHeader>
        <TableHead>
          <TableRow>
            <TableCell>Entity</TableCell>
            <TableCell align="right">Count</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {rows.map(([k, v]) => (
            <TableRow key={k}>
              <TableCell>{k.replace(/_/g, " ")}</TableCell>
              <TableCell align="right">{v}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  );
};

/** Map a thrown API error (message is the machine `error` string) to prose. */
function friendlyError(e: any): string {
  const msg = String(e?.message ?? e);
  if (msg.includes("structural_issues_unacknowledged"))
    return "blocking structural issues were not acknowledged — go back and review them";
  if (msg.includes("import_failed"))
    return "the import ran but some records failed — see the counts";
  if (msg.includes("db_not_found")) return "the database file was not found";
  return msg;
}
