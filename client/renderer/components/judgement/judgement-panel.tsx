"use client";
/**
 * The Judgement tab: the task's judgement applied to this job, drawn.
 *
 * The verdict comes from GET jobs/<id>/judgement/, which judges a finished job
 * once and keeps it (judgement.json in the job directory); the judgement file
 * itself (traps, task title) from GET agent/tasks/<task>/. Each next step is a
 * follow-on job or a rerun of this one: "Create job" / "Rerun" makes it with
 * its inputs set (POST jobs/<id>/apply_next/) and opens it pending, for the
 * user to check and run. See server/ccp4i2/agent/ and docs/agentic-knowledge.md.
 */
import { useCallback, useState } from "react";
import {
  Alert,
  Box,
  Button,
  Card,
  CardActions,
  CardContent,
  Chip,
  CircularProgress,
  LinearProgress,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  CheckCircle,
  Cancel,
  HelpOutline,
  PlayArrow,
  Replay,
} from "@mui/icons-material";
import { useRouter } from "next/navigation";
import { useApi } from "../../api";
import { usePopcorn } from "../../providers/popcorn-provider";
import { Job } from "../../types/models";

export interface Clause {
  text: string;
  holds: boolean | null;
  values: Record<string, unknown>;
  name?: string;
  op?: string;
  threshold?: number;
  value?: number | null;
}

export interface NextStep {
  when: string;
  task?: string;
  rerun?: boolean;
  inputs?: Record<string, unknown>;
  advice?: string;
}

export interface Verdict {
  task: string;
  outcome: string | null;
  tone?: "good" | "caution" | "bad" | "unknown";
  because?: string;
  basis?: string;
  clauses?: Clause[];
  results?: Record<string, unknown>;
  meanings?: Record<string, string>;
  missing?: string[];
  next?: NextStep[];
  note?: string;
  judgement_version?: string | null;
  status?: string;
}

const TONE_COLOUR = {
  good: "success",
  caution: "warning",
  bad: "error",
  unknown: "default",
} as const;

const words = (text?: string) => (text || "").split(/\s+/).join(" ").trim();

const showValue = (value: unknown) => {
  if (value === null || value === undefined) return "—";
  if (typeof value === "number") return Number.isInteger(value) ? String(value) : value.toPrecision(4).replace(/\.?0+$/, "");
  if (typeof value === "object") return JSON.stringify(value);
  return String(value);
};

/**
 * One threshold the verdict used, drawn as a bar: the threshold marked, the
 * value as a dot, the side of the line where the condition holds shaded.
 * The general form of the meter the REFMAC verdict draws for one score.
 */
export const ThresholdGauge = ({ clause }: { clause: Clause }) => {
  const { threshold, value, op } = clause;
  if (threshold === undefined || value === null || value === undefined || typeof value !== "number")
    return null;
  const span = Math.max(Math.abs(threshold), Math.abs(value), 1e-6);
  const lo = Math.min(threshold, value) - 0.35 * span;
  const hi = Math.max(threshold, value) + 0.35 * span;
  const x = (v: number) => (100 * (v - lo)) / (hi - lo);
  const holdsBelow = op === "<" || op === "<=";
  const shade = clause.holds ? "rgba(76,175,80,0.25)" : "rgba(244,67,54,0.2)";
  return (
    <svg
      width="100%"
      height="22"
      viewBox="0 0 100 22"
      preserveAspectRatio="none"
      role="img"
      aria-label={`${clause.name} ${showValue(value)} against ${op} ${threshold}`}
    >
      <rect x="0" y="9" width="100" height="4" fill="rgba(128,128,128,0.25)" />
      <rect
        x={holdsBelow ? 0 : x(threshold)}
        y="9"
        width={holdsBelow ? x(threshold) : 100 - x(threshold)}
        height="4"
        fill={shade}
      />
      <line x1={x(threshold)} x2={x(threshold)} y1="3" y2="19" stroke="currentColor" strokeWidth="0.6" />
      <circle cx={x(value)} cy="11" r="3" fill={clause.holds ? "#2e7d32" : "#c62828"} />
    </svg>
  );
};

const ClauseRow = ({ clause }: { clause: Clause }) => {
  const icon =
    clause.holds === true ? (
      <CheckCircle color="success" fontSize="small" />
    ) : clause.holds === false ? (
      <Cancel color="error" fontSize="small" />
    ) : (
      <HelpOutline color="disabled" fontSize="small" />
    );
  const values = Object.entries(clause.values || {})
    .map(([name, value]) => `${name} = ${showValue(value)}`)
    .join(", ");
  return (
    <Stack direction="row" spacing={1} alignItems="center" sx={{ py: 0.5 }}>
      {icon}
      <Box sx={{ minWidth: 0, flex: 1 }}>
        <Typography variant="body2" sx={{ fontFamily: "monospace" }}>
          {clause.text}
        </Typography>
        <Typography variant="caption" color="text.secondary">
          {clause.holds === null ? `not known: ${values}` : values}
        </Typography>
      </Box>
      <Box sx={{ width: 160, flexShrink: 0 }}>
        <ThresholdGauge clause={clause} />
      </Box>
    </Stack>
  );
};

const StepCard = ({
  step,
  index,
  job,
  title,
  onApplied,
}: {
  step: NextStep;
  index: number;
  job: Job;
  title: string;
  onApplied: (newJob: Job) => void;
}) => {
  const api = useApi();
  const { setMessage } = usePopcorn();
  const [busy, setBusy] = useState(false);
  const apply = useCallback(async () => {
    setBusy(true);
    try {
      const result: any = await api.post(`jobs/${job.id}/apply_next/`, { index });
      if (!result?.success) {
        setMessage(`Could not make the job: ${result?.error || "unknown error"}`, "error");
        return;
      }
      const failed = (result.data.inputs || []).filter((i: any) => !i.ok);
      setMessage(
        failed.length
          ? `Job ${result.data.job.number} made; not set: ${failed.map((i: any) => i.name).join(", ")}`
          : `Job ${result.data.job.number} made, with ${result.data.inputs.length} input(s) set. Check it, then run.`,
        failed.length ? "warning" : "success"
      );
      onApplied(result.data.job);
    } catch (error) {
      setMessage(`Could not make the job: ${error instanceof Error ? error.message : String(error)}`, "error");
    } finally {
      setBusy(false);
    }
  }, [api, job.id, index, setMessage, onApplied]);

  const inputs = Object.entries(step.inputs || {});
  return (
    <Card variant="outlined">
      <CardContent sx={{ pb: 1 }}>
        <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 0.5 }}>
          {step.rerun ? <Replay fontSize="small" /> : <PlayArrow fontSize="small" />}
          <Typography variant="subtitle2">{step.task ? title : "Advice"}</Typography>
        </Stack>
        {step.advice && (
          <Typography variant="body2" sx={{ mb: inputs.length ? 1 : 0 }}>
            {words(step.advice)}
          </Typography>
        )}
        {inputs.length > 0 && (
          <Stack direction="row" spacing={0.5} useFlexGap flexWrap="wrap">
            {inputs.map(([name, value]) => (
              <Chip key={name} size="small" variant="outlined" label={`${name} = ${showValue(value)}`} />
            ))}
          </Stack>
        )}
      </CardContent>
      {step.task && (
        <CardActions sx={{ pt: 0 }}>
          <Button
            size="small"
            variant="contained"
            disabled={busy}
            startIcon={busy ? <CircularProgress size={14} /> : step.rerun ? <Replay /> : <PlayArrow />}
            onClick={apply}
          >
            {step.rerun ? "Rerun with these changes" : "Create job"}
          </Button>
        </CardActions>
      )}
    </Card>
  );
};

export const JudgementPanel = ({ job, onJobCreated }: { job: Job; onJobCreated?: () => void }) => {
  const api = useApi();
  const router = useRouter();
  const { data: response, isLoading } = api.get_endpoint<any>({
    type: "jobs",
    id: job.id,
    endpoint: "judgement",
  });
  const { data: taskInfo } = api.get<any>(`agent/tasks/${job.task_name}/`);
  const { data: taskLookup } = api.get<any>(`task_lookup/`);

  const onApplied = useCallback(
    (newJob: Job) => {
      onJobCreated?.();
      router.push(`/ccp4i2/project/${job.project}/job/${newJob.id}`);
    },
    [onJobCreated, router, job.project]
  );

  if (isLoading) return <LinearProgress />;
  const verdict: Verdict | undefined = response?.success ? response.data : undefined;
  if (!verdict) return <Alert severity="error">The judgement could not be read: {response?.error || "no response"}</Alert>;
  if (!verdict.outcome && !verdict.results)
    return <Alert severity="info">{verdict.note || "This task has no judgement."}</Alert>;

  const judgementFile = taskInfo?.data?.judgement || taskInfo?.judgement;
  const traps: string[] = (judgementFile?.traps || []).map((t: unknown) => words(String(t)));
  const used = new Set((verdict.clauses || []).flatMap((c) => Object.keys(c.values || {})));
  const titleOf = (step: NextStep) =>
    step.rerun
      ? `Rerun ${job.task_name}`
      : taskLookup?.[step.task || ""]?.TASKTITLE || step.task || "";

  return (
    <Stack spacing={2} sx={{ p: 2 }}>
      <Card variant="outlined">
        <CardContent>
          <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 1 }}>
            <Chip
              label={verdict.outcome ?? "not judged"}
              color={TONE_COLOUR[verdict.tone || "unknown"]}
              sx={{ fontWeight: 600 }}
            />
            {verdict.because && (
              <Typography variant="body2" sx={{ fontFamily: "monospace" }} color="text.secondary">
                because {verdict.because}
              </Typography>
            )}
          </Stack>
          {verdict.basis && <Typography variant="body2">{words(verdict.basis)}</Typography>}
          {(verdict.clauses || []).length > 0 && (
            <Box sx={{ mt: 1.5 }}>
              {verdict.clauses!.map((clause, i) => (
                <ClauseRow key={i} clause={clause} />
              ))}
            </Box>
          )}
          {verdict.note && (
            <Alert severity="info" sx={{ mt: 1.5 }}>
              {verdict.note}
              {verdict.judgement_version && ` (judgement version ${verdict.judgement_version})`}
            </Alert>
          )}
        </CardContent>
      </Card>

      {(verdict.next || []).length > 0 && (
        <Box>
          <Typography variant="h6" sx={{ mb: 1 }}>
            What next
          </Typography>
          <Stack spacing={1}>
            {verdict.next!.map((step, i) => (
              <StepCard key={i} step={step} index={i} job={job} title={titleOf(step)} onApplied={onApplied} />
            ))}
          </Stack>
        </Box>
      )}

      <Box>
        <Typography variant="h6" sx={{ mb: 1 }}>
          The numbers
        </Typography>
        <Table size="small">
          <TableHead>
            <TableRow>
              <TableCell>Result</TableCell>
              <TableCell>Value</TableCell>
              <TableCell>What it is</TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {Object.entries(verdict.results || {}).map(([name, value]) => (
              <TableRow key={name} selected={used.has(name)}>
                <TableCell sx={{ fontFamily: "monospace", whiteSpace: "nowrap" }}>{name}</TableCell>
                <TableCell sx={{ whiteSpace: "nowrap" }}>
                  {value === null || value === undefined ? (
                    <Typography variant="caption" color="text.secondary">
                      {verdict.missing?.includes(name) ? "not read" : "absent"}
                    </Typography>
                  ) : (
                    showValue(value)
                  )}
                </TableCell>
                <TableCell>
                  <Tooltip title={verdict.meanings?.[name] || ""}>
                    <Typography variant="caption" sx={{ display: "-webkit-box", WebkitLineClamp: 2, WebkitBoxOrient: "vertical", overflow: "hidden" }}>
                      {verdict.meanings?.[name]}
                    </Typography>
                  </Tooltip>
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </Box>

      {traps.length > 0 && (
        <Box>
          <Typography variant="h6" sx={{ mb: 1 }}>
            Traps
          </Typography>
          <Stack spacing={1}>
            {traps.map((trap, i) => (
              <Alert key={i} severity="warning" variant="outlined">
                {trap}
              </Alert>
            ))}
          </Stack>
        </Box>
      )}
    </Stack>
  );
};
