import { useCallback, useState } from "react";
import {
  Alert,
  Box,
  Button,
  Paper,
  Stack,
  Tooltip,
  Typography,
} from "@mui/material";
import { PlaylistAdd as FillIcon } from "@mui/icons-material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { useApi } from "../../../api";
import { apiPost } from "../../../api-fetch";

/**
 * pandda_campaign — run PanDDA 2 over a declared list of datasets.
 *
 * The task never asks which datasets a campaign has (design note §3.1): the
 * DATASETS list is the whole input, the record of what was submitted, and
 * the lever for batching. Campaign-awareness lives here, in job
 * construction: "Fill from campaign" asks the plugin (through the generic
 * object_method endpoint) for every member with a finished dimple job and
 * appends the ones not already listed. The page knows whether this project
 * is a campaign parent and why nothing can be filled, so the button is
 * disabled with that reason rather than hidden.
 */

interface Candidate {
  label: string;
  project_uuid: string | null;
  listed: boolean;
}
interface Skipped {
  project: string;
  reason: string;
}
interface Preview {
  campaigns: { name: string; members: number }[];
  candidates: Candidate[];
  skipped: Skipped[];
  reason: string | null;
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const api = useApi();
  const { mutateContainer, mutateValidation, useTaskItem } = useJob(job.id);
  const { value: runMode } = useTaskItem("RUN_MODE");
  const [busy, setBusy] = useState(false);
  const [message, setMessage] = useState<string | null>(null);

  const editable = job.status === 1 || job.status === 0; // pending / unknown
  const { data: previewResponse, mutate: mutatePreview } =
    api.objectMethod<{ data?: { result?: Preview } }>(
      job.id,
      "pandda_campaign",
      "campaignCandidates",
      {},
      [job.status],
      editable,
    );
  const preview = previewResponse?.data?.result;
  const toAdd = (preview?.candidates ?? []).filter((c) => !c.listed);
  const campaignName = preview?.campaigns?.[0]?.name;

  let disabledReason: string | null = null;
  if (!editable) disabledReason = "The job has already been run.";
  else if (!preview) disabledReason = "Checking the campaign...";
  else if (preview.reason) disabledReason = preview.reason;
  else if (toAdd.length === 0) disabledReason = "Every dataset with a finished dimple job is already listed.";

  const fill = useCallback(async () => {
    setBusy(true);
    setMessage(null);
    try {
      const response: any = await apiPost(`jobs/${job.id}/object_method/`, {
        object_path: "pandda_campaign",
        method_name: "fillDatasetsFromCampaign",
        args: [],
        kwargs: {},
      });
      const result = response?.data?.result;
      if (result?.success === false) {
        setMessage(result.error ?? "Could not fill the list.");
      } else {
        const added: string[] = result?.data?.added ?? [];
        setMessage(
          added.length
            ? `Added ${added.length} dataset${added.length === 1 ? "" : "s"}: ${added.join(", ")}`
            : "Nothing to add.",
        );
      }
      // The list changed underneath the cached validation: refresh it too, or
      // the confirm dialog keeps saying "add at least one dataset".
      await mutateContainer();
      await mutateValidation();
      await mutatePreview();
    } catch (err: any) {
      setMessage(err?.message ?? String(err));
    } finally {
      setBusy(false);
    }
  }, [job.id, mutateContainer, mutateValidation, mutatePreview]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Datasets">
        <Paper sx={{ p: 1, mb: 1 }}>
          <Stack direction="row" spacing={2} alignItems="center" flexWrap="wrap">
            <Tooltip title={disabledReason ?? `Add every member of ${campaignName} with a finished dimple job`}>
              <span>
                <Button
                  variant="contained"
                  startIcon={<FillIcon />}
                  disabled={busy || disabledReason !== null}
                  onClick={fill}
                >
                  Fill from campaign
                  {toAdd.length > 0 ? ` (${toAdd.length})` : ""}
                </Button>
              </span>
            </Tooltip>
            {preview?.skipped?.length ? (
              <Typography variant="caption" color="text.secondary">
                {preview.skipped.length} member{preview.skipped.length === 1 ? "" : "s"} without a
                usable dimple result:{" "}
                {preview.skipped.map((s) => `${s.project} (${s.reason})`).join("; ")}
              </Typography>
            ) : null}
          </Stack>
          {message ? (
            <Alert severity="info" sx={{ mt: 1 }} onClose={() => setMessage(null)}>
              {message}
            </Alert>
          ) : null}
          <Typography variant="caption" color="text.secondary" sx={{ display: "block", mt: 1 }}>
            PanDDA needs at least 25 datasets to characterise a ground state. The list as
            submitted is the record of the run; add a dataset the campaign rule missed with +.
          </Typography>
        </Paper>
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Datasets" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement {...props} itemName="DATASETS" qualifiers={{ guiLabel: "Datasets" }} />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
      <CCP4i2Tab label="Run">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "How to run" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement {...props} itemName="RUN_MODE" qualifiers={{ guiLabel: "Run mode" }} />
          <CCP4i2TaskElement
            {...props}
            itemName="MIN_CHARACTERISATION_DATASETS"
            qualifiers={{ guiLabel: "Minimum datasets to characterise a ground state (PanDDA default 25)" }}
          />
          {runMode === "local" ? (
            <Box sx={{ display: "flex", flexDirection: "column", gap: 1 }}>
              <CCP4i2TaskElement {...props} itemName="LOCAL_CPUS" qualifiers={{ guiLabel: "CPUs" }} />
              <CCP4i2TaskElement
                {...props}
                itemName="SCRATCH_DIR"
                qualifiers={{ guiLabel: "Scratch directory (RAY_TMPDIR)" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="PANDDA_EXECUTABLE"
                qualifiers={{ guiLabel: "pandda2.analyse (optional explicit path)" }}
              />
            </Box>
          ) : (
            <Typography variant="body2" color="text.secondary">
              Stage only: the input tree is staged and the job finishes. Run PanDDA elsewhere
              with the command the report shows, then fan out from its pandda2_out with this
              job&apos;s manifest (manage.py pandda_fanout).
            </Typography>
          )}
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};

export default TaskInterface;
