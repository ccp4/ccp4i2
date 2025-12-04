import {
  Button,
  Stack,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Menu as MuiMenu,
  MenuItem,
  IconButton,
} from "@mui/material";
import {
  Code,
  ContentCopy,
  Description,
  DirectionsRun,
  Help,
  Menu,
  MenuBook,
  SystemUpdateAlt,
  MoreVert,
} from "@mui/icons-material";
import React, { useRef, useEffect, useState, useMemo } from "react";
import { useApi } from "../api";
import { Job } from "../types/models";
import { useCCP4i2Window } from "../app-context";
import { apiGet } from "../api-fetch";
import { useRouter } from "next/navigation";
import { HelpIframe } from "./help_iframe";
import { usePopcorn } from "../providers/popcorn-provider";
import { useRunCheck } from "../providers/run-check-provider";
import { useJobTab } from "../providers/job-tab-provider";
import { I2RunDialog } from "./i2run-dialog";

interface ToolbarButton {
  label: string;
  icon: React.ReactNode;
  onClick: () => void;
  disabled?: boolean;
  show: boolean;
}

export default function ToolBar() {
  const panelRef = useRef<HTMLDivElement>(null);
  const [panelWidth, setPanelWidth] = useState<number>(0);
  const { projectId, jobId } = useCCP4i2Window();
  const api = useApi();
  const { data: job, mutate: mutateJob } = api.get_endpoint<Job>({
    type: "jobs",
    id: jobId,
    endpoint: "",
  });
  const { mutate: mutateJobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: projectId,
    endpoint: "jobs",
  });
  const router = useRouter();
  const [showHelpPanel, setShowHelpPanel] = useState(false);
  const [showI2RunDialog, setShowI2RunDialog] = useState(false);
  const [i2RunCommand, setI2RunCommand] = useState<string>("");
  const { setMessage } = usePopcorn();
  const { confirmTaskRun } = useRunCheck();
  const { setJobTabValue } = useJobTab();
  const [menuAnchor, setMenuAnchor] = useState<null | HTMLElement>(null);

  useEffect(() => {
    if (!panelRef.current) return;
    const observer = new window.ResizeObserver((entries) => {
      for (let entry of entries) {
        setPanelWidth(entry.contentRect.width);
      }
    });
    observer.observe(panelRef.current);
    return () => observer.disconnect();
  }, []);

  const handleClone = async () => {
    if (job) {
      try {
        const cloneResult: any = await api.post(`jobs/${job?.id}/clone/`);
        if (cloneResult?.success === false) {
          setMessage(`Failed to clone job: ${cloneResult?.error || "Unknown error"}`, "error");
          return;
        }
        if (cloneResult?.id) {
          mutateJob();
          mutateJobs();
          router.push(`/project/${projectId}/job/${cloneResult.id}`);
        }
      } catch (error) {
        setMessage(`Error cloning job: ${error instanceof Error ? error.message : String(error)}`, "error");
      }
    }
  };

  const handleRun = async () => {
    if (job) {
      const confirmed = await confirmTaskRun(job.id);
      if (!confirmed) return;
      try {
        const runResult: any = await api.post(`jobs/${job.id}/run/`);
        if (runResult?.success === false) {
          setMessage(`Failed to run job: ${runResult?.error || "Unknown error"}`, "error");
          return;
        }
        setMessage(`Submitted job ${runResult?.number}: ${runResult?.task_name}`, "success");
        if (runResult?.id) {
          setTimeout(() => {
            mutateJob();
            mutateJobs();
          }, 1000);
        }
      } catch (error) {
        setMessage(`Error running job: ${error instanceof Error ? error.message : String(error)}`, "error");
      }
    }
  };

  const handleI2Run = async () => {
    if (job) {
      const result: { status: string; command: string } = await apiGet(
        `jobs/${job.id}/i2run_command`
      );
      if (result?.command) {
        navigator.clipboard.writeText(result.command);
        setMessage("i2run command copied to clipboard");
        setI2RunCommand(result.command);
        setShowI2RunDialog(true);
      }
    }
  };

  const handleLog = () => setJobTabValue(10);

  // Button definitions with breakpoints
  const toolbarButtons: ToolbarButton[] = useMemo(
    () => [
      {
        label: "Task menu",
        icon: <Menu />,
        onClick: () => router.push(`/project/${projectId}`),
        show: true,
      },
      {
        label: "Run",
        icon: <DirectionsRun />,
        onClick: handleRun,
        disabled: job?.status !== 1,
        show: true,
      },
      {
        label: "Clone job",
        icon: <ContentCopy />,
        onClick: handleClone,
        show: panelWidth > 550,
      },
      {
        label: "Help",
        icon: <Help />,
        onClick: () => {
          if (window?.open) {
            window.open(
              `https://ccp4i2.gitlab.io/rstdocs/tasks/${job?.task_name}/index.html`
            );
          }
        },
        show: panelWidth > 650,
      },
      {
        label: "Bibliography",
        icon: <MenuBook />,
        onClick: () => {},
        show: panelWidth > 750,
      },
      {
        label: "Export MTZ",
        icon: <SystemUpdateAlt />,
        onClick: () => {},
        disabled: job?.status !== 6,
        show: panelWidth > 950,
      },
      {
        label: "Show log file",
        icon: <Description />,
        onClick: handleLog,
        disabled: job?.status !== 6,
        show: panelWidth > 1100,
      },
      {
        label: "i2run command",
        icon: <Code />,
        onClick: handleI2Run,
        show: panelWidth > 1200,
      },
    ],
    [panelWidth, job, projectId, router]
  );

  const visibleButtons = toolbarButtons.filter((btn) => btn.show);
  const hiddenButtons = toolbarButtons.filter((btn) => !btn.show);

  return (
    <>
      <div ref={panelRef}>
        <Stack
          direction="row"
          spacing={2}
          useFlexGap
          sx={{ flexWrap: "wrap", justifyContent: "center", px: 2, mb: 1 }}
        >
          {visibleButtons.map((btn) => (
            <Button
              key={btn.label}
              variant="outlined"
              startIcon={btn.icon}
              onClick={btn.onClick}
              disabled={btn.disabled}
            >
              {btn.label}
            </Button>
          ))}
          {hiddenButtons.length > 0 && (
            <>
              <IconButton
                aria-label="More"
                onClick={(e) => setMenuAnchor(e.currentTarget)}
                sx={{ ml: 1 }}
              >
                <MoreVert />
              </IconButton>
              <MuiMenu
                anchorEl={menuAnchor}
                open={Boolean(menuAnchor)}
                onClose={() => setMenuAnchor(null)}
              >
                {hiddenButtons.map((btn) => (
                  <MenuItem
                    key={btn.label}
                    onClick={() => {
                      btn.onClick();
                      setMenuAnchor(null);
                    }}
                    disabled={btn.disabled}
                  >
                    {btn.icon}
                    <span style={{ marginLeft: 8 }}>{btn.label}</span>
                  </MenuItem>
                ))}
              </MuiMenu>
            </>
          )}
          <HelpIframe
            url={`/help/html/tasks/${job?.task_name}/index.html`}
            open={showHelpPanel}
            handleClose={() => setShowHelpPanel(false)}
          />
        </Stack>
      </div>
      <I2RunDialog
        open={showI2RunDialog}
        command={i2RunCommand}
        onClose={() => setShowI2RunDialog(false)}
      />
    </>
  );
}
