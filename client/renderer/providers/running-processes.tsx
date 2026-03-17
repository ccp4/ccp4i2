import {
  createContext,
  PropsWithChildren,
  useCallback,
  useContext,
  useMemo,
  useState,
} from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  Table,
  TableHead,
  TableRow,
  TableBody,
  TableCell,
  Avatar,
  IconButton,
  Tooltip,
} from "@mui/material";
import CancelIcon from "@mui/icons-material/Cancel";
import { useApi } from "../api";
import { apiPost } from "../api-fetch";
import { Project } from "../types/models";

interface RunningProcessesProps {
  jobsAndProcessesDialogOpen: boolean;
  setJobsAndProcessesDialogOpen: (open: boolean) => void;
}

export const RunningProcessesContext = createContext<RunningProcessesProps>({
  jobsAndProcessesDialogOpen: false,
  setJobsAndProcessesDialogOpen: () => {},
});

interface RunningProcess {
  project: string;
  job_id: number;
  job_task_name: string;
  job_uuid: string;
  job_number: string;
  status: string;
  pid: number | null;
  cpu_percent: number | null;
  memory_usage_bytes: string | null;
}
export const RunningProcessesProvider: React.FC<PropsWithChildren> = (
  props
) => {
  const [jobsAndProcessesDialogOpen, setJobsAndProcessesDialogOpen] =
    useState<boolean>(false);
  const [cancellingJobs, setCancellingJobs] = useState<Set<number>>(new Set());
  const api = useApi();
  // Only poll active_jobs when the dialog is open
  const { data: activeJobs, mutate: mutateActiveJobs } = api.get<any>(
    jobsAndProcessesDialogOpen ? "active_jobs/" : null,
    jobsAndProcessesDialogOpen ? 5000 : 0
  );
  const { data: projects } = api.get<Project[]>("projects/");
  // Handle new API response format: {success: true, data: {active_jobs: [...]}}
  const runningProcesses =
    activeJobs?.success ? activeJobs.data?.active_jobs : [];

  const onSelectRow = useCallback(
    (row: RunningProcess) => {
      const project = projects?.find((project) => project.name === row.project);
      if (project) {
        const url = `/ccp4i2/project/${project.id}/job/${row.job_id}`;
        window.open(url, "_blank");
      }
    },
    [projects]
  );

  const onCancelJob = useCallback(
    async (e: React.MouseEvent, jobId: number) => {
      e.stopPropagation(); // Don't trigger row click
      setCancellingJobs((prev) => new Set(prev).add(jobId));
      try {
        await apiPost(`jobs/${jobId}/cancel/`);
        mutateActiveJobs(); // Refresh the list
      } catch (err) {
        console.error("Failed to cancel job:", err);
      } finally {
        setCancellingJobs((prev) => {
          const next = new Set(prev);
          next.delete(jobId);
          return next;
        });
      }
    },
    [mutateActiveJobs]
  );

  const contextValue = useMemo(
    () => ({ jobsAndProcessesDialogOpen, setJobsAndProcessesDialogOpen }),
    [jobsAndProcessesDialogOpen]
  );

  return (
    <RunningProcessesContext.Provider value={contextValue}>
      {props.children}
      <Dialog
        open={jobsAndProcessesDialogOpen}
        onClose={() => {
          setJobsAndProcessesDialogOpen(false);
        }}
        fullWidth
        maxWidth="lg"
      >
        <DialogTitle>Running Processes</DialogTitle>
        <DialogContent>
          {runningProcesses && (
            <Table>
              <TableHead>
                <TableRow>
                  <TableCell variant="head">Project Name</TableCell>
                  <TableCell variant="head">Id</TableCell>
                  <TableCell variant="head">Task name</TableCell>
                  <TableCell variant="head">Status</TableCell>
                  <TableCell variant="head">number</TableCell>
                  <TableCell variant="head">Pid</TableCell>
                  <TableCell variant="head">%CPU</TableCell>
                  <TableCell variant="head">Memory usage (MB)</TableCell>
                  <TableCell variant="head">Actions</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {runningProcesses.map((process: RunningProcess) => (
                  <TableRow
                    key={process.job_id}
                    onClick={() => onSelectRow(process)}
                    sx={{ cursor: "pointer" }}
                  >
                    <TableCell variant="body">{process.project}</TableCell>
                    <TableCell variant="body">{process.job_id}</TableCell>
                    <TableCell variant="body">
                      <Avatar
                        src={`/svgicons/${process.job_task_name}.svg`}
                        alt={process.job_task_name}
                        imgProps={{
                          onError: (e: any) => {
                            e.target.src = `/qticons/${process.job_task_name}.png`;
                          },
                        }}
                      >
                        {process.job_task_name?.[0]?.toUpperCase()}
                      </Avatar>{" "}
                      {process.job_task_name}
                    </TableCell>
                    <TableCell variant="body">{process.status}</TableCell>
                    <TableCell variant="body">{process.job_number}</TableCell>
                    <TableCell variant="body">
                      {process.pid ?? "-"}
                    </TableCell>
                    <TableCell variant="body">
                      {process.cpu_percent != null
                        ? process.cpu_percent
                        : "-"}
                    </TableCell>
                    <TableCell variant="body">
                      {process.memory_usage_bytes != null
                        ? (
                            Number(process.memory_usage_bytes) /
                            (1024 * 1024)
                          ).toFixed(2)
                        : "-"}
                    </TableCell>
                    <TableCell variant="body">
                      <Tooltip title="Cancel job">
                        <span>
                          <IconButton
                            size="small"
                            color="error"
                            disabled={cancellingJobs.has(process.job_id)}
                            onClick={(e) => onCancelJob(e, process.job_id)}
                          >
                            <CancelIcon />
                          </IconButton>
                        </span>
                      </Tooltip>
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          )}
        </DialogContent>
      </Dialog>
    </RunningProcessesContext.Provider>
  );
};

export const useRunningProcesses = () => {
  const context = useContext(RunningProcessesContext);
  if (!context) {
    throw new Error(
      "useRunningProcesses must be used within a RunningProcessesProvider"
    );
  }
  return context;
};
