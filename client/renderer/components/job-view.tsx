"use client";
import { useEffect, useMemo } from "react";
import { Box, Container, LinearProgress, Tab, Tabs } from "@mui/material";
import { Editor } from "@monaco-editor/react";
import { JobHeader } from "../components/job-header";
import { CCP4i2ReportXMLView } from "../components/report/CCP4i2ReportXMLView";
import { useCCP4i2Window } from "../app-context";
import { TaskContainer } from "../components/task/task-interfaces/task-container";
import { prettifyXml, useJob, usePrevious, useProject } from "../utils";
import ToolBar from "../components/tool-bar";
import { JobCommentEditor } from "../components/job-comment-editor";
import { JobMenu } from "../providers/job-context-menu";
import { JobDirectoryView } from "../components/job-directory-view";
import useSWR from "swr";
import $ from "jquery";
import Diagnostic from "../components/diagnostic";
import { swrFetcher } from "../api-fetch";
import { JobLogViewer } from "../components/job-log-viewer";
import { TaskProvider } from "../providers/task-provider";
import { ValidationViewer } from "../components/validation-viewer";
import { useRunCheck } from "../providers/run-check-provider";
import { useJobTab } from "../providers/job-tab-provider";
import { useTheme } from "../theme/theme-provider";
import { useIsJobEffectivelyActive } from "../providers/recently-started-jobs-context";

export interface JobViewProps {
  jobid: number;
}
export const JobView: React.FC<JobViewProps> = ({ jobid }) => {
  const { devMode, jobId, setJobId, projectId, setProjectId } =
    useCCP4i2Window();

  // Get project and jobs list first (from job_tree - always up to date)
  const { project, jobs, mutateJobs } = useProject(projectId || 0);

  // Find job status from jobs array (consistent with jobs list icons)
  const jobFromTree = useMemo(() => {
    return jobs?.find((j) => j.id === jobid);
  }, [jobs, jobid]);

  // Get detailed job data (params_xml, container, etc.) - status may lag behind
  const { job, params_xml, validation, diagnostic_xml, def_xml, container } =
    useJob(jobid);

  const {
    setExtraDialogActions,
    setProcessedErrors,
    extraDialogActions,
    processedErrors,
  } = useRunCheck();
  const { mode } = useTheme();

  // Use job_tree status (jobFromTree) for polling - it's always current
  // Fall back to useJob's status if job_tree hasn't loaded yet
  const currentStatus = jobFromTree?.status ?? job?.status;
  const isJobActive = useIsJobEffectivelyActive(jobid, currentStatus);

  // Create merged job with current status from job_tree for consistent UI
  // This ensures JobHeader shows the same status as the jobs list
  const jobWithCurrentStatus = useMemo(() => {
    if (!job) return undefined;
    if (jobFromTree?.status !== undefined && jobFromTree.status !== job.status) {
      return { ...job, status: jobFromTree.status };
    }
    return job;
  }, [job, jobFromTree?.status]);

  // Debug: log status sources
  console.log(`[JobView] jobid=${jobid}, treeStatus=${jobFromTree?.status}, jobStatus=${job?.status}, mergedStatus=${jobWithCurrentStatus?.status}, isJobActive=${isJobActive}`);

  const previousJob = usePrevious(job);
  const { jobTabValue: tabValue, setJobTabValue: setTabValue } = useJobTab();

  // This SWR is for the raw XML editor view (tabValue == 2)
  // Uses same key as CCP4i2ReportXMLView so SWR deduplicates - keep polling logic consistent
  const { data: report_xml_json } = useSWR<any>(
    job ? `/api/proxy/jobs/${job.id}/report_xml/` : null,
    swrFetcher,
    { refreshInterval: isJobActive ? 5000 : 0 }
  );

  const report_xml: XMLDocument | null = useMemo(() => {
    if (!report_xml_json) return null;
    // Handle both wrapped response {success: true, data: {xml: ...}} and direct {xml: ...}
    const xmlString = report_xml_json.data?.xml || report_xml_json.xml;
    console.log("[JobView] report_xml_json:", report_xml_json, "xmlString length:", xmlString?.length);
    if (!xmlString) return null;
    return $.parseXML(xmlString);
  }, [report_xml_json]);

  const handleTabChange = (event: React.SyntheticEvent, value: number) => {
    setTabValue(value);
  };

  useEffect(() => {
    const asyncFunc = async () => {
      if (job && setJobId && job.id !== jobId) {
        setJobId(job.id);
      }
      if (job && setProjectId && job.project !== projectId) {
        setProjectId(job.project);
      }
      // Use currentStatus (from job_tree) for tab selection to be consistent with UI
      // Status mapping: 1=Pending, 2=Queued, 3=Running, 4=Failed, 5=Unsatisfactory, 6=Finished, 7=Running remotely
      // - Pending (1) → Task interface (tab 0)
      // - Queued/Running/Finished/Running remotely (2,3,6,7) → Report (tab 3) with polling
      // - Failed/Unsatisfactory (4,5) → Diagnostics (tab 4)
      if (job && job != previousJob && currentStatus !== undefined) {
        setTabValue(currentStatus == 1 ? 0 : [2, 3, 6, 7].includes(currentStatus) ? 3 : 4);
      }
    };
    asyncFunc();
  }, [job, setJobId, currentStatus]);

  //Here a useEffect that will clear processedErrors and extraDialogActions when job changes
  useEffect(() => {
    if (!setExtraDialogActions || !setProcessedErrors) return;
    if (extraDialogActions) setExtraDialogActions(null);
    if (processedErrors) setProcessedErrors(null);
  }, [jobid, setExtraDialogActions, setProcessedErrors]);

  return !project || !jobs || !jobWithCurrentStatus ? (
    <LinearProgress />
  ) : (
    <>
      <ToolBar />
      <Container>
        <JobHeader job={jobWithCurrentStatus} mutateJobs={mutateJobs} />
        <Tabs value={tabValue} onChange={handleTabChange} variant="fullWidth">
          <Tab value={0} label="Task interface" />
          {devMode && <Tab value={1} label="Params as xml" />}
          {devMode && <Tab value={2} label="Report as xml" />}
          {(devMode || [3, 4, 6, 7, 9, 10].includes(jobWithCurrentStatus.status)) && (
            <Tab value={3} label="Report" />
          )}
          {(devMode || jobWithCurrentStatus.status === 5) && (
            <Tab value={4} label="Diagnostics" />
          )}
          {devMode && <Tab value={5} label="Def xml" />}
          {(devMode || jobWithCurrentStatus.status === 1) && (
            <Tab value={6} label="Validation" />
          )}
          {devMode && <Tab value={7} label="Job container" />}
          <Tab value={8} label="Comments" />
          <Tab value={9} label="Directory" />
          <Tab value={10} label="Logs" />
        </Tabs>
        {/* Scroll container with fixed height calc - TODO: make responsive to header height */}
        <Box
          sx={(theme) => ({
            height: "calc(100vh - 20rem)",
            overflowY: "auto",
            // Theme-aware scrollbar styling
            scrollbarColor: `${theme.palette.action.disabled} transparent`,
            scrollbarWidth: "thin",
            "&::-webkit-scrollbar": {
              width: 8,
            },
            "&::-webkit-scrollbar-track": {
              background: "transparent",
            },
            "&::-webkit-scrollbar-thumb": {
              backgroundColor: theme.palette.action.disabled,
              borderRadius: 4,
            },
          })}
        >
          {tabValue == 0 && (
            <TaskProvider>
              <TaskContainer jobId={jobid} />
            </TaskProvider>
          )}
          {devMode && tabValue == 1 && params_xml && (
            <Editor
              height="100%"
              value={params_xml}
              language="xml"
              theme={mode === "dark" ? "vs-dark" : "light"}
            />
          )}
          {devMode && tabValue == 2 && report_xml && (
            <Editor
              height="100%"
              value={prettifyXml(report_xml)}
              language="xml"
              theme={mode === "dark" ? "vs-dark" : "light"}
            />
          )}
          {tabValue == 3 && jobid && <CCP4i2ReportXMLView />}
          {(devMode || jobWithCurrentStatus.status === 5) && tabValue == 4 && diagnostic_xml && (
            <Diagnostic xmlDocument={diagnostic_xml} />
          )}
          {devMode && tabValue == 5 && def_xml && (
            <Editor
              height="100%"
              value={def_xml}
              language="xml"
              theme={mode === "dark" ? "vs-dark" : "light"}
            />
          )}
          {(devMode || jobWithCurrentStatus.status === 1) && tabValue == 6 && validation && (
            <ValidationViewer job={jobWithCurrentStatus} />
          )}
          {tabValue == 7 && container && (
            <Editor
              height="100%"
              value={JSON.stringify(container.container, null, 2)}
              language="json"
              theme={mode === "dark" ? "vs-dark" : "light"}
            />
          )}
          {tabValue == 8 && container && (
            <JobCommentEditor jobId={jobWithCurrentStatus.id} />
          )}
          {tabValue == 9 && project && (
            <JobDirectoryView job={jobWithCurrentStatus} project={project} />
          )}
          {tabValue == 10 && project && (
            <JobLogViewer job={jobWithCurrentStatus} project={project} />
          )}
        </Box>
        <JobMenu />
      </Container>
    </>
  );
};
