"use client";
import { useCallback, useEffect, useMemo, useState } from "react";
import {
  Box,
  Button,
  Container,
  LinearProgress,
  Stack,
  Tab,
  Tabs,
} from "@mui/material";
import { Save as SaveIcon, Restore as RestoreIcon } from "@mui/icons-material";
import { Editor } from "@monaco-editor/react";
import { JobHeader } from "../components/job-header";
import { CCP4i2ReportXMLView } from "../components/report/CCP4i2ReportXMLView";
import { CCP4i2WhatNext } from "../components/report/CCP4i2WhatNext";
import { useCCP4i2Window } from "../app-context";
import { TaskContainer } from "../components/task/task-interfaces/task-container";
import { prettifyXml, useJob, usePrevious, useProject } from "../utils";
import ToolBar from "../components/tool-bar";
import { JobCommentEditor } from "../components/job-comment-editor";
import { JobMenu } from "../providers/job-context-menu";
import { JobDirectoryView } from "../components/job-directory-view";
import { useApi } from "../api";
import { apiPut } from "../api-fetch";
import $ from "jquery";
import Diagnostic from "../components/diagnostic";
import { JobLogViewer } from "../components/job-log-viewer";
import { TaskProvider } from "../providers/task-provider";
import { ValidationViewer } from "../components/validation-viewer";
import { useJobTab } from "../providers/job-tab-provider";
import { JobStatus } from "../types/models";
import {
  TAB,
  REPORT_STATUSES,
  landingTab,
  visibleTabs,
} from "./job-view-tabs";
import { useTheme } from "../theme/theme-provider";
import { useIsJobEffectivelyActive } from "../providers/recently-started-jobs-context";

export interface JobViewProps {
  jobid: number;
}
export const JobView: React.FC<JobViewProps> = ({ jobid }) => {
  const { devMode, jobId, setJobId, projectId, setProjectId } =
    useCCP4i2Window();
  const api = useApi();

  const { project, jobs, mutateJobs } = useProject(projectId);

  // Get detailed job data (params_xml, container, etc.)
  const {
    job,
    params_xml,
    mutateParams_xml,
    validation,
    diagnostic_xml,
    mutateDiagnosticXml,
    def_xml,
    container,
  } = useJob(jobid);

  const { mode } = useTheme();

  // Status from job_tree (shared SWR key, polled by ClassicJobsList every 3-30s).
  // Fall back to useJob's status during initial load.
  const jobFromTree = useMemo(() => jobs?.find((j) => j.id === jobid), [jobs, jobid]);
  const currentStatus = jobFromTree?.status ?? job?.status;
  const isJobActive = useIsJobEffectivelyActive(jobid, currentStatus);

  // Merge: prefer job_tree status (more current), with useJob for detail fields
  const jobWithCurrentStatus = useMemo(() => {
    if (!job) return jobFromTree;
    if (jobFromTree && jobFromTree.status !== job.status) {
      return { ...job, status: jobFromTree.status };
    }
    return job;
  }, [job, jobFromTree]);

  const previousJob = usePrevious(job);
  const previousStatus = usePrevious(currentStatus);
  const { jobTabValue: rawTabValue, setJobTabValue: setTabValue } = useJobTab();

  // State for editable params XML (only used for pending jobs)
  const [editedParamsXml, setEditedParamsXml] = useState<string | null>(null);
  const [isSavingXml, setIsSavingXml] = useState(false);
  const [xmlSaveError, setXmlSaveError] = useState<string | null>(null);

  // Track if the XML has been modified
  const isXmlModified = useMemo(() => {
    return editedParamsXml !== null && editedParamsXml !== params_xml;
  }, [editedParamsXml, params_xml]);

  // Check if job is pending (editable)
  const isPending = currentStatus === 1;

  // Reset edited XML when job changes or params_xml changes
  useEffect(() => {
    setEditedParamsXml(null);
    setXmlSaveError(null);
  }, [jobid, params_xml]);

  // Handle XML editor changes
  const handleXmlChange = useCallback((value: string | undefined) => {
    if (value !== undefined) {
      setEditedParamsXml(value);
      setXmlSaveError(null);
    }
  }, []);

  // Save XML to server
  const handleSaveXml = useCallback(async () => {
    if (!editedParamsXml || !job) return;

    setIsSavingXml(true);
    setXmlSaveError(null);

    try {
      const response = await apiPut(`jobs/${job.id}/params_xml/`, {
        xml: editedParamsXml,
      });

      if (response.success) {
        // Clear edited state on successful save
        setEditedParamsXml(null);
        // Trigger a refresh of params_xml from server
        mutateParams_xml();
      } else {
        setXmlSaveError(response.error || "Failed to save XML");
      }
    } catch (error: any) {
      console.error("[JobView] Failed to save XML:", error);
      setXmlSaveError(error.message || "Failed to save XML");
    } finally {
      setIsSavingXml(false);
    }
  }, [editedParamsXml, job, mutateParams_xml]);

  // Reset XML to original
  const handleResetXml = useCallback(() => {
    setEditedParamsXml(null);
    setXmlSaveError(null);
  }, []);

  // This is for the raw XML editor view (tabValue == 2)
  // Uses same key as CCP4i2ReportXMLView so SWR deduplicates - keep polling logic consistent
  // Don't fetch report_xml for pending jobs - no meaningful report exists yet
  const shouldFetchReport = currentStatus !== undefined && currentStatus > 1;
  const { data: report_xml_json } = api.jobReportXml(
    shouldFetchReport ? job?.id : null,
    isJobActive
  );

  const report_xml: XMLDocument | null = useMemo(() => {
    if (!report_xml_json) return null;
    // Handle both wrapped response {success: true, data: {xml: ...}} and direct {xml: ...}
    const xmlString = report_xml_json.data?.xml || report_xml_json.xml;
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
      // Land on the tab that answers the question this job raises. Written
      // against JobStatus rather than bare numbers: the comment here used to
      // read "4=Failed, 5=Unsatisfactory", which is 4=Interrupted, 5=Failed,
      // 10=Unsatisfactory, and Interrupted and Unsatisfactory jobs were sent
      // to a Diagnostics tab that is shown only for Failed --- so the clamp
      // above bounced them to the task interface.
      if (job && job != previousJob && currentStatus !== undefined) {
        setTabValue(landingTab(currentStatus));
      }
    };
    asyncFunc();
  }, [job, setJobId, currentStatus]);

  // Invalidate diagnostic_xml cache when job transitions to a terminal state
  // so that the Diagnostics tab shows fresh data immediately
  useEffect(() => {
    if (previousStatus === undefined || currentStatus === undefined) return;
    const wasActive = [2, 3, 7].includes(previousStatus);
    const isTerminal = [4, 5, 6].includes(currentStatus);
    if (wasActive && isTerminal) {
      mutateDiagnosticXml();
    }
  }, [currentStatus, previousStatus, mutateDiagnosticXml]);

  // Clamp synchronously so MUI never receives a value for a hidden tab
  const status = jobWithCurrentStatus?.status;
  const tabValue = useMemo(() => {
    const visible = visibleTabs(status, devMode);
    return visible.has(rawTabValue) ? rawTabValue : TAB.TASK_INTERFACE;
  }, [devMode, status, rawTabValue]);

  return !project || !jobs || !jobWithCurrentStatus ? (
    <LinearProgress />
  ) : (
    <>
      <Stack sx={{ height: "100%" }}>
        <ToolBar />
        <JobHeader job={jobWithCurrentStatus} mutateJobs={mutateJobs} />
        {/* Kept in step with the `visible` set above: a tab rendered here and
            missing there (or the reverse) is how a job lands on a tab that is
            not on screen. */}
        <Tabs value={tabValue} onChange={handleTabChange} variant="fullWidth">
          <Tab value={TAB.TASK_INTERFACE} label="Task interface" />
          {devMode && <Tab value={TAB.PARAMS_XML} label="Params as xml" />}
          {devMode && <Tab value={TAB.REPORT_XML} label="Report as xml" />}
          {(devMode ||
            REPORT_STATUSES.includes(jobWithCurrentStatus.status)) && (
            <Tab value={TAB.REPORT} label="Report" />
          )}
          {(devMode || jobWithCurrentStatus.status === JobStatus.FAILED) && (
            <Tab value={TAB.DIAGNOSTICS} label="Diagnostics" />
          )}
          {devMode && <Tab value={TAB.DEF_XML} label="Def xml" />}
          {(devMode || jobWithCurrentStatus.status === JobStatus.PENDING) && (
            <Tab value={TAB.VALIDATION} label="Validation" />
          )}
          {devMode && <Tab value={TAB.JOB_CONTAINER} label="Job container" />}
          <Tab value={TAB.COMMENTS} label="Comments" />
          <Tab value={TAB.DIRECTORY} label="Directory" />
          <Tab value={TAB.LOGS} label="Logs" />
        </Tabs>
        <Box
          sx={{
            flex: "auto",
            minHeight: 0,
            overflowY: "auto",
            scrollbarWidth: "thin",
            pb: 4,
          }}
        >
          {tabValue == 0 && (
            <TaskProvider>
              <TaskContainer jobId={jobid} />
            </TaskProvider>
          )}
          {devMode && tabValue == 1 && params_xml && (
            <Box sx={{ height: "100%", display: "flex", flexDirection: "column" }}>
              {/* Toolbar with Save/Reset buttons for pending jobs */}
              {isPending && (
                <Stack
                  direction="row"
                  spacing={1}
                  sx={{ p: 1, borderBottom: 1, borderColor: "divider" }}
                >
                  <Button
                    variant="contained"
                    size="small"
                    startIcon={<SaveIcon />}
                    onClick={handleSaveXml}
                    disabled={!isXmlModified || isSavingXml}
                  >
                    {isSavingXml ? "Saving..." : "Save"}
                  </Button>
                  <Button
                    variant="outlined"
                    size="small"
                    startIcon={<RestoreIcon />}
                    onClick={handleResetXml}
                    disabled={!isXmlModified || isSavingXml}
                  >
                    Reset
                  </Button>
                  {isXmlModified && (
                    <Box
                      component="span"
                      sx={{ color: "warning.main", alignSelf: "center", ml: 1 }}
                    >
                      (unsaved changes)
                    </Box>
                  )}
                  {xmlSaveError && (
                    <Box
                      component="span"
                      sx={{ color: "error.main", alignSelf: "center", ml: 1 }}
                    >
                      {xmlSaveError}
                    </Box>
                  )}
                </Stack>
              )}
              <Box sx={{ flex: 1 }}>
                <Editor
                  height="100%"
                  value={editedParamsXml ?? params_xml}
                  language="xml"
                  theme={mode === "dark" ? "vs-dark" : "light"}
                  onChange={isPending ? handleXmlChange : undefined}
                  options={{
                    readOnly: !isPending,
                    minimap: { enabled: false },
                  }}
                />
              </Box>
            </Box>
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
        {tabValue == 3 && jobid && <CCP4i2WhatNext />}
        <JobMenu />
      </Stack>
    </>
  );
};
