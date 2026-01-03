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
  const api = useApi();

  // Get project and jobs list first (from job_tree - always up to date)
  const { project, jobs, mutateJobs } = useProject(projectId);

  // Find job status from jobs array (consistent with jobs list icons)
  const jobFromTree = useMemo(() => {
    return jobs?.find((j) => j.id === jobid);
  }, [jobs, jobid]);

  // Get detailed job data (params_xml, container, etc.) - status may lag behind
  const {
    job,
    params_xml,
    mutateParams_xml,
    validation,
    diagnostic_xml,
    def_xml,
    container,
  } = useJob(jobid);

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
        console.log("[JobView] XML saved successfully:", response.data?.message);
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
  const { data: report_xml_json } = api.jobReportXml(job?.id, isJobActive);

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
        <JobMenu />
      </Container>
    </>
  );
};
