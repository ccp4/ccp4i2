"use client";
import { use, useEffect, useMemo, useState } from "react";
import { Box, Container, LinearProgress, Tab, Tabs } from "@mui/material";
import { useApi } from "../api";
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

export interface JobViewProps {
  jobid: number;
}
export const JobView: React.FC<JobViewProps> = ({ jobid }) => {
  const { job, params_xml, validation, diagnostic_xml, def_xml, container } =
    useJob(jobid);
  const { project, jobs, mutateJobs } = useProject(job?.project || 0);
  const { devMode, jobId, setJobId, projectId, setProjectId } =
    useCCP4i2Window();
  const {
    setExtraDialogActions,
    setProcessedErrors,
    extraDialogActions,
    processedErrors,
  } = useRunCheck();
  const { mode } = useTheme();

  const previousJob = usePrevious(job);
  const { jobTabValue: tabValue, setJobTabValue: setTabValue } = useJobTab();
  //const [tabValue, setTabValue] = useState<Number>(job?.status == 1 ? 0 : 3);
  const { data: report_xml_json, mutate: mutateReportXml } = useSWR<any>(
    job ? `/api/proxy/jobs/${job.id}/report_xml/` : null,
    swrFetcher,
    { refreshInterval: job?.status == 3 || job?.status == 2 ? 5000 : 0 }
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
      if (job && job != previousJob) {
        setTabValue(job.status == 1 ? 0 : [3, 6].includes(job.status) ? 3 : 4);
      }
    };
    asyncFunc();
  }, [job, setJobId]);

  //Here a useEffect that will clear processedErrors and extraDialogActions when job changes
  useEffect(() => {
    if (!setExtraDialogActions || !setProcessedErrors) return;
    if (extraDialogActions) setExtraDialogActions(null);
    if (processedErrors) setProcessedErrors(null);
  }, [jobid, setExtraDialogActions, setProcessedErrors]);

  return !project || !jobs || !job ? (
    <LinearProgress />
  ) : (
    <Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <ToolBar />
      <Container
        sx={{
          flex: 1,
          display: "flex",
          flexDirection: "column",
          minHeight: 0, // Critical: allows flex children to shrink below content size
          pb: 1,
        }}
      >
        <Box sx={{ flexShrink: 0 }}>
          <JobHeader job={job} mutateJobs={mutateJobs} />
        </Box>
        <Tabs
          value={tabValue}
          onChange={handleTabChange}
          variant="fullWidth"
          sx={{ flexShrink: 0 }}
        >
          <Tab value={0} label="Task interface" />
          {devMode && <Tab value={1} label="Params as xml" />}
          {devMode && <Tab value={2} label="Report as xml" />}
          {(devMode || [3, 4, 6, 7, 9, 10].includes(job?.status)) && (
            <Tab value={3} label="Report" />
          )}
          {(devMode || job?.status === 5) && (
            <Tab value={4} label="Diagnostics" />
          )}
          {devMode && <Tab value={5} label="Def xml" />}
          {(devMode || job?.status === 1) && (
            <Tab value={6} label="Validation report" />
          )}
          {devMode && <Tab value={7} label="Job container" />}
          <Tab value={8} label="Comments" />
          <Tab value={9} label="Directory" />
          <Tab value={10} label="Logs" />
        </Tabs>
        {/* Scroll container fills remaining space via flexbox */}
        <Box
          sx={(theme) => ({
            flex: 1,
            minHeight: 0, // Critical: allows this flex child to shrink and scroll
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
              <TaskContainer />
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
          {(devMode || job?.status === 5) && tabValue == 4 && diagnostic_xml && (
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
          {(devMode || job?.status === 1) && tabValue == 6 && validation && (
            <ValidationViewer job={job} />
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
            <JobCommentEditor jobId={job.id} />
          )}
          {tabValue == 9 && job && project && (
            <JobDirectoryView job={job} project={project} />
          )}
          {tabValue == 10 && job && project && (
            <JobLogViewer job={job} project={project} />
          )}
        </Box>
        <JobMenu />
      </Container>
    </Box>
  );
};
