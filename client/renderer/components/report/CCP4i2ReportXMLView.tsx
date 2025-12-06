import { ReactNode, useCallback, useEffect, useMemo } from "react";
import $ from "jquery";
import {
  Avatar,
  Button,
  Paper,
  Skeleton,
  Stack,
  Typography,
} from "@mui/material";
import { Job } from "../../types/models";
import { CCP4i2ReportElement } from "./CCP4i2ReportElement";
import { useApi } from "../../api";
import { useCCP4i2Window } from "../../app-context";
import { useJob, usePrevious } from "../../utils";
import { useRouter } from "next/navigation";
import { usePopcorn } from "../../providers/popcorn-provider";
import useSWR from "swr";
import { swrFetcher } from "../../api-fetch";
import { useTheme } from "../../theme/theme-provider";
import { CCP4i2WhatNext } from "./CCP4i2WhatNext";
import { useIsJobEffectivelyActive } from "../../providers/recently-started-jobs-context";

export const CCP4i2ReportXMLView = () => {
  const { customColors } = useTheme();
  const api = useApi();
  const { jobId } = useCCP4i2Window();
  const { job } = useJob(jobId);
  const { mutate: mutateJobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: job?.project,
    endpoint: "jobs",
  });

  // Use grace period-aware hook to handle DB latency after job submission
  // This ensures we keep polling even if the job status hasn't updated yet
  const isJobActive = useIsJobEffectivelyActive(job?.id, job?.status);

  const { data: report_xml_json, error: fetchError, mutate: mutateReportXml } = useSWR<any>(
    job ? `/api/proxy/jobs/${job.id}/report_xml/` : null,
    swrFetcher,
    { refreshInterval: isJobActive ? 5000 : 0 }
  );

  // Debug logging
  useEffect(() => {
    if (job) {
      console.log(`[Report] Fetching report for job ${job.id}, status=${job.status}`);
    }
    if (fetchError) {
      console.error(`[Report] Fetch error for job ${job?.id}:`, fetchError);
    }
    if (report_xml_json) {
      // Handle both wrapped response {success: true, data: {xml: ...}} and direct {xml: ...}
      const xmlString = report_xml_json.data?.xml || report_xml_json.xml;
      console.log(`[Report] Received response for job ${job?.id}:`,
        xmlString ? `${xmlString.length} chars` : 'no xml field',
        report_xml_json);
    }
  }, [job, fetchError, report_xml_json]);

  const report_xml: XMLDocument | null = useMemo(() => {
    if (!report_xml_json) return null;
    // Handle both wrapped response {success: true, data: {xml: ...}} and direct {xml: ...}
    const xmlString = report_xml_json.data?.xml || report_xml_json.xml;
    if (!xmlString) return null;
    return $.parseXML(xmlString);
  }, [report_xml_json]);

  const oldJob = usePrevious(job);

  const router = useRouter();

  const { setMessage } = usePopcorn();

  // Detect job status transitions from active to inactive
  useEffect(() => {
    if (!job || !oldJob || job.id !== oldJob.id) return;
    if (job.status === oldJob.status) return;

    // Check if job transitioned from active (running/queued) to inactive (finished/failed)
    const wasActive = [2, 3, 7].includes(oldJob.status);
    const isNowInactive = ![2, 3, 7].includes(job.status);

    if (wasActive && isNowInactive) {
      console.log(`[Report] Job ${job.id} finished (${oldJob.status} -> ${job.status}), refreshing report`);
      setMessage(`Job finished with status ${job.status}`);
      mutateReportXml(); // Force re-fetch
    }
  }, [job, oldJob, setMessage, mutateReportXml]);

  const reportContent = useMemo<ReactNode[] | null>(() => {
    if (!report_xml) return null;
    if (!job) return null;
    return $(report_xml)
      .children()
      .children()
      .map((iItem: number, item: any) => {
        return (
          <CCP4i2ReportElement
            key={`${iItem}`}
            iItem={iItem}
            item={item}
            job={job}
          />
        );
      })
      .toArray();
  }, [report_xml, job]);

  // Show error state
  if (fetchError) {
    return (
      <Paper sx={{ p: 2, m: 2 }}>
        <Typography color="error" variant="h6">
          Failed to load report
        </Typography>
        <Typography variant="body2" color="text.secondary">
          {fetchError.message || String(fetchError)}
        </Typography>
      </Paper>
    );
  }

  return !reportContent ? (
    <Skeleton />
  ) : (
    <>
      <Paper
        sx={{
          width: "100%",
          height: "calc(100vh - 22rem)",
          overflowY: "auto",
        }}
      >
        {reportContent}
        <CCP4i2WhatNext />
      </Paper>
    </>
  );
};
