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

  const { data: report_xml_json, error: fetchError, mutate: mutateReportXml } = useSWR<any>(
    job ? `/api/proxy/jobs/${job.id}/report_xml/` : null,
    swrFetcher,
    { refreshInterval: job?.status == 3 || job?.status == 2 ? 5000 : 0 }
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
      console.log(`[Report] Received response for job ${job?.id}:`,
        report_xml_json.xml ? `${report_xml_json.xml.length} chars` : 'no xml field',
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

  useEffect(() => {
    if (job && oldJob && job.status !== oldJob.status) {
      if (job.status > 3 && job.id === oldJob.id) {
        setMessage(`Job finished with status ${job.status}`);
        mutateReportXml(() => null); // Force re-fetch
      }
    }
  }, [job, oldJob]);

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
