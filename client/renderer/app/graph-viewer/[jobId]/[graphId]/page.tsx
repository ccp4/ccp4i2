"use client";

import { useMemo, useEffect, useState } from "react";
import { useParams } from "next/navigation";
import { Box, Container, LinearProgress, Typography } from "@mui/material";
import $ from "jquery";
import { useApi } from "../../../../api";
import { CCP4i2ApplicationOutputView } from "../../../../components/report/CCP4i2ApplicationOutputView";

/**
 * Standalone graph viewer page.
 *
 * Opens a graph in its own window, allowing users to view and interact
 * with graphs independently of the main report view.
 *
 * URL: /graph-viewer/[jobId]/[graphId]
 */
export default function GraphViewerPage() {
  const params = useParams();
  const { jobId, graphId } = params as { jobId: string; graphId: string };
  const api = useApi();
  const [title, setTitle] = useState<string>("Graph Viewer");

  // Fetch the report XML for this job
  const { data: report_xml_json, error: fetchError } = api.jobReportXml(
    jobId ? parseInt(jobId) : null,
    false // no polling needed for standalone viewer
  );

  // Parse the XML and find the specific graph element
  const graphElement = useMemo<Element | null>(() => {
    if (!report_xml_json || !graphId) return null;

    // Handle both wrapped response {success: true, data: {xml: ...}} and direct {xml: ...}
    const xmlString = report_xml_json.data?.xml || report_xml_json.xml;
    if (!xmlString) return null;

    try {
      const xmlDoc = $.parseXML(xmlString);
      // Find the graph element by its key attribute (which contains the internalId)
      const $graph = $(xmlDoc).find(`CCP4i2ReportFlotGraph[key="${graphId}"]`);

      if ($graph.length > 0) {
        // Extract title for the window
        const graphTitle = $graph.attr("title");
        if (graphTitle) {
          setTitle(graphTitle);
        }
        // Find the ccp4_data element within the graph
        const dataElement = $graph.find("ccp4\\:ccp4_data, ccp4_data, ns0\\:ccp4_data").get(0);
        return dataElement ?? $graph.get(0) ?? null;
      }
      return null;
    } catch (err) {
      console.error("Failed to parse report XML:", err);
      return null;
    }
  }, [report_xml_json, graphId]);

  // Update document title
  useEffect(() => {
    document.title = title;
  }, [title]);

  if (fetchError) {
    return (
      <Container sx={{ mt: 4 }}>
        <Typography color="error">
          Failed to load graph: {fetchError.message}
        </Typography>
      </Container>
    );
  }

  if (!report_xml_json) {
    return (
      <Box sx={{ width: "100%", mt: 2 }}>
        <LinearProgress />
        <Typography sx={{ mt: 2, textAlign: "center" }}>
          Loading graph...
        </Typography>
      </Box>
    );
  }

  if (!graphElement) {
    return (
      <Container sx={{ mt: 4 }}>
        <Typography color="warning.main">
          Graph not found: {graphId}
        </Typography>
      </Container>
    );
  }

  return (
    <Container maxWidth="lg" sx={{ py: 2 }}>
      <Typography variant="h5" sx={{ mb: 2 }}>
        {title}
      </Typography>
      <Box sx={{ height: "calc(100vh - 120px)" }}>
        <CCP4i2ApplicationOutputView
          output={graphElement}
          height="100%"
        />
      </Box>
    </Container>
  );
}
